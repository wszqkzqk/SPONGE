#include "../third_party/device_backend/vulkan_api.h"

#include <vulkan/vulkan.h>

#include <glslang/Public/ResourceLimits.h>
#include <glslang/Public/ShaderLang.h>
#include <glslang/SPIRV/GlslangToSpv.h>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "vulkan_internal.h"

namespace sponge_vk
{
const char* CommonSource();
int KernelCount();
const char* KernelName(int id);
const char* KernelSource(int id);
int KernelBindingCount(int id);
}  // namespace sponge_vk

namespace
{
struct VulkanState
{
    VkInstance instance = VK_NULL_HANDLE;
    std::vector<VkPhysicalDevice> physical_devices;
    VkPhysicalDevice physical_device = VK_NULL_HANDLE;
    VkDevice device = VK_NULL_HANDLE;
    VkQueue queue = VK_NULL_HANDLE;
    uint32_t queue_family = 0;
    VkCommandPool command_pool = VK_NULL_HANDLE;

    bool atomic_float = false;
    bool float64 = false;
    uint32_t subgroup_size = 32;
    uint32_t max_workgroup_invocations = 1024;
    uint32_t max_push_constants_size = 128;
    uint32_t api_version = 0;

    VkBuffer staging = VK_NULL_HANDLE;
    VkDeviceMemory staging_memory = VK_NULL_HANDLE;
    void* staging_mapped = nullptr;
    VkDeviceSize staging_size = 0;
    VkDeviceSize staging_offset = 0;

    std::unordered_map<const void*, sponge_vk::Allocation> allocations;
    std::unordered_map<std::string, int> kernel_ids;
    std::unordered_map<uint64_t, VkPipeline> pipelines;
    std::unordered_map<int, VkShaderModule> shader_modules;
    std::unordered_map<int, VkDescriptorSetLayout> ds_layouts;
    std::unordered_map<uint64_t, VkPipelineLayout> pipeline_layouts;
    std::string last_error;
    std::recursive_mutex mutex;
};

VulkanState& S()
{
    static VulkanState state;
    return state;
}

static bool SyncEach()
{
    static const bool sync_each = getenv("SPONGE_VK_SYNC_EACH") != nullptr;
    return sync_each;
}

struct StreamState
{
    VkCommandBuffer cb = VK_NULL_HANDLE;
    VkFence fence = VK_NULL_HANDLE;
    VkDescriptorPool pool = VK_NULL_HANDLE;
    bool recording = false;
};

StreamState g_default_stream;
std::unordered_map<void*, StreamState*> g_streams;
std::unordered_map<const void*, sponge_vk::SerialPoolInfo> g_serial_pools;

struct GraphState
{
    VkCommandBuffer cb = VK_NULL_HANDLE;
    VkDescriptorPool pool = VK_NULL_HANDLE;
};

GraphState* g_recording = nullptr;

StreamState* StreamOf(void* stream)
{
    if (stream == nullptr) return &g_default_stream;
    auto it = g_streams.find(stream);
    return it != g_streams.end() ? it->second : &g_default_stream;
}

#define VK_CHECK(call)                                              \
    do                                                              \
    {                                                               \
        VkResult res_ = (call);                                     \
        if (res_ != VK_SUCCESS)                                     \
        {                                                           \
            char buf[256];                                          \
            snprintf(buf, sizeof(buf), "%s failed with VkResult %d", \
                     #call, (int)res_);                             \
            Fail(buf);                                              \
        }                                                           \
    } while (0)

[[noreturn]] void Fail(const std::string& message)
{
    S().last_error = message;
    fprintf(stderr, "[SPONGE Vulkan] %s\n", message.c_str());
    abort();
}

struct StreamState;
void SubmitAndWaitImpl(StreamState* stream);

uint32_t FindMemoryType(uint32_t type_bits, VkMemoryPropertyFlags required,
                        VkMemoryPropertyFlags preferred)
{
    VkPhysicalDeviceMemoryProperties props;
    vkGetPhysicalDeviceMemoryProperties(S().physical_device, &props);
    uint32_t fallback = UINT32_MAX;
    for (uint32_t i = 0; i < props.memoryTypeCount; i++)
    {
        if (!(type_bits & (1u << i))) continue;
        const VkMemoryPropertyFlags flags = props.memoryTypes[i].propertyFlags;
        if ((flags & required) != required) continue;
        if ((flags & preferred) == preferred) return i;
        if (fallback == UINT32_MAX) fallback = i;
    }
    if (fallback != UINT32_MAX) return fallback;
    Fail("no suitable memory type");
}

void EnsureStaging(VkDeviceSize size)
{
    if (size <= S().staging_size) return;
    if (S().staging != VK_NULL_HANDLE)
    {
        SubmitAndWaitImpl(&g_default_stream);
        for (auto& [_, stream] : g_streams) SubmitAndWaitImpl(stream);
        vkDestroyBuffer(S().device, S().staging, nullptr);
        vkFreeMemory(S().device, S().staging_memory, nullptr);
    }
    VkDeviceSize new_size = S().staging_size > 0 ? S().staging_size : 16 << 20;
    while (new_size < size) new_size *= 2;
    VkBufferCreateInfo buffer_info{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    buffer_info.size = new_size;
    buffer_info.usage = VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                        VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    buffer_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    VK_CHECK(vkCreateBuffer(S().device, &buffer_info, nullptr, &S().staging));
    VkMemoryRequirements reqs;
    vkGetBufferMemoryRequirements(S().device, S().staging, &reqs);
    VkMemoryAllocateInfo alloc_info{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    alloc_info.allocationSize = reqs.size;
    alloc_info.memoryTypeIndex =
        FindMemoryType(reqs.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                                VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                       0);
    VK_CHECK(vkAllocateMemory(S().device, &alloc_info, nullptr, &S().staging_memory));
    VK_CHECK(vkBindBufferMemory(S().device, S().staging, S().staging_memory, 0));
    VK_CHECK(vkMapMemory(S().device, S().staging_memory, 0, new_size, 0,
                         &S().staging_mapped));
    S().staging_size = new_size;
    S().staging_offset = 0;
}

void EnsureRecordingImpl(StreamState* stream)
{
    if (stream->recording) return;
    VK_CHECK(vkWaitForFences(S().device, 1, &stream->fence, VK_TRUE, UINT64_MAX));
    VK_CHECK(vkResetFences(S().device, 1, &stream->fence));
    VK_CHECK(vkResetCommandBuffer(stream->cb, 0));
    if (stream->pool != VK_NULL_HANDLE)
        VK_CHECK(vkResetDescriptorPool(S().device, stream->pool, 0));
    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    VK_CHECK(vkBeginCommandBuffer(stream->cb, &begin));
    stream->recording = true;
}

void SubmitAndWaitImpl(StreamState* stream)
{
    if (!stream->recording) return;
    VK_CHECK(vkEndCommandBuffer(stream->cb));
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &stream->cb;
    VK_CHECK(vkQueueSubmit(S().queue, 1, &submit, stream->fence));
    VK_CHECK(vkWaitForFences(S().device, 1, &stream->fence, VK_TRUE, UINT64_MAX));
    stream->recording = false;
    S().staging_offset = 0;
}

void PipelineBarrierCompute(VkCommandBuffer cb)
{
    VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask =
        VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    vkCmdPipelineBarrier(cb, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier,
                         0, nullptr, 0, nullptr);
}

void PipelineBarrierTransferToCompute(StreamState* stream)
{
    VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    barrier.dstAccessMask =
        VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    vkCmdPipelineBarrier(stream->cb, VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier, 0,
                         nullptr, 0, nullptr);
}

void PipelineBarrierComputeToTransfer(StreamState* stream)
{
    VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    vkCmdPipelineBarrier(stream->cb, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 1, &barrier, 0,
                         nullptr, 0, nullptr);
}

std::string BuildShaderDefines()
{
    std::string defines;
    if (S().atomic_float)
        defines += "#extension GL_EXT_shader_atomic_float : enable\n";
    defines +=
        "#extension GL_KHR_shader_subgroup_basic : enable\n"
        "#extension GL_KHR_shader_subgroup_arithmetic : enable\n"
        "#extension GL_KHR_shader_subgroup_ballot : enable\n"
        "#extension GL_KHR_shader_subgroup_shuffle : enable\n";
    if (S().float64) defines += "#extension GL_ARB_gpu_shader_fp64 : enable\n";
    defines += "#define SPONGE_VK_ATOMIC_FLOAT " +
               std::to_string(S().atomic_float ? 1 : 0) + "\n";
    defines += "#define SPONGE_VK_FP64 " +
               std::to_string(S().float64 ? 1 : 0) + "\n";
    defines += "#define SPONGE_VK_SUBGROUP_SIZE " +
               std::to_string(S().subgroup_size) + "\n";
    return defines;
}

VkShaderModule CompileKernel(int kernel_id)
{
    static std::once_flag init_once;
    std::call_once(init_once, []() { glslang::InitializeProcess(); });

    const std::string source = std::string("#version 450\n") +
                               BuildShaderDefines() +
                               sponge_vk::CommonSource() + "\n" +
                               sponge_vk::KernelSource(kernel_id);

    glslang::TShader shader(EShLangCompute);
    const char* strings[] = {source.c_str()};
    shader.setStrings(strings, 1);
    shader.setEnvInput(glslang::EShSourceGlsl, EShLangCompute,
                       glslang::EShClientVulkan, 100);
    shader.setEnvClient(glslang::EShClientVulkan,
                        glslang::EShTargetVulkan_1_1);
    shader.setEnvTarget(glslang::EShTargetSpv, glslang::EShTargetSpv_1_3);
    shader.setEntryPoint("main");
    shader.setSourceEntryPoint("main");

    const TBuiltInResource* resources = GetDefaultResources();
    const EShMessages messages =
        (EShMessages)(EShMsgSpvRules | EShMsgVulkanRules);
    if (!shader.parse(resources, 450, false, messages))
    {
        Fail(std::string("glslang parse failed for ") +
             sponge_vk::KernelName(kernel_id) + ":\n" + shader.getInfoLog() +
             "\n" + shader.getInfoDebugLog());
    }
    glslang::TProgram program;
    program.addShader(&shader);
    if (!program.link(messages))
    {
        Fail(std::string("glslang link failed for ") +
             sponge_vk::KernelName(kernel_id) + ":\n" + program.getInfoLog());
    }
    std::vector<uint32_t> spirv;
    glslang::GlslangToSpv(*program.getIntermediate(EShLangCompute), spirv);

    VkShaderModuleCreateInfo module_info{
        VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    module_info.codeSize = spirv.size() * sizeof(uint32_t);
    module_info.pCode = spirv.data();
    VkShaderModule module;
    VK_CHECK(vkCreateShaderModule(S().device, &module_info, nullptr, &module));
    return module;
}

VkDescriptorSetLayout DescriptorLayoutOf(int kernel_id)
{
    auto it = S().ds_layouts.find(kernel_id);
    if (it != S().ds_layouts.end()) return it->second;
    const int binding_count = sponge_vk::KernelBindingCount(kernel_id);
    std::vector<VkDescriptorSetLayoutBinding> bindings(binding_count);
    for (int i = 0; i < binding_count; i++)
    {
        bindings[i].binding = i;
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo layout_info{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    layout_info.bindingCount = binding_count;
    layout_info.pBindings = bindings.data();
    VkDescriptorSetLayout layout;
    VK_CHECK(vkCreateDescriptorSetLayout(S().device, &layout_info, nullptr,
                                         &layout));
    S().ds_layouts[kernel_id] = layout;
    return layout;
}

VkPipelineLayout PipelineLayoutOf(int kernel_id, uint32_t push_size)
{
    const uint64_t key = (uint64_t)kernel_id << 32 | push_size;
    auto it = S().pipeline_layouts.find(key);
    if (it != S().pipeline_layouts.end()) return it->second;
    VkPipelineLayout layout;
    VkDescriptorSetLayout ds_layout = DescriptorLayoutOf(kernel_id);
    VkPipelineLayoutCreateInfo layout_info{
        VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layout_info.setLayoutCount = 1;
    layout_info.pSetLayouts = &ds_layout;
    VkPushConstantRange range{VK_SHADER_STAGE_COMPUTE_BIT, 0, push_size};
    if (push_size > 0)
    {
        layout_info.pushConstantRangeCount = 1;
        layout_info.pPushConstantRanges = &range;
    }
    VK_CHECK(vkCreatePipelineLayout(S().device, &layout_info, nullptr, &layout));
    S().pipeline_layouts[key] = layout;
    return layout;
}

VkPipeline PipelineOf(int kernel_id, unsigned int bx, unsigned int by,
                      uint32_t push_size)
{
    const uint64_t key = (uint64_t)kernel_id << 40 | (uint64_t)bx << 20 |
                         (uint64_t)by << 4 | (push_size > 0 ? 1ull : 0ull);
    auto it = S().pipelines.find(key);
    if (it != S().pipelines.end()) return it->second;

    auto module_it = S().shader_modules.find(kernel_id);
    if (module_it == S().shader_modules.end())
    {
        VkShaderModule module = CompileKernel(kernel_id);
        module_it = S().shader_modules.emplace(kernel_id, module).first;
    }
    VkPipelineLayout layout = PipelineLayoutOf(kernel_id, push_size);

    const uint32_t spec_data[3] = {bx, by, 1u};
    const VkSpecializationMapEntry spec_entries[3] = {
        {0, 0, sizeof(uint32_t)}, {1, sizeof(uint32_t), sizeof(uint32_t)},
        {2, 2 * sizeof(uint32_t), sizeof(uint32_t)}};
    VkSpecializationInfo spec_info;
    spec_info.mapEntryCount = 3;
    spec_info.pMapEntries = spec_entries;
    spec_info.dataSize = sizeof(spec_data);
    spec_info.pData = spec_data;

    VkComputePipelineCreateInfo pipeline_info{
        VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    pipeline_info.stage.sType =
        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    pipeline_info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    pipeline_info.stage.module = module_it->second;
    pipeline_info.stage.pName = "main";
    pipeline_info.stage.pSpecializationInfo = &spec_info;
    pipeline_info.layout = layout;
    VkPipeline pipeline;
    VK_CHECK(vkCreateComputePipelines(S().device, VK_NULL_HANDLE, 1,
                                      &pipeline_info, nullptr, &pipeline));
    S().pipelines[key] = pipeline;
    return pipeline;
}

StreamState* CreateStream()
{
    auto* stream = new StreamState();
    VkCommandBufferAllocateInfo cb_info{
        VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cb_info.commandPool = S().command_pool;
    cb_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cb_info.commandBufferCount = 1;
    VK_CHECK(vkAllocateCommandBuffers(S().device, &cb_info, &stream->cb));
    VkFenceCreateInfo fence_info{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
    VK_CHECK(vkCreateFence(S().device, &fence_info, nullptr, &stream->fence));
    VkDescriptorPoolSize pool_size{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 4096};
    VkDescriptorPoolCreateInfo pool_info{
        VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    pool_info.maxSets = 1024;
    pool_info.poolSizeCount = 1;
    pool_info.pPoolSizes = &pool_size;
    VK_CHECK(vkCreateDescriptorPool(S().device, &pool_info, nullptr,
                                    &stream->pool));
    return stream;
}

void DestroyStream(StreamState* stream)
{
    if (stream == nullptr) return;
    SubmitAndWaitImpl(stream);
    vkDestroyDescriptorPool(S().device, stream->pool, nullptr);
    vkDestroyFence(S().device, stream->fence, nullptr);
    vkFreeCommandBuffers(S().device, S().command_pool, 1, &stream->cb);
    delete stream;
}

void CopyHostToDevice(void* dst, const void* src, size_t size,
                      StreamState* stream)
{
    EnsureStaging(size);
    if (S().staging_offset + size > S().staging_size)
    {
        SubmitAndWaitImpl(stream);
    }
    memcpy((char*)S().staging_mapped + S().staging_offset, src, size);
    EnsureRecordingImpl(stream);
    VkBufferCopy region{S().staging_offset, 0, size};
    vkCmdCopyBuffer(stream->cb, S().staging,
                    sponge_vk::AllocationOf(dst)->buffer, 1, &region);
    PipelineBarrierTransferToCompute(stream);
    S().staging_offset += (size + 255) & ~VkDeviceSize(255);
}

void CopyDeviceToHost(void* dst, const void* src, size_t size,
                      StreamState* stream)
{
    SubmitAndWaitImpl(&g_default_stream);
    for (auto& [_, other] : g_streams) SubmitAndWaitImpl(other);
    const auto* alloc = sponge_vk::AllocationOf(src);
    EnsureStaging(size);
    EnsureRecordingImpl(&g_default_stream);
    PipelineBarrierComputeToTransfer(&g_default_stream);
    VkBufferCopy region{0, 0, size};
    vkCmdCopyBuffer(g_default_stream.cb, alloc->buffer, S().staging, 1,
                    &region);
    PipelineBarrierTransferToCompute(&g_default_stream);
    SubmitAndWaitImpl(&g_default_stream);
    memcpy(dst, S().staging_mapped, size);
}

void CopyDeviceToHostAsync(void* dst, const void* src, size_t size,
                           StreamState* stream)
{
    SubmitAndWaitImpl(&g_default_stream);
    for (auto& [_, other] : g_streams) SubmitAndWaitImpl(other);
    const auto* alloc = sponge_vk::AllocationOf(src);
    EnsureStaging(size);
    EnsureRecordingImpl(stream);
    PipelineBarrierComputeToTransfer(stream);
    VkBufferCopy region{0, 0, size};
    vkCmdCopyBuffer(stream->cb, alloc->buffer, S().staging, 1, &region);
    PipelineBarrierTransferToCompute(stream);
    SubmitAndWaitImpl(stream);
    memcpy(dst, S().staging_mapped, size);
}

void CopyDeviceToDevice(void* dst, const void* src, size_t size,
                        StreamState* stream)
{
    EnsureRecordingImpl(stream);
    PipelineBarrierComputeToTransfer(stream);
    VkBufferCopy region{0, 0, size};
    vkCmdCopyBuffer(stream->cb, sponge_vk::AllocationOf(src)->buffer,
                    sponge_vk::AllocationOf(dst)->buffer, 1, &region);
    PipelineBarrierTransferToCompute(stream);
}

}  // namespace

namespace sponge_vk
{
VkInstance Instance() { return S().instance; }
VkPhysicalDevice PhysicalDevice() { return S().physical_device; }
VkDevice Device() { return S().device; }
VkQueue Queue() { return S().queue; }
uint32_t QueueFamilyIndex() { return S().queue_family; }
VkCommandPool CommandPool() { return S().command_pool; }
bool HasAtomicFloat() { return S().atomic_float; }
bool HasFloat64() { return S().float64; }
uint32_t SubgroupSize() { return S().subgroup_size; }
const char* LastError() { return S().last_error.c_str(); }

Allocation* AllocationOf(const void* ptr)
{
    auto it = S().allocations.find(ptr);
    if (it == S().allocations.end())
        Fail("unknown device pointer passed to Vulkan backend");
    return &it->second;
}

void RegisterSerialPool(const void* nl, const void* pool, int stride)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    g_serial_pools[nl] = SerialPoolInfo{pool, stride};
}

SerialPoolInfo SerialPoolOf(const void* nl)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    auto it = g_serial_pools.find(nl);
    if (it == g_serial_pools.end())
        Fail("neighbor list has no registered serial pool");
    return it->second;
}

size_t AllocationSize(const void* ptr)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    return AllocationOf(ptr)->size;
}

void EnsureRecording(void* stream) { EnsureRecordingImpl(StreamOf(stream)); }
void SubmitAndWait(void* stream) { SubmitAndWaitImpl(StreamOf(stream)); }
VkCommandBuffer CurrentCommandBuffer() { return g_default_stream.cb; }

void HostBarrier()
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    SubmitAndWaitImpl(&g_default_stream);
    for (auto& [_, stream] : g_streams) SubmitAndWaitImpl(stream);
}

int KernelId(const char* name)
{
    auto it = S().kernel_ids.find(name);
    if (it != S().kernel_ids.end()) return it->second;
    for (int i = 0; i < KernelCount(); i++)
    {
        if (strcmp(KernelName(i), name) == 0)
        {
            S().kernel_ids[name] = i;
            return i;
        }
    }
    Fail(std::string("unknown Vulkan kernel: ") + name);
}

void Launch(int kernel_id, unsigned int grid_x, unsigned int grid_y,
            unsigned int block_x, unsigned int block_y,
            const void* const* buffers, int buffer_count, const void* params,
            size_t params_size, deviceStream_t stream_handle)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    if (params_size > S().max_push_constants_size)
        Fail(std::string("push constant size exceeded for ") +
             KernelName(kernel_id));
    StreamState* stream = StreamOf(stream_handle);
    VkCommandBuffer cb = stream->cb;
    VkDescriptorPool pool = stream->pool;
    if (g_recording != nullptr)
    {
        cb = g_recording->cb;
        pool = g_recording->pool;
    }
    else
    {
        EnsureRecordingImpl(stream);
    }

    const uint32_t push_size = (uint32_t)params_size;
    VkPipeline pipeline = PipelineOf(kernel_id, block_x, block_y, push_size);
    VkDescriptorSetLayout ds_layout = DescriptorLayoutOf(kernel_id);
    VkPipelineLayout pipeline_layout = PipelineLayoutOf(kernel_id, push_size);

    VkDescriptorSetAllocateInfo set_info{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    set_info.descriptorPool = pool;
    set_info.descriptorSetCount = 1;
    set_info.pSetLayouts = &ds_layout;
    VkDescriptorSet descriptor_set;
    VK_CHECK(vkAllocateDescriptorSets(S().device, &set_info, &descriptor_set));

    std::vector<VkWriteDescriptorSet> writes(buffer_count);
    std::vector<VkDescriptorBufferInfo> buffer_infos(buffer_count);
    for (int i = 0; i < buffer_count; i++)
    {
        const Allocation* alloc = AllocationOf(buffers[i]);
        buffer_infos[i].buffer = alloc->buffer;
        buffer_infos[i].offset = 0;
        buffer_infos[i].range = alloc->size;
        writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[i].dstSet = descriptor_set;
        writes[i].dstBinding = i;
        writes[i].descriptorCount = 1;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo = &buffer_infos[i];
    }
    vkUpdateDescriptorSets(S().device, buffer_count, writes.data(), 0, nullptr);

    vkCmdBindPipeline(cb, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
    vkCmdBindDescriptorSets(cb, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout,
                            0, 1, &descriptor_set, 0, nullptr);
    if (params_size > 0)
        vkCmdPushConstants(cb, pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0,
                           (uint32_t)params_size, params);
    vkCmdDispatch(cb, grid_x, grid_y, 1);
    PipelineBarrierCompute(cb);
    if (g_recording == nullptr && SyncEach()) SubmitAndWaitImpl(stream);
}

void* GraphCreate()
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    if (SyncEach() || getenv("SPONGE_VK_NO_GRAPH") != nullptr) return nullptr;
    auto* graph = new GraphState();
    VkCommandBufferAllocateInfo cb_info{
        VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cb_info.commandPool = S().command_pool;
    cb_info.level = VK_COMMAND_BUFFER_LEVEL_SECONDARY;
    cb_info.commandBufferCount = 1;
    VK_CHECK(vkAllocateCommandBuffers(S().device, &cb_info, &graph->cb));
    VkDescriptorPoolSize pool_size{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 4096};
    VkDescriptorPoolCreateInfo pool_info{
        VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    pool_info.maxSets = 1024;
    pool_info.poolSizeCount = 1;
    pool_info.pPoolSizes = &pool_size;
    VK_CHECK(vkCreateDescriptorPool(S().device, &pool_info, nullptr,
                                    &graph->pool));
    return graph;
}

void GraphBeginRecord(void* graph_ptr, deviceStream_t stream_handle)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    auto* graph = static_cast<GraphState*>(graph_ptr);
    SubmitAndWaitImpl(StreamOf(stream_handle));
    VK_CHECK(vkResetCommandBuffer(graph->cb, 0));
    VK_CHECK(vkResetDescriptorPool(S().device, graph->pool, 0));
    VkCommandBufferInheritanceInfo inheritance{
        VK_STRUCTURE_TYPE_COMMAND_BUFFER_INHERITANCE_INFO};
    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    begin.pInheritanceInfo = &inheritance;
    VK_CHECK(vkBeginCommandBuffer(graph->cb, &begin));
    g_recording = graph;
}

void GraphEndRecord(void* graph_ptr)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    VK_CHECK(vkEndCommandBuffer(static_cast<GraphState*>(graph_ptr)->cb));
    g_recording = nullptr;
}

void GraphExecute(void* graph_ptr, deviceStream_t stream_handle)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    auto* graph = static_cast<GraphState*>(graph_ptr);
    StreamState* stream = StreamOf(stream_handle);
    EnsureRecordingImpl(stream);
    vkCmdExecuteCommands(stream->cb, 1, &graph->cb);
    PipelineBarrierCompute(stream->cb);
    if (SyncEach()) SubmitAndWaitImpl(stream);
}

}  // namespace sponge_vk

int deviceInit(unsigned int)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    if (S().instance != VK_NULL_HANDLE) return 0;

    uint32_t api_version = VK_API_VERSION_1_1;
    if (vkEnumerateInstanceVersion != nullptr)
        vkEnumerateInstanceVersion(&api_version);

    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "SPONGE";
    app.apiVersion = api_version;

    std::vector<const char*> extensions;
#ifdef __APPLE__
    extensions.push_back(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
#endif
    VkInstanceCreateInfo instance_info{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
#ifdef __APPLE__
    instance_info.flags = VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
#endif
    instance_info.pApplicationInfo = &app;
    instance_info.enabledExtensionCount = (uint32_t)extensions.size();
    instance_info.ppEnabledExtensionNames = extensions.data();
    if (vkCreateInstance(&instance_info, nullptr, &S().instance) != VK_SUCCESS)
    {
        S().last_error = "vkCreateInstance failed";
        return -1;
    }

    uint32_t count = 0;
    vkEnumeratePhysicalDevices(S().instance, &count, nullptr);
    if (count == 0)
    {
        S().last_error = "no Vulkan physical device found";
        return -1;
    }
    S().physical_devices.resize(count);
    vkEnumeratePhysicalDevices(S().instance, &count,
                               S().physical_devices.data());
    S().api_version = api_version;
    return 0;
}

int deviceGetDeviceCount(int* count)
{
    *count = (int)S().physical_devices.size();
    return 0;
}

int getDeviceProperties(deviceProp* prop, int device_index)
{
    VkPhysicalDevice physical = S().physical_devices[device_index];
    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(physical, &props);
    memset(prop, 0, sizeof(*prop));
    snprintf(prop->name, sizeof(prop->name), "%s", props.deviceName);
    prop->major = VK_API_VERSION_MAJOR(props.apiVersion);
    prop->minor = VK_API_VERSION_MINOR(props.apiVersion);
    prop->maxThreadsPerBlock =
        (int)props.limits.maxComputeWorkGroupInvocations;
    VkPhysicalDeviceMemoryProperties mem_props;
    vkGetPhysicalDeviceMemoryProperties(physical, &mem_props);
    size_t total = 0;
    for (uint32_t i = 0; i < mem_props.memoryHeapCount; i++)
        if (mem_props.memoryHeaps[i].flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT)
            total += mem_props.memoryHeaps[i].size;
    prop->totalGlobalMem = total > 0 ? total : mem_props.memoryHeaps[0].size;

    VkPhysicalDeviceSubgroupProperties subgroup{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
    VkPhysicalDeviceProperties2 props2{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
    props2.pNext = &subgroup;
    vkGetPhysicalDeviceProperties2(physical, &props2);
    prop->warp_size = (int)subgroup.subgroupSize;
    return 0;
}

int setWorkingDevice(int device_index)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    S().physical_device = S().physical_devices[device_index];

    uint32_t family_count = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(S().physical_device, &family_count,
                                             nullptr);
    std::vector<VkQueueFamilyProperties> families(family_count);
    vkGetPhysicalDeviceQueueFamilyProperties(S().physical_device, &family_count,
                                             families.data());
    uint32_t family = UINT32_MAX;
    for (uint32_t i = 0; i < family_count; i++)
    {
        if (families[i].queueFlags & VK_QUEUE_COMPUTE_BIT)
        {
            family = i;
            break;
        }
    }
    if (family == UINT32_MAX) Fail("no compute queue family");
    S().queue_family = family;

    uint32_t extension_count = 0;
    vkEnumerateDeviceExtensionProperties(S().physical_device, nullptr,
                                         &extension_count, nullptr);
    std::vector<VkExtensionProperties> available(extension_count);
    vkEnumerateDeviceExtensionProperties(S().physical_device, nullptr,
                                         &extension_count, available.data());
    auto has_extension = [&](const char* name)
    {
        for (const auto& ext : available)
            if (strcmp(ext.extensionName, name) == 0) return true;
        return false;
    };
    std::vector<const char*> extensions;
    const bool portability =
        has_extension("VK_KHR_portability_subset");
    if (portability) extensions.push_back("VK_KHR_portability_subset");
    S().atomic_float = has_extension(VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME);
    if (S().atomic_float)
        extensions.push_back(VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME);

    VkPhysicalDeviceFeatures2 features2{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    VkPhysicalDeviceShaderAtomicFloatFeaturesEXT atomic_features{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT};
    VkPhysicalDeviceSubgroupProperties subgroup{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
    VkPhysicalDeviceProperties2 props2{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
    props2.pNext = &subgroup;
    vkGetPhysicalDeviceProperties2(S().physical_device, &props2);
    features2.pNext = S().atomic_float ? &atomic_features : nullptr;
    vkGetPhysicalDeviceFeatures2(S().physical_device, &features2);
    if (S().atomic_float && !atomic_features.shaderBufferFloat32AtomicAdd)
        S().atomic_float = false;
    if (getenv("SPONGE_VK_NO_ATOMIC_FLOAT") != nullptr) S().atomic_float = false;
    if (!(subgroup.supportedStages & VK_SHADER_STAGE_COMPUTE_BIT) ||
        !(subgroup.supportedOperations & VK_SUBGROUP_FEATURE_BASIC_BIT) ||
        !(subgroup.supportedOperations & VK_SUBGROUP_FEATURE_ARITHMETIC_BIT) ||
        !(subgroup.supportedOperations & VK_SUBGROUP_FEATURE_BALLOT_BIT) ||
        !(subgroup.supportedOperations & VK_SUBGROUP_FEATURE_SHUFFLE_BIT))
        Fail("device lacks required subgroup operations");

    S().float64 = features2.features.shaderFloat64 == VK_TRUE;
    S().subgroup_size = subgroup.subgroupSize;
    S().max_workgroup_invocations =
        props2.properties.limits.maxComputeWorkGroupInvocations;
    S().max_push_constants_size = props2.properties.limits.maxPushConstantsSize;

    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queue_info{
        VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    queue_info.queueFamilyIndex = family;
    queue_info.queueCount = 1;
    queue_info.pQueuePriorities = &priority;

    VkPhysicalDeviceFeatures enabled_features{};
    enabled_features.shaderFloat64 = S().float64 ? VK_TRUE : VK_FALSE;
    enabled_features.robustBufferAccess = VK_TRUE;
    VkPhysicalDeviceShaderAtomicFloatFeaturesEXT enabled_atomic{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT};
    enabled_atomic.shaderBufferFloat32AtomicAdd = VK_TRUE;

    VkDeviceCreateInfo device_info{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    device_info.queueCreateInfoCount = 1;
    device_info.pQueueCreateInfos = &queue_info;
    device_info.pEnabledFeatures = &enabled_features;
    device_info.enabledExtensionCount = (uint32_t)extensions.size();
    device_info.ppEnabledExtensionNames = extensions.data();
    if (S().atomic_float) device_info.pNext = &enabled_atomic;
    VK_CHECK(vkCreateDevice(S().physical_device, &device_info, nullptr,
                            &S().device));
    vkGetDeviceQueue(S().device, family, 0, &S().queue);

    VkCommandPoolCreateInfo pool_info{
        VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pool_info.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pool_info.queueFamilyIndex = family;
    VK_CHECK(vkCreateCommandPool(S().device, &pool_info, nullptr,
                                 &S().command_pool));

    StreamState* default_stream = CreateStream();
    g_default_stream = *default_stream;
    delete default_stream;
    return 0;
}

deviceError_t deviceMalloc(void** ptr, size_t size)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    *ptr = nullptr;
    if (size == 0) size = 4;
    VkBufferCreateInfo buffer_info{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    buffer_info.size = size;
    buffer_info.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
                        VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                        VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    buffer_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    VkBuffer buffer;
    if (vkCreateBuffer(S().device, &buffer_info, nullptr, &buffer) != VK_SUCCESS)
        return 1;
    VkMemoryRequirements reqs;
    vkGetBufferMemoryRequirements(S().device, buffer, &reqs);
    VkMemoryAllocateInfo alloc_info{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    alloc_info.allocationSize = reqs.size;
    alloc_info.memoryTypeIndex =
        FindMemoryType(reqs.memoryTypeBits,
                       VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 0);
    VkDeviceMemory memory;
    if (vkAllocateMemory(S().device, &alloc_info, nullptr, &memory) !=
        VK_SUCCESS)
    {
        vkDestroyBuffer(S().device, buffer, nullptr);
        return 1;
    }
    VK_CHECK(vkBindBufferMemory(S().device, buffer, memory, 0));
    void* token = (void*)buffer;
    S().allocations.emplace(token,
                            sponge_vk::Allocation{buffer, memory, reqs.size});
    *ptr = token;
    return 0;
}

void deviceFree(void* ptr)
{
    if (ptr == nullptr) return;
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    sponge_vk::HostBarrier();
    g_serial_pools.erase(ptr);
    auto it = S().allocations.find(ptr);
    if (it == S().allocations.end()) return;
    vkDestroyBuffer(S().device, it->second.buffer, nullptr);
    vkFreeMemory(S().device, it->second.memory, nullptr);
    S().allocations.erase(it);
}

deviceError_t deviceMemcpy(void* to, const void* from, size_t size,
                           deviceMemcpyKind kind)
{
    if (size == 0) return 0;
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    StreamState* stream = &g_default_stream;
    switch (kind)
    {
        case deviceMemcpyHostToDevice:
            CopyHostToDevice(to, from, size, stream);
            break;
        case deviceMemcpyDeviceToHost:
            CopyDeviceToHost(to, from, size, stream);
            break;
        case deviceMemcpyDeviceToDevice:
            CopyDeviceToDevice(to, from, size, stream);
            break;
        default:
            memcpy(to, from, size);
            break;
    }
    return 0;
}

deviceError_t deviceMemcpyAsync(void* to, const void* from, size_t size,
                                deviceMemcpyKind kind,
                                deviceStream_t stream_handle)
{
    if (size == 0) return 0;
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    StreamState* stream = StreamOf(stream_handle);
    switch (kind)
    {
        case deviceMemcpyHostToDevice:
            CopyHostToDevice(to, from, size, stream);
            break;
        case deviceMemcpyDeviceToHost:
            CopyDeviceToHostAsync(to, from, size, stream);
            break;
        case deviceMemcpyDeviceToDevice:
            CopyDeviceToDevice(to, from, size, stream);
            break;
        default:
            memcpy(to, from, size);
            break;
    }
    return 0;
}

void deviceMemset(void* to, int val, size_t size)
{
    if (size == 0) return;
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    StreamState* stream = &g_default_stream;
    EnsureRecordingImpl(stream);
    PipelineBarrierComputeToTransfer(stream);
    vkCmdFillBuffer(stream->cb, sponge_vk::AllocationOf(to)->buffer, 0, size,
                    (uint32_t)val);
    PipelineBarrierTransferToCompute(stream);
    if (SyncEach()) SubmitAndWaitImpl(stream);
}

void deviceStreamCreate(deviceStream_t* stream)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    StreamState* state = CreateStream();
    g_streams[state] = state;
    *stream = state;
}

void deviceStreamDestroy(deviceStream_t stream)
{
    if (stream == nullptr) return;
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    auto it = g_streams.find(stream);
    if (it == g_streams.end()) return;
    DestroyStream(it->second);
    g_streams.erase(it);
}

void deviceStreamSynchronize(deviceStream_t stream)
{
    std::lock_guard<std::recursive_mutex> lock(S().mutex);
    SubmitAndWaitImpl(StreamOf(stream));
}

void hostDeviceSynchronize() { sponge_vk::HostBarrier(); }

const char* deviceGetErrorName(deviceError_t) { return "vulkanError"; }
const char* deviceGetErrorString(deviceError_t)
{
    return S().last_error.c_str();
}
deviceError_t deviceGetLastError() { return 0; }
