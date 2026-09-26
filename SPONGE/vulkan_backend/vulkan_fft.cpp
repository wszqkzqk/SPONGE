#include "../third_party/device_backend/vulkan_api.h"

#ifdef warpSize
#undef warpSize
#endif

#include <cstdlib>
#include <cstring>
#include <vector>

#include "vulkan_internal.h"
#include "vkFFT/vkFFT.h"

namespace
{
struct SpongeFftPlan
{
    VkFFTApplication app;
    VkFFTConfiguration config;
    VkFence fence = VK_NULL_HANDLE;
    int type = 0;
    int batch = 1;
    uint64_t real_bytes = 0;
    uint64_t complex_bytes = 0;
    uint64_t real_batch_stride = 0;
    uint64_t complex_batch_stride = 0;
};

void FftFail(VkFFTResult result, const char* what)
{
    if (result == VKFFT_SUCCESS) return;
    fprintf(stderr, "[SPONGE Vulkan] %s failed with VkFFTResult %d\n", what,
            (int)result);
    abort();
}
}  // namespace

int deviceFFTPlanMany(FFT_HANDLE* handle, int dimension, FFT_SIZE_t* length,
                      void*, int, int, void*, int, int, FFT_TYPE type,
                      int batch)
{
    auto* plan = new SpongeFftPlan();
    plan->type = type;
    plan->batch = batch;

    uint64_t real_elems = 1, complex_elems = 1;
    for (int i = 0; i < dimension; i++) real_elems *= (uint64_t)length[i];
    complex_elems = real_elems / (uint64_t)length[dimension - 1] *
                    (uint64_t)(length[dimension - 1] / 2 + 1);
    plan->real_bytes = real_elems * sizeof(float);
    plan->complex_bytes = complex_elems * 2 * sizeof(float);
    plan->real_batch_stride = plan->real_bytes;
    plan->complex_batch_stride = plan->complex_bytes;

    VkFFTConfiguration& config = plan->config;
    memset(&config, 0, sizeof(config));
    config.FFTdim = dimension;
    for (int i = 0; i < dimension; i++)
        config.size[i] = (uint64_t)length[dimension - 1 - i];

    VkFenceCreateInfo fence_info{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    if (vkCreateFence(sponge_vk::Device(), &fence_info, nullptr,
                      &plan->fence) != VK_SUCCESS)
    {
        delete plan;
        return 1;
    }

    static VkPhysicalDevice physical_device = sponge_vk::PhysicalDevice();
    static VkDevice device = sponge_vk::Device();
    static VkQueue queue = sponge_vk::Queue();
    static VkCommandPool command_pool = sponge_vk::CommandPool();
    config.physicalDevice = &physical_device;
    config.device = &device;
    config.queue = &queue;
    config.commandPool = &command_pool;
    config.fence = &plan->fence;
    config.performR2C = 1;
    config.normalize = 0;
    config.makeForwardPlanOnly = type == FFT_R2C ? 1 : 0;
    config.makeInversePlanOnly = type == FFT_C2R ? 1 : 0;
    config.numberBatches = 1;
    config.specifyOffsetsAtLaunch = batch > 1 ? 1 : 0;
    config.bufferSize = &plan->complex_bytes;
    config.inputBufferSize = &plan->real_bytes;
    config.isInputFormatted = 1;
    config.inverseReturnToInputBuffer = 1;

    VkFFTResult result = initializeVkFFT(&plan->app, config);
    if (result != VKFFT_SUCCESS)
    {
        vkDestroyFence(device, plan->fence, nullptr);
        delete plan;
        FftFail(result, "initializeVkFFT");
        return 1;
    }
    *handle = plan;
    return 0;
}

static void FftAppend(SpongeFftPlan* plan, void* real_ptr, void* complex_ptr,
                      int inverse)
{
    sponge_vk::EnsureRecording(nullptr);
    VkCommandBuffer command_buffer = sponge_vk::CurrentCommandBuffer();
    VkBuffer real_buffer = sponge_vk::AllocationOf(real_ptr)->buffer;
    VkBuffer complex_buffer = sponge_vk::AllocationOf(complex_ptr)->buffer;
    for (int b = 0; b < plan->batch; b++)
    {
        VkFFTLaunchParams launch_params{};
        launch_params.commandBuffer = &command_buffer;
        launch_params.buffer = &complex_buffer;
        launch_params.inputBuffer = &real_buffer;
        launch_params.bufferOffset = b * plan->complex_batch_stride;
        launch_params.inputBufferOffset = b * plan->real_batch_stride;
        FftFail(VkFFTAppend(&plan->app, inverse, &launch_params),
                "VkFFTAppend");
    }
}

int deviceFFTExecR2C(FFT_HANDLE handle, float* input, FFT_COMPLEX* output)
{
    auto* plan = (SpongeFftPlan*)handle;
    FftAppend(plan, input, output, 0);
    return 0;
}

int deviceFFTExecC2R(FFT_HANDLE handle, FFT_COMPLEX* input, float* output)
{
    auto* plan = (SpongeFftPlan*)handle;
    FftAppend(plan, output, input, 1);
    return 0;
}

void deviceFFTDestroy(FFT_HANDLE handle)
{
    auto* plan = (SpongeFftPlan*)handle;
    if (plan == nullptr) return;
    sponge_vk::SubmitAndWait(nullptr);
    deleteVkFFT(&plan->app);
    vkDestroyFence(sponge_vk::Device(), plan->fence, nullptr);
    delete plan;
}
