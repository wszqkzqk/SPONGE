#pragma once

#include <vulkan/vulkan.h>

#include <cstddef>
#include <cstdint>

namespace sponge_vk
{
struct Allocation
{
    VkBuffer buffer;
    VkDeviceMemory memory;
    VkDeviceSize size;
};

VkInstance Instance();
VkPhysicalDevice PhysicalDevice();
VkDevice Device();
VkQueue Queue();
uint32_t QueueFamilyIndex();
VkCommandPool CommandPool();

bool HasAtomicFloat();
bool HasFloat64();
uint32_t SubgroupSize();

Allocation* AllocationOf(const void* ptr);
void SubmitAndWait(void* stream);
void EnsureRecording(void* stream);
VkCommandBuffer CurrentCommandBuffer();
const char* LastError();
}  // namespace sponge_vk
