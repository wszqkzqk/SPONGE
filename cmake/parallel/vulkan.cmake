set(CPP_DIALECT "CXX")

add_definitions(-DUSE_VULKAN)

find_package(Vulkan REQUIRED)

find_path(
  GLSLANG_INCLUDE_DIR
  NAMES "glslang/Public/ShaderLang.h"
  HINTS "$ENV{CONDA_PREFIX}"
  PATH_SUFFIXES "include")
find_library(
  GLSLANG_LIBRARY
  NAMES "glslang"
  HINTS "$ENV{CONDA_PREFIX}"
  PATH_SUFFIXES "lib")
find_library(
  GLSLANG_SPIRV_LIBRARY
  NAMES "SPIRV"
  HINTS "$ENV{CONDA_PREFIX}"
  PATH_SUFFIXES "lib")
find_library(
  GLSLANG_RESOURCE_LIMITS_LIBRARY
  NAMES "glslang-default-resource-limits"
  HINTS "$ENV{CONDA_PREFIX}"
  PATH_SUFFIXES "lib")
foreach(REQUIRED_VAR GLSLANG_INCLUDE_DIR GLSLANG_LIBRARY GLSLANG_SPIRV_LIBRARY
                     GLSLANG_RESOURCE_LIMITS_LIBRARY)
  if(NOT ${REQUIRED_VAR})
    message(FATAL_ERROR "glslang not found: ${REQUIRED_VAR}")
  endif()
endforeach()

target_include_directories(common_libraries INTERFACE ${GLSLANG_INCLUDE_DIR})
target_link_libraries(
  common_libraries
  INTERFACE Vulkan::Vulkan ${GLSLANG_LIBRARY} ${GLSLANG_SPIRV_LIBRARY}
            ${GLSLANG_RESOURCE_LIMITS_LIBRARY})

find_package(ZLIB REQUIRED)
find_package(LLVM CONFIG REQUIRED)
find_package(Clang CONFIG QUIET)
target_include_directories(common_libraries INTERFACE ${LLVM_INCLUDE_DIRS})
if(WIN32)
  llvm_map_components_to_libnames(
    SPONGE_LLVM_LIBS
    support
    core
    executionengine
    native
    nativecodegen
    orcjit
    runtimedyld
    targetparser)
  target_link_libraries(common_libraries INTERFACE ${SPONGE_LLVM_LIBS})
  target_link_libraries(common_libraries INTERFACE clangFrontendTool)
else()
  target_link_libraries(common_libraries INTERFACE LLVM clang-cpp)
endif()

if(ON_ARM)
  message(STATUS "Use Open Source Math Libraries")
  include("${PROJECT_ROOT_DIR}/cmake/math/open_source.cmake")
else()
  message(STATUS "Use MKL as Math Library")
  include("${PROJECT_ROOT_DIR}/cmake/math/mkl.cmake")
endif()

set(SPONGE_VULKAN_GLSL_DIR
    "${PROJECT_ROOT_DIR}/SPONGE/third_party/vulkan_backend/glsl")
set(SPONGE_VULKAN_GENERATED "${CMAKE_BINARY_DIR}/vulkan_shader_sources.cpp")
file(GLOB SPONGE_VULKAN_SHADER_FILES CONFIGURE_DEPENDS
     "${SPONGE_VULKAN_GLSL_DIR}/*.comp" "${SPONGE_VULKAN_GLSL_DIR}/*.glsl")
add_custom_command(
  OUTPUT ${SPONGE_VULKAN_GENERATED}
  COMMAND
    ${CMAKE_COMMAND} -DGLSL_DIR=${SPONGE_VULKAN_GLSL_DIR}
    -DOUT=${SPONGE_VULKAN_GENERATED} -P
    ${PROJECT_ROOT_DIR}/cmake/utils/embed_vulkan_shaders.cmake
  DEPENDS ${SPONGE_VULKAN_SHADER_FILES}
          ${PROJECT_ROOT_DIR}/cmake/utils/embed_vulkan_shaders.cmake
  COMMENT "Embedding Vulkan GLSL kernels")
set(SPONGE_VULKAN_BACKEND_SOURCES
    ${PROJECT_ROOT_DIR}/SPONGE/third_party/vulkan_backend/vulkan_runtime.cpp
    ${PROJECT_ROOT_DIR}/SPONGE/third_party/vulkan_backend/vulkan_fft.cpp
    ${SPONGE_VULKAN_GENERATED})
set_source_files_properties(${SPONGE_VULKAN_BACKEND_SOURCES}
                            PROPERTIES LANGUAGE CXX)
find_path(
  GLSLANG_C_INTERFACE_INCLUDE_DIR
  NAMES "glslang_c_interface.h"
  HINTS "$ENV{CONDA_PREFIX}"
  PATH_SUFFIXES "include/glslang/Include" "include")
set_source_files_properties(
  ${PROJECT_ROOT_DIR}/SPONGE/third_party/vulkan_backend/vulkan_fft.cpp
  PROPERTIES COMPILE_DEFINITIONS "VKFFT_BACKEND=0"
             INCLUDE_DIRECTORIES
             "${PROJECT_ROOT_DIR}/SPONGE/third_party/vkfft;${PROJECT_ROOT_DIR}/SPONGE/third_party/vkfft/vkFFT;${GLSLANG_C_INTERFACE_INCLUDE_DIR}"
)

if(ON_ARM)
  target_compile_definitions(common_libraries INTERFACE USE_NEON)
else()
  target_compile_options(common_libraries INTERFACE -mavx2)
  target_compile_definitions(common_libraries INTERFACE USE_AVX2)
endif()
