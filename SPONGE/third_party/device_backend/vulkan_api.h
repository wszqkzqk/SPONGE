#ifndef BASIC_BACKEND_H
#define BASIC_BACKEND_H
#define VULKAN_ARCH_NAME "Vulkan"

#include <cstddef>
#include <cstdint>

#include "../philox.hpp"
#define Philox4_32_10_t Philox4x32_10
#define device_rand_init(seed, id, offset, state_ptr) \
    (state_ptr)[0] = Philox4_32_10_t(seed, id, offset)
#define device_get_4_normal_distributed_random_numbers(rand_float4,   \
                                                       rand_state, i) \
    rand_state[i].normal((float*)rand_float4, i)
#define __device__
#define __host__
#define __global__
#define __forceinline__ inline
#define __noinline__
#define __launch_bounds__(THREAD)
#if defined(_MSC_VER) && !defined(__restrict__)
#define __restrict__ __restrict
#endif

struct float4
{
    float x, y, z, w;
};
float rnorm3df(float, float, float);
float norm3df(float, float, float);
float erfcxf(float);
float atomicAdd(float*, float);
double atomicAdd(double*, double);
int atomicAdd(int*, int);
int atomicExch(int* address, int val);

struct dim3
{
    unsigned int x;
    unsigned int y;
    unsigned int z;
    dim3(unsigned int ux, unsigned int uy = 1u, unsigned int uz = 1u)
        : x(ux), y(uy), z(uz) {};
};

#define warpSize 0

#define DEVICE_INIT_SUCCESS 0
#define DEVICE_MALLOC_SUCCESS 0

struct deviceProp
{
    char name[256];
    size_t totalGlobalMem;
    int maxThreadsPerBlock;
    int warp_size;
    int major;
    int minor;
};

int deviceInit(unsigned int flags);
int deviceGetDeviceCount(int* count);
int getDeviceProperties(deviceProp* prop, int device);
int setWorkingDevice(int device);

typedef int deviceError_t;
const char* deviceGetErrorName(deviceError_t error);
const char* deviceGetErrorString(deviceError_t error);
deviceError_t deviceGetLastError();
#define deviceErrorLaunchOutOfResources 7
#define deviceErrorInvalidValue 1
#define deviceErrorInvalidConfiguration 2

deviceError_t deviceMalloc(void** ptr, size_t size);
void deviceFree(void* ptr);

enum deviceMemcpyKind
{
    deviceMemcpyHostToHost,
    deviceMemcpyHostToDevice,
    deviceMemcpyDeviceToHost,
    deviceMemcpyDeviceToDevice,
    deviceMemcpyDefault
};
deviceError_t deviceMemcpy(void* to, const void* from, size_t size,
                           deviceMemcpyKind kind);

typedef void* deviceStream_t;
void deviceStreamCreate(deviceStream_t* stream);
void deviceStreamDestroy(deviceStream_t stream);
void deviceStreamSynchronize(deviceStream_t stream);
deviceError_t deviceMemcpyAsync(void* to, const void* from, size_t size,
                                deviceMemcpyKind kind, deviceStream_t stream);
void deviceMemset(void* to, int val, size_t size);
void hostDeviceSynchronize();

namespace sponge_vk
{
void HostBarrier();
void Launch(int kernel_id, unsigned int grid_x, unsigned int grid_y,
            unsigned int block_x, unsigned int block_y,
            const void* const* buffers, int buffer_count, const void* params,
            size_t params_size, deviceStream_t stream);
int KernelId(const char* name);
void* GraphCreate();
void GraphBeginRecord(void* graph, deviceStream_t stream);
void GraphEndRecord(void* graph);
void GraphExecute(void* graph, deviceStream_t stream);
}  // namespace sponge_vk

#define Launch_Device_Kernel(kernel, grid, block, sm_memory, stream, ...) \
    (sponge_vk::HostBarrier(), kernel(__VA_ARGS__))

#define VK_LAUNCH(kernel_name, grid_x, grid_y, block_x, block_y, buffers, \
                  params, stream)                                        \
    sponge_vk::Launch(sponge_vk::KernelId(#kernel_name), grid_x, grid_y, \
                      block_x, block_y, buffers, sizeof(buffers) / sizeof(void*), \
                      params, sizeof(*(params)), stream)

#endif  // BASIC_BACKEND_H

#ifndef FFT_BACKEND_H
#define FFT_BACKEND_H

struct _fft_complex
{
    float x, y;
};
#define FFT_COMPLEX _fft_complex
#define REAL(c) c.x
#define IMAGINARY(c) c.y
#define FFT_LIBRARY_NAME "VkFFT"
typedef void* FFT_HANDLE;
#define FFT_SUCCESS 0
#define FFT_RESULT int
#define FFT_SIZE_t int

enum FFT_TYPE
{
    FFT_R2C,
    FFT_C2R
};

int deviceFFTPlanMany(FFT_HANDLE* handle, int dimension, FFT_SIZE_t* length,
                      void* inembed, int istride, int idist, void* onembed,
                      int ostride, int odist, FFT_TYPE type, int batch);
int deviceFFTExecR2C(FFT_HANDLE handle, float* input, FFT_COMPLEX* output);
int deviceFFTExecC2R(FFT_HANDLE handle, FFT_COMPLEX* input, float* output);
void deviceFFTDestroy(FFT_HANDLE handle);

#endif  // FFT_BACKEND_H

#ifndef BLAS_BACKEND_H
#define BLAS_BACKEND_H

#ifdef USE_MKL
#include <mkl.h>
#define BLAS_LIBRARY_NAME "MKL-BLAS"
#elif defined(USE_OPENBLAS)
#include <cblas.h>
#include <lapacke.h>
#define BLAS_LIBRARY_NAME "OpenBLAS"
#endif

#define BLAS_HANDLE int
#define BLAS_SUCCESS 0

enum deviceBlasOperation_t
{
    DEVICE_BLAS_OP_N,
    DEVICE_BLAS_OP_T,
    DEVICE_BLAS_OP_C
};

enum deviceFillMode_t
{
    DEVICE_FILL_MODE_UPPER
};

enum deviceEigMode_t
{
    DEVICE_EIG_MODE_VECTOR
};

#define deviceBlasCreate(handle)
#define deviceBlasDestroy(handle)

#if defined(USE_MKL) || defined(USE_OPENBLAS)
#define deviceBlasSgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, \
                        ldb, C, ldc)                                          \
    do                                                                        \
    {                                                                         \
        for (int i = 0; i < (m) * (n); ++i)                                   \
            (C)[i] = (*(alpha)) * (A)[i] + (*(beta)) * (B)[i];                \
    } while (0)

#define deviceBlasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B,   \
                        ldb, beta, C, ldc)                                   \
    cblas_sgemm(CblasColMajor,                                               \
                (transa == DEVICE_BLAS_OP_N ? CblasNoTrans : CblasTrans),    \
                (transb == DEVICE_BLAS_OP_N ? CblasNoTrans : CblasTrans), m, \
                n, k, *(alpha), A, lda, B, ldb, *(beta), C, ldc)

#define deviceBlasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B,   \
                        ldb, beta, C, ldc)                                   \
    cblas_dgemm(CblasColMajor,                                               \
                (transa == DEVICE_BLAS_OP_N ? CblasNoTrans : CblasTrans),    \
                (transb == DEVICE_BLAS_OP_N ? CblasNoTrans : CblasTrans), m, \
                n, k, *(alpha), A, lda, B, ldb, *(beta), C, ldc)

#define deviceBlasDdot(handle, n, x, incx, y, incy, result) \
    (*(result) = cblas_ddot(n, x, incx, y, incy), 0)
#endif

#endif  // BLAS_BACKEND_H

#ifndef SOLVER_BACKEND_H
#define SOLVER_BACKEND_H

#ifdef USE_MKL
#include <mkl.h>
#define SOLVER_LIBRARY_NAME "MKL-SOLVER"
#elif defined(USE_OPENBLAS)
#include <lapacke.h>
#define SOLVER_LIBRARY_NAME "LAPACKE"
#endif

#define SOLVER_HANDLE int
#define SOLVER_SUCCESS 0

#define deviceSolverCreate(handle)
#define deviceSolverDestroy(handle)

#if defined(USE_MKL) || defined(USE_OPENBLAS)
#define deviceSolverDsyevdBufferSize(handle, jobz, uplo, n, A, lda, W, lwork) \
    [&]() -> int                                                              \
    {                                                                         \
        double wq;                                                            \
        lapack_int iwq;                                                       \
        LAPACKE_dsyevd_work(LAPACK_COL_MAJOR, 'V', 'U', (lapack_int)(n), (A), \
                            (lapack_int)(lda), (W), &wq, -1, &iwq, -1);       \
        *(lwork) = (int)(wq + 0.5);                                           \
        return 0;                                                             \
    }()

#define deviceSolverDsyevd(handle, jobz, uplo, n, A, lda, W, work, lwork,     \
                           info)                                              \
    do                                                                        \
    {                                                                         \
        lapack_int _liw = 0;                                                  \
        double _wq;                                                           \
        lapack_int _iwq;                                                      \
        LAPACKE_dsyevd_work(LAPACK_COL_MAJOR, 'V', 'U', (lapack_int)(n), (A), \
                            (lapack_int)(lda), (W), &_wq, -1, &_iwq, -1);     \
        _liw = _iwq;                                                          \
        std::vector<lapack_int> _iwork(_liw);                                 \
        *(info) = (int)LAPACKE_dsyevd_work(                                   \
            LAPACK_COL_MAJOR, 'V', 'U', (lapack_int)(n), (A),                 \
            (lapack_int)(lda), (W), (work), (lapack_int)(lwork),              \
            _iwork.data(), (lapack_int)_liw);                                 \
    } while (0)
#endif

#endif  // SOLVER_BACKEND_H
