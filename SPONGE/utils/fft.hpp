/*
CUFFT for cuda backend - license: Nvidia SDK License
MKL for x86 backend - license: Intel Simplified Software License
KML for kunpeng backend - license: 鲲鹏应用使能套件BoostKit用户许可协议 2.0
ARMPL for general arm backend - license: Arm Simplified End User License
Agreement HCFFT for hpcc backend -  cooperated with 沐曦MetaX
*/

#pragma once
#include <string.h>

#ifdef USE_HIP
#include "../third_party/device_backend/hip_api.h"
#elif defined(USE_CUDA)
#include "../third_party/device_backend/cuda_api.h"
#elif defined(USE_VULKAN)
#include "../third_party/device_backend/vulkan_api.h"
#ifdef USE_VULKAN
#include <fftw3.h>
#endif
#else
#include "../third_party/device_backend/cpu_api.h"
#endif

__device__ __host__ __forceinline__ FFT_COMPLEX expc(FFT_COMPLEX z);
__device__ __host__ __forceinline__ FFT_COMPLEX divc(FFT_COMPLEX a,
                                                     FFT_COMPLEX b);

struct SPONGE_FFT_WRAPPER
{
    friend __device__ __host__ __forceinline__ FFT_COMPLEX expc(FFT_COMPLEX z)
    {
        FFT_COMPLEX res;
        float t = expf(REAL(z));
#ifdef _WIN32
        IMAGINARY(res) = sinf(IMAGINARY(z));
        REAL(res) = cosf(IMAGINARY(z));
#else
        sincosf(IMAGINARY(z), &IMAGINARY(res), &REAL(res));
#endif
        REAL(res) *= t;
        IMAGINARY(res) *= t;
        return res;
    }

    friend __device__ __host__ __forceinline__ FFT_COMPLEX divc(FFT_COMPLEX a,
                                                                FFT_COMPLEX b)
    {
        FFT_COMPLEX result;
        float denom = REAL(b) * REAL(b) + IMAGINARY(b) * IMAGINARY(b);
        REAL(result) =
            (REAL(a) * REAL(b) + IMAGINARY(a) * IMAGINARY(b)) / denom;
        IMAGINARY(result) =
            (IMAGINARY(a) * REAL(b) - REAL(a) * IMAGINARY(b)) / denom;
        return result;
    }

    static bool Use_Host_FFT()
    {
#if defined(USE_VULKAN) && !defined(USE_GPU)
        static const bool host_fft = getenv("SPONGE_VK_HOST_FFT") != nullptr;
        return host_fft;
#else
        return false;
#endif
    }

    static FFT_RESULT Make_FFT_Plan(FFT_HANDLE* handle, int batch,
                                    int dimension, FFT_SIZE_t* length,
                                    FFT_TYPE type)
    {
#if defined(USE_GPU) || defined(USE_VULKAN)
        if (!Use_Host_FFT())
            return deviceFFTPlanMany(handle, dimension, length, NULL, 0, 0,
                                     NULL, 0, 0, type, batch);
#endif
#if !defined(USE_GPU) || defined(USE_VULKAN)
        int* c_length = (int*)malloc(sizeof(int) * dimension);
        memcpy(c_length, length, sizeof(int) * dimension);
        c_length[dimension - 1] = c_length[dimension - 1] / 2 + 1;
        int r_dim = 1, c_dim = 1;
        for (int i = 0; i < dimension; i++)
        {
            r_dim *= length[i];
            c_dim *= c_length[i];
        }
        float* tmp_in = (float*)fftwf_malloc(sizeof(float) * r_dim * batch);
        fftwf_complex* tmp_out =
            (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * c_dim * batch);

        if (type == FFT_R2C)
        {
            handle[0] = fftwf_plan_many_dft_r2c(
                dimension, length, batch, tmp_in, length, 1, r_dim, tmp_out,
                c_length, 1, c_dim, FFTW_ESTIMATE);
        }
        else
        {
            handle[0] = fftwf_plan_many_dft_c2r(
                dimension, length, batch, tmp_out, c_length, 1, c_dim, tmp_in,
                length, 1, r_dim, FFTW_ESTIMATE);
        }
        fftwf_free(tmp_in);
        fftwf_free(tmp_out);
        free(c_length);
        return handle[0] == NULL;
#endif
    }

    static void R2C(FFT_HANDLE handle, float* input, FFT_COMPLEX* output)
    {
#if defined(USE_GPU) || defined(USE_VULKAN)
        if (!Use_Host_FFT())
        {
            deviceFFTExecR2C(handle, input, output);
            return;
        }
#endif
        fftwf_execute_dft_r2c((fftwf_plan)handle, input, (fftwf_complex*)output);
    }

    static void C2R(FFT_HANDLE handle, FFT_COMPLEX* input, float* output)
    {
#if defined(USE_GPU) || defined(USE_VULKAN)
        if (!Use_Host_FFT())
        {
            deviceFFTExecC2R(handle, input, output);
            return;
        }
#endif
        fftwf_execute_dft_c2r((fftwf_plan)handle, (fftwf_complex*)input, output);
    }

    static void Destroy_FFT_Plan(FFT_HANDLE* handle)
    {
#if defined(USE_GPU) || defined(USE_VULKAN)
        if (!Use_Host_FFT())
        {
            deviceFFTDestroy(handle[0]);
            return;
        }
#endif
        fftwf_destroy_plan((fftwf_plan)handle[0]);
    }
};
