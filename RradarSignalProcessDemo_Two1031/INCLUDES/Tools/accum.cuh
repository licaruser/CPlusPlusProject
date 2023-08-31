#ifndef _ACCUM_CUH__
#define _ACCUM_CUH__


#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>
#include <curand_kernel.h>
#include <time.h>
#include "CudaArray.cuh"


// 前缀和
FCMat accum(FCMat& Data, cudaStream_t stream);

// 前缀和
void accum(FCMat& Data, FCMat& Result, cudaStream_t stream);


#endif // !_ACCUM_CUH__
