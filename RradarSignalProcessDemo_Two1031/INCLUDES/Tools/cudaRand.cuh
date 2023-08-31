#ifndef __CUDARAND_CUH__
#define __CUDARAND_CUH__

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>
#include <curand_kernel.h>
#include <time.h>
#include "CudaArray.cuh"

//// 高斯分布， 一维数据
//void cuGuassRand(float* GRand, unsigned int N);

// 生成高斯分布随机数序列
CCMat Randn(unsigned int Rows, cudaStream_t stream);
// 生成高斯分布随机数序列
void Randn(CCMat& data, unsigned int Rows, cudaStream_t stream);

// 生成高斯分布随机数序列
CCMat Randn(unsigned int Rows, unsigned int Cols, cudaStream_t stream);
// 生成高斯分布随机数序列
void Randn(CCMat& data, unsigned int Rows, unsigned int Cols, cudaStream_t stream);

// 生成高斯分布随机数序列
CCMat Randn(unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

// 生成高斯分布随机数序列
void Randn(CCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

//生成实数高斯分布随机数序列
void Randn(FCMat& result, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

//服从(bias, float+bias)的均匀分布
FCMat Randu(unsigned int Rows, float scale, float bias, cudaStream_t stream);
//服从(bias, float+bias)的均匀分布
void Randu(FCMat& data, unsigned int Rows, float scale, float bias, cudaStream_t stream);

//服从(bias, float+bias)的均匀分布
FCMat Randu(unsigned int Rows, unsigned int Cols, float scale, float bias, cudaStream_t stream);

//服从(bias, float+bias)的均匀分布
void Randu(FCMat& data, unsigned int Rows, unsigned int Cols, float scale, float bias, cudaStream_t stream);


//服从(bias, float+bias)的均匀分布
FCMat Randu(unsigned int Rows, unsigned int Cols, unsigned int Bands, float scale, float bias, cudaStream_t stream);

//服从(bias, float+bias)的均匀分布
void Randu(FCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, float scale, float bias, cudaStream_t stream);


#endif 
