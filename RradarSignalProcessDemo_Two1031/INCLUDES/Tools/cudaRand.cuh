#ifndef __CUDARAND_CUH__
#define __CUDARAND_CUH__

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>
#include <curand_kernel.h>
#include <time.h>
#include "CudaArray.cuh"

//// ��˹�ֲ��� һά����
//void cuGuassRand(float* GRand, unsigned int N);

// ���ɸ�˹�ֲ����������
CCMat Randn(unsigned int Rows, cudaStream_t stream);
// ���ɸ�˹�ֲ����������
void Randn(CCMat& data, unsigned int Rows, cudaStream_t stream);

// ���ɸ�˹�ֲ����������
CCMat Randn(unsigned int Rows, unsigned int Cols, cudaStream_t stream);
// ���ɸ�˹�ֲ����������
void Randn(CCMat& data, unsigned int Rows, unsigned int Cols, cudaStream_t stream);

// ���ɸ�˹�ֲ����������
CCMat Randn(unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

// ���ɸ�˹�ֲ����������
void Randn(CCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

//����ʵ����˹�ֲ����������
void Randn(FCMat& result, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

//����(bias, float+bias)�ľ��ȷֲ�
FCMat Randu(unsigned int Rows, float scale, float bias, cudaStream_t stream);
//����(bias, float+bias)�ľ��ȷֲ�
void Randu(FCMat& data, unsigned int Rows, float scale, float bias, cudaStream_t stream);

//����(bias, float+bias)�ľ��ȷֲ�
FCMat Randu(unsigned int Rows, unsigned int Cols, float scale, float bias, cudaStream_t stream);

//����(bias, float+bias)�ľ��ȷֲ�
void Randu(FCMat& data, unsigned int Rows, unsigned int Cols, float scale, float bias, cudaStream_t stream);


//����(bias, float+bias)�ľ��ȷֲ�
FCMat Randu(unsigned int Rows, unsigned int Cols, unsigned int Bands, float scale, float bias, cudaStream_t stream);

//����(bias, float+bias)�ľ��ȷֲ�
void Randu(FCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, float scale, float bias, cudaStream_t stream);


#endif 
