#include "cuPeriodicExtension.cuh"
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"





__global__ void cuPeriodicExtensionKernel(cuComplex* Data, cuComplex* Result, int Period, int Th, int DataLen, int ResultLen)
{
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < ResultLen)
	{
		unsigned int index = idx % DataLen;
		unsigned int index2 = idx % Period;

		if (index2 > Th)
		{
			Result[idx] = Data[index];
		}
		else
		{
			Result[idx].x = 0.0;
			Result[idx].y = 0.0;
		}
	}
}

/*
��������
data: �������һ������
result:�������غ������
Period:ѭ������������Ӧ�շ����ڵĵ���
Th:������Ч������Ӧ����ʱ��ĵ���

example:
Period: 4
Th: 2
data:		1 2 3 4 5 6 1 2 3 4 5 6
��������:	1 1 0 0 1 1 0 0 1 1 0 0 
result:		0 0 3 4 0 0 1 2 0 0 5 6 
*/
void PeriodicExtension(CCMat& Data, CCMat& Result, int Period, int Th, cudaStream_t Stream)
{

	auto DataLen = Data.elements();
	auto ResultLen = Result.elements();


	//// �����߳̾���ָ��GPU�̷߳���
	int threadsPerBlock = 1024;
	unsigned int blocksPerGrid = (ResultLen + threadsPerBlock - 1) / threadsPerBlock;

	cuPeriodicExtensionKernel << < blocksPerGrid, threadsPerBlock, 0, Stream >> >(Data.FirstAddr(), Result.FirstAddr(), Period, Th, DataLen, ResultLen);

}




