#ifndef _CUDISTRIBUTEKERNEL_CUH_
#define _CUDISTRIBUTEKERNEL_CUH_

#include <cmath>
#include <cuComplex.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void GaussKernel(float sigmaf, float *hf, int elements);
__global__ void ExponentialKernel(float sigmaf, float *hf, int elements);
__global__ void CauchyKernel(float sigmaf, float *hf, int elements);
__global__ void FullSpectrumKernel(float sigmaf, float *hf, int SpectrumPara, int elements);
__global__ void Swerlling12Kernel(float *Uniformx, float* y, double sigmac, int elements);

//ʵ��pow(a-tile(meam(a),N),2)
//__global__ void cuPowATileMean(float *data, float *atilemean, float *sum_mean, float *powatilemean, int N);
////ʵ����һ���ͺ���
//__global__ void SqrtTileDiv(float *sum, float *alitemean, float *result, int N);
////����pow(a-tile(meam(a),N),2)
//__global__ void cuComplexPowATileMean(cuComplex *data, cuComplex *atilemean, cuComplex *sum_mean, float *powatilemean, int N);
////������һ���ͺ���
//__global__ void SqrtTileComplexDiv(float *sum, cuComplex *alitemean, cuComplex *result, int N);

//pow
__global__ void cuPow(float *data, float *result, int N, float base);
//// ����ֵ
__global__ void cuAbs(cuComplex *Data, float *Result, int Len);
//ʵ��+�鲿=����
__global__ void floatToComplex(float *real, float *imag, cuComplex *tar, int len);
//����=ʵ��+�鲿
__global__ void ComplexTofloat(float *real, float *imag, cuComplex *tar, int len);
//�����ݽ������Ա任 y=scale * x + bias
__global__ void LinearTrans(float *result, float *data, float scale, float bias, int N);
//������ָ��  ������Ƕȣ�ȡֵ��ΧΪ(-pi, pi)
__global__ void ComplexIndex(cuComplex *data, int N, float base);
/*
���
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(float *Data, float* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);
/*
���
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(cuComplex* Data, cuComplex* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);



/*  float / ��float���� + float�� */
__global__ void MatDivKernel(float *data, float *result, float scale, float bias, int N);
//����float���У�start:step:start+(N-1)*step
__global__ void seqKernel(float *result, float start, float step, int N);
//�����������Ƶ�һ��������
__global__ void vec2mat(float * data, float *result, int rows, int cols);
//�ɽǶȸ���ŷ����ʽ�õ����ź�
__global__ void EulerFormula(float * sata, cuComplex *result, int N);
__global__ void Complexvec2mat(cuComplex * data, cuComplex *result, int rows, int cols);
//����ת��
__global__ void transposeKernel(float* in, float* out, int Rows, int Cols);

#endif