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

//实数pow(a-tile(meam(a),N),2)
//__global__ void cuPowATileMean(float *data, float *atilemean, float *sum_mean, float *powatilemean, int N);
////实数归一化和函数
//__global__ void SqrtTileDiv(float *sum, float *alitemean, float *result, int N);
////复数pow(a-tile(meam(a),N),2)
//__global__ void cuComplexPowATileMean(cuComplex *data, cuComplex *atilemean, cuComplex *sum_mean, float *powatilemean, int N);
////复数归一化和函数
//__global__ void SqrtTileComplexDiv(float *sum, cuComplex *alitemean, cuComplex *result, int N);

//pow
__global__ void cuPow(float *data, float *result, int N, float base);
//// 绝对值
__global__ void cuAbs(cuComplex *Data, float *Result, int Len);
//实部+虚部=复数
__global__ void floatToComplex(float *real, float *imag, cuComplex *tar, int len);
//复数=实部+虚部
__global__ void ComplexTofloat(float *real, float *imag, cuComplex *tar, int len);
//对数据进行线性变换 y=scale * x + bias
__global__ void LinearTrans(float *result, float *data, float scale, float bias, int N);
//复数求指数  复数求角度，取值范围为(-pi, pi)
__global__ void ComplexIndex(cuComplex *data, int N, float base);
/*
求和
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(float *Data, float* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);
/*
求和
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(cuComplex* Data, cuComplex* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);



/*  float / （float矩阵 + float） */
__global__ void MatDivKernel(float *data, float *result, float scale, float bias, int N);
//产生float序列，start:step:start+(N-1)*step
__global__ void seqKernel(float *result, float start, float step, int N);
//将列向量复制到一个矩阵中
__global__ void vec2mat(float * data, float *result, int rows, int cols);
//由角度根据欧拉公式得到复信号
__global__ void EulerFormula(float * sata, cuComplex *result, int N);
__global__ void Complexvec2mat(cuComplex * data, cuComplex *result, int rows, int cols);
//矩阵转置
__global__ void transposeKernel(float* in, float* out, int Rows, int Cols);

#endif