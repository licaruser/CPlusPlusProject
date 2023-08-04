#ifndef _TOOLS_CUH_
#define _TOOLS_CUH_

#include <iostream>
#include <fstream>
#include <string>

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "CudaArray.cuh"
#include "SigGeneratePublicDefinition.h"
#include "cufft.h"
#include "accum.cuh"
#include "Error.cuh"

using namespace std;

namespace tools
{
	/*保存主机中的复数矩阵*/
	void SaveHostComplexMatrix(cuComplex *h_data, int rows, int cols, string path);

	/*保存设备中的复数矩阵*/
	void SaveDeviceComplexMatrix(cuComplex *d_data, int rows, int cols, string path);

	/*保存设备中的实数矩阵*/
	void SaveDeviceFloatMatrix(float* d_data, int rows, int cols, string path);

	/*保存设备中的整数矩阵*/
	void SaveDeviceIntMatrix(int *d_data, int rows, int cols, string path);

	/*保存设备中的向量*/
	void SaveDeviceFloatVector(float *d_data, int elements, std::string path);

	/*保存主机中的short型向量*/
	void SaveHostShortVector(short* d_data, int elements, std::string path);

	/*保存主机中的复数向量*/
	void SaveHostComplexVector(cuComplex *h_data, int elements, std::string path);

	/*保存设备中的复数向量*/
	void SaveDeviceComplexVector(cuComplex *d_data, int elements, std::string path);

	/*保存设备中的整数向量*/
	void SaveDeviceIntVector(int* d_data, int elements, std::string path);	

	/*保存设备中的复数矩阵*/
	template <typename Type>
	void SaveDeviceVector(Type *d_data, int elements, string path);

	/*保存主机中的复数矩阵*/
	template <typename Type>
	void SaveHostVector(Type *d_data, int elements, string path);

	void Zeros(CudaArray<cuComplex> &data, unsigned int N, cudaStream_t stream);
	void Zeros(CudaArray<float> &data, unsigned int N, cudaStream_t stream);

	/*dB Result = 20*log10(Data)*/
	FCMat dB(FCMat& Data, cudaStream_t stream);

	/*idB Result = pow(10, Data / 20)*/
	FCMat idB(FCMat& Data, cudaStream_t stream);



	/*将实数矩阵转为虚部为0的复数矩阵*/
	CCMat Complex(FCMat& Data, cudaStream_t stream);
	/*将实数矩阵转为虚部为0的复数矩阵*/
	void Complex(FCMat& Data, CCMat& Result, cudaStream_t stream);

	/*将两个实数矩阵转为复数矩阵*/
	CCMat Complex(FCMat& Real, FCMat& Imag, cudaStream_t stream);
	/*将两个实数矩阵转为复数矩阵*/
	void Complex(FCMat& Real, FCMat& Imag, CCMat& Result, cudaStream_t stream);

	/*复数矩阵取实部*/
	FCMat Real(CCMat& Data, cudaStream_t Stream);
	void Real(CCMat& Data, FCMat& Result, cudaStream_t Stream);

	/*复数矩阵取虚部*/
	FCMat Imag(CCMat& Data, cudaStream_t Stream);
	void Imag(CCMat& Data, FCMat& Result, cudaStream_t Stream);

	/*复数矩阵点乘*/
	CCMat DotMul(CCMat& Data1, CCMat& Data2, cudaStream_t stream);
	/*复数矩阵点乘, 结果为Data1*/
	void StillDotMul(CCMat& Data1, CCMat& Data2, cudaStream_t stream);

	/*复数矩阵加法*/
	CCMat MatAdd(CCMat& Data1, CCMat& Data2, cudaStream_t stream);
	/*复数矩阵加法*/
	void MatAdd(CCMat& Data1, cuComplex Data2, cudaStream_t stream);

	/*复数矩阵加法*/
	void MatAdd(CCMat& Data1, CCMat& Data2, CCMat& Result, cudaStream_t stream);
	/*复数矩阵加法*/
	void MatAdd(CCMat& Data1, CCMat& Data2, CCMat& Data3, CCMat& Result, cudaStream_t stream);

	/*float矩阵乘法*/
	FCMat MatMul(FCMat& A, FCMat& B, cudaStream_t stream);
	/*复数矩阵乘法*/
	CCMat MatMul(CCMat& A, CCMat& B, cudaStream_t stream);
	/*复数矩阵乘法*/
	void MatMul(CCMat& A, CCMat& B, CCMat& C, cudaStream_t stream);

	/*复数矩阵乘float*/
	void StillMatFMul(CCMat& Data1, float Data2, cudaStream_t stream);

	/*复数矩阵除float*/
	CCMat CFDiv(CCMat& Data1, float Data2, cudaStream_t stream);
	/*原地复数矩阵除float*/
	void StillCFDiv(CCMat& Data1, float Data2, cudaStream_t stream);

	/*复数矩阵乘float*/
	void MatFMul(CCMat& Data1, FCMat& Data2, cudaStream_t stream);
	/*复数矩阵乘float*/
	CCMat MatFMul(CCMat& Data1, float Data2, cudaStream_t stream);
	/*float矩阵乘float*/
	FCMat MatFMul(FCMat& Data1, float Data2, cudaStream_t stream);
	/*复数矩阵乘float*/
	void MatFMul(CCMat& Data1, CCMat& Result, float Data2, cudaStream_t stream);

	/*截取深拷贝（切片）*/
	CCMat Slice(CCMat& Data, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream);

	/*截取深拷贝（切片）*/
	void Slice(CCMat& Data, CCMat& Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream);
	/*截取深拷贝(FCMat)*/
	void Slice(FCMat& Data, FCMat& Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream);


	//// 累加和
	//FCMat cuAccum(FCMat Data);


	/*矩阵按列扩展
	描述：将Data扩展为Result长度（按列扩展）
	Data : (Rows, Cols, Bands)
	Result : (Rows + abs(AddLength), Cols, Bands)
	AddLength : 扩展量（正右负左）
	Value : 扩展的值
	*/
	FCMat cuExtend(FCMat& Data, int AddLength, float Value, cudaStream_t stream);
	void cuExtend(CCMat& Data, CCMat& Result, int AddLength, cuComplex Value, cudaStream_t stream);
	CCMat cuExtend(CCMat& Data, int AddLength, cuComplex Value, cudaStream_t stream);

	// 转置
	void Transpose(FCMat& data, cudaStream_t stream);
	void Transpose(CCMat& data, cudaStream_t stream);
	void Transpose(CCMat& data, CCMat& result, cudaStream_t stream);
	// 共轭转置
	void CTranspose(CCMat& data, cudaStream_t stream);

	// 转置
	FCMat cuTranspose(FCMat& data, cudaStream_t stream);
	CCMat cuTranspose(CCMat& data, cudaStream_t stream);

	// 共轭转置
	CCMat cuCTranspose(CCMat& data, cudaStream_t stream);

	/*矩阵按列求和*/
	FCMat sum(FCMat& Data, cudaStream_t stream);
	/*float矩阵按列求和*/
	void sum(FCMat& Data, FCMat& result, cudaStream_t stream);

	CCMat sum(CCMat& Data, cudaStream_t stream);
	/*复数矩阵按列求和*/
	void sum(CCMat& Data, CCMat& result, cudaStream_t stream);

	ICMat sum(ICMat& Data, cudaStream_t stream);
	// 求和
	float cuSum(FCMat& Data, cudaStream_t stream);
	// 复数abs求和
	float cuSum(CCMat& Data, cudaStream_t stream);

	// 归一化
	void Z_ScoreStandardization(FCMat &data, FCMat &z_score, int N, cudaStream_t stream);
	void Z_ScoreStandardization(CCMat &data, CCMat &z_score, int N, cudaStream_t stream);

	// 复数归一化
	void Z_ScoreStandardization_con(CCMat &data, CCMat &z_score, int N, cudaStream_t stream);

	// 获取最大值
	float cuMax(FCMat& Data, cudaStream_t stream);

	// 获取最大值result及下标index
	void cuMax(FCMat& Data, float& result, dim3& index, cudaStream_t stream);

	// 获取最大值result及下标index
	void cuMax(FCMat& Data, FCMat& result, dim3& index, cudaStream_t stream);

	// FFT
	/*扩展到FFTNum后进行FFT 按列*/
	CCMat cuFFT(CCMat &Data, unsigned int FFTNum, cudaStream_t stream);
	void cuFFT(CCMat &Data, CCMat& DataF, unsigned int FFTNum, cudaStream_t stream);

	// IFFT
	/*IFFT后再截断到OriginNum 按列*/
	void cuIFFT(CCMat &Data, unsigned int OriginNum, cudaStream_t stream);
	void cuIFFT(CCMat &Data, CCMat &Data1, unsigned int OriginNum, int Cols, int Bands, cudaStream_t stream);
	// 求面积中心（CPU函数）
	// 矩阵共轭
	void cuConj(CCMat& Data, cudaStream_t stream);
	double CenterOfMass(Vector<double> envelope);

	// 按列循环移动，shiftNum大于0，则右移，否则左移
	void cuShiftCol(CCMat& Data, CCMat& Result, int shiftNum, cudaStream_t stream);
	// 按列循环移动，shiftNum大于0，则右移，否则左移
	void cuShiftCol(FCMat& Data, FCMat& Result, int shiftNum, cudaStream_t stream);
	
	template <typename Type>
	void DisDim(CudaArray<Type>& Data)
	{
		printf("(%d, %d, %d)\n", Data.dims(0), Data.dims(1), Data.dims(2));
	}

}

#endif // _TOOLS_CUH_



