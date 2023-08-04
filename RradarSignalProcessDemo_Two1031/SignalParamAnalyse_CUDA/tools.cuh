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
	/*���������еĸ�������*/
	void SaveHostComplexMatrix(cuComplex *h_data, int rows, int cols, string path);

	/*�����豸�еĸ�������*/
	void SaveDeviceComplexMatrix(cuComplex *d_data, int rows, int cols, string path);

	/*�����豸�е�ʵ������*/
	void SaveDeviceFloatMatrix(float* d_data, int rows, int cols, string path);

	/*�����豸�е���������*/
	void SaveDeviceIntMatrix(int *d_data, int rows, int cols, string path);

	/*�����豸�е�����*/
	void SaveDeviceFloatVector(float *d_data, int elements, std::string path);

	/*���������е�short������*/
	void SaveHostShortVector(short* d_data, int elements, std::string path);

	/*���������еĸ�������*/
	void SaveHostComplexVector(cuComplex *h_data, int elements, std::string path);

	/*�����豸�еĸ�������*/
	void SaveDeviceComplexVector(cuComplex *d_data, int elements, std::string path);

	/*�����豸�е���������*/
	void SaveDeviceIntVector(int* d_data, int elements, std::string path);	

	/*�����豸�еĸ�������*/
	template <typename Type>
	void SaveDeviceVector(Type *d_data, int elements, string path);

	/*���������еĸ�������*/
	template <typename Type>
	void SaveHostVector(Type *d_data, int elements, string path);

	void Zeros(CudaArray<cuComplex> &data, unsigned int N, cudaStream_t stream);
	void Zeros(CudaArray<float> &data, unsigned int N, cudaStream_t stream);

	/*dB Result = 20*log10(Data)*/
	FCMat dB(FCMat& Data, cudaStream_t stream);

	/*idB Result = pow(10, Data / 20)*/
	FCMat idB(FCMat& Data, cudaStream_t stream);



	/*��ʵ������תΪ�鲿Ϊ0�ĸ�������*/
	CCMat Complex(FCMat& Data, cudaStream_t stream);
	/*��ʵ������תΪ�鲿Ϊ0�ĸ�������*/
	void Complex(FCMat& Data, CCMat& Result, cudaStream_t stream);

	/*������ʵ������תΪ��������*/
	CCMat Complex(FCMat& Real, FCMat& Imag, cudaStream_t stream);
	/*������ʵ������תΪ��������*/
	void Complex(FCMat& Real, FCMat& Imag, CCMat& Result, cudaStream_t stream);

	/*��������ȡʵ��*/
	FCMat Real(CCMat& Data, cudaStream_t Stream);
	void Real(CCMat& Data, FCMat& Result, cudaStream_t Stream);

	/*��������ȡ�鲿*/
	FCMat Imag(CCMat& Data, cudaStream_t Stream);
	void Imag(CCMat& Data, FCMat& Result, cudaStream_t Stream);

	/*����������*/
	CCMat DotMul(CCMat& Data1, CCMat& Data2, cudaStream_t stream);
	/*����������, ���ΪData1*/
	void StillDotMul(CCMat& Data1, CCMat& Data2, cudaStream_t stream);

	/*��������ӷ�*/
	CCMat MatAdd(CCMat& Data1, CCMat& Data2, cudaStream_t stream);
	/*��������ӷ�*/
	void MatAdd(CCMat& Data1, cuComplex Data2, cudaStream_t stream);

	/*��������ӷ�*/
	void MatAdd(CCMat& Data1, CCMat& Data2, CCMat& Result, cudaStream_t stream);
	/*��������ӷ�*/
	void MatAdd(CCMat& Data1, CCMat& Data2, CCMat& Data3, CCMat& Result, cudaStream_t stream);

	/*float����˷�*/
	FCMat MatMul(FCMat& A, FCMat& B, cudaStream_t stream);
	/*��������˷�*/
	CCMat MatMul(CCMat& A, CCMat& B, cudaStream_t stream);
	/*��������˷�*/
	void MatMul(CCMat& A, CCMat& B, CCMat& C, cudaStream_t stream);

	/*���������float*/
	void StillMatFMul(CCMat& Data1, float Data2, cudaStream_t stream);

	/*���������float*/
	CCMat CFDiv(CCMat& Data1, float Data2, cudaStream_t stream);
	/*ԭ�ظ��������float*/
	void StillCFDiv(CCMat& Data1, float Data2, cudaStream_t stream);

	/*���������float*/
	void MatFMul(CCMat& Data1, FCMat& Data2, cudaStream_t stream);
	/*���������float*/
	CCMat MatFMul(CCMat& Data1, float Data2, cudaStream_t stream);
	/*float�����float*/
	FCMat MatFMul(FCMat& Data1, float Data2, cudaStream_t stream);
	/*���������float*/
	void MatFMul(CCMat& Data1, CCMat& Result, float Data2, cudaStream_t stream);

	/*��ȡ�������Ƭ��*/
	CCMat Slice(CCMat& Data, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream);

	/*��ȡ�������Ƭ��*/
	void Slice(CCMat& Data, CCMat& Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream);
	/*��ȡ���(FCMat)*/
	void Slice(FCMat& Data, FCMat& Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream);


	//// �ۼӺ�
	//FCMat cuAccum(FCMat Data);


	/*��������չ
	��������Data��չΪResult���ȣ�������չ��
	Data : (Rows, Cols, Bands)
	Result : (Rows + abs(AddLength), Cols, Bands)
	AddLength : ��չ�������Ҹ���
	Value : ��չ��ֵ
	*/
	FCMat cuExtend(FCMat& Data, int AddLength, float Value, cudaStream_t stream);
	void cuExtend(CCMat& Data, CCMat& Result, int AddLength, cuComplex Value, cudaStream_t stream);
	CCMat cuExtend(CCMat& Data, int AddLength, cuComplex Value, cudaStream_t stream);

	// ת��
	void Transpose(FCMat& data, cudaStream_t stream);
	void Transpose(CCMat& data, cudaStream_t stream);
	void Transpose(CCMat& data, CCMat& result, cudaStream_t stream);
	// ����ת��
	void CTranspose(CCMat& data, cudaStream_t stream);

	// ת��
	FCMat cuTranspose(FCMat& data, cudaStream_t stream);
	CCMat cuTranspose(CCMat& data, cudaStream_t stream);

	// ����ת��
	CCMat cuCTranspose(CCMat& data, cudaStream_t stream);

	/*���������*/
	FCMat sum(FCMat& Data, cudaStream_t stream);
	/*float���������*/
	void sum(FCMat& Data, FCMat& result, cudaStream_t stream);

	CCMat sum(CCMat& Data, cudaStream_t stream);
	/*�������������*/
	void sum(CCMat& Data, CCMat& result, cudaStream_t stream);

	ICMat sum(ICMat& Data, cudaStream_t stream);
	// ���
	float cuSum(FCMat& Data, cudaStream_t stream);
	// ����abs���
	float cuSum(CCMat& Data, cudaStream_t stream);

	// ��һ��
	void Z_ScoreStandardization(FCMat &data, FCMat &z_score, int N, cudaStream_t stream);
	void Z_ScoreStandardization(CCMat &data, CCMat &z_score, int N, cudaStream_t stream);

	// ������һ��
	void Z_ScoreStandardization_con(CCMat &data, CCMat &z_score, int N, cudaStream_t stream);

	// ��ȡ���ֵ
	float cuMax(FCMat& Data, cudaStream_t stream);

	// ��ȡ���ֵresult���±�index
	void cuMax(FCMat& Data, float& result, dim3& index, cudaStream_t stream);

	// ��ȡ���ֵresult���±�index
	void cuMax(FCMat& Data, FCMat& result, dim3& index, cudaStream_t stream);

	// FFT
	/*��չ��FFTNum�����FFT ����*/
	CCMat cuFFT(CCMat &Data, unsigned int FFTNum, cudaStream_t stream);
	void cuFFT(CCMat &Data, CCMat& DataF, unsigned int FFTNum, cudaStream_t stream);

	// IFFT
	/*IFFT���ٽضϵ�OriginNum ����*/
	void cuIFFT(CCMat &Data, unsigned int OriginNum, cudaStream_t stream);
	void cuIFFT(CCMat &Data, CCMat &Data1, unsigned int OriginNum, int Cols, int Bands, cudaStream_t stream);
	// ��������ģ�CPU������
	// ������
	void cuConj(CCMat& Data, cudaStream_t stream);
	double CenterOfMass(Vector<double> envelope);

	// ����ѭ���ƶ���shiftNum����0�������ƣ���������
	void cuShiftCol(CCMat& Data, CCMat& Result, int shiftNum, cudaStream_t stream);
	// ����ѭ���ƶ���shiftNum����0�������ƣ���������
	void cuShiftCol(FCMat& Data, FCMat& Result, int shiftNum, cudaStream_t stream);
	
	template <typename Type>
	void DisDim(CudaArray<Type>& Data)
	{
		printf("(%d, %d, %d)\n", Data.dims(0), Data.dims(1), Data.dims(2));
	}

}

#endif // _TOOLS_CUH_



