#ifndef __COMMONKERNEL_CUH__
#define __COMMONKERNEL_CUH__


#include <cuComplex.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cufft.h"

/*������ʱ����������*/
__global__ void TBaseGen(double *t, double fs, double t0, int elements);

/*�������Ӻ˺���
������*c = *a + *b
*/
__global__ void ComplexAddKernel(const cuComplex *a, const cuComplex *b, cuComplex *c, int elements);

__global__ void MatAddComplexKernal(cuComplex *Data1, cuComplex Data2, int elements);


/*���������˺���
������*b = *a
*/
__global__ void ComplexCopyKernel(const cuComplex *a, cuComplex *b, int elements);

/*dB
������20*log10(*a)
*/
__global__ void dBKernel(const float *a, float* b, int elements);

/*idB
������pow(10, *a / 20)
*/
__global__ void idBKernel(const float *a, float* b, int elements);

/*��ʵ������תΪ�鲿Ϊ0�ĸ�������
������*b = *(a+0i)
*/
__global__ void ComplexKernel(const float *a, cuComplex* b, int elements);

/*����
������*a = *(a+0i)
*/
__global__ void ConjKernel(cuComplex *data, int elements);

/*������ʵ��ƴ��һ������
������*Res = *(Real+Imagi)
*/
__global__ void ComplexMat(cuComplex *Res, float *Real, float *Imag, int elements);

/*����ת��
����:����ת�ã�����ά����
*/
__global__ void TransposeKernel(cuComplex* in, cuComplex* out, int Rows, int Cols, int Bands);
__global__ void TransposeKernel(float* in, float* out, int Rows, int Cols, int Bands);
__global__ void TransposeKernel(int* in, int* out, int Rows, int Cols, int Bands);

//����ת��
__global__ void CTransposeKernel(cuComplex* in, cuComplex* out, int Rows, int Cols, int Bands);


/*ȡʵ��
������*Real = real(*data��
*/
__global__ void Realkernel(cuComplex *Data, float *Real, int elements);

/*ȡ�鲿
������*Imag = imag(*data��
*/
__global__ void Imagkernel(cuComplex *Data, float *Imag, int elements);

/*��CCMat��abs
*/
__global__ void Abskernel(cuComplex *DataIn, float *DataOut, int elements);

// ���Ե�ͨ�ź���
__global__ void LowPass(cuComplex *Res, int StartPoint, int EndPoint, int elements);

/*��������������
������*Res = *(Real+Imagi)
*/
__global__ void DotMulKernal(cuComplex *Data1, cuComplex *Data2, cuComplex *Res, int elements);
__global__ void DotMul2Kernal(cuComplex *Data1, cuComplex *Data2, int elements);

/*���������������
������*Res = *Data1 + *Data2
*/
__global__ void MatAddKernal(cuComplex *Data1, cuComplex *Data2, cuComplex *Res, int elements);

/*���������������
������*Res = *Data1 + *Data2 + *Data3
*/
__global__ void MatAddKernal(cuComplex* Data1, cuComplex* Data2, cuComplex* Data3, cuComplex* Res, int elements);

/*����������ʵ������
������*Data1.x = *Data1.x * *Data2
*Data1.y = *Data1.y * *Data2
*/
__global__ void MatMulKernal(cuComplex *Data1, float *Data2, int elements);

/*���������float
������*Res = *Data1 / Data2
*/
__global__ void C2FDivKernal(cuComplex *Data1, float Data2, cuComplex *Res, int elements);

__global__ void C2FDiv2Kernal(cuComplex *Data1, float Data2, int elements);

/*���������float
������*Res = *Data1 * Data2
*/
__global__ void C2FMulKernal(cuComplex *Data1, float Data2, cuComplex *Res, int elements);



/*float�����float
������*Res = *Data1 * Data2
*/
__global__ void F2FMulKernal(float *Data1, float Data2, float *Res, int elements);

/*��ȡ���
��������Data
*/
__global__ void SliceKernel(cuComplex *Data, cuComplex *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void SliceKernel(float *Data, float *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void SliceKernel(int *Data, int *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands);


/*��������չ
��������Data��չΪResult���ȣ�������չ��
Data : (Rows, Cols, Bands)
Result : (Rows + abs(AddLength), Cols, Bands)
AddLength : ��չ�������Ҹ���
Value : ��չ��ֵ
*/
__global__ void ExtendKernel(float* Data, float* Result, int AddLength, float Value, unsigned int Rows, unsigned int Cols, unsigned int Bands);

/*��������չ
��������Data��չΪResult���ȣ�������չ��
Data : (Rows, Cols, Bands)
Result : (Rows + abs(AddLength), Cols, Bands)
AddLength : ��չ�������Ҹ���
Value : ��չ��ֵ
*/
__global__ void ExtendKernel(cuComplex* Data, cuComplex* Result, int AddLength, cuComplex Value, unsigned int Rows, unsigned int Cols, unsigned int Bands);


/*��ֵɸѡ
������
Data		: (Rows, Cols, Bands)
Result		: (Rows + abs(AddLength), Cols, Bands)
Th			: ��ֵ
CompareFlag	: �Ƚ�����{'<','>','<=','>=','=='}
����
Th = 3; CompareFlag = '>'
Data	:  4 5 6 1 8 3 2 4
Result	:  1 1 1 0 1 0 0 1
*/
template <typename Type>
__global__ void CompareKernel(Type* Data, bool* Result, Type Th, char CompareFlag, unsigned int elements);
__global__ void CompareKernel(cuComplex* Data, bool* Result, float Th, char CompareFlag, unsigned int elements);
__global__ void CompareKernel(float* Data, int* Result, float Th, char CompareFlag, unsigned int elements);
__global__ void CompareKernel(cuComplex* Data, int* Result, float Th, char CompareFlag, unsigned int elements);

/*
���
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(float *Data, float* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void cuSumKernel(cuComplex* Data, cuComplex* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void cuSumKernel(int *Data, int* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);

//ʵ��pow(a-tile(meam(a),N),2)
__global__ void cuPowATileMean(float *data, float *atilemean, float *sum_mean, float *powatilemean, int N);
//ʵ����һ���ͺ���
__global__ void SqrtTileDiv(float *sum, float *alitemean, float *result, int N);
//����pow(a-tile(meam(a),N),2)
__global__ void cuComplexPowATileMean(cuComplex *data, cuComplex *atilemean, cuComplex *sum_mean, float *powatilemean, int N);
//������һ���ͺ���
__global__ void SqrtTileComplexDiv(float *sum, cuComplex *alitemean, cuComplex *result, int N);

// ������һ���˺���
__global__ void z_score_kernel(cuComplex *data, cuComplex *atiledata, float* PowerResult, unsigned int Rows, unsigned int Cols, unsigned int Bands, cuComplex *Sum_first, float *Sum_second, cuComplex* z_score);


// ���н���3άѭ����λ
__global__ void CirculShiftCol(float *src, float *tar, int shift_num, int Rows, int Cols, int Bands);
__global__ void CirculShiftCol(cuComplex* Data, cuComplex* Result, int shift_num, int Rows, int Cols, int Bands);

// ����2λѭ����λ
__global__ void Circcushift(cuComplex *src, cuComplex *tar, int  Row, int Col, int shift_num);
__global__ void Circcushift(float *src, float *tar, int  Row, int Col, int shift_num);




/*�ź�����ר��*/
__global__ void SetPerTarHn(cuComplex *Data, cuComplex *Temp, int colIndex, int elements);

#endif // !__COMMONKERNEL_CUH__
