#ifndef __COMMONKERNEL_CUH__
#define __COMMONKERNEL_CUH__


#include <cuComplex.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cufft.h"

/*采样段时间序列生成*/
__global__ void TBaseGen(double *t, double fs, double t0, int elements);

/*复数叠加核函数
描述：*c = *a + *b
*/
__global__ void ComplexAddKernel(const cuComplex *a, const cuComplex *b, cuComplex *c, int elements);

__global__ void MatAddComplexKernal(cuComplex *Data1, cuComplex Data2, int elements);


/*复数拷贝核函数
描述：*b = *a
*/
__global__ void ComplexCopyKernel(const cuComplex *a, cuComplex *b, int elements);

/*dB
描述：20*log10(*a)
*/
__global__ void dBKernel(const float *a, float* b, int elements);

/*idB
描述：pow(10, *a / 20)
*/
__global__ void idBKernel(const float *a, float* b, int elements);

/*将实数矩阵转为虚部为0的复数矩阵
描述：*b = *(a+0i)
*/
__global__ void ComplexKernel(const float *a, cuComplex* b, int elements);

/*共轭
描述：*a = *(a+0i)
*/
__global__ void ConjKernel(cuComplex *data, int elements);

/*将两个实数拼成一个复数
描述：*Res = *(Real+Imagi)
*/
__global__ void ComplexMat(cuComplex *Res, float *Real, float *Imag, int elements);

/*矩阵转置
描述:行列转置，第三维不变
*/
__global__ void TransposeKernel(cuComplex* in, cuComplex* out, int Rows, int Cols, int Bands);
__global__ void TransposeKernel(float* in, float* out, int Rows, int Cols, int Bands);
__global__ void TransposeKernel(int* in, int* out, int Rows, int Cols, int Bands);

//共轭转置
__global__ void CTransposeKernel(cuComplex* in, cuComplex* out, int Rows, int Cols, int Bands);


/*取实部
描述：*Real = real(*data）
*/
__global__ void Realkernel(cuComplex *Data, float *Real, int elements);

/*取虚部
描述：*Imag = imag(*data）
*/
__global__ void Imagkernel(cuComplex *Data, float *Imag, int elements);

/*求CCMat的abs
*/
__global__ void Abskernel(cuComplex *DataIn, float *DataOut, int elements);

// 乘以低通门函数
__global__ void LowPass(cuComplex *Res, int StartPoint, int EndPoint, int elements);

/*两个复数矩阵点乘
描述：*Res = *(Real+Imagi)
*/
__global__ void DotMulKernal(cuComplex *Data1, cuComplex *Data2, cuComplex *Res, int elements);
__global__ void DotMul2Kernal(cuComplex *Data1, cuComplex *Data2, int elements);

/*两个复数矩阵相加
描述：*Res = *Data1 + *Data2
*/
__global__ void MatAddKernal(cuComplex *Data1, cuComplex *Data2, cuComplex *Res, int elements);

/*三个复数矩阵相加
描述：*Res = *Data1 + *Data2 + *Data3
*/
__global__ void MatAddKernal(cuComplex* Data1, cuComplex* Data2, cuComplex* Data3, cuComplex* Res, int elements);

/*复数矩阵点乘实数矩阵
描述：*Data1.x = *Data1.x * *Data2
*Data1.y = *Data1.y * *Data2
*/
__global__ void MatMulKernal(cuComplex *Data1, float *Data2, int elements);

/*复数矩阵除float
描述：*Res = *Data1 / Data2
*/
__global__ void C2FDivKernal(cuComplex *Data1, float Data2, cuComplex *Res, int elements);

__global__ void C2FDiv2Kernal(cuComplex *Data1, float Data2, int elements);

/*复数矩阵乘float
描述：*Res = *Data1 * Data2
*/
__global__ void C2FMulKernal(cuComplex *Data1, float Data2, cuComplex *Res, int elements);



/*float矩阵乘float
描述：*Res = *Data1 * Data2
*/
__global__ void F2FMulKernal(float *Data1, float Data2, float *Res, int elements);

/*截取深拷贝
描述：将Data
*/
__global__ void SliceKernel(cuComplex *Data, cuComplex *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void SliceKernel(float *Data, float *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void SliceKernel(int *Data, int *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands);


/*矩阵按列扩展
描述：将Data扩展为Result长度（按列扩展）
Data : (Rows, Cols, Bands)
Result : (Rows + abs(AddLength), Cols, Bands)
AddLength : 扩展量（正右负左）
Value : 扩展的值
*/
__global__ void ExtendKernel(float* Data, float* Result, int AddLength, float Value, unsigned int Rows, unsigned int Cols, unsigned int Bands);

/*矩阵按列扩展
描述：将Data扩展为Result长度（按列扩展）
Data : (Rows, Cols, Bands)
Result : (Rows + abs(AddLength), Cols, Bands)
AddLength : 扩展量（正右负左）
Value : 扩展的值
*/
__global__ void ExtendKernel(cuComplex* Data, cuComplex* Result, int AddLength, cuComplex Value, unsigned int Rows, unsigned int Cols, unsigned int Bands);


/*阈值筛选
描述：
Data		: (Rows, Cols, Bands)
Result		: (Rows + abs(AddLength), Cols, Bands)
Th			: 阈值
CompareFlag	: 比较类型{'<','>','<=','>=','=='}
例：
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
求和
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(float *Data, float* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void cuSumKernel(cuComplex* Data, cuComplex* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);
__global__ void cuSumKernel(int *Data, int* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands);

//实数pow(a-tile(meam(a),N),2)
__global__ void cuPowATileMean(float *data, float *atilemean, float *sum_mean, float *powatilemean, int N);
//实数归一化和函数
__global__ void SqrtTileDiv(float *sum, float *alitemean, float *result, int N);
//复数pow(a-tile(meam(a),N),2)
__global__ void cuComplexPowATileMean(cuComplex *data, cuComplex *atilemean, cuComplex *sum_mean, float *powatilemean, int N);
//复数归一化和函数
__global__ void SqrtTileComplexDiv(float *sum, cuComplex *alitemean, cuComplex *result, int N);

// 复数归一化核函数
__global__ void z_score_kernel(cuComplex *data, cuComplex *atiledata, float* PowerResult, unsigned int Rows, unsigned int Cols, unsigned int Bands, cuComplex *Sum_first, float *Sum_second, cuComplex* z_score);


// 按列进行3维循环移位
__global__ void CirculShiftCol(float *src, float *tar, int shift_num, int Rows, int Cols, int Bands);
__global__ void CirculShiftCol(cuComplex* Data, cuComplex* Result, int shift_num, int Rows, int Cols, int Bands);

// 按行2位循环移位
__global__ void Circcushift(cuComplex *src, cuComplex *tar, int  Row, int Col, int shift_num);
__global__ void Circcushift(float *src, float *tar, int  Row, int Col, int shift_num);




/*信号生成专属*/
__global__ void SetPerTarHn(cuComplex *Data, cuComplex *Temp, int colIndex, int elements);

#endif // !__COMMONKERNEL_CUH__
