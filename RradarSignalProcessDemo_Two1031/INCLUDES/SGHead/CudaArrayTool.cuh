#ifndef _CUDAARRAYTOOLS_CUH_
#define _CUDAARRAYTOOLS_CUH_

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/*获取矩阵中某个位置的值*/
template <typename Type>
__global__ void getRowKernel(Type *a, Type *result, int index, int elements);

/*设置矩阵中的一行*/
__global__ void setRowKernel(cuComplex *a, cuComplex *Value, int OneRow, int CurrentBand, int Row, int Col, int Band);
__global__ void setRowKernel(float *a, float *Value, int OneRow, int CurrentBand, int Row, int Col, int Band);
__global__ void setRowKernel(int *a, int *Value, int OneRow, int CurrentBand, int Row, int Col, int Band);

/*设置矩阵中的一列*/
__global__ void setColKernel(cuComplex *a, cuComplex *Value, int OneCol, int CurrentBand, int Row, int Col, int Band);
__global__ void setColKernel(float *a, float *Value, int OneCol, int CurrentBand, int Row, int Col, int Band);
__global__ void setColKernel(int *a, int *Value, int OneCol, int CurrentBand, int Row, int Col, int Band);

/*设置矩阵中的一个阵*/
template <typename Type>
__global__ void setBandKernel(cuComplex *a, cuComplex *Value, int OneBand, int Row, int Col, int Band);
__global__ void setBandKernel(float *a, float *Value, int OneBand, int Row, int Col, int Band);
__global__ void setBandKernel(int *a, int *Value, int OneBand, int Row, int Col, int Band);

/*赋值核函数
描述：*a = elements
*/
__global__ void ValuateKernel(cuComplex *a, cuComplex Value, unsigned int Elements);
__global__ void ValuateKernel(float *a, float Value, unsigned int Elements);
__global__ void ValuateKernel(int *a, int Value, unsigned int Elements);
__global__ void ValuateKernel(bool *a, bool Value, unsigned int Elements);


/*设置某行，不开辟新的空间*/
void cuSetRow(cuComplex *Data, cuComplex* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneRow, cudaStream_t stream);
void cuSetRow(float *Data, float* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneRow, cudaStream_t stream);
void cuSetRow(int *Data, int* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneRow, cudaStream_t stream);

/*设置某列，不开辟新的空间*/
void cuSetCol(cuComplex* Data, cuComplex* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneCol, cudaStream_t stream);
void cuSetCol(float* Data, float* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneCol, cudaStream_t stream);
void cuSetCol(int* Data, int* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneCol, cudaStream_t stream);

/*设置某阵，不开辟新的空间*/
void cuSetBand(cuComplex* Data, cuComplex* Value, int OneBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);
void cuSetBand(float* Data, float* Value, int OneBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);
void cuSetBand(int* Data, int* Value, int OneBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream);

/*将data赋值为Value*/
void Valuate(cuComplex* data, cuComplex Value, unsigned int Size, cudaStream_t stream);
void Valuate(float* data, float Value, unsigned int Size, cudaStream_t stream);
void Valuate(int* data, int Value, unsigned int Size, cudaStream_t stream);
void Valuate(bool* data, bool Value, unsigned int Size, cudaStream_t stream);

#endif // _CUDAARRAYTOOLS_CUH_
