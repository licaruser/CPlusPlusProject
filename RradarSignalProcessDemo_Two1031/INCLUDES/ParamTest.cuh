#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "sm_20_atomic_functions.h"
#include <cstdio>

#include <vector_types.h>
#include <vector>
//#include <complex>
#include <cuda.h>
#include <cuComplex.h>
#include <CudaArray.cuh>
#include <tools.cuh>
#include <cuda_runtime.h>


__global__ void hello_from_gpu();



// GPU数据核心处理
void ParamGpuProgressHe(vector<vector<double>>& lep, const vector<complex<double>>& Source_Data);

void Test(vector<complex<double>>& temp);

__global__ void vectorAverage(float* d_vector, int width, int height, float *d_average);

__global__ void CalculateNN_Array(float* d_vector, int width, int heigth, float h_average, unsigned int* NN);
//__global__ void ComplexMat(cuComplex* Res, float* Real, float* Imag, int elements);

// 凭借脉冲上升沿下降沿进行脉内带宽识别
//__global__ void WithPulseCalculate_GPU(unsigned int* NN, int UpVec, int DownVec, int* Array);