#include "CudaArrayTool.cuh"

/*获取矩阵中某个位置的值*/
template <typename Type>
__global__ void getRowKernel(Type *a, Type *result, int index, int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < elements)
	{
		if (index == idx)
			result = a;
	}
}

/*设置矩阵中的一行*/
__global__ void setRowKernel(cuComplex *a, cuComplex *Value, int OneRow, int CurrentBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneRow == R_bias)
			a[CurrentBand * (Row * Col) + C_bias * Row + R_bias] = Value[C_bias];
	}
}
__global__ void setRowKernel(float *a, float *Value, int OneRow, int CurrentBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneRow == R_bias)
			a[CurrentBand * (Row * Col) + C_bias * Row + R_bias] = Value[C_bias];
	}
}
__global__ void setRowKernel(int *a, int *Value, int OneRow, int CurrentBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneRow == R_bias)
			a[CurrentBand * (Row * Col) + C_bias * Row + R_bias] = Value[C_bias];
	}
}

/*设置矩阵中的一列*/
__global__ void setColKernel(cuComplex *a, cuComplex *Value, int OneCol, int CurrentBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneCol == C_bias)
			a[CurrentBand * (Row * Col) + C_bias * Row + R_bias] = Value[R_bias];
	}
}
__global__ void setColKernel(float *a, float *Value, int OneCol, int CurrentBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneCol == C_bias)
			a[CurrentBand * (Row * Col) + C_bias * Row + R_bias] = Value[R_bias];
	}
}
__global__ void setColKernel(int *a, int *Value, int OneCol, int CurrentBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneCol == C_bias)
			a[CurrentBand * (Row * Col) + C_bias * Row + R_bias] = Value[R_bias];
	}
}

/*设置矩阵中的一个阵*/
__global__ void setBandKernel(cuComplex *a, cuComplex *Value, int OneBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneBand == B_bias)
			a[B_bias * (Row * Col) + C_bias * Row + R_bias] = Value[C_bias * Row + R_bias];
	}
}
__global__ void setBandKernel(float *a, float *Value, int OneBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneBand == B_bias)
			a[B_bias * (Row * Col) + C_bias * Row + R_bias] = Value[C_bias * Row + R_bias];
	}
}
__global__ void setBandKernel(int *a, int *Value, int OneBand, int Row, int Col, int Band)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < Row * Col * Band)
	{
		unsigned int B_bias = idx / (Row * Col);
		unsigned int C_bias = (idx % (Row * Col)) / Row;
		unsigned int R_bias = (idx % (Row * Col)) % Row;
		if (OneBand == B_bias)
			a[B_bias * (Row * Col) + C_bias * Row + R_bias] = Value[C_bias * Row + R_bias];
	}
}


/*赋值核函数
描述：*a = elements
*/
__global__ void ValuateKernel(cuComplex *a, cuComplex Value, unsigned int Elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Elements)
	{
		a[idx].x = Value.x;
		a[idx].y = Value.y;
	}
}
__global__ void ValuateKernel(float *a, float Value, unsigned int Elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Elements)
	{
		a[idx] = Value;
	}
}
__global__ void ValuateKernel(int *a, int Value, unsigned int Elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Elements)
	{
		a[idx] = Value;
	}
}
__global__ void ValuateKernel(bool *a, bool Value, unsigned int Elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Elements)
	{
		a[idx] = Value;
	}
}

//*获取某个阵
//data 为输入数据
//Index 为需要获取的阵列索引
//result 为输出数据
//*/
//template <typename Type>
//void getBands(Type* data, int* Index, Type* result, unsigned int Rows, unsigned int Cols, unsigned int Bands)
//{
//	result = nullptr;
//
//	int threadsPerBlock = 1024;
//	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
//	getBandsKernel << <blocksPerGrid, threadsPerBlock >> >(data, result, OneCol * Rows, Rows * Cols * Bands);
//}


/*设置某行，不开辟新的空间*/
void cuSetRow(cuComplex* Data, cuComplex* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneRow, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setRowKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneRow, CurrentBand, Rows, Cols, Bands);
}
void cuSetRow(float* Data, float* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneRow, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setRowKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneRow, CurrentBand, Rows, Cols, Bands);
}
void cuSetRow(int* Data, int* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneRow, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setRowKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneRow, CurrentBand, Rows, Cols, Bands);
}

/*设置某列，不开辟新的空间*/
void cuSetCol(cuComplex* Data, cuComplex* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneCol, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setColKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneCol, CurrentBand, Rows, Cols, Bands);
}
void cuSetCol(float* Data, float* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneCol, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setColKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneCol, CurrentBand, Rows, Cols, Bands);
}
void cuSetCol(int* Data, int* Value, int CurrentBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, unsigned int OneCol, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setColKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneCol, CurrentBand, Rows, Cols, Bands);
}

/*设置某阵，不开辟新的空间*/
void cuSetBand(cuComplex* Data, cuComplex* Value, int OneBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setBandKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneBand, Rows, Cols, Bands);

}
void cuSetBand(float* Data, float* Value, int OneBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setBandKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneBand, Rows, Cols, Bands);

}
void cuSetBand(int* Data, int* Value, int OneBand, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
	setBandKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data, Value, OneBand, Rows, Cols, Bands);

}


/*将data赋值为Value*/
void Valuate(cuComplex* data, cuComplex Value, unsigned int Size, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
	ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data, Value, Size);
}
void Valuate(float* data, float Value, unsigned int Size, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
	ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data, Value, Size);
}
void Valuate(int* data, int Value, unsigned int Size, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
	ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data, Value, Size);
}
void Valuate(bool* data, bool Value, unsigned int Size, cudaStream_t stream)
{
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
	ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data, Value, Size);
}


//template <typename Type>
//void constants(CudaArray<Type> &data, unsigned int N)
//{
//	data.Resize(N, 1);
//	int threadsPerBlock = 256;
//	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
//	OneKernel << <blocksPerGrid, threadsPerBlock >> >(data.FirstAddr(), 1, N);
//
//
//}


