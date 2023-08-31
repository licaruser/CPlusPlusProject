#include "CudaArrayTool.cuh"

/*��ȡ������ĳ��λ�õ�ֵ*/
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

/*���þ����е�һ��*/
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

/*���þ����е�һ��*/
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

/*���þ����е�һ����*/
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


/*��ֵ�˺���
������*a = elements
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

//*��ȡĳ����
//data Ϊ��������
//Index Ϊ��Ҫ��ȡ����������
//result Ϊ�������
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


/*����ĳ�У��������µĿռ�*/
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

/*����ĳ�У��������µĿռ�*/
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

/*����ĳ�󣬲������µĿռ�*/
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


/*��data��ֵΪValue*/
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


