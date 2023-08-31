#include "accum.cuh"

__global__ void ReducePartSumKernel(float* input, float* part_sum,
	unsigned int SumBlockNum, unsigned int Len,
	unsigned int part_num)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < SumBlockNum)
	{
		for (int i = 0; i < part_num; i++)
		{
			part_sum[idx] = part_sum[idx] + input[idx * part_num];
		}
	}
}

__global__ void ScanWithBaseSum(float* data, float* part_sum,float* output,
	unsigned int SumBlockNum, unsigned int Len,
	unsigned int part_num) 
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Len)
	{
		unsigned SumBlockBias = idx / part_num;
		unsigned XBias = idx % part_num;
		
		output[idx] = 0.0;

		__syncthreads();

		// 各部分求和叠加
		for (int i = 0; i < SumBlockBias; i++)
		{
			output[idx] = output[idx] + part_sum[i];
		}

		// 区域内求和
		for (int i = 1; i < XBias; i++)
		{
			output[idx] = output[idx] + data[SumBlockBias * part_num + i];
		}

	}

}

// 前缀和
FCMat accum(FCMat& Data, cudaStream_t stream)
{
	auto XLen = Data.elements();

	int part_num = 1024;
	unsigned int SumBlockNum = XLen / part_num + 1;
	int threadsPerBlock = 1024;
	int blocksPerGrid = (SumBlockNum + threadsPerBlock - 1) / threadsPerBlock;

	// 计算各part的累加和
	FCMat SumTemp(SumBlockNum, 1, 1);
	ReducePartSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data.FirstAddr(), SumTemp.FirstAddr(),SumBlockNum, XLen, part_num);

	FCMat Result(Data.dims(0), Data.dims(1), Data.dims(2));
	blocksPerGrid = (XLen + threadsPerBlock - 1) / threadsPerBlock;
	ScanWithBaseSum << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data.FirstAddr(), SumTemp.FirstAddr(), Result.FirstAddr(),
		SumBlockNum, XLen, part_num);

	SumTemp.Free();

	return Result;
}

// 前缀和
void accum(FCMat& Data, FCMat& Result, cudaStream_t stream)
{
	Result.Resize(Data.dims(0), Data.dims(1), Data.dims(2));

	auto XLen = Data.elements();

	int part_num = 1024;
	unsigned int SumBlockNum = XLen / part_num + 1;
	int threadsPerBlock = 1024;
	int blocksPerGrid = (SumBlockNum + threadsPerBlock - 1) / threadsPerBlock;

	// 计算各part的累加和
	FCMat SumTemp(SumBlockNum, 1, 1);
	ReducePartSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data.FirstAddr(), SumTemp.FirstAddr(), SumBlockNum, XLen, part_num);

	
	blocksPerGrid = (XLen + threadsPerBlock - 1) / threadsPerBlock;
	ScanWithBaseSum << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data.FirstAddr(), SumTemp.FirstAddr(), Result.FirstAddr(),
		SumBlockNum, XLen, part_num);

	SumTemp.Free();

}

