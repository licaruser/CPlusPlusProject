#include "cudaRand.cuh"

__global__ void setup_kernel(curandState *state, unsigned long seed, unsigned int len)
{
	// ��ȡ�̺߳ţ�һά�ṹ
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < len)
	{
		// ��ʼ�������������
		curand_init(seed, tid, 0, &state[tid]);
	}

}

__global__ void use(curandState *globalState, cuComplex* data, unsigned int len, unsigned int elements)
{
	// ��ȡ�̺߳ţ�һά�ṹ
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < len)
	{
		for (int ii = 0; ii < 128; ii++)
		{
			int idx = ii * len + tid;
			if (idx < elements)
			{
				data[idx] = curand_normal2(globalState + tid);
			}
		}
	}
}

__global__ void use(curandState *globalState, float* data, unsigned int len, unsigned int elements)
{
	// ��ȡ�̺߳ţ�һά�ṹ
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < len)
	{
		for (int ii = 0; ii < 128; ii++)
		{
			int idx = ii * len + tid;
			if (idx < elements)
			{
				data[idx] = curand_normal(globalState + tid);
			}
		}
	}
}

//����(bias, float+bias)�ľ��ȷֲ�
__global__ void cuuseUrand(curandState *globalState, float* data, unsigned int len, float scale, float bias, unsigned int elements)
{
	// ��ȡ�̺߳ţ�һά�ṹ
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < len)
	{
		for (int ii = 0; ii < 128; ii++)
		{
			int idx = ii * len + tid;
			if (idx < elements)
			{
				data[idx] = curand_uniform(globalState + tid) * scale + bias;
			}
		}
	}
}



//void cuGuassRand(float* GRand, unsigned int N)
//{
//	curandGenerator_t gen;
//	// ��������������������������ʽ����ǰΪĬ��
//	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
//	// �������������
//	curandSetPseudoRandomGeneratorSeed(gen, time(0));
//	// ���ɸ�˹����
//	curandGenerateNormal(gen, GRand, N, 0.0, 1.0);
//	// ����������
//	curandDestroyGenerator(gen);
//}

// ���ɸ�˹�ֲ����������
CCMat Randn(unsigned int Rows, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	CCMat data(Rows, 1, 1);
	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows);
	devStates.Free();

	return data;
}

// ���ɸ�˹�ֲ����������
void Randn(CCMat& data, unsigned int Rows, cudaStream_t stream)
{
	data.Resize(Rows, 1, 1);

	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));


	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows);
	devStates.Free();

}

// ���ɸ�˹�ֲ����������
CCMat Randn(unsigned int Rows, unsigned int Cols, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	CCMat data(Rows, Cols, 1);	
	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows*Cols);
	devStates.Free();

	return data;
}

// ���ɸ�˹�ֲ����������
void Randn(CCMat& data, unsigned int Rows, unsigned int Cols, cudaStream_t stream)
{
	data.Resize(Rows, Cols, 1);

	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	
	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows*Cols);

	devStates.Free();
}

// ���ɸ�˹�ֲ����������
CCMat Randn(unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols*Bands + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	CCMat data(Rows, Cols, Bands);
	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows*Cols*Bands);
	devStates.Free();

	return data;
}

// ���ɸ�˹�ֲ����������
void Randn(CCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream)
{
	data.Resize(Rows, Cols, Bands);

	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols*Bands + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	
	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows*Cols*Bands);

	devStates.Free();
}

// ����ʵ����˹�ֲ����������
void Randn(FCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols*Bands + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));
	
	data.Resize(Rows, Cols, Bands, 0.0, stream);
	use << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), Rows*Cols*Bands);

	devStates.Free();
}


//����(bias, float+bias)�ľ��ȷֲ�
FCMat Randu(unsigned int Rows, float scale, float bias, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	FCMat data(Rows, 1, 1);
	cuuseUrand << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), scale, bias, Rows);
	devStates.Free();

	return data;
}

//����(bias, float+bias)�ľ��ȷֲ�
void Randu(FCMat& data, unsigned int Rows, float scale, float bias, cudaStream_t stream)
{
	data.Resize(Rows, 1, 1);

	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	
	cuuseUrand << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), scale, bias, Rows);

	devStates.Free();
}

//����(bias, float+bias)�ľ��ȷֲ�
FCMat Randu(unsigned int Rows, unsigned int Cols, float scale, float bias, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));
	
	FCMat data(Rows, Cols, 1);
	cuuseUrand << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), scale, bias, Rows*Cols);

	devStates.Free();

	return data;
}

//����(bias, float+bias)�ľ��ȷֲ�
void Randu(FCMat& data, unsigned int Rows, unsigned int Cols, float scale, float bias, cudaStream_t stream)
{
	data.Resize(Rows, Cols, 1);
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));
	
	cuuseUrand << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), scale, bias, Rows*Cols);

	devStates.Free();
}

//����(bias, float+bias)�ľ��ȷֲ�
FCMat Randu(unsigned int Rows, unsigned int Cols, unsigned int Bands, float scale, float bias, cudaStream_t stream)
{
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols*Bands + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));

	FCMat data(Rows, Cols, Bands);
	cuuseUrand << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), scale, bias, Rows*Cols*Bands);
	devStates.Free();

	return data;
}


//����(bias, float+bias)�ľ��ȷֲ�
void Randu(FCMat& data, unsigned int Rows, unsigned int Cols, unsigned int Bands, float scale, float bias, cudaStream_t stream)
{
	data.Resize(Rows, Cols, Bands);
	srand(time(0) * 1e4);
	int threadsPerBlock1 = 128;
	CudaArray<curandState> devStates((Rows*Cols*Bands + threadsPerBlock1 - 1) / threadsPerBlock1, 1, 1);
	int blocksPerGrid = (devStates.dims(0) + threadsPerBlock1 - 1) / threadsPerBlock1;
	setup_kernel << <blocksPerGrid, threadsPerBlock1, 0, stream >> >(devStates.FirstAddr(), unsigned long(rand()), devStates.dims(0));
	
	cuuseUrand << <blocksPerGrid, threadsPerBlock1, 0, stream >> > (devStates.FirstAddr(), data.FirstAddr(), devStates.dims(0), scale, bias, Rows*Cols*Bands);

	devStates.Free();
}



