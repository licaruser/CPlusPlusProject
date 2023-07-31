#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <vector_types.h>

__global__ void hello_from_gpu()
{
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	printf("%d,%d.\n", bid, tid);
}

int main()
{
	const dim3 gridSize(2);
	const dim3 blockSize(3);
	printf("start\n");
	hello_from_gpu << <gridSize, blockSize >> > ();
	printf("endn\n");
	cudaDeviceSynchronize();
	return 0;
}