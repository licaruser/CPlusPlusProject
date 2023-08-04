

#ifndef __ERROR_CUH__
#define __ERROR_CUH__

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <fstream>
#include <string>

#define HANDLEERROR(err) (HandleError(err, __FILE__, __LINE__))

static void HandleError(cudaError_t err, const char  *file, int line) {
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Error %d: \"%s\" in %s at line %d\n", int(err), cudaGetErrorString(err), file, line);
		exit(int(err));
	}
}



inline void CheckLastError()
{
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "failed-%d: %s\n", int(cudaStatus), cudaGetErrorString(cudaStatus));
	}
	else
	{
		printf("current cudaStatus is success.\n");
	}
}

inline void DisStreamErr(cudaStream_t stream)
{
	cudaError err = cudaStreamSynchronize(stream);
	printf("µ±Ç°Á÷×´Ì¬:%d\n", int(err));
}

#endif // !__ERROR_CUH__