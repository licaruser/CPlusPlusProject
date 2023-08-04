#pragma once


#ifndef SIGNALPARAMANALYSE
#define SIGNALPARAMANALYSE

#ifdef SIGNALPARAMANALYSE_EXPORTS
#define SIGNALPARAMANALYSE_API __declspec(dllexport)
#else
#define SIGNALPARAMANALYSE_API __declspec(dllexport)
#endif // SIGNALPARAMANALYSE_EXPORTS

#include <iostream>
#include <fstream>
#include <iomanip>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <vector_types.h>
#include <vector>
#include <complex>
#include <cuda.h>
#include <cuComplex.h>
#include <CudaArray.cuh>


using namespace std;



class SIGNALPARAMANALYSE_API CUDASignalParamAnalyse
{
public:
	CUDASignalParamAnalyse();
	~CUDASignalParamAnalyse();
	

	void StepAdvance(const vector<vector<complex<double>>>& trf, const vector<complex<double>>& SourceData);

protected:



private:



};

#endif // !SIGNALPARAMANALYSE
