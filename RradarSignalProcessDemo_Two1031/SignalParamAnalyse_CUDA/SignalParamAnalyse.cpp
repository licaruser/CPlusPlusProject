#include "SignalParamAnalyse.h"

inline CUDASignalParamAnalyse::CUDASignalParamAnalyse()
{

}

CUDASignalParamAnalyse::~CUDASignalParamAnalyse()
{

}

void CUDASignalParamAnalyse::StepAdvance(const vector<vector<complex<double>>>& trf, const vector<complex<double>>& SourceData)
{
	CCMat Data;

	// 数据以矩阵的形式保存
	int cols;
	int rows;
	cols = trf.size();
	if (trf.at(0).size() > 0)
	{
		rows = trf.at(0).size();
	}
	
	complex<double>* Matrix[128][62500];


	//将CPU代码转换到GPU上
	cuComplex *gpuMatrix;
	cudaMalloc((void**)& gpuMatrix, int(trf.size()) * int(trf.at(0).size()) * sizeof(cuComplex));

	cudaMemcpy(d_matrix, matrix, rows * cols * sizeof(Complex), cudaMemcpyHostToDevice);
}
