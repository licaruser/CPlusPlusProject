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
	int cols = 0;
	int rows = 0;
	cols = trf.size();//62500
	if (trf.at(0).size() > 0)
	{
		rows = trf.at(0).size(); //128
	}
	
	Data.Resize(rows, cols, 1);


	cuComplex* Matrix[128][62500];
	for (int ii = 0; ii < cols; ii++)  //62500
	{
		for (int jj = 0; jj < rows; jj++) //128
		{
			Matrix[jj][ii]->x = trf.at(ii).at(jj).real();
			Matrix[jj][ii]->y = trf.at(ii).at(jj).imag();
		}
	}

	//将CPU代码转换到GPU上
	cuComplex *gpuMatrix;
	cudaMalloc((void**)& gpuMatrix, rows * cols * sizeof(cuComplex));

	cudaMemcpy(gpuMatrix, Matrix, rows * cols * sizeof(cuComplex), cudaMemcpyHostToDevice);


}
