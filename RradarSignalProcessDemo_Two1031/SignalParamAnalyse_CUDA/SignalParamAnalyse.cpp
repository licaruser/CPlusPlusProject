#include "SignalParamAnalyse.h"
#include "ParamTest.cuh"
#include "CommonKernel.cuh"
#include "ParamAnalyse.h"

inline CUDASignalParamAnalyse::CUDASignalParamAnalyse()
{

}

CUDASignalParamAnalyse::~CUDASignalParamAnalyse()
{

}



void CUDASignalParamAnalyse::StepAdvance(const vector<vector<complex<double>>>& trf, const vector<complex<double>>& SourceData)
{
	//执行截取及模值计算
	ParamAnalyse* cpuParamAnalyse = new ParamAnalyse;
	vector<vector<double>> SourceDouData = cpuParamAnalyse->ReturnNeedData(trf);
	
	//使用GPU计算二维vector的平均值
	ParamGpuProgressHe(SourceDouData, SourceData);


	// 先试验一段数据，就全当试错了！  试验正确
	vector<complex<double>> temp = trf.at(0);
	//Test(temp);

	// GPU数据处理核心


	//gpu核函数<<<Block,Thread>>>
	//核函数第一个参数表示设备在执行核函数时使用的“并行线程块的数量”，也就是创建核函数的数量，并以并行的方式来运行他；
	//blockIdx是一个内置变量，变量中的值就是当前执行设备代码的线程块的索引；
	//Cuda支持二维线程块数组，并行线程块集合也称为一个线程格（Grid）；
	//启动线程块数组时，数组每一维的最大数量都不能超过65535，这是一种硬件限制。

	//核函数第二个参数表示每一个Block的线程数（Thread）
	//最大线程数量不能超过设备属性结构中maxThreadsPerBlock域的值，每个线程块大概是512；
	//BlockDim 是一个常量，是每个线程块中的线程数量
	//gridDim是线程块的大小
	//（N+127/128，128）
	//核函数不会对越过数组边界的内存进行读取或者写入；






	//CCMat Data;   //GPU上的数据格式

	////// 数据以矩阵的形式保存
	//int cols = 0;
	//int rows = 0;
	//cols = trf.size();//62500
	//if (trf.at(0).size() > 0)
	//{
	//	rows = trf.at(0).size(); //128
	//}
	//int rows = 2;
	//int cols = 1;


	//Data.Resize(rows, cols, 1);


	//cuComplex* Matrix[2][1];
	//Data = float(trf.at(0).at(0).real());


	//for (int ii = 0; ii < cols; ii++)  //62500
	//{
	//	for (int jj = 0; jj < rows; jj++) //256
	//	{
	//		//Matrix[jj][ii]->x = trf.at(ii).at(jj).real();
	//		//Matrix[jj][ii]->y = trf.at(ii).at(jj).imag();
	//		Data.setCol()
	//	}
	//}

	////将CPU代码转换到GPU上
	//cuComplex *gpuMatrix;
	//cudaMalloc((void**)& gpuMatrix, rows * cols * sizeof(cuComplex));

	////cudaMemcpy(gpuMatrix, Matrix, rows * cols * sizeof(cuComplex), cudaMemcpyHostToDevice);

	//tools::SaveDeviceComplexMatrix(gpuMatrix, rows, cols, "save\\gpuMatrix.txt");
}
