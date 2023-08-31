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
	//ִ�н�ȡ��ģֵ����
	ParamAnalyse* cpuParamAnalyse = new ParamAnalyse;
	vector<vector<double>> SourceDouData = cpuParamAnalyse->ReturnNeedData(trf);
	
	//ʹ��GPU�����άvector��ƽ��ֵ
	ParamGpuProgressHe(SourceDouData, SourceData);


	// ������һ�����ݣ���ȫ���Դ��ˣ�  ������ȷ
	vector<complex<double>> temp = trf.at(0);
	//Test(temp);

	// GPU���ݴ������


	//gpu�˺���<<<Block,Thread>>>
	//�˺�����һ��������ʾ�豸��ִ�к˺���ʱʹ�õġ������߳̿����������Ҳ���Ǵ����˺��������������Բ��еķ�ʽ����������
	//blockIdx��һ�����ñ����������е�ֵ���ǵ�ǰִ���豸������߳̿��������
	//Cuda֧�ֶ�ά�߳̿����飬�����߳̿鼯��Ҳ��Ϊһ���̸߳�Grid����
	//�����߳̿�����ʱ������ÿһά��������������ܳ���65535������һ��Ӳ�����ơ�

	//�˺����ڶ���������ʾÿһ��Block���߳�����Thread��
	//����߳��������ܳ����豸���Խṹ��maxThreadsPerBlock���ֵ��ÿ���߳̿�����512��
	//BlockDim ��һ����������ÿ���߳̿��е��߳�����
	//gridDim���߳̿�Ĵ�С
	//��N+127/128��128��
	//�˺��������Խ������߽���ڴ���ж�ȡ����д�룻






	//CCMat Data;   //GPU�ϵ����ݸ�ʽ

	////// �����Ծ������ʽ����
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

	////��CPU����ת����GPU��
	//cuComplex *gpuMatrix;
	//cudaMalloc((void**)& gpuMatrix, rows * cols * sizeof(cuComplex));

	////cudaMemcpy(gpuMatrix, Matrix, rows * cols * sizeof(cuComplex), cudaMemcpyHostToDevice);

	//tools::SaveDeviceComplexMatrix(gpuMatrix, rows, cols, "save\\gpuMatrix.txt");
}
