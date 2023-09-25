#include "ParamTest.cuh"
#include <numeric>
#include "CommonKernel.cuh"




void ParamGpuProgressHe(vector<vector<double>>& lep,const vector<complex<double>>& Source_Data)
{
	//int rows = 0;
	//int cols = 0;
	//if (lep.size() > 0)
	//{
	//	cols = lep.size();
	//	if (lep.at(0).size() > 0)
	//	{
	//		rows = lep.at(0).size();
	//	}
	//}
	//else
	//{
	//	return;
	//}
	//CCMat Matrix_Test(256, 1, 1);
	//Matrix_Test[0].x = 
	//float Matrix_Real[][];
	//float Matrix_Imag[][];
	// begin,开始
	const int rows = 128;
	const int cols = 62500;
	const int sizes = rows * cols;

	static float Matrix_vector[sizes];
	static unsigned int NN_Matrix[sizes];
	int aa_cols = 0;
	int bb_rows = 0;
	// 循环赋值
	for (int ii = 0; ii < sizes; ii++)
	{
		//printf("输出一列的每行数据aa_cols:%d,bb_rows:%d\n", aa_cols, bb_rows);//
		Matrix_vector[ii] = lep.at(aa_cols).at(bb_rows);  //先放的是一列的每行数据，再放下一列
		NN_Matrix[ii] = 0;

		if (bb_rows == rows - 1)
		{
			bb_rows = -1;
			if (aa_cols == cols - 1)
			{
				break;
			}
			aa_cols++;
			//continue;
		}
		bb_rows++;
	}
	//printf("Matrix_vector[8000000]:%f\n", Matrix_vector[7999999]);
	////测试总和使用，测试结论正确
	////float total_he_vector(0.0);
	////for (int ii = 0; ii < rows * cols; ii++)
	////{
	////	total_he_vector = total_he_vector + Matrix_vector[ii];
	////}
	////printf("total_he_vector:%f\n", total_he_vector);

	float* d_vector;
	float* d_average;
	cudaMalloc((void**)& d_vector, sizes * sizeof(float));
	cudaMalloc((void**)& d_average, sizeof(float));
	cudaError_t cudaStatus_vector = cudaMemcpy(d_vector, Matrix_vector, sizes * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus_vector != cudaSuccess)
	{
		std::cout << "CudaMemcpy failed:" << cudaGetErrorString(cudaStatus_vector) << std::endl;
	}
	else
	{
		std::cout << "cudaMemcpy succeeded." << std::endl;
	}
	cudaError_t cudaStatus_daverage = cudaMemset(d_average, 0, sizeof(float));
	if (cudaStatus_daverage != cudaSuccess)
	{
		std::cout << "CudaMemcpy failed:" << cudaGetErrorString(cudaStatus_daverage) << std::endl;
	}
	else
	{
		std::cout << "cudaMemcpy succeeded." << std::endl;
	}

	//float* cudaCols;
	//float* cudaRows;
	//cudaMalloc((void**)& cudaCols, sizeof(float));
	//cudaMemset(cudaCols, cols, sizeof(float));
	//cudaMalloc((void**)& cudaRows, sizeof(float));
	//cudaMemset(cudaRows, rows, sizeof(float));
	//printf("d_vector:%f\n", d_vector);

	dim3 blockSize(128, 62500);//线程块    //128行，62500列
	int grid_x = (rows + blockSize.x - 1) / blockSize.x;//cols是列，rows是行
	int grid_y = (cols + blockSize.y - 1) / blockSize.y;
	//printf("grid_x:%d,grid_y:%d\n", grid_x, grid_y);
	dim3 gridSize(grid_x, grid_y);//1，1
	//printf("d_vector:%f\n", d_vector);

	cudaError_t cudaStatus;//设置cuda核函数状态变量
	//blockSize是并行线程块的数量，blockIdx是一个内置变量;gridSize is 每个线程块对应的线程数，每个线程块的线程被限制在512
	//并行的线程块集合又称为一个线程格(Grid)
	vectorAverage << <blockSize, gridSize >> > (d_vector, rows, cols, d_average);
	cudaStatus = cudaGetLastError(); //获取核函数执行后的状态变量状态
	if (cudaStatus != cudaSuccess)
	{
		//核函数执行失败
		printf("Kernel execution failed:%s\n", cudaGetErrorString(cudaStatus));
	}

	//所有数据的平均值
	float h_average;
	cudaError_t cudaStatus_haverage = cudaMemcpy(&h_average, d_average, sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus_haverage != cudaSuccess)
	{
		std::cout << "CudaMemcpy failed:" << cudaGetErrorString(cudaStatus_haverage) << '\n'/*std::endl*/;
	}
	else
	{
		std::cout << "cudaMemcpy succeeded." << '\n'/*std::endl*/;
	}
	//printf("cudaCols:%f\n", cudaCols);
	printf("h_average:%f\n", h_average);

	h_average /= sizes;
	//第一步利用核函数计算平均值已完成，第二步利用门限值判别超过门限值的所有位置信息
	unsigned int* NN;
	cudaError_t cudaMalloc_NN = cudaMalloc((void**)& NN, sizes * sizeof(unsigned int));
	if (cudaMalloc_NN != cudaSuccess)
	{
		std::cout << "CudaMalloc failed:" << cudaGetErrorString(cudaMalloc_NN) << std::endl;
	}
	//printf("sizeof(unsigned int): %d.\n", sizeof(unsigned int));
	cudaError_t cudaStatus_NN = cudaMemcpy(NN, NN_Matrix, sizes * sizeof(unsigned int), cudaMemcpyHostToDevice);
	if (cudaStatus_NN != cudaSuccess)
	{
		std::cout << "CudaMemcpy failed:" << cudaGetErrorString(cudaStatus_NN) << std::endl;
	}
	//编写一个核函数
	cudaError_t cudaStatus_NNfunction;//设置cuda核函数状态变量
	CalculateNN_Array << <blockSize, gridSize >> > (d_vector, rows, cols, h_average * 30, NN);
	cudaStatus_NNfunction = cudaGetLastError(); //获取核函数执行后的状态变量状态
	if (cudaStatus_NNfunction != cudaSuccess)
	{
		//核函数执行失败
		printf("Kernel execution failed:%s\n", cudaGetErrorString(cudaStatus_NNfunction));
	}

	//完成第二步骤将数据与门限判别之后，进行第三步，识别出每一列和行连续是1的索引，并将其输出
	//Source_Data
	//先识别一维原始数据的脉冲起始时间和结束时间，再根据此时间，计算二维带宽

	////根据时域数据计算脉冲上升沿下降沿时刻
	//vector<int> UpVec;           //上升沿vector
	//vector<int> DownVec;         //下降沿vector
	//vector<int> PulseWidthVec;   //对应脉宽
	//double CurrentTime = 0.0;
	//double m_JamInitSampleFs = 0.0;
	////m_pSignalProcessor->GetJamInitSampleFs(m_JamInitSampleFs);
	//// 参数估计--估计脉冲上升沿、下降沿、脉宽
	//if (Source_Data.size() > 0)
	//{
	//	ParameterAnalysis(Source_Data, UpVec, DownVec, PulseWidthVec, CurrentTime, m_JamInitSampleFs);//输出雷达信号参数单位(毫秒)[上升沿、下降沿、脉宽]
	//}



	//for (int ii = 0;ii<UpVec.size();ii++)
	//{
	//	int begin_data = UpVec.at(ii);  //起始位置
	//	int end_data = DownVec.at(ii);  //结束位置

	//	int ThreadNum = end_data - begin_data;

	//	int* Array;
	//	cudaError_t cudaMalloc_Array = cudaMalloc((void**)& Array, ThreadNum * 2 * sizeof(int));
	//	if (cudaMalloc_Array != cudaSuccess)
	//	{
	//		std::cout << "CudaMalloc failed:" << cudaGetErrorString(cudaMalloc_Array) << std::endl;
	//	}
	//	cudaError_t cudaStatus_Array = cudaMemcpy(Array, 0, ThreadNum * 2 * sizeof(int), cudaMemcpyHostToDevice);
	//	if (cudaStatus_Array != cudaSuccess)
	//	{
	//		std::cout << "CudaMemcpy failed:" << cudaGetErrorString(cudaStatus_Array) << std::endl;
	//	}


	//}




	//printf("sizes:%d\n", sizes);
	//printf("Average: %f\n", h_average);

	cudaDeviceSynchronize(); //在此之后，可以确保核函数已经执行完成
	cudaFree(d_vector);
	cudaFree(d_average);


}

void ParameterAnalysis(const vector<complex<double>>& AllRadarData, vector<int>& UpData, vector<int>& DownData, vector<int>& PulseWidth, double& time, const double Fs)
{
	// 估计AllRadarData的上升沿、下降沿和脉宽，并存入对应的vector中，便于后续进行分选处理
	int PulseUpFlag = 0;
	int PulseDownFlag = 0;
	int PulseUpPos = 0;
	int PulseDownPos = 0;

	// 1、计算abs值
	vector<double> HeDataABS;
	HeDataABS.resize(AllRadarData.size());
	for (int ii = 0; ii < AllRadarData.size(); ii++)
	{
		float RealData = AllRadarData.at(ii).real();
		float ImagData = AllRadarData.at(ii).imag();
		HeDataABS[ii] = sqrt(RealData * RealData + ImagData * ImagData);
	}

	// 2、利用前面纯噪声部分计算噪声功率和方差
	vector<double> DataRead;
	DataRead.resize(HeDataABS.size());
	for (int jj = 0; jj < HeDataABS.size(); jj++)
	{
		DataRead[jj] = HeDataABS[jj] * HeDataABS[jj];   //计算功率
	}
	HeDataABS.clear();

	// 检测所用参数
	int PulseFlag = 0;
	int	CntBegin = 0;
	int CntEnd = 0;
	double NoisePower = 0.0;
	double DataPower = 0.0;
	int BeginFrameThrh = 128;
	int EndFrameThrh = 128;
	int FrameSize = 128; //帧长度
	NoisePower = accumulate(DataRead.begin(), DataRead.begin() + FrameSize, 0.0) / FrameSize; //计算前128个点的噪声功率均值

	vector<double> data_buff;
	data_buff.resize(FrameSize);
	// 3、能量检测方法检测脉冲边沿
	vector <int> UpPoint;
	vector <int> DownPoint;
	vector <int> PulseWidthPoint;
	for (int k = 0; k < DataRead.size(); k++)
	{

		for (int j = 1; j < FrameSize; j++)
		{
			data_buff[FrameSize - j] = data_buff[FrameSize - j - 1];
		}
		data_buff[0] = DataRead[k];

		//if (k == 800)
		//{
		//	int oo = 0;  //测试
		//}
		DataPower = accumulate(data_buff.begin(), data_buff.end(), 0.0) / data_buff.size();  //计算该帧信号平均功率

		if (DataPower > 1.5 * NoisePower)
		{
			CntEnd = 0;
			CntBegin = CntBegin + 1;

			if (CntBegin >= BeginFrameThrh) // 连续有超过BeginFrameThrh帧信号超过检测门限，则认为是脉冲开始
			{
				if (PulseFlag == 0) // 此前还未检测到脉冲
				{
					// 寻找精确的脉冲起始位置
					PulseUpPos = k - FrameSize;
					UpPoint.push_back(PulseUpPos);
				}

				PulseUpFlag = 1;
				PulseFlag = 1;

			}
		}
		else
		{
			if (PulseFlag == 1)
			{
				CntEnd = CntEnd + 1;
				if (CntEnd >= EndFrameThrh)  //有连续超过EndFrameThrh帧信号低于检测门限，则认为是脉冲结束
				{
					// 寻找精确的脉冲结束位置
					PulseDownPos = k - FrameSize - 127;
					DownPoint.push_back(PulseDownPos);
					PulseWidthPoint.push_back(PulseDownPos - PulseUpPos);
					PulseDownFlag = 1;
					PulseFlag = 0;
					CntBegin = 0;
				}
			}
			else
			{
				CntBegin = 0;
				CntEnd = 0;
			}
		}
	}

	// 4、点数换算时间--将点数Pos位置换算到时间上(单位/毫秒)
	for (int aa = 0; aa < UpPoint.size(); aa++)
	{
		//double Tmp_UpTime;
		//Tmp_UpTime = ((time + UpPoint[aa] / Fs) * 1e3);
		UpData.push_back(UpPoint.at(aa));
	}
	for (int bb = 0; bb < DownPoint.size(); bb++)
	{
		//double Tmp_DownTime;
		//Tmp_DownTime = ((time + DownPoint[bb] / Fs) * 1e3);
		DownData.push_back(DownPoint.at(bb));
	}
	for (int cc = 0; cc < PulseWidthPoint.size(); cc++)
	{
		//double Tmp_PulseWidthTime;
		//Tmp_PulseWidthTime = ((time * 1e3 + PulseWidthPoint[cc] / Fs) * 1e3);
		PulseWidth.push_back(PulseWidthPoint.at(cc));
	}



	//// 4、点数换算时间--将点数Pos位置换算到时间上(单位/毫秒)
	//for (int aa = 0; aa < UpPoint.size(); aa++)
	//{
	//	double Tmp_UpTime;
	//	Tmp_UpTime = ((time + UpPoint[aa] / Fs) * 1e3);
	//	UpData.push_back(Tmp_UpTime);
	//}
	//for (int bb = 0; bb < DownPoint.size(); bb++)
	//{
	//	double Tmp_DownTime;
	//	Tmp_DownTime = ((time + DownPoint[bb] / Fs) * 1e3);
	//	DownData.push_back(Tmp_DownTime);
	//}
	//for (int cc = 0; cc < PulseWidthPoint.size(); cc++)
	//{
	//	double Tmp_PulseWidthTime;
	//	Tmp_PulseWidthTime = ((time * 1e3 + PulseWidthPoint[cc] / Fs) * 1e3);
	//	PulseWidth.push_back(Tmp_PulseWidthTime);
	//}
}


void Test(vector<complex<double>>& temp)
{

	float Matrix_Real[256];
	float Matrix_Imag[256];

	for (int ii = 0; ii < temp.size(); ii++)
	{
		Matrix_Real[ii] = temp.at(ii).real();
		Matrix_Imag[ii] = temp.at(ii).imag();
	}

	const int Mem = 256 * sizeof(float);

	//Gpu申请内存
	float* matrix_cuda_real;
	float* matrix_cuda_imag;

	cudaMalloc((void**)& matrix_cuda_real, Mem);
	cudaMalloc((void**)& matrix_cuda_imag, Mem);

	//将cpu数据拷至gpu上
	cudaMemcpy(matrix_cuda_real, Matrix_Real, Mem, cudaMemcpyHostToDevice);
	cudaMemcpy(matrix_cuda_imag, Matrix_Imag, Mem, cudaMemcpyHostToDevice);

	int threadsPerBlock1 = 256;  //每个线程块的线程数为256；
	int blocksPerGrid1 = (256 + threadsPerBlock1 - 1) / threadsPerBlock1;

	//printf("%d,%d\n", blocksPerGrid1, threadsPerBlock1);
	CCMat Matrix_Test(256, 1, 1);
	//cudaStream_t Stream;
	//ComplexMat << <blocksPerGrid1, threadsPerBlock1 >> > (Matrix_Test.FirstAddr(), matrix_cuda_real, matrix_cuda_imag, 256);
	tools::SaveDeviceComplexMatrix(Matrix_Test.FirstAddr(), Matrix_Test.Row(), Matrix_Test.Col(), "save\\gpuMatrix.txt");

	//const dim3 gridSize(2);
	//const dim3 blockSize(3);
	//printf("start\n");
	//hello_from_gpu << <gridSize, blockSize >> > ();
	//printf("endn\n");
	//cudaDeviceSynchronize();

}
//GPU计算的应用前景再很大程度上取决于能否从许多问题中发掘出大规模并行性;
__global__ void vectorAverage(float* d_vector, int width, int height, float* d_average)
{                                              //width is 128; height is 62500;
	                                         
	int idx = blockIdx.x * blockDim.x + threadIdx.x;   //blockDim is 线程块中每一维的线程数量
	//线程块 128 * 62500二维;blockIdx.x是从0~127，blockDim.x是每一个块内的线程数 is 1，threadIdx.x是0;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	//线程块 128 * 62500二维;blockIdx.y是从0~62499，blockDim.y is 1，threadIdx.y是0;
	//printf("%d,%d,%d,%d.\n", blockDim.x, blockDim.y, threadIdx.x, threadIdx.y);
	//printf("%d,%d.\n", blockIdx.x, blockIdx.y);
	if (idx < width && idy < height)
	{
		int index = idy * width + idx;
		//printf("%d.\n", index);
		atomicAdd(d_average, d_vector[index]);
	}
}

__global__ void CalculateNN_Array(float* d_vector, int width, int heigth, float h_average, unsigned int* NN)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	if (idx < width && idy < heigth)
	{
		int index = idy * width + idx;
		if (d_vector[index] > h_average)
		{
			NN[index] = 1;
		}
	}
}

//__global__ void ComplexMat(cuComplex* Res, float* Real, float* Imag, int elements)
//{
//	// 线程块中线程量 * 线程块x + 线程x
//	int i = blockDim.x * blockIdx.x + threadIdx.x;
//
//	if (i < elements)
//	{
//		Res[i].x = Real[i];
//		Res[i].x = Imag[i];
//	}
//	//return __global__ void();
//}

//int main()
//{
//	const dim3 gridSize(2);
//	const dim3 blockSize(3);
//	printf("start\n");
//	hello_from_gpu << <gridSize, blockSize >> > ();
//	printf("endn\n");
//	cudaDeviceSynchronize();
//	return 0;
//}