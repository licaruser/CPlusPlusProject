#include "ParamTest.cuh"
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






	//printf("sizes:%d\n", sizes);
	//printf("Average: %f\n", h_average);

	cudaDeviceSynchronize(); //在此之后，可以确保核函数已经执行完成
	cudaFree(d_vector);
	cudaFree(d_average);


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