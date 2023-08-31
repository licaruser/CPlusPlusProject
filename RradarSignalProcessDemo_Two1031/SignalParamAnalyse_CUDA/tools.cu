#include "tools.cuh"
#include "CommonKernel.cuh"
#include "cublas_v2.h"

namespace tools
{
	/*保存主机中的复数矩阵*/
	void SaveHostComplexMatrix(cuComplex *h_data, int rows, int cols, string path)
	{
		ofstream fp;
		fp.open(path, ios::out | ios::trunc);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				float RData = h_data[int(j * rows + i)].x;
				float IData = h_data[int(j * rows + i)].y;
				// 写入实部
				fp << RData;
				// 写入符号和虚部
				if (IData != -0 || IData != 0)
				{
					fp << (IData > 0 ? "+" : "") << IData << "i";
				}
				// 写入列上的间隔
				if (j != cols - 1)
				{
					fp << " ";
				}
			}
			fp << endl;
		}
		fp << endl;
		fp.close();
		std::cout << path << "	Matrix Save sucess!" << endl;

	}

	/*保存设备中的复数矩阵*/
	void SaveDeviceComplexMatrix(cuComplex *d_data, int rows, int cols, string path)
	{
		cudaDeviceSynchronize();
		cuComplex *h_data = (cuComplex *)malloc(sizeof(cuComplex) * rows * cols);
		cudaMemcpy(h_data, d_data, sizeof(cuComplex) * rows * cols, cudaMemcpyDeviceToHost);

		ofstream fp;
		fp.open(path, ios::out | ios::trunc);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				float RData = h_data[int(j * rows + i)].x;
				float IData = h_data[int(j * rows + i)].y;
				// 写入实部
				fp << RData;
				// 写入符号和虚部
				if (IData != -0 || IData != 0)
				{
					fp << (IData > 0 ? "+" : "") << IData << "i";
				}
				// 写入列上的间隔
				if (j != cols - 1)
				{
					fp << " ";
				}
			}
			fp << endl;
		}
		fp << endl;
		fp.close();
		std::cout << path << "	Matrix Save sucess!" << endl;
	}

	/*保存设备中的实数矩阵*/
	void SaveDeviceFloatMatrix(float* d_data, int rows, int cols, string path)
	{
		cudaDeviceSynchronize();
		float* h_data = (float*)malloc(sizeof(float) * rows * cols);
		cudaMemcpy(h_data, d_data, sizeof(float) * rows * cols, cudaMemcpyDeviceToHost);

		ofstream fp;
		fp.open(path, ios::out | ios::trunc);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				float Data = h_data[int(j * rows + i)];

				fp << Data;

				// 写入列上的间隔
				if (j != cols - 1)
				{
					fp << " ";
				}
			}
			fp << endl;
		}
		fp << endl;
		fp.close();
		std::cout << path << "	Matrix Save sucess!" << endl;
	}

	/*保存设备中的整数矩阵*/
	void SaveDeviceIntMatrix(int *d_data, int rows, int cols, string path)
	{
		cudaDeviceSynchronize();
		int *h_data = (int *)malloc(sizeof(int) * rows * cols);
		cudaMemcpy(h_data, d_data, sizeof(int) * rows * cols, cudaMemcpyDeviceToHost);

		ofstream fp;
		fp.open(path, ios::out | ios::trunc);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				// 写入实部
				fp << h_data[int(j * rows + i)];

				// 写入列上的间隔
				if (j != cols - 1)
				{
					fp << " ";
				}
			}
			fp << endl;
		}
		fp << endl;
		fp.close();
		std::cout << path << "	Matrix Save sucess!" << endl;
	}


	/*保存主机中的复数向量*/
	void SaveHostComplexVector(cuComplex *h_data, int elements, std::string path)
	{
		
		std::fstream fp;
		fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
		for (int i = 0; i < elements; i++)
		{
			float RData = h_data[i].x;
			float IData = h_data[i].y;
			// 写入实部
			fp << RData;
			// 写入符号和虚部
			if (IData != -0 || IData != 0)
			{
				fp << (IData > 0 ? "+" : "") << IData << "i";
			}
			fp << endl;
		}
		fp.close();
		std::cout << path << " Save sucess!" << endl;
	}

	/*保存设备中的复数向量*/
	void SaveDeviceComplexVector(cuComplex *d_data, int elements, std::string path)
	{
		cudaDeviceSynchronize();
		cuComplex *h_data = (cuComplex *)malloc(sizeof(cuComplex) * elements);
		cudaMemcpy(h_data, d_data, sizeof(cuComplex) * elements, cudaMemcpyDeviceToHost);

		std::fstream fp;
		fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
		for (int i = 0; i < elements; i++)
		{
			float RData = h_data[i].x;
			float IData = h_data[i].y;
			// 写入实部
			fp << RData;
			// 写入符号和虚部
			if (IData != -0 || IData != 0)
			{
				fp << (IData > 0 ? "+" : "") << IData << "i";
			}
			fp << endl;
		}
		fp.close();
		std::cout << path << " Save sucess!" << endl;
	}
	/*保存设备中的复数向量*/
	void SaveDeviceIntVector(int* d_data, int elements, std::string path)
	{
		cudaDeviceSynchronize();
		int* h_data = (int*)malloc(sizeof(int) * elements);
		cudaMemcpy(h_data, d_data, sizeof(int) * elements, cudaMemcpyDeviceToHost);

		std::fstream fp;
		fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
		for (int i = 0; i < elements; i++)
		{
			int Data = h_data[i];
			fp << Data;
			fp << endl;
		}
		fp.close();
		std::cout << path << " Save sucess!" << endl;
	}

	/*保存设备中的向量*/
	void SaveDeviceFloatVector(float *d_data, int elements, std::string path)
	{
		cudaDeviceSynchronize();
		float *h_data = (float *)malloc(sizeof(float) * elements);
		cudaMemcpy(h_data, d_data, sizeof(float) * elements, cudaMemcpyDeviceToHost);

		std::fstream fp;
		fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
		for (int i = 0; i < elements; i++)
		{
			float RData = h_data[i];
			// 写入实部
			fp << RData;
			fp << endl;
		}
		fp.close();
		std::cout << path << " Save sucess!" << endl;
	}


	/*保存主机中的short型向量*/
	void SaveHostShortVector(short* h_data, int elements, std::string path)
	{
		cudaDeviceSynchronize();

		std::fstream fp;
		fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
		for (int i = 0; i < elements; i++)
		{
			float RData = h_data[i];
			// 写入实部
			fp << RData;
			fp << endl;
		}
		fp.close();
		std::cout << path << " Save sucess!" << endl;
	}

	/*保存主机中的复数矩阵*/
	template <typename Type>
	void SaveHostVector(Type *d_data, int elements, string path)
	{
		std::fstream fp;
		fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
		for (int i = 0; i < elements; i++)
		{
			// 写入实部
			fp << d_data[i];
			fp << endl;
		}
		fp.close();
		std::cout << path << " Save sucess!" << endl;
	}

	void Zeros(CudaArray<cuComplex> &data, unsigned int N, cudaStream_t stream)
	{
		data.Resize(N, 1, 1);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
		cuComplex Temp;Temp.x = 0.0;Temp.y = 0.0;
		ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data.FirstAddr(), Temp, N);
	}
	void Zeros(CudaArray<float> &data, unsigned int N, cudaStream_t stream)
	{
		data.Resize(N, 1, 1);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
		ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data.FirstAddr(), 0.0, N);
	}

	/*dB Result = 20*log10(Data)*/
	FCMat dB(FCMat& Data, cudaStream_t stream)
	{
		unsigned int Size = Data.elements();
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;

		FCMat Result(Data.dims(0), Data.dims(1), Data.dims(2));
		dBKernel << <blocksPerGrid, threadsPerBlock >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

		return Result;
	}

	/*idB Result = pow(10, Data / 20)*/
	FCMat idB(FCMat& Data, cudaStream_t stream)
	{
		unsigned int Size = Data.elements();
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;

		FCMat Result(Data.dims(0), Data.dims(1), Data.dims(2));
		idBKernel << <blocksPerGrid, threadsPerBlock >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

		return Result;
	}

	/*将实数矩阵转为虚部为0的复数矩阵*/
	CCMat Complex(FCMat& Data, cudaStream_t stream)
	{
		unsigned int Size = Data.elements();
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;

		CCMat Result(Data.dims(0), Data.dims(1), Data.dims(2));
		ComplexKernel << <blocksPerGrid, threadsPerBlock >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

		return Result;
	}
	/*将实数矩阵转为虚部为0的复数矩阵*/
	void Complex(FCMat& Data, CCMat& Result, cudaStream_t stream)
	{
		unsigned int Size = Data.elements();
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;

		Result.Resize(Data.dims(0), Data.dims(1), Data.dims(2));
		ComplexKernel << <blocksPerGrid, threadsPerBlock >> > (Data.FirstAddr(), Result.FirstAddr(), Size);


	}


	/*将两个实数矩阵转为复数矩阵*/
	CCMat Complex(FCMat& Real, FCMat& Imag, cudaStream_t stream)
	{
		unsigned int RealSize = Real.elements();
		unsigned int ImagSize = Imag.elements();
		if (RealSize != ImagSize)
		{
			printf("维度不一致\n");throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (RealSize + threadsPerBlock - 1) / threadsPerBlock;

		CCMat Result(Real.dims(0), Real.dims(1), Real.dims(2));
		ComplexMat << <blocksPerGrid, threadsPerBlock >> > (Result.FirstAddr(), Real.FirstAddr(), Imag.FirstAddr(), RealSize);

		return Result;
	}

	/*将两个实数矩阵转为复数矩阵*/
	void Complex(FCMat& Real, FCMat& Imag, CCMat& Result,  cudaStream_t stream)
	{
		unsigned int RealSize = Real.elements();
		unsigned int ImagSize = Imag.elements();
		if (RealSize != ImagSize)
		{
			printf("维度不一致\n"); throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (RealSize + threadsPerBlock - 1) / threadsPerBlock;

		Result.Resize(Real.dims(0), Real.dims(1), Real.dims(2));
		ComplexMat << <blocksPerGrid, threadsPerBlock >> > (Result.FirstAddr(), Real.FirstAddr(), Imag.FirstAddr(), RealSize);


	}

	FCMat Real(CCMat& Data, cudaStream_t Stream)
	{
		unsigned int Size = Data.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
		FCMat Result(Data.dims(0), Data.dims(1), Data.dims(2));
		Realkernel << <blocksPerGrid, threadsPerBlock, 0, Stream >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

		return Result;
	}

	void Real(CCMat& Data, FCMat& Result, cudaStream_t Stream)
	{
		unsigned int Size = Data.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
		Result.Resize(Data.dims(0), Data.dims(1), Data.dims(2));
		Realkernel << <blocksPerGrid, threadsPerBlock, 0, Stream >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

	}

	FCMat Imag(CCMat& Data, cudaStream_t Stream)
	{
		unsigned int Size = Data.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
		FCMat Result(Data.dims(0), Data.dims(1), Data.dims(2));
		Imagkernel << <blocksPerGrid, threadsPerBlock, 0, Stream >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

		return Result;
	}

	void Imag(CCMat& Data, FCMat& Result, cudaStream_t Stream)
	{
		unsigned int Size = Data.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Size + threadsPerBlock - 1) / threadsPerBlock;
		Result.Resize(Data.dims(0), Data.dims(1), Data.dims(2));
		Imagkernel << <blocksPerGrid, threadsPerBlock, 0, Stream >> > (Data.FirstAddr(), Result.FirstAddr(), Size);

	}


	/*复数矩阵点乘*/
	CCMat DotMul(CCMat& Data1, CCMat& Data2, cudaStream_t stream)
	{
		CCMat Result(Data1.dims(0), Data1.dims(1), Data1.dims(2));
		unsigned int Data1Size = Data1.elements();
		unsigned int Data2Size = Data2.elements();
		if (Data1Size != Data2Size)
		{
			printf("维度不一致, Data1: %d, Data2: %d\n", Data1Size, Data2Size);
			throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		DotMulKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2.FirstAddr(), Result.FirstAddr(), Data1Size);

		return Result;
	}

	/*复数矩阵点乘*/
	void StillDotMul(CCMat& Data1, CCMat& Data2, cudaStream_t stream)
	{
		unsigned int Data1Size = Data1.elements();
		unsigned int Data2Size = Data2.elements();
		if (Data1Size != Data2Size)
		{
			printf("维度不一致, Data1: %d, Data2: %d\n", Data1Size, Data2Size);
			throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		DotMul2Kernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2.FirstAddr(), Data1Size);

	}

	/*复数矩阵加法*/
	CCMat MatAdd(CCMat& Data1, CCMat& Data2, cudaStream_t stream)
	{
		CCMat Result(Data1.dims(0), Data1.dims(1), Data1.dims(2));
		unsigned int Data1Size = Data1.elements();
		unsigned int Data2Size = Data2.elements();
		if (Data1Size != Data2Size)
		{
			printf("复数矩阵加法 时 维度不一致\n"); throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		MatAddKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2.FirstAddr(), Result.FirstAddr(), Data1Size);

		return Result;
	}

	/*复数矩阵加法*/
	void MatAdd(CCMat& Data1, CCMat& Data2, CCMat& Result, cudaStream_t stream)
	{
		Result.Resize(Data1.dims(0), Data1.dims(1), Data1.dims(2));

		unsigned int Data1Size = Data1.elements();
		unsigned int Data2Size = Data2.elements();
		if (Data1Size != Data2Size)
		{
			printf("复数矩阵加法 时 维度不一致\n"); throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		MatAddKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2.FirstAddr(), Result.FirstAddr(), Data1Size);

	}

	/*复数矩阵加法*/
	void MatAdd(CCMat& Data1, CCMat& Data2, CCMat& Data3, CCMat& Result, cudaStream_t stream)
	{
		Result.Resize(Data1.dims(0), Data1.dims(1), Data1.dims(2));

		unsigned int Data1Size = Data1.elements();
		unsigned int Data2Size = Data2.elements();
		if (Data1Size != Data2Size)
		{
			printf("复数矩阵加法 时 维度不一致\n"); throw;
		}

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		MatAddKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2.FirstAddr(), Result.FirstAddr(), Data1Size);

	}

	/*复数矩阵加法*/
	void MatAdd(CCMat& Data1, cuComplex Data2, cudaStream_t stream)
	{

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1.elements() + threadsPerBlock - 1) / threadsPerBlock;

		MatAddComplexKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Data1.elements());

	}

	/*float矩阵乘法*/
	FCMat MatMul(FCMat& A, FCMat& B, cudaStream_t stream)
	{
		if (A.dims(1) != B.dims(0))
		{
			printf("维度不一致");throw;
		}

		FCMat C(A.dims(0), B.dims(1), 1);

		int A_Row = A.dims(0);
		int A_Col = A.dims(1);
		int B_Row = B.dims(0);
		int B_Col = B.dims(1);

		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		float a = 1, b = 0;
		cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, A_Row, B_Col, A_Col, &a, A.FirstAddr(), A_Row, B.FirstAddr(), B_Row, &b, C.FirstAddr(), A_Row);

		cublasDestroy(handle);

		return C;
	}

	
	/*复数矩阵乘法*/
	CCMat MatMul(CCMat& A, CCMat& B, cudaStream_t stream)
	{
		if (A.dims(1) != B.dims(0))
		{
			printf("维度不一致");throw;
		}

		CCMat C(A.dims(0), B.dims(1), 1);

		int A_Row = A.dims(0);
		int A_Col = A.dims(1);
		int B_Row = B.dims(0);
		int B_Col = B.dims(1);

		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cuComplex a,b;
		a.x = 1;a.y = 0;
		b.x = 0;b.y = 0;
		cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, A_Row, B_Col, A_Col, &a, A.FirstAddr(), A_Row, B.FirstAddr(), B_Row, &b, C.FirstAddr(), A_Row);

		cublasDestroy(handle);
		return C;
	}

	/*复数矩阵乘法*/
	void MatMul(CCMat& A, CCMat& B, CCMat& C, cudaStream_t stream)
	{
		if (A.dims(1) != B.dims(0))
		{
			printf("维度不一致"); throw;
		}

		C.Resize(A.dims(0), B.dims(1), 1);

		int A_Row = A.dims(0);
		int A_Col = A.dims(1);
		int B_Row = B.dims(0);
		int B_Col = B.dims(1);

		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cuComplex a, b;
		a.x = 1; a.y = 0;
		b.x = 0; b.y = 0;
		cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, A_Row, B_Col, A_Col, &a, A.FirstAddr(), A_Row, B.FirstAddr(), B_Row, &b, C.FirstAddr(), A_Row);

		cublasDestroy(handle);

	}

	/*复数矩阵除float*/
	CCMat CFDiv(CCMat& Data1, float Data2, cudaStream_t stream)
	{
		CCMat Result(Data1.dims(0), Data1.dims(1), Data1.dims(2));
		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		C2FDivKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Result.FirstAddr(), Data1Size);

		return Result;
	}

	/*原地复数矩阵除float*/
	void StillCFDiv(CCMat& Data1, float Data2, cudaStream_t stream)
	{
		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		C2FDiv2Kernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Data1Size);

	}



	/*复数矩阵乘float*/
	void MatFMul(CCMat& Data1, FCMat& Data2, cudaStream_t stream)
	{
		if (Data1.elements() != Data2.elements()) return;

		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		MatMulKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2.FirstAddr(), Data1Size);

	}

	/*复数矩阵乘float*/
	CCMat MatFMul(CCMat& Data1, float Data2, cudaStream_t stream)
	{
		CCMat Result(Data1.dims(0), Data1.dims(1), Data1.dims(2));
		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		C2FMulKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Result.FirstAddr(), Data1Size);

		return Result;
	}

	/*复数矩阵乘float*/
	void MatFMul(CCMat& Data1, CCMat& Result, float Data2, cudaStream_t stream)
	{
		Result.Resize(Data1.dims(0), Data1.dims(1), Data1.dims(2));
		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		C2FMulKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Result.FirstAddr(), Data1Size);

	}

	/*复数矩阵乘float*/
	void StillMatFMul(CCMat& Data1, float Data2, cudaStream_t stream)
	{
		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		C2FMulKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Data1.FirstAddr(), Data1Size);

	}


	/*float矩阵乘float*/
	FCMat MatFMul(FCMat& Data1, float Data2, cudaStream_t stream)
	{
		FCMat Result(Data1.dims(0), Data1.dims(1), Data1.dims(2));
		unsigned int Data1Size = Data1.elements();

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data1Size + threadsPerBlock - 1) / threadsPerBlock;

		F2FMulKernal << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data1.FirstAddr(), Data2, Result.FirstAddr(), Data1Size);

		return Result;
	}


	/*截取深拷贝*/
	CCMat Slice(CCMat& Data, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		CCMat Result(REnd - RStart + 1, CEnd - CStart + 1, BEnd - BStart + 1);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows*Cols*Bands + threadsPerBlock - 1) / threadsPerBlock;

		SliceKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Result.FirstAddr(), RStart, REnd, CStart, CEnd, BStart, BEnd, Rows, Cols, Bands);

		return Result;
	}

	/*截取深拷贝(FCMat)*/
	FCMat Slice(FCMat& Data, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		FCMat Result(REnd - RStart + 1, CEnd - CStart + 1, BEnd - BStart + 1);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows*Cols*Bands + threadsPerBlock - 1) / threadsPerBlock;

		SliceKernel << <blocksPerGrid, threadsPerBlock >> > (Data.FirstAddr(), Result.FirstAddr(), RStart, REnd, CStart, CEnd, BStart, BEnd, Rows, Cols, Bands);

		return Result;
	}

	/*截取深拷贝(FCMat)*/
	void Slice(FCMat& Data, FCMat& Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		Result.Resize(REnd - RStart + 1, CEnd - CStart + 1, BEnd - BStart + 1);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows*Cols*Bands + threadsPerBlock - 1) / threadsPerBlock;

		SliceKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Result.FirstAddr(), RStart, REnd, CStart, CEnd, BStart, BEnd, Rows, Cols, Bands);

	}

	/*截取深拷贝*/
	void Slice(CCMat& Data, CCMat& Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		Result.Resize(REnd - RStart + 1, CEnd - CStart + 1, BEnd - BStart + 1);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows*Cols*Bands + threadsPerBlock - 1) / threadsPerBlock;

		SliceKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Result.FirstAddr(), RStart, REnd, CStart, CEnd, BStart, BEnd, Rows, Cols, Bands);

	}


	/*矩阵按列扩展
		描述：将Data扩展为Result长度（按列扩展）
		Data : (Rows, Cols, Bands)
		Result : (Rows + abs(AddLength), Cols, Bands)
		AddLength : 扩展量（正右负左）
		Value : 扩展的值
	*/
	FCMat cuExtend(FCMat& Data, int AddLength, float Value, cudaStream_t stream)
	{
		if (AddLength == 0) return Data;

		FCMat Result(Data.dims(0) + abs(AddLength), Data.dims(1), Data.dims(2));
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Result.elements() + threadsPerBlock - 1) / threadsPerBlock;

		ExtendKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Result.FirstAddr(), AddLength, Value, Data.dims(0), Data.dims(1), Data.dims(2));

		return Result;
	}

	/*矩阵按列扩展
	描述：将Data扩展为Result长度（按列扩展）
	Data : (Rows, Cols, Bands)
	Result : (Rows + abs(AddLength), Cols, Bands)
	AddLength : 扩展量（正右负左）
	Value : 扩展的值
	*/
	CCMat cuExtend(CCMat& Data, int AddLength, cuComplex Value, cudaStream_t stream)
	{
		if (AddLength == 0) return Data;

		CCMat Result(Data.dims(0) + abs(AddLength), Data.dims(1), Data.dims(2));
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Result.elements() + threadsPerBlock - 1) / threadsPerBlock;

		ExtendKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Result.FirstAddr(), AddLength, Value, Data.dims(0), Data.dims(1), Data.dims(2));

		return Result;
	}

	void cuExtend(CCMat& Data, CCMat& Result, int AddLength, cuComplex Value, cudaStream_t stream)
	{
		if (AddLength == 0)
		{
			Result.CopyFromDevice(Data);
			return;
		}

		Result.Resize(Data.dims(0) + abs(AddLength), Data.dims(1), Data.dims(2));
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Result.elements() + threadsPerBlock - 1) / threadsPerBlock;

		ExtendKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Result.FirstAddr(), AddLength, Value, Data.dims(0), Data.dims(1), Data.dims(2));

	}

	// 转置
	void Transpose(FCMat& data, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		FCMat result(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		TransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		cudaDeviceSynchronize();
		data = result;

		result.Free();
	};

	// 转置
	FCMat cuTranspose(FCMat& data, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		FCMat result(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		TransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		
		return result;
	};

	// 转置
	void Transpose(CCMat& data, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		CCMat result(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		TransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		
		data.CopyFromDevice(result);

		result.Free();
	};

	// 转置
	CCMat cuTranspose(CCMat& data, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		CCMat result(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		TransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		
		return result;
	};
void Transpose(CCMat& data, CCMat& result, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		result.Resize(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		TransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);

		/*data.CopyFromDevice(result);

		result.Free();*/
	};

	// 共轭转置
	void CTranspose(CCMat& data, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		CCMat result(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		CTransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		cudaDeviceSynchronize();
		data = result;

		result.Free();
	};

	// 共轭转置
	CCMat cuCTranspose(CCMat& data, cudaStream_t stream)
	{
		int Rows = data.dims(0);
		int Cols = data.dims(1);
		int Bands = data.dims(2);

		CCMat result(Cols, Rows, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;

		CTransposeKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);

		return result;
	};

	// 求和
	float cuSum(FCMat& Data, cudaStream_t stream)
	{
		float result;
		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cublasSasum(handle, Data.elements(), Data.FirstAddr(), 1, &result);
		cublasDestroy(handle);
		return result;
	}

	// 复数abs求和
	float cuSum(CCMat& Data, cudaStream_t stream)
	{
		float result;
		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cublasScasum(handle, Data.elements(), Data.FirstAddr(), 1, &result);
		cublasDestroy(handle);
		return result;
	}

	// 获取最大值
	float cuMax(FCMat& Data, cudaStream_t stream)
	{
		int idx;
		int Rows =  Data.dims(0);
		int Cols =  Data.dims(1);
		int Bands = Data.dims(2);

		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cublasIsamax(handle, Data.elements(), Data.FirstAddr(), 1, &idx);
		cublasDestroy(handle);

		// 阵索引
		unsigned int B_bias = idx / (Rows * Cols);
		// 列索引
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// 行索引
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;
		
		float result;
		// 获取单个元素的值
		cudaMemcpy(&result, Data.FirstAddr() + idx, sizeof(float), cudaMemcpyDeviceToHost);

		return result;
	}

	// 获取最大值及下标
	void cuMax(FCMat& Data, float& result, dim3& index, cudaStream_t stream)
	{
		int idx;
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cublasIsamax(handle, Data.elements(), Data.FirstAddr(), 1, &idx);
		cublasDestroy(handle);

		idx = idx - 1;

		// 阵索引
		index.z = idx / (Rows * Cols);
		// 列索引
		index.y = (idx % (Rows * Cols)) / Rows;
		// 行索引
		index.x = (idx % (Rows * Cols)) % Rows;

		// 获取单个元素的值
		cudaMemcpy(&result, Data.FirstAddr() + idx, sizeof(float), cudaMemcpyDeviceToHost);
	}

	// 获取最大值及下标
	void cuMax(FCMat& Data, FCMat& result, dim3& index, cudaStream_t stream)
	{
		int idx;
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		cublasHandle_t handle;
		cublasCreate(&handle);
		cublasSetStream(handle, stream);
		cublasIsamax(handle, Data.elements(), Data.FirstAddr(), 1, &idx);
		cublasDestroy(handle);

		idx = idx - 1;

		// 阵索引
		index.z = idx / (Rows * Cols);
		// 列索引
		index.y = (idx % (Rows * Cols)) / Rows;
		// 行索引
		index.x = (idx % (Rows * Cols)) % Rows;

		// 获取单个元素的值
		result.CopyFromDevicePoint(Data.FirstAddr() + idx, 1, 1, 1);
	}

	/*float矩阵按列求和*/
	FCMat sum(FCMat& Data, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		FCMat result(1, Cols, Bands, 0.0);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		cuSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		return result;
	}

	/*float矩阵按列求和*/
	void sum(FCMat& Data, FCMat& result, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		result.Resize(1, Cols, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(result.FirstAddr(), 0.0, Cols * Bands);

		blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		cuSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);

	}

	/*复数矩阵按列求和*/
	CCMat sum(CCMat& Data, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		cuComplex temp; temp.x = 0.0; temp.y = 0.0;
		CCMat result(1, Cols, Bands, temp);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		cuSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		return result;
	}
	/*复数矩阵按列求和*/
	void sum(CCMat& Data, CCMat& result, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		cuComplex temp; temp.x = 0.0; temp.y = 0.0;
		result.Resize(1, Cols, Bands, temp, stream);
		
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		ValuateKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> >(result.FirstAddr(), temp, Cols * Bands);

		blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		cuSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
	}
	/*实数矩阵按列求和*/
	ICMat sum(ICMat& Data, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		ICMat result(1, Cols, Bands);

		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		cuSumKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), result.FirstAddr(), Rows, Cols, Bands);
		return result;
	}

	// 归一化
	void Z_ScoreStandardization(FCMat &data, FCMat &z_score, int N, cudaStream_t stream)
	{
		FCMat SumValue(1, 1, 1, 0.0, stream);
		sum(data, SumValue, stream);
		FCMat atiledata(N, 1, 1);

		int threadsPerBlock = 1024;
		unsigned int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
		FCMat powatilemean(N, 1, 1);
		cuPowATileMean << <blocksPerGrid, threadsPerBlock, 0, stream>> >(data.FirstAddr(), atiledata.FirstAddr(), SumValue.FirstAddr(), powatilemean.FirstAddr(), N);

		FCMat Sum(1, 1, 1, 0.0, stream);
		sum(powatilemean, Sum, stream);

		z_score.Resize(N, 1, 1);
		//以上为了获取Sum，以便后面求sigma
		SqrtTileDiv << <blocksPerGrid, threadsPerBlock, 0, stream>> >(Sum.FirstAddr(), atiledata.FirstAddr(), z_score.FirstAddr(), N);

		//释放空间
		SumValue.Free();
		atiledata.Free();
		powatilemean.Free();
		Sum.Free();
	}

	void Z_ScoreStandardization(CCMat &data, CCMat &z_score, int N, cudaStream_t stream)
	{
		//复数求和
		cuComplex initValue;	initValue.x = 0.0;	initValue.y = 0.0;
		CCMat SumValue(1,1,1, initValue, stream);
		sum(data, SumValue, stream);

		CCMat atiledata(N, 1, 1);
		int threadsPerBlock = 1024;
		unsigned int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
		FCMat PowerResult(N, 1, 1);

		cuComplexPowATileMean << <blocksPerGrid, threadsPerBlock, 0, stream >> >(data.FirstAddr(), atiledata.FirstAddr(), SumValue.FirstAddr(), PowerResult.FirstAddr(), N);

		FCMat Sum(1,1,1,0.0,stream);
		sum(PowerResult, Sum, stream);

		z_score.Resize(N, 1, 1);
		//以上为了获取Sum，以便后面求sigma
		SqrtTileComplexDiv << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Sum.FirstAddr(), atiledata.FirstAddr(),
					z_score.FirstAddr(), N);

		//释放空间
		SumValue.Free();
		atiledata.Free();
		PowerResult.Free();
		Sum.Free();
	
		/******************************************************/
	}

	void Z_ScoreStandardization_con(CCMat &data, CCMat &z_score, int N, cudaStream_t stream)
	{
		//cuComplex initialValue;	initialValue.x = 0.0; initialValue.y = 0.0;
		z_score.Resize(data.dims(0), data.dims(1), data.dims(2));
		unsigned int Rows = data.dims(0);
		unsigned int Cols = data.dims(1);
		unsigned int Bands = data.dims(2);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Rows * Cols * Bands + threadsPerBlock - 1) / threadsPerBlock;
		cuComplex temp; temp.x = 0.0, temp.y = 0.0;
		CCMat Sum_first(1, 1, 1, temp, stream);
		CCMat atiledata(N, 1, 1);
		FCMat PowerResult(N, 1, 1);
		FCMat Sum_second(1, 1, 1, temp.x, stream);
		z_score_kernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (data.FirstAddr(), atiledata.FirstAddr(), PowerResult.FirstAddr(), Rows, Cols, Bands, Sum_first.FirstAddr(), Sum_second.FirstAddr(), z_score.FirstAddr());

		Sum_first.Free();
		atiledata.Free();
		PowerResult.Free();
		Sum_second.Free();
	}

	// FFT
	/*扩展到FFTNum后进行FFT 按列*/
	CCMat cuFFT(CCMat &Data, unsigned int FFTNum, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		CCMat DataF;
		cuComplex Value; Value.x = 0.0; Value.y = 0.0;
		// 需要扩展再扩展
		if (FFTNum > Rows)
		{
			// 按fft点数扩展（按列）
			cuExtend(Data, DataF, FFTNum - Rows, Value, stream);

		}
		else
		{
			DataF.CopyFromDevice(Data);
			FFTNum = Rows;
		}

		cufftHandle Handle;
		cufftPlan1d(&Handle, FFTNum, CUFFT_C2C, Cols);
		cufftSetStream(Handle, stream);
		for (int i = 0; i < Bands; i++)
		{
			cufftExecC2C(Handle, DataF.FirstAddr() + i * (FFTNum * Cols), DataF.FirstAddr() + i * (FFTNum * Cols), CUFFT_FORWARD);
		}
		cufftDestroy(Handle);

		return DataF;
	}
	void cuFFT(CCMat &Data, CCMat& DataF, unsigned int FFTNum, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		
		cuComplex Value; Value.x = 0.0; Value.y = 0.0;
		// 需要扩展再扩展
		//if (FFTNum > Rows)
		//{
		//	// 按fft点数扩展（按列）
		//	cuExtend(Data, DataF, FFTNum - Rows, Value, stream);

		//}
		//else
		//{
		//	DataF.CopyFromDevice(Data);
		//	FFTNum = Rows;
		//}

		cuExtend(Data, DataF, FFTNum - Rows, Value, stream);

		//CCMat DataTemp(DataF.dims(0), DataF.dims(1), DataF.dims(2));

		cufftHandle Handle;
		cufftPlan1d(&Handle, FFTNum, CUFFT_C2C, Cols);
		cufftSetStream(Handle, stream);
		for (int i = 0; i < Bands; i++)
		{
			//cufftExecC2C(Handle, DataF.FirstAddr() + i * (FFTNum * Cols), DataTemp.FirstAddr() + i * (FFTNum * Cols), CUFFT_FORWARD);
			cufftExecC2C(Handle, DataF.FirstAddr() + i * (FFTNum * Cols), DataF.FirstAddr() + i * (FFTNum * Cols), CUFFT_FORWARD);
		}
		cufftDestroy(Handle);

		//DataF.CopyFromDevice(DataTemp);

		//DataTemp.Free();


	}

	// IFFT
	/*IFFT后再截断到OriginNum 按列*/
	void cuIFFT(CCMat &Data, unsigned int OriginNum, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);

		CCMat DataTemp(Data.dims(0), Data.dims(1), Data.dims(2));

		cufftHandle plan;
		cufftPlan1d(&plan, Rows, CUFFT_C2C, Cols);
		cufftSetStream(plan, stream);
		for (int i = 0; i < Bands; i++)
		{
			cufftExecC2C(plan, Data.FirstAddr() + i * (Rows * Cols), DataTemp.FirstAddr() + i * (Rows * Cols), CUFFT_INVERSE);
		}
		cufftDestroy(plan);

		//cuComplex Value; Value.x = 0.0; Value.y = 0.0;

		StillMatFMul(DataTemp, 1.0 / float(Rows), stream);

		
		// 需要截断再截断	
		//if (OriginNum < Rows)
		//{
		//	// 按OriginNum点数截断（按列）
		//	Slice(DataTemp, Data, 0, OriginNum-1, 0, Data.dims(1)-1, 0, Data.dims(2)-1, stream);

		//}
		//else
		//{
		//	Data.CopyFromDevice(DataTemp);
		//	OriginNum = Rows;
		//}		

		Slice(DataTemp, Data, 0, OriginNum - 1, 0, Data.dims(1) - 1, 0, Data.dims(2) - 1, stream);

		DataTemp.Free();

	}
	//输入输出不为同一变量的
	// IFFT
	/*IFFT后再截断到OriginNum 按列*/
	void cuIFFT(CCMat &Data, CCMat &Data1, unsigned int OriginNum, int Cols, int Bands,  cudaStream_t stream)
	{
		//int Rows = Data.dims(0);
		//int Cols = Data.dims(1);
		//int Bands = Data.dims(2);

		//CCMat DataTemp(Data.dims(0), Data.dims(1), Data.dims(2));

		cufftHandle plan;
		cufftPlan1d(&plan, OriginNum, CUFFT_C2C, Cols);
		cufftSetStream(plan, stream);
		for (int i = 0; i < Bands; i++)
		{
			cufftExecC2C(plan, Data.FirstAddr() + i * (OriginNum * Cols), Data1.FirstAddr() + i * (OriginNum * Cols), CUFFT_INVERSE);
		}
		cufftDestroy(plan);

		StillMatFMul(Data1, 1.0 / float(OriginNum), stream);

		Slice(Data1, Data,  0, OriginNum - 1, 0, Data.dims(1) - 1, 0, Data.dims(2) - 1, stream);

		//DataTemp.Free();

	}
	// 矩阵共轭
	void cuConj(CCMat& Data, cudaStream_t stream)
	{
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data.elements() + threadsPerBlock - 1) / threadsPerBlock;

		ConjKernel << <blocksPerGrid, threadsPerBlock, 0, stream >> > (Data.FirstAddr(), Data.elements());
	}

	// 按列循环移动，shiftNum大于0，则右移，否则左移
	void cuShiftCol(CCMat& Data, CCMat& Result, int shiftNum, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data.elements() + threadsPerBlock - 1) / threadsPerBlock;
		CirculShiftCol << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data.FirstAddr(), Result.FirstAddr(), shiftNum, Rows, Cols, Bands);
	}

	// 按列循环移动，shiftNum大于0，则右移，否则左移
	void cuShiftCol(FCMat& Data, FCMat& Result, int shiftNum, cudaStream_t stream)
	{
		int Rows = Data.dims(0);
		int Cols = Data.dims(1);
		int Bands = Data.dims(2);
		int threadsPerBlock = 1024;
		int blocksPerGrid = (Data.elements() + threadsPerBlock - 1) / threadsPerBlock;
		CirculShiftCol << <blocksPerGrid, threadsPerBlock, 0, stream >> >(Data.FirstAddr(), Result.FirstAddr(), shiftNum, Rows, Cols, Bands);
	}

	//// 求面积中心
	//double CenterOfMass(Vector<double> envelope) {
	//	//求向量面积中心
	//	double Mass = 0;
	//	double area = 0;
	//	double areaLast = 0;
	//	double center = 0;
	//	int flag = 0;
	//	//cout << envelope[0] << " " << envelope[1] << " " << envelope[2] << " " << envelope[3] << " " << envelope[4] << endl;
	//	//envelope[0] = envelope[0] - 0.1;
	//	if (envelope.size() == 1)
	//		center = 0;
	//	else {
	//		Mass = sum(envelope) - envelope[0] / 2 - envelope[envelope.size() - 1] / 2; // 总面积
	//																					//cout << Mass << endl;
	//		for (int ii = 0; (ii < envelope.size() - 1) && (flag == 0); ii++)
	//		{
	//			area += (envelope[ii] + envelope[ii + 1]) / 2;
	//			areaLast = area - (envelope[ii] + envelope[ii + 1]) / 2;
	//			if (area == Mass / 2) {
	//				//面积中心在envelop[ii+1]处
	//				center = ii + 1;
	//			}
	//			// 			cout << Mass / 2 << endl;
	//			// 			cout << Mass / 2 - envelope[ii] << endl;
	//			if (area > Mass / 2) {
	//				//面积中心在envelop[ii] 与 envelop[ii+1]之间
	//				if (abs(envelope[ii + 1] - envelope[ii]) < 1e-6)
	//				{
	//					center = ii + 0.5;
	//				}
	//				else
	//				{
	//					/*cout << "L2-L1:" << envelope[ii + 1] - envelope[ii] << endl;
	//					cout << "L1^2:"<< pow(envelope[ii], 2) << endl;
	//					cout << "P:" << (Mass / 2 - envelope[ii]) << endl;
	//					cout << "根号下:" << sqrt(pow(envelope[ii], 2) + 2 * (envelope[ii + 1] - envelope[ii])*(Mass / 2 - envelope[ii])) << endl;*/
	//					center = ii + (sqrt(pow(envelope[ii], 2) + 2 * (envelope[ii + 1] - envelope[ii])*(Mass / 2 - areaLast)) - envelope[ii]) / (envelope[ii + 1] - envelope[ii]);
	//				}
	//				//ii = envelope.size()-1;
	//				flag = 1;
	//			}

	//		}
	//	}
	//	return center;
	//}


}

