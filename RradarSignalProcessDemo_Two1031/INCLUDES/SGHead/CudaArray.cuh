#ifndef  _CudaArray_H_
#define  _CudaArray_H_

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>
#include <string>
#include "CommonKernel.cuh"
#include <cublas.h>
#include "CudaArrayTool.cuh"
#include "Error.cuh"

template <typename Type>
class SigHostMatrix;
template <typename Type>
class CudaArray;

typedef CudaArray<cuComplex>   CCMat;	// 复数类型矩阵
typedef CudaArray<float>	   FCMat;	// 浮点类型矩阵
typedef CudaArray<bool>		   BCMat;	// 布尔类型矩阵
typedef CudaArray<int>		   ICMat;	// 整型矩阵

// 信号设置 GPU设备中的信号
template <typename Type>
class CudaArray
{
public:
	CudaArray()
	{
		m_Rows = 0;
		m_Cols = 0;
		m_Bands = 0;
		d_Data = NULL;
		m_stream = NULL;
	}
	 // 非深度拷贝过程，仅传递首地址指针和尺度参数
	CudaArray(CudaArray<Type> &m)
	{
		err = CopyFromDevicePoint(m.FirstAddr(), m.Row(), m.Col(), m.Band());
	}

	// 定义大小构造
	CudaArray(unsigned int Rows)
	{
		m_Rows = Rows;
		m_Cols = 1;
		m_Bands = 1;
		d_Data = NULL;
		m_stream = NULL;
		if (d_Data == NULL)
			err = Malloc(Rows, m_Cols, m_Bands);
		else
		{
			Free();
			err = Malloc(Rows, m_Cols, m_Bands);
		}
	}

	// 定义大小构造
	CudaArray(unsigned int Rows, unsigned int Cols)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = 1;
		d_Data = NULL;
		m_stream = NULL;
		if (d_Data == NULL)
			err = Malloc(Rows, Cols, m_Bands);
		else
		{
			Free();
			err = Malloc(Rows, Cols, m_Bands);
		}
	}

	// 定义大小构造
	CudaArray(unsigned int Rows, unsigned int Cols, unsigned int Bands)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = Bands;
		d_Data = NULL;
		m_stream = NULL;
		if (d_Data == NULL)
			err = Malloc(Rows, Cols, Bands);
		else
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
		}
	}


	// 定义大小+float赋值构造
	CudaArray(unsigned int Rows, unsigned int Cols, unsigned int Bands, float Value)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = Bands;
		d_Data = NULL;
		m_stream = NULL;

		if (d_Data == NULL)
			err = Malloc(Rows, Cols, Bands);
		else
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
		}
		
		// 将值赋为Value
		Valuate(this->d_Data, Value, m_Rows * m_Cols * m_Bands, m_stream);
	}

	// 定义大小+float赋值构造
	CudaArray(unsigned int Rows, unsigned int Cols, unsigned int Bands, float Value, cudaStream_t stream)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = Bands;
		d_Data = NULL;
		m_stream = stream;

		if (d_Data == NULL)
			err = Malloc(Rows, Cols, Bands);
		else
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
		}

		// 将值赋为Value
		Valuate(this->d_Data, Value, m_Rows * m_Cols * m_Bands, stream);
	}

	// 定义大小+cuComplex赋值构造
	CudaArray(unsigned int Rows, unsigned int Cols, unsigned int Bands, cuComplex Value)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = Bands;
		d_Data = NULL;
		m_stream = NULL;
		if (d_Data == NULL)
			err = Malloc(Rows, Cols, Bands);
		else
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
		}

		// 将值赋为Value
		Valuate(this->d_Data, Value, m_Rows * m_Cols * m_Bands, m_stream);
	}

	// 定义大小+cuComplex赋值构造
	CudaArray(unsigned int Rows, unsigned int Cols, unsigned int Bands, cuComplex Value, cudaStream_t stream)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = Bands;
		d_Data = NULL;
		m_stream = stream;
		if (d_Data == NULL)
			err = Malloc(Rows, Cols, Bands);
		else
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
		}

		// 将值赋为Value
		Valuate(this->d_Data, Value, m_Rows * m_Cols * m_Bands, m_stream);
	}


	~CudaArray(){ }

	// 释放内存
	void Free()
	{
		m_Rows = 0;
		m_Cols = 0;
		m_Bands = 0;
		if (d_Data != NULL)
		{
			err = cudaFree(d_Data);
			if (err != cudaSuccess)
			{
				printf("释放内存不成功!错误为:%s\n", cudaGetErrorString(err));
			}
		}
		d_Data = NULL;
	}

	// 清除原来的数据，重新分配
	cudaError Resize(unsigned int Rows, unsigned int Cols, unsigned int Bands)
	{
		err = cudaSuccess;
		if (m_Rows != Rows || m_Cols != Cols || m_Bands != Bands)
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
		}
		return err;
	}
	// 清除原来的数据，重新分配
	cudaError Resize(unsigned int Rows, unsigned int Cols, unsigned int Bands, Type Value, cudaStream_t stream)
	{
		err = cudaSuccess;
		if (m_Rows != Rows || m_Cols != Cols || m_Bands != Bands)
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
			if (err == cudaSuccess)
			{
				// 将值赋为Value
				Valuate(this->d_Data, Value, Rows * Cols * Bands, stream);
			}
		}
		return err;
	}


	// 获取首地址
	Type* FirstAddr() { return d_Data; };
	// 获取行数
	unsigned int Row() { return m_Rows; };
	// 获取列数
	unsigned int Col() { return m_Cols; };
	// 获取阵数
	unsigned int Band() { return m_Bands; };
	// 获取元素个数
	unsigned int elements() { return m_Rows * m_Cols * m_Bands; };

	unsigned int dims(int dims)
	{
		if (dims == 0) return m_Rows;
		else if (dims == 1) return m_Cols;
		else if (dims == 2) return m_Bands;
		else return -1;
	}

	// 从主机拷贝数据到设备中
	cudaError CopyFromHost(Type * h_data, unsigned int Rows, unsigned int Cols, unsigned int Bands)
	{
		if (m_Rows != Rows || m_Cols != Cols || m_Bands != Bands)
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
			if(err == cudaSuccess) 
				err = cudaMemcpy(d_Data, h_data, sizeof(Type) * Rows * Cols * Bands, cudaMemcpyHostToDevice);
		}
		else
			err = cudaMemcpy(d_Data, h_data, sizeof(Type) * Rows * Cols * Bands, cudaMemcpyHostToDevice);

		return err;
	}


	// 从主机拷贝数据到设备中
	cudaError CopyFromHost(SigHostMatrix<Type> &h_data)
	{
		err = CopyFromDevicePoint(h_data.FirstAddr(), h_data.Row(), h_data.Col(), h_data.Bands());
		return err;
	}

	// 从设备拷贝数据到设备中
	cudaError CopyFromDevicePoint(Type* d_data, unsigned int Rows, unsigned int Cols, unsigned int Bands)
	{
		if (m_Rows != Rows || m_Cols != Cols || m_Bands != Bands)
		{
			Free();
			err = Malloc(Rows, Cols, Bands);
			if (err == cudaSuccess)
				err = cudaMemcpy(d_Data, d_data, sizeof(Type) * Rows * Cols * Bands, cudaMemcpyDeviceToDevice);
		}
		else
			err = cudaMemcpy(d_Data, d_data, sizeof(Type)  * Rows * Cols * Bands, cudaMemcpyDeviceToDevice);

		return err;
	}
	// 从设备拷贝数据到设备中
	cudaError CopyFromDevice(CudaArray<Type> &m)
	{
		//unsigned int R = m.Row();
		//unsigned int C = m.Col();
		//Type* d_data = m.FirstAddr();
		err = CopyFromDevicePoint(m.FirstAddr(), m.Row(), m.Col(), m.Band());
		return err;
	}

	cudaError FromArrayToDevice(Type* DataAddr, unsigned int Rows, unsigned int Cols, unsigned int Bands)
	{
		Free();
		this->d_Data = DataAddr;
		this->m_Rows = Rows;
		this->m_Cols = Cols;
		this->m_Bands = Bands;
		err = cudaSuccess;
		return err;
	}
	// 深度拷贝
	CudaArray<Type> operator=(CudaArray<Type> &m)
	{
		if (m.FirstAddr() == FirstAddr())
			return *this;
		
		// this->d_Data = NULL;
		err = CopyFromDevicePoint(m.FirstAddr(), m.Row(), m.Col(), m.Band());
		m.Free();
		return *this;
	}



	//// 复数点乘
	//CudaArray<Type> operator*(CudaArray<Type> &m)
	//{
	//	if (m.FirstAddr() == FirstAddr())
	//		return *this;

	//	err = CopyFromDevicePoint(m.FirstAddr(), m.Row(), m.Col(), m.Band());
	//	return *this;
	//}


	// 获取某一列, 非深度拷贝
	CudaArray<Type> getCol(unsigned int Col)
	{
		Type* result = this->d_Data + m_Rows * Col;
		CudaArray<Type> resultM(m_Rows, 1, 1);
		resultM.FromArrayToDevice(result, m_Rows, 1, 1);

		return resultM;
	}

	// 获取某一列, 深度拷贝
	CudaArray<Type> getColCopy(unsigned int Col)
	{
		Type* result = this->d_Data + m_Rows * Col;
		CudaArray<Type> resultM(m_Rows, 1, 1);
		resultM.CopyFromDevicePoint(result, m_Rows, 1, 1);

		return resultM;
	}

	// 获取某一阵, 非深度拷贝
	void getBand(CudaArray<Type>& resultM, unsigned int Band)
	{
		Type* result = this->d_Data + m_Rows * m_Cols * Band;
		resultM.Resize(m_Rows, m_Cols, 1);
		resultM.FromArrayToDevice(result, m_Rows, m_Cols, 1);

	}

	void getBandCopy(CudaArray<Type>& resultM, unsigned int Band)
	{
		Type* result = this->d_Data + m_Rows * m_Cols * Band;
		resultM.Resize(m_Rows, m_Cols, 1);
		resultM.CopyFromDevicePoint(result, m_Rows, m_Cols, 1);

	}

	/* 获取某个元素(Row, Col, Band)到主机 */
	Type getValueToHost(unsigned int Row, unsigned int Col, unsigned int Band, cudaStream_t stream)
	{
		if (Row >= m_Rows || Col >= m_Cols || Band >= m_Bands)
		{
			printf("超出维度\n");throw;
		}

		Type result;
		// 获取偏移
		Type* bias = this->d_Data + Band * m_Rows * m_Cols + Col * m_Rows + Row;
		// 获取单个元素的值
		cudaMemcpy(&result, bias, sizeof(Type), cudaMemcpyDeviceToHost);

		return result;
	}

	/* 获取某个元素(Row, Col, Band)到GPU*/
	CudaArray<Type> getValue(unsigned int Row, unsigned int Col, unsigned int Band, cudaStream_t stream)
	{
		if (Row >= m_Rows || Col >= m_Cols || Band >= m_Bands)
		{
			printf("超出维度\n");throw;
			return;
		}

		CudaArray<Type> Result;
		Result.FromArrayToDevice(this->d_Data + Band * m_Rows * m_Cols + Col * m_Rows + Row, 1, 1, 1);

		return Result;
	}

	// 设置某一行
	void setRow(unsigned int Row, unsigned int CurrentBand, CudaArray<Type>& Value, cudaStream_t stream)
	{
		if (Value.elements() != m_Cols || CurrentBand >= m_Bands)
		{
			printf("维度不一致\n");throw;
			return;
		}
		cuSetRow(this->d_Data, Value.FirstAddr(), CurrentBand, m_Rows, m_Cols, m_Bands, Row, stream);
	}

	// 设置某几行
	void setRow(CudaArray<int>& Row, unsigned int CurrentBand, CudaArray<Type>& Value, cudaStream_t stream)
	{
		if (Value.dims(1) != m_Cols || Value.dims(0) >= m_Rows || CurrentBand >= m_bands)
		{
			printf("维度不一致\n");throw;
			return;
		}
		cuSetRow(this->d_Data, Value.FirstAddr(), CurrentBand, m_Rows, m_Cols, m_Bands, Row, stream);
	}
	
	// 设置某一列
	void setCol(unsigned int Col, unsigned int CurrentBand, CudaArray<Type>& Value, cudaStream_t stream)
	{
		if (Value.elements() != m_Rows || CurrentBand >= m_bands)
		{
			printf("维度不一致\n");throw;
			return;
		}
		cuSetCol(this->d_Data, Value.FirstAddr(), CurrentBand, m_Rows, m_Cols, m_Bands, Col, stream);
	}

	// 设置某一列
	void setBand(unsigned int OneBand, CudaArray<Type>& Value, cudaStream_t stream)
	{
		if (Value.elements() != m_Rows * m_Cols || OneBand >= m_Bands)
		{
			printf("维度不一致\n");throw;
			return;
		}
		cuSetBand(this->d_Data, Value.FirstAddr(), OneBand, m_Rows, m_Cols, m_Bands, stream);
	}

	// 设置流, 该流已在外部创建
	void setStream(cudaStream_t stream)
	{
		this->m_stream = stream;
	}



	void Reshape(unsigned int R, unsigned int C, unsigned int B)
	{
		if (R*C*B != this->m_Rows*this->m_Cols*this->m_Bands)
		{
			printf("reshape大小不一致");throw;
		}
		else
		{
			this->m_Rows = R;
			this->m_Cols = C;
			this->m_Bands = B;
		}
	}

private:
	unsigned int m_Rows;
	unsigned int m_Cols;
	unsigned int m_Bands;

	Type* d_Data; // 设备中的和路信号指针
	cudaError err;
	cudaStream_t m_stream; // 流编号

	// 分配内存
	cudaError Malloc(unsigned int Rows, unsigned int Cols, unsigned int Bands)
	{
		m_Rows = Rows;
		m_Cols = Cols;
		m_Bands = Bands;
		err = cudaMalloc((void**)&d_Data, sizeof(Type)*m_Rows * m_Cols * m_Bands);
		if (err != cudaSuccess)
		{
			m_Rows = 0;
			m_Cols = 0;
			m_Bands = 0;
			printf("GPU内存分配失败!错误为:%s\n", cudaGetErrorString(err));
		}
		return err;
	}

};

// 信号设置 主机中的矩阵
template <typename Type>
class SigHostMatrix
{
public:
	SigHostMatrix()
	{
		Rows = 0;
		Cols = 0;
		Bands = 0;
	}
	SigHostMatrix(unsigned int R, unsigned int C, unsigned int B)
	{
		Malloc(R, C, B);
	}
	SigHostMatrix(SigHostMatrix<Type> &h)
	{
		Rows = 0;
		Cols = 0;
		Bands = 0;
		CopyFromHost(h);
	}
	~SigHostMatrix() { }
	// 分配内存
	bool Malloc(unsigned int R, unsigned int C, unsigned int B)
	{
		Rows = R;
		Cols = C;
		Bands = B;
		err = true;
		h_Data = (Type*)malloc(sizeof(Type) * Rows * Cols * Bands);
		if (h_Data == nullptr)
		{
			err = false;
			Rows = 0;
			Cols = 0;
			Bands = 0;
		}
		return err;
	}

	// 释放内存
	void Free()
	{
		Rows = 0;
		Cols = 0;
		Bands = 0;
		if (h_Data != nullptr) free(h_Data);
		h_Data = nullptr;
	}

	// 获取首地址
	Type* FirstAddr() { return h_Data; };
	// 获取行数
	unsigned int Row() { return Rows; };
	// 获取列数
	unsigned int Col() { return Cols; };
	// 获取阵数
	unsigned int Band() { return Bands; };
	// 获取元素个数
	unsigned int elements() { return Rows * Cols * Bands; };

	unsigned int dims(int dims)
	{
		if (dims == 0) return m_Rows;
		else if (dims == 1) return m_Cols;
		else if (dims == 2) return m_Bands;
		else return -1;
	}

	// 从主机拷贝数据到主机中
	bool CopyFromHostPoint(Type* h_data, unsigned int R, unsigned int C, unsigned int B)
	{
		if (Rows != R || Cols != C || Bands != B)
		{
			free(d_Data);
			err = Malloc(R, C, B);
			if (err)
				err =memcpy(h_Data, h_data, R * C * B);
		}
		else
			err = memcpy(h_Data, h_data, R * C * B);

		return err;
	}
	// 从主机拷贝数据到主机中
	bool CopyFromHost(SigHostMatrix<Type> &h_data)
	{
		err = CopyFromHostPoint((Type*)h_data.FirstAddr(), h_data.Row(), h_data.Col(), h_data.Band());
	}

	// 从设备拷贝数据到主机中
	bool CopyFromDevice(Type* d_data, unsigned int R, unsigned int C, unsigned int B)
	{
		if (Rows != R || Cols != C || Bands != B)
		{
			free(d_data);
			err = Malloc(R, C, B);
			if (err)
				err = cudaMemcpy(h_Data, d_data, sizeof(Type) * R * C * B, cudaMemcpyDeviceToHost);
		}
		else
			err = cudaMemcpy(h_Data, d_data, sizeof(Type) * R * C * B, cudaMemcpyDeviceToHost);

		return err;
	}
	// 从设备拷贝数据到主机中
	bool CopyFromDevice(CudaArray<Type> &m)
	{
		err = CopyFromDevice(m.FirstAddr(), m.Row(), m.Col(), m.Band());
		return err;
	}

	SigHostMatrix<Type> operator=(SigHostMatrix<Type> &m)
	{
		CopyFromHostPoint(m.h_Data, m.Row(), m.Col(), m.Band());
		return *this;
	}

private:
	bool err;
	unsigned int Rows;
	unsigned int Cols;
	unsigned int Bands;
	Type* h_Data;

};

#endif // ! _CudaArray_H_
