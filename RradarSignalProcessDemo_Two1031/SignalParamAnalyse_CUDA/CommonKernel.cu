#include "CommonKernel.cuh"

__global__ void hello_from_gpu()
{

	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	printf("%d,%d.\n", bid, tid);

}

/*������ʱ����������*/
__global__ void TBaseGen(double *t, double fs, double t0, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		t[i] = t0 + i / fs;
	}
}

/*�������Ӻ˺���
������*c = *a + *b
*/
__global__ void ComplexAddKernel(const cuComplex *a, const cuComplex *b, cuComplex *c, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		c[i].x = a[i].x + b[i].x;
		c[i].y = a[i].y + b[i].y;
	}
}


/*���������˺���
������*b = *a
*/
__global__ void ComplexCopyKernel(const cuComplex *a, cuComplex *b, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		b[i].x = a[i].x;
		b[i].y = a[i].y;
	}
}

/*dB
������*b = 20*log10(*a)
*/
__global__ void dBKernel(const float *a, float* b, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		b[i] = 20 * log10(a[i]);
	}
}

/*idB
������*b = pow(10, *a / 20)
*/
__global__ void idBKernel(const float *a, float* b, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		b[i] = pow(10.0, a[i] / 20.0);
	}
}

/*��ʵ������תΪ�鲿Ϊ0�ĸ�������
������*a = *(a+0i)
*/
__global__ void ComplexKernel(const float *a, cuComplex* b, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		b[i].x = a[i];
		b[i].y = 0.0;
	}

}


/*����
������*a = *(a+0i)
*/
__global__ void ConjKernel(cuComplex *data, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		data[i].x = data[i].x;
		data[i].y = -data[i].y;
	}

}
/*������ʵ��ƴ��һ������
������*Res = *(Real+Imagi)
*/
__global__ void ComplexMat(cuComplex *Res, float *Real, float *Imag, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i].x = Real[i];
		Res[i].y = Imag[i];
	}
}


//����ת��
__global__ void TransposeKernel(cuComplex* in, cuComplex* out, int Rows, int Cols, int Bands)
{
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		out[B_bias * (Rows * Cols) + R_bias * Cols + C_bias] = in[idx];

	}
}
__global__ void TransposeKernel(float* in, float* out, int Rows, int Cols, int Bands)
{
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		out[B_bias * (Rows * Cols) + R_bias * Cols + C_bias] = in[idx];

	}
}
__global__ void TransposeKernel(int* in, int* out, int Rows, int Cols, int Bands)
{
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		out[B_bias * (Rows * Cols) + R_bias * Cols + C_bias] = in[idx];

	}
}


//����ת��
__global__ void CTransposeKernel(cuComplex* in, cuComplex* out, int Rows, int Cols, int Bands)
{
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		out[B_bias * (Rows * Cols) + R_bias * Cols + C_bias].x = in[idx].x;
		out[B_bias * (Rows * Cols) + R_bias * Cols + C_bias].y = -in[idx].y;

	}
}


/*ȡʵ��
������*Real = real(*data��
*/
__global__ void Realkernel(cuComplex *Data, float *Real, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Real[i] = Data[i].x;
	}
}

/*ȡ�鲿
������*Imag = imag(*data��
*/
__global__ void Imagkernel(cuComplex *Data, float *Imag, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Imag[i] = Data[i].y;
	}
}

/*��CCMat��abs
*/
__global__ void Abskernel(cuComplex *DataIn, float *DataOut, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		DataOut[i] = sqrt(DataIn[i].x * DataIn[i].x + DataIn[i].y * DataIn[i].y);
	}
}


// ���Ե�ͨ�ź���
__global__ void LowPass(cuComplex *Res, int StartPoint, int EndPoint, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		if (i >= StartPoint && i < EndPoint)
		{
			Res[i].x = 0;
			Res[i].y = 0;
		}
	}
}

/*��������������
������*Res = *Data1 .* *Data2
*/
__global__ void DotMulKernal(cuComplex *Data1, cuComplex *Data2, cuComplex *Res, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i] = cuCmulf(Data1[i], Data2[i]);
	}
}

__global__ void DotMul2Kernal(cuComplex *Data1, cuComplex *Data2, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Data1[i] = cuCmulf(Data1[i], Data2[i]);
	}
}

/*���������������
������*Res = *Data1 + *Data2
*/
__global__ void MatAddKernal(cuComplex *Data1, cuComplex *Data2, cuComplex *Res, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i].x = Data1[i].x + Data2[i].x;
		Res[i].y = Data1[i].y + Data2[i].y;
	}
}

/*���������������
������*Res = *Data1 + *Data2 + *Data3
*/
__global__ void MatAddKernal(cuComplex *Data1, cuComplex *Data2, cuComplex* Data3, cuComplex *Res, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i].x = Data1[i].x + Data2[i].x + Data3[i].x;
		Res[i].y = Data1[i].y + Data2[i].y + Data3[i].y;
	}
}

__global__ void MatAddComplexKernal(cuComplex *Data1, cuComplex Data2, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Data1[i].x = Data1[i].x + Data2.x;
		Data1[i].y = Data1[i].y + Data2.y;
	}
}

/*����������ʵ������
������*Data1.x = *Data1.x * *Data2
*Data1.y = *Data1.y * *Data2
*/
__global__ void MatMulKernal(cuComplex *Data1, float *Data2, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Data1[i].x = Data1[i].x * Data2[i];
		Data1[i].y = Data1[i].y * Data2[i];
	}
}


/*���������float
������*Res = *Data1 / Data2
*/
__global__ void C2FDivKernal(cuComplex *Data1, float Data2, cuComplex *Res, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i].x = Data1[i].x / Data2;
		Res[i].y = Data1[i].y / Data2;
	}
}

__global__ void C2FDiv2Kernal(cuComplex *Data1, float Data2,int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Data1[i].x = Data1[i].x / Data2;
		Data1[i].y = Data1[i].y / Data2;
	}
}

/*���������float
������*Res = *Data1 * Data2
*/
__global__ void C2FMulKernal(cuComplex *Data1, float Data2, cuComplex *Res, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i].x = Data1[i].x * Data2;
		Res[i].y = Data1[i].y * Data2;
	}
}

/*float�����float
������*Res = *Data1 * Data2
*/
__global__ void F2FMulKernal(float *Data1, float Data2, float *Res, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Res[i] = Data1[i] * Data2;
	}
}

/*��ȡ���
��������Data
*/
__global__ void SliceKernel(cuComplex *Data, cuComplex *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		if ((R_bias >= RStart && R_bias <= REnd) &&
			(C_bias >= CStart && C_bias <= CEnd) &&
			(B_bias >= BStart && B_bias <= BEnd))
		{
			unsigned int ResultRowNum = (REnd - RStart + 1);
			unsigned int ResultColNum = (CEnd - CStart + 1);
			Result[(B_bias - BStart) * ResultRowNum * ResultColNum + (C_bias - CStart) * ResultRowNum + (R_bias - RStart)] = Data[idx];

		}

	}


}
__global__ void SliceKernel(float *Data, float *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		if ((R_bias >= RStart && R_bias <= REnd) &&
			(C_bias >= CStart && C_bias <= CEnd) &&
			(B_bias >= BStart && B_bias <= BEnd))
		{
			unsigned int ResultRowNum = (REnd - RStart + 1);
			unsigned int ResultColNum = (CEnd - CStart + 1);
			Result[(B_bias - BStart) * ResultRowNum * ResultColNum + (C_bias - CStart) * ResultRowNum + (R_bias - RStart)] = Data[idx];

		}

	}


}
__global__ void SliceKernel(int *Data, int *Result, int RStart, int REnd, int CStart, int CEnd, int BStart, int BEnd, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows*Cols*Bands)
	{
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		if ((R_bias >= RStart && R_bias <= REnd) &&
			(C_bias >= CStart && C_bias <= CEnd) &&
			(B_bias >= BStart && B_bias <= BEnd))
		{
			unsigned int ResultRowNum = (REnd - RStart + 1);
			unsigned int ResultColNum = (CEnd - CStart + 1);
			Result[(B_bias - BStart) * ResultRowNum * ResultColNum + (C_bias - CStart) * ResultRowNum + (R_bias - RStart)] = Data[idx];

		}

	}


}

/*��������չ
��������Data��չΪResult���ȣ�������չ��
Data : (Rows, Cols, Bands)
Result : (Rows + abs(AddLength), Cols, Bands)
AddLength : ��չ�������Ҹ���
Value : ��չ��ֵ
*/
__global__ void ExtendKernel(float* Data, float* Result, int AddLength, float Value, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int ResultRows = Rows + abs(AddLength);
	if (idx < ResultRows * Cols * Bands)
	{
		/*resultƫ��*/
		// ������
		unsigned int B_bias = idx / (ResultRows * Cols);
		// ������
		unsigned int C_bias = (idx % (ResultRows * Cols)) / ResultRows;
		// ������
		unsigned int R_bias = (idx % (ResultRows * Cols)) % ResultRows;

		// ����չ
		if (AddLength > 0)
		{
			// ��չ����
			if (R_bias >= Rows)
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Value;
			}
			// ����չ����
			else
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Data[B_bias * Rows * Cols + C_bias * Rows + R_bias];
			}
		}
		// ����չ
		else
		{
			// ��չ����
			if (R_bias < -AddLength)
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Value;
			}
			// ����չ����
			else
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Data[B_bias * Rows * Cols + C_bias * Rows + R_bias + AddLength];
			}
		}
	}
}


/*��������չ
��������Data��չΪResult���ȣ�������չ��
Data : (Rows, Cols, Bands)
Result : (Rows + abs(AddLength), Cols, Bands)
AddLength : ��չ�������Ҹ���
Value : ��չ��ֵ
*/
__global__ void ExtendKernel(cuComplex* Data, cuComplex* Result, int AddLength, cuComplex Value, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	int ResultRows = Rows + abs(AddLength);
	if (idx < ResultRows * Cols * Bands)
	{
		/*resultƫ��*/
		// ������
		unsigned int B_bias = idx / (ResultRows * Cols);
		// ������
		unsigned int C_bias = (idx % (ResultRows * Cols)) / ResultRows;
		// ������
		unsigned int R_bias = (idx % (ResultRows * Cols)) % ResultRows;

		// ����չ
		if (AddLength > 0)
		{
			// ��չ����
			if (R_bias >= Rows)
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Value;
			}
			// ����չ����
			else
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Data[B_bias * Rows * Cols + C_bias * Rows + R_bias];
			}
		}
		// ����չ
		else
		{
			// ��չ����
			if (R_bias < -AddLength)
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Value;
			}
			// ����չ����
			else
			{
				Result[B_bias * ResultRows * Cols + C_bias * ResultRows + R_bias] = Data[B_bias * Rows * Cols + C_bias * Rows + R_bias + AddLength];
			}
		}
	}
}

/*��ֵɸѡ
������
Data		: (Rows, Cols, Bands)
Result		: (Rows + abs(AddLength), Cols, Bands)
Th			: ��ֵ
CompareFlag	: �Ƚ�����{'<','>','<=','>=','=='}
����
Th = 3; CompareFlag = '>'
Data	:  4 5 6 1 8 3 2 4
Result	:  1 1 1 0 1 0 0 1
*/
__global__ void CompareKernel(float* Data, bool* Result, float Th, char CompareFlag, unsigned int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		if (CompareFlag == '<')
		{
			if (Data[idx] < Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '>')
		{
			if (Data[idx] > Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '<=')
		{
			if (Data[idx] <= Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '>=')
		{
			if (Data[idx] >= Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '==')
		{
			if (Data[idx] == Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
	}
}
__global__ void CompareKernel(cuComplex* Data, bool* Result, float Th, char CompareFlag, unsigned int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		if (CompareFlag == '<')
		{
			if (cuCabsf(Data[idx]) < Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '>')
		{
			if (cuCabsf(Data[idx]) > Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '<=')
		{
			if (cuCabsf(Data[idx]) <= Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '>=')
		{
			if (cuCabsf(Data[idx]) >= Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
		else if (CompareFlag == '==')
		{
			if (cuCabsf(Data[idx]) == Th)
			{
				Result[idx] = true;
			}
			else
			{
				Result[idx] = false;
			}
		}
	}
}
__global__ void CompareKernel(float* Data, int* Result, float Th, char CompareFlag, unsigned int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		if (CompareFlag == '<')
		{
			if (Data[idx] < Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '>')
		{
			if (Data[idx] > Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '<=')
		{
			if (Data[idx] <= Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '>=')
		{
			if (Data[idx] >= Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '==')
		{
			if (Data[idx] == Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
	}
}
__global__ void CompareKernel(cuComplex* Data, int* Result, float Th, char CompareFlag, unsigned int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		if (CompareFlag == '<')
		{
			if (cuCabsf(Data[idx]) < Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '>')
		{
			if (cuCabsf(Data[idx]) > Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '<=')
		{
			if (cuCabsf(Data[idx]) <= Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '>=')
		{
			if (cuCabsf(Data[idx]) >= Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
		else if (CompareFlag == '==')
		{
			if (cuCabsf(Data[idx]) == Th)
			{
				Result[idx] = 1;
			}
			else
			{
				Result[idx] = 0;
			}
		}
	}
}

/*
���
Data: [R C B]
Result: [1 C B]
*/
__global__ void cuSumKernel(float *Data, float* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows * Cols * Bands)
	{
		/*resultƫ��*/
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		atomicAdd(&Result[B_bias * Cols + C_bias], Data[idx]);
	}
}
__global__ void cuSumKernel(cuComplex* Data, cuComplex* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows * Cols * Bands)
	{
		/*resultƫ��*/
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		atomicAdd(&Result[B_bias * Cols + C_bias].x, Data[idx].x);
		atomicAdd(&Result[B_bias * Cols + C_bias].y, Data[idx].y);
	}

}
__global__ void cuSumKernel(int *Data, int* Result, unsigned int Rows, unsigned int Cols, unsigned int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows * Cols * Bands)
	{
		/*resultƫ��*/
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		atomicAdd(&Result[B_bias * Cols + C_bias], Data[idx]);
	}
}

__global__ void cuComplexPowATileMean(cuComplex *data, cuComplex *atilemean, cuComplex *sum_mean, float *powatilemean, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		float mean1 = sum_mean[0].x / N;
		float mean2 = sum_mean[0].y / N;
		float dataI = data[idx].x - mean1;//�ظ�
		float dataQ = data[idx].y - mean2;//�ظ�
		atilemean[idx].x = dataI;//�ظ�
		atilemean[idx].y = dataQ;//�ظ�
		powatilemean[idx] = dataI * dataI + dataQ * dataQ;
	}
}

__global__ void SqrtTileComplexDiv(float *sum, cuComplex *alitemean, cuComplex *result, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		float temp = sum[0] / (N - 1);
		float sigmma = sqrt(temp);
		result[idx].x = alitemean[idx].x / sigmma;
		result[idx].y = alitemean[idx].y / sigmma;
	}
}
__global__ void SqrtTileDiv(float *sum, float *alitemean, float *result, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		float temp = sum[0] / (N - 1);
		float sigmma = sqrt(temp);
		result[idx] = alitemean[idx] / sigmma;
	}
}
__global__ void cuPowATileMean(float *data, float *atilemean, float *sum_mean, float *powatilemean, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		float mean = sum_mean[0] / N;
		float tmp = data[idx] - mean;//�ظ�
		atilemean[idx] = tmp;
		powatilemean[idx] = tmp * tmp;
	}
}

__global__ void z_score_kernel(cuComplex *data, cuComplex *atiledata, float* PowerResult, unsigned int Rows, unsigned int Cols, unsigned int Bands, cuComplex *Sum_first, float *Sum_second, cuComplex* z_score)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Rows * Cols * Bands)
	{
		unsigned int N = Rows*Cols*Bands;
		/*Sum_firstƫ��*/
		// ������
		unsigned int B_bias = idx / (Rows * Cols);
		// ������
		unsigned int C_bias = (idx % (Rows * Cols)) / Rows;
		// ������
		unsigned int R_bias = (idx % (Rows * Cols)) % Rows;

		//ԭ�Ӳ������
		atomicAdd(&Sum_first[B_bias * Cols + C_bias].x, data[idx].x);
		atomicAdd(&Sum_first[B_bias * Cols + C_bias].y, data[idx].y);

		
		float mean1 = Sum_first[0].x / N;
		float mean2 = Sum_first[0].y / N;
		atiledata[idx].x = data[idx].x - mean1;//�ظ�
		atiledata[idx].y = data[idx].y - mean2;
		PowerResult[idx] = (atiledata[idx].x) *(atiledata[idx].x) + (atiledata[idx].y)*(atiledata[idx].y);
		atomicAdd(&Sum_second[B_bias * Cols + C_bias], PowerResult[idx]);

		float temp = Sum_second[0] / (N - 1);
		float sigmma = sqrt(temp);
		z_score[idx].x = (atiledata[idx].x) / sigmma;
		z_score[idx].y = (atiledata[idx].y) / sigmma;
	}
}

// ���н���3άѭ����λ
__global__ void CirculShiftCol(float* Data, float* Result, int shift_num, int Rows, int Cols, int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	unsigned int RCMul = Rows * Cols;
	if (idx < RCMul * Bands)
	{
		// ��ά����
		unsigned int indx = (idx % RCMul) % Rows;
		unsigned int indy = (idx % RCMul) / Rows;
		unsigned int indz = idx / RCMul;

		int indx_shift = (indx + Rows + shift_num) % Rows;

		int index_new = indz * RCMul + indy * Rows + indx_shift;

		Result[idx] = Data[index_new];
	}
}
__global__ void CirculShiftCol(cuComplex* Data, cuComplex* Result, int shift_num, int Rows, int Cols, int Bands)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	unsigned int RCMul = Rows * Cols;
	if (idx < RCMul * Bands)
	{
		// ��ά����
		unsigned int indx = (idx % RCMul) % Rows;
		unsigned int indy = (idx % RCMul) / Rows;
		unsigned int indz = idx / RCMul;

		int indx_shift = (indx + Rows + shift_num) % Rows;

		int index_new = indz * RCMul + indy * Rows + indx_shift;

		Result[idx] = Data[index_new];
	}
}

// ����2λѭ����λ
__global__ void Circcushift(cuComplex *src, cuComplex *tar, int  Row, int Col, int shift_num)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int shift = (idx + shift_num * Row) % (Row * Col);

	if (idx < Row * Col)
	{
		tar[idx] = src[shift];
	}
}
__global__ void Circcushift(float *src, float *tar, int  Row, int Col, int shift_num)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int shift = (idx + shift_num * Row) % (Row * Col);

	if (idx < Row * Col)
	{
		tar[idx] = src[shift];
	}
}




/*�ź�����ר��*/
__global__ void SetPerTarHn(cuComplex *Data, cuComplex *Temp, int colIndex, int elements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < elements)
	{
		Data[colIndex * elements + i] = Temp[i];
	}
}