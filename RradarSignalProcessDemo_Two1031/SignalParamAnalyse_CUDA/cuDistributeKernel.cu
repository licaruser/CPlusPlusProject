
#include"cuDistributeKernel.cuh"

#define PI	double(3.141592653589793)
// ���ɸ�˹�ֲ��Ӳ�����
__global__ void GaussKernel(float sigmaf, float *hf, int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		float f = -6e3 + idx * 12e3 / elements;

		hf[idx] = exp( - f * f / 2 / sigmaf / sigmaf);
	}
}
// ����ָ���ֲ��Ӳ�����
__global__ void ExponentialKernel(float sigmaf, float *hf, int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		float f = -6e3 + idx * 12e3 / elements;

		hf[idx] = exp(-abs(f) / sigmaf);
	}
}

__global__ void CauchyKernel(float sigmaf, float *hf, int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		float f = -6e3 + idx * 12e3 / elements;

		float tmp = abs(f / sigmaf);
		hf[idx] = 1. / (1.0 + pow(tmp, (float)2.0));
	}
}

__global__ void FullSpectrumKernel(float sigmaf, float *hf, int SpectrumPara, int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		float f = -6e3 + idx * 12e3 / elements;

		hf[idx] = 1. / (1.0 + pow(abs(f / sigmaf), (float)SpectrumPara));
	}
}
// ˹ά����� 1��2 ��
__global__ void Swerlling12Kernel(float *Uniformx, float* y, double sigmac, int elements)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < elements)
	{
		y[idx] = -sigmac * log(1.0 - Uniformx[idx]);
	}

}
__global__ void LinearTrans(float *result, float *data, float scale, float bias, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		result[idx] = scale * data[idx] + bias;
	}
}
//��pow
__global__ void cuPow(float *data, float *result, int N, float base)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		result[idx] = std::pow(data[idx], base);
	}
}
//ʵ��+�鲿=����
__global__ void floatToComplex(float *real, float *imag, cuComplex *tar, int len)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < len)
	{
		tar[idx].x = real[idx];
		tar[idx].y = imag[idx];
	}
}
//����=ʵ��+�鲿
__global__ void ComplexTofloat(float *real, float *imag, cuComplex *tar, int len)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < len)
	{
		real[idx] = tar[idx].x;
		imag[idx] = tar[idx].y;
	}
}
// ����ֵ
__global__ void cuAbs(cuComplex *Data, float *Result, int Len)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < Len)
	{
		Result[idx] = cuCabsf(Data[idx]);

	}
}
//������ָ��
__global__ void ComplexIndex(cuComplex *data, int N, float base)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < N)
	{
		//��Ƕ�
		float eps = 10e-10;
		float AngleResult = 0.0;
		if (abs(data[idx].x) <= eps)
			AngleResult = 0.0;
		else if (data[idx].x > 0)
		{
			float angleTmp = data[idx].y / data[idx].x;
			AngleResult = atan(angleTmp);
		}
		else if (data[idx].x <= 0 && data[idx].y >0)
		{
			float angleTmp = data[idx].y / data[idx].x;
			AngleResult = PI + atan(angleTmp);
		}
		else
		{
			float angleTmp = data[idx].y / data[idx].x;
			AngleResult = atan(angleTmp) - PI;
		}
		float Amp = pow(cuCabsf(data[idx]), base);

		//������ָ���������˷�ֵ�����Ǻ�ָ��
		data[idx].x = Amp*cosf(AngleResult * base);
		data[idx].y = Amp*sinf(AngleResult * base);
	}
}

//__global__ void SqrtTileComplexDiv(float *sum, cuComplex *alitemean, cuComplex *result, int N)
//{
//	int idx = blockDim.x * blockIdx.x + threadIdx.x;
//	if (idx < N)
//	{
//		float temp = sum[0] / (N - 1);
//		float sigmma = sqrt(temp);
//		result[idx].x = alitemean[idx].x / sigmma;
//		result[idx].y = alitemean[idx].y / sigmma;
//	}
//}
//__global__ void cuPowATileMean(float *data, float *atilemean, float *sum_mean, float *powatilemean, int N)
//{
//	int idx = blockDim.x * blockIdx.x + threadIdx.x;
//	if (idx < N)
//	{
//		float mean = sum_mean[0] / N;
//		atilemean[idx] = data[idx] - mean;//�ظ�
//		powatilemean[idx] = atilemean[idx] * atilemean[idx];
//	}
//}
//__global__ void cuComplexPowATileMean(cuComplex *data, cuComplex *atilemean, cuComplex *sum_mean, float *powatilemean, int N)
//{
//	int idx = blockDim.x * blockIdx.x + threadIdx.x;
//	if (idx < N)
//	{
//		float mean1 = sum_mean[0].x / N;
//		float mean2 = sum_mean[0].y / N;
//		atilemean[idx].x = data[idx].x - mean1;//�ظ�
//		atilemean[idx].y = data[idx].y - mean2;//�ظ�
//		powatilemean[idx] = atilemean[idx].x * atilemean[idx].x + atilemean[idx].y * atilemean[idx].y;
//	}
//}

/*  float / ��float���� + float�� */
__global__ void MatDivKernel(float *data, float *result, float scale, float bias, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		result[idx] = scale / (data[idx] + bias);
	}
}

//����float���У�start:step:start+(N-1)*step
__global__ void seqKernel(float *result, float start, float step, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		result[idx] = start + step * (float)idx;
	}
}

//�����������Ƶ�һ��������
__global__ void vec2mat(float * data, float *result, int rows, int cols)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int N = rows * cols;
	int bias = idx % rows;
	if (idx < N)
	{
		result[idx] = data[bias];
	}
}

__global__ void Complexvec2mat(cuComplex * data, cuComplex *result, int rows, int cols)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int N = rows * cols;
	int bias = idx % rows;
	if (idx < N)
	{
		result[idx].x = data[bias].x;
		result[idx].y = data[bias].y;
	}
}

//�ɽǶȸ���ŷ����ʽ�õ����ź�
__global__ void EulerFormula(float * sata, cuComplex *result, int N)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < N)
	{
		result[idx].x = cos(sata[idx]);
		result[idx].y = sin(sata[idx]);
	}

}

//����ת��
__global__ void transposeKernel(float* in, float* out, int Rows, int Cols)
{
	long idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Rows*Cols)
	{
		int R = int(idx / Cols);
		int C = idx % Cols;
		out[C*Rows + R] = in[idx];

	}
}

