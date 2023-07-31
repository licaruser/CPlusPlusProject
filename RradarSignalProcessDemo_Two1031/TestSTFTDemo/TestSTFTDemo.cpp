#include "TestSTFTDemo.h"



ParamAnalyse* m_TmpParamEsta = new(ParamAnalyse);       //参数估计


int main()
{
	LARGE_INTEGER Code_start;
	LARGE_INTEGER Code_stop;
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);

	// 读入信号
	STFT_lgy Tmp_stft;
	vector<complex<double>> SignalData;
	Tmp_stft.Read_Csv_Data(SignalData, 1, 1);


	vector<vector<complex<double>>> trf;
	double f;
	vector<int> t;

	for (int ii = 1; ii < SignalData.size() + 1; ii++)
	{
		t.push_back(ii);
	}
	int N = 256;
	//代码开始，计时开始
	QueryPerformanceCounter(&Code_start);
	Tmp_stft.tfrstft(trf, f, SignalData, t, N);
	//代码结束，计时结束
	QueryPerformanceCounter(&Code_stop);
	double time_sec = (unsigned long long)(Code_stop.QuadPart - Code_start.QuadPart) / (double)freq.QuadPart;
	cout << "stft程序运行时间："<< time_sec << endl;

	//ParamAnalyse TestClass;     //直接使用类定义申明
	//TestClass.StepAdvance(trf);
	m_TmpParamEsta->StepAdvance(trf, SignalData);   //使用new申请对象，在堆中分配了内存，堆上的内存分配，亦称动态内存分配
	
	
	//录取短时傅里叶变换后的数据
	//Tmp_stft.SaveData_Csv(trf);


	//int signal_len = SignalData.size();
	//fftw_complex* input, * output;
	//fftw_plan plan;
	//int N = 256; // 窗口大小
	//int hop = N / 2; // 滑动步长
	//double* x = new double[N];
	//double* win = new double[N]; // 窗口函数
	//// 初始化窗口函数
	//for (int i = 0; i < N; i++) {
	//	win[i] = 0.54 - 0.46 * cos(2.0 * M_PI * i / N);
	//}
	//// 实现STFT
	//int L = signal_len / hop; // 帧数
	//input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	//output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	//plan = fftw_plan_dft_1d(N, input, output, FFTW_FORWARD, FFTW_ESTIMATE);
	//for (int i = 0; i < L; i++) 
	//{
	//	for (int j = 0; j < N; j++) 
	//	{
	//		x[j] = x_signal_buf[i * hop + j] * win[j];
	//		input[j][0] = x[j];
	//		input[j][1] = 0.0;
	//	}
	//	fftw_execute(plan);
	//	// 处理FFT结果
	//	// ...
	//}
	//fftw_free(input);
	//fftw_free(output);
	//fftw_destroy_plan(plan);
	return 0;
}


STFT_lgy::STFT_lgy()
	//:m_TmpParamEsta(new ParamAnalyse)
{

}

STFT_lgy::~STFT_lgy()
{
}


void STFT_lgy::Read_Excel_Data(int init_num)
{
	Book* book = xlCreateXMLBook();//创建一个Excel Book
	if (book)
	{
		if (book->load("")) // 打开xlsv的
		{
			//加载excel文件
			Sheet* sheet = book->getSheet(0);//获取第一个工作表
			if (sheet)
			{
				int row = sheet->lastRow(), col = sheet->lastCol();
				for (int ii = init_num; ii < row; ++ii)  // 行
				{
					for (int jj = 0; jj < col; ++jj) // 列
					{
						const char* text = sheet->readStr(ii, jj);
						if (text)
						{
							std::cout << text << "\t";
						}
					}
					std::cout << std::endl;
				}
			}
		}
		book->release(); //释放Excel Book
	}

}

void STFT_lgy::Read_Csv_Data(vector<complex<double>>& SignalData, int Init_row, int Init_col)
{
	SignalData.clear();
	int start_row = Init_row;
	int start_col = Init_col;

	std::ifstream file("..\\Data\\4.3GHz_400MHz_Pw1us_PRI10us_000.csv");
	std::vector<vector<std::string>> data;
	std::string line;

	int row_index = 0;
	while (getline(file, line))
	{
		// 跳过起始行之前的所有行
		if (row_index < start_row - 1)
		{
			row_index++;
			continue;
		}
		// 解析当前行并跳过起始列之前的所有列
		std::vector<std::string> row;
		std::stringstream ss(line);
		int col_index = 0;
		std::string cell;
		while (getline(ss, cell, ','))
		{
			if (col_index < start_col - 1)
			{
				col_index++;
				continue;
			}
			row.push_back(cell);
			col_index++;
		}
		data.push_back(row);
		row_index++;
	}
	//vector<complex<double>>  RadarSignalData;
	complex<double> Tmp_complex_cell;

	for (const auto& row : data) {
		//Tmp_complex_cell.real(std::atof(row.at(0).c_str())); //时间
		//Tmp_complex_cell.imag(std::atof(row.at(1).c_str())); //信号

		//读取只有实部的信号
		Tmp_complex_cell.real(std::atof(row.at(1).c_str())); //时间
		Tmp_complex_cell.imag(0.0); //信号
		SignalData.push_back(Tmp_complex_cell);
	}
	
}

//template<typename T>
//T min(T a, T b)
//{
//	return(a < b) ? a : b;
//}
//template<typename T, typename Args>
//T min(T a, T b, Args c)
//{
//	return min(a < b ? a : b, c);
//}



void STFT_lgy::tfrstft(vector<vector<complex<double>>>& trf, double& f, const vector<complex<double>>& x, vector<int>& t, int& N)
{

	vector<vector<complex<double>>> Tmp_trf;
	vector<complex<double>> Tmp_SourceData;
	vector<int> Time;
	int FsPoint;
	Tmp_SourceData = x;
	Time = t;
	FsPoint = N;

	int xrow = Tmp_SourceData.size();
	int hlength = floor(FsPoint / 4);
	hlength = hlength + 1 - (hlength % 2);//取余数

	Tmp_trf.resize(xrow, std::vector<complex<double>>(FsPoint,0));

	// 创建窗口函数
	std::vector<double> window;
	for (int ii = 1; ii < hlength + 1; ii++)
	{
		window.push_back(0.54 - 0.46 * (std::cos(2.0 * M_PI * ii / (hlength + 1))));
	}
	int hrow = window.size();
	int Lh = (hrow - 1) / 2;

	VectorXd v(hrow);
	for (int ii = 0; ii < window.size(); ii++)
	{
		v[ii] =  window.at(ii);
	}
	double n = v.norm();
	for (int ii = 0; ii < window.size(); ii++)
	{
		window.at(ii) = window.at(ii) / n;
	}

	int tcol = Time.size();

	//短时傅里叶变换的核心
	for (int icol = 0; icol < tcol; icol++)
	{
		int ti = Time.at(icol);
		int tau_FuMin = -min_member(round(FsPoint / 2), Lh, ti - 1);
		int tau_ZenMin = min_member(round(FsPoint / 2) - 1, Lh, xrow - ti);
		vector<int> tau;

		//保护tau_FuMin = tau_ZenMin+1 = 0；
		if ((tau_ZenMin + 1) == 0)
		{
			if (tau_FuMin == 0)
			{
				tau.push_back(0);
			}
		}

		for (int ii_tau = tau_FuMin; ii_tau < tau_ZenMin + 1; ii_tau++)
		{
			tau.push_back(ii_tau);
		}

		vector<int> indices;
		complex<double> Tempcomplex;
		vector<complex<double>> vectorTempComplex;
		for (int indi = 0; indi < tau.size(); indi++)
		{
			indices.push_back(((tau.at(indi) + FsPoint) % FsPoint));

			double window_real = window.at(Lh + tau.at(indi));
			//stft核心计算
			Tempcomplex.real(Tmp_SourceData.at(ti + tau.at(indi) - 1).real() * window_real);
			Tempcomplex.imag(Tmp_SourceData.at(ti + tau.at(indi) - 1).imag() * window_real);

			vectorTempComplex.push_back(Tempcomplex);
			//Tmp_trf.at(icol).at(((tau.at(indi) + FsPoint) % FsPoint) + 1).real(Tmp_SourceData.at(ti + tau.at(indi)).real() * window.at(Lh + 1 + tau.at(indi)));
			//Tmp_trf.at(icol).at(((tau.at(indi) + FsPoint) % FsPoint) + 1).imag(Tmp_SourceData.at(ti + tau.at(indi)).imag() * window.at(Lh + 1 + tau.at(indi)));
		}

		for (int pos = 0; pos < indices.size(); pos++)
		{
			int pos_indices = indices.at(pos);
			Tmp_trf.at(icol).at(pos_indices).real(vectorTempComplex.at(pos).real());
			Tmp_trf.at(icol).at(pos_indices).imag(vectorTempComplex.at(pos).imag());
		}
		//int oo = 1;
	}


	//二维FFt失败
	//int N1 = xrow;   //62500
	//int N2 = FsPoint;  // 64
	//fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N1 * N2);
	//fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N1 * N2);
	////初始化输入数据
	////for (int fft_ii = 0; fft_ii < N1 * N2; fft_ii++)
	////{
	//int fft_ii = 0;
	////for (const auto& col_ii : Tmp_trf)
	////{
	////	for (const auto& row_ii : col_ii)
	////	{
	////		in[fft_ii][0] = row_ii.real();
	////		in[fft_ii][1] = row_ii.imag();
	////		fft_ii++;
	////	}
	////}
	//// 先放64行，再放62500
	//for (int ii_row = 0; ii_row < FsPoint; ii_row++)  //64
	//{
	//	for (int ii_col = 0; ii_col < xrow; ii_col++)  //62500 
	//	{
	//		in[fft_ii][0] = Tmp_trf.at(ii_col).at(ii_row).real();
	//		in[fft_ii][1] = Tmp_trf.at(ii_col).at(ii_row).imag();
	//		fft_ii++;
	//	}
	//}
	////}
	//// 创建 FFT 变换的计划对象
	//fftw_plan  p = fftw_plan_dft_2d(N1, N2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//// 执行 FFT 变换
	//fftw_execute(p);
	//// 输出结果
	//for (int out_ii = 0; out_ii < N1; out_ii++)  // 62500
	//{
	//	for (int out_jj = 0; out_jj < N2; out_jj++)  // 64
	//	{
	//		Tmp_trf.at(out_ii).at(out_jj).real(out[out_ii * N2 + out_jj][0]);//实部
	//		Tmp_trf.at(out_ii).at(out_jj).imag(out[out_ii * N2 + out_jj][1]);//虚部
	//	}
	//}


	vector<complex<double>> Tmp_in;
	// 一维FFT
	int NN = FsPoint;  // 64
	fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * NN);
	fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * NN);
	fftw_plan  plan = fftw_plan_dft_1d(NN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for (int M_col = 0; M_col < Tmp_trf.size(); M_col++)
	{
		Tmp_in = Tmp_trf.at(M_col);

		// 创建 FFT 变换的计划对象

		// 存储输入数据到数组 in
		for (int ii = 0; ii < NN; ii++)
		{
			in[ii][0] = Tmp_in.at(ii).real();
			in[ii][1] = Tmp_in.at(ii).imag();
		}
		// 执行DFT计算
		fftw_execute(plan);

		//// fftshift
		//int kk = NN / 2;
		//for (int i = 0; i < kk; i++)
		//{
		//	double tmp1 = out[i][0];
		//	double tmp2 = out[i][1];
		//	out[i][0] = out[i + kk][0];
		//	out[i][1] = out[i + kk][1];
		//	out[i + kk][0] = tmp1;
		//	out[i + kk][1] = tmp2;
		//}

		// 输出DFT结果
		for (int ii = 0; ii < NN; ii++)
		{
			Tmp_trf.at(M_col).at(ii).real(out[ii][0]);
			Tmp_trf.at(M_col).at(ii).imag(out[ii][1]);
		}

	}
	// 释放内存和销毁声明 
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan);
	



	trf = Tmp_trf;

}

int STFT_lgy::min_member(int a, int b, int c)
{
	int Temp = (a < b) ? a : b;
	return (Temp < c) ? Temp : c;
}

void STFT_lgy::SaveData_Csv(const vector<vector<complex<double>>>& Data)
{
	//vector<complex<double>> Tmp_data;
	//for (int ii = 0; ii < Data.size(); ii++)
	//{
	//	string path = "Data\\HonstCompleData.csv";
	//	Tmp_data = Data.at(ii);
	//	std::fstream fp;
	//	fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
	//	for (int jj = 0; jj < Tmp_data.size(); jj++)
	//	{
	//		double RData = Tmp_data.at(jj).real();
	//		double IData = Tmp_data.at(jj).imag();
	//		// 写入实部
	//		fp << RData;
	//		// 写入符号和虚部
	//		if (IData != -0 || IData != 0)
	//		{
	//			fp << (IData > 0 ? "+" : "") << IData << "i";
	//		}
	//		fp << endl;
	//	}
	//	fp.close();
	//	std::cout << path << " 保存成功! " << endl;
	//	Tmp_data.clear();
	//}

	string path = "..\\Data\\HonstCompleData.csv";
	std::fstream fp;
	fp.open(path, std::ios::in | std::ios::out | std::ios::trunc);
	int row = Data.at(0).size();  //行 64
	int col = Data.size();        //列 2500

	for (int jj = 0; jj < col; jj++)  // 2500
	{
		for (int ii = 0; ii < row; ii++)  // 64
		{
			double RData = Data.at(jj).at(ii).real();
			double IData = Data.at(jj).at(ii).imag();
			// 写入实部
			fp << RData;
			// 写入符号和虚部
			if (IData != -0 || IData != 0)
			{
				fp << (IData > 0 ? "+" : "") << IData << "i";
			}
			fp << ",";
		}
		fp << endl;
	}
	fp.close();
	std::cout << path << " 保存成功! " << endl;

	//// 保存二维复数矩阵 data 为二进制文件 data.bin
	//ofstream ofs("..\\Data\\HonstCompleData.bin", ios::binary);
	//size_t nrows = Data.size();
	//size_t ncols = Data[0].size();
	//ofs.write(reinterpret_cast<char*>(&nrows), sizeof(nrows));
	//ofs.write(reinterpret_cast<char*>(&ncols), sizeof(ncols));
	//for (size_t i = 0; i < nrows; ++i) {
	//	for (size_t j = 0; j < ncols; ++j) {
	//		complex<double> value = Data[i][j];
	//		ofs.write(reinterpret_cast<char*>(&value), sizeof(value));
	//	}
	//}
	//ofs.close();
}



//std::vector<std::vector<double>> STFT_lgy::tfrstft(const std::vector<double>& signal, int window_size, int overlap_size)
//{
//	// 创建窗口函数
//	std::vector<double> window(window_size);
//	for (int i = 0; i < window_size; i++) {
//		window[i] = 0.5 * (1 - std::cos(2 * M_PI * i / (window_size - 1)));
//	}
//
//	// 创建FFT计算对象
//	fftw_plan fft = fftw_plan_r2r_1d(window_size, FFTW_R2HC, FFTW_MEASURE);
//
//	// 计算窗口数
//	int num_windows = (signal.size() - window_size) / overlap_size + 1;
//
//	// 创建存储结果的二维数据结构
//	std::vector<std::vector<double>> result(num_windows, std::vector<double>(window_size / 2 + 1));
//
//	// 对每个窗口进行FFT变换
//	for (int i = 0; i < num_windows; i++) {
//		// 计算窗口的起始位置
//		int start_index = i * overlap_size;
//
//		// 应用窗口函数
//		std::vector<double> windowed_signal(window.begin(), window.end());
//		for (int j = 0; j < window_size; j++) {
//			windowed_signal[j] *= signal[start_index + j];
//		}
//
//		// 执行FFT计算
//		fftw_execute(fft);
//
//		// 将结果存储到二维数据结构中
//		for (int j = 0; j < window_size / 2 + 1; j++) {
//			result[i][j] = std::abs(signal[j]);
//		}
//	}
//
//	// 释放FFTW计算对象
//	fftw_destroy_plan(fft);
//
//	return result;
//	//return std::vector<std::vector<double>>();
//}



