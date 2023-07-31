#pragma once
#ifndef TEST_STFT_DEMO_h__
#define TEST_STFT_DEMO_h__

#ifdef	TEST_STFT_DEMO_EXPORTS
#define TEST_STFT_DEMO_API __declspec(dllexport)
#else		  
#define TEST_STFT_DEMO_API __declspec(dllimport)
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <string.h>
#include <vector>
#include <Dense>
#include <windows.h>//精度时间
#include "libxl_412\include_cpp\libxl.h"
#include "ParamAnalyse.h"
#include "PublicDefinition.h"


#define M_PI 3.1415926535897932384626433

using namespace std;
using namespace libxl;
using Eigen::VectorXd;


class TEST_STFT_DEMO_API STFT_lgy
{
public:
	STFT_lgy();
	~STFT_lgy();

	void Read_Excel_Data(int num);

	void Read_Csv_Data(std::vector<complex<double>>& SignalData, int Init_row, int Init_col);

	std::vector<std::vector<double>> tfrstft(const std::vector<double>& signal, int window_size, int overlap_size);

	//Matlab中tfrstft复现，trf 输出二维矩阵，f 输出频率， x 输入原始数据， t 输入时间， N 输入频率分辨率点数
	void tfrstft(vector<vector<complex<double>>>& trf, double& f, const vector<complex<double>>& x, vector<int>& t, int& N);

	int min_member(int a, int b, int c);

	void SaveData_Csv(const vector<vector<complex<double>>>& Data);

protected:





private:
	//ParamAnalyse *m_TmpParamEsta; //参数估计的类



};




#endif