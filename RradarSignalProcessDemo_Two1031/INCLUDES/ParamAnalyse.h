#pragma once

#ifndef PARAMANALYSE_H
#define PARAMANALYSE_H

#ifdef PARAMANALYSE_EXPORTS
#define PARAMANALYSE_API __declspec(dllexport)
#else
#define PARAMANALYSE_API __declspec(dllexport)
#endif // PARAMANALYSE_EXPORTS

using namespace std;

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
#include <numeric>
#include "PublicDefinition.h"


//参数估计类
class PARAMANALYSE_API ParamAnalyse
{
public:
	ParamAnalyse();
	~ParamAnalyse();

	vector<PulseBand> ReturnPulseBandData();

	//以BIN文件形式保存，以便MATLAB画图查看
	void RelyOnBINDou(const vector<vector<int>>& NN);

	//小型化的滚动搜索信号拾取
	void AnalyseWinPulse(const vector<int>& DuanYiPulse);

	void WithPulseCalculate(const vector<vector<int>>& SourceData,const vector<int>& UpData, const vector<int>& DownData);


	void ParameterAnalysis(const vector<complex<double>>& AllRadarData, vector<int>& UpData, vector<int>& DownData, vector<int>& PulseWidth, double& time, const double Fs);

	//量化门限函数
	vector<vector<int>> QuanfyThreshold(const vector<vector<double>>& aveVec,const double& AveData);

	//计算vector<vector<double>>结构体的平均值
	double AverageVec(const vector<vector<double>>& aveVec);

	vector<vector<double>> ReturnNeedData(const vector<vector<complex<double>>>& trf);


	//参数估计类的主要执行函数
	void StepAdvance(const vector<vector<complex<double>>>& trf, const vector<complex<double>>& SourceData);

protected:





private:
	//vector<int>      m_InnerDownData{ vector<int>(5,0) };
	//vector<int>      m_InnerUpData{ vector<int>(5,0) };
	//vector<string> name = vector<string>(5);
	//vector<int> val{ vector<int>(5,0) };
	//

	vector<int>      m_InnerDownData;
	vector<int>      m_InnerUpData;

	vector<PulseBand> m_PulseBandVec;



};




#endif