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
#include <windows.h>//����ʱ��
#include "libxl_412\include_cpp\libxl.h"
#include "ParamAnalyse.h"
#include <numeric>
#include "PublicDefinition.h"


//����������
class PARAMANALYSE_API ParamAnalyse
{
public:
	ParamAnalyse();
	~ParamAnalyse();

	vector<PulseBand> ReturnPulseBandData();

	//��BIN�ļ���ʽ���棬�Ա�MATLAB��ͼ�鿴
	void RelyOnBINDou(const vector<vector<int>>& NN);

	//С�ͻ��Ĺ��������ź�ʰȡ
	void AnalyseWinPulse(const vector<int>& DuanYiPulse);

	void WithPulseCalculate(const vector<vector<int>>& SourceData,const vector<int>& UpData, const vector<int>& DownData);


	void ParameterAnalysis(const vector<complex<double>>& AllRadarData, vector<int>& UpData, vector<int>& DownData, vector<int>& PulseWidth, double& time, const double Fs);

	//�������޺���
	vector<vector<int>> QuanfyThreshold(const vector<vector<double>>& aveVec,const double& AveData);

	//����vector<vector<double>>�ṹ���ƽ��ֵ
	double AverageVec(const vector<vector<double>>& aveVec);

	vector<vector<double>> ReturnNeedData(const vector<vector<complex<double>>>& trf);


	//�������������Ҫִ�к���
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