#include "ParamAnalyse.h"




ParamAnalyse::ParamAnalyse()
	//:m_InnerUpData(5,0)
	//,m_InnerDownData(5,0)
{
	//int aaa3 = 0;
}

ParamAnalyse::~ParamAnalyse()
{
}

vector<PulseBand> ParamAnalyse::ReturnPulseBandData()
{
	return m_PulseBandVec;
}

void ParamAnalyse::RelyOnBINDou(const vector<vector<int>>& NN)
{
	//vector<int> bbc;
	//保存数据到BIN文件中
	std::ofstream file("..\\Data\\data.bin", std::ios::binary);
	if (file.is_open()) {
		// 将数据写入文件
		for (const auto& iter : NN)
		{
			/*cout << sizeof(int) << endl;
			bbc = iter;*/
			file.write(reinterpret_cast<const char*>(iter.data()), iter.size() * sizeof(int));
			//std::cout << (reinterpret_cast<const char*>(bbc.data()), bbc.size() * sizeof(int)) << std::endl;
		}
		file.close();
		std::cout << "数据已保存为BIN文件" << std::endl;
	}
	else {
		std::cout << "无法打开文件" << std::endl;
	}
}

void ParamAnalyse::AnalyseWinPulse(const vector<int>& DuanYiPulse)
{
	// 估计AllRadarData的上升沿、下降沿和脉宽，并存入对应的vector中，便于后续进行分选处理
	int PulseUpFlag = 0;
	int PulseDownFlag = 0;
	int PulseUpPos = 0;
	int PulseDownPos = 0;

	vector<double> DataRead;
	DataRead.resize(DuanYiPulse.size());
	for (int jj = 0; jj < DuanYiPulse.size(); jj++)
	{
		DataRead[jj] = DuanYiPulse[jj] * DuanYiPulse[jj];   //计算功率
	}
	
	// 检测所用参数
	int PulseFlag = 0;
	int	CntBegin = 0;
	int CntEnd = 0;
	double NoisePower = 0.0;
	double DataPower = 0.0;
	int BeginFrameThrh = 4;
	int EndFrameThrh = 4;
	int FrameSize = 4; //帧长度
	int NoiseSize = DuanYiPulse.size();
	NoisePower = accumulate(DataRead.begin(), DataRead.begin() + NoiseSize, 0.0) / NoiseSize; //计算前全部点的噪声功率均值

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
					PulseUpPos = k - FrameSize + 1;
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
					PulseDownPos = k - FrameSize - FrameSize + 1;
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
		m_InnerUpData.push_back(UpPoint.at(aa));
	}
	for (int bb = 0; bb < DownPoint.size(); bb++)
	{
		//double Tmp_DownTime;
		//Tmp_DownTime = ((time + DownPoint[bb] / Fs) * 1e3);
		m_InnerDownData.push_back(DownPoint.at(bb));
	}


}

void ParamAnalyse::WithPulseCalculate(const vector<vector<int>>& SourceInnerData, const vector<int>& UpData, const vector<int>& DownData)
{
	vector<int> InnerData;
	//vector<int> FunInnerTmpUpData;
	//vector<int> funInnerTmpDownData;
	//计算过幅值的数据
	for (int ii = 0; ii < UpData.size(); ii++)
	{
		//FunInnerTmpUpData.clear();
		//funInnerTmpDownData.clear();
		m_InnerUpData.clear();
		m_InnerDownData.clear();
		for (int jj = UpData.at(ii); jj < DownData.at(ii); jj++)
		{
			if (jj == 8723)
			{
				int ooc = 0;
			}
			
			InnerData = SourceInnerData.at(jj);
			AnalyseWinPulse(InnerData);
			////第一次输出的是一个点时间下的傅里叶变换结果
			//FunInnerTmpUpData.push_back(*min_element(m_InnerUpData.begin(), m_InnerUpData.end()));
			//funInnerTmpDownData.push_back(*max_element(m_InnerDownData.begin(), m_InnerDownData.end()));
		}

		//私有成员已被赋值
		PulseBand TmpPulseData;
		TmpPulseData.minGle = *min_element(m_InnerUpData.begin(), m_InnerUpData.end());
		TmpPulseData.maxGle = *max_element(m_InnerDownData.begin(), m_InnerDownData.end());
		TmpPulseData.BandWidth = TmpPulseData.maxGle - TmpPulseData.minGle;

		m_PulseBandVec.push_back(TmpPulseData);
	}
}




void ParamAnalyse::ParameterAnalysis(const vector<complex<double>>& AllRadarData, vector<int>& UpData, vector<int>& DownData, vector<int>& PulseWidth, double& time, const double Fs)
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

vector<vector<int>> ParamAnalyse::QuanfyThreshold(const vector<vector<double>>& aveVec, const double& AveData)
{
	vector<vector<int> > NN(aveVec.size(), vector<int>(aveVec.at(0).size()));
	for (int ii = 0; ii < NN.size(); ii++)
	{
		for (int jj = 0; jj < NN.at(0).size(); jj++)
		{
			NN.at(ii).at(jj) = 0;
		}
	}

	for (int ii = 0; ii < aveVec.size(); ii++)
	{
		vector<double> TempSourData;
		TempSourData = aveVec.at(ii);
		for (int jj = 0; jj < TempSourData.size(); jj++)
		{
			if (TempSourData.at(jj) > AveData * 30)
			{
				NN.at(ii).at(jj) = 1;
			}
		}
	}

	return NN;
}

double ParamAnalyse::AverageVec(const vector<vector<double>>& aveVec)
{
	vector<double> TempAveData;
	double sumValue = 0.0;
	double aveVecData = 0.0;
	for (int ii = 0; ii < aveVec.size(); ii++)
	{
		sumValue = accumulate(begin(aveVec.at(ii)), end(aveVec.at(ii)), 0.0);
		aveVecData = sumValue / aveVec.at(ii).size();
		TempAveData.push_back(aveVecData);
	}

	double sumTotalValue = accumulate(begin(TempAveData), end(TempAveData), 0.0);
	double aveVecTotalData = sumTotalValue / TempAveData.size();

	return aveVecTotalData;
}

vector<vector<double>> ParamAnalyse::ReturnNeedData(const vector<vector<complex<double>>>& trf)
{
	vector<vector<double>> tmpData;
	vector<double> OtherTmpData;
	for (int ii = 0; ii < trf.size(); ii++)
	{
		for (int jj = 0; jj < trf.at(0).size()/2; jj++)  //默认子集中size相等
		{
			OtherTmpData.push_back(sqrt(pow(trf.at(ii).at(jj).real(), 2) + pow(trf.at(ii).at(jj).imag(), 2)));
		}
		tmpData.push_back(OtherTmpData);
		OtherTmpData.clear();
	}

	return tmpData;
}

void ParamAnalyse::StepAdvance(const vector<vector<complex<double>>> &trf,const vector<complex<double>>& SourceData)
{
	//执行截取及模值计算
	vector<vector<double>> SourceDouData = ReturnNeedData(trf);

	//计算所有值的平均值，找到合适的门限
	double AveVectorDatal = AverageVec(SourceDouData);

	//量化超过门限的时频分析值
	vector<vector<int>> NN = QuanfyThreshold(SourceDouData, AveVectorDatal);

	//保存量化之后的结果，确保门限值设置符合预期，使用MATLAB画图
	RelyOnBINDou(NN);


	//根据时域数据计算脉冲上升沿下降沿时刻
	vector<int> UpVec;           //上升沿vector
	vector<int> DownVec;         //下降沿vector
	vector<int> PulseWidthVec;   //对应脉宽
	double CurrentTime = 0.0;
	double m_JamInitSampleFs = 0.0;
	//m_pSignalProcessor->GetJamInitSampleFs(m_JamInitSampleFs);
	// 参数估计--估计脉冲上升沿、下降沿、脉宽
	if (SourceData.size() > 0)
	{
		ParameterAnalysis(SourceData, UpVec, DownVec, PulseWidthVec, CurrentTime, m_JamInitSampleFs);//输出雷达信号参数单位(毫秒)[上升沿、下降沿、脉宽]
	}

	//根据脉冲所在时刻计算其带宽,将结果存储于私有成员变量中
	WithPulseCalculate(NN, UpVec, DownVec);


	//根据计算的私有成员
	vector<PulseBand> DataPulseBand;
	
	DataPulseBand = ReturnPulseBandData();


	int ada = 0;
}
