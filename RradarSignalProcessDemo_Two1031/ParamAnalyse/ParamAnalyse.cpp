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
	//�������ݵ�BIN�ļ���
	std::ofstream file("..\\Data\\data.bin", std::ios::binary);
	if (file.is_open()) {
		// ������д���ļ�
		for (const auto& iter : NN)
		{
			/*cout << sizeof(int) << endl;
			bbc = iter;*/
			file.write(reinterpret_cast<const char*>(iter.data()), iter.size() * sizeof(int));
			//std::cout << (reinterpret_cast<const char*>(bbc.data()), bbc.size() * sizeof(int)) << std::endl;
		}
		file.close();
		std::cout << "�����ѱ���ΪBIN�ļ�" << std::endl;
	}
	else {
		std::cout << "�޷����ļ�" << std::endl;
	}
}

void ParamAnalyse::AnalyseWinPulse(const vector<int>& DuanYiPulse)
{
	// ����AllRadarData�������ء��½��غ������������Ӧ��vector�У����ں������з�ѡ����
	int PulseUpFlag = 0;
	int PulseDownFlag = 0;
	int PulseUpPos = 0;
	int PulseDownPos = 0;

	vector<double> DataRead;
	DataRead.resize(DuanYiPulse.size());
	for (int jj = 0; jj < DuanYiPulse.size(); jj++)
	{
		DataRead[jj] = DuanYiPulse[jj] * DuanYiPulse[jj];   //���㹦��
	}
	
	// ������ò���
	int PulseFlag = 0;
	int	CntBegin = 0;
	int CntEnd = 0;
	double NoisePower = 0.0;
	double DataPower = 0.0;
	int BeginFrameThrh = 4;
	int EndFrameThrh = 4;
	int FrameSize = 4; //֡����
	int NoiseSize = DuanYiPulse.size();
	NoisePower = accumulate(DataRead.begin(), DataRead.begin() + NoiseSize, 0.0) / NoiseSize; //����ǰȫ������������ʾ�ֵ

	vector<double> data_buff;
	data_buff.resize(FrameSize);
	// 3��������ⷽ������������
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
		//	int oo = 0;  //����
		//}
		DataPower = accumulate(data_buff.begin(), data_buff.end(), 0.0) / data_buff.size();  //�����֡�ź�ƽ������

		if (DataPower > 1.5 * NoisePower)
		{
			CntEnd = 0;
			CntBegin = CntBegin + 1;

			if (CntBegin >= BeginFrameThrh) // �����г���BeginFrameThrh֡�źų���������ޣ�����Ϊ�����忪ʼ
			{
				if (PulseFlag == 0) // ��ǰ��δ��⵽����
				{
					// Ѱ�Ҿ�ȷ��������ʼλ��
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
				if (CntEnd >= EndFrameThrh)  //����������EndFrameThrh֡�źŵ��ڼ�����ޣ�����Ϊ���������
				{
					// Ѱ�Ҿ�ȷ���������λ��
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

	// 4����������ʱ��--������Posλ�û��㵽ʱ����(��λ/����)
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
	//�������ֵ������
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
			////��һ���������һ����ʱ���µĸ���Ҷ�任���
			//FunInnerTmpUpData.push_back(*min_element(m_InnerUpData.begin(), m_InnerUpData.end()));
			//funInnerTmpDownData.push_back(*max_element(m_InnerDownData.begin(), m_InnerDownData.end()));
		}

		//˽�г�Ա�ѱ���ֵ
		PulseBand TmpPulseData;
		TmpPulseData.minGle = *min_element(m_InnerUpData.begin(), m_InnerUpData.end());
		TmpPulseData.maxGle = *max_element(m_InnerDownData.begin(), m_InnerDownData.end());
		TmpPulseData.BandWidth = TmpPulseData.maxGle - TmpPulseData.minGle;

		m_PulseBandVec.push_back(TmpPulseData);
	}
}




void ParamAnalyse::ParameterAnalysis(const vector<complex<double>>& AllRadarData, vector<int>& UpData, vector<int>& DownData, vector<int>& PulseWidth, double& time, const double Fs)
{
	// ����AllRadarData�������ء��½��غ������������Ӧ��vector�У����ں������з�ѡ����
	int PulseUpFlag = 0;
	int PulseDownFlag = 0;
	int PulseUpPos = 0;
	int PulseDownPos = 0;

	// 1������absֵ
	vector<double> HeDataABS;
	HeDataABS.resize(AllRadarData.size());
	for (int ii = 0; ii < AllRadarData.size(); ii++)
	{
		float RealData = AllRadarData.at(ii).real();
		float ImagData = AllRadarData.at(ii).imag();
		HeDataABS[ii] = sqrt(RealData * RealData + ImagData * ImagData);
	}

	// 2������ǰ�洿�������ּ����������ʺͷ���
	vector<double> DataRead;
	DataRead.resize(HeDataABS.size());
	for (int jj = 0; jj < HeDataABS.size(); jj++)
	{
		DataRead[jj] = HeDataABS[jj] * HeDataABS[jj];   //���㹦��
	}
	HeDataABS.clear();

	// ������ò���
	int PulseFlag = 0;
	int	CntBegin = 0;
	int CntEnd = 0;
	double NoisePower = 0.0;
	double DataPower = 0.0;
	int BeginFrameThrh = 128;
	int EndFrameThrh = 128;
	int FrameSize = 128; //֡����
	NoisePower = accumulate(DataRead.begin(), DataRead.begin() + FrameSize, 0.0) / FrameSize; //����ǰ128������������ʾ�ֵ

	vector<double> data_buff;
	data_buff.resize(FrameSize);
	// 3��������ⷽ������������
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
		//	int oo = 0;  //����
		//}
		DataPower = accumulate(data_buff.begin(), data_buff.end(), 0.0) / data_buff.size();  //�����֡�ź�ƽ������

		if (DataPower > 1.5 * NoisePower)
		{
			CntEnd = 0;
			CntBegin = CntBegin + 1;

			if (CntBegin >= BeginFrameThrh) // �����г���BeginFrameThrh֡�źų���������ޣ�����Ϊ�����忪ʼ
			{
				if (PulseFlag == 0) // ��ǰ��δ��⵽����
				{
					// Ѱ�Ҿ�ȷ��������ʼλ��
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
				if (CntEnd >= EndFrameThrh)  //����������EndFrameThrh֡�źŵ��ڼ�����ޣ�����Ϊ���������
				{
					// Ѱ�Ҿ�ȷ���������λ��
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

	// 4����������ʱ��--������Posλ�û��㵽ʱ����(��λ/����)
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



	//// 4����������ʱ��--������Posλ�û��㵽ʱ����(��λ/����)
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
		for (int jj = 0; jj < trf.at(0).size()/2; jj++)  //Ĭ���Ӽ���size���
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
	//ִ�н�ȡ��ģֵ����
	vector<vector<double>> SourceDouData = ReturnNeedData(trf);

	//��������ֵ��ƽ��ֵ���ҵ����ʵ�����
	double AveVectorDatal = AverageVec(SourceDouData);

	//�����������޵�ʱƵ����ֵ
	vector<vector<int>> NN = QuanfyThreshold(SourceDouData, AveVectorDatal);

	//��������֮��Ľ����ȷ������ֵ���÷���Ԥ�ڣ�ʹ��MATLAB��ͼ
	RelyOnBINDou(NN);


	//����ʱ�����ݼ��������������½���ʱ��
	vector<int> UpVec;           //������vector
	vector<int> DownVec;         //�½���vector
	vector<int> PulseWidthVec;   //��Ӧ����
	double CurrentTime = 0.0;
	double m_JamInitSampleFs = 0.0;
	//m_pSignalProcessor->GetJamInitSampleFs(m_JamInitSampleFs);
	// ��������--�������������ء��½��ء�����
	if (SourceData.size() > 0)
	{
		ParameterAnalysis(SourceData, UpVec, DownVec, PulseWidthVec, CurrentTime, m_JamInitSampleFs);//����״��źŲ�����λ(����)[�����ء��½��ء�����]
	}

	//������������ʱ�̼��������,������洢��˽�г�Ա������
	WithPulseCalculate(NN, UpVec, DownVec);


	//���ݼ����˽�г�Ա
	vector<PulseBand> DataPulseBand;
	
	DataPulseBand = ReturnPulseBandData();


	int ada = 0;
}
