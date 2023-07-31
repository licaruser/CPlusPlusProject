#include "PublicDefinition.h"
using namespace std;



//���ܣ�ģ������״�����ǰ�ص���ʱ��TOA�����ڸ���Ϊ���ݽ���
//���룺RadarNum���״���Ŀ��RadarPRI,ÿ���״ﵽ��ʱ����,RadarSigBegin,�״�TOA��ʼֵ��RadarTOANum���״�TOA����
//������ã�TOATotal���������״��źŵ���ʱ��
void TOA_RadarSignalProduce(int RadarNum, const vector<double> RadarPRI, const vector<double> RadarSigBegin, const vector<int> RadarTOANum, vector<double>& TOATotal);

//���ܣ�ģ���״ﵽ��ʱ��Ӷ�������
//���룺ToaTotal���״�TOA���У�JitterRata��������Χ
//���أ�ToaTotal���״�TOA����
vector<double> RadarToaAddJitter(const vector<double> ToaTotal, const double JitterRata);

//���ܣ�����BeginData��EndData֮���Numλ�����������ΪС�����Precision��0<=BeginData<EndData
//���룺BeginData��EndData��Num��Precision
//�����RandArray
void GetRand(double* RandArray, double BeginData, double EndData, int Num, int Precision);

//���ܣ�����PRI�źŷ�ѡ����
//���룺PRITotal,�״ﵽ��ʱ������ Taumin,Taumax,����������ϢԤ��PRI���·�Χ��K��С�и���,Jitter,����ϵ��
//���أ��������޵�Ǳ��PRIֵ
vector<double> RadarPRISort(const vector<double> PRITotal, double Taumin, double Taumax, int K, double Jitter);

//���ܣ�����PRI��Ϣ�����״��źż�����ɸѡ������ÿ������Դ�ź�
//���룺PotentialPRI����ѡ����PRI�źţ�ToaToatal��������������
//�����UniqueToa����������Դ����ʱ��
void SearchRadarSignal(const vector<double> PotentialPRI, const vector<double> ToaToatal, vector<double> UniqueToa);

//���ܣ�ɾ��������ָ��������Ԫ��
//���룺*dat������ָ�룻*len ���鳤��ָ�룻idx��ɾ��Ԫ������
//�����dat�����º������
void Remove_ArrayElement(double* Data, int* len, int idx)
{
	cout << "���������Ԫ�ظ�����" << *len << "   Ҫɾ�������������ǣ�" << Data[idx] << "   ���������һ�������ǣ�" << Data[*len - 1] << endl;
	if (idx < 0 || idx >= *len) { return; }
	for (int ii = idx; ii < *len; ii++)
	{
		Data[ii] = Data[ii + 1];
	}
	(*len)--;
}




int main()
{
	//LARGE_INTEGER Code_start;
	//LARGE_INTEGER Code_stop;
	//LARGE_INTEGER freq;
	//QueryPerformanceFrequency(&freq);

	//���뿪ʼ����ʱ��ʼ
	//QueryPerformanceCounter(&Code_start);

	//���嵽��ʱ��ģ�⣬ģ��1000������
	vector<double> InitRadarPRI;
	InitRadarPRI.push_back(1);
	InitRadarPRI.push_back(sqrt(2));
	InitRadarPRI.push_back(sqrt(5));

	vector<double> InitRadarSigBegin;
	InitRadarSigBegin.push_back(0);
	InitRadarSigBegin.push_back(0.1);
	InitRadarSigBegin.push_back(0.2);

	vector<int> InitRadarTOANum;
	InitRadarTOANum.push_back(333);
	InitRadarTOANum.push_back(333);
	InitRadarTOANum.push_back(334);

	vector<double> TOATotal;
	TOA_RadarSignalProduce(3, InitRadarPRI, InitRadarSigBegin, InitRadarTOANum, TOATotal);
	//�״ﵽ��ʱ����������
	double JitterToa = 0.0;//��ÿһ���״��źŵ���ʱ������0-10%�����������Χ
	TOATotal = RadarToaAddJitter(TOATotal, JitterToa);
	//�źŲ������֮�󣬽��з�ѡ����
	//����PRI��ѡ
	vector<double> PotentialPRI;//Ǳ��PRIֵ
	double taumin = 0.0;
	double taumax = 10.0;
	int K = 501;
	PotentialPRI = RadarPRISort(TOATotal, taumin, taumax, K, JitterToa);

	//����������ɵ�PRI��Ϣ���������壬ɸѡ����������Դ
	vector<double> UniqueRadatToa;
	SearchRadarSignal(PotentialPRI, TOATotal, UniqueRadatToa);//�����״��źź���

	//�����������ʱ����
	//QueryPerformanceCounter(&Code_stop);
	//double time_sec = (unsigned long long)(Code_stop.QuadPart - Code_start.QuadPart) / (double)freq.QuadPart;
	//cout << "������������ʱ�䣺"<< time_sec << endl;
	return 0;

}


//ģ������״�����ǰ�ص���ʱ��TOA�����ڸ���Ϊ���ݽ���
//���룺RadarNum���״���Ŀ��RadarPRI,ÿ���״ﵽ��ʱ����,RadarSigBegin,�״�TOA��ʼֵ��RadarTOANum���״�TOA����
//�����TOATotal���������״��źŵ���ʱ��
void TOA_RadarSignalProduce(int RadarNum, const vector<double> RadarPRI, const vector<double> RadarSigBegin, const vector<int> RadarTOANum, vector<double>& TOATotal)
{
	//ģ������̶�PRI
	for (int ii = 0; ii < RadarNum; ii++)
	{
		double TemVarPRI = RadarPRI.at(ii);
		double TemVarSigBegin = RadarSigBegin.at(ii);
		int TemVarToaNum = RadarTOANum.at(ii);
		for (int jj = 0; jj < TemVarToaNum; jj++)
		{
			double TemVarData = TemVarSigBegin + jj * TemVarPRI;
			TOATotal.push_back(TemVarData);
		}
	}
	//�������״��źŵ���ʱ������
	sort(TOATotal.begin(), TOATotal.end());


	//for (vector<double>::iterator it = TOATotal.begin(); it != TOATotal.end(); it++)//������
	//{
	//	cout << *it << endl;
	//}
}

//ģ���״ﵽ��ʱ��Ӷ�������
//���룺ToaTotal���״�TOA���У�JitterRata��������Χ
//���أ�ToaTotal���״�TOA����
vector<double> RadarToaAddJitter(const vector<double> ToaTotal, const double JitterRata)
{
	vector<double> TemVarToaTotal;
	TemVarToaTotal = ToaTotal;
	int TotalNumber = TemVarToaTotal.size();
	double* RandArray = new double[TotalNumber];
	GetRand(RandArray, -1, 1, TotalNumber, 3);
	for (int ii = 0; ii < TotalNumber; ii++)
	{
		RandArray[ii] = RandArray[ii] * JitterRata;
		//cout << TemVarToaTotal[ii] << "," << RandArray[ii] << endl;//δ�Ӷ���ǰ���״ﵽ��ʱ��
		TemVarToaTotal[ii] = TemVarToaTotal[ii] + TemVarToaTotal[ii] * RandArray[ii];//�Ӷ���֮����״ﵽ��ʱ��
		//cout << TemVarToaTotal[ii] << endl;
	}

	//for (vector<double>::iterator it = TemVarToaTotal.begin(); it != TemVarToaTotal.end(); it++)//������
	//{
	//	cout << *it << endl;
	//}
	delete[] RandArray;//�ͷ��ڴ�
	return TemVarToaTotal;
}

//���ܣ�����BeginData��EndData֮���Numλ�����������ΪС�����Precision��0<=BeginData<EndData
//���룺BeginData��EndData��Num��Precision
//�����RandArray
void GetRand(double* RandArray, double BeginData, double EndData, int Num, int Precision)
{
	int N;
	if (Precision == 0)
	{
		N = 0;
	}
	else if (Precision == 1)
	{
		N = 9;
	}
	else if (Precision == 2)
	{
		N = 99;
	}
	else if (Precision == 3)
	{
		N = 999;
	}
	else if (Precision == 4)
	{
		N = 9999;
	}
	else
	{
		N = 9999;
		cout << "��λ�����������ƣ�Ĭ��Ϊ4λ" << endl;
	}

	srand(time(NULL));//������������ӣ�ʹÿ�λ�ȡ����������в�ͬ
	double DifferenceData = EndData - BeginData;//�����ֵ
	for (int ii = 0; ii < Num; ii++)
	{
		//����0-1֮��������Num��
		RandArray[ii] = rand() % (N + 1) / (float)(N + 1);
		//ͨ��0-1֮���������������BeginData��EndData֮��������
		RandArray[ii] = BeginData + RandArray[ii] * DifferenceData;
	}
}


//���ܣ�����PRI�źŷ�ѡ����
//���룺PRITotal,�״ﵽ��ʱ������ Taumin,Taumax,����������ϢԤ��PRI���·�Χ��K��С�и���
//���أ��������޵�Ǳ��PRIֵ
vector<double> RadarPRISort(const vector<double> PRITotal, double Taumin, double Taumax, int K, double JitterRata)
{
	double zetazero = 0.03;
	int ToaNum = PRITotal.size();//�״ﵽ��������ܸ���
	double* ToaKCenter = new double[K];//PRIС������ֵ
	double* Bk = new double[K];
	int* Flag = new int[K];
	double* PhaseTimeBegin = new double[K];
	vector<double> PotentialPRI;
	vector<_complex> Drange;//PRI��ָ���ķ����ܺ�
	vector<int> Crange;//PRIδ��ָ���ķ����ܺ�

	_complex TemComData;
	TemComData.x = 0;
	TemComData.y = 0;
	int TemComCrange = 0;
	for (int zz = 0; zz < K; zz++)
	{
		Drange.push_back(TemComData);
		Crange.push_back(TemComCrange);
	}
	//����С������ֵ
	for (int ii = 0; ii < K; ii++)
	{
		ToaKCenter[ii] = (ii + 0.5) * (Taumax - Taumin) / K + Taumin;
		Flag[ii] = 1;
		Bk[ii] = 2 * JitterRata * ToaKCenter[ii];
		//cout << ToaKCenter[ii] << "        " << Flag[ii] << endl;
	}

	int n = 1;
	//��ѡ
	while (n < ToaNum)
	{
		int m;
		m = n - 1;
		while (m >= 0)
		{
			double tau = PRITotal[n] - PRITotal[m];
			if ((tau > (1 - JitterRata) * Taumin) & (tau <= (1 + JitterRata) * Taumax))
			{
				int K1 = trunc((tau / (1 + JitterRata) - Taumin) * K / (Taumax - Taumin) + 1);//����ȡ��
				int K2 = trunc((tau / (1 - JitterRata) - Taumin) * K / (Taumax - Taumin) + 1);
				if (K2 > K)
				{
					break;
				}
				for (int kk = K1 - 1; kk < K2; kk++)
				{
					double Phase;//��λ
					int* Tempflag = &Flag[kk];
					if (*Tempflag == 1)
					{
						double temp = PRITotal[n];
						PhaseTimeBegin[kk] = temp;
					}
					double* PhaseTimeBegina = &PhaseTimeBegin[kk];//�൱��O��k��
					double* ToaKCentera = &ToaKCenter[kk];
					Phase = (PRITotal[n] - *PhaseTimeBegina) / *ToaKCentera;
					double nunu = Phase + 0.4999999;
					double zeta = Phase / nunu - 1;
					int nu = trunc(nunu);
					//cout << "��λ��ʼֵ��" << PhaseTimeBegin[kk] << endl;
					if ((nu == 1) & (PRITotal[m] == PhaseTimeBegin[kk]) | (nu >= 2) & (abs(zeta) <= zetazero))
					{
						PhaseTimeBegin[kk] = PRITotal[n];
					}
					double eta = (PRITotal[n] - PhaseTimeBegin[kk]) / ToaKCenter[kk];
					//cout << "TOAʱ�䣺" << PRITotal[n] << "   " << "��λ��ʼʱ�䣺" << PhaseTimeBegin[kk] << "   " << "С������ֵ��" << ToaKCenter[kk] << "   " << "�������λֵ��" << eta << endl;
					//ָ������
					_complex Phase_Eta;
					Phase_Eta.x = cos(2 * PI * eta);
					Phase_Eta.y = sin(2 * PI * eta);
					Drange[kk].x = Drange[kk].x + Phase_Eta.x;
					Drange[kk].y = Drange[kk].y + Phase_Eta.y;
					Crange[kk] = Crange[kk] + 1;
					Flag[kk] = 0;
				}
			}
			else if (tau > Taumax * (1 + JitterRata))
			{
				break;
			}
			m = m - 1;
		}
		n = n + 1;
	}
	//��ʱ��ѡ����ɣ��������������б���ѡ���������޵�PRI
	vector<double> Model;
	double ModelPhase;
	double ThresholdValueA;//����ֵA
	double ThresholdValueB;//����ֵB
	double ThresholdValueC;//����ֵC
	double ThresholdValueEnd;//��������ֵ

	vector<double> ThresholdValueVector;
	for (int vv = 0; vv < Drange.size(); vv++)
	{
		ModelPhase = sqrt((Drange.at(vv).x) * (Drange.at(vv).x) + (Drange.at(vv).y) * (Drange.at(vv).y));//��ģֵ
		Model.push_back(ModelPhase);

		ThresholdValueA = 225 / *(&ToaKCenter[vv]);//ȡ��ַ֮����ȡֵ�����û�������ˣ�
		ThresholdValueB = 0.15 * Crange[vv];
		ThresholdValueC = 4 * sqrt(pow(ToaNum, 2) * *(&Bk[vv]) / 750);
		ThresholdValueEnd = max(ThresholdValueA, ThresholdValueB);
		ThresholdValueEnd = max(ThresholdValueEnd, ThresholdValueC);

		ThresholdValueVector.push_back(ThresholdValueEnd);

		if (ModelPhase > ThresholdValueEnd)
		{
			PotentialPRI.push_back(*(&ToaKCenter[vv]));//���ô������޵�PRIֵ
		}
	}
	delete[] ToaKCenter;//�ͷ��ڴ�
	delete[] Bk;
	delete[] Flag;
	delete[] PhaseTimeBegin;

	return PotentialPRI;
}



//���ܣ�����PRI��Ϣ�����״��źż�����ɸѡ������ÿ������Դ�ź�
//���룺PotentialPRI����ѡ����PRI�źţ�ToaToatal��������������
//�����UniqueToa����������Դ����ʱ��
void SearchRadarSignal(const vector<double> PotentialPRI, const vector<double> ToaToatal, vector<double> UniqueToa)
{
	int Radar_PRINum = PotentialPRI.size();//��ѡ������PRI����Ŀ
	int Radar_ToaNum = ToaToatal.size();//�״�ȫ������TOA����Ŀ
	int GataNum = 0;//�ж�PRI��Ч����
	int Sig_Total = 0;//����Դ�źŶ�λ��
	double * Data = new double[Radar_ToaNum];//���ڴ�Ÿ�������Դ�ź�
	if (Radar_PRINum == 0)
	{
		//PRI��ѡ�޽��
		cout << "��������PRI�źţ����ѡʧ��" << endl;
	}
	else
	{
		//PRI��ѡ�н���������������м�����Ҳ����������������Դ�ź�
		//������ҵ���PRI���д��������������У����ι������δ�滮���������塱�͡����嶪ʧ�����
		for (int ii = 0; ii < Radar_PRINum; ii++)
		{
			//����ÿһ������
			for (int jj = 0; jj < Radar_ToaNum; jj++)
			{
				double Tem_Toa = PotentialPRI.at(ii) + ToaToatal.at(jj);//��Toa�ĵ�0λ��ʼ����֮ǰ�������PRI
				for (int zz = jj + 1; zz < Radar_ToaNum; zz++)//����һλ��ʼ�����Ƿ��е������������Tem_Toa
				{
					if (abs(Tem_Toa - ToaToatal.at(zz)) < 1e-2)//Tem_Toa = Toatoatal.at(zz)ʱ��
					{
						GataNum = GataNum + 1;
						break;
					}
				}
				if (GataNum > 5) //˵�����PRIȷʵ��Ч�������޷�Ӧ�����嶪ʧ�����嶶��
				{
					//���м���ǰ����
					for (int dd = jj; dd > 0; dd--)
					{
						Tem_Toa = ToaToatal.at(dd) - PotentialPRI.at(ii);
						//�������˵�һ�������ǵڶ���ֵ,Ҳ����0
						/*if (abs(Tem_Toa) < 1e-2)
						{
							Data[Sig_Total] = ToaToatal.at(dd);
							cout << "�����������Դ�ź�TOA��" << Data[Sig_Total] << endl;
							Sig_Total = Sig_Total + 1;
							break;
						}*/
						//��ǰ����ʱ��������ǰ�ж�
						for (int ww = dd - 1; ww > -1; ww--)
						{
							if (abs(Tem_Toa - ToaToatal.at(ww)) < 1e-2)
							{
								Data[Sig_Total] = ToaToatal.at(ww);
								//cout << "�����������Դ�ź�TOA��" << Data[Sig_Total] << endl;
								Sig_Total = Sig_Total + 1;
								break;
							}
						}
					}
					//����ȡ��
					int CenterData = trunc(Sig_Total / 2);
					//����ǰ��Ե�
					for (int ppi = 0; ppi < CenterData; ppi++)
					{
						double Tem_Toa = Data[Sig_Total - ppi - 1];
						Data[Sig_Total - ppi - 1] = Data[ppi];
						Data[ppi] = Tem_Toa;						
					}
					Data[Sig_Total] = ToaToatal.at(jj);
					Sig_Total = Sig_Total + 1;
					//������
					for (int ppt = jj; ppt < Radar_ToaNum - 2; ppt++)
					{
						Tem_Toa = ToaToatal.at(ppt) + PotentialPRI.at(ii);
						for (int wwi = ppt + 1; wwi < Radar_ToaNum - 1; wwi++)
						{
							if (abs(Tem_Toa - ToaToatal.at(wwi)) < 1e-2)
							{
								Data[Sig_Total] = ToaToatal.at(wwi);
								Sig_Total = Sig_Total + 1;
							}
						}
					}

					//�������
					for (int ll = 0; ll < Sig_Total; ll++)
					{
						cout << "�����������Դ�ź�TOA��" << Data[ll] << endl;
					}
					int aa = 1;


					//��ģ����������ź�TOA���Ӻ�PRI��ѡ�������������������������׼ȷ���ʲ���ǰ��ȶ��ټ������ҳ����ϱ�PRI���״��ź�
					double DifferenceDataA_B;
					for (int Diff = 0; Diff < Sig_Total; Diff++)
					{
						DifferenceDataA_B = Data[Diff + 1] - Data[Diff];//�������������ݺ���һ����ȥǰ��һ��
						if (DifferenceDataA_B > PotentialPRI.at(ii)*1.03)
						{
							//������ж��ڣ�˵��ǰ�����������⣬����130%��˵���������ݴ��󣬻����Ƕ�ʧ���壡
							//�������嶪ʧ���ʣ��жϲ�ֵ�Ƿ��ǵ�ǰPRI����������������������˵������ʱ���嶪ʧ����������������˵���������ֲ��
							int Flag_Error = 0;
							for (int zfii = 2; zfii < 6; zfii++)//�ݶ�2����5�����
							{
								if (abs(PotentialPRI.at(ii) * zfii - DifferenceDataA_B) < 1e-2)//0.01
								{
									cout << "���������ǰ���ֵΪ��ǰPRI��" << zfii << "�����м䶪ʧ������Ϊ" << zfii - 1 << endl;
									//����������ɾ����Ӧ����ȥԭʼ�����в��Ҹ�����ķ���Դ���ݣ�ȷ���Ƿ��Ǵ�����������ﲻ��ִ���������裡
									Flag_Error = 1;
								}
							}
							if (Flag_Error == 1)
							{
								//˵���������嶪ʧ�����Ǵ���
							}
							else
							{
								//˵���������嶪ʧ���Ǽ��������ˣ���������������
								Remove_ArrayElement(Data, &Sig_Total, Diff + 1);//ɾ��֮��õ����µ�����
								//Sig_Total�Ǹ÷���Դ������Diff+1��Ҫɾ�����±������ݴ��㿪ʼ���±�Ҫ��������1
								if (Sig_Total > Diff + 1)
								{
									Diff = Diff - 1;//ɾ��֮�������·����ж�һ��
								}
							}
						}
						if (DifferenceDataA_B < PotentialPRI.at(ii) * 0.97)
						{
							//������ж��ڣ�˵��ǰ�����������⣬С��70%��˵���϶�����֮ǰ�ļ��������ˣ�
							//��Ҫ��һ���жϵ������ĸ�����������
							double DataDiff_DiffMinus = Data[Diff] - Data[Diff - 1];
							double DataDiffAdd_DiffMinus = Data[Diff + 1] - Data[Diff - 1];
							if (abs(DataDiff_DiffMinus - PotentialPRI.at(ii)) < abs(DataDiffAdd_DiffMinus - PotentialPRI.at(ii)))//Data[Diff]��Data[Diff-1]��PRI�����PRI�Ĳ�ֵ   С��  ��Data[Diff + 1]��Data[Diff - 1]��ǰǰһ��PRI�����PRI�Ĳ�ֵ����˵��ǰ������ݸ��ӽ�����PRI
							{
								//�Ѻ��������ֵ�ϴ������ɾ����
								Remove_ArrayElement(Data, &Sig_Total, Diff + 1);//ɾ��֮��õ����µ�����
								//Sig_Total�Ǹ÷���Դ������Diff+1��Ҫɾ�����±������ݴ��㿪ʼ���±�Ҫ��������1
								if (Sig_Total > Diff + 1)
								{
									Diff = Diff - 1;//ɾ��֮�������·����ж�һ��
								}
							}
							else
							{
								//���Ǿ���ǰ��������ڵ��ں��������˵������������ӽ���ѡ��PRI��ɾ��ǰ�����
								Remove_ArrayElement(Data, &Sig_Total, Diff);//ɾ��֮��õ����µ�����
								//Sig_Total�Ǹ÷���Դ������Diff+1��Ҫɾ�����±������ݴ��㿪ʼ���±�Ҫ��������1
								if (Sig_Total > Diff + 1)
								{
									Diff = Diff - 1;//ɾ��֮�������·����ж�һ��
								}
							}
						}
					}
					

					//�������
					for (int ll = 0; ll < Sig_Total ; ll++)
					{
						cout << "�����������Դ�ź�TOA��" << Data[ll] << endl;
					}
					int aab = 1;
				}

			}
		}
	}
}