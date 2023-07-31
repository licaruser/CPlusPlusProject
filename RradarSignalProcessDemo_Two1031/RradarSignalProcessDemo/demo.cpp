#include "PublicDefinition.h"
using namespace std;



//功能：模拟产生雷达脉冲前沿到达时间TOA，后期更改为数据接收
//输入：RadarNum，雷达数目，RadarPRI,每部雷达到达时间间隔,RadarSigBegin,雷达TOA起始值，RadarTOANum，雷达TOA数量
//输出引用：TOATotal，多脉冲雷达信号到达时间
void TOA_RadarSignalProduce(int RadarNum, const vector<double> RadarPRI, const vector<double> RadarSigBegin, const vector<int> RadarTOANum, vector<double>& TOATotal);

//功能：模拟雷达到达时间加抖动处理
//输入：ToaTotal，雷达TOA序列，JitterRata，抖动范围
//返回：ToaTotal，雷达TOA序列
vector<double> RadarToaAddJitter(const vector<double> ToaTotal, const double JitterRata);

//功能：产生BeginData到EndData之间的Num位随机数，精度为小数点后Precision，0<=BeginData<EndData
//输入：BeginData，EndData，Num，Precision
//输出：RandArray
void GetRand(double* RandArray, double BeginData, double EndData, int Num, int Precision);

//功能：进行PRI信号分选工作
//输入：PRITotal,雷达到达时间序列 Taumin,Taumax,根据先验信息预设PRI大致范围，K，小盒个数,Jitter,抖动系数
//返回：超过门限的潜在PRI值
vector<double> RadarPRISort(const vector<double> PRITotal, double Taumin, double Taumax, int K, double Jitter);

//功能：根据PRI信息进行雷达信号检索，筛选出具体每部辐射源信号
//输入：PotentialPRI，分选出的PRI信号，ToaToatal，多脉冲数据流
//输出：UniqueToa，独立辐射源到达时间
void SearchRadarSignal(const vector<double> PotentialPRI, const vector<double> ToaToatal, vector<double> UniqueToa);

//功能：删除数组中指定索引的元素
//输入：*dat，数组指针；*len 数组长度指针；idx，删除元素索引
//输出：dat，更新后的数组
void Remove_ArrayElement(double* Data, int* len, int idx)
{
	cout << "输出数组中元素个数：" << *len << "   要删除的索引数据是：" << Data[idx] << "   数组中最后一个数据是：" << Data[*len - 1] << endl;
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

	//代码开始，计时开始
	//QueryPerformanceCounter(&Code_start);

	//脉冲到达时间模拟，模拟1000个脉冲
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
	//雷达到达时间加随机抖动
	double JitterToa = 0.0;//给每一个雷达信号到达时间增加0-10%的随机抖动范围
	TOATotal = RadarToaAddJitter(TOATotal, JitterToa);
	//信号产生完成之后，进行分选工作
	//进行PRI分选
	vector<double> PotentialPRI;//潜在PRI值
	double taumin = 0.0;
	double taumax = 10.0;
	int K = 501;
	PotentialPRI = RadarPRISort(TOATotal, taumin, taumax, K, JitterToa);

	//根据搜索完成的PRI信息，检索脉冲，筛选出各个辐射源
	vector<double> UniqueRadatToa;
	SearchRadarSignal(PotentialPRI, TOATotal, UniqueRadatToa);//检索雷达信号函数

	//代码结束，计时结束
	//QueryPerformanceCounter(&Code_stop);
	//double time_sec = (unsigned long long)(Code_stop.QuadPart - Code_start.QuadPart) / (double)freq.QuadPart;
	//cout << "整个程序运行时间："<< time_sec << endl;
	return 0;

}


//模拟产生雷达脉冲前沿到达时间TOA，后期更改为数据接收
//输入：RadarNum，雷达数目，RadarPRI,每部雷达到达时间间隔,RadarSigBegin,雷达TOA起始值，RadarTOANum，雷达TOA数量
//输出：TOATotal，多脉冲雷达信号到达时间
void TOA_RadarSignalProduce(int RadarNum, const vector<double> RadarPRI, const vector<double> RadarSigBegin, const vector<int> RadarTOANum, vector<double>& TOATotal)
{
	//模拟产生固定PRI
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
	//多脉冲雷达信号到达时间排序
	sort(TOATotal.begin(), TOATotal.end());


	//for (vector<double>::iterator it = TOATotal.begin(); it != TOATotal.end(); it++)//迭代器
	//{
	//	cout << *it << endl;
	//}
}

//模拟雷达到达时间加抖动处理
//输入：ToaTotal，雷达TOA序列，JitterRata，抖动范围
//返回：ToaTotal，雷达TOA序列
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
		//cout << TemVarToaTotal[ii] << "," << RandArray[ii] << endl;//未加抖动前的雷达到达时间
		TemVarToaTotal[ii] = TemVarToaTotal[ii] + TemVarToaTotal[ii] * RandArray[ii];//加抖动之后的雷达到达时间
		//cout << TemVarToaTotal[ii] << endl;
	}

	//for (vector<double>::iterator it = TemVarToaTotal.begin(); it != TemVarToaTotal.end(); it++)//迭代器
	//{
	//	cout << *it << endl;
	//}
	delete[] RandArray;//释放内存
	return TemVarToaTotal;
}

//功能：产生BeginData到EndData之间的Num位随机数，精度为小数点后Precision，0<=BeginData<EndData
//输入：BeginData，EndData，Num，Precision
//输出：RandArray
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
		cout << "数位超过设置限制，默认为4位" << endl;
	}

	srand(time(NULL));//设置随机数种子，使每次获取的随机数序列不同
	double DifferenceData = EndData - BeginData;//计算差值
	for (int ii = 0; ii < Num; ii++)
	{
		//产生0-1之间的随机数Num个
		RandArray[ii] = rand() % (N + 1) / (float)(N + 1);
		//通过0-1之间的随机数，计算出BeginData到EndData之间的随机数
		RandArray[ii] = BeginData + RandArray[ii] * DifferenceData;
	}
}


//功能：进行PRI信号分选工作
//输入：PRITotal,雷达到达时间序列 Taumin,Taumax,根据先验信息预设PRI大致范围，K，小盒个数
//返回：超过门限的潜在PRI值
vector<double> RadarPRISort(const vector<double> PRITotal, double Taumin, double Taumax, int K, double JitterRata)
{
	double zetazero = 0.03;
	int ToaNum = PRITotal.size();//雷达到达脉冲的总个数
	double* ToaKCenter = new double[K];//PRI小盒中心值
	double* Bk = new double[K];
	int* Flag = new int[K];
	double* PhaseTimeBegin = new double[K];
	vector<double> PotentialPRI;
	vector<_complex> Drange;//PRI加指数的幅度总和
	vector<int> Crange;//PRI未加指数的幅度总和

	_complex TemComData;
	TemComData.x = 0;
	TemComData.y = 0;
	int TemComCrange = 0;
	for (int zz = 0; zz < K; zz++)
	{
		Drange.push_back(TemComData);
		Crange.push_back(TemComCrange);
	}
	//计算小盒中心值
	for (int ii = 0; ii < K; ii++)
	{
		ToaKCenter[ii] = (ii + 0.5) * (Taumax - Taumin) / K + Taumin;
		Flag[ii] = 1;
		Bk[ii] = 2 * JitterRata * ToaKCenter[ii];
		//cout << ToaKCenter[ii] << "        " << Flag[ii] << endl;
	}

	int n = 1;
	//分选
	while (n < ToaNum)
	{
		int m;
		m = n - 1;
		while (m >= 0)
		{
			double tau = PRITotal[n] - PRITotal[m];
			if ((tau > (1 - JitterRata) * Taumin) & (tau <= (1 + JitterRata) * Taumax))
			{
				int K1 = trunc((tau / (1 + JitterRata) - Taumin) * K / (Taumax - Taumin) + 1);//向零取整
				int K2 = trunc((tau / (1 - JitterRata) - Taumin) * K / (Taumax - Taumin) + 1);
				if (K2 > K)
				{
					break;
				}
				for (int kk = K1 - 1; kk < K2; kk++)
				{
					double Phase;//相位
					int* Tempflag = &Flag[kk];
					if (*Tempflag == 1)
					{
						double temp = PRITotal[n];
						PhaseTimeBegin[kk] = temp;
					}
					double* PhaseTimeBegina = &PhaseTimeBegin[kk];//相当于O（k）
					double* ToaKCentera = &ToaKCenter[kk];
					Phase = (PRITotal[n] - *PhaseTimeBegina) / *ToaKCentera;
					double nunu = Phase + 0.4999999;
					double zeta = Phase / nunu - 1;
					int nu = trunc(nunu);
					//cout << "相位初始值：" << PhaseTimeBegin[kk] << endl;
					if ((nu == 1) & (PRITotal[m] == PhaseTimeBegin[kk]) | (nu >= 2) & (abs(zeta) <= zetazero))
					{
						PhaseTimeBegin[kk] = PRITotal[n];
					}
					double eta = (PRITotal[n] - PhaseTimeBegin[kk]) / ToaKCenter[kk];
					//cout << "TOA时间：" << PRITotal[n] << "   " << "相位起始时间：" << PhaseTimeBegin[kk] << "   " << "小盒中心值：" << ToaKCenter[kk] << "   " << "计算的相位值：" << eta << endl;
					//指数计算
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
	//此时分选已完成，后续进行门限判别，挑选出超过门限的PRI
	vector<double> Model;
	double ModelPhase;
	double ThresholdValueA;//门限值A
	double ThresholdValueB;//门限值B
	double ThresholdValueC;//门限值C
	double ThresholdValueEnd;//最后的门限值

	vector<double> ThresholdValueVector;
	for (int vv = 0; vv < Drange.size(); vv++)
	{
		ModelPhase = sqrt((Drange.at(vv).x) * (Drange.at(vv).x) + (Drange.at(vv).y) * (Drange.at(vv).y));//求模值
		Model.push_back(ModelPhase);

		ThresholdValueA = 225 / *(&ToaKCenter[vv]);//取地址之后再取值好像就没有提醒了？
		ThresholdValueB = 0.15 * Crange[vv];
		ThresholdValueC = 4 * sqrt(pow(ToaNum, 2) * *(&Bk[vv]) / 750);
		ThresholdValueEnd = max(ThresholdValueA, ThresholdValueB);
		ThresholdValueEnd = max(ThresholdValueEnd, ThresholdValueC);

		ThresholdValueVector.push_back(ThresholdValueEnd);

		if (ModelPhase > ThresholdValueEnd)
		{
			PotentialPRI.push_back(*(&ToaKCenter[vv]));//放置大于门限的PRI值
		}
	}
	delete[] ToaKCenter;//释放内存
	delete[] Bk;
	delete[] Flag;
	delete[] PhaseTimeBegin;

	return PotentialPRI;
}



//功能：根据PRI信息进行雷达信号检索，筛选出具体每部辐射源信号
//输入：PotentialPRI，分选出的PRI信号，ToaToatal，多脉冲数据流
//输出：UniqueToa，独立辐射源到达时间
void SearchRadarSignal(const vector<double> PotentialPRI, const vector<double> ToaToatal, vector<double> UniqueToa)
{
	int Radar_PRINum = PotentialPRI.size();//分选出来的PRI总数目
	int Radar_ToaNum = ToaToatal.size();//雷达全脉冲流TOA总数目
	int GataNum = 0;//判定PRI有效门限
	int Sig_Total = 0;//辐射源信号定位数
	double * Data = new double[Radar_ToaNum];//用于存放各个辐射源信号
	if (Radar_PRINum == 0)
	{
		//PRI分选无结果
		cout << "脉冲流无PRI信号，或分选失败" << endl;
	}
	else
	{
		//PRI分选有结果，进行脉冲序列检索，也就是搜索出各辐射源信号
		//逐个对找到的PRI进行带入搜索脉冲序列，本次工程设计未规划“抖动脉冲”和“脉冲丢失”情况
		for (int ii = 0; ii < Radar_PRINum; ii++)
		{
			//遍历每一个脉冲
			for (int jj = 0; jj < Radar_ToaNum; jj++)
			{
				double Tem_Toa = PotentialPRI.at(ii) + ToaToatal.at(jj);//从Toa的第0位开始加上之前计算出的PRI
				for (int zz = jj + 1; zz < Radar_ToaNum; zz++)//从下一位开始查找是否有等于上面算出的Tem_Toa
				{
					if (abs(Tem_Toa - ToaToatal.at(zz)) < 1e-2)//Tem_Toa = Toatoatal.at(zz)时候
					{
						GataNum = GataNum + 1;
						break;
					}
				}
				if (GataNum > 5) //说明这个PRI确实有效，可能无法应对脉冲丢失和脉冲抖动
				{
					//从中间向前检索
					for (int dd = jj; dd > 0; dd--)
					{
						Tem_Toa = ToaToatal.at(dd) - PotentialPRI.at(ii);
						//检索到了第一个或者是第二个值,也就是0
						/*if (abs(Tem_Toa) < 1e-2)
						{
							Data[Sig_Total] = ToaToatal.at(dd);
							cout << "检索结果辐射源信号TOA：" << Data[Sig_Total] << endl;
							Sig_Total = Sig_Total + 1;
							break;
						}*/
						//向前检索时，依次向前判断
						for (int ww = dd - 1; ww > -1; ww--)
						{
							if (abs(Tem_Toa - ToaToatal.at(ww)) < 1e-2)
							{
								Data[Sig_Total] = ToaToatal.at(ww);
								//cout << "检索结果辐射源信号TOA：" << Data[Sig_Total] << endl;
								Sig_Total = Sig_Total + 1;
								break;
							}
						}
					}
					//向下取整
					int CenterData = trunc(Sig_Total / 2);
					//数据前后对调
					for (int ppi = 0; ppi < CenterData; ppi++)
					{
						double Tem_Toa = Data[Sig_Total - ppi - 1];
						Data[Sig_Total - ppi - 1] = Data[ppi];
						Data[ppi] = Tem_Toa;						
					}
					Data[Sig_Total] = ToaToatal.at(jj);
					Sig_Total = Sig_Total + 1;
					//向后检索
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

					//调试输出
					for (int ll = 0; ll < Sig_Total; ll++)
					{
						cout << "检索结果辐射源信号TOA：" << Data[ll] << endl;
					}
					int aa = 1;


					//因模拟产生脉冲信号TOA叠加和PRI分选结果存在误差，导致脉冲检索结果不准确，故采用前后比对再检索，找出符合本PRI的雷达信号
					double DifferenceDataA_B;
					for (int Diff = 0; Diff < Sig_Total; Diff++)
					{
						DifferenceDataA_B = Data[Diff + 1] - Data[Diff];//检索出来的数据后面一个减去前面一个
						if (DifferenceDataA_B > PotentialPRI.at(ii)*1.03)
						{
							//进入此判断内，说明前后数据有问题，大于130%，说明后面数据错误，或者是丢失脉冲！
							//根据脉冲丢失概率，判断差值是否是当前PRI的整数倍，若是整数倍，说明检索时脉冲丢失，若不是整数倍，说明检索出现差错
							int Flag_Error = 0;
							for (int zfii = 2; zfii < 6; zfii++)//暂定2倍到5倍检测
							{
								if (abs(PotentialPRI.at(ii) * zfii - DifferenceDataA_B) < 1e-2)//0.01
								{
									cout << "检索结果中前后差值为当前PRI的" << zfii << "倍，中间丢失脉冲数为" << zfii - 1 << endl;
									//整数倍数则不删除，应该再去原始数据中查找该区域的辐射源数据，确定是否是错误检索，这里不再执行其他步骤！
									Flag_Error = 1;
								}
							}
							if (Flag_Error == 1)
							{
								//说明属于脉冲丢失，不是错误
							}
							else
							{
								//说明不是脉冲丢失，是检索出错了，则抛弃错误脉冲
								Remove_ArrayElement(Data, &Sig_Total, Diff + 1);//删除之后得到更新的数组
								//Sig_Total是该辐射源总数，Diff+1是要删除的下表，因数据从零开始，下标要比总数少1
								if (Sig_Total > Diff + 1)
								{
									Diff = Diff - 1;//删除之后再重新返回判断一次
								}
							}
						}
						if (DifferenceDataA_B < PotentialPRI.at(ii) * 0.97)
						{
							//进入此判断内，说明前后数据有问题，小于70%，说明肯定错误，之前的检索出错了！
							//需要进一步判断到底是哪个数据有问题
							double DataDiff_DiffMinus = Data[Diff] - Data[Diff - 1];
							double DataDiffAdd_DiffMinus = Data[Diff + 1] - Data[Diff - 1];
							if (abs(DataDiff_DiffMinus - PotentialPRI.at(ii)) < abs(DataDiffAdd_DiffMinus - PotentialPRI.at(ii)))//Data[Diff]与Data[Diff-1]处PRI与检索PRI的差值   小于  （Data[Diff + 1]与Data[Diff - 1]与前前一个PRI与检索PRI的差值），说明前面的数据更接近检索PRI
							{
								//把后面这个差值较大的数据删除掉
								Remove_ArrayElement(Data, &Sig_Total, Diff + 1);//删除之后得到更新的数组
								//Sig_Total是该辐射源总数，Diff+1是要删除的下表，因数据从零开始，下标要比总数少1
								if (Sig_Total > Diff + 1)
								{
									Diff = Diff - 1;//删除之后再重新返回判断一次
								}
							}
							else
							{
								//否是就是前面的数大于等于后面的数，说明后面的数更接近分选的PRI，删除前面的数
								Remove_ArrayElement(Data, &Sig_Total, Diff);//删除之后得到更新的数组
								//Sig_Total是该辐射源总数，Diff+1是要删除的下表，因数据从零开始，下标要比总数少1
								if (Sig_Total > Diff + 1)
								{
									Diff = Diff - 1;//删除之后再重新返回判断一次
								}
							}
						}
					}
					

					//调试输出
					for (int ll = 0; ll < Sig_Total ; ll++)
					{
						cout << "检索结果辐射源信号TOA：" << Data[ll] << endl;
					}
					int aab = 1;
				}

			}
		}
	}
}