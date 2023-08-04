#pragma once
#ifndef __SignalGenerateDefinition_h__
#define __SignalGenerateDefinition_h__

#include <vector>
#include <vector.h>
#include <matrix.h>
#include "PublicDefinition.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "CudaArray.cuh"
#include "CommonKernel.cuh"

using namespace std;
using namespace splab;

//const   double  PI = 3.141592653589793;         //Բ����
//const   double  LIGHTSPEED = 3e8;               //����

typedef Vector<std::complex<double>> DCVec;//˫���Ȳ����ź�����
typedef Vector<std::complex<float>>  SCVec;//�����Ȳ����ź�����
typedef Vector<double>	DVec;          //˫���Ȳ����ź�����
typedef Vector<float>	SVec;          //�����Ȳ����ź�����
typedef Matrix<double>	DMat;          //˫���Ȳ����źž���
typedef Matrix<float>	SMat;          //�����Ȳ����źž���
typedef Matrix<std::complex<double>> DCMat;//˫���Ȳ����źž���
typedef Matrix<std::complex<float>>  SCMat;//�����Ȳ����źž���

typedef std::complex<double>		DCSig;//˫���ȵ�����������ֵ��(����I,Q˫ͨ��)
typedef std::complex<float>			SCSig;//�����ȵ�����������ֵ��(����I,Q˫ͨ��)

// ========================  �źŲ����̶���������ʼ��������  =============================

struct SGSGAPInitParaStruct {
	// ������������
	int JamSamplePointNum;

	//¼ȡ��־λ
	bool SPDataRecord;

	//NLFM��Ƶб��
	vector<double> NLFM_K;

	// ��������
	double NoisePower;

	// ���Ӳ���FFT����
	int ComplexFFTNum;

	// ���LFM FFT����
	int WBLFMFFTNum;

	// ���һά���������
	int HRRPPoint;

	// ԭʼ��Ϣ����
	int OriginalInfoLength;

	// �ٶ�һ����Ϣ����
	int VelocityOnceInfoLength;

	// ���һά��ظ���
	int MaxCorrNum;

	// ���۵�Ԫ����
	int CFARN;

	// �����޼������
	double SlowThresholdDetectFactor;

	// �����޼������
	double FastThresholdDetectFactor;

	// ȫ����������ޣ�fNoise�ı�����
	int	AllJamDetectTh;

	// ȫ��������������Ÿ���
	int DetectAllJamNumMax;

	// �������������ޣ�fNoise�ı�����
	int	MainLobeJamDetectTh;

	// �����������������Ÿ���
	int DetectMainLobeJamNumMax;

	// ��������
	int MaxDetectNum;

	// ���Ŀ�곤��
	int MaxTarLength;

	// ���ջ�������ϵ��(dB)
	double RecNoiseFactor;

	// ��ֵ���书��(W)
	double TransPower;

	// ������ۺ����ֵ(dB)
	double TransLoss;

	// �״��������������������(dB)
	double MainAnteGainMax;

	SGSGAPInitParaStruct() : JamSamplePointNum(32), SPDataRecord(false), NoisePower(1.0),
		ComplexFFTNum(512), WBLFMFFTNum(512), HRRPPoint(536),
		OriginalInfoLength(536), VelocityOnceInfoLength(64),
		MaxCorrNum(20), CFARN(32), SlowThresholdDetectFactor(4),
		FastThresholdDetectFactor(8), AllJamDetectTh(5),
		DetectAllJamNumMax(10), MainLobeJamDetectTh(5),
		DetectMainLobeJamNumMax(10), MaxDetectNum(10), MaxTarLength(30),
		RecNoiseFactor(0.0), TransPower(100e3), TransLoss(6), MainAnteGainMax(35) {
		NLFM_K = {
			-0.11417607723306,
			0.03960138311910,
			-0.02048549632323,
			0.01253307329411,
			-0.00840992355201,
			0.00598620805378
		};
	}
};


// ���ջ������źŴ���ṹ��
struct SGReceiver2SigProcessStruct
{
	// ����������ȡ
	DCMat NoiseSideLobeJamSample;

	// �������������ȡ
	DCMat NoiseMainLobeJamSample;

	// ��ȡĿ��ز����ֵ
	double SrMax;

	// ��������
	int DeletePluseNum;

	// ���ջ���������ƽ
	double fNoise;

	SGReceiver2SigProcessStruct() : SrMax(0.0), DeletePluseNum(0), fNoise(0.0) {}
};


// ========================  ������Ϣ  ==============================
// Ŀ�곡������  ==============================
struct SGTarParaStruct
{
	double			TarRangeReceive;		// Ŀ�����
	double			TarVelocityReceive;		// Ŀ���ٶ�
	DCSig			AlphaErroWeight;		// Ŀ�귽λ�����Ȩֵ
	DCSig			BetaErroWeight;			// Ŀ�긩�������Ȩֵ
	double			TarSNRReceive;			// Ŀ�������
	vector<float>	SNRA;					// ��·����� 
	double			TarSigmacReceive;		// Ŀ��ƽ��RCS
	int				SwerllingType;          // Ŀ��RCS���ģ��

	SGTarParaStruct() : TarRangeReceive(0.0), TarVelocityReceive(0.0), TarSNRReceive(0.0), TarSigmacReceive(0.0), SwerllingType(0) {}
};
struct SGTargetInsightStruct
{
	vector<SGTarParaStruct>		TarPara;
	double						ArrayPara;   // ��������

	SGTargetInsightStruct() : ArrayPara(0.0) {}
};
// ���ų�����Ϣ  ================================
// ѹ�Ƹ��ŷ�ʽ
enum SGJamMannEnum
{
	SGNoiseAM,		//��������
	SGNoisePM,		//��������
	SGNoiseFM,		//������Ƶ
	SGTriggerNoise,	//ת��ʽ
};
// ѹ�Ƹ��Ų���
struct SGNoiseParaStruct
{
	SGJamMannEnum	JamMann;			// ���ŷ�ʽ�������������������ࣻ������Ƶ��ת��ʽ
	double			Fjam;				// ��������Ƶ��
	double			ReceivTime;			// ֡��ʱ��
	double			ReceivTransTime;	// ת������
	double			JamBegRange;		// ���ſ�ʼ���룺����ǰ�صľ���
	double			JamBandwidth;		// ���Ŵ���
	
	SGNoiseParaStruct() : Fjam(0.0), ReceivTime(0.0), ReceivTransTime(0.0), JamBegRange(0.0), JamBandwidth(0.0) {}
};
//ѹ�Ƹ���
struct SGNoiseJamParaStruct
{
	SGNoiseParaStruct	NoisePara;			// ѹ�Ƹ��Ų���
	double				JNRHe;				// ��·�����
	DVec				JNRA;				// ��·�����
	DCSig				AlphaErroWeightJam;	// ��λ��ͨ������Ȩֵ����Ⱥ�·��
	DCSig				BetaErroWeightJam;	// ������ͨ������Ȩֵ
	int					SpaceFeat;			// �������� 0-������ţ�1-�԰����

	SGNoiseJamParaStruct() : JNRHe(0.0), SpaceFeat(0) {}
};
//��Ŀ����Ų���
struct SGCheatTarParaStruct
{
	double Range;		// ��Ŀ�����
	double Velocity;	// ��Ŀ�꾶��

	SGCheatTarParaStruct() : Range(0.0), Velocity(0.0) {}
};

//��ƭ����
struct SGCheatJamParaStruct
{
	vector<SGCheatTarParaStruct>	CheatTarPara;		// ��Ŀ�����
	double							JNRHe;				// ��·�����
	vector<float>					JNRA;				// ��·�����
	DCSig							AlphaErroWeightJam;	// ��1ͨ������ϵ������Ⱥ�·��
	DCSig							BetaErroWeightJam;	// ��2ͨ������ϵ��
	int								SpaceFeat;			// �������� 0-������ţ�1-�԰����

	SGCheatJamParaStruct() : JNRHe(0.0), SpaceFeat(0) {}
};

// ʵװ���Ų���
/*
	ʾ��1����ͨ����������D:�ɼ����� N:������ |:פ����ʼ

	ͨ��1	N N	N N | D D D D
	ͨ��2		  N	| D D D D D D D
	ͨ��3			|   D D D D D D D D
	ͨ��4		N N	| D D D D D D
					  1 2 3 4 5 6 7 8 9
	ʵ�ʲɼ����ݳ��ȣ�JamSampleLen{9��9��9��9}
	��ͨ�������������ȣ�JamFillLen{4��1��0��2}
	��ͨ��ͷ��������Ҫ����ȡ�ĵ�����JamCutBeg{0��0��1��0}
	��ͨ��ǰ�����ݽ�ȡ���ȣ�JamCutLen{4��7��8��6}


	ʾ��2�����崮��������- 1 1 -:����һ��PRT  |:פ����ʼ  D:�ɼ�����  N:������

	�������崮�ɼ���			- 1 1 -	- 1 | 1 - - 1 1 - - 1 1 - -
	��ͨPDʵ�ʲɼ���			N N N N N N | D D D D D D D D D D D
	�׸�Prt������ʼ��λ�ã�PrtSampleStartLen 1
	ÿ������Ĳ���������PrtSampleLen 2


	ʾ��3��ISAR+���崮SAR���崮������

	�������崮�ɼ���			- 1 1 -	- 1 | 1 - - 1 1 - - 1 1 - -
	ISAR+���崮SARʵ�ʲɼ�		- N N - - N | N - -	D D	- -	D D - -
	ʵ�ʲɼ����ݳ��ȣ�JamSampleLen 2 * 2
	��ͨ�������������ȣ�JamFillLen 2 * 2
	�׸�Prt������ʼ��λ�ã�PrtSampleStartLen 1 ��δʹ�ã�
	ÿ������Ĳ���������PrtSampleLen 2
*/
struct SGTrueJamParaStruct
{
	int				ChannelState;       // ͨ��״̬
	unsigned short	JamID;				// ���Ż�ID
	int				JamSampleLen;		// �ɼ������źų���
	int				JamCutLen;			// �����źŽ�ȡ����
	int				JamCutBeg;			// ��ͨ��ͷ��������Ҫ����ȡ�ĵ��� 
	int				JamFillLen;			// �����źŲ�����������
	short*			BegAddress;			// �����ź���ʼ��ַ
	double			ActualJNR;			// ʵ�ʸ����(dB)
	vector<double>	ExpJNRHCC;			// ÿ���ź���������ȡ��Ͳ��(dB)
	vector<double>  ExpJNRASS;			// ÿ���ź����������-����ͨ��(dB)
	double			ActualNoisePower;   // ͨ��ʵ���������ʣ�dBm��
	int				PrtSampleStartLen;	// �׸�Prt������ʼ��λ��
	int				PrtSampleLen;		// ÿ������Ĳ�������Nwid

	SGTrueJamParaStruct() : ChannelState(0), JamID(0), JamSampleLen(0), JamCutLen(0), JamCutBeg(0), JamFillLen(0), BegAddress(0), ActualJNR(0.0), ActualNoisePower(0.0), PrtSampleStartLen(0), PrtSampleLen(0) {}
};

// ���ų�����Ϣ���ܣ�
struct SGJamInsightStruct
{
	vector<SGNoiseJamParaStruct> NoiseJamPara;
	vector<SGCheatJamParaStruct> CheatJamPara;
	vector<SGTrueJamParaStruct> TrueJamPara;
};

// ��������  ===================================
struct SGMissileParaStruct
{
	double	MisRangeReceive;	// ��������
	DCSig	MisAlphaErroWeight;	// ������λ�����Ȩֵ
	DCSig	MisBetaErroWeight;	// �������������Ȩֵ
	double	MisSNRReceive;		// ���������

	SGMissileParaStruct() : MisRangeReceive(0.0), MisAlphaErroWeight(0.0, 0.0), MisBetaErroWeight(0.0, 0.0), MisSNRReceive(0.0) {}
};
struct SGMissileInsightStruct
{
	int					MisID;				// ����ID
	int					DispChannNum;		// ͨ����
	int					FreNum;				//	Ƶ��
	double				MisTrackRCenter;	// �������벨������
	SGMissileParaStruct MissilePara;

	SGMissileInsightStruct() : MisID(0), DispChannNum(0), FreNum(0), MisTrackRCenter(0.0) {}
};
// �Ӳ�����ö����
enum SGClutterPowerSpectrumEnum
{
	SGGauss,			// ��˹��
	SGExponential,		// ָ����
	SGCauchy,			// ������
	SGNthPower,			// N�η���
	SGFullSpectrum,		// ȫ����
};
// �Ӳ�����  =================================
struct SGClutterParaStruct
{
	// �Ӳ�����ͳ�Ʋ���
	// �����ײ���
	SGClutterPowerSpectrumEnum	PowerSpectrumDistrib;	// �Ӳ������׷ֲ�
	int							SpectrumPara;           // N�η���ȫ���ף������׷ֲ������в��� ����Χ 2 ~ 5��
	// �����ײ���
	int							AmpDistrib;				// ���ȷֲ� 0-�����ֲ���1-������̬�ֲ���2-�������ֲ���3-K�ֲ�
	double						ShapeP;					// ��״���� (���������ȷֲ���ʽ���в���)
	double						ScalQ;					// �߶Ȳ���

	// �źŲ���
	double						ClutterRangeBeg;		// �Ӳ�������ʼ
	double						ClutterRangeEnd;		// �Ӳ�������ֹ
	double						ClutterVelocity;		// �Ӳ�����
	DCSig						AlphaErroWeight;		// �Ӳ���λ��Ȩֵ
	DCSig						BetaErroWeight;			// �Ӳ�������Ȩֵ
	double						CNRMax;					// �������� 
	
	SGClutterParaStruct() : SpectrumPara(0), AmpDistrib(0), ShapeP(0.0), ScalQ(0.0), ClutterRangeBeg(0.0), ClutterRangeEnd(0.0), ClutterVelocity(0.0), CNRMax(0.0) {}
};
// �Ӳ�����
struct SGClutterInsightStruct
{
	vector<SGClutterParaStruct> ClutterPara;	// �Ӳ�����
};

// �źŴ������ ===================================
// �������
struct SGFillZerosParaStruct
{
	int		FillZerosFlag;				// �Ƿ������־
	double	FillZerosRange;				// �������
	double	FillZerosRangeBeforTrans;	// ����ǰ�������
	double	FillZerosRangeAfterTrans;	// ������������
	int		DeletePulseNum;				// ���������
	double	VelocityTh;					// �ٶ���������

	SGFillZerosParaStruct() : FillZerosFlag(1), FillZerosRange(0.0), FillZerosRangeBeforTrans(0.0), FillZerosRangeAfterTrans(0.0), DeletePulseNum(0), VelocityTh(0.0){}
};

// ����������
enum SGWinTypeEnum
{
	SGNonWinType = 0,
	SGHamming,
	SGHanning,
	SGChebywin,
};

// ���������Ʋ���
struct SGWinCtrlStruct
{
	SGWinTypeEnum		WinTypeR;
	SGWinTypeEnum		WinTypeV;

	SGWinCtrlStruct() : WinTypeR(SGHamming), WinTypeV(SGChebywin) {}
};

// ������
struct SGDetectStruct
{
	int ThType; //������ͣ�0 - �����ޣ�1 - ������
	
	SGDetectStruct() : ThType(0) {}
};

// ��ѹǰ�������
struct SGBeforePCParaStruct
{
	SGFillZerosParaStruct	FillZerosPara;		// �������
	bool					EchoplexFlag;		// ������Ϣ�����֣�false-ԭʼ�ز���Ϣ true-һ����Ϣ
	bool					AntiMainLobeFlag;	// �����꿪����־λ
	bool					SLBAfterPCFlag;		// ��ѹ���԰�����
	bool					ADBFFlag;			// ADBF������־λ
	Vector<int>				AChannelFlag;		// ��ͨ�����ر�־
	Vector<int>				SLCFlag;			// �������ر�־
	Vector<int>				SLBFlag;			// �������ر�־
	int						OpenedANum;			// �ѿ����ĸ�������

	SGBeforePCParaStruct() : EchoplexFlag(true), AntiMainLobeFlag(true), SLBAfterPCFlag(true), ADBFFlag(true), OpenedANum(0) {}
};

struct SGTrackGateStruct
{
	int		BatchNum;				// ����Ŀ������
	double	TrackGateRCenter;		// ���벨������
	double	TrackGateRSize;			// ���벨�ſ��
	double	TrackGateVCenter;		// �ٶȲ�������
	double	TrackGateVSize;			// �ٶȲ��ſ��

	SGTrackGateStruct() : BatchNum(0), TrackGateRCenter(0.0), TrackGateRSize(0.0), TrackGateVCenter(0.0), TrackGateVSize(0.0) {}
};

// ����㼣����ʱ��������ݽṹ 
struct SGTargetMesStr
{
	int		RgateS;			// �������������ʼ 
	int		RgateE;			/* ����������Ž��� */
	int		VgateS;			/* ������ٶ�����ʼ */
	int		VgateE;			/* ������ٶ��Ž��� */
	float	Amax;			/* ������ź����ֵ */
	int		AmaxRgate;		/* ������ź����ֵ��Ӧ�ľ����� */
	int		AmaxVgate;		/* ������ź����ֵ��Ӧ���ٶ��� */
	bool	Flag;			/* �������־��һά���ʱ��Flag=1����ά���ʱ�������һά�����������Flag=0 */
	bool	VgateFlag;		/* �ٶ����������־ */
	bool	TarValidFlag;	/* ������Ч��־*/

	SGTargetMesStr() : RgateS(0), RgateE(0), VgateS(0), VgateE(0), Amax(0.0), AmaxRgate(0), AmaxVgate(0), Flag(true), VgateFlag(false), TarValidFlag(false) {}
};

// ���ٲ�������
struct SGMeasureStruct
{
	int							DispChannNum;	// Ŀ��ͨ��
	int							TrackGateNum;	// ���Ÿ���
	vector<SGTrackGateStruct>	TrackGate;		// ���Ų���

	SGMeasureStruct() : DispChannNum(0), TrackGateNum(0) {}
};
// �źŴ���������ܣ�
struct SGSignalProcessStruct
{
	SGBeforePCParaStruct	BeforePCPara;	// ��ѹǰ�������
	SGWinCtrlStruct			WinCtrl;		// ���������Ʋ���
	SGDetectStruct			Detect;			// ������
	SGMeasureStruct			Measure;		// ��������
};

// ������Ϣ���ܣ�
struct SGInsightStruct
{
	//����Ŀ�����
	SGTargetInsightStruct TargetInsight;
	//������������
	SGMissileInsightStruct MissileInsight;
	// ���Ų���
	SGJamInsightStruct JamInsight;
	//�Ӳ�����
	SGClutterInsightStruct ClutterInsight;
};

// ===========================  ���Ʋ���  ==============================
// ״̬����  ==================================
//�״��¼�����(��Ϊ��Ҫ�ı�ʶ��ִ�к�������Ĺؼ�����)
//[ע��]:���ﰴ�����ȼ�˳������������Խ��к������ۺ����ȼ�����
enum SGEventTypeEnum
{
	SGSearch,		//����
	SGTrack,		//����
	SGPassive,		//��������
	SGHRRP,			//����
	SGCatch,    	//�ػ�
	SGMissile,		//����
	SGLead,			//��ʾ
	SGError,		//����״̬
};

struct SGStateStruct
{
	SGEventTypeEnum		Mode;		// ģʽ
	int					DwellNum;	// פ������
	bool				DEBUGON;	// ���Ա�־λ
	float				DwellTime;	//פ��ʱ��
	
	SGStateStruct() : DwellNum(0), DEBUGON(false), DwellTime(0.0) {}
};

// ���β���  ==================================
// ��������
enum SGWaveTypeEnum
{
	SGErrorWaveType = -1,	//δָ���ź���ʽ
	SGSingleFrequency,		//����Ƶ
	SGLFM,					//���Ե�Ƶ
	SGNLFM,					//�����Ե�Ƶ
	SGPCM,					//��λ����
	SGMTI,					//MTI
	SGPD,					//���������
	SGComplexWave,			//���Ӳ���
	SGWBLFM,				//������Ե�Ƶȥб
};
// ���Ӳ��β���
struct SGModulOuterPulseStruct
{
	int				StepMode;			// �ݱ䷽ʽ 0-�Զ���1-������2-�����
	double			StepFreInterval;	// ��Ƶ���
	int				StepNum;			// ��Ƶ�������������壩
	int				TeamNum;			// ������

	int				SubPulseNum;		// ���������
	double			SubStepFreInterval; // ��������Ƶ���

	vector<float>	FreGroup;			// ����Ƶ�ʼ��ϣ� ���ڽݱ�Ƶ��ת��ʽ�����ź�����

	SGModulOuterPulseStruct() : StepMode(0), StepFreInterval(0.0), StepNum(0), TeamNum(0), SubPulseNum(0), SubStepFreInterval(0.0) {}
};
//// ����Ƶ����
//struct StepFreParaStruct
//{
//	double StepFreInterval; // ��Ƶ���
//	int StepNum; // ��Ƶ����
//	int TeamNum;// ������
//};
////���������
//struct SubPulseParaStruct
//{
//	int SubPulseModulType;//��������Ʒ�ʽ��0-����Ƶ��1-LFM��2-PCM
//	double SubPulseFre;//��������Ƶ����λMHz���������Ƶ��
//	double SubPulsePulseWidth;//������������λ��us
//	double SubPulseBandWidth;//�����������λ��MHz
//	splab::Vector<double> SubPulsePCMCode;//���������ڵ��Ʊ���
//};
//// ���Ӳ����������
//struct ModulInterPulseStruct
//{
//	double PulseCenterFre;//��������Ƶ��
//	int SubPulseNum; //���������
//	vector<SubPulseParaStruct> SubPulsePara;//���������
//};

// �򵥲��β���
//���ڵ��Ʒ�ʽ
enum SGModulTypeEnum
{
	SGNonModulType = -1,			// δ֪���ڵ��Ʒ�ʽ
	SGModulLFM,						// ���Ե�Ƶ
	SGModulComplexCode,				// ���Ӳ���
	SGModulSingleFre,				// ����Ƶ
	SGModulPCM,						// ��λ����
};

struct SGWaveStruct
{
	SGWaveTypeEnum			WaveType;			// ��������
	SGModulTypeEnum			ModulType;			// ���ڵ��Ʒ�ʽ
	double					PulseWidth;			// ����
	double					fc;					// ��Ƶ
	double					PRT;				// �����ظ�����
	int						PulseNum;			// �������
	DVec					MTIPara;			// MTI�����ظ�����
	double					BandWidth;			// ���ƴ��� Hz
	double					CodeWidth;			// ��Ԫ��� us
	int						PCMCodeIndex;		// ��λ�������г�������

	SGModulOuterPulseStruct	ModulOuterPulse;	// �ݱ�Ƶ����
	
	SGWaveStruct() : PulseWidth(0.0), fc(0.0), PRT(0.0), PulseNum(0), BandWidth(0.0), CodeWidth(0.0), PCMCodeIndex(0.0) {}
};

// ���ջ�����  ===========================
// ���䲨��ָ��
struct SGTransBeamPointingStruct
{
	double	TransAlpha;
	double	TransBeta;

	SGTransBeamPointingStruct() : TransAlpha(0.0), TransBeta(0.0) {}
};

// ���ղ���ָ��
struct SGReceivBeamPointingStruct
{
	double	ReceivAlpha;
	double	ReceivBeta;

	SGReceivBeamPointingStruct() : ReceivAlpha(0.0), ReceivBeta(0.0) {}
};
//��������
struct SGSampleStruct
{
	double	fs;				// ������
	double	SampleBegRange; // ������ʼ���� ��λ:m
	double	SampleRange;	// �������� ��λ:m
	int		SampleBegPoint;	// ������ʼ�㣨�����פ����ʼ��
	int		SamplePoint;	// ��������

	SGSampleStruct() : fs(0.0), SampleBegRange(0.0), SampleRange(0.0), SampleBegPoint(0), SamplePoint(0) {}
};
//���ջ����Ͳ���
struct SGReceivTypeStruct
{
	int		Constitute;	// 0-ֱ�� 1-ȥб
	double	ReceivB;	// ���ջ�����

	SGReceivTypeStruct() : Constitute(0), ReceivB(0.0) {}
};
//���ջ�����
struct SGReceiveStruct
{
	splab::Vector<int>			AChannelFlag;			// ��ͨ�����ر�־
	int							OpenedANum;				// �ѿ����ĸ�������
	SGTransBeamPointingStruct	TransBeamPointing;		// ���䲨��ָ��
	SGReceivBeamPointingStruct	ReceivBeamPointing;		// ���ܲ���ָ��
	double						ReceiverVelocityCompen;	// ���ջ������ٶ�
	double						ReferenceRange;			// �ο�����
	SGReceivTypeStruct			ReceivType;
	SGSampleStruct				Sample;					// ����
	double						AGC;					// �Զ��������

	SGReceiveStruct() : OpenedANum(0), ReceiverVelocityCompen(0.0), ReferenceRange(0.0), AGC(0.0) {}
};

struct SGCWStruct
{
	//״̬����
	SGStateStruct State;
	//���β���
	SGWaveStruct Wave;
	//���ջ�����
	SGReceiveStruct Receive;
	//�źŴ������
	SGSignalProcessStruct SignalProcess;
};


//// ==========================  ������Ϣ  =============================
//
//// Ŀ�������Ϣ  ============================
//struct TargetRWStruct
//{
//	int BatchNum;				// Ŀ������
//	double TarRange;			// Ŀ�����
//	double TarVelocity;			// Ŀ���ٶ�
//	double TarAAngle;			// Ŀ�귽λ�����
//	double TarEAngle;			// Ŀ�긩�������
//	double TarSNR;				// Ŀ�������
//	double TarQuality;			//�ز�����
//};
//
//// ������������Ϣ  =============================
//struct SufJamStruct
//{
//	// �����ܸ������
//	int JamNum;        //���Ÿ���
//	int SideJamNum;    //�԰���Ÿ���
//	int MainJamNum;    //������Ÿ���
//};
//
//// ���������������Ϣ
//
////���Ż�����
//enum SPJamTypeEnum
//{
//	SPNoiseJam,			//���Ƹ���
//	SPCheatJam,			//��ƭ����
//	SPNoiseAndCheatJam,	//ѹ��+��ƭ����
//	SPNonJam,			//�޸���
//};
//
//struct MainJamFeatStruct
//{
//	SPJamTypeEnum Type; // ����������ͣ�0-ѹ�ƣ�1-��ƭ��2-ѹ�Ƽ���ƭ
//	double Bandwidth;   // ���Ŵ���
//	double DutyRate;    // �������ռ�ձ�
//	double JNR;         // ��������
//};
//
//// ��������
//struct JamFeatStruct
//{
//	SPJamTypeEnum JamType;	// ��������
//	double DutyRate;    // ����ռ�ձ�
//	double Bandwidth;   // ���Ŵ���
//	double JNR;         // �����
//	MainJamFeatStruct MainJamFeat; // �����������
//};
//
//// ������������Ϣ����
//struct JamDetectParaRWStruct
//{
//	SufJamStruct SufJam;    //�����ܸ������
//	JamFeatStruct JamFeat;  //��/�����ϸ�������
//	double JamBegRange;
//};
//
//// �������ٻ���  ====================================
//struct PassiveTrackRWStruct
//{
//	double JNR;//�����
//	double AlphaJamAngle; //���ŷ�λ��
//	double BetaJamAngle; //���Ÿ�����
//	int JamQuality;//�������� 0:<12.5; 1:12.5-18; 2:18-30; 3:>30
//};
//
//// ����������Ϣ
//struct MissileRWStruct {
//	int ID;					// ����ID
//	int FreNum;				// ����Ƶ��
//	double RangeErro;		// �������
//	double AAngleErro;		// ��λ�����
//	double EAngleErro;		// ���������
//	double SNR;				// �����
//};
//
//// ������Ϣ���ܣ�
//struct ReturnWordStruct // ������Ϣ����
//{
//	StateStruct State;
//	int TarNum; // Ŀ�����
//	vector<TargetRWStruct> TargetRW;//Ŀ�����
//	JamDetectParaRWStruct JamDetectPara;//����������
//	PassiveTrackRWStruct PassiveTrackRW;//��������
//	MissileRWStruct MissileRW;//����
//};
//
//// һ����Ϣ����  =============================
//////��������һ����Ϣ����
////struct PassiveOnceInfoStruct
////{
////	double JamRangeOnceInfo[536] = { 0 };//���ž���һ����Ϣ
////	double JamVelocityOnceInfo[64] = { 0 };//�����ٶ�һ����Ϣ
////};
//
//// Ŀ�����һ����Ϣ���Ͳ���
//// ���ٲ��β���
//struct TrackWaveStruct
//{
//	int DispChannNum;			// Ŀ��ͨ����
//	int BatchNum;				// Ŀ������
//	SPWaveTypeEnum WaveType;	// ��������
//	double BandWidth;			// ���ƴ���
//	double PulseWidth;			// ����
//	double PRT;					// �����ظ�����
//	int FreNum;					// Ƶ���
//};
//
//// Ŀ����ٲ���һ����Ϣ���Ͳ���
//struct MeasureParaStruct
//{
//	double TarRangeErro;				// Ŀ��������
//	double TarVelocityErro;				// Ŀ���ٶ����
//	double TarAAngleErro;				// Ŀ�귽λ�����
//	double TarEAngleErro;				// Ŀ�긩�������
//	double TarSNR;						// Ŀ�������
//	double TarAmp;						// Ŀ�����
//};
//
//// һ����Ϣ����
//struct OnceInfoStruct
//{
//	TrackWaveStruct TrackWave;
//	MeasureParaStruct MeasurePara;
//
//	int TrackGateRCenter;					//���벨������
//	int TrackGateRSize;						//���벨�Ŵ�С
//	int TrackGateVCenter;					//�ٶȲ�������
//	int TrackGateVSize;						//�ٶȲ��Ŵ�С
//	double RangeOnceInfo[536] = { 0 };		//����һ����Ϣ����
//	double VelocityOnceInfo[64] = { 0 };	//�ٶ�һ����Ϣ����
//};


// �ṹ�����He��C1��C2��A(�ศ·)�źž���
struct SGChannelSignalMatrix
{
	DCMat He;
	DCMat C1;
	Matrix<std::complex<double>> C2;
	Vector<Matrix<std::complex<double>>> A;

	SGChannelSignalMatrix()
	{
		memset(this, 0, sizeof(SGChannelSignalMatrix));
	}
};

//// �������
//struct MeasureOutStruct
//{
//	double Range;		// ����Ŀ�����
//	double Velocity;	// ����Ŀ���ٶ�
//	double AlphaAngle;	// ����Ŀ�귽λ���鲿
//	double BetaAngle;	// ����Ŀ�긩�����鲿
//	double Amp;			// ����Ŀ�����
//	double SNR;			// ����Ŀ��SNR
//	int Quality;		// ����Ŀ������
//};



/*���׶εĴ���ʱ����Ϣ*/
struct SGTimeInfo
{
	float	TarTimeInfo;				// Ŀ�����ɺ�ʱ
	float	ClutterTimeInfo;			// �Ӳ����ɺ�ʱ
	float	JamTimeInfo;				// �������ɺ�ʱ
	float	ReceiverTimeInfo;			// ���ջ��ź����ɺ�ʱ
	float	BeforePCTimeInfo;			// ��ѹǰ�����ʱ
	float	PCTimeInfo;					// ��ѹ��ʱ
	float	TargetDetectTimeInfo;		// Ŀ�����ʱ
	float	ParaMeasureTimeInfo;		// ����������ʱ
	float	SignalProcessTimeInfo;		// �źŴ����ʱ
	float	AllProcessTimeInfo;			// ȫ����ʱ��

	SGTimeInfo():TarTimeInfo(0.0), ClutterTimeInfo(0.0), JamTimeInfo(0.0), ReceiverTimeInfo(0.0), BeforePCTimeInfo(0.0), PCTimeInfo(0.0), TargetDetectTimeInfo(0.0), ParaMeasureTimeInfo(0.0), SignalProcessTimeInfo(0.0), AllProcessTimeInfo(0.0)
	{}

};

/*�ڴ���Ϣ*/
struct SGMemInfo
{
	size_t TarMemInfo;					// Ŀ�������ڴ�
	size_t ClutterMemInfo;				// �Ӳ������ڴ�
	size_t JamMemInfo;					// ���������ڴ�
	size_t ReceiverMemInfo;				// ���ջ��ź������ڴ�
	size_t BeforePCMemInfo;				// ��ѹǰ�����ڴ�
	size_t PCMemInfo;					// ��ѹ�ڴ�
	size_t TargetDetectMemInfo;			// Ŀ�����ڴ�
	size_t ParaMeasureMemInfo;			// ���������ڴ�
	size_t SignalProcessMemInfo;		// �źŴ����ڴ�
	size_t AllProcessMemInfo;			// ȫ�����ڴ�

	size_t TarBufInfo;					// Ŀ�����ɻ���
	size_t ClutterBufInfo;				// �Ӳ����ɻ���
	size_t JamBufInfo;					// �������ɻ���
	size_t ReceiverBufInfo;				// ���ջ��ź����ɻ���
	size_t BeforePCBufInfo;				// ��ѹǰ������
	size_t PCBufInfo;					// ��ѹ����
	size_t TargetDetectBufInfo;			// Ŀ���⻺��
	size_t ParaMeasureBufInfo;			// ������������
	size_t SignalProcessBufInfo;		// �źŴ�����
	size_t AllProcessBufInfo;			// ȫ���̻���


	SGMemInfo():TarMemInfo(0), ClutterMemInfo(0), JamMemInfo(0), ReceiverMemInfo(0), PCMemInfo(0), TargetDetectMemInfo(0), ParaMeasureMemInfo(0), SignalProcessMemInfo(0), AllProcessMemInfo(0),TarBufInfo(0), ClutterBufInfo(0), JamBufInfo(0), ReceiverBufInfo(0), PCBufInfo(0), TargetDetectBufInfo(0), ParaMeasureBufInfo(0), SignalProcessBufInfo(0), AllProcessBufInfo(0)
	{}
};


#endif