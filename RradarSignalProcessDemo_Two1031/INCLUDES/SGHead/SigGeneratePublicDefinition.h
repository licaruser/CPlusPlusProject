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

//const   double  PI = 3.141592653589793;         //圆周率
//const   double  LIGHTSPEED = 3e8;               //光速

typedef Vector<std::complex<double>> DCVec;//双精度采样信号序列
typedef Vector<std::complex<float>>  SCVec;//单精度采样信号序列
typedef Vector<double>	DVec;          //双精度采样信号序列
typedef Vector<float>	SVec;          //单精度采样信号序列
typedef Matrix<double>	DMat;          //双精度采样信号矩阵
typedef Matrix<float>	SMat;          //单精度采样信号矩阵
typedef Matrix<std::complex<double>> DCMat;//双精度采样信号矩阵
typedef Matrix<std::complex<float>>  SCMat;//单精度采样信号矩阵

typedef std::complex<double>		DCSig;//双精度单个复数采样值点(包括I,Q双通道)
typedef std::complex<float>			SCSig;//单精度单个复数采样值点(包括I,Q双通道)

// ========================  信号产生固定参数（初始化参数）  =============================

struct SGSGAPInitParaStruct {
	// 干扰样本点数
	int JamSamplePointNum;

	//录取标志位
	bool SPDataRecord;

	//NLFM调频斜率
	vector<double> NLFM_K;

	// 噪声功率
	double NoisePower;

	// 复杂波形FFT点数
	int ComplexFFTNum;

	// 宽带LFM FFT点数
	int WBLFMFFTNum;

	// 宽带一维距离像点数
	int HRRPPoint;

	// 原始信息长度
	int OriginalInfoLength;

	// 速度一次信息场景
	int VelocityOnceInfoLength;

	// 最大一维相关个数
	int MaxCorrNum;

	// 积累单元个数
	int CFARN;

	// 慢门限检测因子
	double SlowThresholdDetectFactor;

	// 快门限检测因子
	double FastThresholdDetectFactor;

	// 全干扰侦察门限（fNoise的倍数）
	int	AllJamDetectTh;

	// 全干扰侦察最大检测干扰个数
	int DetectAllJamNumMax;

	// 主瓣干扰侦察门限（fNoise的倍数）
	int	MainLobeJamDetectTh;

	// 主瓣干扰侦察最大检测干扰个数
	int DetectMainLobeJamNumMax;

	// 最大检测个数
	int MaxDetectNum;

	// 最大目标长度
	int MaxTarLength;

	// 接收机热噪声系数(dB)
	double RecNoiseFactor;

	// 峰值发射功率(W)
	double TransPower;

	// 发射机综合损耗值(dB)
	double TransLoss;

	// 雷达天线主阵主瓣最大增益(dB)
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


// 接收机传给信号处理结构体
struct SGReceiver2SigProcessStruct
{
	// 干扰样本获取
	DCMat NoiseSideLobeJamSample;

	// 主瓣干扰样本获取
	DCMat NoiseMainLobeJamSample;

	// 获取目标回波最大值
	double SrMax;

	// 舍脉冲数
	int DeletePluseNum;

	// 接收机热噪声电平
	double fNoise;

	SGReceiver2SigProcessStruct() : SrMax(0.0), DeletePluseNum(0), fNoise(0.0) {}
};


// ========================  场景信息  ==============================
// 目标场景参数  ==============================
struct SGTarParaStruct
{
	double			TarRangeReceive;		// 目标距离
	double			TarVelocityReceive;		// 目标速度
	DCSig			AlphaErroWeight;		// 目标方位角误差权值
	DCSig			BetaErroWeight;			// 目标俯仰角误差权值
	double			TarSNRReceive;			// 目标信噪比
	vector<float>	SNRA;					// 辅路信噪比 
	double			TarSigmacReceive;		// 目标平均RCS
	int				SwerllingType;          // 目标RCS起伏模型

	SGTarParaStruct() : TarRangeReceive(0.0), TarVelocityReceive(0.0), TarSNRReceive(0.0), TarSigmacReceive(0.0), SwerllingType(0) {}
};
struct SGTargetInsightStruct
{
	vector<SGTarParaStruct>		TarPara;
	double						ArrayPara;   // 阵列因子

	SGTargetInsightStruct() : ArrayPara(0.0) {}
};
// 干扰场景信息  ================================
// 压制干扰方式
enum SGJamMannEnum
{
	SGNoiseAM,		//噪声调幅
	SGNoisePM,		//噪声调相
	SGNoiseFM,		//噪声调频
	SGTriggerNoise,	//转发式
};
// 压制干扰参数
struct SGNoiseParaStruct
{
	SGJamMannEnum	JamMann;			// 干扰方式：噪声调幅；噪声调相；噪声调频；转发式
	double			Fjam;				// 干扰中心频率
	double			ReceivTime;			// 帧收时间
	double			ReceivTransTime;	// 转发周期
	double			JamBegRange;		// 干扰开始距离：干扰前沿的距离
	double			JamBandwidth;		// 干扰带宽
	
	SGNoiseParaStruct() : Fjam(0.0), ReceivTime(0.0), ReceivTransTime(0.0), JamBegRange(0.0), JamBandwidth(0.0) {}
};
//压制干扰
struct SGNoiseJamParaStruct
{
	SGNoiseParaStruct	NoisePara;			// 压制干扰参数
	double				JNRHe;				// 和路干噪比
	DVec				JNRA;				// 辅路干噪比
	DCSig				AlphaErroWeightJam;	// 方位差通道干扰权值（相比和路）
	DCSig				BetaErroWeightJam;	// 俯仰差通道干扰权值
	int					SpaceFeat;			// 空域特征 0-主瓣干扰，1-旁瓣干扰

	SGNoiseJamParaStruct() : JNRHe(0.0), SpaceFeat(0) {}
};
//假目标干扰参数
struct SGCheatTarParaStruct
{
	double Range;		// 假目标距离
	double Velocity;	// 假目标径速

	SGCheatTarParaStruct() : Range(0.0), Velocity(0.0) {}
};

//欺骗干扰
struct SGCheatJamParaStruct
{
	vector<SGCheatTarParaStruct>	CheatTarPara;		// 假目标参数
	double							JNRHe;				// 和路干噪比
	vector<float>					JNRA;				// 辅路干噪比
	DCSig							AlphaErroWeightJam;	// 差1通道干扰系数（相比和路）
	DCSig							BetaErroWeightJam;	// 差2通道干扰系数
	int								SpaceFeat;			// 空域特征 0-主瓣干扰，1-旁瓣干扰

	SGCheatJamParaStruct() : JNRHe(0.0), SpaceFeat(0) {}
};

// 实装干扰参数
/*
	示例1（各通道参数）：D:采集数据 N:补噪声 |:驻留起始

	通道1	N N	N N | D D D D
	通道2		  N	| D D D D D D D
	通道3			|   D D D D D D D D
	通道4		N N	| D D D D D D
					  1 2 3 4 5 6 7 8 9
	实际采集数据长度：JamSampleLen{9，9，9，9}
	各通道补充噪声长度：JamFillLen{4，1，0，2}
	各通道头部数据需要被截取的点数：JamCutBeg{0，0，1，0}
	各通道前端数据截取长度：JamCutLen{4，7，8，6}


	示例2（脉冲串参数）：- 1 1 -:表征一个PRT  |:驻留起始  D:采集数据  N:补噪声

	理论脉冲串采集：			- 1 1 -	- 1 | 1 - - 1 1 - - 1 1 - -
	普通PD实际采集：			N N N N N N | D D D D D D D D D D D
	首个Prt采样起始点位置：PrtSampleStartLen 1
	每个脉冲的采样点数：PrtSampleLen 2


	示例3（ISAR+脉冲串SAR脉冲串参数）

	理论脉冲串采集：			- 1 1 -	- 1 | 1 - - 1 1 - - 1 1 - -
	ISAR+脉冲串SAR实际采集		- N N - - N | N - -	D D	- -	D D - -
	实际采集数据长度：JamSampleLen 2 * 2
	各通道补充噪声长度：JamFillLen 2 * 2
	首个Prt采样起始点位置：PrtSampleStartLen 1 （未使用）
	每个脉冲的采样点数：PrtSampleLen 2
*/
struct SGTrueJamParaStruct
{
	int				ChannelState;       // 通道状态
	unsigned short	JamID;				// 干扰机ID
	int				JamSampleLen;		// 采集干扰信号长度
	int				JamCutLen;			// 干扰信号截取长度
	int				JamCutBeg;			// 各通道头部数据需要被截取的点数 
	int				JamFillLen;			// 干扰信号补充噪声长度
	short*			BegAddress;			// 干扰信号起始地址
	double			ActualJNR;			// 实际干噪比(dB)
	vector<double>	ExpJNRHCC;			// 每阵信号期望信噪比—和差差(dB)
	vector<double>  ExpJNRASS;			// 每阵信号期望信噪比-辅助通道(dB)
	double			ActualNoisePower;   // 通道实际噪声功率（dBm）
	int				PrtSampleStartLen;	// 首个Prt采样起始点位置
	int				PrtSampleLen;		// 每个脉冲的采样点数Nwid

	SGTrueJamParaStruct() : ChannelState(0), JamID(0), JamSampleLen(0), JamCutLen(0), JamCutBeg(0), JamFillLen(0), BegAddress(0), ActualJNR(0.0), ActualNoisePower(0.0), PrtSampleStartLen(0), PrtSampleLen(0) {}
};

// 干扰场景信息（总）
struct SGJamInsightStruct
{
	vector<SGNoiseJamParaStruct> NoiseJamPara;
	vector<SGCheatJamParaStruct> CheatJamPara;
	vector<SGTrueJamParaStruct> TrueJamPara;
};

// 导弹场景  ===================================
struct SGMissileParaStruct
{
	double	MisRangeReceive;	// 导弹距离
	DCSig	MisAlphaErroWeight;	// 导弹方位角误差权值
	DCSig	MisBetaErroWeight;	// 导弹俯仰角误差权值
	double	MisSNRReceive;		// 导弹信噪比

	SGMissileParaStruct() : MisRangeReceive(0.0), MisAlphaErroWeight(0.0, 0.0), MisBetaErroWeight(0.0, 0.0), MisSNRReceive(0.0) {}
};
struct SGMissileInsightStruct
{
	int					MisID;				// 导弹ID
	int					DispChannNum;		// 通道号
	int					FreNum;				//	频点
	double				MisTrackRCenter;	// 导弹距离波门中心
	SGMissileParaStruct MissilePara;

	SGMissileInsightStruct() : MisID(0), DispChannNum(0), FreNum(0), MisTrackRCenter(0.0) {}
};
// 杂波功率枚举类
enum SGClutterPowerSpectrumEnum
{
	SGGauss,			// 高斯谱
	SGExponential,		// 指数谱
	SGCauchy,			// 柯西谱
	SGNthPower,			// N次方谱
	SGFullSpectrum,		// 全极谱
};
// 杂波场景  =================================
struct SGClutterParaStruct
{
	// 杂波类型统计参数
	// 功率谱参数
	SGClutterPowerSpectrumEnum	PowerSpectrumDistrib;	// 杂波功率谱分布
	int							SpectrumPara;           // N次方与全极谱（功率谱分布）特有参数 （范围 2 ~ 5）
	// 幅度谱参数
	int							AmpDistrib;				// 幅度分布 0-瑞利分布，1-对数正态分布，2-威布尔分布，3-K分布
	double						ShapeP;					// 形状参数 (威布尔幅度分布公式特有参数)
	double						ScalQ;					// 尺度参数

	// 信号参数
	double						ClutterRangeBeg;		// 杂波距离起始
	double						ClutterRangeEnd;		// 杂波距离终止
	double						ClutterVelocity;		// 杂波径速
	DCSig						AlphaErroWeight;		// 杂波方位角权值
	DCSig						BetaErroWeight;			// 杂波俯仰角权值
	double						CNRMax;					// 最大杂噪比 
	
	SGClutterParaStruct() : SpectrumPara(0), AmpDistrib(0), ShapeP(0.0), ScalQ(0.0), ClutterRangeBeg(0.0), ClutterRangeEnd(0.0), ClutterVelocity(0.0), CNRMax(0.0) {}
};
// 杂波场景
struct SGClutterInsightStruct
{
	vector<SGClutterParaStruct> ClutterPara;	// 杂波参数
};

// 信号处理参数 ===================================
// 舍近参数
struct SGFillZerosParaStruct
{
	int		FillZerosFlag;				// 是否填零标志
	double	FillZerosRange;				// 舍近距离
	double	FillZerosRangeBeforTrans;	// 发射前置零距离
	double	FillZerosRangeAfterTrans;	// 发射后置零距离
	int		DeletePulseNum;				// 舍脉冲个数
	double	VelocityTh;					// 速度填零门限

	SGFillZerosParaStruct() : FillZerosFlag(1), FillZerosRange(0.0), FillZerosRangeBeforTrans(0.0), FillZerosRangeAfterTrans(0.0), DeletePulseNum(0), VelocityTh(0.0){}
};

// 窗函数类型
enum SGWinTypeEnum
{
	SGNonWinType = 0,
	SGHamming,
	SGHanning,
	SGChebywin,
};

// 窗函数控制参数
struct SGWinCtrlStruct
{
	SGWinTypeEnum		WinTypeR;
	SGWinTypeEnum		WinTypeV;

	SGWinCtrlStruct() : WinTypeR(SGHamming), WinTypeV(SGChebywin) {}
};

// 检测参数
struct SGDetectStruct
{
	int ThType; //检测类型：0 - 慢门限；1 - 快门限
	
	SGDetectStruct() : ThType(0) {}
};

// 脉压前处理参数
struct SGBeforePCParaStruct
{
	SGFillZerosParaStruct	FillZerosPara;		// 填零参数
	bool					EchoplexFlag;		// 回送信息控制字：false-原始回波信息 true-一次信息
	bool					AntiMainLobeFlag;	// 抗主瓣开启标志位
	bool					SLBAfterPCFlag;		// 脉压后旁瓣消隐
	bool					ADBFFlag;			// ADBF开启标志位
	Vector<int>				AChannelFlag;		// 辅通道开关标志
	Vector<int>				SLCFlag;			// 对消开关标志
	Vector<int>				SLBFlag;			// 消隐开关标志
	int						OpenedANum;			// 已开启的辅助个数

	SGBeforePCParaStruct() : EchoplexFlag(true), AntiMainLobeFlag(true), SLBAfterPCFlag(true), ADBFFlag(true), OpenedANum(0) {}
};

struct SGTrackGateStruct
{
	int		BatchNum;				// 跟踪目标批号
	double	TrackGateRCenter;		// 距离波门中心
	double	TrackGateRSize;			// 距离波门宽度
	double	TrackGateVCenter;		// 速度波门中心
	double	TrackGateVSize;			// 速度波门宽度

	SGTrackGateStruct() : BatchNum(0), TrackGateRCenter(0.0), TrackGateRSize(0.0), TrackGateVCenter(0.0), TrackGateVSize(0.0) {}
};

// 定义点迹搜索时相关区数据结构 
struct SGTargetMesStr
{
	int		RgateS;			// 相关区距离门起始 
	int		RgateE;			/* 相关区距离门结束 */
	int		VgateS;			/* 相关区速度门起始 */
	int		VgateE;			/* 相关区速度门结束 */
	float	Amax;			/* 相关区信号最大值 */
	int		AmaxRgate;		/* 相关区信号最大值对应的距离门 */
	int		AmaxVgate;		/* 相关区信号最大值对应的速度门 */
	bool	Flag;			/* 相关区标志，一维相关时置Flag=1，二维相关时，如果该一维相关区被关联Flag=0 */
	bool	VgateFlag;		/* 速度区相关区标志 */
	bool	TarValidFlag;	/* 包络有效标志*/

	SGTargetMesStr() : RgateS(0), RgateE(0), VgateS(0), VgateE(0), Amax(0.0), AmaxRgate(0), AmaxVgate(0), Flag(true), VgateFlag(false), TarValidFlag(false) {}
};

// 跟踪测量参数
struct SGMeasureStruct
{
	int							DispChannNum;	// 目标通道
	int							TrackGateNum;	// 波门个数
	vector<SGTrackGateStruct>	TrackGate;		// 波门参数

	SGMeasureStruct() : DispChannNum(0), TrackGateNum(0) {}
};
// 信号处理参数（总）
struct SGSignalProcessStruct
{
	SGBeforePCParaStruct	BeforePCPara;	// 脉压前处理参数
	SGWinCtrlStruct			WinCtrl;		// 窗函数控制参数
	SGDetectStruct			Detect;			// 检测参数
	SGMeasureStruct			Measure;		// 测量参数
};

// 场景信息（总）
struct SGInsightStruct
{
	//场景目标参数
	SGTargetInsightStruct TargetInsight;
	//场景导弹参数
	SGMissileInsightStruct MissileInsight;
	// 干扰参数
	SGJamInsightStruct JamInsight;
	//杂波参数
	SGClutterInsightStruct ClutterInsight;
};

// ===========================  控制参数  ==============================
// 状态参数  ==================================
//雷达事件类型(极为重要的标识，执行后续处理的关键依据)
//[注意]:这里按照优先级顺序进行了排序，以进行后续的综合优先级计算
enum SGEventTypeEnum
{
	SGSearch,		//搜索
	SGTrack,		//跟踪
	SGPassive,		//被动跟踪
	SGHRRP,			//成像
	SGCatch,    	//截获
	SGMissile,		//导弹
	SGLead,			//领示
	SGError,		//错误状态
};

struct SGStateStruct
{
	SGEventTypeEnum		Mode;		// 模式
	int					DwellNum;	// 驻留计数
	bool				DEBUGON;	// 调试标志位
	float				DwellTime;	//驻留时长
	
	SGStateStruct() : DwellNum(0), DEBUGON(false), DwellTime(0.0) {}
};

// 波形参数  ==================================
// 波形类型
enum SGWaveTypeEnum
{
	SGErrorWaveType = -1,	//未指定信号形式
	SGSingleFrequency,		//单载频
	SGLFM,					//线性调频
	SGNLFM,					//非线性调频
	SGPCM,					//相位编码
	SGMTI,					//MTI
	SGPD,					//脉冲多普勒
	SGComplexWave,			//复杂波形
	SGWBLFM,				//宽带线性调频去斜
};
// 复杂波形参数
struct SGModulOuterPulseStruct
{
	int				StepMode;			// 捷变方式 0-自动；1-步进；2-随机；
	double			StepFreInterval;	// 跳频间隔
	int				StepNum;			// 跳频脉冲数（大脉冲）
	int				TeamNum;			// 分组数

	int				SubPulseNum;		// 子脉冲个数
	double			SubStepFreInterval; // 子脉冲跳频间隔

	vector<float>	FreGroup;			// 中心频率集合， 用于捷变频的转发式干扰信号生成

	SGModulOuterPulseStruct() : StepMode(0), StepFreInterval(0.0), StepNum(0), TeamNum(0), SubPulseNum(0), SubStepFreInterval(0.0) {}
};
//// 步进频参数
//struct StepFreParaStruct
//{
//	double StepFreInterval; // 跳频间隔
//	int StepNum; // 跳频脉冲
//	int TeamNum;// 分组数
//};
////子脉冲参数
//struct SubPulseParaStruct
//{
//	int SubPulseModulType;//子脉冲调制方式：0-单载频；1-LFM；2-PCM
//	double SubPulseFre;//子脉冲载频，单位MHz，相对中心频率
//	double SubPulsePulseWidth;//子脉冲脉宽，单位：us
//	double SubPulseBandWidth;//子脉冲带宽，单位：MHz
//	splab::Vector<double> SubPulsePCMCode;//子脉冲脉内调制编码
//};
//// 复杂波形脉间参数
//struct ModulInterPulseStruct
//{
//	double PulseCenterFre;//脉冲中心频率
//	int SubPulseNum; //子脉冲个数
//	vector<SubPulseParaStruct> SubPulsePara;//子脉冲参数
//};

// 简单波形参数
//脉内调制方式
enum SGModulTypeEnum
{
	SGNonModulType = -1,			// 未知脉内调制方式
	SGModulLFM,						// 线性调频
	SGModulComplexCode,				// 复杂波形
	SGModulSingleFre,				// 单载频
	SGModulPCM,						// 相位编码
};

struct SGWaveStruct
{
	SGWaveTypeEnum			WaveType;			// 波形类型
	SGModulTypeEnum			ModulType;			// 脉内调制方式
	double					PulseWidth;			// 脉宽
	double					fc;					// 载频
	double					PRT;				// 脉冲重复周期
	int						PulseNum;			// 脉冲个数
	DVec					MTIPara;			// MTI脉冲重复周期
	double					BandWidth;			// 调制带宽 Hz
	double					CodeWidth;			// 码元宽度 us
	int						PCMCodeIndex;		// 相位编码序列长度索引

	SGModulOuterPulseStruct	ModulOuterPulse;	// 捷变频参数
	
	SGWaveStruct() : PulseWidth(0.0), fc(0.0), PRT(0.0), PulseNum(0), BandWidth(0.0), CodeWidth(0.0), PCMCodeIndex(0.0) {}
};

// 接收机参数  ===========================
// 发射波束指向
struct SGTransBeamPointingStruct
{
	double	TransAlpha;
	double	TransBeta;

	SGTransBeamPointingStruct() : TransAlpha(0.0), TransBeta(0.0) {}
};

// 接收波束指向
struct SGReceivBeamPointingStruct
{
	double	ReceivAlpha;
	double	ReceivBeta;

	SGReceivBeamPointingStruct() : ReceivAlpha(0.0), ReceivBeta(0.0) {}
};
//采样参数
struct SGSampleStruct
{
	double	fs;				// 采样率
	double	SampleBegRange; // 采样起始距离 单位:m
	double	SampleRange;	// 采样距离 单位:m
	int		SampleBegPoint;	// 采样起始点（相对于驻留起始）
	int		SamplePoint;	// 采样点数

	SGSampleStruct() : fs(0.0), SampleBegRange(0.0), SampleRange(0.0), SampleBegPoint(0), SamplePoint(0) {}
};
//接收机类型参数
struct SGReceivTypeStruct
{
	int		Constitute;	// 0-直采 1-去斜
	double	ReceivB;	// 接收机带宽

	SGReceivTypeStruct() : Constitute(0), ReceivB(0.0) {}
};
//接收机参数
struct SGReceiveStruct
{
	splab::Vector<int>			AChannelFlag;			// 辅通道开关标志
	int							OpenedANum;				// 已开启的辅助个数
	SGTransBeamPointingStruct	TransBeamPointing;		// 发射波束指向
	SGReceivBeamPointingStruct	ReceivBeamPointing;		// 接受波束指向
	double						ReceiverVelocityCompen;	// 接收机补偿速度
	double						ReferenceRange;			// 参考距离
	SGReceivTypeStruct			ReceivType;
	SGSampleStruct				Sample;					// 采样
	double						AGC;					// 自动增益控制

	SGReceiveStruct() : OpenedANum(0), ReceiverVelocityCompen(0.0), ReferenceRange(0.0), AGC(0.0) {}
};

struct SGCWStruct
{
	//状态参数
	SGStateStruct State;
	//波形参数
	SGWaveStruct Wave;
	//接收机参数
	SGReceiveStruct Receive;
	//信号处理参数
	SGSignalProcessStruct SignalProcess;
};


//// ==========================  回送信息  =============================
//
//// 目标回送信息  ============================
//struct TargetRWStruct
//{
//	int BatchNum;				// 目标批号
//	double TarRange;			// 目标距离
//	double TarVelocity;			// 目标速度
//	double TarAAngle;			// 目标方位角误差
//	double TarEAngle;			// 目标俯仰角误差
//	double TarSNR;				// 目标信噪比
//	double TarQuality;			//回波质量
//};
//
//// 干扰侦查回送信息  =============================
//struct SufJamStruct
//{
//	// 整体受干扰情况
//	int JamNum;        //干扰个数
//	int SideJamNum;    //旁瓣干扰个数
//	int MainJamNum;    //主瓣干扰个数
//};
//
//// 主瓣干扰侦查回送信息
//
////干扰机类型
//enum SPJamTypeEnum
//{
//	SPNoiseJam,			//抑制干扰
//	SPCheatJam,			//欺骗干扰
//	SPNoiseAndCheatJam,	//压制+欺骗干扰
//	SPNonJam,			//无干扰
//};
//
//struct MainJamFeatStruct
//{
//	SPJamTypeEnum Type; // 主瓣干扰类型，0-压制，1-欺骗，2-压制加欺骗
//	double Bandwidth;   // 干扰带宽
//	double DutyRate;    // 主瓣干扰占空比
//	double JNR;         // 主瓣干噪比
//};
//
//// 干扰特征
//struct JamFeatStruct
//{
//	SPJamTypeEnum JamType;	// 干扰类型
//	double DutyRate;    // 干扰占空比
//	double Bandwidth;   // 干扰带宽
//	double JNR;         // 干噪比
//	MainJamFeatStruct MainJamFeat; // 主瓣干扰特征
//};
//
//// 干扰侦查回送信息参数
//struct JamDetectParaRWStruct
//{
//	SufJamStruct SufJam;    //整体受干扰情况
//	JamFeatStruct JamFeat;  //旁/主瓣混合干扰特征
//	double JamBegRange;
//};
//
//// 被动跟踪回送  ====================================
//struct PassiveTrackRWStruct
//{
//	double JNR;//干噪比
//	double AlphaJamAngle; //干扰方位角
//	double BetaJamAngle; //干扰俯仰角
//	int JamQuality;//干扰质量 0:<12.5; 1:12.5-18; 2:18-30; 3:>30
//};
//
//// 导弹回送信息
//struct MissileRWStruct {
//	int ID;					// 导弹ID
//	int FreNum;				// 导弹频点
//	double RangeErro;		// 距离误差
//	double AAngleErro;		// 方位角误差
//	double EAngleErro;		// 俯仰角误差
//	double SNR;				// 信噪比
//};
//
//// 回送信息（总）
//struct ReturnWordStruct // 回送信息参数
//{
//	StateStruct State;
//	int TarNum; // 目标个数
//	vector<TargetRWStruct> TargetRW;//目标回送
//	JamDetectParaRWStruct JamDetectPara;//干扰侦察回送
//	PassiveTrackRWStruct PassiveTrackRW;//被动跟踪
//	MissileRWStruct MissileRW;//导弹
//};
//
//// 一次信息回送  =============================
//////被动跟踪一次信息回送
////struct PassiveOnceInfoStruct
////{
////	double JamRangeOnceInfo[536] = { 0 };//干扰距离一次信息
////	double JamVelocityOnceInfo[64] = { 0 };//干扰速度一次信息
////};
//
//// 目标跟踪一次信息回送参数
//// 跟踪波形参数
//struct TrackWaveStruct
//{
//	int DispChannNum;			// 目标通道号
//	int BatchNum;				// 目标批号
//	SPWaveTypeEnum WaveType;	// 波形类型
//	double BandWidth;			// 调制带宽
//	double PulseWidth;			// 脉宽
//	double PRT;					// 脉冲重复周期
//	int FreNum;					// 频点号
//};
//
//// 目标跟踪测量一次信息回送参数
//struct MeasureParaStruct
//{
//	double TarRangeErro;				// 目标距离误差
//	double TarVelocityErro;				// 目标速度误差
//	double TarAAngleErro;				// 目标方位角误差
//	double TarEAngleErro;				// 目标俯仰角误差
//	double TarSNR;						// 目标信噪比
//	double TarAmp;						// 目标幅度
//};
//
//// 一次信息回送
//struct OnceInfoStruct
//{
//	TrackWaveStruct TrackWave;
//	MeasureParaStruct MeasurePara;
//
//	int TrackGateRCenter;					//距离波门中心
//	int TrackGateRSize;						//距离波门大小
//	int TrackGateVCenter;					//速度波门中心
//	int TrackGateVSize;						//速度波门大小
//	double RangeOnceInfo[536] = { 0 };		//距离一次信息数据
//	double VelocityOnceInfo[64] = { 0 };	//速度一次信息数据
//};


// 结构体包含He、C1、C2、A(多辅路)信号矩阵
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

//// 测量结果
//struct MeasureOutStruct
//{
//	double Range;		// 测量目标距离
//	double Velocity;	// 测量目标速度
//	double AlphaAngle;	// 测量目标方位角虚部
//	double BetaAngle;	// 测量目标俯仰角虚部
//	double Amp;			// 测量目标幅度
//	double SNR;			// 测量目标SNR
//	int Quality;		// 测量目标质量
//};



/*各阶段的处理时间信息*/
struct SGTimeInfo
{
	float	TarTimeInfo;				// 目标生成耗时
	float	ClutterTimeInfo;			// 杂波生成耗时
	float	JamTimeInfo;				// 干扰生成耗时
	float	ReceiverTimeInfo;			// 接收机信号生成耗时
	float	BeforePCTimeInfo;			// 脉压前处理耗时
	float	PCTimeInfo;					// 脉压耗时
	float	TargetDetectTimeInfo;		// 目标检测耗时
	float	ParaMeasureTimeInfo;		// 参数测量耗时
	float	SignalProcessTimeInfo;		// 信号处理耗时
	float	AllProcessTimeInfo;			// 全流程时间

	SGTimeInfo():TarTimeInfo(0.0), ClutterTimeInfo(0.0), JamTimeInfo(0.0), ReceiverTimeInfo(0.0), BeforePCTimeInfo(0.0), PCTimeInfo(0.0), TargetDetectTimeInfo(0.0), ParaMeasureTimeInfo(0.0), SignalProcessTimeInfo(0.0), AllProcessTimeInfo(0.0)
	{}

};

/*内存信息*/
struct SGMemInfo
{
	size_t TarMemInfo;					// 目标生成内存
	size_t ClutterMemInfo;				// 杂波生成内存
	size_t JamMemInfo;					// 干扰生成内存
	size_t ReceiverMemInfo;				// 接收机信号生成内存
	size_t BeforePCMemInfo;				// 脉压前处理内存
	size_t PCMemInfo;					// 脉压内存
	size_t TargetDetectMemInfo;			// 目标检测内存
	size_t ParaMeasureMemInfo;			// 参数测量内存
	size_t SignalProcessMemInfo;		// 信号处理内存
	size_t AllProcessMemInfo;			// 全流程内存

	size_t TarBufInfo;					// 目标生成缓存
	size_t ClutterBufInfo;				// 杂波生成缓存
	size_t JamBufInfo;					// 干扰生成缓存
	size_t ReceiverBufInfo;				// 接收机信号生成缓存
	size_t BeforePCBufInfo;				// 脉压前处理缓存
	size_t PCBufInfo;					// 脉压缓存
	size_t TargetDetectBufInfo;			// 目标检测缓存
	size_t ParaMeasureBufInfo;			// 参数测量缓存
	size_t SignalProcessBufInfo;		// 信号处理缓存
	size_t AllProcessBufInfo;			// 全流程缓存


	SGMemInfo():TarMemInfo(0), ClutterMemInfo(0), JamMemInfo(0), ReceiverMemInfo(0), PCMemInfo(0), TargetDetectMemInfo(0), ParaMeasureMemInfo(0), SignalProcessMemInfo(0), AllProcessMemInfo(0),TarBufInfo(0), ClutterBufInfo(0), JamBufInfo(0), ReceiverBufInfo(0), PCBufInfo(0), TargetDetectBufInfo(0), ParaMeasureBufInfo(0), SignalProcessBufInfo(0), AllProcessBufInfo(0)
	{}
};


#endif