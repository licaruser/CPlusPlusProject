#pragma once
#include "SigGeneratePublicDefinition.h"
#include "CudaArray.cuh"

#ifdef GPUSIGNALGENERATE_EXPORTS
#define GPUSIGNALGENERATE_API __declspec(dllexport)
#else
#define GPUSIGNALGENERATE_API __declspec(dllimport)
#endif

// 回波信号生成
GPUSIGNALGENERATE_API void ReceiverEchoGenerate(CCMat &Data, SGCWStruct &CW, SGInsightStruct &Insight, SGSGAPInitParaStruct &SGAPInitPara, cudaStream_t Stream);

// 回波信号生成，包含数据录取
GPUSIGNALGENERATE_API void ReceiverEchoGenerate(CCMat& Data, SGCWStruct& CW, SGInsightStruct& Insight, SGSGAPInitParaStruct& SGAPInitPara, unsigned char* addrToRecord, cudaStream_t Stream);

// 回波信号生成
GPUSIGNALGENERATE_API void ReceiverEchoGenerateInfo(CCMat &Data, SGCWStruct &CW, SGInsightStruct &Insight, SGSGAPInitParaStruct &SGAPInitPara, SGTimeInfo &TimeInfo, SGMemInfo &MemInfo, cudaStream_t Stream);

// 回波信号生成，包含数据录取
GPUSIGNALGENERATE_API void ReceiverEchoGenerateInfo(CCMat& Data, SGCWStruct& CW, SGInsightStruct& Insight, SGSGAPInitParaStruct& SGAPInitPara, SGTimeInfo& TimeInfo, SGMemInfo& MemInfo, unsigned char* addrToRecord, cudaStream_t Stream);



