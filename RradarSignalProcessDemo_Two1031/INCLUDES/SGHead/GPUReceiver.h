#pragma once
#include "SigGeneratePublicDefinition.h"
#include "CudaArray.cuh"

#ifdef GPUSIGNALGENERATE_EXPORTS
#define GPUSIGNALGENERATE_API __declspec(dllexport)
#else
#define GPUSIGNALGENERATE_API __declspec(dllimport)
#endif

// �ز��ź�����
GPUSIGNALGENERATE_API void ReceiverEchoGenerate(CCMat &Data, SGCWStruct &CW, SGInsightStruct &Insight, SGSGAPInitParaStruct &SGAPInitPara, cudaStream_t Stream);

// �ز��ź����ɣ���������¼ȡ
GPUSIGNALGENERATE_API void ReceiverEchoGenerate(CCMat& Data, SGCWStruct& CW, SGInsightStruct& Insight, SGSGAPInitParaStruct& SGAPInitPara, unsigned char* addrToRecord, cudaStream_t Stream);

// �ز��ź�����
GPUSIGNALGENERATE_API void ReceiverEchoGenerateInfo(CCMat &Data, SGCWStruct &CW, SGInsightStruct &Insight, SGSGAPInitParaStruct &SGAPInitPara, SGTimeInfo &TimeInfo, SGMemInfo &MemInfo, cudaStream_t Stream);

// �ز��ź����ɣ���������¼ȡ
GPUSIGNALGENERATE_API void ReceiverEchoGenerateInfo(CCMat& Data, SGCWStruct& CW, SGInsightStruct& Insight, SGSGAPInitParaStruct& SGAPInitPara, SGTimeInfo& TimeInfo, SGMemInfo& MemInfo, unsigned char* addrToRecord, cudaStream_t Stream);



