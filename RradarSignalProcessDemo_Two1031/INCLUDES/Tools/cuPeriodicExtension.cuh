#ifndef __PERIODICEXTENSION_CUH__
#define __PERIODICEXTENSION_CUH__

#include <cuComplex.h>
//#include "SigGeneratePublicDefinition.h"
#include "tools.cuh"


/*
周期延拓
data: 待处理的一列数据
result:周期延拓后的数据
Period:循环周期数，对应收发周期的点数
Th:延拓有效数，对应侦收时间的点数

example:
Period: 4
Th: 2
data:		1 2 3 4 5 6 1 2 3 4 5 6
周期序列:	1 1 0 0 1 1 0 0 1 1 0 0
result:		0 0 3 4 0 0 1 2 0 0 5 6
*/
void PeriodicExtension(CCMat& Data, CCMat& Result, int Period, int Th, cudaStream_t Stream);




#endif

