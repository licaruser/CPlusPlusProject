#ifndef __PERIODICEXTENSION_CUH__
#define __PERIODICEXTENSION_CUH__

#include <cuComplex.h>
//#include "SigGeneratePublicDefinition.h"
#include "tools.cuh"


/*
��������
data: �������һ������
result:�������غ������
Period:ѭ������������Ӧ�շ����ڵĵ���
Th:������Ч������Ӧ����ʱ��ĵ���

example:
Period: 4
Th: 2
data:		1 2 3 4 5 6 1 2 3 4 5 6
��������:	1 1 0 0 1 1 0 0 1 1 0 0
result:		0 0 3 4 0 0 1 2 0 0 5 6
*/
void PeriodicExtension(CCMat& Data, CCMat& Result, int Period, int Th, cudaStream_t Stream);




#endif

