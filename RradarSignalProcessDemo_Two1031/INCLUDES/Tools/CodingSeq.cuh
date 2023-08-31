#ifndef _CODINGSEQ_H_
#define _CODINGSEQ_H_
#include <vector>
#include "CudaArray.cuh"

using namespace std;

void getCodeM3(ICMat &code);
void getCodeM4(ICMat &code);
void getCodeM5(ICMat &code);
void getCodeM6(ICMat &code);
void getCodeM7(ICMat &code);
void getCodeM8(ICMat &code);
void getCodeM9(ICMat &code);
void getCodeM10(ICMat &code);
void getCodeM11(ICMat &code);
void getCodeM12(ICMat &code);
void getCodeM13(ICMat &code);
void getCodeM14(ICMat &code);
void getCodeM15(ICMat &code);

void getComplexCode5(FCMat& PCMCode_5);
void getComplexCode7(FCMat& PCMCode_7);
/* 选择不同的M码
r: 可选 m码: {3,4,5,...,14} 巴克码: {107,,111,113}
*/
void SelectMCode(int r, ICMat &Code);
void SelectComplexCode(int r, FCMat& result);

#endif // _CODINGSEQ_H_
