#pragma once
#ifndef TESTINPUTPOINT_h__
#define TESTINPUTPOINT_h__

#ifdef	TESTINPUTPOINT_EXPORTS
#define TESTINPUTPOINT_API __declspec(dllexport)
#else		  
#define TESTINPUTPOINT_API __declspec(dllimport)
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
//#include <fftw3.h>
#include <string.h>
#include <vector>
//#include <Dense>
#include <windows.h>//精度时间
#include "libxl_412\include_cpp\libxl.h"


using namespace std;

class TESTINPUTPOINT_API InputTransPoint
{


public:
	//InputTransPoint();
	//~InputTransPoint();


	//void InputTransFunction(vector<double>& InputData,vector<double>& OutputData, double & Span);


protected:



private:




};
#endif