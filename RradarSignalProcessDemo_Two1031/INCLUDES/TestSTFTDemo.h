#pragma once
#ifndef TEST_STFT_DEMO_h__
#define TEST_STFT_DEMO_h__

#ifdef	TEST_STFT_DEMO_EXPORTS
#define TEST_STFT_DEMO_API __declspec(dllexport)
#else		  
#define TEST_STFT_DEMO_API __declspec(dllimport)
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <string.h>
#include <vector>
#include "libxl_412/include_cpp/libxl.h"



#define M_PI 3.1415926535897932384626433


class TEST_STFT_DEMO_API STFT_lgy
{
public:
	STFT_lgy();
	~STFT_lgy();
	void Read_Excel_Data(int num);

	std::vector<complex<double>> Read_Csv_Data(int Init_row,int Init_col);


protected:





private:




};




#endif