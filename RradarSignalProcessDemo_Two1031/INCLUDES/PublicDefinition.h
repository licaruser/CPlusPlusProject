#ifndef PUBLICDEFINITION_H
#define PUBLICDEFINITION_H

//#include "Exports.h"
#include <vector>
#include <map>
#include <stack>
#include <algorithm>
#include <complex>



using namespace std;

struct PulseBand
{
	int minGle;
	int maxGle;
	int BandWidth;

	PulseBand()
	{
		memset(this, 0, sizeof(PulseBand));
	};

	//PulseBand operator = (const PulseBand& it)
	//{

	//	minGle = it.minGle;
	//	maxGle = it.maxGle;
	//	BandWidth = it.BandWidth;
	//};

};


#endif // !PUBLICDEFINITION_H



#pragma once
