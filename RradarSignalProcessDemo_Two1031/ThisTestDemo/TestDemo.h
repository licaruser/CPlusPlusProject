#pragma once


#ifndef TESTDEMO_h__
#define TESTDEMO_h__

#ifdef	TESTDEMO_EXPORTS
#define TESTDEMO_API __declspec(dllexport)
#else		  
#define TESTDEMO_API __declspec(dllimport)
#endif
#include <malloc.h>
#include <stdio.h>
#include <windows.h>
#include <iostream>



using namespace std;

class TESTDEMO_API TestDemo
{
public:
	TestDemo();
	virtual ~TestDemo();

	static DWORD WINAPI hMutexThreadFun(LPVOID lpParam);

	void StepAdvance();

protected:





private:





};




#endif

