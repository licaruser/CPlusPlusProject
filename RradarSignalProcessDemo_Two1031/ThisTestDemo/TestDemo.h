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
#include <string.h>
#include <vector>



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




class Animal {
public:
	virtual void eat() const
	{
		std::cout << "I eat like a generic Animal. " << std::endl;
	}
	virtual ~Animal()
	{
	}

	virtual void runFun() const = 0;//纯虚函数，必须重新定义；
};
//狼
class Wolf : public Animal
{
public:
	void eat() const
	{
		std::cout << "I eat like a wolf!" << std::endl;
	}

	void runFun() const;
};


class Fish : public Animal {
public:
	void eat() const
	{
		std::cout << "I eat like a fish!" << std::endl;
	}
	void runFun() const;

};


class GoldFish : public Fish
{
public:
	void eat() const
	{
		std::cout << "I eat like a goldfish!" << std::endl;
	}

	void runFun() const;

};


class OtherAnimal : public Animal
{
	void runFun() const;
};





#endif

