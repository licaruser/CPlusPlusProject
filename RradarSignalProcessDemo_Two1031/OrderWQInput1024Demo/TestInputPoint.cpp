#include "TestInputPoint.h"
#include <iostream>

using namespace std;

class classA {
public:
	classA() 
	{ 
		cout << " >> classA Ctor"; 
	}
	virtual ~classA() { cout << " >> classA Dtor"; }
};

class interfaceB :public classA 
{
public:
	virtual void methodS() = 0;
};

class classB :public interfaceB 
{
public:
	classB() { cout << " >> classB Ctor"; }
	~classB() { cout << " >> classB Dtor"; }
	virtual void methodS() {};
};


int main()
{
	classA a;
	classB* b = new classB();
	delete b;

	return 0;
	//vector<double> InputData;
	//vector<double> OutpuData;
	//double span;

	//


	//InputTransPoint TmpClassTranPoint;
	//TmpClassTranPoint.InputTransFunction(InputData, OutpuData, span);


	//return 0;
}

//InputTransPoint::InputTransPoint()
//{
//}
//
//InputTransPoint::~InputTransPoint()
//{
//}
//
//void InputTransPoint::InputTransFunction(vector<double>& InputData, vector<double>& OutputData, double& Span)
//{
//	// 已知输入数据，输出的数间隔，输入数据的间隔也知道
//	// 进来的所有的值都是点
//
//}
