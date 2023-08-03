#include "TestDemo.h"


//int* remove_element(int* array, int sizeOfArray, int indexToRemove)
//{
//	int* temp = (int*)malloc((sizeOfArray - 1) * sizeof(int)); // allocate an array with a size 1 less than the current one if (indexToRemove != 0) 
//	memcpy(temp, array, indexToRemove * sizeof(int)); // copy everything BEFORE the index if (indexToRemove != (sizeOfArray - 1)) memcpy(temp+indexToRemove, array+indexToRemove+1, 
//	(sizeOfArray - indexToRemove - 1) * sizeof(int)); // copy everything AFTER the index 
//	free(array);
//	return temp;
//}

// 创建一个线程，专用于输出信号量被释放的值
HANDLE hMutexThread;
// 创建一个互斥量
HANDLE hMutex;

int main()
{
	//int howMany = 20;
	//int* test = (int*)malloc(howMany * sizeof(int*));
	//for (int i = 0; i < howMany; ++i)
	//	(test[i]) = i;
	//printf("%dn", test[16]);
	//remove_element(test, howMany, 16);
	//--howMany;
	//printf("%dn", test[16]); return 0;

	/************************************
	Title：自学测试多线程代码
	Time：2023.08.01
	Author：李广阳
	NUmber:001
	************************************/
	TestDemo *MyClass = new TestDemo;
	MyClass->StepAdvance();





	return 0;
}

//DWORD WINAPI hMutexThreadFun(LPVOID lpParam)
//{
//	while (1)
//	{
//		// 等待互斥量
//		WaitForSingleObject(hMutex, INFINITE);
//		//DWORD dwRet = WaitForSingleObject(hMutex, INFINITE);
//		//if (dwRet == WAIT_FAILED)
//		//{
//		//	cout << "WaitForSingleObject error: " << GetLastError() << endl;
//		//	return 1;
//		//}
//
//		cout << "WaitForSingleObject return value: " << "Successful" << endl;
//
//		// 释放互斥量
//		if (!ReleaseMutex(hMutex))
//		{
//			cout << "ReleaseMutex error: " << GetLastError() << endl;
//			return 1;
//		}
//
//		// 关闭互斥量的句柄
//		if (!CloseHandle(hMutex))
//		{
//			cout << "CloseHandle error: " << GetLastError() << endl;
//			return 1;
//		}
//	}
//
//	return 0;
//}

TestDemo::TestDemo()
{
}

TestDemo::~TestDemo()
{
}

DWORD __stdcall TestDemo::hMutexThreadFun(LPVOID lpParam)
{
	while (1)
	{
		// 等待互斥量
		WaitForSingleObject(hMutex, INFINITE);
		//DWORD dwRet = WaitForSingleObject(hMutex, INFINITE);
		//if (dwRet == WAIT_FAILED)
		//{
		//	cout << "WaitForSingleObject error: " << GetLastError() << endl;
		//	return 1;
		//}

		cout << "WaitForSingleObject return value: " << "Successful" << endl;
	}

	// 释放互斥量
	if (!ReleaseMutex(hMutex))
	{
		cout << "ReleaseMutex error: " << GetLastError() << endl;
		return 1;
	}

	// 关闭互斥量的句柄
	if (!CloseHandle(hMutex))
	{
		cout << "CloseHandle error: " << GetLastError() << endl;
		return 1;
	}

	return 0;
}

void TestDemo::StepAdvance()
{

	hMutexThread = CreateThread(NULL, NULL, hMutexThreadFun, this, 0, NULL);    //创建线程
	//hMutex = CreateMutex(NULL, FALSE, NULL);
	hMutex = CreateSemaphore(NULL, 0, 100, NULL);  //队列中最多阻塞100个信号

	if (hMutex == NULL)
	{
		cout << "CreateMutex error: " << GetLastError() << endl;
		//return 1;
	}

	int XunHuan = 0;

	while (1)
	{
		if (XunHuan == 6)
		{
			break;
		}
		Sleep(1000);//1000ms

		//当前是一个信号量执行一个线程，信号量的实现也是靠条件变量和互斥锁
		ReleaseSemaphore(hMutex, 1, NULL);

		XunHuan++;
	}

}
