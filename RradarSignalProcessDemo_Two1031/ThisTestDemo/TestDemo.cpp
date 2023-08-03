#include "TestDemo.h"


//int* remove_element(int* array, int sizeOfArray, int indexToRemove)
//{
//	int* temp = (int*)malloc((sizeOfArray - 1) * sizeof(int)); // allocate an array with a size 1 less than the current one if (indexToRemove != 0) 
//	memcpy(temp, array, indexToRemove * sizeof(int)); // copy everything BEFORE the index if (indexToRemove != (sizeOfArray - 1)) memcpy(temp+indexToRemove, array+indexToRemove+1, 
//	(sizeOfArray - indexToRemove - 1) * sizeof(int)); // copy everything AFTER the index 
//	free(array);
//	return temp;
//}

// ����һ���̣߳�ר��������ź������ͷŵ�ֵ
HANDLE hMutexThread;
// ����һ��������
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
	Title����ѧ���Զ��̴߳���
	Time��2023.08.01
	Author�������
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
//		// �ȴ�������
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
//		// �ͷŻ�����
//		if (!ReleaseMutex(hMutex))
//		{
//			cout << "ReleaseMutex error: " << GetLastError() << endl;
//			return 1;
//		}
//
//		// �رջ������ľ��
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
		// �ȴ�������
		WaitForSingleObject(hMutex, INFINITE);
		//DWORD dwRet = WaitForSingleObject(hMutex, INFINITE);
		//if (dwRet == WAIT_FAILED)
		//{
		//	cout << "WaitForSingleObject error: " << GetLastError() << endl;
		//	return 1;
		//}

		cout << "WaitForSingleObject return value: " << "Successful" << endl;
	}

	// �ͷŻ�����
	if (!ReleaseMutex(hMutex))
	{
		cout << "ReleaseMutex error: " << GetLastError() << endl;
		return 1;
	}

	// �رջ������ľ��
	if (!CloseHandle(hMutex))
	{
		cout << "CloseHandle error: " << GetLastError() << endl;
		return 1;
	}

	return 0;
}

void TestDemo::StepAdvance()
{

	hMutexThread = CreateThread(NULL, NULL, hMutexThreadFun, this, 0, NULL);    //�����߳�
	//hMutex = CreateMutex(NULL, FALSE, NULL);
	hMutex = CreateSemaphore(NULL, 0, 100, NULL);  //�������������100���ź�

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

		//��ǰ��һ���ź���ִ��һ���̣߳��ź�����ʵ��Ҳ�ǿ����������ͻ�����
		ReleaseSemaphore(hMutex, 1, NULL);

		XunHuan++;
	}

}
