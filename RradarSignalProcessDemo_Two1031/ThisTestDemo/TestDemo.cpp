//#include "TestDemo.h"
//#include <vector>
//
////int* remove_element(int* array, int sizeOfArray, int indexToRemove)
////{
////	int* temp = (int*)malloc((sizeOfArray - 1) * sizeof(int)); // allocate an array with a size 1 less than the current one if (indexToRemove != 0) 
////	memcpy(temp, array, indexToRemove * sizeof(int)); // copy everything BEFORE the index if (indexToRemove != (sizeOfArray - 1)) memcpy(temp+indexToRemove, array+indexToRemove+1, 
////	(sizeOfArray - indexToRemove - 1) * sizeof(int)); // copy everything AFTER the index 
////	free(array);
////	return temp;
////}
//
//// 创建一个线程，专用于输出信号量被释放的值
//HANDLE hMutexThread;
//// 创建一个互斥量
//HANDLE hMutex;
//
//class BaseA
//{
//public:
//	BaseA()
//	{
//		cout << "构造A类" << endl;
//	}
//	BaseA(int a)
//	{
//		cout << "构造A类，重载" << endl;
//	}
//	~BaseA()
//	{
//		std::cout << "析构父类A类" << std::endl;
//	}
//};
//
//class SonB :public BaseA
//{
//public:
//	SonB()
//	{
//		std::cout << "构造派生类B类" << std::endl;
//	}
//
//	~SonB()
//	{
//		std::cout << "析构B类" << std::endl;
//	}
//};
//
//void PrinfFunc(int Pri)
//{
//	std::cout << "vaule:" << Pri << std::endl;
//}
//
//void Func(vector<int>& Vaule, void(*func)(int))
//{
//	for (int Pri : Vaule)
//		func(Pri);
//		//std::cout << "vaule:" << Pri << std::endl;
//}
//
//class B 
//{
//	int m_data;
//public:
//	B() { cout << "B constructed\n"; }
//	~B() { cout << "B destructed\n"; }
//	B(int i) : m_data(i) 
//	{ 
//		cout << "B constructed with " << i << "\n"; 
//	}
//};
//
//B PlayB(B bb) 
//{ 
//	return bb;
//}

//#include <iostream>
//#include <exception>
//using namespace std;
//class CBase { virtual void dummy() {} };
//class CDerived : public CBase { int a; };
//int main() {
//	try 
//	{
//		CBase* pba = new CDerived;
//		CBase* pbb = new CBase;
//		CDerived* pd;
//		pd = dynamic_cast<CDerived*>(pba);
//		if (pd == 0) cout << "Null pointer on first type-cast " << endl;
//		pd = dynamic_cast<CDerived*>(pbb);
//		if (pd == 0) cout << "Null pointer on second type-cast " << endl;
//	}
//	catch (exception& e) 
//	{ 
//		cout << "Exception: " << e.what(); 
//	}
//	return 0;
//}

//
//class Base 
//{
//public:
//	Base();
//	virtual ~Base();
//	void func();
//private:
//	int m_DataBase;
//};
//
//class SonAA :public Base
//{
//public:
//	SonAA();
//	virtual ~SonAA();
//	void func();
//private:
//	int m_DataAA;
//};
//
//Base::Base()
//{
//	std::cout << "The function is a Base Class!" << std::endl;
//}
//Base::~Base()
//{
//	std::cout << "The function is a Destructor!" << std::endl;
//}
//void Base::func()
//{
//	std::cout << "func_Base" << std::endl;
//}
//
//SonAA::SonAA()
//{
//	std::cout << "The function is SonAA constructor!" << std::endl;
//}
//SonAA::~SonAA()
//{
//	std::cout << "The function is a SonAA destructor!" << std::endl;
//}
//void SonAA::func()
//{
//	std::cout << "func_SonAA" << std::endl;
//}


//int main() 
//{
//	//std::vector<Animal*> animals;
//	////animals.push_back(new Animal());//含有纯虚函数的抽象类不能实例化
//	//animals.push_back(new Wolf());
//	//animals.push_back(new Fish());
//	//animals.push_back(new GoldFish());
//	//animals.push_back(new OtherAnimal());
//	//for (std::vector<Animal*>::const_iterator it = animals.begin(); it != animals.end(); ++it) {
//	//	(*it)->eat();
//	//	(*it)->runFun();
//	//	delete* it;
//	//}
//
//	Base* aa = new SonAA;
//	aa->func();
//
//	delete(aa);
//
//	return 0;
//}

//#include <iostream>
//#include <vector>
//#include <cmath>
//
////using namespace std;
//# define MY_PI 3.1415926
//// 基类 Shape
//class Shape {
//public:
//	virtual double CalculateArea() = 0;
//};
//
//// 矩形类 Rectangle
//class Rectangle : public Shape 
//{
//public:
//	Rectangle() 
//	{
//
//	}
//	double CalculateArea() override 
//	{
//		return width * height;
//	}
//protected:
//	double width;
//	double height;
//};
//
////// 正方形类 Square
////class Square : public Rectangle {
////public:
////	Square(double side) : Rectangle(side, side) {}
////};
//
//// 圆形类 Circle
//class Circle : public Shape {
//protected:
//	double radius;
//public:
//	Circle(double r) : radius(r) {}
//	double CalculateArea() override 
//	{
//		return MY_PI * radius * radius;
//	}
//};
//
//int main() {
//	vector<Shape*> shapes;
//	Shape* aa = new Rectangle;
//	//shapes.push_back(new Square(5.0));
//	shapes.push_back(new Circle(2.8));
//
//	for (Shape* shape : shapes) {
//		cout << "Area: " << shape->CalculateArea() << endl;
//	}
//
//	// 释放内存
//	for (Shape* shape : shapes) {
//		delete shape;
//	}
//
//	return 0;
//}


//void Animal::runFun() const {}

//#include<stdio.h>          
//int main()
//{
//	char a[] = "hello,world";
//	char* ptr = a;
//	printf("%c\n", *(ptr + 4));
//	printf("%c\n", ptr[4]);
//	printf("%c\n", a[4]);
//	printf("%c\n", *(a + 4));
//	*(ptr + 4) += 1;
//	printf("%s\n", a);
//	return 0;
//}

//#include <iostream>
//#include <list>
//using namespace std;
//int main()
//{
//	list<int> list1;
//	for (int i = 0; i < 8; i++)
//	{
//		list1.push_back(i);
//
//	}
//	for (list<int>::iterator it = list1.begin(); it != list1.end(); ++it)
//	{
//		if (*it % 2 == 0)
//		{
//			list1.erase(it);
//		}
//	}
//	return 0;
//}

//void swap(int *p1, int* p2)
//{
//	int* p = new int;
//	*p = *p1;
//	*p1 = *p2;
//	*p2 = *p;
//}
//
//
//
//
//void main()
//{
//	int* pp = new int;
//	int* cc = new int;
//	*pp = 2;
//	*cc = 3;
//	swap(pp, cc);
//
//	cout << *pp << "," << *cc << endl;
//}



//int main()
//{
//	PlayB(5);
//	{
//		B bb = PlayB(5);
//		int aa;
//	}
//	
//	return 0;
//}
























//
//char* stcpy_my(char* charDest, const char* charSrc)
//{
//	char* charTemp = charDest;
//	while (*charSrc != '\0') 
//	{
//		*charDest = *charSrc;
//		charDest++;
//		charSrc++;
//	}
//	*charDest = '\0';
//	return charTemp;
//}

//class My_String
//{
//public:
//	//My_String();
//	My_String(const char* str = 0 );
//	My_String(const My_String& str);
//	~My_String();
//
//private:
//	char* m_str;
//
//};
//


////My_String::My_String()
////{
////	std::cout << "执行无参构造函数" << std::endl;
////}
//
//My_String::My_String(const char* str)
//{
//	if (str)
//	{
//		m_str = new char[strlen(str) + 1];
//		stcpy_my(m_str, str);
//	}
//	else
//	{
//		std::cout << "执行有参构造函数的else语句" << std::endl;
//		m_str = new char[1];
//		*m_str = '\0';
//	}
//
//	std::cout << m_str << std::endl;
//}
//
//My_String::My_String(const My_String& str)
//{
//
//}
//
//My_String::~My_String()
//{
//}
//
//class A
//{
//public:
//	static int value;
//};
//int value = 5;
//int A::value;

//int main()
//{
//	//int value = 10;
//	//A::value++;
//	//value++;
//	//std::cout << ::value << "" << A::value << "" << value << std::endl;
//	//////int howMany = 20;
//	//////int* test = (int*)malloc(howMany * sizeof(int*));
//	//////for (int i = 0; i < howMany; ++i)
//	//////	(test[i]) = i;
//	//////printf("%dn", test[16]);
//	//////remove_element(test, howMany, 16);
//	//////--howMany;
//	//////printf("%dn", test[16]); return 0;
//	/////************************************
//	////Title：自学测试多线程代码
//	////Time：2023.08.01
//	////Author：李广阳
//	////NUmber:001
//	////************************************/
//	////TestDemo *MyClass = new TestDemo;
//	////MyClass->StepAdvance();
//
//
//	//函数指针
//	//vector<int> Vaules = { 1,3,5,8,9,10 };
//	//Func(Vaules,PrinfFunc);   //传入函数中一个函数，可能也就是将指针传进去了
//
//	//继承以及虚函数使用
//	//BaseA Father = SonB EricB;
//	//My_String s1();
//	//BaseA a[5];
//	//std::vector<int> data = { 1,2,3,4 };
//	//for (std::size_t i = 0; i < data.size(); i++)
//	//{
//	//	std::cout << data[i] << std::endl;
//	//}
//	
//
//	//char src[] = "hello,world!";
//	//char dest[20];
//	//stcpy_my(dest, src);
//	//printf("%s\n", dest);
//
//
//	//std::cout << (7 & 3 + 12) << std::endl;
//
//
//	PlayB();
//	B(4);
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	return 0;
//}
//







































































//
//#include <chrono>
//#include <iostream>
//#include <iostream>
//#include <string>
//#include <chrono>
//
//using namespace std;
//using namespace std::chrono;
//
//
//// 判断是否为闰年
//bool isLeapYear(int year) {
//	if (year % 4 == 0 && year % 100 != 0)
//		return true;
//	if (year % 400 == 0)
//		return true;
//	return false;
//}
//
//// 获取指定年份的天数
//int getDaysInYear(int year) {
//	if (isLeapYear(year))
//		return 366;
//	return 365;
//}
//
//// 获取指定月份的天数
//int getDaysInMonth(int year, int month) {
//	int days[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
//	if (month == 2 && isLeapYear(year))
//		return 29;
//	return days[month - 1];
//}
//
//// 计算两个日期之间相差的天数
//int getDaysBetweenDates(int startYear, int startMonth, int startDay, int endYear, int endMonth, int endDay) {
//	int days = 0;
//
//	// 计算起始日期之后的天数
//	for (int month = startMonth; month <= 12; month++)
//		days += getDaysInMonth(startYear, month);
//	days -= startDay;
//
//	// 计算中间整年的天数
//	for (int year = startYear + 1; year < endYear; year++)
//		days += getDaysInYear(year);
//
//	// 计算结束日期之前的天数
//	for (int month = 1; month < endMonth; month++)
//		days += getDaysInMonth(endYear, month);
//	days += endDay - 1;
//
//	return days;
//}
//
//// 根据身份证号提取出生日期
//string getBirthdayFromId(const string& id) {
//	return id.substr(6, 8);
//}
//
//
//int main() {
//	string id;
//	int year, month, day;
//
//	cin >> id >> year >> month >> day;
//
//	string birthday = getBirthdayFromId(id);
//
//	int currYear = year;
//	int currMonth = month;
//	int currDay = day;
//
//	int birthYear = stoi(birthday.substr(0, 4));
//	int birthMonth = stoi(birthday.substr(4, 2));
//	int birthDay = stoi(birthday.substr(6, 2));
//
//	// 如果当前月日在生日之前,则以去年生日计算
//	if (month > birthMonth || (month == birthMonth && day > birthDay)) 
//	{
//		birthYear = currYear + 1;
//	}
//
//	int days = getDaysBetweenDates(birthYear, birthMonth, birthDay,
//		currYear, currMonth, currDay);
//	
//
//	cout << "Days to next birthday: " << days << endl;
//
//	return 0;
//}
//#include <Stdint.h>
//
//struct A
//{
//	uint8_t  a : 3;
//	uint8_t  b : 5;
//	uint8_t  c : 2;
//	uint8_t  d : 6;
//};
//
//struct A a = { 1,2,3,4 };
//
//int main()
//{
//	uint8_t  b = *((uint8_t*)& a + 1);
//}

//
//#include <stdio.h>
//#include <string>
//
//using namespace std;
//
//void fun(char** m)
//{
//	m++;
//	printf("%s", *m);
//}
//
//void main()
//{
//	static char* a[] = {"MORNING","B","C"};
//
//	char** n;
//	n = a;
//	fun(n);
//
//}



//#include <string> 
//int main()                     //mian函数缺少返回值
//{
//	char* src = "helloworld";
//	char* dest = NULL;
//	int len = strlen(src);
//	dest = (char*)malloc(len);
//	char* d = dest;
//	char* s = &src[len-1];
//	*s--;
//	*s--;
//	printf("%d", *s--);
//	while (len-- != 0)
//	{
//		*d++ = *s--;
//	}
//	//*d = '\0';
//	printf("%s",dest);
//	free(dest);
//	return 0;
//}



//#include <stdio.h>
//void test(void)
//{
//	int zero = 0;
//	int b = 2;
//	bool flag;
//	flag = zero = (b - 2);
//	if (zero = (b - 2))
//	{
//		printf("eq zero");
//	}
//	else
//	{
//		printf("ne zero");
//	}
//}
//
//void main()
//{
//	test();
//}



































//
//
//
//
//
//
//
//
//
//
////DWORD WINAPI hMutexThreadFun(LPVOID lpParam)
////{
////	while (1)
////	{
////		// 等待互斥量
////		WaitForSingleObject(hMutex, INFINITE);
////		//DWORD dwRet = WaitForSingleObject(hMutex, INFINITE);
////		//if (dwRet == WAIT_FAILED)
////		//{
////		//	cout << "WaitForSingleObject error: " << GetLastError() << endl;
////		//	return 1;
////		//}
////
////		cout << "WaitForSingleObject return value: " << "Successful" << endl;
////
////		// 释放互斥量
////		if (!ReleaseMutex(hMutex))
////		{
////			cout << "ReleaseMutex error: " << GetLastError() << endl;
////			return 1;
////		}
////
////		// 关闭互斥量的句柄
////		if (!CloseHandle(hMutex))
////		{
////			cout << "CloseHandle error: " << GetLastError() << endl;
////			return 1;
////		}
////	}
////
////	return 0;
////}
//
//TestDemo::TestDemo()
//{
//}
//
//TestDemo::~TestDemo()
//{
//}
//
//DWORD __stdcall TestDemo::hMutexThreadFun(LPVOID lpParam)
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
//	}
//
//	// 释放互斥量
//	if (!ReleaseMutex(hMutex))
//	{
//		cout << "ReleaseMutex error: " << GetLastError() << endl;
//		return 1;
//	}
//
//	// 关闭互斥量的句柄
//	if (!CloseHandle(hMutex))
//	{
//		cout << "CloseHandle error: " << GetLastError() << endl;
//		return 1;
//	}
//
//	return 0;
//}
//
//void TestDemo::StepAdvance()
//{
//
//	hMutexThread = CreateThread(NULL, NULL, hMutexThreadFun, this, 0, NULL);    //创建线程
//	//hMutex = CreateMutex(NULL, FALSE, NULL);
//	hMutex = CreateSemaphore(NULL, 0, 100, NULL);  //队列中最多阻塞100个信号
//
//	if (hMutex == NULL)
//	{
//		cout << "CreateMutex error: " << GetLastError() << endl;
//		//return 1;
//	}
//
//	int XunHuan = 0;
//
//	while (1)
//	{
//		if (XunHuan == 6)
//		{
//			break;
//		}
//		Sleep(1000);//1000ms
//
//		//当前是一个信号量执行一个线程，信号量的实现也是靠条件变量和互斥锁
//		ReleaseSemaphore(hMutex, 1, NULL);
//
//		XunHuan++;
//	}
//
//}
//
//void Wolf::runFun() const
//{
//	std::cout << "狼儿跑！" << std::endl;
//}
//
//void Fish::runFun() const
//{
//	std::cout << "鱼儿游！" << std::endl;
//}
//
//void GoldFish::runFun() const
//{
//	std::cout << "金鱼儿跳龙门！" << std::endl;
//}
//
//void OtherAnimal::runFun() const
//{
//	std::cout << "其他动物悠哉游哉走！" << std::endl;
//}
//
//
