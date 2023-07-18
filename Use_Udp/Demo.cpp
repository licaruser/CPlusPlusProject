#include <stdio.h>
#include <WinSock2.h> //windows socket的头文件
#include <ws2tcpip.h> 
#include <Windows.h>
#include <iostream>
#include <thread>
#include <process.h>
#pragma comment(lib, "ws2_32.lib") //连接winsock2.h的静态库文件
using namespace std;

#define COMSOCKET_UDP      //UDP通信_服务器代码控制标志
#define COMSOCKET_TCP1       //TCP通信_服务器代码控制标志

#ifdef COMSOCKET_UDP   //UDP通信测试代码块
int main() {
	SOCKET m_Socket;
	SOCKADDR_IN m_RemoteAddress; //远程地址
	int m_RemoteAddressLen;

	// socket环境
	WSADATA  wsaData;
	if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
		cout << "WSAStartup error:" << GetLastError() << endl;
		return false;
	}

	// socket对象
	m_Socket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (m_Socket == INVALID_SOCKET)
	{
		closesocket(m_Socket);
		m_Socket = INVALID_SOCKET;
		return false;
	}

	// 远端地址
	const char* ip = "192.168.1.103";    //客户端远端IP
	int port = 9700;
	m_RemoteAddress.sin_family = AF_INET;
	m_RemoteAddress.sin_port = htons(port);
	m_RemoteAddressLen = sizeof(m_RemoteAddress);
	inet_pton(AF_INET, ip, &m_RemoteAddress.sin_addr); //将十进制IP转换到sin_addr（二进制）中

	// 接收和发送
	char recvBuf[1024] = { 0 };
	char sendBuf[1024] = "Nice to meet you!我的";

	while (1) {
		int sendLen = sendto(m_Socket, sendBuf, strlen(sendBuf), 0, (sockaddr*)& m_RemoteAddress, m_RemoteAddressLen);
		if (sendLen > 0) {
			std::printf("发送到远程端连接, 其ip: %s, port: %d\n", inet_ntoa(m_RemoteAddress.sin_addr), ntohs(m_RemoteAddress.sin_port));
			cout << "发送到远程端的信息： " << sendBuf << endl;
		}

		int recvLen = recvfrom(m_Socket, recvBuf, 1024, 0, NULL, NULL);  //建立连接之后，无所谓地址，只需要接收就行
		if (recvLen > 0) {
			std::printf("接收到一个连接, 其ip: %s, port: %d\n", inet_ntoa(m_RemoteAddress.sin_addr), ntohs(m_RemoteAddress.sin_port));
			cout << "接收到一个信息： " << recvBuf << endl;
		}
	}

	//没定义客户端发射的地址、端口，只定义了服务端的接收地址、端口
	//服务端定义了，服务端绑定的IP和端口，发送到的目的地无所谓。

	closesocket(m_Socket);
	WSACleanup();
	return true;
}
#endif

#ifdef COMSOCKET_TCP       //TCP通信_服务器代码控制标志
int main()
{
	//加载winsock库
	WSADATA wsadata;
	WSAStartup(MAKEWORD(2, 3), &wsadata);

	//客户端socket
	SOCKET clientSock = socket(PF_INET, SOCK_STREAM, 0);  // clientSock 也是句柄，是linux内核的句柄

	//初始化socket信息
	//memset:作用是在一段内存块中填充某个给定的值，它对较大的结构体或数组进行清零操作的一种最快方法。
	sockaddr_in clientAddr;
	memset(&clientAddr, 0, sizeof(SOCKADDR));

	//设置Socket的连接地址、方式和端口
	clientAddr.sin_addr.s_addr = inet_addr("192.168.1.103");
	clientAddr.sin_family = PF_INET;
	clientAddr.sin_port = htons(9700);

	//建立连接
	connect(clientSock, (SOCKADDR*)& clientAddr, sizeof(SOCKADDR));
	cout << "已建立连接。" << endl;

	//发送消息
	char* s = new char[100];
	cout << "请输入你要发送的文字消息: ";
	cin >> s;
	send(clientSock, s, strlen(s) * sizeof(char) + 1, NULL);
	cout << "已发送:" << s << endl;

	//接收消息
	system("pause");
	char Buffer[MAXBYTE] = { 0 };
	recv(clientSock, Buffer, MAXBYTE, 0);
	cout << "通过端口:" << ntohs(clientAddr.sin_port) << "接收到:" << Buffer << endl;

	//关闭连接
	closesocket(clientSock);
	WSACleanup();
	cout << "客户端连接已关闭。" << endl;

	system("pause");
	return 0;
}
#endif