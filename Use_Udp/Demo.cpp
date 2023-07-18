#include <stdio.h>
#include <WinSock2.h> //windows socket��ͷ�ļ�
#include <ws2tcpip.h> 
#include <Windows.h>
#include <iostream>
#include <thread>
#include <process.h>
#pragma comment(lib, "ws2_32.lib") //����winsock2.h�ľ�̬���ļ�
using namespace std;

#define COMSOCKET_UDP      //UDPͨ��_������������Ʊ�־
#define COMSOCKET_TCP1       //TCPͨ��_������������Ʊ�־

#ifdef COMSOCKET_UDP   //UDPͨ�Ų��Դ����
int main() {
	SOCKET m_Socket;
	SOCKADDR_IN m_RemoteAddress; //Զ�̵�ַ
	int m_RemoteAddressLen;

	// socket����
	WSADATA  wsaData;
	if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
		cout << "WSAStartup error:" << GetLastError() << endl;
		return false;
	}

	// socket����
	m_Socket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (m_Socket == INVALID_SOCKET)
	{
		closesocket(m_Socket);
		m_Socket = INVALID_SOCKET;
		return false;
	}

	// Զ�˵�ַ
	const char* ip = "192.168.1.103";    //�ͻ���Զ��IP
	int port = 9700;
	m_RemoteAddress.sin_family = AF_INET;
	m_RemoteAddress.sin_port = htons(port);
	m_RemoteAddressLen = sizeof(m_RemoteAddress);
	inet_pton(AF_INET, ip, &m_RemoteAddress.sin_addr); //��ʮ����IPת����sin_addr�������ƣ���

	// ���պͷ���
	char recvBuf[1024] = { 0 };
	char sendBuf[1024] = "Nice to meet you!�ҵ�";

	while (1) {
		int sendLen = sendto(m_Socket, sendBuf, strlen(sendBuf), 0, (sockaddr*)& m_RemoteAddress, m_RemoteAddressLen);
		if (sendLen > 0) {
			std::printf("���͵�Զ�̶�����, ��ip: %s, port: %d\n", inet_ntoa(m_RemoteAddress.sin_addr), ntohs(m_RemoteAddress.sin_port));
			cout << "���͵�Զ�̶˵���Ϣ�� " << sendBuf << endl;
		}

		int recvLen = recvfrom(m_Socket, recvBuf, 1024, 0, NULL, NULL);  //��������֮������ν��ַ��ֻ��Ҫ���վ���
		if (recvLen > 0) {
			std::printf("���յ�һ������, ��ip: %s, port: %d\n", inet_ntoa(m_RemoteAddress.sin_addr), ntohs(m_RemoteAddress.sin_port));
			cout << "���յ�һ����Ϣ�� " << recvBuf << endl;
		}
	}

	//û����ͻ��˷���ĵ�ַ���˿ڣ�ֻ�����˷���˵Ľ��յ�ַ���˿�
	//����˶����ˣ�����˰󶨵�IP�Ͷ˿ڣ����͵���Ŀ�ĵ�����ν��

	closesocket(m_Socket);
	WSACleanup();
	return true;
}
#endif

#ifdef COMSOCKET_TCP       //TCPͨ��_������������Ʊ�־
int main()
{
	//����winsock��
	WSADATA wsadata;
	WSAStartup(MAKEWORD(2, 3), &wsadata);

	//�ͻ���socket
	SOCKET clientSock = socket(PF_INET, SOCK_STREAM, 0);  // clientSock Ҳ�Ǿ������linux�ں˵ľ��

	//��ʼ��socket��Ϣ
	//memset:��������һ���ڴ�������ĳ��������ֵ�����Խϴ�Ľṹ�������������������һ����췽����
	sockaddr_in clientAddr;
	memset(&clientAddr, 0, sizeof(SOCKADDR));

	//����Socket�����ӵ�ַ����ʽ�Ͷ˿�
	clientAddr.sin_addr.s_addr = inet_addr("192.168.1.103");
	clientAddr.sin_family = PF_INET;
	clientAddr.sin_port = htons(9700);

	//��������
	connect(clientSock, (SOCKADDR*)& clientAddr, sizeof(SOCKADDR));
	cout << "�ѽ������ӡ�" << endl;

	//������Ϣ
	char* s = new char[100];
	cout << "��������Ҫ���͵�������Ϣ: ";
	cin >> s;
	send(clientSock, s, strlen(s) * sizeof(char) + 1, NULL);
	cout << "�ѷ���:" << s << endl;

	//������Ϣ
	system("pause");
	char Buffer[MAXBYTE] = { 0 };
	recv(clientSock, Buffer, MAXBYTE, 0);
	cout << "ͨ���˿�:" << ntohs(clientAddr.sin_port) << "���յ�:" << Buffer << endl;

	//�ر�����
	closesocket(clientSock);
	WSACleanup();
	cout << "�ͻ��������ѹرա�" << endl;

	system("pause");
	return 0;
}
#endif