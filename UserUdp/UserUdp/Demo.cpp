#include <stdio.h>
#include <WinSock2.h> //windows socket��ͷ�ļ�
#include <Windows.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <process.h>
#pragma comment(lib, "ws2_32.lib") //����winsock2.h�ľ�̬���ļ�
using namespace std;


#define COMSOCKET_UDP      //UDPͨ��_������������Ʊ�־
#define COMSOCKET_TCP1       //TCPͨ��_������������Ʊ�־

#ifdef COMSOCKET_UDP   //UDPͨ�Ų��Դ����
int main() {
	SOCKET m_Socket;
	SOCKADDR_IN m_BindAddress;   //�󶨵�ַ
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

	// ��ռ��<ip, port>
	const char* ip = "192.168.1.103";
	int port = 9700;
	m_BindAddress.sin_family = AF_INET;
	m_BindAddress.sin_addr.S_un.S_addr = inet_addr(ip);
	m_BindAddress.sin_port = htons(port);
	auto ret = bind(m_Socket, (sockaddr*)& m_BindAddress, sizeof(SOCKADDR));  //��ô�о�socketֻ��Ҫ���ӵ��Լ��ı���IP�Ϳ����أ�
	if (ret == SOCKET_ERROR)
	{
		closesocket(m_Socket);
		m_Socket = INVALID_SOCKET;
		return false;
	}

	//����ip��port����
	const char* sendip = "192.168.1.200";
	int sendport = 9700;
	struct sockaddr_in form;
	form.sin_family = AF_INET;
	form.sin_addr.S_un.S_addr = inet_addr(sendip);
	form.sin_port = htons(sendport);
	int form_len = sizeof(form);


	// ���պͷ���
	char recvBuf[1024] = { 0 };
	char sendBuf[1024] = "Nice to meet you, too!";
	//m_RemoteAddressLen = sizeof(m_RemoteAddress);

	std::printf("�����ð�ռ�õ�����, ��ip: %s, port: %d\n", inet_ntoa(m_BindAddress.sin_addr), ntohs(m_BindAddress.sin_port));

	while (1) {
		int recvLen = recvfrom(m_Socket, recvBuf, 1024, 0, (sockaddr*)& form, &form_len); //���յĵ�ַ�о�����ν������ĺ���Ҳ���ܣ�������Լ�����
		if (recvLen > 0) {
			std::printf("���յ�һ������, ��ip: %s, port: %d\n", inet_ntoa(form.sin_addr), ntohs(form.sin_port));
			cout << "���յ�һ����Ϣ�� " << recvBuf << endl;


		}
		int sendLen = sendto(m_Socket, sendBuf, strlen(sendBuf), 0, (sockaddr*)& form, form_len);
		if (sendLen > 0) {
			cout << "���͵�Զ�̶˵���Ϣ�� " << sendBuf << endl;
		}
	}

	closesocket(m_Socket);
	WSACleanup();
	return true;
}

#endif // COMSOCKET_UDP

#ifdef COMSOCKET_TCP
	//mutex ÿ���߳��ڶ���Դ����ǰ�������ȼ������ɹ��������ܲ�������������������
	//ͬһʱ�̣�ֻ����һ���̳߳��и�����
	mutex m;
	//����ṹ����������
	typedef struct my_file
	{
		SOCKET clientSocket; //�ļ��ڲ�������һ��SOCKET ���ںͿͻ��˽���ͨ��
		sockaddr_in clientAddr; //���ڱ���ͻ��˵�socket��ַ
		int id; //�ļ�������
	}F;

	DWORD WINAPI transmmit(const LPVOID arg)
	{
		//ʵ��������Ϊ��׷�󲢷��Բ�Ӧ�ü�����������Ϊ�˷��㿴���
		m.lock();

		F* temp = (F*)arg;
		//��ȡ�ļ������
		//int file_id = temp->id;
		//��ȡ�ͻ����Ķ˿ں�
		//ntohs(temp -> clientAddr.sin_port); 
		cout << "���Կ�ʼ,�ȴ��ͻ��˷�����Ϣ..." << endl;
		//�ӿͻ��˴���������
		char Buffer[MAXBYTE] = { 0 }; //������
		recv(temp->clientSocket, Buffer, MAXBYTE, 0); //recv���� �ӿͻ���ͨ��clientScocket����
		cout << "�߳�" << temp->id << "�ӿͻ��˵�" << ntohs(temp->clientAddr.sin_port) << "�Ŷ˿��յ�:" << Buffer << endl;

		//���ͼ򵥵��ַ������ͻ���
		const char* s = "Server file";
		send(temp->clientSocket, s, strlen(s) * sizeof(char) + 1, NULL);
		cout << "�߳�" << temp->id << "ͨ���ͻ��˵�" << ntohs(temp->clientAddr.sin_port) << "�Ŷ˿ڷ���:" << s << endl;

		m.unlock();

		return 0;
	}

	int main()
	{
		//����winsock��,��һ��������winsocket load�İ汾�ţ�2.3��
		WSADATA wsaData;
		WSAStartup(MAKEWORD(2, 3), &wsaData);

		//�����������˵�socket��Э���壬 sokcet���ͣ�
		SOCKET servSocket = socket(AF_INET, SOCK_STREAM, 0);//����ĳ�SOCK_DGRAM��ʹ��UDP

		// ��ʼ��socket��Ϣ
		sockaddr_in servAddr; //��������socket��ַ������sin_addr��ʾIP��ַ��sin_port���ֶ˿ںź�sin_zero����ֽ�
		memset(&servAddr, 0, sizeof(SOCKADDR)); //��ʼ��socket��ַ

		//����Socket�����ӵ�ַ����ʽ�Ͷ˿ڣ�����
		servAddr.sin_family = PF_INET; //����ʹ�õ�Э����
		servAddr.sin_port = htons(9700); //����ʹ�õĶ˿�
		servAddr.sin_addr.s_addr = inet_addr("192.168.1.103"); //���ð󶨵�IP��ַ
		::bind(servSocket, (SOCKADDR*)& servAddr, sizeof(SOCKADDR)); //��֮ǰ������servSocket�Ͷ˿ڣ�IP��ַ��

		HANDLE hThread[20]; //��ȡ���
		listen(servSocket, 20); //�����������˿ڣ����20������
		for (int i = 0; i < 20; i++)
		{
			F* temp = new F; //�����µĴ���ṹ��
			sockaddr_in clntAddr;
			int nSize = sizeof(SOCKADDR);
			SOCKET clientSock = accept(servSocket, (SOCKADDR*)& clntAddr, &nSize);
			//temp���ݳ�Ա��ֵ
			temp->clientSocket = clientSock;
			temp->id = i + 1;
			temp->clientAddr = clntAddr;
			//ͨ������������߳�
			hThread[i] = CreateThread(NULL, 0, &transmmit, temp, 0, NULL);
		}

		//�ȴ����߳����
		WaitForMultipleObjects(20, hThread, TRUE, INFINITE);
		cout << WSAGetLastError() << endl; //�鿴������Ϣ

		//�ر�socket���ͷ�winsock
		closesocket(servSocket);
		WSACleanup();

		cout << "�����������ѹرա�" << endl;
		system("pause");

		return 0;
}
#endif // COMSOCKET_TCP

