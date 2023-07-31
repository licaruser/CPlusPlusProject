% �����˲����
% �Ծ���2020/06/18

clear; close all;

load Motion.txt;

d2r = pi/180.0;

Nsim = size(Motion,1); % �������
t = Motion(:,1); % ʱ���(s)

% �״�վֱ������ϵ�ڵ�TBMĿ��λ�ú��ٶȷ���
RealPosX = Motion(:,2);
RealPosY = Motion(:,3);
RealPosZ = Motion(:,4);
RealVelX = Motion(:,5);
RealVelY = Motion(:,6);
RealVelZ = Motion(:,7);

% �����˲���Ľ��
FilterPosX = Motion(:,8);
FilterPosY = Motion(:,9);
FilterPosZ = Motion(:,10);
FilterVelX = Motion(:,11);
FilterVelY = Motion(:,12);
FilterVelZ = Motion(:,13);
CovPosX = Motion(:,14);CovPosY = Motion(:,15);CovPosZ = Motion(:,16);
CovVelX = Motion(:,17);CovVelY = Motion(:,18);CovVelZ = Motion(:,19);

figure,plot(t,RealPosX,'r',t,FilterPosX,'k');