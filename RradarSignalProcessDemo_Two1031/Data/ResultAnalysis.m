% 分析滤波结果
% 赵晶，2020/06/18

clear; close all;

load Motion.txt;

d2r = pi/180.0;

Nsim = size(Motion,1); % 仿真点数
t = Motion(:,1); % 时间戳(s)

% 雷达站直角坐标系内的TBM目标位置和速度分量
RealPosX = Motion(:,2);
RealPosY = Motion(:,3);
RealPosZ = Motion(:,4);
RealVelX = Motion(:,5);
RealVelY = Motion(:,6);
RealVelZ = Motion(:,7);

% 经过滤波后的结果
FilterPosX = Motion(:,8);
FilterPosY = Motion(:,9);
FilterPosZ = Motion(:,10);
FilterVelX = Motion(:,11);
FilterVelY = Motion(:,12);
FilterVelZ = Motion(:,13);
CovPosX = Motion(:,14);CovPosY = Motion(:,15);CovPosZ = Motion(:,16);
CovVelX = Motion(:,17);CovVelY = Motion(:,18);CovVelZ = Motion(:,19);

figure,plot(t,RealPosX,'r',t,FilterPosX,'k');