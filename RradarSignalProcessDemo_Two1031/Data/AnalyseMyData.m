clc;clear all; close all;

xxData = [16 32 64 128 256 512 1024];
yy1Data = [0.326831 0.60012 1.193 2.386 4.92375 10.3315 22.8668];
yy2Data = [1.00272 1.80688 3.57798 7.44231 15.7472 33.0173 71.2359];

figure
plot(xxData,yy1Data,'-*','color','#0072BD');
hold on;
plot(xxData,yy2Data,'>-','color','#7E2F8E');
grid on;
xlabel('频率分辨点数/个');ylabel('耗时/sec');legend('采样率：625MHz','采样率：2GHz')
title('基于CPU的时频耗时统计图');


zp = BaseZoom();
zp.plot;
% x = 50;
% y = 300;
% size = 100;
% 
% axes('Position',[0.2 0.55 0.3 0.3]);
% plot((y:y+size,x:x+size),:)

