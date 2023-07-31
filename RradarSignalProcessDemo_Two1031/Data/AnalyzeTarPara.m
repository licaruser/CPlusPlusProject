clc;clear all; close all;

Path1 = "lfm256point.csv";

%SouData_Num1 = dlmread(path1);
%SouData_Num1 = xlsread(Path1);

SouData_Num1 = readlfm256fspoint(Path1,[1,inf]);
SouData_Num1 = SouData_Num1(:,1:256);
SouData_Num1 = table2array(SouData_Num1);
SouData_Num1 = SouData_Num1.';

%SouData_Num1 = hilbert(SouData_Num1);%希尔伯特变换（Hilbert变换）

Fs = 2e9;

xxDouble = linspace(0,100e-6,200000);
M = 256;
f = [0:M/2-1 -M/2:-1]'/M;
%f = linspace(0,-1,M);
figure;
mesh(xxDouble,f(1:M)*Fs,abs(SouData_Num1(1:M,:)))
view(2)
xlabel("时间/s");ylabel("频率/Hz");title("模拟信号时频分析结果");zlabel("幅度")
xlim([0 100e-6]);

% 读取二进制文件 data.bin
fid = fopen('HonstCompleData.bin', 'rb');
nrows = fread(fid, 1, 'size_t');
ncols = fread(fid, 1, 'size_t');

% 将读取的数据转换为二维复数矩阵
data = zeros(nrows, ncols);
for i = 1:nrows
    for j = 1:ncols
        real_part = fread(fid, 1, 'double');
        imag_part = fread(fid, 1, 'double');
        data(i, j) = complex(real_part, imag_part);
    end
end
fclose(fid);


%  AllSignalData = dlmread('HonstCompleData1.txt');
%  AllSignalData = AllSignalData.';
%  
%  SiganlData1 = dlmread('HonstCompleData0.txt');
%  SiganlData1 = SiganlData1.';
 
%  AllSignalData = dlmread('3.txt');
%  AllSignalData = AllSignalData.';
%  cols = length(AllSignalData);
%  figure
%  plot(abs(AllSignalData(1:cols)));
%  
% %  hold on;
% %  plot(real(SiganlData1(1:cols)),'blue');
% 
% %Tar11 = importfile_Tar11('Tar11.txt',[1,inf]);
% Tar11 = dlmread('ModelData.txt');
% Threshold = dlmread('ThresholdValueVectorData.txt');
% 
% %Tar11 = dlmread('expAll.txt');
% %Tar11 = table2array(Tar11);
% Tar11 = Tar11.';
% Threshold = Threshold.';
% 
% figure
% plot(abs(Tar11(1:501)))
% hold on 
% plot(abs(Threshold(1:501)),'red')