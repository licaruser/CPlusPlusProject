clc;clear all;close all;

Para = "data.bin";
fid = fopen(Para,'r');
[rows,ElementNums] = fread(fid,'int');
data = zeros(128,62500);
row = 1;
col = 1;
for ii = 1:ElementNums
    data(row,col) = rows(ii);
    row = row + 1;
    
    if row > 128
        col = col + 1;
        row = 1;
    end
end

figure
xtable = linspace(1,128,128);
ytable = linspace(1,62500,62500);
mesh(ytable,xtable,data(:,:));

aa = 11; 