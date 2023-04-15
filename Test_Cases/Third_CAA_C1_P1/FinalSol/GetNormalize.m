clc
clear all
clf

figure(1)
datfiles = dir('*.dat');

for k = 1 : 1: length(datfiles)
    data = load(datfiles(k).name); %load just this file
    figure(1)
        plot(data(:,1),data(:,2),'-o','MarkerIndices',1:10:length(data),'LineWidth',2.0)
        hold on
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        pause(0.0001)
end


