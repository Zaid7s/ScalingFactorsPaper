clc
clear all
clf

figure(1)
datfiles = dir('*.dat');
% datfiles = dir('*RDRP*');
% datfiles = dir('*Fou*');
% datfiles = dir('Unsteady_TriDi*');
% datfiles = dir('Unsteady_LUSGS*');

for k = 1 : length(datfiles)
    data = load(datfiles(k).name); %load just this file
    (datfiles(k).name)
    figure(1)
    subplot(3, 1, 1)
        plot(data(:,1),data(:,2),'-o','MarkerIndices',1:10:length(data),'LineWidth',2.0)
        hold on
        plot(data(:,1),data(:,6),'LineWidth',2.0)
        hold off
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        ylim([0.7 1.1])
    subplot(3, 1, 2)
        plot(data(:,1),data(:,3),'-x','MarkerIndices',1:10:length(data),'LineWidth',2.0)
        hold on
        plot(data(:,1),data(:,7),'LineWidth',2.0)
        hold off
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        ylim([0.1 1.0])
    subplot(3, 1, 3)
        plot(data(:,1),data(:,4),'-*','MarkerIndices',1:10:length(data),'LineWidth',2.0)
        hold on
        plot(data(:,1),data(:,8),'LineWidth',2.0)
        hold off
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        ylim([0.4 0.8])
        ylabel('\rho')
        pause(0.0001)
end
