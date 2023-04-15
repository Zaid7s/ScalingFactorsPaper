clc
clear all
clf

figure(1)
datfiles = dir('*_Sec*');
% datfiles = dir('*_Fou*');
% datfiles = dir('*_Six*');
% datfiles = dir('*_RDRP*');
% datfiles = dir('*_DRP*');
% datfiles = dir('*_Compact*');

for k = 1 : 4: length(datfiles)
    data = load(datfiles(k).name); %load just this file
    figure(1)
        plot(data(:,1),data(:,4), 'color', [0.25 0.25 0.25],'LineWidth',2.0)
%         hold on
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        pause(0.0001)
end

