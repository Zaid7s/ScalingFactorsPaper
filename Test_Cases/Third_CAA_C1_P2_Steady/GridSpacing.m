clc
clear all
clf

Area    = importdata('Area_Compact4_NPTS00401.dat');

figure(1)
    plot(Area(:, 2), Area(:, 3), 'LineWidth', 2.0)
    hold on
    ylabel('Area')
    xlabel('Domain, x')
    grid on
    ylim([0.2 1.2])
    grid minor
    hold off
    print(['Area'], '-depsc', '-r900')
figure(2)
    plot(Area(:, 2), Area(:, 5), 'LineWidth', 2.0)
    hold on
    ylabel('Grid Distribution, \Delta{x}')
    xlabel('Domain, x')
    ylim([0 0.2])
    grid on
    grid minor
    print(['DX'], '-depsc', '-r900')


