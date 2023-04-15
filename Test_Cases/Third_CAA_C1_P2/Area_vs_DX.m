clc
clear all
clf

Area    = importdata('Area_SecOrder_NPTS00401.dat');

figure(1)
    yyaxis right
        plot(Area(:, 2), Area(:, 3), 'LineWidth', 1.0)
        hold on
        ylabel('Area')
        xlabel('Domain, x')
        grid on
        ylim([0.2 1.2])
        xlim([-10 10])
        grid minor
        hold off
    yyaxis left
        plot(Area(:, 2), Area(:, 5), '-o','MarkerIndices',1:20:length(Area(:, 5)),'LineWidth', 1.0)
        hold on
        ylabel('Grid Distribution, \Delta{x}')
        xlabel('Domain, x')
        ylim([0 0.2])
        grid on
        grid minor
        xlim([-10 10])
%         print(['Area_dX'], '-depsc', '-r900')

