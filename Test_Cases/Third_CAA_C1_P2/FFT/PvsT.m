clc
clear all
datfiles  = dir('*.dat');
for k = 1: length(datfiles)
    data  = load(datfiles(k).name);
    figure(1)
        plot(data(:, 1), data(:, 2), 'LineWidth', 1.0)
        hold on
        grid on
        grid minor
        xlabel('Time')
        ylabel('Pressure Perturbation')
        title('Exit Pressure')
        hold off
        ylim([-0.000015 0.000015])
        pause(0.001)
end
