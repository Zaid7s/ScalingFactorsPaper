clc
clear all


clc
clear all
clf

%% LUSGS
figure(1)
datfiles = dir('LUSGS*');

for k = 1 :  length(datfiles) - 2
    data = load(datfiles(k).name);
    datfiles(k).name
    Length_k = [1; 5; 11; 21; 31; 41; 51; 55; 61; 71; 81; 91; 101];
    all_marks = {   'o','+','*','.','x','s',...
                    'd','^','v','>','<','p',...
                    'h','o','+','*','o','x',...
                    's','d'};
        plot(data(:, 1), data(:, 4),'LineWidth', 1.0, 'Marker',all_marks{k},...
                            'MarkerIndices',Length_k(k):20:length(data))
        hold on
        ylabel('Mean Pressure', 'Fontsize', 14)
        xlabel('Domain', 'Fontsize', 14)
        xlim([-10 10])
        ylim([0.4 0.8])
    pause(0.001)
end
    plot(data(:,1),data(:,8), 'color', [0.25 0.25 0.25], 'LineWidth',1.5)
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', 'Exact', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    print(['C1P1_SteadyState'], '-depsc', '-r900')

