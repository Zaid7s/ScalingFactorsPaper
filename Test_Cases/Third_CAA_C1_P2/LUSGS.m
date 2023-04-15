clc
clear all
clf

Second_TriDi    = importdata('LUSGS_SecOrder_D_10_LHS_Diss_OffMean.dat');
Fourth_TriDi    = importdata('LUSGS_FouOrder_D_10_LHS_Diss_OffMean.dat');
Sixth_TriDi     = importdata('LUSGS_SixOrder_D_10_LHS_Diss_OffMean.dat');
RDRP_TriDi      = importdata('LUSGS_RDRPSten_D_10_LHS_Diss_OffMean.dat');

Exact           = importdata('Exact_SixOrder.dat');

%% LUSGS
figure(1)
    plot(Second_TriDi(:, 1), Second_TriDi(:, 4), 'o', 'MarkerIndices', 1:20:length(Second_TriDi), 'LineWidth', 1.0)
    hold on
    plot(Fourth_TriDi(:, 1), Fourth_TriDi(:, 4), 'x', 'MarkerIndices', 5:20:length(Second_TriDi), 'LineWidth', 1.0)
    hold on
    plot(Sixth_TriDi(:, 1), Sixth_TriDi(:, 4), 's', 'MarkerIndices', 10:20:length(Second_TriDi), 'LineWidth', 1.0)
    hold on
    plot(RDRP_TriDi(:, 1), RDRP_TriDi(:, 4), '^', 'MarkerIndices', 15:20:length(Second_TriDi), 'LineWidth', 1.0)
    hold on
    plot(Exact(:, 1), Exact(:, 4), 'color', [0.25 0.25 0.25], 'LineWidth', 1.5)
    hold on
    legend('E2','E4','E6','RDRP','Exact','Location','Best')
    grid on
    grid minor
    xlabel('Domain, x')
    ylabel('Mean Pressure')
    xlim([-10 10])
%     print(['LUSGS'], '-depsc', '-r900')
    
