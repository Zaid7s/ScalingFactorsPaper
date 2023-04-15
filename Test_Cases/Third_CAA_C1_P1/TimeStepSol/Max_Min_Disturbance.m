clc
clear all

Second_LUSGS_Off01      = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_036_Max_Min_Pres.dat');
Fourth_LUSGS_Off01      = importdata('LUSGS_LHS_Diss_Off_FouOrder_D_10_036_Max_Min_Pres.dat');
Sixth_LUSGS_Off01       = importdata('LUSGS_LHS_Diss_Off_SixOrder_D_10_036_Max_Min_Pres.dat');
RDRP_LUSGS_Off01        = importdata('LUSGS_LHS_Diss_Off_RDRPSten_D_10_036_Max_Min_Pres.dat');
DRP_LUSGS_Off01         = importdata('LUSGS_LHS_Diss_Off__DRPSten_D_10_036_Max_Min_Pres.dat');
Compact_LUSGS_Off01     = importdata('LUSGS_LHS_Diss_Off_Compact4_D_10_036_Max_Min_Pres.dat');
% PreFac4_LUSGS_Off01     = importdata('LUSGS_LHS_Diss_Off_PreFac_4_D_10_036_Max_Min_Pres.dat');
% PreFac6_LUSGS_Off01     = importdata('LUSGS_LHS_Diss_Off_PreFac_6_D_10_036_Max_Min_Pres.dat');

Exact       = importdata('CAT1P1.FIG1.DATA'); 

%% LUSGS Scaling Factors
figure(1)
    plot(Second_LUSGS_Off01(:, 1), Second_LUSGS_Off01(:, 5),'-o','MarkerIndices',1:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Fourth_LUSGS_Off01(:, 1), Fourth_LUSGS_Off01(:, 5),'-x','MarkerIndices',5:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Sixth_LUSGS_Off01(:, 1), Sixth_LUSGS_Off01(:, 5),'-*','MarkerIndices',10:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(RDRP_LUSGS_Off01(:, 1), RDRP_LUSGS_Off01(:, 5),'-s','MarkerIndices',15:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(DRP_LUSGS_Off01(:, 1), DRP_LUSGS_Off01(:, 5),'-h','MarkerIndices',20:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Compact_LUSGS_Off01(:, 1), Compact_LUSGS_Off01(:, 5),'-d','MarkerIndices',25:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
%     plot(PreFac4_LUSGS_Off01(:, 1), PreFac4_LUSGS_Off01(:, 5),'-.','MarkerIndices',30:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
%     hold on
%     plot(PreFac6_LUSGS_Off01(:, 1), PreFac6_LUSGS_Off01(:, 5),'-p','MarkerIndices',35:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
%     hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:40:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
    xlim([-10 10])
    xlabel('Domain', 'Fontsize', 14)
    ylabel('Maximum Pressure Perturbation', 'Fontsize', 14)
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', 'Exact', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    print(['C1P1_MaxDisturbance'], '-depsc', '-r900')
%% Zoomed In
figure(2)
subplot(1, 2, 1)
    plot(Second_LUSGS_Off01(:, 1), Second_LUSGS_Off01(:, 5),'-o','MarkerIndices',1:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Fourth_LUSGS_Off01(:, 1), Fourth_LUSGS_Off01(:, 5),'-x','MarkerIndices',5:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Sixth_LUSGS_Off01(:, 1), Sixth_LUSGS_Off01(:, 5),'-*','MarkerIndices',10:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(RDRP_LUSGS_Off01(:, 1), RDRP_LUSGS_Off01(:, 5),'-s','MarkerIndices',15:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(DRP_LUSGS_Off01(:, 1), DRP_LUSGS_Off01(:, 5),'-h','MarkerIndices',20:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Compact_LUSGS_Off01(:, 1), Compact_LUSGS_Off01(:, 5),'-d','MarkerIndices',25:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
%     plot(PreFac4_LUSGS_Off01(:, 1), PreFac4_LUSGS_Off01(:, 5),'-|','MarkerIndices',30:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
%     hold on
%     plot(PreFac6_LUSGS_Off01(:, 1), PreFac6_LUSGS_Off01(:, 5),'-p','MarkerIndices',35:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
%     hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:40:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
    xlim([-10 10])
    ylim([5*10^-6 12*10^-6])
    xlabel('Domain', 'Fontsize', 14)
    ylabel('Maximum Pressure Perturbation', 'Fontsize', 14)
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', 'Exact', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
subplot(1, 2, [2])
    plot(Second_LUSGS_Off01(:, 1), Second_LUSGS_Off01(:, 5),'-o','MarkerIndices',1:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Fourth_LUSGS_Off01(:, 1), Fourth_LUSGS_Off01(:, 5),'-x','MarkerIndices',5:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Sixth_LUSGS_Off01(:, 1), Sixth_LUSGS_Off01(:, 5),'-*','MarkerIndices',10:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(RDRP_LUSGS_Off01(:, 1), RDRP_LUSGS_Off01(:, 5),'-s','MarkerIndices',15:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(DRP_LUSGS_Off01(:, 1), DRP_LUSGS_Off01(:, 5),'-h','MarkerIndices',20:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
    plot(Compact_LUSGS_Off01(:, 1), Compact_LUSGS_Off01(:, 5),'-d','MarkerIndices',25:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
    hold on
%     plot(PreFac4_LUSGS_Off01(:, 1), PreFac4_LUSGS_Off01(:, 5),'-|','MarkerIndices',30:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
%     hold on
%     plot(PreFac6_LUSGS_Off01(:, 1), PreFac6_LUSGS_Off01(:, 5),'-p','MarkerIndices',35:40:length(Second_LUSGS_Off01),'LineWidth', 1.0)
%     hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:40:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
    xlim([-1 1])
    ylim([0 1.2*10^-4])
    xlabel('Domain', 'Fontsize', 14)
    set(gcf, 'Position', [372 527 749 420])
    ylabel('Maximum Pressure Perturbation', 'Fontsize', 14)
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', 'Exact', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    print(['C1P1_MaxDisturbance_zoom'], '-depsc', '-r900')
