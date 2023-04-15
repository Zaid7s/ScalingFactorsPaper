clc
clear all

Second_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_PvsTime.dat');
Fourth_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_FouOrder_D_10_PvsTime.dat');
Sixth_TriDi_LUSGS     = importdata('LUSGS_LHS_Diss_Off_SixOrder_D_10_PvsTime.dat');
RDRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_RDRPSten_D_10_PvsTime.dat');
DRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off__DRPSten_D_10_PvsTime.dat');
Compact_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_Compact4_D_10_PvsTime.dat');
% PreFac4_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_Compact4_D_10_PvsTime.dat');
% PreFac6_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_PreFac_6_D_10_PvsTime.dat');

Exact       = importdata('CAT1P2.FIG4.DATA'); 
 
%% LUSGS Scaling Factors
figure(3)
    plot(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'-o','MarkerIndices',1:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Fourth_TriDi_LUSGS(:, 1), Fourth_TriDi_LUSGS(:, 2),'-x','MarkerIndices',2:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Sixth_TriDi_LUSGS(:, 1), Sixth_TriDi_LUSGS(:, 2),'-*','MarkerIndices',3:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'-s','MarkerIndices',4:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(DRP_TriDi_LUSGS(:, 1), DRP_TriDi_LUSGS(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Compact_TriDi_LUSGS(:, 1), Compact_TriDi_LUSGS(:, 2),'-d','MarkerIndices',1:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:5:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', 'Exact', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    xlabel('Time', 'Fontsize', 14)
    ylabel('Perturbation Pressure', 'Fontsize', 14)
    print(['C1P2_P_vs_Time'], '-depsc', '-r900')  
    
