clc
clear all

Second_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_072_Max_Min_Pres.dat');
Fourth_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_FouOrder_D_10_072_Max_Min_Pres.dat');
Sixth_TriDi_LUSGS       = importdata('LUSGS_LHS_Diss_Off_SixOrder_D_10_072_Max_Min_Pres.dat');
RDRP_TriDi_LUSGS        = importdata('LUSGS_LHS_Diss_Off_RDRPSten_D_10_072_Max_Min_Pres.dat');
DRP_TriDi_LUSGS         = importdata('LUSGS_LHS_Diss_Off__DRPSten_D_10_072_Max_Min_Pres.dat');
Compact_TriDi_LUSGS     = importdata('LUSGS_LHS_Diss_Off_Compact4_D_10_072_Max_Min_Pres.dat');

Exact       = importdata('CAT1P2.FIG3.DATA'); 
 
%% LUSGS Scaling Factors
figure(3)
    plot(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 6),'-o','MarkerIndices',1:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Fourth_TriDi_LUSGS(:, 1), Fourth_TriDi_LUSGS(:, 6),'-x','MarkerIndices',2:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Sixth_TriDi_LUSGS(:, 1), Sixth_TriDi_LUSGS(:, 6),'-*','MarkerIndices',3:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 6),'-s','MarkerIndices',4:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(DRP_TriDi_LUSGS(:, 1), DRP_TriDi_LUSGS(:, 6),'-h','MarkerIndices',5:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Compact_TriDi_LUSGS(:, 1), Compact_TriDi_LUSGS(:, 6),'-d','MarkerIndices',1:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:5:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
    xlim([-10 10])
    ylim([-2e-5 2e-5])
    legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP',  ...
            'Prefactored Fourth-Order Compact', 'Exact', ...
            'NumColumns',1, 'Location', 'best') 
    xlabel('Domain')
    ylabel('Perturbation Pressure')
%     print(['C1P2_P_vs_Time'], '-depsc', '-r900')
  
