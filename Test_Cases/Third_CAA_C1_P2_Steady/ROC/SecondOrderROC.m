clc
clear all

Second_TriDi_Off01      = importdata('TriDi_SecOrder_D_10LHS_Diss_Off_Converge0001.dat');
Fourth_TriDi_Off01      = importdata('TriDi_FouOrder_D_10LHS_Diss_Off_Converge0001.dat');
Sixth_TriDi_Off01       = importdata('TriDi_SixOrder_D_10LHS_Diss_Off_Converge0001.dat');
RDRP_TriDi_Off01        = importdata('TriDi_RDRPSten_D_10LHS_Diss_Off_Converge0001.dat');
DRP_TriDi_Off01         = importdata('TriDi__DRPSten_D_10LHS_Diss_Off_Converge0001.dat');
Compact_TriDi_Off01     = importdata('TriDi_Compact4_D_10LHS_Diss_Off_Converge0001.dat');
PreFac4Compact_TriDi_Off01     = importdata('TriDi_PreFac_4_D_10LHS_Diss_Off_Converge0001.dat');
PreFac6Compact_TriDi_Off01     = importdata('TriDi_PreFac_6_D_10LHS_Diss_Off_Converge0001.dat');

Second_TriDi_On_01      = importdata('TriDi_SecOrder_D_10LHS_Diss_On__Converge0001.dat');
Fourth_TriDi_On_01      = importdata('TriDi_FouOrder_D_10LHS_Diss_On__Converge0001.dat');
Sixth_TriDi_On_01       = importdata('TriDi_SixOrder_D_10LHS_Diss_On__Converge0001.dat');
RDRP_TriDi_On_01        = importdata('TriDi_RDRPSten_D_10LHS_Diss_On__Converge0001.dat');
DRP_TriDi_On_01         = importdata('TriDi__DRPSten_D_10LHS_Diss_On__Converge0001.dat');
Compact_TriDi_On_01     = importdata('TriDi_Compact4_D_10LHS_Diss_On__Converge0001.dat');
PreFac4Compact_TriDi_On_01     = importdata('TriDi_PreFac_4_D_10LHS_Diss_On__Converge0001.dat');
PreFac6Compact_TriDi_On_01     = importdata('TriDi_PreFac_6_D_10LHS_Diss_On__Converge0001.dat');

Second_TriDi_LUSGS      = importdata('LUSGS_SecOrder_D_10LHS_Diss_Off_Converge0001.dat');
Fourth_TriDi_LUSGS      = importdata('LUSGS_FouOrder_D_10LHS_Diss_Off_Converge0001.dat');
Sixth_TriDi_LUSGS       = importdata('LUSGS_SixOrder_D_10LHS_Diss_Off_Converge0001.dat');
RDRP_TriDi_LUSGS        = importdata('LUSGS_RDRPSten_D_10LHS_Diss_Off_Converge0001.dat');
DRP_TriDi_LUSGS         = importdata('LUSGS__DRPSten_D_10LHS_Diss_Off_Converge0001.dat');
Compact_TriDi_LUSGS     = importdata('LUSGS_Compact4_D_10LHS_Diss_Off_Converge0001.dat');
PreFac4Compact_LUSGS_Off01     = importdata('LUSGS_PreFac_4_D_10LHS_Diss_Off_Converge0001.dat');
PreFac6Compact_LUSGS_Off01     = importdata('LUSGS_PreFac_6_D_10LHS_Diss_Off_Converge0001.dat');

%% LUSGS
figure(1)
    semilogy(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Fourth_TriDi_LUSGS(:, 1), Fourth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Sixth_TriDi_LUSGS(:, 1), Sixth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(DRP_TriDi_LUSGS(:, 1), DRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(PreFac4Compact_LUSGS_Off01(:, 1), PreFac4Compact_LUSGS_Off01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP',  ...
            'Prefactored Fourth-Order Compact ', ...
             'NumColumns',1, 'Location', 'best')  
        xlabel('Iteration')
    ylabel('Residual')
   print(['C1P2_LUSGS_ROC'], '-depsc', '-r900')
    
%% TriDi
figure(2)
    semilogy(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Fourth_TriDi_Off01(:, 1), Fourth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Sixth_TriDi_Off01(:, 1), Sixth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(DRP_TriDi_Off01(:, 1), DRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(PreFac4Compact_TriDi_Off01(:, 1), PreFac4Compact_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    xlabel('Iteration', 'Fontsize', 14)
    ylabel('Residual', 'Fontsize', 14)
   print(['C1P2_TriDi_Scaling_ROC'], '-depsc', '-r900')
    
   
%% TriDi
figure(3)
    semilogy(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Fourth_TriDi_On_01(:, 1), Fourth_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Sixth_TriDi_On_01(:, 1), Sixth_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(DRP_TriDi_On_01(:, 1), DRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(PreFac4Compact_TriDi_On_01(:, 1), PreFac4Compact_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12; 
    legend('E2', 'E4', 'E6', 'RDRP', 'DRP',  ...
            'C4', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    xlabel('Iteration', 'Fontsize', 14)
    ylabel('Residual', 'Fontsize', 14)
   print(['C1P2_TriDi_Dissipation_ROC'], '-depsc', '-r900')
    
% Second vs RDRP
figure(5)
    semilogy(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12; 
    legend('E2 Scaling Factors', 'E2 Dissipation', 'E2 LU-SGS', ...
           'RDRP Scaling Factors', 'RDRP Dissipation', 'RDRP LU-SGS', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    xlabel('Iteration', 'Fontsize', 14)
    ylabel('Residual', 'Fontsize', 14)
   print(['C1P2_ROC_E2_vs_RDRP'], '-depsc', '-r900')
  
%% TriDi
figure(4)
    semilogy(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12; 
    legend('Scaling Factors', 'LU-SGS', 'Dissipation', ...
            'NumColumns',1, 'Location', 'best', 'Fontsize', 12) 
    xlabel('Iteration', 'Fontsize', 14)
    ylabel('Residual', 'Fontsize', 14)
    print(['RDRP_ROC_C1P2'], '-depsc', '-r900')
     
