clc
clear all

Second_TriDi_Off01    = importdata('TriDi_SecOrder_D_10LHS_Diss_Off_Converge.dat');
Fourth_TriDi_Off01    = importdata('TriDi_FouOrder_D_10LHS_Diss_Off_Converge.dat');
Sixth_TriDi_Off01     = importdata('TriDi_SixOrder_D_10LHS_Diss_Off_Converge.dat');
RDRP_TriDi_Off01      = importdata('TriDi_RDRPSten_D_10LHS_Diss_Off_Converge.dat');
DRP_TriDi_Off01      = importdata('TriDi__DRPSten_D_10LHS_Diss_Off_Converge.dat');

Second_TriDi_On_01    = importdata('TriDi_SecOrder_D_10LHS_Diss_On__Converge.dat');
Fourth_TriDi_On_01    = importdata('TriDi_FouOrder_D_10LHS_Diss_On__Converge.dat');
Sixth_TriDi_On_01     = importdata('TriDi_SixOrder_D_10LHS_Diss_On__Converge.dat');
RDRP_TriDi_On_01      = importdata('TriDi_RDRPSten_D_10LHS_Diss_On__Converge.dat');
DRP_TriDi_On_01      = importdata('TriDi__DRPSten_D_10LHS_Diss_On__Converge.dat');

Second_TriDi_LUSGS    = importdata('LUSGS_SecOrder_D_10LHS_Diss_Off_Converge.dat');
Fourth_TriDi_LUSGS    = importdata('LUSGS_FouOrder_D_10LHS_Diss_Off_Converge.dat');
Sixth_TriDi_LUSGS     = importdata('LUSGS_SixOrder_D_10LHS_Diss_Off_Converge.dat');
RDRP_TriDi_LUSGS      = importdata('LUSGS_RDRPSten_D_10LHS_Diss_Off_Converge.dat');
DRP_TriDi_LUSGS      = importdata('LUSGS__DRPSten_D_10LHS_Diss_Off_Converge.dat');

%% LUSGS
figure(1)
    loglog(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    loglog(Fourth_TriDi_LUSGS(:, 1), Fourth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    loglog(Sixth_TriDi_LUSGS(:, 1), Sixth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    loglog(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    loglog(DRP_TriDi_LUSGS(:, 1), DRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    title('LUSGS')
    legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP',  'DRP', 'Location', 'Best')
    xlabel('Total # of Iterations')
    ylabel('Residual')
    print(['LUSGS_C1P1'], '-depsc', '-r900')
    
%% TriDi
figure(2)
    loglog(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(Fourth_TriDi_Off01(:, 1), Fourth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(Sixth_TriDi_Off01(:, 1), Sixth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(DRP_TriDi_Off01(:, 1), DRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    title('TriDi Scaling Factors')
    legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP', 'Location', 'Best')
    xlabel('Total # of Iterations')
    ylabel('Residual')
    print(['TriDi_C1P1'], '-depsc', '-r900')
    
%% TriDi + Dissipation
figure(3)
    loglog(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(Fourth_TriDi_On_01(:, 1), Fourth_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(Sixth_TriDi_On_01(:, 1), Sixth_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    loglog(DRP_TriDi_On_01(:, 1), DRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    title('TriDi Dissipation')
    legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP', 'Location', 'Best')
    xlabel('Total # of Iterations')
    ylabel('Residual')
    print(['TriDi_Dissipation_C1P1'], '-depsc', '-r900')
