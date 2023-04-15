clc
clear all

Second_TriDi_Off    = importdata('TriDi_SecOrder_D_10LHS_Diss_Off_Converge.dat');
Fourth_TriDi_Off    = importdata('TriDi_FouOrder_D_10LHS_Diss_Off_Converge.dat');
Sixth_TriDi_Off     = importdata('TriDi_SixOrder_D_10LHS_Diss_Off_Converge.dat');
RDRP_TriDi_Off      = importdata('TriDi_RDRPSten_D_10LHS_Diss_Off_Converge.dat');


Second_TriDi_On    = importdata('TriDi_SecOrder_D_10LHS_Diss_On__Converge.dat');
Fourth_TriDi_On    = importdata('TriDi_FouOrder_D_10LHS_Diss_On__Converge.dat');
Sixth_TriDi_On     = importdata('TriDi_SixOrder_D_10LHS_Diss_On__Converge.dat');
RDRP_TriDi_On      = importdata('TriDi_RDRPSten_D_10LHS_Diss_On__Converge.dat');

Second_LUSGS    = importdata('LUSGS_SecOrder_D_10LHS_Diss_Off_Converge.dat');
Fourth_LUSGS    = importdata('LUSGS_FouOrder_D_10LHS_Diss_Off_Converge.dat');
Sixth_LUSGS     = importdata('LUSGS_SixOrder_D_10LHS_Diss_Off_Converge.dat');
RDRP_LUSGS      = importdata('LUSGS_RDRPSten_D_10LHS_Diss_Off_Converge.dat');

%% Second
figure(1)
    semilogy(Second_TriDi_On(:, 3), Second_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    semilogy(Second_TriDi_Off(:, 3), Second_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    semilogy(Second_LUSGS(:, 3), Second_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Computational Work')
    ylabel('Residual')
    print(['Second'], '-depsc', '-r900')
    
%% Fourth
figure(2)
    semilogy(Fourth_TriDi_On(:, 3), Fourth_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    semilogy(Fourth_TriDi_Off(:, 3), Fourth_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    semilogy(Fourth_LUSGS(:, 3), Fourth_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Computational Work')
    ylabel('Residual')
    print(['Fourth'], '-depsc', '-r900')
    
%% Sixth
figure(3)
    semilogy(Sixth_TriDi_On(:, 3), Sixth_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    semilogy(Sixth_TriDi_Off(:, 3), Sixth_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    semilogy(Sixth_LUSGS(:, 3), Sixth_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Computational Work')
    ylabel('Residual')
    print(['Sixth'], '-depsc', '-r900')
  
%% RDRP
figure(4)
    semilogy(RDRP_TriDi_On(:, 3), RDRP_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    semilogy(RDRP_TriDi_Off(:, 3), RDRP_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    semilogy(RDRP_LUSGS(:, 3), RDRP_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Computational Work')
    ylabel('Residual')
    print(['RDRP'], '-depsc', '-r900')
