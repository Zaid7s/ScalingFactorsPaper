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
    loglog(Second_TriDi_On(:, 3)-Second_TriDi_On(1, 3)+1, Second_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    loglog(Second_TriDi_Off(:, 3)-Second_TriDi_Off(1, 3)+1, Second_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    loglog(Second_LUSGS(:, 3)-Second_LUSGS(1, 3)+1, Second_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Total # of Iterations') 
    ylabel('Residual')
%     print(['Second_LAE'], '-depsc', '-r900')
    
%% Fourth
figure(2)
    loglog(Fourth_TriDi_On(:, 3)-Fourth_TriDi_On(1, 3)+1, Fourth_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    loglog(Fourth_TriDi_Off(:, 3)-Fourth_TriDi_Off(1, 3)+1, Fourth_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    loglog(Fourth_LUSGS(:, 3)-Fourth_LUSGS(1, 3)+1, Fourth_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Total # of Iterations') 
    ylabel('Residual')
%     print(['Fourth_LAE'], '-depsc', '-r900')
    
%% Sixth
figure(3)
    loglog(Sixth_TriDi_On(:, 3)-Sixth_TriDi_On(1, 3)+1, Sixth_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    loglog(Sixth_TriDi_Off(:, 3)-Sixth_TriDi_Off(1, 3)+1, Sixth_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    loglog(Sixth_LUSGS(:, 3) - Sixth_LUSGS(1, 3) + 1, Sixth_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Total # of Iterations') 
    ylabel('Residual')
%     print(['Sixth_LAE'], '-depsc', '-r900')
  
%% RDRP
figure(4)
    loglog(RDRP_TriDi_On(:, 3) - RDRP_TriDi_On(1, 3) + 1, RDRP_TriDi_On(:, 4), 'LineWidth', 2.0)
    hold on
    loglog(RDRP_TriDi_Off(:, 3) - RDRP_TriDi_Off(1, 3) + 1, RDRP_TriDi_Off(:, 4),'LineWidth', 2.0)
    hold on
    loglog(RDRP_LUSGS(:, 3) - RDRP_LUSGS(1, 3) + 1, RDRP_LUSGS(:, 4),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    legend('TriDi Dissipation', 'TriDi Scaling Factors', 'LU-SGS', 'Location', 'Best')
    xlabel('Total # of Iterations') 
    ylabel('Residual')
%     print(['RDRP_LAE'], '-depsc', '-r900')
% 
