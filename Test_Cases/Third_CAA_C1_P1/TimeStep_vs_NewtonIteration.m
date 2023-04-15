clc
clear all
clf 

Second_TriDi_Off01    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
Fourth_TriDi_Off01    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
Sixth_TriDi_Off01     = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
RDRP_TriDi_Off01      = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
DRP_TriDi_Off01      = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');

Second_TriDi_On_01    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
Fourth_TriDi_On_01    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
Sixth_TriDi_On_01     = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
RDRP_TriDi_On_01      = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
DRP_TriDi_On_01      = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');

Second_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_Time.dat');
Fourth_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_FouOrder_D_10_Time.dat');
Sixth_TriDi_LUSGS     = importdata('LUSGS_LHS_Diss_Off_SixOrder_D_10_Time.dat');
RDRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_RDRPSten_D_10_Time.dat');
DRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off__DRPSten_D_10_Time.dat');

%% Maximum
max1(1, 1)  = max(Second_TriDi_Off01(:, 2));
max1(2, 1)  = max(Fourth_TriDi_Off01(:, 2));
max1(3, 1)  = max(Sixth_TriDi_Off01(:, 2));
max1(4, 1)  = max(RDRP_TriDi_Off01(:, 2));
max1(5, 1)  = max(DRP_TriDi_Off01(:, 2));

max2(1, 1) = max(Second_TriDi_On_01(:, 2));
max2(2, 1) = max(Fourth_TriDi_On_01(:, 2));
max2(3, 1) = max(Sixth_TriDi_On_01(:, 2));
max2(4, 1) = max(RDRP_TriDi_On_01(:, 2));
max2(5, 1) = max(DRP_TriDi_On_01(:, 2));

max3(1, 1) = max(Second_TriDi_LUSGS(:, 2));
max3(2, 1) = max(Fourth_TriDi_LUSGS(:, 2));
max3(3, 1) = max(Sixth_TriDi_LUSGS(:, 2));
max3(4, 1) = max(RDRP_TriDi_LUSGS(:, 2));
max3(5, 1) = max(DRP_TriDi_LUSGS(:, 2));

%% Minumum
min1(1, 1)  = min(Second_TriDi_Off01(:, 2));
min1(2, 1)  = min(Fourth_TriDi_Off01(:, 2));
min1(3, 1)  = min(Sixth_TriDi_Off01(:, 2));
min1(4, 1)  = min(RDRP_TriDi_Off01(:, 2));
min1(5, 1)  = min(DRP_TriDi_Off01(:, 2));

min2(1, 1) = min(Second_TriDi_On_01(:, 2));
min2(2, 1) = min(Fourth_TriDi_On_01(:, 2));
min2(3, 1) = min(Sixth_TriDi_On_01(:, 2));
min2(4, 1) = min(RDRP_TriDi_On_01(:, 2));
min2(5, 1) = min(DRP_TriDi_On_01(:, 2));

min3(1, 1) = min(Second_TriDi_LUSGS(:, 2));
min3(2, 1) = min(Fourth_TriDi_LUSGS(:, 2));
min3(3, 1) = min(Sixth_TriDi_LUSGS(:, 2));
min3(4, 1) = min(RDRP_TriDi_LUSGS(:, 2));
min3(5, 1) = min(DRP_TriDi_LUSGS(:, 2));

%% Average
mean1(1, 1)  = mean(Second_TriDi_Off01(:, 2));
mean1(2, 1)  = mean(Fourth_TriDi_Off01(:, 2));
mean1(3, 1)  = mean(Sixth_TriDi_Off01(:, 2));
mean1(4, 1)  = mean(RDRP_TriDi_Off01(:, 2));
mean1(5, 1)  = mean(DRP_TriDi_Off01(:, 2));

mean2(1, 1) = mean(Second_TriDi_On_01(:, 2));
mean2(2, 1) = mean(Fourth_TriDi_On_01(:, 2));
mean2(3, 1) = mean(Sixth_TriDi_On_01(:, 2));
mean2(4, 1) = mean(RDRP_TriDi_On_01(:, 2));
mean2(5, 1) = mean(DRP_TriDi_On_01(:, 2));

mean3(1, 1) = mean(Second_TriDi_LUSGS(:, 2));
mean3(2, 1) = mean(Fourth_TriDi_LUSGS(:, 2));
mean3(3, 1) = mean(Sixth_TriDi_LUSGS(:, 2));
mean3(4, 1) = mean(RDRP_TriDi_LUSGS(:, 2));
mean3(5, 1) = mean(DRP_TriDi_LUSGS(:, 2));

%% Sum
sum1(1, 1)  = sum(Second_TriDi_Off01(:, 2));
sum1(2, 1)  = sum(Fourth_TriDi_Off01(:, 2));
sum1(3, 1)  = sum(Sixth_TriDi_Off01(:, 2));
sum1(4, 1)  = sum(RDRP_TriDi_Off01(:, 2));
sum1(5, 1)  = sum(DRP_TriDi_Off01(:, 2));

sum2(1, 1) = sum(Second_TriDi_On_01(:, 2));
sum2(2, 1) = sum(Fourth_TriDi_On_01(:, 2));
sum2(3, 1) = sum(Sixth_TriDi_On_01(:, 2));
sum2(4, 1) = sum(RDRP_TriDi_On_01(:, 2));
sum2(5, 1) = sum(DRP_TriDi_On_01(:, 2));

sum3(1, 1) = sum(Second_TriDi_LUSGS(:, 2));
sum3(2, 1) = sum(Fourth_TriDi_LUSGS(:, 2));
sum3(3, 1) = sum(Sixth_TriDi_LUSGS(:, 2));
sum3(4, 1) = sum(RDRP_TriDi_LUSGS(:, 2));
sum3(5, 1) = sum(DRP_TriDi_LUSGS(:, 2));

%% Create Table
% TriDi Scale
table(max1, min1, mean1, sum1)

% TriDi Dissipation
table(max2, min2, mean2, sum2)

% LUSGS
table(max3, min3, mean3, sum3)

%% TriDi Scaling Factors
% figure(1)
    subplot(3, 1,  1)
    semilogx(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Fourth_TriDi_Off01(:, 1), Fourth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Sixth_TriDi_Off01(:, 1), Sixth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(DRP_TriDi_Off01(:, 1), DRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    %ylim([0 65])
    title('TriDi Scaling Factors')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP','Location', 'Best')
    xlabel('Time Step')
    ylabel('Newton Iteration')
%     print(['TriDi_ScalingFactors_C1P1_Time'], '-depsc', '-r900')
 
%% TriDi Dissipation
% figure(2)
    subplot(3, 1,  2)
    semilogx(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Fourth_TriDi_On_01(:, 1), Fourth_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Sixth_TriDi_On_01(:, 1), Sixth_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(DRP_TriDi_On_01(:, 1), DRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    %ylim([0 65])
    title('TriDi Dissipation')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP','Location', 'Best')
    xlabel('Time Step')
    ylabel('Newton Iteration')
%     print(['TriDi_Dissipation_Time'], '-depsc', '-r900')
 
%% LUSGS
% figure(3)
    subplot(3, 1,  3)
    semilogx(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Fourth_TriDi_LUSGS(:, 1), Fourth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Sixth_TriDi_LUSGS(:, 1), Sixth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(DRP_TriDi_LUSGS(:, 1), DRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    %ylim([0 65])
    title('LUSGS')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP','Location', 'Best')
    xlabel('Time Step')
    ylabel('Newton Iteration')
    set(gcf, 'Position', [1466 237 560 697])
    
    print(['C1P1_Iteration'], '-depsc', '-r900')
%% All
figure(2)
    semilogx(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
    hold on
    semilogx(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    %ylim([0 65])
    title('LUSGS')
    legend('TriDi Scaling Factor, E2', 'TriDi Dissipation, E2', 'LU-SGS, E2', 'TriDi Scaling Factor, RDRP', 'TriDi Dissipation, RDRP', 'LU-SGS, RDRP','Exact','Location', 'Best')
    xlabel('Time Step')
    ylabel('Newton Iteration')
%     set(gcf, 'Position', [1466 237 560 697])
    
    print(['C1P1_Iteration_All'], '-depsc', '-r900')
    

