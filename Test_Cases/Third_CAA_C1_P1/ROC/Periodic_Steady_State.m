clc
clear all

Second_TriDi_Off01    = importdata('TriDi_LHS_Diss_Off_SecOrder_D_10_128_Periodic_Steady_State.dat');
Fourth_TriDi_Off01    = importdata('TriDi_LHS_Diss_Off_FouOrder_D_10_128_Periodic_Steady_State.dat');
Sixth_TriDi_Off01     = importdata('TriDi_LHS_Diss_Off_SixOrder_D_10_128_Periodic_Steady_State.dat');
RDRP_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off_RDRPSten_D_10_128_Periodic_Steady_State.dat');
DRP_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off__DRPSten_D_10_128_Periodic_Steady_State.dat');
Compact_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off_Compact4_D_10_128_Periodic_Steady_State.dat');

Second_TriDi_On_01    = importdata('TriDi_LHS_Diss_On__SecOrder_D_10_128_Periodic_Steady_State.dat');
Fourth_TriDi_On_01    = importdata('TriDi_LHS_Diss_On__FouOrder_D_10_128_Periodic_Steady_State.dat');
Sixth_TriDi_On_01     = importdata('TriDi_LHS_Diss_On__SixOrder_D_10_128_Periodic_Steady_State.dat');
RDRP_TriDi_On_01      = importdata('TriDi_LHS_Diss_On__RDRPSten_D_10_128_Periodic_Steady_State.dat');
DRP_TriDi_On_01      = importdata('TriDi_LHS_Diss_On___DRPSten_D_10_128_Periodic_Steady_State.dat');
Compact_TriDi_On_01      = importdata('TriDi_LHS_Diss_Off_Compact4_D_10_128_Periodic_Steady_State.dat');

Second_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_128_Periodic_Steady_State.dat');
Fourth_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_FouOrder_D_10_128_Periodic_Steady_State.dat');
Sixth_TriDi_LUSGS     = importdata('LUSGS_LHS_Diss_Off_SixOrder_D_10_128_Periodic_Steady_State.dat');
RDRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_RDRPSten_D_10_128_Periodic_Steady_State.dat');
DRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off__DRPSten_D_10_128_Periodic_Steady_State.dat');
Compact_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_Compact4_D_10_128_Periodic_Steady_State.dat');

%% TriDi Scaling Factors
figure(1)
    semilogy(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'-o','MarkerIndices',1:5:length(Second_TriDi_Off01),'LineWidth', 2.0)
    hold on
    semilogy(Fourth_TriDi_Off01(:, 1), Fourth_TriDi_Off01(:, 2),'-x','MarkerIndices',2:5:length(Second_TriDi_Off01),'LineWidth', 2.0)
    hold on
    semilogy(Sixth_TriDi_Off01(:, 1), Sixth_TriDi_Off01(:, 2),'-*','MarkerIndices',3:5:length(Second_TriDi_Off01),'LineWidth', 2.0)
    hold on
    semilogy(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'-s','MarkerIndices',4:5:length(Second_TriDi_Off01),'LineWidth', 2.0)
    hold on
    semilogy(DRP_TriDi_Off01(:, 1), DRP_TriDi_Off01(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_Off01),'LineWidth', 2.0)
    hold on
    semilogy(Compact_TriDi_Off01(:, 1), Compact_TriDi_Off01(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_Off01),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    title('TriDi Scaling Factors')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP','4^{th} Compact', 'Location', 'Best')
    xlabel('Cycle Number')
    ylabel('Residual')
    print(['TriDi_ScalingFactors_C1P2_PeriodicSteadyState'], '-depsc', '-r900')
 
%% TriDi Dissipation
figure(2)
    semilogy(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'-o','MarkerIndices',1:5:length(Second_TriDi_On_01),'LineWidth', 2.0)
    hold on
    semilogy(Fourth_TriDi_On_01(:, 1), Fourth_TriDi_On_01(:, 2),'-x','MarkerIndices',2:5:length(Second_TriDi_On_01),'LineWidth', 2.0)
    hold on
    semilogy(Sixth_TriDi_On_01(:, 1), Sixth_TriDi_On_01(:, 2),'-*','MarkerIndices',3:5:length(Second_TriDi_On_01),'LineWidth', 2.0)
    hold on
    semilogy(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'-s','MarkerIndices',4:5:length(Second_TriDi_On_01),'LineWidth', 2.0)
    hold on
    semilogy(DRP_TriDi_On_01(:, 1), DRP_TriDi_On_01(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_On_01),'LineWidth', 2.0)
    hold on
    semilogy(Compact_TriDi_On_01(:, 1), Compact_TriDi_On_01(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_On_01),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    title('TriDi Dissipation')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP','4^{th} Compact', 'Location', 'Best')
    xlabel('Cycle Number')
    ylabel('Residual')
    print(['TriDi_Dissipation_C1P2_PeriodicSteadyState'], '-depsc', '-r900')
 
%% LUSGS
figure(3)
    semilogy(Second_TriDi_LUSGS(:, 1), Second_TriDi_LUSGS(:, 2),'-o','MarkerIndices',1:5:length(Second_TriDi_LUSGS),'LineWidth', 2.0)
    hold on
    semilogy(Fourth_TriDi_LUSGS(:, 1), Fourth_TriDi_LUSGS(:, 2),'-x','MarkerIndices',2:5:length(Second_TriDi_LUSGS),'LineWidth', 2.0)
    hold on
    semilogy(Sixth_TriDi_LUSGS(:, 1), Sixth_TriDi_LUSGS(:, 2),'-*','MarkerIndices',3:5:length(Second_TriDi_LUSGS),'LineWidth', 2.0)
    hold on
    semilogy(RDRP_TriDi_LUSGS(:, 1), RDRP_TriDi_LUSGS(:, 2),'-s','MarkerIndices',4:5:length(Second_TriDi_LUSGS),'LineWidth', 2.0)
    hold on
    semilogy(DRP_TriDi_LUSGS(:, 1), DRP_TriDi_LUSGS(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_LUSGS),'LineWidth', 2.0)
    hold on
    semilogy(Compact_TriDi_LUSGS(:, 1), Compact_TriDi_LUSGS(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_LUSGS),'LineWidth', 2.0)
    hold on
    grid on
    grid minor
    title('LUSGS')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP','4^{th} Compact', 'Location', 'Best')
    xlabel('Cycle Number')
    ylabel('Residual')
    print(['LUSGS_C1P2_PeriodicSteadyState'], '-depsc', '-r900')
