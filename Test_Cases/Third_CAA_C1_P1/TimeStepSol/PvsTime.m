clc
clear all

Second_TriDi_Off01    = importdata('TriDi_LHS_Diss_Off_SecOrder_D_10_PvsTime.dat');
Fourth_TriDi_Off01    = importdata('TriDi_LHS_Diss_Off_FouOrder_D_10_PvsTime.dat');
Sixth_TriDi_Off01     = importdata('TriDi_LHS_Diss_Off_SixOrder_D_10_PvsTime.dat');
RDRP_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off_RDRPSten_D_10_PvsTime.dat');
DRP_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off__DRPSten_D_10_PvsTime.dat');
Compact_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off_Compact4_D_10_PvsTime.dat');
PreFac4_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off_PreFac_4_D_10_PvsTime.dat');
PreFac6_TriDi_Off01      = importdata('TriDi_LHS_Diss_Off_PreFac_6_D_10_PvsTime.dat');

Second_TriDi_On_01    = importdata('TriDi_LHS_Diss_On__SecOrder_D_10_PvsTime.dat');
Fourth_TriDi_On_01    = importdata('TriDi_LHS_Diss_On__FouOrder_D_10_PvsTime.dat');
Sixth_TriDi_On_01     = importdata('TriDi_LHS_Diss_On__SixOrder_D_10_PvsTime.dat');
RDRP_TriDi_On_01      = importdata('TriDi_LHS_Diss_On__RDRPSten_D_10_PvsTime.dat');
DRP_TriDi_On_01      = importdata('TriDi_LHS_Diss_On___DRPSten_D_10_PvsTime.dat');
Compact_TriDi_On_01      = importdata('TriDi_LHS_Diss_On__Compact4_D_10_PvsTime.dat');
PreFac4_TriDi_On_01      = importdata('TriDi_LHS_Diss_On__PreFac_4_D_10_PvsTime.dat');
PreFac6_TriDi_On_01      = importdata('TriDi_LHS_Diss_On__PreFac_6_D_10_PvsTime.dat');

Second_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_SecOrder_D_10_PvsTime.dat');
Fourth_TriDi_LUSGS    = importdata('LUSGS_LHS_Diss_Off_FouOrder_D_10_PvsTime.dat');
Sixth_TriDi_LUSGS     = importdata('LUSGS_LHS_Diss_Off_SixOrder_D_10_PvsTime.dat');
RDRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_RDRPSten_D_10_PvsTime.dat');
DRP_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off__DRPSten_D_10_PvsTime.dat');
Compact_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_Compact4_D_10_PvsTime.dat');
PreFac4_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_PreFac_4_D_10_PvsTime.dat');
PreFac6_TriDi_LUSGS      = importdata('LUSGS_LHS_Diss_Off_PreFac_6_D_10_PvsTime.dat');

Exact       = importdata('CAT1P2.FIG4.DATA'); 

%% TriDi Scaling Factors
figure(1)
    plot(Second_TriDi_Off01(:, 1), Second_TriDi_Off01(:, 2),'-o','MarkerIndices',1:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(Fourth_TriDi_Off01(:, 1), Fourth_TriDi_Off01(:, 2),'-x','MarkerIndices',2:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(Sixth_TriDi_Off01(:, 1), Sixth_TriDi_Off01(:, 2),'-*','MarkerIndices',3:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(RDRP_TriDi_Off01(:, 1), RDRP_TriDi_Off01(:, 2),'-s','MarkerIndices',4:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(DRP_TriDi_Off01(:, 1), DRP_TriDi_Off01(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(Compact_TriDi_Off01(:, 1), Compact_TriDi_Off01(:, 2),'-d','MarkerIndices',1:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(PreFac4_TriDi_Off01(:, 1), PreFac4_TriDi_Off01(:, 2),'-.','MarkerIndices',2:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(PreFac6_TriDi_Off01(:, 1), PreFac6_TriDi_Off01(:, 2),'-p','MarkerIndices',3:5:length(Second_TriDi_Off01),'LineWidth', 1.0)
    hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:5:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
%     title('TriDi Scaling Factors')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP', ...
            '4^{th} Compact',  '4^{th} Prefactored Compact',  ...
            '6^{th} Prefactored Compact', 'Exact','Location', 'Best')
    xlabel('Time')
    ylabel('Perturbation Pressure')
    print(['TriDi_ScalingFactors_C1P2'], '-depsc', '-r900')
 
%% TriDi Dissipation
figure(2)
    plot(Second_TriDi_On_01(:, 1), Second_TriDi_On_01(:, 2),'-o','MarkerIndices',1:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(Fourth_TriDi_On_01(:, 1), Fourth_TriDi_On_01(:, 2),'-x','MarkerIndices',2:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(Sixth_TriDi_On_01(:, 1), Sixth_TriDi_On_01(:, 2),'-*','MarkerIndices',3:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(RDRP_TriDi_On_01(:, 1), RDRP_TriDi_On_01(:, 2),'-s','MarkerIndices',4:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(DRP_TriDi_On_01(:, 1), DRP_TriDi_On_01(:, 2),'-h','MarkerIndices',5:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(Compact_TriDi_On_01(:, 1), Compact_TriDi_On_01(:, 2),'-d','MarkerIndices',1:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(PreFac4_TriDi_On_01(:, 1), PreFac4_TriDi_On_01(:, 2),'-.','MarkerIndices',2:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(PreFac6_TriDi_On_01(:, 1), PreFac6_TriDi_On_01(:, 2),'-p','MarkerIndices',3:5:length(Second_TriDi_On_01),'LineWidth', 1.0)
    hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:5:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
%     title('TriDi Scaling Factors')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP', ...
            '4^{th} Compact',  '4^{th} Prefactored Compact',  ...
            '6^{th} Prefactored Compact', 'Exact','Location', 'Best')
    xlabel('Time')
    ylabel('Perturbation Pressure')
    print(['TriDi_Dissipation_C1P2'], '-depsc', '-r900')
    
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
    plot(PreFac4_TriDi_LUSGS(:, 1), PreFac4_TriDi_LUSGS(:, 2),'-.','MarkerIndices',2:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(PreFac6_TriDi_LUSGS(:, 1), PreFac6_TriDi_LUSGS(:, 2),'-p','MarkerIndices',3:5:length(Second_TriDi_LUSGS),'LineWidth', 1.0)
    hold on
    plot(Exact(:, 1), Exact(:, 2),'-','MarkerIndices',20:5:length(Exact),'LineWidth', 1.5, 'Color', [0.25 0.25 0.25])
    hold on
    grid on
    grid minor
%     title('TriDi Scaling Factors')
    legend('2^{nd} Order', '4^{th} Order', '6^{th} Order', 'RDRP', 'DRP', ...
            '4^{th} Compact',  '4^{th} Prefactored Compact',  ...
            '6^{th} Prefactored Compact', 'Exact','Location', 'Best')
    xlabel('Time')
    ylabel('Perturbation Pressure')
    print(['LUSGS_C1P2'], '-depsc', '-r900')
    
