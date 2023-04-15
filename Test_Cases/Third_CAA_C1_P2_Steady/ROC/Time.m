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

% Second_TriDi_On_01      = importdata('TriDi_SecOrder_D_10LHS_Diss_On__Converge0001.dat');
% Fourth_TriDi_On_01      = importdata('TriDi_FouOrder_D_10LHS_Diss_On__Converge0001.dat');
% Sixth_TriDi_On_01       = importdata('TriDi_SixOrder_D_10LHS_Diss_On__Converge0001.dat');
% RDRP_TriDi_On_01        = importdata('TriDi_RDRPSten_D_10LHS_Diss_On__Converge0001.dat');
% DRP_TriDi_On_01         = importdata('TriDi__DRPSten_D_10LHS_Diss_On__Converge0001.dat');
% Compact_TriDi_On_01     = importdata('TriDi_Compact4_D_10LHS_Diss_On__Converge0001.dat');
% PreFac4Compact_TriDi_On_01     = importdata('TriDi_PreFac_4_D_10LHS_Diss_On__Converge0001.dat');
% PreFac6Compact_TriDi_On_01     = importdata('TriDi_PreFac_6_D_10LHS_Diss_On__Converge0001.dat');
% 
% Second_TriDi_LUSGS      = importdata('LUSGS_SecOrder_D_10LHS_Diss_Off_Converge0001.dat');
% Fourth_TriDi_LUSGS      = importdata('LUSGS_FouOrder_D_10LHS_Diss_Off_Converge0001.dat');
% Sixth_TriDi_LUSGS       = importdata('LUSGS_SixOrder_D_10LHS_Diss_Off_Converge0001.dat');
% RDRP_TriDi_LUSGS        = importdata('LUSGS_RDRPSten_D_10LHS_Diss_Off_Converge0001.dat');
% DRP_TriDi_LUSGS         = importdata('LUSGS__DRPSten_D_10LHS_Diss_Off_Converge0001.dat');
% Compact_TriDi_LUSGS     = importdata('LUSGS_Compact4_D_10LHS_Diss_Off_Converge0001.dat');
% PreFac4Compact_LUSGS_Off01     = importdata('LUSGS_PreFac_4_D_10LHS_Diss_Off_Converge0001.dat');
% PreFac6Compact_LUSGS_Off01     = importdata('LUSGS_PreFac_6_D_10LHS_Diss_Off_Converge0001.dat');

% %% LUSGS
% figure(1)
%     semilogy(Second_TriDi_LUSGS(:, 3), Second_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(Fourth_TriDi_LUSGS(:, 3), Fourth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(Sixth_TriDi_LUSGS(:, 3), Sixth_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(RDRP_TriDi_LUSGS(:, 3), RDRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(DRP_TriDi_LUSGS(:, 3), DRP_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(Compact_TriDi_LUSGS(:, 3), Compact_TriDi_LUSGS(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(PreFac4Compact_LUSGS_Off01(:, 3), PreFac4Compact_LUSGS_Off01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(PreFac6Compact_LUSGS_Off01(:, 3), PreFac6Compact_LUSGS_Off01(:, 2), 'color', [0.25 0.25 0.25], 'LineWidth', 1.0)
%     hold on
%     grid on
%     grid minor
%     title('LU-SGS')
%     legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP',  ...
%             '4^{th} Compact',  '4^{th} PreFactored Compact', ...
%             '6^{th} PreFactored Compact', 'NumColumns',1, 'Location', 'best')  
%         xlabel('Computational Time (s)')
%     ylabel('Residual')
%     print(['LUSGS_Time_C1P2'], '-depsc', '-r900')
    
%% TriDi
figure(2)
    semilogy(Second_TriDi_Off01(:, 3), Second_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Fourth_TriDi_Off01(:, 3), Fourth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Sixth_TriDi_Off01(:, 3), Sixth_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(RDRP_TriDi_Off01(:, 3), RDRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(DRP_TriDi_Off01(:, 3), DRP_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(Compact_TriDi_Off01(:, 3), Compact_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(PreFac4Compact_TriDi_Off01(:, 3), PreFac4Compact_TriDi_Off01(:, 2),'LineWidth', 1.0)
    hold on
    semilogy(PreFac6Compact_TriDi_Off01(:, 3), PreFac6Compact_TriDi_Off01(:, 2), 'color', [0.25 0.25 0.25], 'LineWidth', 1.0)
    hold on
    grid on
    grid minor
    title('TriDi with Scaling Factors')
    legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP',  ...
            '4^{th} Compact',  '4^{th} PreFactored Compact', ...
            '6^{th} PreFactored Compact', 'NumColumns',1, 'Location', 'best')  
        xlabel('Computational Time (s)')
    ylabel('Residual')
    print(['TriDi_Scaling_Time_C1P2'], '-depsc', '-r900')
    
%    
% %% TriDi
% figure(3)
%     semilogy(Second_TriDi_On_01(:, 3), Second_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(Fourth_TriDi_On_01(:, 3), Fourth_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(Sixth_TriDi_On_01(:, 3), Sixth_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(RDRP_TriDi_On_01(:, 3), RDRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(DRP_TriDi_On_01(:, 3), DRP_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(Compact_TriDi_On_01(:, 3), Compact_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(PreFac4Compact_TriDi_On_01(:, 3), PreFac4Compact_TriDi_On_01(:, 2),'LineWidth', 1.0)
%     hold on
%     semilogy(PreFac6Compact_TriDi_On_01(:, 3), PreFac6Compact_TriDi_On_01(:, 2), 'color', [0.25 0.25 0.25], 'LineWidth', 1.0)
%     hold on
%     grid on
%     grid minor
%     title('TriDi with Dissipation')
%     legend('Second Order', 'Fourth Order', 'Sixth Order', 'RDRP', 'DRP',  ...
%             '4^{th} Compact',  '4^{th} PreFactored Compact', ...
%             '6^{th} PreFactored Compact', 'NumColumns',1, 'Location', 'best')  
%         xlabel('Computational Time (s)')
%     ylabel('Residual')
%     print(['TriDi_Dissipation_Time_C1P2'], '-depsc', '-r900')
    
