MODULE GetLUSGS
!!! Working Test Case
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetRHSVector
    USE GetRungeKutta
    USE GetLMatrix
    USE GetUMatrix
    USE GetDMatrix
    USE GetSolveLUSGS
    USE GetLDU_Vectors
    USE ForwardSweep
    USE BackwardSweep

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: LUSGS

    INTERFACE LUSGS
        MODULE PROCEDURE LUSGS_Mean
        MODULE PROCEDURE LUSGS_Mean_Pet
    END INTERFACE LUSGS

    CONTAINS

    SUBROUTINE LUSGS_Mean(  iMin          ,&
                            iMax          ,&
                            bMin          ,&
                            bMax          ,&
                            Q_n           ,&
                            Q_np1_k       ,&
                            Q_np1         ,&
                            Q_np1_kp1     ,&
                            RHS           ,&
                            time          ,&
                            dQdTau        ,&
                            delta_x       ,&
                            delta_tau     ,&
                            delta_t       ,&
                            DS            ,&
                            nD            ,&
                            Dis           ,&
                            L_LowerDiag   ,&
                            L_Diag        ,&
                            U_UpperDiag   ,&
                            U_Diag        ,&
                            D_Diag        ,&
                            A_Plus        ,&
                            A_Minus       ,&
                            S_Plus        ,&
                            S_Minus       ,&
                            S_Jac         ,&
                            Delta_Q_star2, &
                            AA            ,&
                            BB            ,&
                            EE            ,&
                            Delta_Q_star  ,&
                            Delta_Q       ,&
                            bet           ,&
                            nTau          ,&
                            nStages       ,&
                            nI            ,&
                            Newtonl2Res   ,&
                            Source_Fac    ,&
                            Mach          ,&
                            u             ,&
                            p             ,&
                            c             ,&
                            gm1           ,&
                            rho           ,&
                            rho_Static    ,&
                            P_Static      ,&
                            Q_Exact       ,&
                            Jac_Curv      ,&
                            Inv_Jac_Curv  ,&
                            dZidX         ,&
                            dZidt         ,&
                            Area          ,&    
                            nDis) 

    USE GetGlobalParameter

    INCLUDE 'Inc_GetLUSGS_Mean.f90'
    
    END SUBROUTINE LUSGS_Mean

    SUBROUTINE LUSGS_Mean_Pet(  iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                Q_n             ,&
                                Q_np1_k         ,&
                                Q_np1           ,&
                                Q_nm1           ,&
                                Q_nm2           ,&
                                Q_np1_kp1       ,&
                                RHS             ,&
                                time            ,&
                                dQdTau          ,&
                                delta_x         ,&
                                delta_tau       ,&
                                delta_t         ,&
                                DS              ,&
                                nD              ,&
                                Dis             ,&
                                L_LowerDiag     ,&
                                L_Diag          ,&
                                U_UpperDiag     ,&
                                U_Diag          ,&
                                D_Diag          ,&
                                A_Plus          ,&
                                A_Minus         ,&
                                S_Plus          ,&
                                S_Minus         ,&
                                S_Jac           ,&
                                Delta_Q_star2   ,&
                                AA              ,&
                                BB              ,&
                                EE              ,&
                                Delta_Q_star    ,&
                                Delta_Q         ,&
                                bet             ,&
                                nTau            ,&
                                nStages         ,&
                                nI              ,&
                                Newtonl2Res     ,&
                                Source_Fac      ,&
                                Mach            ,&
                                u               ,&
                                p               ,&
                                c               ,&
                                gm1             ,&
                                rho             ,&
                                rho_Static      ,&
                                P_Static        ,&
                                Q_Exact         ,&
                                Jac_Curv        ,&
                                Inv_Jac_Curv    ,&
                                dZidX           ,&
                                dZidt           ,&
                                Area            ,&    
                                A1_Inflow       ,&
                                A2_Inflow       ,&
                                A3_Inflow       ,&
                                A1_Outflow      ,&
                                A2_Outflow      ,&
                                A3_Outflow      ,&
                                nDis            ,&
                                nBound) 

    USE GetGlobalParameter

    INCLUDE 'Inc_GetLUSGS_Mean_Pet.f90'
    
    END SUBROUTINE LUSGS_Mean_Pet

END MODULE GetLUSGS
