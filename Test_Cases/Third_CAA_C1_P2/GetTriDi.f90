MODULE GetTriDi

    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetRHSVector
    USE GetRungeKutta
    USE GetDiag
    USE BlockTriDiSolver1D
    USE GetFixBoundaries

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: TriDi

    INTERFACE TriDi
        MODULE PROCEDURE TriDi_Mean
        MODULE PROCEDURE TriDi_Mean_Pet
    END INTERFACE TriDi

    CONTAINS

    SUBROUTINE TriDi_Mean(  iMin            ,&
                            iMax            ,&
                            bMin            ,&
                            bMax            ,&
                            Delta_Q         ,&
                            A_Jac           ,&
                            S_Jac           ,&
                            Q_n             ,&
                            Q_np1_k         ,&
                            Q_np1_kp1       ,&
                            Q_np1           ,&
                            Q_nm1           ,&
                            Q_nm2           ,&
                            RHS             ,&
                            time            ,&
                            delta_x         ,&
                            delta_tau       ,&
                            delta_t         ,&
                            DS              ,&
                            nD              ,&
                            Dis             ,&
                            bet             ,&
                            nTau            ,&
                            nStages         ,&
                            Scaling_Fac     ,&
                            E               ,&
                            E_GPT           ,&
                            Newtonl2Res     ,&
                            dQdTau          ,&
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
                            nDis            ,&
                            Area) 

    USE GetGlobalParameter

    INCLUDE 'Inc_GetTriDi_Mean.f90'

    END SUBROUTINE TriDi_Mean

    SUBROUTINE TriDi_Mean_Pet(  iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                Delta_Q         ,&
                                A_Jac           ,&
                                S_Jac           ,&
                                Q_n             ,&
                                Q_np1_k         ,&
                                Q_np1_kp1       ,&
                                Q_np1           ,&
                                Q_nm1           ,&
                                Q_nm2           ,&
                                RHS             ,&
                                time            ,&
                                delta_x         ,&
                                delta_tau       ,&
                                delta_t         ,&
                                DS              ,&
                                nD              ,&
                                Dis             ,&
                                bet             ,&
                                nTau            ,&
                                nStages         ,&
                                Scaling_Fac     ,&
                                E               ,&
                                E_GPT           ,&
                                Newtonl2Res     ,&
                                dQdTau          ,&
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
                                nDis            ,&
                                Area            ,&    
                                A1_Inflow       ,&
                                A2_Inflow       ,&
                                A3_Inflow       ,&
                                A1_Outflow      ,&
                                A2_Outflow      ,&
                                A3_Outflow      ,&
                                nBound) 

    USE GetGlobalParameter

    INCLUDE 'Inc_GetTriDi_Mean_Pet.f90'

    END SUBROUTINE TriDi_Mean_Pet
!!!!!!!!!!!!!!!!!!!!
END MODULE GetTriDi
