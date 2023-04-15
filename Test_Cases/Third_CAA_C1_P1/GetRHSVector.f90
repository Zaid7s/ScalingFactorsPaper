MODULE GetRHSVector
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetdEdXVector
    USE GetDissipation
    USE GetdQdT_BDF

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: RHSVector

    INTERFACE RHSVector
        MODULE PROCEDURE RHSVector_Mean
        MODULE PROCEDURE RHSVector_Mean_Pet
    END INTERFACE RHSVector

    CONTAINS

    SUBROUTINE RHSVector_Mean(  iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                Q_n             ,&
                                Q_np1_k         ,&
                                dQdTau          ,&
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
                                nDis            ,&
                                A3_Outflow)
    USE GetGlobalParameter

    INCLUDE 'Inc_GetRHSVector_Mean.f90'

    END SUBROUTINE RHSVector_Mean

    SUBROUTINE RHSVector_Mean_Pet(  iMin            ,&
                                    iMax            ,&
                                    bMin            ,&
                                    bMax            ,&
                                    Q_n             ,&
                                    Q_np1_k         ,&
                                    Q_nm1           ,&
                                    Q_nm2           ,&
                                    dQdTau          ,&
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
                                    nDis            ,&
                                    A3_Outflow      ,&
                                    nBound)
    USE GetGlobalParameter

    INCLUDE 'Inc_GetRHSVector_Mean_Pet.f90'

    END SUBROUTINE RHSVector_Mean_Pet

END MODULE GetRHSVector
