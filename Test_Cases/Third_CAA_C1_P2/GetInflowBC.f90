MODULE GetInflowBC
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetDeltaA
    USE BlockTriDiSolver1D
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: InflowBC

    INTERFACE InflowBC
        MODULE PROCEDURE InflowBC_Mean
        MODULE PROCEDURE InflowBC_Mean_Pet
    END INTERFACE InflowBC

    CONTAINS
         
    SUBROUTINE InflowBC_Mean(   iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                E               ,&
                                time            ,&
                                dEdX            ,&
                                delta_x         ,&
                                delta_t         ,&
                                delta_tau       ,&
                                DS              ,&
                                Mean_Mach       ,&
                                Mean_u          ,&
                                Mean_p          ,&
                                Mean_c          ,&
                                gm1             ,&
                                Mean_rho        ,&
                                rho_Static      ,&
                                P_Static        ,&
                                Q_np1_k         ,&
                                Q_Exact         ,&
                                Jac_Curv        ,&
                                Inv_Jac_Curv    ,&
                                dZidX)

    USE GetGlobalParameter

    INCLUDE 'Inc_GetInflowBC_Mean.f90'

    END SUBROUTINE InflowBC_Mean

    SUBROUTINE InflowBC_Mean_Pet(   iMin            ,&
                                    iMax            ,&
                                    bMin            ,&
                                    bMax            ,&
                                    E               ,&
                                    time            ,&
                                    dEdX            ,&
                                    delta_x         ,&
                                    delta_t         ,&
                                    delta_tau       ,&
                                    DS              ,&
                                    Mean_Mach       ,&
                                    Mean_u          ,&
                                    Mean_p          ,&
                                    Mean_c          ,&
                                    gm1             ,&
                                    Mean_rho        ,&
                                    rho_Static      ,&
                                    P_Static        ,&
                                    Q_np1_k         ,&
                                    Q_Exact         ,&
                                    Jac_Curv        ,&
                                    Inv_Jac_Curv    ,&
                                    dZidX           ,&
                                    A_S_BC          ,&
                                    A_Plus_BC       ,&
                                    A_Minus_BC      ,&
                                    nBound)

    USE GetGlobalParameter

    INCLUDE 'Inc_GetInflowBC_Mean_Pet.f90'

    END SUBROUTINE InflowBC_Mean_Pet

END MODULE GetInflowBC
