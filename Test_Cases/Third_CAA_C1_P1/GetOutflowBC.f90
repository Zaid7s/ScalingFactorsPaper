MODULE GetOutflowBC
!! New Boundaries

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE BlockTriDiSolver1D
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: OutflowBC

    INTERFACE OutflowBC
        MODULE PROCEDURE OutflowBC_Mean
        MODULE PROCEDURE OutflowBC_Mean_Pet
    END INTERFACE OutflowBC

    CONTAINS
         

    SUBROUTINE OutflowBC_Mean(  iMin            ,&
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

    INCLUDE 'Inc_GetOutflowBC_Mean.f90'

    END SUBROUTINE OutflowBC_Mean

    SUBROUTINE OutflowBC_Mean_Pet(  iMin            ,&
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

    INCLUDE 'Inc_GetOutflowBC_Mean_Pet.f90'

    END SUBROUTINE OutflowBC_Mean_Pet

END MODULE GetOutflowBC
