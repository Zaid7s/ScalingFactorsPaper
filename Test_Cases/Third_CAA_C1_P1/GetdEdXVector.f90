MODULE GetdEdXVector

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetSpatialDerivative
    USE GetInflowBC
    USE GetOutflowBC
    USE GetCompactScheme
    
    IMPLICIT NONE
    PRIVATE

    PUBLIC:: dEdXVector

    INTERFACE dEdXVector
        MODULE PROCEDURE dEdXVector_Mean
        MODULE PROCEDURE dEdXVector_Mean_Pet
    END INTERFACE dEdXVector

    CONTAINS
         
    SUBROUTINE dEdXVector_Mean( iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                E               ,&
                                time            ,&
                                dEdX            ,&
                                delta_X         ,&
                                delta_tau       ,&
                                delta_t         ,&
                                DS              ,&
                                Dis             ,&
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
                                dZidt)

    USE GetGlobalParameter

    INCLUDE 'Inc_GetdEdXVector_Mean.f90'

    END SUBROUTINE dEdXVector_Mean
         
    SUBROUTINE dEdXVector_Mean_Pet( iMin            ,&
                                    iMax            ,&
                                    bMin            ,&
                                    bMax            ,&
                                    E               ,&
                                    time            ,&
                                    dEdX            ,&
                                    delta_X         ,&
                                    delta_tau       ,&
                                    delta_t         ,&
                                    DS              ,&
                                    Dis             ,&
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
                                    dZidt           ,&
                                    A1_Inflow       ,&
                                    A2_Inflow       ,&
                                    A3_Inflow       ,&
                                    A1_Outflow      ,&
                                    A2_Outflow      ,&
                                    A3_Outflow      ,&
                                    nBound)

    USE GetGlobalParameter

    INCLUDE 'Inc_GetdEdXVector_Mean_Pet.f90'

    END SUBROUTINE dEdXVector_Mean_Pet

END MODULE GetdEdXVector
