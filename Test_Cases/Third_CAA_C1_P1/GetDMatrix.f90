MODULE GetDMatrix
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: DMatrix

    CONTAINS
         

    SUBROUTINE DMatrix( iMin            ,&
                        iMax            ,&
                        bMin            ,&
                        bMax            ,&
                        delta_x         ,&
                        delta_t         ,&
                        A_Plus          ,&
                        A_Minus         ,&
                        S_Plus          ,&
                        S_Minus         ,&
                        S_Jac           ,&
                        Source_Fac      ,&
                        D_Diag          ,&
                        Jac_Curv        ,&
                        Inv_Jac_Curv    ,&
                        dZidX           ,&
                        dZidt)

    USE GetGlobalParameter

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin

    REAL(KIND = rDef), INTENT(IN) ::    delta_x    ,&
                                        delta_t

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt
    
    REAL(KIND = rDef), INTENT(IN), DIMENSION(:, :, :) ::    A_Plus      ,&
                                                            A_Minus     ,&
                                                            S_Jac       ,&
                                                            S_Plus      ,&
                                                            S_Minus

    REAL(KIND = rDef), INTENT(INOUT), DIMENSION(:, :, :) ::   D_Diag

    INTEGER:: i, j
    REAL(KIND = rDef):: alpha

    alpha   = 1.0_rDef/100000000.0_rDef


    DO i = iMin, iMax
        D_Diag(1, 1, i) = alpha+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(1, 1, i) - A_Minus(1, 1, i)) -     &
                            1.0_rDef*S_Jac(1, 1, i)*a1   
        D_Diag(1, 2, i) = 0.0_rDef+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(1, 2, i) - A_Minus(1, 2, i)) -     &
                            1.0_rDef*S_Jac(1, 2, i)*a1    
        D_Diag(1, 3, i) = 0.0_rDef+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(1, 3, i) - A_Minus(1, 3, i)) -     &
                            1.0_rDef*S_Jac(1, 3, i)*a1    
        D_Diag(2, 1, i) = 0.0_rDef+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(2, 1, i) - A_Minus(2, 1, i)) -     &
                            1.0_rDef*S_Jac(2, 1, i)*a1    
        D_Diag(2, 2, i) = alpha+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(2, 2, i) - A_Minus(2, 2, i)) -     &
                            1.0_rDef*S_Jac(2, 2, i)*a1    
        D_Diag(2, 3, i) = 0.0_rDef+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(2, 3, i) - A_Minus(2, 3, i)) -     &
                            1.0_rDef*S_Jac(2, 3, i)*a1 
        D_Diag(3, 1, i) = 0.0_rDef+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(3, 1, i) - A_Minus(3, 1, i)) -     &
                            1.0_rDef*S_Jac(3, 1, i)*a1 
        D_Diag(3, 2, i) = 0.0_rDef+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(3, 2, i) - A_Minus(3, 2, i)) -     &
                            1.0_rDef*S_Jac(3, 2, i)*a1 
        D_Diag(3, 3, i) = alpha+(1.0_rDef/delta_Zi)*a1*    &
                             (A_Plus(3, 3, i) - A_Minus(3, 3, i)) -     &
                            1.0_rDef*S_Jac(3, 3, i)*a1 
    END DO

    END SUBROUTINE DMatrix

END MODULE GetDMatrix
