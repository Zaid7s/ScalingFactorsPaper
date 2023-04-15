MODULE GetLDU_Vectors
!!! Working Test Case
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: LDU_Vectors

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE LDU_Vectors( iMin            ,&
                            iMax            ,&
                            bMin            ,&
                            bMax            ,&
                            delta_x         ,&
                            delta_Zi        ,&
                            delta_t         ,&
                            A_Plus          ,&
                            A_Minus         ,&
                            S_Jac           ,&
                            Source_Fac      ,&
                            L_LowerDiag     ,&
                            U_UpperDiag     ,&
                            Tau             ,&
                            Jac_Curv        ,&
                            Inv_Jac_Curv    ,&
                            dZidX           ,&
                            dZidt)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin

    REAL(KIND = rDef), INTENT(IN) ::    delta_x     ,&
                                        delta_Zi    ,&
                                        delta_t

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt
    
    REAL(KIND = rDef), INTENT(IN), DIMENSION(:, :, :) ::    A_Plus      ,&
                                                            A_Minus     ,&
                                                            S_Jac

    REAL(KIND = rDef), INTENT(INOUT), DIMENSION(:, :, :) :: L_LowerDiag ,&
                                                            Tau         ,&
                                                            U_UpperDiag

    INTEGER:: i, j

    REAL(KIND = rDef):: alpha   ,&
                        a1

    alpha   = 1.0_rDef/1000000000.0_rDef
    a1      = 1.0_rDef/delta_Zi

    i            = iMin
    INCLUDE 'Tau.f90'
    L_LowerDiag(1, 1, i) = 0.0_rDef
    L_LowerDiag(1, 2, i) = 0.0_rDef
    L_LowerDiag(1, 3, i) = 0.0_rDef
    L_LowerDiag(2, 1, i) = 0.0_rDef
    L_LowerDiag(2, 2, i) = 0.0_rDef
    L_LowerDiag(2, 3, i) = 0.0_rDef
    L_LowerDiag(3, 1, i) = 0.0_rDef
    L_LowerDiag(3, 2, i) = 0.0_rDef
    L_LowerDiag(3, 3, i) = 0.0_rDef
    U_UpperDiag(1, 1, i) = a1*(A_Minus(1, 1, i + 1))
    U_UpperDiag(1, 2, i) = a1*(A_Minus(1, 2, i + 1))
    U_UpperDiag(1, 3, i) = a1*(A_Minus(1, 3, i + 1))
    U_UpperDiag(2, 1, i) = a1*(A_Minus(2, 1, i + 1))
    U_UpperDiag(2, 2, i) = a1*(A_Minus(2, 2, i + 1))
    U_UpperDiag(2, 3, i) = a1*(A_Minus(2, 3, i + 1))
    U_UpperDiag(3, 1, i) = a1*(A_Minus(3, 1, i + 1))
    U_UpperDiag(3, 2, i) = a1*(A_Minus(3, 2, i + 1))
    U_UpperDiag(3, 3, i) = a1*(A_Minus(3, 3, i + 1))

    DO i = iMin + 1, iMax - 1
        INCLUDE 'Tau.f90'
        L_LowerDiag(1, 1, i) = -a1*(A_Plus(1, 1, i - 1))
        L_LowerDiag(1, 2, i) = -a1*(A_Plus(1, 2, i - 1))
        L_LowerDiag(1, 3, i) = -a1*(A_Plus(1, 3, i - 1))
        L_LowerDiag(2, 1, i) = -a1*(A_Plus(2, 1, i - 1))
        L_LowerDiag(2, 2, i) = -a1*(A_Plus(2, 2, i - 1))
        L_LowerDiag(2, 3, i) = -a1*(A_Plus(2, 3, i - 1))
        L_LowerDiag(3, 1, i) = -a1*(A_Plus(3, 1, i - 1))
        L_LowerDiag(3, 2, i) = -a1*(A_Plus(3, 2, i - 1))
        L_LowerDiag(3, 3, i) = -a1*(A_Plus(3, 3, i - 1))

        U_UpperDiag(1, 1, i) = a1*(A_Minus(1, 1, i + 1))
        U_UpperDiag(1, 2, i) = a1*(A_Minus(1, 2, i + 1))
        U_UpperDiag(1, 3, i) = a1*(A_Minus(1, 3, i + 1))
        U_UpperDiag(2, 1, i) = a1*(A_Minus(2, 1, i + 1))
        U_UpperDiag(2, 2, i) = a1*(A_Minus(2, 2, i + 1))
        U_UpperDiag(2, 3, i) = a1*(A_Minus(2, 3, i + 1))
        U_UpperDiag(3, 1, i) = a1*(A_Minus(3, 1, i + 1))
        U_UpperDiag(3, 2, i) = a1*(A_Minus(3, 2, i + 1))
        U_UpperDiag(3, 3, i) = a1*(A_Minus(3, 3, i + 1))
    END DO    

    i            = iMax
    INCLUDE 'Tau.f90'
    U_UpperDiag(1, 1, i) = 0.0_rDef
    U_UpperDiag(1, 2, i) = 0.0_rDef
    U_UpperDiag(1, 3, i) = 0.0_rDef
    U_UpperDiag(2, 1, i) = 0.0_rDef
    U_UpperDiag(2, 2, i) = 0.0_rDef
    U_UpperDiag(2, 3, i) = 0.0_rDef
    U_UpperDiag(3, 1, i) = 0.0_rDef
    U_UpperDiag(3, 2, i) = 0.0_rDef
    U_UpperDiag(3, 3, i) = 0.0_rDef
    L_LowerDiag(1, 1, i) = -a1*(A_Plus(1, 1, i - 1))
    L_LowerDiag(1, 2, i) = -a1*(A_Plus(1, 2, i - 1))
    L_LowerDiag(1, 3, i) = -a1*(A_Plus(1, 3, i - 1))
    L_LowerDiag(2, 1, i) = -a1*(A_Plus(2, 1, i - 1))
    L_LowerDiag(2, 2, i) = -a1*(A_Plus(2, 2, i - 1))
    L_LowerDiag(2, 3, i) = -a1*(A_Plus(2, 3, i - 1))
    L_LowerDiag(3, 1, i) = -a1*(A_Plus(3, 1, i - 1))
    L_LowerDiag(3, 2, i) = -a1*(A_Plus(3, 2, i - 1))
    L_LowerDiag(3, 3, i) = -a1*(A_Plus(3, 3, i - 1))

    END SUBROUTINE LDU_Vectors

END MODULE GetLDU_Vectors
