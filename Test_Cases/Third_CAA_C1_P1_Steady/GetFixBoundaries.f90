MODULE GetFixBoundaries

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetInverse

    IMPLICIT NONE
    PRIVATE
    PUBLIC:: FixBoundaries

    INTEGER, PARAMETER:: rDef = REAL64

    INTERFACE FixBoundaries
        MODULE PROCEDURE FixBoundaries_64Bit
    END INTERFACE FixBoundaries

    CONTAINS

    SUBROUTINE FixBoundaries_64Bit( iMin                ,& 
                                    iMax                ,& 
                                    bMin                ,& 
                                    bMax                ,&
                                    LowerDiag           ,&
                                    Diag                ,&
                                    UpperDiag           ,&
                                    RHS                 ,&
                                    LowerDiag_Bound)

    INTEGER, INTENT(IN)::   iMin         ,&
                            iMax         ,&
                            bMin         ,&
                            bMax
    
    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(INOUT)::  LowerDiag     ,&
                                                            Diag          ,&
                                                            UpperDiag

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) :: LowerDiag_Bound

    REAL(KIND = rDef), DIMENSION(bMax, bMax) :: Inv_LowerDiag       ,&
                                                LowerDiagFac        ,&
                                                Dummy1              ,&
                                                Dummy2

    REAL(KIND = rDef), DIMENSION(bMax, 1) ::    Dummy3  ,&
                                                Dummy4  ,&
                                                Dummy5

    REAL(KIND = rDef) :: DetInv

    LOGICAL :: OK_FLAG

    INTEGER :: i, j, k

    LowerDiagFac(1, 1) = LowerDiag(iMax - 1, 1, 1)
    LowerDiagFac(1, 2) = LowerDiag(iMax - 1, 1, 2)
    LowerDiagFac(1, 3) = LowerDiag(iMax - 1, 1, 3)
    LowerDiagFac(2, 1) = LowerDiag(iMax - 1, 2, 1)
    LowerDiagFac(2, 2) = LowerDiag(iMax - 1, 2, 2)
    LowerDiagFac(2, 3) = LowerDiag(iMax - 1, 2, 3)
    LowerDiagFac(3, 1) = LowerDiag(iMax - 1, 3, 1)
    LowerDiagFac(3, 2) = LowerDiag(iMax - 1, 3, 2)
    LowerDiagFac(3, 3) = LowerDiag(iMax - 1, 3, 3)

    CALL Inverse(   A       = LowerDiagFac  ,&
                    AINV    = Inv_LowerDiag ,&
                    OK_FLAG = OK_FLAG       ,&
                    bMax    = bMax)

    Dummy1(1, 1)  = Diag(iMax - 1, 1, 1)
    Dummy1(1, 2)  = Diag(iMax - 1, 1, 2)
    Dummy1(1, 3)  = Diag(iMax - 1, 1, 3)
    Dummy1(2, 1)  = Diag(iMax - 1, 2, 1)
    Dummy1(2, 2)  = Diag(iMax - 1, 2, 2)
    Dummy1(2, 3)  = Diag(iMax - 1, 2, 3)
    Dummy1(3, 1)  = Diag(iMax - 1, 3, 1)
    Dummy1(3, 2)  = Diag(iMax - 1, 3, 2)
    Dummy1(3, 3)  = Diag(iMax - 1, 3, 3)

    Dummy2(1, 1)  = UpperDiag(iMax - 1, 1, 1)
    Dummy2(1, 2)  = UpperDiag(iMax - 1, 1, 2)
    Dummy2(1, 3)  = UpperDiag(iMax - 1, 1, 3)
    Dummy2(2, 1)  = UpperDiag(iMax - 1, 2, 1)
    Dummy2(2, 2)  = UpperDiag(iMax - 1, 2, 2)
    Dummy2(2, 3)  = UpperDiag(iMax - 1, 2, 3)
    Dummy2(3, 1)  = UpperDiag(iMax - 1, 3, 1)
    Dummy2(3, 2)  = UpperDiag(iMax - 1, 3, 2)
    Dummy2(3, 3)  = UpperDiag(iMax - 1, 3, 3)

    Dummy3(1, 1)  = RHS(iMax - 1, 1)
    Dummy3(2, 1)  = RHS(iMax - 1, 2)
    Dummy3(3, 1)  = RHS(iMax - 1, 3)

    Dummy4(1, 1)  = RHS(iMax, 1)
    Dummy4(2, 1)  = RHS(iMax, 2)
    Dummy4(3, 1)  = RHS(iMax, 3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    LowerDiag(iMax, :, :) = LowerDiag(iMax, :, :) - &
                            MATMUL(MATMUL(LowerDiag_Bound,Inv_LowerDiag),Dummy1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Diag(iMax, :, :) = Diag(iMax, :, :) -           &
                        MATMUL(MATMUL(LowerDiag_Bound,Inv_LowerDiag),Dummy2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Dummy5 = Dummy4 - MATMUL((MATMUL(LowerDiag_Bound,Inv_LowerDiag)), Dummy3)

    RHS(iMax, 1) = Dummy5(1, 1)
    RHS(iMax, 2) = Dummy5(2, 1)
    RHS(iMax, 3) = Dummy5(3, 1)

    END SUBROUTINE FixBoundaries_64Bit

END MODULE GetFixBoundaries      
