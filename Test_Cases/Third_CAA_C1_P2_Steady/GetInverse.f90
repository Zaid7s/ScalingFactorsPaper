MODULE GetInverse
    
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Inverse
    
    INTEGER, PARAMETER:: rDef = REAL64

    INTERFACE Inverse
        MODULE PROCEDURE Inverse3X3
    END INTERFACE Inverse

    CONTAINS

    SUBROUTINE Inverse3X3(  A           ,&
                            AINV        ,&
                            OK_FLAG     ,&
                            bMax)
    
    REAL(KIND = rDef), PARAMETER :: EPS = 5E-13_rDef
    INTEGER, INTENT(IN) :: bMax
    REAL(KIND = rDef), DIMENSION(bMax, bMax), INTENT(IN)  :: A
    REAL(KIND = rDef), DIMENSION(bMax, bMax), INTENT(OUT) :: AINV
    
    LOGICAL, INTENT(OUT) :: OK_FLAG
    
    REAL(KIND = rDef) :: DET
    
    REAL(KIND = rDef), DIMENSION(bMax, bMax) :: COFACTOR
    
    
    DET =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)
    
    IF (ABS(DET) .LE. EPS) THEN
        AINV = 0.0_rDef
        OK_FLAG = .FALSE.
        RETURN
    END IF
    
    COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    
    AINV = TRANSPOSE(COFACTOR) / DET
    
    OK_FLAG = .TRUE.
    
    END SUBROUTINE Inverse3X3

END MODULE GetInverse
