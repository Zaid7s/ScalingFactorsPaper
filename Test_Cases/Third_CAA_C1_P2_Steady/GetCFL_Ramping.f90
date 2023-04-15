!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Switched Evolution Relaxation (SER) Method By Mulder and van Leer           ! 
! SER is a CFL evolution scheme that ties CFL                                 !
! ramping to the nonlinear residual                                           !
!           CFL^{n+1}   = epsi*CFL^{n}                                        !
!           epsi        = Residual^{n-1}/Residual^{n}                         !
! There are also CFL max and min limits                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu            ----> Is the CFL which is used both as an input and output    !
!                       (CFL is modified after it leaves)                     !
! Residual     ----> Referes to the Residual                                  !
! Fac2         ----> Stores the Reisudal at the previous time level
! Fac1         ----> Gets the residual of the cirrent time level
!
! ***   Note, given that the CFL Ramping algorithm needs theResidual of the ***
! ***   previous time level, the algorithm isn't self starting. Hence an    ***
! ***   IF statement is needed such that R^{n-1} at n = 0 is equal to     ***
! ***   R^{n}. Otherwise code should perform nornally over                  ***
! ***   CFL_Ramping_Part1 and CFL_Ramping_Part2                             ***
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE GetCFL_Ramping
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    IMPLICIT NONE
    PRIVATE
    PUBLIC:: CFLRampObject,  CreateObject, CFL_Ramping

    INTEGER, PARAMETER:: rDef = REAL64

    INTERFACE CreateObject
        MODULE PROCEDURE CreateRampingObject
    END INTERFACE CreateObject

    INTERFACE CFL_Ramping
        MODULE PROCEDURE CFL_Ramp_Part1
        MODULE PROCEDURE CFL_Ramp_Part2
    END INTERFACE CFL_Ramping

    TYPE CFLRampObject
        PRIVATE

        REAL(KIND = rDef) ::    Residual    ,&
                                nu          ,&
                                Fac2

END TYPE CFLRampObject

    !!! Define Local Variables !!!

    REAL(KIND = rDef):: epsi    ,&
                        Fac1
CONTAINS

    SUBROUTINE CreateRampingObject( object      ,&
                                    Residual    ,&
                                    Fac2)

        TYPE(CFLRampObject), INTENT(INOUT):: object

        REAL(KIND = rDef), INTENT(IN)    :: Residual   ,&
                                            Fac2

        object%Residual = Residual
        object%Fac2     = Fac2

    END SUBROUTINE CreateRampingObject

    SUBROUTINE CFL_Ramp_Part1(  object, &
                                Factor2)

        TYPE(CFLRampObject), INTENT(INOUT):: object

        REAL(KIND = rDef), INTENT(OUT) ::  Factor2

        object%Fac2 = Factor2
        Fac1        = object%Residual

        epsi        = Fac1/object%Residual
        Factor2     = Fac1

    END SUBROUTINE CFL_Ramp_Part1

    SUBROUTINE CFL_Ramp_Part2(  object  ,&
                                CFL     ,&
                                Factor2)

        TYPE(CFLRampObject), INTENT(INOUT):: object

        REAL(KIND = rDef), INTENT(INOUT):: CFL     ,&
                                            Factor2

        object%nu   = CFL
        object%Fac2 = Factor2
        epsi        = Factor2/object%Residual
        CFL         = epsi*object%nu
        Factor2     = object%Residual

        !!!! This sets the minimum CFL !!!!!
        IF (CFL < 1.0_rDef) THEN
            CFL = 1.0_rDef
        END IF
        !!!! This sets the maximum CFL !!!!!
        IF (CFL > 7.0_rDef) THEN
            CFL = 7.0_rDef
        END IF

    END SUBROUTINE CFL_Ramp_Part2

END MODULE GetCFL_Ramping
