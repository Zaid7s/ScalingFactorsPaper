MODULE GetGlobalParameter
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE

    INTEGER, PARAMETER ::   rDef    = REAL64    ,&
                            xMin    = -10       ,&
                            xMax    = 10

    REAL(KIND = rDef), PARAMETER :: gam         = 1.40_rDef                 ,&
                                    delta_Zi    = 1.0_rDef                  ,&
                                    tol         = 5E-12_rDef                ,&
                                    pi          = 4.0_rDef*ATAN(1.0_rDef)   ,&
                                    p_Exit      = 0.6071752_rDef            ,&
                                    p0In        = 0.734620030346377_rDef    ,&
                                    rho0In      = 1.020252612017331_rDef    ,&
                                    t0In        = 1.008052349360178_rDef 


    !INTEGER, PARAMETER :: ndQdT = 1
    !REAL(KIND = rDef), PARAMETER :: a1  = 1.0_rDef

    !INTEGER, PARAMETER :: ndQdT = 2
    !REAL(KIND = rDef), PARAMETER :: a1  = 2.0_rDef/3.0_rDef

    INTEGER, PARAMETER :: ndQdT = 3
    REAL(KIND = rDef), PARAMETER :: a1  = 3.0_rDef/5.0_rDef

END MODULE GetGlobalParameter
