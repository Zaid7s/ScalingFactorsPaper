MODULE GetDX

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: DX

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         
    SUBROUTINE DX(  xMin            ,&
                    xMax            ,&
                    iMin            ,&
                    iMax            ,&
                    delta_xNew      ,& 
                    xNew            ,&
                    x0              ,&
                    dX_Left         ,&
                    dX_Mid          ,&
                    dX_Right        ,&
                    nG)

    INTEGER, INTENT(IN) ::  iMin        ,&
                            iMax        ,&
                            x0          ,&
                            xMin        ,&
                            xMax        ,&
                            nG

    REAL(KIND = rDef), INTENT(IN) ::    dX_Left     ,&
                                        dX_Mid      ,&
                                        dX_Right

    REAL(KIND = rDef), DIMENSION(:), INTENT(INOUT) ::   delta_xNew  ,&
                                                        xNew
    
    INTEGER ::  i
    REAL(KIND = rDef) ::    dXX         ,&
                            MMin        ,&
                            nFac        ,&
                            Fac         ,&
                            xxx
    IF (nG == 1) THEN
        ! 401 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.010867959652098_rDef
        Fac     = 0.048922224510661_rDef
        dXX     = 0.108679596520982_rDef
    ELSEIF (nG == 2) THEN
        ! 351 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.011920491998808_rDef
        Fac     = 0.053665022391179_rDef
        dXX     = 0.119204919988084_rDef
    ELSEIF (nG == 3) THEN
        ! 301 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.013199691209066_rDef
        Fac     = 0.059418536456994_rDef
        dXX     = 0.131996912090662_rDef
    ELSEIF (nG == 4) THEN
        ! 251 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.014785850480053_rDef
        Fac     = 0.066558647613910_rDef
        dXX     = 0.147858504800525_rDef
    ELSEIF (nG == 5) THEN
        ! 211 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.016358437753843_rDef
        Fac     = 0.073637664295401_rDef
        dXX     = 0.163584377538427_rDef
    ELSEIF (nG == 6) THEN
        ! 181 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.017775404594133_rDef
        Fac     = 0.080040209759216_rDef
        dXX     = 0.177754045941335_rDef
    ELSEIF (nG == 7) THEN
        ! 151 Grid Points
        xxx     = 8.0_rDef
        MMin    = -0.999999999924497_rDef
        nFac    = 0.018865457118764_rDef
        Fac     = 0.084957938109137_rDef
        dXX     = 0.188654571187638_rDef
    END IF

    IF (nG == 8) THEN
        delta_xNew = delta_xNew
    ELSEIF (nG == 9) THEN
        OPEN(15, FILE = 'Area_Matlab.txt', FORM = 'FORMATTED', STATUS = 'OLD')
            DO i = iMin, iMax
                READ(15, *) xNew(i)                       ,&
                            delta_xNew(i)
            END DO
        CLOSE(15)
    ELSE
        DO i = iMin, iMax
            IF (xNew(i) < REAL(-x0, rDef)) THEN
                delta_xNew(i) = dXX
            ELSEIF ((xNew(i) >= REAL(-x0, rDef)).AND.&
                                    (xNew(i) < 0.0_rDef))                   THEN
                delta_xNew(i)   = ((TANH(-xxx*(xNew(i) + &
                                    (REAL(x0, rDef) - 0.5_rDef)))-MMin))&
                                    *Fac+nFac
            ELSEIF ((xNew(i) >= 0.0_rDef).AND.(xNew(i) <= REAL(x0, rDef)))  THEN
                delta_xNew(i)   = ((TANH(-xxx*(-xNew(i) + &
                                    (REAL(x0, rDef) - 0.5_rDef)))-MMin))&
                                    *Fac+nFac
            ELSEIF ((xNew(i) > REAL(x0, rDef)).AND.&
                                (xNew(i) <= REAL(xMax, rDef)))              THEN
                delta_xNew(i) = dXX
            ELSE
                delta_xNew(i) = 0.0_rDef
            END IF
        END DO

        DO i = iMin, iMax
            xNew(i) = 0.0_rDef
        END DO

        DO i = iMin, iMax
            IF (i == iMin) THEN
                xNew(i) = REAL(xMin, rDef)
            ELSE
                xNew(i) = xNew(i-1) + delta_xNew(i)
            END IF
        END DO
    END IF

    END SUBROUTINE DX

END MODULE GetDX
