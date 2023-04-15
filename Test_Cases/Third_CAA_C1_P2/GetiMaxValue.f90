MODULE GetiMaxValue

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: iMaxValue


    CONTAINS
         
    SUBROUTINE iMaxValue(   x0          ,&
                            iMin        ,&
                            iMax        ,&
                            dX_Left     ,&
                            dX_Mid      ,&
                            dX_Right    ,&
                            xNew        ,&
                            deltaX_New  ,&
                            nG)
    USE GetGlobalParameter

    INTEGER, INTENT(IN) ::  nG

    INTEGER, INTENT(INOUT)::    iMax    ,&
                                iMin    ,&
                                x0

    REAL(KIND = rDef), INTENT(INOUT) ::   dX_Left     ,&
                                        dX_Mid      ,&
                                        dX_Right

    REAL(KIND = rDef), INTENT(INOUT), DIMENSION(:), ALLOCATABLE ::    xNew    ,&
                                                                    deltaX_New

    INTEGER ::  iMin_Left   ,& 
                iMax_Left   ,&
                iMax_Left1  ,&
                iMin_Mid    ,&
                iMax_Mid    ,&
                iMax_Mid1   ,&
                iMin_Right  ,&
                iMax_Right  ,&
                iMax_Right1, &
                i

    CHARACTER(LEN = 255) :: Current_Directory 
    CALL GET_ENVIRONMENT_VARIABLE('PWD',Current_Directory)

    x0          = 2
    
    iMax_Left1  = 101
    iMax_Right1 = 101

    IF      (nG == 1) THEN
        iMax_Mid1   = 1601     ! 401 Grid Points
    ELSEIF  (nG == 2) THEN
        iMax_Mid1   = 1351     ! 351 Grid Points
    ELSEIF  (nG == 3) THEN
        iMax_Mid1   = 1101     ! 301 Grid Points
    ELSEIF  (nG == 4) THEN
        iMax_Mid1   = 851      ! 251 Grid Points
    ELSEIF  (nG == 5) THEN
        iMax_Mid1   = 651      ! 211 Grid Points
    ELSEIF  (nG == 6) THEN
        iMax_Mid1   = 501      ! 181 Grid Points
    ELSEIF  (nG == 7) THEN
        iMax_Mid1   = 401      ! 151 Grid Points
    END IF

    dX_Left     = (REAL(xMax-xMin, rDef))/(REAL(iMax_Left1-1, rDef))
    dX_Mid      = (REAL(xMax-xMin, rDef))/(REAL(iMax_Mid1-1, rDef))
    dX_Right    = (REAL(xMax-xMin, rDef))/(REAL(iMax_Right1-1, rDef))

    iMin_Left   = 1
    iMax_Left   = (xMax-x0)/(dX_Right) + 1
    
    iMin_Mid    = iMax_Left
    iMax_Mid    = iMin_Mid+NINT((x0+x0)/(dX_Mid))

    iMin_Right  = iMax_Mid
    iMax_Right  = iMax_Mid + (-x0-xMin)/(dX_Left)


    IF (nG == 8) THEN
        iMax = 6401
    ELSE
        iMax        = iMax_Right
    END IF

    IF (nG == 9) THEN
        OPEN(14, FILE = 'Area_Matlab_Grid_Point.txt', FORM = 'FORMATTED')
            READ(14, *) iMax
        CLOSE(14)
    END IF
    
    ALLOCATE(   deltaX_New(iMax)  ,&
                xNew(iMax))

    DO i = iMin, iMax
        xNew(i)         = 0.0_rDef
        deltaX_New(i)   = 0.0_rDef
    END DO

    DO i = iMin, iMax
        IF ((i >= iMin_Left).AND.(i <= iMax_Left)) THEN
            deltaX_New(i) = dX_Right 
        ELSEIF ((i > iMin_Mid).AND.(i <= iMax_Mid)) THEN 
            deltaX_New(i) = dX_Mid
        ELSEIF ((i > iMin_Right).AND.(i <= iMax_Right)) THEN 
            deltaX_New(i) = dX_Left
        END IF
    END DO

    IF (nG == 8) THEN
        deltaX_New = REAL(xMax-xMin, rDef)/REAL(iMax-1, rDef)
    ELSEIF (nG == 9) THEN
        OPEN(15, FILE = 'Area_Matlab.txt', FORM = 'FORMATTED', STATUS = 'OLD')
            DO i = iMin, iMax
                READ(15, *) xNew(i)                       ,&
                            deltaX_New(i)
            END DO
        CLOSE(15)
    ELSE
        deltaX_New = deltaX_New
    END IF

    DO i = iMin, iMax
        IF (i == iMin) THEN
            xNew(i) = REAL(xMin, rDef)
        ELSE
            xNew(i) = xNew(i-1) + deltaX_New(i)
        END IF
    END DO

    END SUBROUTINE iMaxValue

END MODULE GetiMaxValue
