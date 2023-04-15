MODULE GetAreaDerivative

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE BlockTriDiSolver1D
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: AreaDerivative

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE AreaDerivative(  iMin      ,&
                                iMax      ,&
                                A         ,&
                                dAdX      ,&
                                delta_X   ,&
                                DS)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           DS

    REAL(KIND = rDef), INTENT(IN):: delta_x
                                        
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN):: A

    REAL(KIND = rDef), DIMENSION(:), INTENT(INOUT):: dAdX

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   LowerDiag       ,&
                                                            UpperDiag       ,&
                                                            MainDiag

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  RHS      ,&
                                                        dAdX_TriDi


    INTEGER:: i

    INTEGER, DIMENSION(1):: iiMin, iiMax

    iiMin = iMin
    iiMax = iMax
    
! dAdX for Interior Domain

    IF      (DS == 1) THEN

        i = iMin
        dAdX(i) = ( -3.0_rDef*A(i)      &
                   + 4.0_rDef*A(i+1)  &
                   -          A(i+2)) &
                       /(2.0_rDef*delta_x)

        DO i = iMin+1, iMax-1 
            dAdX(i)   = (1.0_rDef/(delta_x*2.0_rDef))&
                        *(A(i+1) - A(i-1))
        END DO

        i = iMax
        dAdX(i) = -(-3.0_rDef*A(i)       &
                   + 4.0_rDef*A(i-1)   &
                   -          A(i-2))  &
                       /(2.0_rDef*delta_x)
        
    ELSE IF (DS == 2) THEN
        i       = iMin
        dAdX(i) = (-25.0_rDef*A(i-0)   &
                   +48.0_rDef*A(i+1)   &
                   -36.0_rDef*A(i+2)   &
                   +16.0_rDef*A(i+3)   &
                    -3.0_rDef*A(i+4))  &
                    /(12.0_rDef*delta_x)
        
        i       = iMin+1
        dAdX(i) = (-3.0_rDef*A(i-1)   &
                  -10.0_rDef*A(i+0)   &
                  +18.0_rDef*A(i+1)   &
                   -6.0_rDef*A(i+2)   &
                   +1.0_rDef*A(i+3))  &
                    /(12.0_rDef*delta_x)
               
        DO i = iMin+2, iMax-2 
            dAdX(i)   = (1.0_rDef/(delta_x*12.0_rDef))&
                                *(        A(i-2)&
                                -8.0_rDef*A(i-1)&
                                +8.0_rDef*A(i+1)&
                                -1.0_rDef*A(i+2))
        END DO
        
        i       = iMax-1
        dAdX(i) = -(-3.0_rDef*A(i+1)   &
                   -10.0_rDef*A(i-0)   &
                   +18.0_rDef*A(i-1)   &
                    -6.0_rDef*A(i-2)   &
                    +1.0_rDef*A(i-3))  &
                     /(12.0_rDef*delta_x)

        i       = iMax
        dAdX(i) = -(-25.0_rDef*A(i+0)   &
                    +48.0_rDef*A(i-1)   &
                    -36.0_rDef*A(i-2)   &
                    +16.0_rDef*A(i-3)   &
                     -3.0_rDef*A(i-4))  &
                     /(12.0_rDef*delta_x)

    ELSE IF (DS == 3) THEN
        
        i       = iMin
        dAdX(i) = (     -147.0_rDef*A(i+0)    &
                        +360.0_rDef*A(i+1)    &
                        -450.0_rDef*A(i+2)    &
                        +400.0_rDef*A(i+3)    &
                        -225.0_rDef*A(i+4)    &
                         +72.0_rDef*A(i+5)    &
                         -10.0_rDef*A(i+6))   &
                            /(60.0_rDef*delta_x)

        i       = iMin+1
        dAdX(i) = (     -10.0_rDef*A(i-1)     &
                        -77.0_rDef*A(i+0)     &
                       +150.0_rDef*A(i+1)     &
                       -100.0_rDef*A(i+2)     &
                        +50.0_rDef*A(i+3)     &
                        -15.0_rDef*A(i+4)     &
                         +2.0_rDef*A(i+5))    &
                            /(60.0_rDef*delta_x)
        
        i       = iMin+2
        dAdX(i) = (      2.0_rDef*A(i-2)      &
                       -24.0_rDef*A(i-1)      &
                       -35.0_rDef*A(i+0)      &
                       +80.0_rDef*A(i+1)      &    
                       -30.0_rDef*A(i+2)      &
                        +8.0_rDef*A(i+3)      &
                        -1.0_rDef*A(i+4))     &
                            /(60.0_rDef*delta_x)
        
        DO i = iMin+3, iMax-3
            dAdX(i)   = (1.0_rDef/(delta_x*60.0_rDef))&
                                *(   -1.0_rDef*A(i-3)&
                                     +9.0_rDef*A(i-2)&
                                    -45.0_rDef*A(i-1)&
                                    +45.0_rDef*A(i+1)&
                                     -9.0_rDef*A(i+2)&
                                     +1.0_rDef*A(i+3))
        END DO
        
        i       = iMax-2
        dAdX(i) = -(  2.0_rDef*A(i+2)      &
                    -24.0_rDef*A(i+1)      &
                    -35.0_rDef*A(i+0)      &
                    +80.0_rDef*A(i-1)      &    
                    -30.0_rDef*A(i-2)      &
                     +8.0_rDef*A(i-3)      &
                     -1.0_rDef*A(i-4))     &
                            /(60.0_rDef*delta_x)
        
        i       = iMax-1
        dAdX(i) = -( -10.0_rDef*A(i+1)     &
                     -77.0_rDef*A(i+0)     &
                    +150.0_rDef*A(i-1)     &
                    -100.0_rDef*A(i-2)     &
                     +50.0_rDef*A(i-3)     &
                     -15.0_rDef*A(i-4)     &
                      +2.0_rDef*A(i-5))    &
                            /(60.0_rDef*delta_x)
        
        i       = iMax
        dAdX(i) = -(-147.0_rDef*A(i-0)    &
                    +360.0_rDef*A(i-1)    &
                    -450.0_rDef*A(i-2)    &
                    +400.0_rDef*A(i-3)    &
                    -225.0_rDef*A(i-4)    &
                     +72.0_rDef*A(i-5)    &
                     -10.0_rDef*A(i-6))   &
                            /(60.0_rDef*delta_x)
    ELSE IF (DS == 4) THEN

        i = iMin
        dAdX(i) = (-119.0_rDef*A(i+0)   &
                   +296.0_rDef*A(i+1)   &
                   -379.0_rDef*A(i+2)   &
                   +344.0_rDef*A(i+3)   &
                   -197.0_rDef*A(i+4)   &
                    +64.0_rDef*A(i+5)   &
                     -9.0_rDef*A(i+6))&
                            /(48.0_rDef*delta_x)
        i = iMin+1
        dAdX(i) =  ( -9.0_rDef*A(i-1)   &
                    -56.0_rDef*A(i+0)   &
                   +107.0_rDef*A(i+1)   &
                    -64.0_rDef*A(i+2)   &
                    +29.0_rDef*A(i+3)   &
                     -8.0_rDef*A(i+4)   &
                     +1.0_rDef*A(i+5))&
                            /(48.0_rDef*delta_x)
        i = iMin+2
        dAdX(i) = (  1.0_rDef*A(i-2)   &
                   -16.0_rDef*A(i-1)   &
                   -35.0_rDef*A(i+0)   &
                   +72.0_rDef*A(i+1)   &
                   -29.0_rDef*A(i+2)   &
                    +8.0_rDef*A(i+3)   &
                    -1.0_rDef*A(i+4))&
                         /(48.0_rDef*delta_x)

        DO i = iMin+3, iMax-3
            dAdX(i)   = (1.0_rDef/(delta_x*48.0_rDef))&
                                *(   -1.0_rDef*A(i-3)&
                                     +8.0_rDef*A(i-2)&
                                    -37.0_rDef*A(i-1)&
                                    +37.0_rDef*A(i+1)&
                                     -8.0_rDef*A(i+2)&
                                     +1.0_rDef*A(i+3))
        END DO
    
        i       = iMax-2
        dAdX(i) = -( 1.0_rDef*A(i+2)   &
                   -16.0_rDef*A(i+1)   &
                   -35.0_rDef*A(i+0)   &
                   +72.0_rDef*A(i-1)   &
                   -29.0_rDef*A(i-2)   &
                    +8.0_rDef*A(i-3)   &
                    -1.0_rDef*A(i-4))&
                         /(48.0_rDef*delta_x)
        
        i       = iMax-1
        dAdX(i) =  -(-9.0_rDef*A(i+1)   &
                    -56.0_rDef*A(i+0)   &
                   +107.0_rDef*A(i-1)   &
                    -64.0_rDef*A(i-2)   &
                    +29.0_rDef*A(i-3)   &
                     -8.0_rDef*A(i-4)   &
                     +1.0_rDef*A(i-5))&
                            /(48.0_rDef*delta_x)
        
        i       = iMax
        dAdX(i) = -(-119.0_rDef*A(i-0)   &
                    +296.0_rDef*A(i-1)   &
                    -379.0_rDef*A(i-2)   &
                    +344.0_rDef*A(i-3)   &
                    -197.0_rDef*A(i-4)   &
                     +64.0_rDef*A(i-5)   &
                      -9.0_rDef*A(i-6))&
                            /(48.0_rDef*delta_x)

    ELSE IF (DS == 5) THEN

        i = iMin
        dAdX(i) = (-2.192280339_rDef*A(i+0)    &
                   +4.748611401_rDef*A(i+1)    &
                   -5.108851915_rDef*A(i+2)    &
                   +4.461567104_rDef*A(i+3)    &
                   -2.833498741_rDef*A(i+4)    &
                   +1.128328861_rDef*A(i+5)    &
                   -0.203876371_rDef*A(i+6))/(delta_x)

        i = iMin + 1
        dAdX(i) = ( -0.20940637445552380_rDef*A(i-1)   &
                    -1.08406399497536900_rDef*A(i-0)   &
                    +2.14474893675524100_rDef*A(i+1)   &
                    -1.38356163287698400_rDef*A(i+2)   &
                    +0.76392684656023610_rDef*A(i+3)   &
                    -0.27940632043194790_rDef*A(i+4)   &
                    +0.04776253942434876_rDef*A(i+5))/(delta_x)

        i       = iMin + 2
        dAdX(i) = ( 0.056444818639386500_rDef*A(i-2)    &
                    -0.51192471192015800_rDef*A(i-1)    &
                    -0.37038205460987970_rDef*A(i-0)    &
                    +1.13854563013693000_rDef*A(i+1)    &
                    -0.42076972392956480_rDef*A(i+2)    &
                    +0.12838542345539630_rDef*A(i+3)    &
                    -0.02029938177210999_rDef*A(i+4))/(delta_x)

        DO i = iMin+3, iMax-3
            dAdX(i) = (-0.0208431427703_rDef*A(i-3)&
                       +0.1667059044150_rDef*A(i-2)&
                       -0.7708823805180_rDef*A(i-1)&
                       +0.7708823805180_rDef*A(i+1)&
                       -0.1667059044150_rDef*A(i+2)&
                       +0.0208431427703_rDef*A(i+3))/delta_x
        END DO

        i       = iMax - 2
        dAdX(i) = ( 0.056444818639386500_rDef*A(i+2)    &
                    -0.51192471192015800_rDef*A(i+1)    &
                    -0.37038205460987970_rDef*A(i-0)    &
                    +1.13854563013693000_rDef*A(i-1)    &
                    -0.42076972392956480_rDef*A(i-2)    &
                    +0.12838542345539630_rDef*A(i-3)    &
                    -0.02029938177210999_rDef*A(i-4))/(delta_x)

        i = iMax - 1
        dAdX(i) = ( -0.20940637445552380_rDef*A(i+1)   &
                    -1.08406399497536900_rDef*A(i-0)   &
                    +2.14474893675524100_rDef*A(i-1)   &
                    -1.38356163287698400_rDef*A(i-2)   &
                    +0.76392684656023610_rDef*A(i-3)   &
                    -0.27940632043194790_rDef*A(i-4)   &
                    +0.04776253942434876_rDef*A(i-5))/(delta_x)
    
        i = iMax
        dAdX(i) = (-2.192280339_rDef*A(i-0)    &
                   +4.748611401_rDef*A(i-1)    &
                   -5.108851915_rDef*A(i-2)    &
                   +4.461567104_rDef*A(i-3)    &
                   -2.833498741_rDef*A(i-4)    &
                   +1.128328861_rDef*A(i-5)    &
                   -0.203876371_rDef*A(i-6))/(delta_x)

    ELSE IF (DS == 6) THEN

        ALLOCATE(   LowerDiag(iMax, 1, 1)   ,&
                    UpperDiag(iMax, 1, 1)   ,&
                    MainDiag(iMax, 1, 1)    ,&
                    RHS(iMax, 1)            ,&
                    dAdX_TriDi(iMax, 1))

        i = iMin
        LowerDiag(i, 1, 1)  = 0.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 3.0_rDef
        RHS(i, 1)           =   ((-17.0_rDef/6.0_rDef)*A(i)     +       &
                                1.5_rDef*A(i + 1)               +       &
                                1.5_rDef*A(i + 2)               -       &
                                1.0_rDef/6.0_rDef*A(i + 3))/delta_x

        DO i = iMin + 1, iMax - 1
            LowerDiag(i, 1, 1)  = 0.25_rDef
            MainDiag(i, 1, 1)   = 1.0_rDef
            UpperDiag(i, 1, 1)  = 0.25_rDef
            RHS(i, 1)           = 0.75_rDef*(A(i + 1) - A(i - 1))/delta_x
        END DO

        i = iMax
        LowerDiag(i, 1, 1)  = 3.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 0.0_rDef
        RHS(i, 1)           =   ((17.0_rDef/6.0_rDef)*A(i)      -       &
                                1.5_rDef*A(i - 1)              -       &
                                1.5_rDef*A(i - 2)              +       &
                                1.0_rDef/6.0_rDef*A(i - 3))/delta_x

    
        CALL BlockTriMatrix(  numVar       = 1,         & 
                              iSolveStart  = iiMin,     &
                              iSolveEnd    = iiMax,     &
                              aStart       = iiMin,     &
                              A            = LowerDiag, &
                              bStart       = iiMin,     &
                              B            = MainDiag,  &
                              cStart       = iiMin,     &
                              C            = UpperDiag, &
                              xStart       = iiMin,     &
                              X            = dAdX_TriDi,&
                              rStart       = iiMin,     &
                              R            = RHS)

        DO i = iMin, iMax
            dAdX(i) = dAdX_TriDi(i, 1)
        END DO

    ELSE IF (DS == 7) THEN

        ALLOCATE(   LowerDiag(iMax, 1, 1)   ,&
                    UpperDiag(iMax, 1, 1)   ,&
                    MainDiag(iMax, 1, 1)    ,&
                    RHS(iMax, 1)            ,&
                    dAdX_TriDi(iMax, 1))

        i = iMin
        LowerDiag(i, 1, 1)  = 0.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 3.0_rDef
        RHS(i, 1)           =   ((-17.0_rDef/6.0_rDef)*A(i)     +       &
                                1.5_rDef*A(i + 1)               +       &
                                1.5_rDef*A(i + 2)               -       &
                                1.0_rDef/6.0_rDef*A(i + 3))/delta_x

        DO i = iMin + 1, iMax - 1
            LowerDiag(i, 1, 1)  = 0.25_rDef
            MainDiag(i, 1, 1)   = 1.0_rDef
            UpperDiag(i, 1, 1)  = 0.25_rDef
            RHS(i, 1)           = 0.75_rDef*(A(i + 1) - A(i - 1))/delta_x
        END DO

        i = iMax
        LowerDiag(i, 1, 1)  = 3.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 0.0_rDef
        RHS(i, 1)           =   ((17.0_rDef/6.0_rDef)*A(i)      -       &
                                1.5_rDef*A(i - 1)              -       &
                                1.5_rDef*A(i - 2)              +       &
                                1.0_rDef/6.0_rDef*A(i - 3))/delta_x

    
        CALL BlockTriMatrix(  numVar       = 1,         & 
                              iSolveStart  = iiMin,     &
                              iSolveEnd    = iiMax,     &
                              aStart       = iiMin,     &
                              A            = LowerDiag, &
                              bStart       = iiMin,     &
                              B            = MainDiag,  &
                              cStart       = iiMin,     &
                              C            = UpperDiag, &
                              xStart       = iiMin,     &
                              X            = dAdX_TriDi,&
                              rStart       = iiMin,     &
                              R            = RHS)

        DO i = iMin, iMax
            dAdX(i) = dAdX_TriDi(i, 1)
        END DO
    ELSE IF (DS == 8) THEN

        ALLOCATE(   LowerDiag(iMax, 1, 1)   ,&
                    UpperDiag(iMax, 1, 1)   ,&
                    MainDiag(iMax, 1, 1)    ,&
                    RHS(iMax, 1)            ,&
                    dAdX_TriDi(iMax, 1))

        i = iMin
        LowerDiag(i, 1, 1)  = 0.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 5.0_rDef
        RHS(i, 1)           =   ((-197.0_rDef/60.0_rDef)*A(i)   -   &
                                (5.0_rDef/12.0_rDef)*A(i + 1)   +   &
                                5.0_rDef*A(i + 2)               -   &
                                5.0_rDef/3.0_rDef*A(i + 3)      +   &
                                5.0_rDef/12.0_rDef*A(i + 4)     -   &
                                1.0_rDef/20.0_rDef*A(i + 5))/delta_x

        i = iMin + 1
        LowerDiag(i, 1, 1)  = 2.0_rDef/11.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 2.0_rDef/11.0_rDef
        RHS(i, 1)           =   ((-20.0_rDef/33.0_rDef)*A(i)    -   &
                                (35.0_rDef/132.0_rDef)*A(i + 1) +   &
                                (34.0_rDef/33.0_rDef)*A(i + 2)  -   &
                                (7.0_rDef/33.0_rDef)*A(i + 3)   +   &
                                (2.0_rDef/33.0_rDef)*A(i + 4)   -   &
                                (1.0_rDef/132.0_rDef)*A(i + 5))/delta_x

        DO i = iMin + 2, iMax - 2
            LowerDiag(i, 1, 1)  = 0.3333333333333_rDef
            MainDiag(i, 1, 1)   = 1.0_rDef
            UpperDiag(i, 1, 1)  = 0.3333333333333_rDef
            RHS(i, 1)           =   14.0_rDef/(9.0_rDef*2.0_rDef)       &
                                    *(A(i + 1) - A(i - 1))/delta_x +    &
                                    1.0_rDef/(9.0_rDef*4.0_rDef)        &
                                    *(A(i + 2) - A(i - 2))/delta_x
        END DO

        i = iMax - 1
        LowerDiag(i, 1, 1)  = 2.0_rDef/11.0_rDef
        MainDiag(i, 1, 1)   = 2.0_rDef/11.0_rDef
        UpperDiag(i, 1, 1)  = 0.0_rDef
        RHS(i, 1)           =  -((-20.0_rDef/33.0_rDef)*A(i)    -   &
                                (35.0_rDef/132.0_rDef)*A(i - 1) +   &
                                (34.0_rDef/33.0_rDef)*A(i - 2)  -   &
                                (7.0_rDef/33.0_rDef)*A(i - 3)   +   &
                                (2.0_rDef/33.0_rDef)*A(i - 4)   -   &
                                (1.0_rDef/132.0_rDef)*A(i - 5))/delta_x
        i = iMax
        LowerDiag(i, 1, 1)  = 5.0_rDef
        MainDiag(i, 1, 1)   = 1.0_rDef
        UpperDiag(i, 1, 1)  = 0.0_rDef
        RHS(i, 1)           =   (( 197.0_rDef/60.0_rDef)*A(i)   +   &
                                (5.0_rDef/12.0_rDef)*A(i - 1)   -   &
                                5.0_rDef*A(i - 2)               +   &
                                5.0_rDef/3.0_rDef*A(i - 3)      -   &
                                5.0_rDef/12.0_rDef*A(i - 4)     +   &
                                1.0_rDef/20.0_rDef*A(i - 5))/delta_x

    
        CALL BlockTriMatrix(  numVar       = 1,         & 
                              iSolveStart  = iiMin,     &
                              iSolveEnd    = iiMax,     &
                              aStart       = iiMin,     &
                              A            = LowerDiag, &
                              bStart       = iiMin,     &
                              B            = MainDiag,  &
                              cStart       = iiMin,     &
                              C            = UpperDiag, &
                              xStart       = iiMin,     &
                              X            = dAdX_TriDi,&
                              rStart       = iiMin,     &
                              R            = RHS)

        DO i = iMin, iMax
            dAdX(i) = dAdX_TriDi(i, 1)
        END DO

    END IF

    END SUBROUTINE AreaDerivative

END MODULE GetAreaDerivative
