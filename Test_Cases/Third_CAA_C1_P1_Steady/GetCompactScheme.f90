MODULE GetCompactScheme
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE BlockTriDiSolver1D
    USE GetSolvePreFac
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: CompactScheme

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE CompactScheme(   iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                E               ,&
                                dEdX            ,&
                                delta_x         ,&
                                delta_tau       ,&
                                DS              ,&
                                dZidX)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS

    REAL(KIND = rDef), INTENT(IN):: delta_tau  ,&
                                     delta_x

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN):: dZidX

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: E

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: dEdX

    INTEGER:: i, j, k

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   LowerDiag       ,&
                                                            UpperDiag       ,&
                                                            MainDiag

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  RHS     ,&
                                                        RHS_B   ,&
                                                        RHS_F   ,&
                                                        dEdX_B     ,&
                                                        dEdX_F     ,&
                                                        Unknown 

    REAL(KIND = rDef) ::    a          ,&
                            b

    INTEGER, DIMENSION(1):: iiMin, iiMax

    iiMin = iMin
    iiMax = iMax
    
    ALLOCATE(   LowerDiag(iMax, bMax, bMax) ,&
                UpperDiag(iMax, bMax, bMax) ,&
                MainDiag(iMax, bMax, bMax)  ,&
                RHS(iMax, bMax)             ,&
                dEdX_B(iMax, bMax)          ,&
                dEdX_F(iMax, bMax)          ,&
                RHS_B(iMax, bMax)           ,&
                RHS_F(iMax, bMax)           ,&
                Unknown(iMax, bMax))

    IF (DS == 6) THEN
        i           = iMin
        LowerDiag(i, 1,  1)  = 0.0_rDef 
        LowerDiag(i, 1,  2)  = 0.0_rDef 
        LowerDiag(i, 1,  3)  = 0.0_rDef 
        LowerDiag(i, 2,  1)  = 0.0_rDef 
        LowerDiag(i, 2,  2)  = 0.0_rDef 
        LowerDiag(i, 2,  3)  = 0.0_rDef 
        LowerDiag(i, 3,  1)  = 0.0_rDef 
        LowerDiag(i, 3,  2)  = 0.0_rDef 
        LowerDiag(i, 3,  3)  = 0.0_rDef 

        MainDiag(i, 1,  1)  = 1.0_rDef 
        MainDiag(i, 1,  2)  = 0.0_rDef 
        MainDiag(i, 1,  3)  = 0.0_rDef 
        MainDiag(i, 2,  1)  = 0.0_rDef 
        MainDiag(i, 2,  2)  = 1.0_rDef 
        MainDiag(i, 2,  3)  = 0.0_rDef 
        MainDiag(i, 3,  1)  = 0.0_rDef 
        MainDiag(i, 3,  2)  = 0.0_rDef 
        MainDiag(i, 3,  3)  = 1.0_rDef 

        UpperDiag(i, 1,  1)  = 0.0_rDef 
        UpperDiag(i, 1,  2)  = 0.0_rDef 
        UpperDiag(i, 1,  3)  = 0.0_rDef 
        UpperDiag(i, 2,  1)  = 0.0_rDef 
        UpperDiag(i, 2,  2)  = 0.0_rDef 
        UpperDiag(i, 2,  3)  = 0.0_rDef 
        UpperDiag(i, 3,  1)  = 0.0_rDef 
        UpperDiag(i, 3,  2)  = 0.0_rDef 
        UpperDiag(i, 3,  3)  = 0.0_rDef 

        DO i = iMin + 1, iMax - 1
            LowerDiag(i, 1,  1)  = 0.25_rDef 
            LowerDiag(i, 1,  2)  = 0.0_rDef 
            LowerDiag(i, 1,  3)  = 0.0_rDef 
            LowerDiag(i, 2,  1)  = 0.0_rDef 
            LowerDiag(i, 2,  2)  = 0.25_rDef 
            LowerDiag(i, 2,  3)  = 0.0_rDef 
            LowerDiag(i, 3,  1)  = 0.0_rDef 
            LowerDiag(i, 3,  2)  = 0.0_rDef 
            LowerDiag(i, 3,  3)  = 0.25_rDef 

            MainDiag(i, 1,  1)  = 1.0_rDef 
            MainDiag(i, 1,  2)  = 0.0_rDef 
            MainDiag(i, 1,  3)  = 0.0_rDef 
            MainDiag(i, 2,  1)  = 0.0_rDef 
            MainDiag(i, 2,  2)  = 1.0_rDef 
            MainDiag(i, 2,  3)  = 0.0_rDef 
            MainDiag(i, 3,  1)  = 0.0_rDef 
            MainDiag(i, 3,  2)  = 0.0_rDef 
            MainDiag(i, 3,  3)  = 1.0_rDef 

            UpperDiag(i, 1,  1)  = 0.25_rDef 
            UpperDiag(i, 1,  2)  = 0.0_rDef 
            UpperDiag(i, 1,  3)  = 0.0_rDef 
            UpperDiag(i, 2,  1)  = 0.0_rDef 
            UpperDiag(i, 2,  2)  = 0.25_rDef 
            UpperDiag(i, 2,  3)  = 0.0_rDef 
            UpperDiag(i, 3,  1)  = 0.0_rDef 
            UpperDiag(i, 3,  2)  = 0.0_rDef 
            UpperDiag(i, 3,  3)  = 0.25_rDef 
        END DO

        i           = iMax
        LowerDiag(i, 1,  1)  = 0.0_rDef 
        LowerDiag(i, 1,  2)  = 0.0_rDef 
        LowerDiag(i, 1,  3)  = 0.0_rDef 
        LowerDiag(i, 2,  1)  = 0.0_rDef 
        LowerDiag(i, 2,  2)  = 0.0_rDef 
        LowerDiag(i, 2,  3)  = 0.0_rDef 
        LowerDiag(i, 3,  1)  = 0.0_rDef 
        LowerDiag(i, 3,  2)  = 0.0_rDef 
        LowerDiag(i, 3,  3)  = 0.0_rDef 

        MainDiag(i, 1,  1)  = 1.0_rDef 
        MainDiag(i, 1,  2)  = 0.0_rDef 
        MainDiag(i, 1,  3)  = 0.0_rDef 
        MainDiag(i, 2,  1)  = 0.0_rDef 
        MainDiag(i, 2,  2)  = 1.0_rDef 
        MainDiag(i, 2,  3)  = 0.0_rDef 
        MainDiag(i, 3,  1)  = 0.0_rDef 
        MainDiag(i, 3,  2)  = 0.0_rDef 
        MainDiag(i, 3,  3)  = 1.0_rDef 

        UpperDiag(i, 1,  1)  = 0.0_rDef 
        UpperDiag(i, 1,  2)  = 0.0_rDef 
        UpperDiag(i, 1,  3)  = 0.0_rDef 
        UpperDiag(i, 2,  1)  = 0.0_rDef 
        UpperDiag(i, 2,  2)  = 0.0_rDef 
        UpperDiag(i, 2,  3)  = 0.0_rDef 
        UpperDiag(i, 3,  1)  = 0.0_rDef 
        UpperDiag(i, 3,  2)  = 0.0_rDef 
        UpperDiag(i, 3,  3)  = 0.0_rDef 

        DO j = bMin, bMax 
            i           = iMin
            RHS(i, j)   =   dEdX(i, j)
            DO i = iMin + 1, iMax - 1
                RHS(i, j) = 0.75_rDef*(E(i + 1, j) - E(i - 1, j))/delta_x
            END DO
            i           = iMax
            RHS(i, j)   =   dEdX(i, j)
        END DO

        CALL BlockTriMatrix(  numVar       = bMax,      & 
                              iSolveStart  = iiMin,     &
                              iSolveEnd    = iiMax,     &
                              aStart       = iiMin,     &
                              A            = LowerDiag, &
                              bStart       = iiMin,     &
                              B            = MainDiag,  &
                              cStart       = iiMin,     &
                              C            = UpperDiag, &
                              xStart       = iiMin,     &
                              X            = Unknown,   &
                              rStart       = iiMin,     &
                              R            = RHS)

        DO j = bMin, bMax 
            DO i = iMin + 1, iMax - 1
                dEdX(i, j) = Unknown(i, j)
            END DO
        END DO

    ELSEIF (DS == 7) THEN
        
        a   = 0.5_rDef - 0.288675134594813_rDef 
        b   = 1.0_rDef

        i           = iMin
        LowerDiag(i, 1,  1)  = 0.0_rDef 
        LowerDiag(i, 1,  2)  = 0.0_rDef 
        LowerDiag(i, 1,  3)  = 0.0_rDef 
        LowerDiag(i, 2,  1)  = 0.0_rDef 
        LowerDiag(i, 2,  2)  = 0.0_rDef 
        LowerDiag(i, 2,  3)  = 0.0_rDef 
        LowerDiag(i, 3,  1)  = 0.0_rDef 
        LowerDiag(i, 3,  2)  = 0.0_rDef 
        LowerDiag(i, 3,  3)  = 0.0_rDef 

        MainDiag(i, 1,  1)  = 0.5_rDef 
        MainDiag(i, 1,  2)  = 0.0_rDef 
        MainDiag(i, 1,  3)  = 0.0_rDef 
        MainDiag(i, 2,  1)  = 0.0_rDef 
        MainDiag(i, 2,  2)  = 0.5_rDef 
        MainDiag(i, 2,  3)  = 0.0_rDef 
        MainDiag(i, 3,  1)  = 0.0_rDef 
        MainDiag(i, 3,  2)  = 0.0_rDef 
        MainDiag(i, 3,  3)  = 0.5_rDef 

        UpperDiag(i, 1,  1)  = 0.0_rDef 
        UpperDiag(i, 1,  2)  = 0.0_rDef 
        UpperDiag(i, 1,  3)  = 0.0_rDef 
        UpperDiag(i, 2,  1)  = 0.0_rDef 
        UpperDiag(i, 2,  2)  = 0.0_rDef 
        UpperDiag(i, 2,  3)  = 0.0_rDef 
        UpperDiag(i, 3,  1)  = 0.0_rDef 
        UpperDiag(i, 3,  2)  = 0.0_rDef 
        UpperDiag(i, 3,  3)  = 0.0_rDef 

        DO i = iMin + 1, iMax - 1
            LowerDiag(i, 1,  1)  = 0.5_rDef*(a/(1.0_rDef - a))  
            LowerDiag(i, 1,  2)  = 0.0_rDef                     
            LowerDiag(i, 1,  3)  = 0.0_rDef                     
            LowerDiag(i, 2,  1)  = 0.0_rDef                     
            LowerDiag(i, 2,  2)  = 0.5_rDef*(a/(1.0_rDef - a))  
            LowerDiag(i, 2,  3)  = 0.0_rDef                     
            LowerDiag(i, 3,  1)  = 0.0_rDef                     
            LowerDiag(i, 3,  2)  = 0.0_rDef                     
            LowerDiag(i, 3,  3)  = 0.5_rDef*(a/(1.0_rDef - a))  

            MainDiag(i, 1,  1)  = 0.5_rDef 
            MainDiag(i, 1,  2)  = 0.0_rDef 
            MainDiag(i, 1,  3)  = 0.0_rDef 
            MainDiag(i, 2,  1)  = 0.0_rDef 
            MainDiag(i, 2,  2)  = 0.5_rDef 
            MainDiag(i, 2,  3)  = 0.0_rDef 
            MainDiag(i, 3,  1)  = 0.0_rDef 
            MainDiag(i, 3,  2)  = 0.0_rDef 
            MainDiag(i, 3,  3)  = 0.5_rDef 

            UpperDiag(i, 1,  1)  = 0.5_rDef*(a/(1.0_rDef - a))
            UpperDiag(i, 1,  2)  = 0.0_rDef 
            UpperDiag(i, 1,  3)  = 0.0_rDef 
            UpperDiag(i, 2,  1)  = 0.0_rDef 
            UpperDiag(i, 2,  2)  = 0.5_rDef*(a/(1.0_rDef - a)) 
            UpperDiag(i, 2,  3)  = 0.0_rDef 
            UpperDiag(i, 3,  1)  = 0.0_rDef 
            UpperDiag(i, 3,  2)  = 0.0_rDef 
            UpperDiag(i, 3,  3)  = 0.5_rDef*(a/(1.0_rDef - a)) 
        END DO

        i           = iMax
        LowerDiag(i, 1,  1)  = 0.0_rDef 
        LowerDiag(i, 1,  2)  = 0.0_rDef 
        LowerDiag(i, 1,  3)  = 0.0_rDef 
        LowerDiag(i, 2,  1)  = 0.0_rDef 
        LowerDiag(i, 2,  2)  = 0.0_rDef 
        LowerDiag(i, 2,  3)  = 0.0_rDef 
        LowerDiag(i, 3,  1)  = 0.0_rDef 
        LowerDiag(i, 3,  2)  = 0.0_rDef 
        LowerDiag(i, 3,  3)  = 0.0_rDef 

        MainDiag(i, 1,  1)  = 0.5_rDef 
        MainDiag(i, 1,  2)  = 0.0_rDef 
        MainDiag(i, 1,  3)  = 0.0_rDef 
        MainDiag(i, 2,  1)  = 0.0_rDef 
        MainDiag(i, 2,  2)  = 0.5_rDef 
        MainDiag(i, 2,  3)  = 0.0_rDef 
        MainDiag(i, 3,  1)  = 0.0_rDef 
        MainDiag(i, 3,  2)  = 0.0_rDef 
        MainDiag(i, 3,  3)  = 0.5_rDef 

        UpperDiag(i, 1,  1)  = 0.0_rDef 
        UpperDiag(i, 1,  2)  = 0.0_rDef 
        UpperDiag(i, 1,  3)  = 0.0_rDef 
        UpperDiag(i, 2,  1)  = 0.0_rDef 
        UpperDiag(i, 2,  2)  = 0.0_rDef 
        UpperDiag(i, 2,  3)  = 0.0_rDef 
        UpperDiag(i, 3,  1)  = 0.0_rDef 
        UpperDiag(i, 3,  2)  = 0.0_rDef 
        UpperDiag(i, 3,  3)  = 0.0_rDef 

        DO j = bMin, bMax 
            i           = iMin
            RHS_B(i, j) = 0.5_rDef*dEdX(i, j)
            RHS_F(i, j) = 0.5_rDef*dEdX(i, j)
            
            DO i = iMin + 1, iMax - 1
                RHS_B(i, j) = &
                    (b/(2.0_rDef*delta_x*(1.0_rDef - a)))*              &
                    (E(i + 1, j) - E(i, j))                         +   &
                    ((1.0_rDef - b)/(2.0_rDef*delta_x*(1.0_rDef - a)))* &
                    (E(i, j) - E(i - 1, j))

                RHS_F(i, j) = &
                    (b/(2.0_rDef*delta_x*(1.0_rDef - a)))*              &
                    (E(i, j) - E(i - 1, j))                         +   &
                    ((1.0_rDef - b)/(2.0_rDef*delta_x*(1.0_rDef - a)))* &
                    (E(i + 1, j) - E(i, j))
            END DO

            i           = iMax
            RHS_B(i, j) = 0.5_rDef*dEdX(i, j)
            RHS_F(i, j) = 0.5_rDef*dEdX(i, j)
        END DO

        CALL SolvePreFac(   L_LowerDiag = LowerDiag     ,&
                            U_UpperDiag = UpperDiag     ,&
                            Tau         = MainDiag      ,&
                            iMin        = iMin          ,&
                            iMax        = iMax          ,&
                            bMin        = bMin          ,&
                            bMax        = bMax          ,&
                            RHS_B       = RHS_B         ,&
                            RHS_F       = RHS_F         ,&
                            D_B         = dEdX_B        ,&
                            D_F         = dEdX_F)

        DO j = bMin, bMax 
            DO i = iMin + 1, iMax - 1
                dEdX(i, j) = (dEdX_B(i, j) + dEdX_F(i, j))*0.5_rDef
            END DO
        END DO

    ELSEIF (DS == 8) THEN
        
        a   = 0.5_rDef - 0.223606797749979_rDef 
        b   = 1.0_rDef - 1.0_rDef/(30.0_rDef*a)

        i           = iMin
        LowerDiag(i, 1,  1)  = 0.0_rDef 
        LowerDiag(i, 1,  2)  = 0.0_rDef 
        LowerDiag(i, 1,  3)  = 0.0_rDef 
        LowerDiag(i, 2,  1)  = 0.0_rDef 
        LowerDiag(i, 2,  2)  = 0.0_rDef 
        LowerDiag(i, 2,  3)  = 0.0_rDef 
        LowerDiag(i, 3,  1)  = 0.0_rDef 
        LowerDiag(i, 3,  2)  = 0.0_rDef 
        LowerDiag(i, 3,  3)  = 0.0_rDef 

        MainDiag(i, 1,  1)  = 0.5_rDef 
        MainDiag(i, 1,  2)  = 0.0_rDef 
        MainDiag(i, 1,  3)  = 0.0_rDef 
        MainDiag(i, 2,  1)  = 0.0_rDef 
        MainDiag(i, 2,  2)  = 0.5_rDef 
        MainDiag(i, 2,  3)  = 0.0_rDef 
        MainDiag(i, 3,  1)  = 0.0_rDef 
        MainDiag(i, 3,  2)  = 0.0_rDef 
        MainDiag(i, 3,  3)  = 0.5_rDef 

        UpperDiag(i, 1,  1)  = 0.0_rDef 
        UpperDiag(i, 1,  2)  = 0.0_rDef 
        UpperDiag(i, 1,  3)  = 0.0_rDef 
        UpperDiag(i, 2,  1)  = 0.0_rDef 
        UpperDiag(i, 2,  2)  = 0.0_rDef 
        UpperDiag(i, 2,  3)  = 0.0_rDef 
        UpperDiag(i, 3,  1)  = 0.0_rDef 
        UpperDiag(i, 3,  2)  = 0.0_rDef 
        UpperDiag(i, 3,  3)  = 0.0_rDef 

        DO i = iMin + 1, iMax - 1
            LowerDiag(i, 1,  1)  = 0.5_rDef*(a/(1.0_rDef - a))  
            LowerDiag(i, 1,  2)  = 0.0_rDef                     
            LowerDiag(i, 1,  3)  = 0.0_rDef                     
            LowerDiag(i, 2,  1)  = 0.0_rDef                     
            LowerDiag(i, 2,  2)  = 0.5_rDef*(a/(1.0_rDef - a))  
            LowerDiag(i, 2,  3)  = 0.0_rDef                     
            LowerDiag(i, 3,  1)  = 0.0_rDef                     
            LowerDiag(i, 3,  2)  = 0.0_rDef                     
            LowerDiag(i, 3,  3)  = 0.5_rDef*(a/(1.0_rDef - a))  

            MainDiag(i, 1,  1)  = 0.5_rDef 
            MainDiag(i, 1,  2)  = 0.0_rDef 
            MainDiag(i, 1,  3)  = 0.0_rDef 
            MainDiag(i, 2,  1)  = 0.0_rDef 
            MainDiag(i, 2,  2)  = 0.5_rDef 
            MainDiag(i, 2,  3)  = 0.0_rDef 
            MainDiag(i, 3,  1)  = 0.0_rDef 
            MainDiag(i, 3,  2)  = 0.0_rDef 
            MainDiag(i, 3,  3)  = 0.5_rDef 

            UpperDiag(i, 1,  1)  = 0.5_rDef*(a/(1.0_rDef - a))
            UpperDiag(i, 1,  2)  = 0.0_rDef 
            UpperDiag(i, 1,  3)  = 0.0_rDef 
            UpperDiag(i, 2,  1)  = 0.0_rDef 
            UpperDiag(i, 2,  2)  = 0.5_rDef*(a/(1.0_rDef - a)) 
            UpperDiag(i, 2,  3)  = 0.0_rDef 
            UpperDiag(i, 3,  1)  = 0.0_rDef 
            UpperDiag(i, 3,  2)  = 0.0_rDef 
            UpperDiag(i, 3,  3)  = 0.5_rDef*(a/(1.0_rDef - a)) 
        END DO

        i           = iMax
        LowerDiag(i, 1,  1)  = 0.0_rDef 
        LowerDiag(i, 1,  2)  = 0.0_rDef 
        LowerDiag(i, 1,  3)  = 0.0_rDef 
        LowerDiag(i, 2,  1)  = 0.0_rDef 
        LowerDiag(i, 2,  2)  = 0.0_rDef 
        LowerDiag(i, 2,  3)  = 0.0_rDef 
        LowerDiag(i, 3,  1)  = 0.0_rDef 
        LowerDiag(i, 3,  2)  = 0.0_rDef 
        LowerDiag(i, 3,  3)  = 0.0_rDef 

        MainDiag(i, 1,  1)  = 0.5_rDef 
        MainDiag(i, 1,  2)  = 0.0_rDef 
        MainDiag(i, 1,  3)  = 0.0_rDef 
        MainDiag(i, 2,  1)  = 0.0_rDef 
        MainDiag(i, 2,  2)  = 0.5_rDef 
        MainDiag(i, 2,  3)  = 0.0_rDef 
        MainDiag(i, 3,  1)  = 0.0_rDef 
        MainDiag(i, 3,  2)  = 0.0_rDef 
        MainDiag(i, 3,  3)  = 0.5_rDef 

        UpperDiag(i, 1,  1)  = 0.0_rDef 
        UpperDiag(i, 1,  2)  = 0.0_rDef 
        UpperDiag(i, 1,  3)  = 0.0_rDef 
        UpperDiag(i, 2,  1)  = 0.0_rDef 
        UpperDiag(i, 2,  2)  = 0.0_rDef 
        UpperDiag(i, 2,  3)  = 0.0_rDef 
        UpperDiag(i, 3,  1)  = 0.0_rDef 
        UpperDiag(i, 3,  2)  = 0.0_rDef 
        UpperDiag(i, 3,  3)  = 0.0_rDef 

        DO j = bMin, bMax 
            i           = iMin
            RHS_B(i, j) = 0.5_rDef*dEdX(i, j)
            RHS_F(i, j) = 0.5_rDef*dEdX(i, j)
            
            DO i = iMin + 1, iMax - 1
                RHS_B(i, j) = &
                    (b/(2.0_rDef*delta_x*(1.0_rDef - a)))*              &
                    (E(i + 1, j) - E(i, j))                         +   &
                    ((1.0_rDef - b)/(2.0_rDef*delta_x*(1.0_rDef - a)))* &
                    (E(i, j) - E(i - 1, j))

                RHS_F(i, j) = &
                    (b/(2.0_rDef*delta_x*(1.0_rDef - a)))*              &
                    (E(i, j) - E(i - 1, j))                         +   &
                    ((1.0_rDef - b)/(2.0_rDef*delta_x*(1.0_rDef - a)))* &
                    (E(i + 1, j) - E(i, j))
            END DO

            i           = iMax
            RHS_B(i, j) = 0.5_rDef*dEdX(i, j)
            RHS_F(i, j) = 0.5_rDef*dEdX(i, j)
        END DO

        CALL SolvePreFac(   L_LowerDiag = LowerDiag     ,&
                            U_UpperDiag = UpperDiag     ,&
                            Tau         = MainDiag      ,&
                            iMin        = iMin          ,&
                            iMax        = iMax          ,&
                            bMin        = bMin          ,&
                            bMax        = bMax          ,&
                            RHS_B       = RHS_B         ,&
                            RHS_F       = RHS_F         ,&
                            D_B         = dEdX_B        ,&
                            D_F         = dEdX_F)

        DO j = bMin, bMax 
            DO i = iMin + 1, iMax - 1
                dEdX(i, j) = (dEdX_B(i, j) + dEdX_F(i, j))*0.5_rDef
            END DO
        END DO

    END IF
    
    END SUBROUTINE CompactScheme

END MODULE GetCompactScheme
