MODULE GetSolvePreFac
!!! Working Test Case
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: SolvePreFac

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE SolvePreFac( L_LowerDiag     ,&    
                            U_UpperDiag     ,& 
                            Tau             ,&
                            iMin            ,& 
                            iMax            ,& 
                            bMin            ,& 
                            bMax            ,& 
                            RHS_B           ,&
                            RHS_F           ,&
                            D_B             ,&
                            D_F) 

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    L_LowerDiag ,&
                                                            U_UpperDiag ,&
                                                            Tau

    INTEGER, INTENT(IN) ::  iMin    ,&
                            iMax    ,&
                            bMin    ,&
                            bMax

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS_B   ,&
                                                            RHS_F   ,&
                                                            D_B     ,&
                                                            D_F

    INTEGER :: i, j, ii, m, n, l
    
    REAL(KIND = rDef) :: FacL, DetInv, Fac3
    
    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  A,              &
                                                        B,              &
                                                        EE,             &
                                                        LowerFac

    ALLOCATE(LowerFac(bMax, bMin)           ,&
             A(bMax, bMax)                  ,&
             B(bMax, bMin)                  ,&
             EE(bMax, bMin))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! I'm defining B as the RHS
    !! Then I will find the inverse of the block Matrix Tau(:, :, iMax)
    !! Then use matmul to find Delta_Q_Star2(iMin, :)
    i = iMin
    B(1, 1) = RHS_B(i, 1)
    B(2, 1) = RHS_B(i, 2)
    B(3, 1) = RHS_B(i, 3)

    Fac3 = (Tau(i, 1, 1)*Tau(i, 2, 2)*Tau(i, 3, 3) - &
            Tau(i, 1, 1)*Tau(i, 2, 3)*Tau(i, 3, 2) - &
            Tau(i, 1, 2)*Tau(i, 2, 1)*Tau(i, 3, 3) + &
            Tau(i, 1, 2)*Tau(i, 2, 3)*Tau(i, 3, 1) + &
            Tau(i, 1, 3)*Tau(i, 2, 1)*Tau(i, 3, 2) - &
            Tau(i, 1, 3)*Tau(i, 2, 2)*Tau(i, 3, 1))

    !WRITE(6, *) Tau(3, 3, 1),Tau(2, 2, 1),Tau(1, 1, 1), "Got Here"
    !READ(5, *)

    DetInv  = 1.0_rDef/Fac3
    A(1, 1) =  DetInv*(Tau(i, 2, 2)*Tau(i, 3, 3) - Tau(i, 2, 3)*Tau(i, 3, 2)) 
    A(2, 1) = -DetInv*(Tau(i, 2, 1)*Tau(i, 3, 3) - Tau(i, 2, 3)*Tau(i, 3, 1)) 
    A(3, 1) =  DetInv*(Tau(i, 2, 1)*Tau(i, 3, 2) - Tau(i, 2, 2)*Tau(i, 3, 1)) 
    A(1, 2) = -DetInv*(Tau(i, 1, 2)*Tau(i, 3, 3) - Tau(i, 1, 3)*Tau(i, 3, 2)) 
    A(2, 2) =  DetInv*(Tau(i, 1, 1)*Tau(i, 3, 3) - Tau(i, 1, 3)*Tau(i, 3, 1)) 
    A(3, 2) = -DetInv*(Tau(i, 1, 1)*Tau(i, 3, 2) - Tau(i, 1, 2)*Tau(i, 3, 1)) 
    A(1, 3) =  DetInv*(Tau(i, 1, 2)*Tau(i, 2, 3) - Tau(i, 1, 3)*Tau(i, 2, 2)) 
    A(2, 3) = -DetInv*(Tau(i, 1, 1)*Tau(i, 2, 3) - Tau(i, 1, 3)*Tau(i, 2, 1)) 
    A(3, 3) =  DetInv*(Tau(i, 1, 1)*Tau(i, 2, 2) - Tau(i, 1, 2)*Tau(i, 2, 1)) 

    EE = MATMUL(A, B)

    D_B(i, 1) = EE(1, 1)
    D_B(i, 2) = EE(2, 1)
    D_B(i, 3) = EE(3, 1)

    DO i = iMin + 1, iMax
        A(1, 1)     = L_LowerDiag(i, 1, 1)
        A(2, 1)     = L_LowerDiag(i, 2, 1)
        A(3, 1)     = L_LowerDiag(i, 3, 1)
        A(1, 2)     = L_LowerDiag(i, 1, 2)
        A(2, 2)     = L_LowerDiag(i, 2, 2)
        A(3, 2)     = L_LowerDiag(i, 3, 2)
        A(1, 3)     = L_LowerDiag(i, 1, 3)
        A(2, 3)     = L_LowerDiag(i, 2, 3)
        A(3, 3)     = L_LowerDiag(i, 3, 3)

        B(1, 1)     = D_B(i - 1, 1)
        B(2, 1)     = D_B(i - 1, 2)
        B(3, 1)     = D_B(i - 1, 3)

        LowerFac    = MATMUL(A, B)

        B(1, 1) = RHS_B(i, 1) - LowerFac(1, 1)
        B(2, 1) = RHS_B(i, 2) - LowerFac(2, 1)
        B(3, 1) = RHS_B(i, 3) - LowerFac(3, 1)

        !! The Inverse of Tau is Here
        Fac3 = (Tau(i, 1, 1)*Tau(i, 2, 2)*Tau(i, 3, 3) - &
                Tau(i, 1, 1)*Tau(i, 2, 3)*Tau(i, 3, 2) - &
                Tau(i, 1, 2)*Tau(i, 2, 1)*Tau(i, 3, 3) + &
                Tau(i, 1, 2)*Tau(i, 2, 3)*Tau(i, 3, 1) + &
                Tau(i, 1, 3)*Tau(i, 2, 1)*Tau(i, 3, 2) - &
                Tau(i, 1, 3)*Tau(i, 2, 2)*Tau(i, 3, 1))

        DetInv  = 1.0_rDef/Fac3
        
        A(1, 1) =  DetInv*(Tau(i, 2, 2)*Tau(i, 3, 3) - &
                                    Tau(i, 2, 3)*Tau(i, 3, 2)) 
        A(2, 1) = -DetInv*(Tau(i, 2, 1)*Tau(i, 3, 3) - &
                                    Tau(i, 2, 3)*Tau(i, 3, 1)) 
        A(3, 1) =  DetInv*(Tau(i, 2, 1)*Tau(i, 3, 2) - &
                                    Tau(i, 2, 2)*Tau(i, 3, 1)) 
        A(1, 2) = -DetInv*(Tau(i, 1, 2)*Tau(i, 3, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 3, 2)) 
        A(2, 2) =  DetInv*(Tau(i, 1, 1)*Tau(i, 3, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 3, 1)) 
        A(3, 2) = -DetInv*(Tau(i, 1, 1)*Tau(i, 3, 2) - &
                                    Tau(i, 1, 2)*Tau(i, 3, 1)) 
        A(1, 3) =  DetInv*(Tau(i, 1, 2)*Tau(i, 2, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 2, 2)) 
        A(2, 3) = -DetInv*(Tau(i, 1, 1)*Tau(i, 2, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 2, 1)) 
        A(3, 3) =  DetInv*(Tau(i, 1, 1)*Tau(i, 2, 2) - &
                                       Tau(i, 1, 2)*Tau(i, 2, 1)) 
        !!

        EE = MATMUL(A, B)

        D_B(i, 1) = EE(1, 1)
        D_B(i, 2) = EE(2, 1)
        D_B(i, 3) = EE(3, 1)
    END DO


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i = iMax
    B(1, 1) = RHS_F(i, 1)
    B(2, 1) = RHS_F(i, 2)
    B(3, 1) = RHS_F(i, 3)

    Fac3 = (Tau(i, 1, 1)*Tau(i, 2, 2)*Tau(i, 3, 3) - &
            Tau(i, 1, 1)*Tau(i, 2, 3)*Tau(i, 3, 2) - &
            Tau(i, 1, 2)*Tau(i, 2, 1)*Tau(i, 3, 3) + &
            Tau(i, 1, 2)*Tau(i, 2, 3)*Tau(i, 3, 1) + &
            Tau(i, 1, 3)*Tau(i, 2, 1)*Tau(i, 3, 2) - &
            Tau(i, 1, 3)*Tau(i, 2, 2)*Tau(i, 3, 1))

    DetInv  = 1.0_rDef/Fac3
    A(1, 1) =  DetInv*(Tau(i, 2, 2)*Tau(i, 3, 3) - Tau(i, 2, 3)*Tau(i, 3, 2)) 
    A(2, 1) = -DetInv*(Tau(i, 2, 1)*Tau(i, 3, 3) - Tau(i, 2, 3)*Tau(i, 3, 1)) 
    A(3, 1) =  DetInv*(Tau(i, 2, 1)*Tau(i, 3, 2) - Tau(i, 2, 2)*Tau(i, 3, 1)) 
    A(1, 2) = -DetInv*(Tau(i, 1, 2)*Tau(i, 3, 3) - Tau(i, 1, 3)*Tau(i, 3, 2)) 
    A(2, 2) =  DetInv*(Tau(i, 1, 1)*Tau(i, 3, 3) - Tau(i, 1, 3)*Tau(i, 3, 1)) 
    A(3, 2) = -DetInv*(Tau(i, 1, 1)*Tau(i, 3, 2) - Tau(i, 1, 2)*Tau(i, 3, 1)) 
    A(1, 3) =  DetInv*(Tau(i, 1, 2)*Tau(i, 2, 3) - Tau(i, 1, 3)*Tau(i, 2, 2)) 
    A(2, 3) = -DetInv*(Tau(i, 1, 1)*Tau(i, 2, 3) - Tau(i, 1, 3)*Tau(i, 2, 1)) 
    A(3, 3) =  DetInv*(Tau(i, 1, 1)*Tau(i, 2, 2) - Tau(i, 1, 2)*Tau(i, 2, 1)) 

    EE      = MATMUL(A, B)

    D_F(i, 1) = EE(1, 1) 
    D_F(i, 2) = EE(2, 1) 
    D_F(i, 3) = EE(3, 1) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Third Step
    DO i = iMax - 1, iMin, -iMin
        A(1, 1)     = U_UpperDiag(i, 1, 1)
        A(2, 1)     = U_UpperDiag(i, 2, 1)
        A(3, 1)     = U_UpperDiag(i, 3, 1)
        A(1, 2)     = U_UpperDiag(i, 1, 2)
        A(2, 2)     = U_UpperDiag(i, 2, 2)
        A(3, 2)     = U_UpperDiag(i, 3, 2)
        A(1, 3)     = U_UpperDiag(i, 1, 3)
        A(2, 3)     = U_UpperDiag(i, 2, 3)
        A(3, 3)     = U_UpperDiag(i, 3, 3)

        B(1, 1)     = D_F(i + 1, 1)
        B(2, 1)     = D_F(i + 1, 2)
        B(3, 1)     = D_F(i + 1, 3)

        LowerFac    = MATMUL(A, B)

        B(1, 1) = RHS_F(i, 1) - LowerFac(1, 1)
        B(2, 1) = RHS_F(i, 2) - LowerFac(2, 1)
        B(3, 1) = RHS_F(i, 3) - LowerFac(3, 1)

        !! The Inverse of Tau is Here
        Fac3 = (Tau(i, 1, 1)*Tau(i, 2, 2)*Tau(i, 3, 3) - &
                Tau(i, 1, 1)*Tau(i, 2, 3)*Tau(i, 3, 2) - &
                Tau(i, 1, 2)*Tau(i, 2, 1)*Tau(i, 3, 3) + &
                Tau(i, 1, 2)*Tau(i, 2, 3)*Tau(i, 3, 1) + &
                Tau(i, 1, 3)*Tau(i, 2, 1)*Tau(i, 3, 2) - &
                Tau(i, 1, 3)*Tau(i, 2, 2)*Tau(i, 3, 1))

        DetInv  = 1.0_rDef/Fac3
        
        A(1, 1) =  DetInv*(Tau(i, 2, 2)*Tau(i, 3, 3) - &
                                    Tau(i, 2, 3)*Tau(i, 3, 2)) 
        A(2, 1) = -DetInv*(Tau(i, 2, 1)*Tau(i, 3, 3) - &
                                    Tau(i, 2, 3)*Tau(i, 3, 1)) 
        A(3, 1) =  DetInv*(Tau(i, 2, 1)*Tau(i, 3, 2) - &
                                    Tau(i, 2, 2)*Tau(i, 3, 1)) 
        A(1, 2) = -DetInv*(Tau(i, 1, 2)*Tau(i, 3, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 3, 2)) 
        A(2, 2) =  DetInv*(Tau(i, 1, 1)*Tau(i, 3, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 3, 1)) 
        A(3, 2) = -DetInv*(Tau(i, 1, 1)*Tau(i, 3, 2) - &
                                    Tau(i, 1, 2)*Tau(i, 3, 1)) 
        A(1, 3) =  DetInv*(Tau(i, 1, 2)*Tau(i, 2, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 2, 2)) 
        A(2, 3) = -DetInv*(Tau(i, 1, 1)*Tau(i, 2, 3) - &
                                    Tau(i, 1, 3)*Tau(i, 2, 1)) 
        A(3, 3) =  DetInv*(Tau(i, 1, 1)*Tau(i, 2, 2) - &
                                       Tau(i, 1, 2)*Tau(i, 2, 1)) 
        EE  = MATMUL(A, B)

        D_F(i, 1) = EE(1, 1) 
        D_F(i, 2) = EE(2, 1) 
        D_F(i, 3) = EE(3, 1) 
    END DO

    END SUBROUTINE SolvePreFac

END MODULE GetSolvePreFac
