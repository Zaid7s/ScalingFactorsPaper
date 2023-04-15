MODULE GetSolveLUSGS
!!! Working Test Case
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: SolveLUSGS

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE SolveLUSGS(  L_LowerDiag     ,&    
                            U_UpperDiag     ,& 
                            Tau             ,&
                            iMin            ,& 
                            iMax            ,& 
                            bMin            ,& 
                            bMax            ,& 
                            RHS             ,& 
                            Delta_Q_star    ,&
                            Delta_Q_star2   ,&
                            Delta_Q) 

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    L_LowerDiag ,&
                                                            U_UpperDiag ,&
                                                            Tau

    INTEGER, INTENT(IN) ::  iMin    ,&
                            iMax    ,&
                            bMin    ,&
                            bMax

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS             ,&
                                                            Delta_Q_star2   ,&
                                                            Delta_Q_star    ,&
                                                            Delta_Q

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
    B(1, 1) = RHS(i, 1)
    B(2, 1) = RHS(i, 2)
    B(3, 1) = RHS(i, 3)

    Fac3 = (Tau(1, 1, i)*Tau(2, 2, i)*Tau(3, 3, i) - &
            Tau(1, 1, i)*Tau(2, 3, i)*Tau(3, 2, i) - &
            Tau(1, 2, i)*Tau(2, 1, i)*Tau(3, 3, i) + &
            Tau(1, 2, i)*Tau(2, 3, i)*Tau(3, 1, i) + &
            Tau(1, 3, i)*Tau(2, 1, i)*Tau(3, 2, i) - &
            Tau(1, 3, i)*Tau(2, 2, i)*Tau(3, 1, i))

    DetInv  = 1.0_rDef/Fac3
    A(1, 1) =  DetInv*(Tau(2, 2, i)*Tau(3, 3, i) - Tau(2, 3, i)*Tau(3, 2, i)) 
    A(2, 1) = -DetInv*(Tau(2, 1, i)*Tau(3, 3, i) - Tau(2, 3, i)*Tau(3, 1, i)) 
    A(3, 1) =  DetInv*(Tau(2, 1, i)*Tau(3, 2, i) - Tau(2, 2, i)*Tau(3, 1, i)) 
    A(1, 2) = -DetInv*(Tau(1, 2, i)*Tau(3, 3, i) - Tau(1, 3, i)*Tau(3, 2, i)) 
    A(2, 2) =  DetInv*(Tau(1, 1, i)*Tau(3, 3, i) - Tau(1, 3, i)*Tau(3, 1, i)) 
    A(3, 2) = -DetInv*(Tau(1, 1, i)*Tau(3, 2, i) - Tau(1, 2, i)*Tau(3, 1, i)) 
    A(1, 3) =  DetInv*(Tau(1, 2, i)*Tau(2, 3, i) - Tau(1, 3, i)*Tau(2, 2, i)) 
    A(2, 3) = -DetInv*(Tau(1, 1, i)*Tau(2, 3, i) - Tau(1, 3, i)*Tau(2, 1, i)) 
    A(3, 3) =  DetInv*(Tau(1, 1, i)*Tau(2, 2, i) - Tau(1, 2, i)*Tau(2, 1, i)) 

    EE = MATMUL(A, B)

    Delta_Q_star2(i, 1) = EE(1, 1)
    Delta_Q_star2(i, 2) = EE(2, 1)
    Delta_Q_star2(i, 3) = EE(3, 1)

    !!!!!!! Second Step !!!!!!!
    A(1, 1)     = Tau(1, 1, i)
    A(2, 1)     = Tau(2, 1, i)
    A(3, 1)     = Tau(3, 1, i)
    A(1, 2)     = Tau(1, 2, i)
    A(2, 2)     = Tau(2, 2, i)
    A(3, 2)     = Tau(3, 2, i)
    A(1, 3)     = Tau(1, 3, i)
    A(2, 3)     = Tau(2, 3, i)
    A(3, 3)     = Tau(3, 3, i)

    B(1, 1)     = Delta_Q_star2(i, 1)
    B(2, 1)     = Delta_Q_star2(i, 2)
    B(3, 1)     = Delta_Q_star2(i, 3)

    EE          = MATMUL(A, B)

    Delta_Q_star(i, 1) = EE(1, 1)
    Delta_Q_star(i, 2) = EE(2, 1)
    Delta_Q_star(i, 3) = EE(3, 1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO i = iMin + 1, iMax
        A(1, 1)     = L_LowerDiag(1, 1, i)
        A(2, 1)     = L_LowerDiag(2, 1, i)
        A(3, 1)     = L_LowerDiag(3, 1, i)
        A(1, 2)     = L_LowerDiag(1, 2, i)
        A(2, 2)     = L_LowerDiag(2, 2, i)
        A(3, 2)     = L_LowerDiag(3, 2, i)
        A(1, 3)     = L_LowerDiag(1, 3, i)
        A(2, 3)     = L_LowerDiag(2, 3, i)
        A(3, 3)     = L_LowerDiag(3, 3, i)

        B(1, 1)     = Delta_Q_star2(i - 1, 1)
        B(2, 1)     = Delta_Q_star2(i - 1, 2)
        B(3, 1)     = Delta_Q_star2(i - 1, 3)

        LowerFac    = MATMUL(A, B)

        B(1, 1) = RHS(i, 1) - LowerFac(1, 1)
        B(2, 1) = RHS(i, 2) - LowerFac(2, 1)
        B(3, 1) = RHS(i, 3) - LowerFac(3, 1)

        !! The Inverse of Tau is Here
        Fac3 = (Tau(1, 1, i)*Tau(2, 2, i)*Tau(3, 3, i) - &
                Tau(1, 1, i)*Tau(2, 3, i)*Tau(3, 2, i) - &
                Tau(1, 2, i)*Tau(2, 1, i)*Tau(3, 3, i) + &
                Tau(1, 2, i)*Tau(2, 3, i)*Tau(3, 1, i) + &
                Tau(1, 3, i)*Tau(2, 1, i)*Tau(3, 2, i) - &
                Tau(1, 3, i)*Tau(2, 2, i)*Tau(3, 1, i))

        DetInv  = 1.0_rDef/Fac3
        
        A(1, 1) =  DetInv*(Tau(2, 2, i)*Tau(3, 3, i) - &
                                    Tau(2, 3, i)*Tau(3, 2, i)) 
        A(2, 1) = -DetInv*(Tau(2, 1, i)*Tau(3, 3, i) - &
                                    Tau(2, 3, i)*Tau(3, 1, i)) 
        A(3, 1) =  DetInv*(Tau(2, 1, i)*Tau(3, 2, i) - &
                                    Tau(2, 2, i)*Tau(3, 1, i)) 
        A(1, 2) = -DetInv*(Tau(1, 2, i)*Tau(3, 3, i) - &
                                    Tau(1, 3, i)*Tau(3, 2, i)) 
        A(2, 2) =  DetInv*(Tau(1, 1, i)*Tau(3, 3, i) - &
                                    Tau(1, 3, i)*Tau(3, 1, i)) 
        A(3, 2) = -DetInv*(Tau(1, 1, i)*Tau(3, 2, i) - &
                                    Tau(1, 2, i)*Tau(3, 1, i)) 
        A(1, 3) =  DetInv*(Tau(1, 2, i)*Tau(2, 3, i) - &
                                    Tau(1, 3, i)*Tau(2, 2, i)) 
        A(2, 3) = -DetInv*(Tau(1, 1, i)*Tau(2, 3, i) - &
                                    Tau(1, 3, i)*Tau(2, 1, i)) 
        A(3, 3) =  DetInv*(Tau(1, 1, i)*Tau(2, 2, i) - &
                                       Tau(1, 2, i)*Tau(2, 1, i)) 
        !!

        EE = MATMUL(A, B)

        Delta_Q_star2(i, 1) = EE(1, 1)
        Delta_Q_star2(i, 2) = EE(2, 1)
        Delta_Q_star2(i, 3) = EE(3, 1)
    
        !!!!!!! Second Step !!!!!!!
        A(1, 1)     = Tau(1, 1, i)
        A(2, 1)     = Tau(2, 1, i)
        A(3, 1)     = Tau(3, 1, i)
        A(1, 2)     = Tau(1, 2, i)
        A(2, 2)     = Tau(2, 2, i)
        A(3, 2)     = Tau(3, 2, i)
        A(1, 3)     = Tau(1, 3, i)
        A(2, 3)     = Tau(2, 3, i)
        A(3, 3)     = Tau(3, 3, i)

        B(1, 1)     = Delta_Q_star2(i, 1)
        B(2, 1)     = Delta_Q_star2(i, 2)
        B(3, 1)     = Delta_Q_star2(i, 3)

        EE          = MATMUL(A, B)

        Delta_Q_star(i, 1) = EE(1, 1)
        Delta_Q_star(i, 2) = EE(2, 1)
        Delta_Q_star(i, 3) = EE(3, 1)
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i = iMax
    B(1, 1) = Delta_Q_star(i, 1)
    B(2, 1) = Delta_Q_star(i, 2)
    B(3, 1) = Delta_Q_star(i, 3)

    Fac3 = (Tau(1, 1, i)*Tau(2, 2, i)*Tau(3, 3, i) - &
            Tau(1, 1, i)*Tau(2, 3, i)*Tau(3, 2, i) - &
            Tau(1, 2, i)*Tau(2, 1, i)*Tau(3, 3, i) + &
            Tau(1, 2, i)*Tau(2, 3, i)*Tau(3, 1, i) + &
            Tau(1, 3, i)*Tau(2, 1, i)*Tau(3, 2, i) - &
            Tau(1, 3, i)*Tau(2, 2, i)*Tau(3, 1, i))

    DetInv  = 1.0_rDef/Fac3
    A(1, 1) =  DetInv*(Tau(2, 2, i)*Tau(3, 3, i) - Tau(2, 3, i)*Tau(3, 2, i)) 
    A(2, 1) = -DetInv*(Tau(2, 1, i)*Tau(3, 3, i) - Tau(2, 3, i)*Tau(3, 1, i)) 
    A(3, 1) =  DetInv*(Tau(2, 1, i)*Tau(3, 2, i) - Tau(2, 2, i)*Tau(3, 1, i)) 
    A(1, 2) = -DetInv*(Tau(1, 2, i)*Tau(3, 3, i) - Tau(1, 3, i)*Tau(3, 2, i)) 
    A(2, 2) =  DetInv*(Tau(1, 1, i)*Tau(3, 3, i) - Tau(1, 3, i)*Tau(3, 1, i)) 
    A(3, 2) = -DetInv*(Tau(1, 1, i)*Tau(3, 2, i) - Tau(1, 2, i)*Tau(3, 1, i)) 
    A(1, 3) =  DetInv*(Tau(1, 2, i)*Tau(2, 3, i) - Tau(1, 3, i)*Tau(2, 2, i)) 
    A(2, 3) = -DetInv*(Tau(1, 1, i)*Tau(2, 3, i) - Tau(1, 3, i)*Tau(2, 1, i)) 
    A(3, 3) =  DetInv*(Tau(1, 1, i)*Tau(2, 2, i) - Tau(1, 2, i)*Tau(2, 1, i)) 

    EE      = MATMUL(A, B)

    Delta_Q(i, 1) = EE(1, 1) 
    Delta_Q(i, 2) = EE(2, 1) 
    Delta_Q(i, 3) = EE(3, 1) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Third Step
    DO i = iMax - 1, iMin, -iMin
        A(1, 1)     = U_UpperDiag(1, 1, i)
        A(2, 1)     = U_UpperDiag(2, 1, i)
        A(3, 1)     = U_UpperDiag(3, 1, i)
        A(1, 2)     = U_UpperDiag(1, 2, i)
        A(2, 2)     = U_UpperDiag(2, 2, i)
        A(3, 2)     = U_UpperDiag(3, 2, i)
        A(1, 3)     = U_UpperDiag(1, 3, i)
        A(2, 3)     = U_UpperDiag(2, 3, i)
        A(3, 3)     = U_UpperDiag(3, 3, i)

        B(1, 1)     = Delta_Q(i + 1, 1)
        B(2, 1)     = Delta_Q(i + 1, 2)
        B(3, 1)     = Delta_Q(i + 1, 3)

        LowerFac    = MATMUL(A, B)

        B(1, 1) = Delta_Q_star(i, 1) - LowerFac(1, 1)
        B(2, 1) = Delta_Q_star(i, 2) - LowerFac(2, 1)
        B(3, 1) = Delta_Q_star(i, 3) - LowerFac(3, 1)

        !! The Inverse of Tau is Here
        Fac3 = (Tau(1, 1, i)*Tau(2, 2, i)*Tau(3, 3, i) - &
                Tau(1, 1, i)*Tau(2, 3, i)*Tau(3, 2, i) - &
                Tau(1, 2, i)*Tau(2, 1, i)*Tau(3, 3, i) + &
                Tau(1, 2, i)*Tau(2, 3, i)*Tau(3, 1, i) + &
                Tau(1, 3, i)*Tau(2, 1, i)*Tau(3, 2, i) - &
                Tau(1, 3, i)*Tau(2, 2, i)*Tau(3, 1, i))

        DetInv  = 1.0_rDef/Fac3
        
        A(1, 1) =  DetInv*(Tau(2, 2, i)*Tau(3, 3, i) - &
                                    Tau(2, 3, i)*Tau(3, 2, i)) 
        A(2, 1) = -DetInv*(Tau(2, 1, i)*Tau(3, 3, i) - &
                                    Tau(2, 3, i)*Tau(3, 1, i)) 
        A(3, 1) =  DetInv*(Tau(2, 1, i)*Tau(3, 2, i) - &
                                    Tau(2, 2, i)*Tau(3, 1, i)) 
        A(1, 2) = -DetInv*(Tau(1, 2, i)*Tau(3, 3, i) - &
                                    Tau(1, 3, i)*Tau(3, 2, i)) 
        A(2, 2) =  DetInv*(Tau(1, 1, i)*Tau(3, 3, i) - &
                                    Tau(1, 3, i)*Tau(3, 1, i)) 
        A(3, 2) = -DetInv*(Tau(1, 1, i)*Tau(3, 2, i) - &
                                    Tau(1, 2, i)*Tau(3, 1, i)) 
        A(1, 3) =  DetInv*(Tau(1, 2, i)*Tau(2, 3, i) - &
                                    Tau(1, 3, i)*Tau(2, 2, i)) 
        A(2, 3) = -DetInv*(Tau(1, 1, i)*Tau(2, 3, i) - &
                                    Tau(1, 3, i)*Tau(2, 1, i)) 
        A(3, 3) =  DetInv*(Tau(1, 1, i)*Tau(2, 2, i) - &
                                       Tau(1, 2, i)*Tau(2, 1, i)) 
        EE  = MATMUL(A, B)

        Delta_Q(i, 1) = EE(1, 1) 
        Delta_Q(i, 2) = EE(2, 1) 
        Delta_Q(i, 3) = EE(3, 1) 
    END DO

    END SUBROUTINE SolveLUSGS

END MODULE GetSolveLUSGS
