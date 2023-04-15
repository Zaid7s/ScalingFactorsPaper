MODULE GetMeanValues
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetLUSGS
    USE GetJacobians
    USE GetTriDi
    USE GetAreaDerivative
    USE GetDX
    USE GetCFL_Ramping
    USE GetExactSol
    USE GetiMaxValue

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: MeanValues

    TYPE(CFLRampObject):: thisCFLObject

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE MeanValues(  Q_Mean      ,&
                            Q_Exact     ,&
                            nFlow       ,&
                            p_Exit      ,&
                            nL          ,&
                            nD          ,&
                            DS          ,&
                            nG)

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT) ::        &
                                                        Q_Mean              ,&
                                                        Q_Exact

    INTEGER, INTENT(IN)::   nFlow   ,&
                            DS      ,&
                            nG

    INTEGER, INTENT(OUT)::  nL  ,&
                            nD

    REAL(KIND = rDef), INTENT(INOUT):: p_Exit

    INTEGER     ::      xMin            ,&
                        xMax            ,&
                        iMax3           ,&
                        iMax4           ,&
                        iMin            ,&
                        iMax_Original   ,&
                        iMax            ,&
                        nPts            ,&
                        bMin            ,&
                        bMax            ,&
                        i               ,&
                        ii              ,&
                        nT              ,&
                        nI              ,&
                        nI2             ,&
                        numTimeSteps    ,&
                        MaxIteration    ,&
                        k               ,&
                        j               ,&
                        l               ,&
                        m               ,&
                        n               ,&
                        nTau            ,&
                        nStages         ,&
                        nLHS            ,&
                        iQuarter        ,&
                        nDis            ,&
                        xC              ,&
                        nScal

    REAL(KIND = rDef) ::    delta_x     ,&
                            delta_x1    ,&
                            delta_x2    ,&
                            delta_x3    ,&
                            pi          ,&
                            delta_t     ,&
                            delta_Zi    ,&
                            delta_tau   ,&
                            tol         ,&
                            gam         ,&
                            gm1         ,&
                            omega1      ,&
                            omega2      ,&
                            time        ,&
                            L2Fac       ,&
                            NewtonL2Fac, &
                            epsi_A      ,&
                            epsi_S      ,&
                            Scaling_Fac, &
                            kDeltaX     ,&
                            sigma_Inv   ,&
                            nu_physical, &
                            nu          ,&
                            xx          ,&
                            rho         ,&
                            u           ,&
                            p           ,&
                            Mach        ,&
                            c           ,&
                            eig1_A      ,&
                            eig2_A      ,&
                            eig3_A      ,&
                            eig1_S      ,&
                            eig2_S      ,&
                            eig3_S      ,&
                            start_time  ,&
                            end_time    ,&
                            total_time  ,&
                            work        ,&
                            T0In        ,&
                            p0In        ,&
                            rho0In      ,&
                            P_Static    ,&
                            rho_Static  ,&
                            Factor2     ,&
                            Higher_ScaleFac         ,&
                            A1_Inflow   ,&
                            A2_Inflow   ,&
                            A3_Inflow   ,&
                            A1_Outflow  ,&
                            A2_Outflow  ,&
                            A3_Outflow

                    
    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE ::     x               ,&
                                                        x1              ,&
                                                        delta_xNew      ,&
                                                        Area1           ,&
                                                        Jac_Curv        ,&
                                                        Inv_Jac_Curv    ,&
                                                        l2Res           ,&
                                                        Residual        ,&
                                                        MaxWaveSpeed    ,&
                                                        NewtonResidual  ,&
                                                        Newtonl2Res     ,&
                                                        bet             ,&
                                                        Area            ,&
                                                        dAdZi            ,&
                                                        Source_Fac      ,&
                                                        dXdt            ,&
                                                        dZidX           ,&
                                                        dZidt

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_peturbation   ,&
                                                        Q_n             ,&
                                                        Q_np1           ,&
                                                        Q_nm1           ,&
                                                        Q_nm2           ,&
                                                        dQdTau          ,&
                                                        RHS             ,&
                                                        E_GPT           ,&
                                                        E               ,&
                                                        Q_np1_EGPT      ,&
                                                        Q_EGPT          ,&
                                                        Delta_Q_EGPT    ,&
                                                        Q_np1_k         ,&
                                                        Q_np1_kp1       ,&
                                                        Delta_Q_star2   ,&
                                                        Delta_Q_star    ,&
                                                        Delta_Q         ,&
                                                        AA              ,&
                                                        BB              ,&
                                                        EE              ,&
                                                        Dis             ,&
                                                        Exact

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   A_Jac           ,&
                                                            S_Jac           ,&
                                                            L_LowerDiag     ,&
                                                            L_Diag          ,&
                                                            U_UpperDiag     ,&
                                                            U_Diag          ,&
                                                            D_Diag          ,&
                                                            A_Plus          ,&
                                                            A_Minus         ,&
                                                            S_Plus          ,&
                                                            S_Minus         ,&
                                                            Q_Store

                                                            
    CHARACTER(LEN = 8)  :: DifferStencil
    CHARACTER(LEN = 4)  :: Dissipation
    CHARACTER(LEN = 5)  :: LeftHandSide
    CHARACTER(LEN = 12)  :: LHS_Diss

    CHARACTER(LEN = 27):: filename1
    CHARACTER(LEN = 48):: filename2
    CHARACTER(LEN = 27):: filename4
    CHARACTER(LEN = 12):: filename5
    CHARACTER(LEN = 18):: filename6
    CHARACTER(LEN = 28):: filename7
    CHARACTER(LEN = 52):: filename8
    CHARACTER(LEN = 255) :: Current_Directory 

    CALL GET_ENVIRONMENT_VARIABLE('PWD',Current_Directory)
    
    pi              = 4.0_rDef*ATAN(1.0_rDef)
    xMin            = -10                      
    xMax            =  10                    
    delta_Zi        = 1.0_rDef     
    iMin            = 1
    bMin            = 1
    bMax            = 3
    tol             = 5E-13_rDef
    gam             = 1.40_rDef
    gm1             = 1.0_rDef/(gam-1.0_rDef)
    numTimeSteps    = 2000000
    MaxIteration    = 1 
    omega1          = pi/10.0_rDef
    omega2          = pi/50.0_rDef
    time            = 0.0_rDef
    nStages         = 4
    delta_t         = 0.0_rDef

    p_Exit          = 1.0_rDef/gam  ! 0.650_rDef
    p0In            = 0.7975371184_rDef 
    rho0In          = 1.081930199_rDef 
    t0In            = 1.032_rDef 

    nT              = 1

    !call get_environment_variable('PWD',cwd)
    !print *, "The current working directory is: ",trim(cwd)
    
    CALL iMaxValue( xMin        = xMin          ,&
                    xMax        = xMax          ,&
                    x0          = xC            ,&
                    iMin        = iMin          ,&
                    iMax        = iMax          ,&
                    dX_Left     = delta_x1      ,&
                    dX_Mid      = delta_x2      ,&
                    dX_Right    = delta_x3      ,&
                    xNew        = x             ,&
                    deltaX_New  = delta_xNew    ,&
                    nG          = nG)

    nPts            = iMax

    ALLOCATE(   Area(iMax)                          ,&
                Area1(iMax)                         ,&
                x1(iMax)                            ,&
                Source_Fac(iMax)                    ,&
                Q_peturbation(iMax, bMax)           ,&
                Q_Mean(iMax, bMax)                  ,&
                Q_n(iMax, bMax)                     ,&
                Q_Store(numTimeSteps, iMax, bMax)   ,&
                Q_np1(iMax, bMax)                   ,&
                Q_nm1(iMax, bMax)                   ,&
                Q_nm2(iMax, bMax)                   ,&
                dQdTau(iMax, bMax)                  ,&
                dAdZi(iMax)                         ,&
                RHS(iMax, bMax)                     ,&
                l2Res(numTimeSteps)                 ,&
                Residual(numTimeSteps)              ,&
                MaxWaveSpeed(numTimeSteps)          ,&
                Newtonl2Res(MaxIteration)           ,&
                NewtonResidual(MaxIteration)        ,&
                E_GPT(iMax, bMax)                   ,&
                E(iMax, bMax)                       ,&
                Dis(iMax, bMax)                     ,&
                Q_np1_k(iMax, bMax)                 ,&
                Q_Exact(iMax, bMax)                 ,&
                Q_np1_kp1(iMax, bMax)               ,&
                L_LowerDiag(bMax, bMax, iMax)       ,&
                L_Diag(bMax, bMax, iMax)            ,&
                U_UpperDiag(bMax, bMax, iMax)       ,&
                U_Diag(bMax, bMax, iMax)            ,&
                D_Diag(bMax, bMax, iMax)            ,&
                Exact(iMax, bMax)                   ,&
                Delta_Q_star2(iMax, bMax)           ,&
                Delta_Q_star(iMax, bMax)            ,&
                Delta_Q(iMax, bMax)                 ,&
                AA(bMax, bMax)                      ,&
                BB(bMax, bMin)                      ,&
                EE(bMax, bMin)                      ,&
                bet(nStages)                        ,&
                A_Jac(bMax, bMax, iMax)             ,&
                A_Plus(bMax, bMax, iMax)            ,&
                A_Minus(bMax, bMax, iMax)           ,&
                S_Jac(bMax, bMax, iMax)             ,&
                S_Plus(bMax, bMax, iMax)            ,&
                S_Minus(bMax, bMax, iMax)           ,&
                Jac_Curv(iMax)                      ,&
                Inv_Jac_Curv(iMax)                  ,&
                dXdt(iMax)                          ,&
                dZidX(iMax)                         ,&
                dZidt(iMax))
    
    DO i = iMin, iMax
        DO j = bMin, bMax 
            DO k = bMin, bMax 
                Area(i)                         = 0.0_rDef
                Area1(i)                        = 0.0_rDef
                x1(i)                           = 0.0_rDef
                Source_Fac(i)                   = 0.0_rDef
                Q_n(i, j)                       = 0.0_rDef
                Q_np1(i, j)                     = 0.0_rDef
                Q_nm1(i, j)                     = 0.0_rDef
                Q_nm2(i, j)                     = 0.0_rDef
                dQdTau(i, j)                    = 0.0_rDef
                dAdZi(i)                        = 0.0_rDef
                RHS(i, j)                       = 0.0_rDef
                E_GPT(i, j)                     = 0.0_rDef
                E(i, j)                         = 0.0_rDef
                Dis(i, j)                       = 0.0_rDef
                Q_np1_k(i, j)                   = 0.0_rDef
                Q_np1_kp1(i, j)                 = 0.0_rDef   
                L_LowerDiag(k, j, i)            = 0.0_rDef
                L_Diag(k, j, i)                 = 0.0_rDef
                U_UpperDiag(k, j, i)            = 0.0_rDef
                U_Diag(k, j, i)                 = 0.0_rDef
                D_Diag(k, j, i)                 = 0.0_rDef
                Exact(i, j)                     = 0.0_rDef   
                Delta_Q_star2(i, j)             = 0.0_rDef   
                Delta_Q_star(i, j)              = 0.0_rDef   
                Delta_Q(i, j)                   = 0.0_rDef
                AA(k, j)                        = 0.0_rDef
                BB(j, 1)                        = 0.0_rDef
                EE(j, 1)                        = 0.0_rDef
                bet(nStages)                    = 0.0_rDef   
                A_Jac(k, j, i)                  = 0.0_rDef
                A_Plus(k, j, i)                 = 0.0_rDef
                A_Minus(k, j, i)                = 0.0_rDef
                S_Jac(k, j, i)                  = 0.0_rDef
                S_Plus(k, j, i)                 = 0.0_rDef
                S_Minus(k, j, i)                = 0.0_rDef
                Jac_Curv(i)                     = 0.0_rDef
                Inv_Jac_Curv(i)                 = 0.0_rDef
                dXdt(i)                         = 0.0_rDef
                dZidX(i)                        = 0.0_rDef
                dZidt(i)                        = 0.0_rDef
            END DO
        END DO
    END DO

    WRITE(filename5, '(a12)') 'Jacobian.dat'
    
    OPEN (64, FILE = filename5, FORM = 'FORMATTED')

    CALL DX(  xMin          = xMin          ,&
              xMax          = xMax          ,&
              iMin          = iMin          ,&
              iMax          = iMax          ,&
              delta_xNew    = delta_xNew    ,&
              xNew          = x             ,&
              x0            = xC            ,&
              dX_Left       = delta_x1      ,&
              dX_Mid        = delta_x2      ,&
              dX_Right      = delta_x3      ,&
              nG            = nG)

    DO i = iMin, iMax
        WRITE(64, *) i, delta_xNew(i)
    END DO

    CLOSE(64)
    
    !! RK Coefficients
    bet(1)              = (1.0_rDef)/(4.0_rDef) 
    bet(2)              = (1.0_rDef)/(3.0_rDef) 
    bet(3)              = (1.0_rDef)/(2.0_rDef) 
    bet(4)              =  1.0_rDef 
    
    DO nScal = 1, 1
        Higher_ScaleFac = REAL(nScal - 1, rDef)*0.1_rDef + 1.0_rDef

        DO nL = 1, 2
        !! Different LHS. 1--> Uses Lower upper Symmetric Gauss Seidel Method by 
        !!                      Yoon and Jameson
        !!                2--> Uses ADI on the LHS
            IF (nL == 1) THEN
                LeftHandSide = "LUSGS"
            ELSE IF (nL == 2) THEN
                LeftHandSide = "TriDi"
            ELSE
                CYCLE
            END IF

            DO nDis = 1, 2
                IF (nDis == 1) THEN            ! First Order Backward Difference
                    LHS_Diss = 'LHS_Diss_Off'
                ELSEIF (nDis == 2) THEN        ! Second Order Backward Difference
                    LHS_Diss = 'LHS_Diss_On_'
                ELSE
                    WRITE(6, *) "Uninitialized Variable, or pick another nDis"
                    STOP
                END IF

                IF ((nScal > 1).AND.(nL == 1)) THEN
                    CYCLE
                END IF
                IF ((nL == 1).AND.(nDis == 2)) THEN
                    CYCLE
                END IF
                IF (((nDis == 2).AND.(nL == 2)).AND.(nScal > 1)) THEN
                    CYCLE
                END IF

                !! Different RHS Dissipation.   1--> Second Order Dissipation
                !!                              2--> Fourth Order
                !!                              3--> Sixth Order
                !!                              4--> 8th Order
                !!                              5--> Tenth Order
                DO nD = 5, 5
                    IF (nD == 1) THEN
                        Dissipation = 'D_02'
                    ELSE IF (nD == 2) THEN
                        Dissipation = 'D_04'
                    ELSE IF (nD == 3) THEN
                        Dissipation = 'D_06'
                    ELSE IF (nD == 4) THEN
                        Dissipation = 'D_08'
                    ELSE IF (nD == 5) THEN
                        Dissipation = 'D_10'
                    ELSE
                        EXIT
                    END IF
                    !! Different RHS Differncing Stencil.   1--> Second Order
                    !!                                      2--> Fourth Order
                    !!                                      3--> Sixth Order
                    !!                                      4--> RDRP
                    IF (DS == 1) THEN
                        Scaling_Fac     = 1.0_rDef 
                        DifferStencil   = 'SecOrder'
                        kDeltaX         = 1.0_rDef
                    ELSE IF (DS == 2) THEN
                        Scaling_Fac     = 2.12179602298_rDef
                        DifferStencil   = 'FouOrder'
                        kDeltaX         = 1.400_rDef
                    ELSE IF (DS == 3) THEN
                        Scaling_Fac     = 6.352796985_rDef
                        DifferStencil   = 'SixOrder'
                        kDeltaX         = 1.586_rDef
                    ELSE IF (DS == 4) THEN
                        Scaling_Fac     = 6.36206870932_rDef
                        DifferStencil   = 'RDRPSten'
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 5) THEN
                        Scaling_Fac     = 4.97487147161_rDef
                        DifferStencil   = '_DRPSten'
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 6) THEN
                        Scaling_Fac     = 3.800180_rDef
                        DifferStencil   = 'Compact4'
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 7) THEN
                        Scaling_Fac     = 3.800180_rDef
                        DifferStencil   = 'PreFac_4'
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 8) THEN
                        Scaling_Fac     = 10.04379_rDef
                        DifferStencil   = 'PreFac_6'
                        kDeltaX         = 1.664_rDef
                    ELSE
                        EXIT
                    END IF
                            
                    Scaling_Fac = Scaling_Fac*Higher_ScaleFac

                    IF (nScal == 7) THEN
                        Scaling_Fac = 6.36206870932_rDef
                    END IF

                    !IF ((DS == 1).AND.(nD < 2)) THEN
                    !    CYCLE
                    !END IF

                    !IF ((DS == 2).AND.(nD < 4)) THEN
                    !    CYCLE
                    !END IF

                    !IF ((DS == 3).AND.(nD < 5)) THEN
                    !    CYCLE
                    !END IF
                    !
                    !IF ((DS == 4).AND.(nD < 5)) THEN
                    !    CYCLE
                    !END IF
                    
                    WRITE(filename2, '(a5, a1, a8, a1, a4, a12, a9, i4.4, a4)') &
                            LeftHandSide, '_', DifferStencil, '_', &
                            Dissipation, LHS_Diss, '_Converge' ,nScal, '.dat'
                    
                    OPEN(56, FILE = TRIM(Current_Directory)//&
                                &'/ROC/'&
                                &//filename2, FORM = 'FORMATTED')

                    !! Define Area and Domain

                    DO i = iMin, iMax
                        IF (x(i) > 0.0_rDef) THEN
                            Area(i) = 0.536572_rDef-0.198086_rDef*&
                                      EXP(-1.0_rDef*LOG(2.0_rDef)*((x(i)/0.6_rDef)&
                                      *(x(i)/0.6_rDef)))
                        ELSEIF (x(i) <= 0.0_rDef) THEN
                            Area(i) = 1.0_rDef-0.661514_rDef*&
                                      EXP(-1.0_rDef*LOG(2.0_rDef)*((x(i)/0.6_rDef)&
                                      *(x(i)/0.6_rDef)))
                        END IF
                    END DO

                    DO i = iMin, iMax
                        dZidX(i)        = 1.0_rDef/delta_xNew(i)
                        Jac_Curv(i)     = dZidX(i) 
                        Inv_Jac_Curv(i) = 1.0_rDef/(Jac_Curv(i))
                        dXdt(i)         = 0.0_rDef 
                        dZidt(i)        = Jac_Curv(i)*dXdt(i)
                    END DO

                    !! Get dAdZi
                    CALL AreaDerivative(    iMin    =   iMin        ,&
                                            iMax    =   iMax        ,&
                                            A       =   Area        ,&
                                            dAdX    =   dAdZi       ,&
                                            delta_x =   delta_Zi    ,&
                                            DS      =   DS)

                    !! Write Out The area to get Exact Solution
                    WRITE(filename1, '(a5, a8, a5, i5.5, a4)') &
                            "Area_", DifferStencil, '_NPTS', nPts, '.dat'
                    OPEN (59, FILE = filename1, FORM = 'FORMATTED')

                    !WRITE(59, *) nPts
                    DO i = iMin, iMax
                        Source_Fac(i)  = -Jac_Curv(i)*(dAdZi(i))
                        WRITE(59, *) i, x(i), Area(i), Jac_Curv(i), delta_xNew(i)
                    END DO
                    CLOSE(59)

                    CALL ExactSol(  iMin    = iMin          ,&
                                    iMax    = iMax          ,&
                                    p0In    = p0In          ,&
                                    t0In    = T0In          ,&
                                    p_Exit  = p_Exit        ,&
                                    tol     = tol           ,&
                                    ii      = i             ,&
                                    x       = x             ,&
                                    A       = Area          ,&
                                    p       = Exact(:, 1)   ,&
                                    rho     = Exact(:, 2)   ,&
                                    u       = Exact(:, 3))

                    !!! NOTE THE EXACT SOLUTION FILE WRITES OUT IN THE FOLLOWING
                    !!! FORMAT : X(i), PRESSURE(i), DENSITY(i), VELOCITY(i). 
                    !!! Why?  WHO THE FUCK KNOWS!!
                    !!! Hence, Q_Exact has a weird format. 
                    DO i = iMin, iMax
                        Q_Exact(i, 1) = Exact(i, 2)
                        Q_Exact(i, 2) = Exact(i, 2)*Exact(i, 3)
                        Q_Exact(i, 3) = Exact(i, 1)/(gam-1.0_rDef) +  &
                                        Exact(i, 2)*Exact(i, 3)*&
                                        Exact(i, 3)*0.5_rDef
                    END DO

                    WRITE(filename6, '(a6, a8, a4)') &
                            "Exact_", DifferStencil, '.dat'
                    OPEN (67, FILE = filename6, FORM = 'FORMATTED')

                    !WRITE(59, *) nPts
                    DO i = iMin, iMax
                        WRITE(67, *) x(i), Q_Exact(i, 1)                        ,&
                                           !!!!!!!!!!!!Exact Velocity!!!!!!!!
                                           Q_Exact(i, 2)/Q_Exact(i, 1)          ,&
                                           !!!!!!!!!!!!Exact Pressure!!!!!!!!
                                           ((gam-1.0_rDef)*(Q_Exact(i, 3) &
                                           - 0.5_rDef*Q_Exact(i, 2)*&
                                           Q_Exact(i, 2)/Q_Exact(i, 1)))        ,&
                                           !!!!!!!!!!!!Exact Mach Number!!!!!
                                           (Q_Exact(i, 2)/Q_Exact(i, 1))/&
                                           (SQRT(gam*(((gam-1.0_rDef)*&
                                           (Q_Exact(i, 3) - 0.5_rDef*&
                                           Q_Exact(i, 2)*Q_Exact(i, 2)&
                                           /Q_Exact(i, 1))))/(Q_Exact(i, 1))))

                    END DO
                    CLOSE(67)

                    !! Initial Condition
                    DO i = iMin, iMax
                        rho         = 1.0_rDef          !Exact(i, 2)
                        u           = 0.2_rDef          !Exact(i, 3)
                        p           = p_Exit            !Exact(i, 1)
                        Q_n(i, 1)   = rho*Inv_Jac_Curv(i)
                        Q_n(i, 2)   = rho*u*Inv_Jac_Curv(i)
                        Q_n(i, 3)   = (p/(gam-1.0_rDef))*Inv_Jac_Curv(i) +    &
                                        0.5_rDef*rho*u*u*Inv_Jac_Curv(i)
                    END DO 

                    DO j = bMin, bMax
                        DO i = iMin, iMax
                            Q_Exact(i, j)   = Q_Exact(i, j)*Inv_Jac_Curv(i)
                        END DO
                    END DO

                    !! Needed For Outflow Boundary Calculation
                    P_Static    = p_Exit*Inv_Jac_Curv(1)
                    rho_Static  = Q_n(1, 1)


                    IF (nDis == 2) THEN
                        Scaling_Fac = 2.0_rDef
                    END IF

                    nI2 = 0
                    total_time  = 0.0_rDef
                    work        = 0.0_rDef
                    
                    DO nT = 1, numTimeSteps

                        CALL CPU_TIME(start_time)

                        INCLUDE 'InitializeQ.f90'

                        DO nI = 1, MaxIteration
                            CALL Jacobians(  iMin           = iMin          ,&
                                             iMax           = iMax          ,&
                                             bMin           = bMin          ,&
                                             bMax           = bMax          ,&
                                             A_Jac          = A_Jac         ,&
                                             u              = u             ,&
                                             p              = p             ,&
                                             c              = c             ,&
                                             rho            = rho           ,&
                                             Mach           = Mach          ,&
                                             Q_n            = Q_np1_k       ,&
                                             gm1            = gm1           ,&
                                             gam            = gam           ,&
                                             eig1_A         = eig1_A        ,&
                                             eig2_A         = eig2_A        ,&
                                             eig3_A         = eig3_A        ,&
                                             epsi_A         = epsi_A        ,&
                                             Source_Fac     = Source_Fac    ,&
                                             A_Plus         = A_Plus        ,&
                                             A_Minus        = A_Minus       ,&
                                             eig1_S         = eig1_S        ,&
                                             eig2_S         = eig2_S        ,&
                                             eig3_S         = eig3_S        ,&
                                             epsi_S         = epsi_S        ,&
                                             S_Plus         = S_Plus        ,&
                                             S_Minus        = S_Minus       ,&
                                             S_Jac          = S_Jac         ,&
                                             Jac_Curv       = Jac_Curv      ,&
                                             Inv_Jac_Curv   = Inv_Jac_Curv  ,&
                                             dZidX          = dZidX         ,&
                                             dZidt          = dZidt         ,&
                                             Area           = Area)          

                            !INCLUDE 'CFL_Calculations.f90'

                            nu_physical = 4.0_rDef
                            sigma_Inv   = epsi_A*kDeltaX
                            delta_tau   = 10000000000.0_rDef 
                            delta_t     = nu_physical/SQRT(sigma_Inv*sigma_Inv)

                            IF (nL == 1) THEN
                                CALL LUSGS( iMin          = iMin            ,&
                                            iMax          = iMax            ,&
                                            bMin          = bMin            ,&
                                            bMax          = bMax            ,&
                                            Q_n           = Q_n             ,&
                                            Q_np1_k       = Q_np1_k         ,&
                                            Q_np1         = Q_np1           ,&
                                            Q_nm1         = Q_nm1           ,&
                                            Q_nm2         = Q_nm2           ,&
                                            Q_np1_kp1     = Q_np1_kp1       ,&
                                            RHS           = RHS             ,&
                                            time          = time            ,&
                                            dQdTau        = dQdTau          ,&
                                            delta_x       = delta_x         ,&
                                            delta_Zi      = delta_Zi        ,&
                                            delta_tau     = delta_tau       ,&
                                            delta_t       = delta_t         ,&
                                            DS            = DS              ,&
                                            nD            = nD              ,&
                                            Dis           = Dis             ,&
                                            L_LowerDiag   = L_LowerDiag     ,&
                                            U_UpperDiag   = U_UpperDiag     ,&
                                            A_Plus        = A_Plus          ,&
                                            A_Minus       = A_Minus         ,&
                                            S_Jac         = S_Jac           ,&
                                            Delta_Q_star2 = Delta_Q_star2   ,&
                                            AA            = AA              ,&
                                            BB            = BB              ,&
                                            EE            = EE              ,&
                                            Delta_Q_star  = Delta_Q_star    ,&
                                            Delta_Q       = Delta_Q         ,&
                                            bet           = bet             ,&
                                            nTau          = nTau            ,&
                                            nStages       = nStages         ,&
                                            nI            = nI              ,&
                                            Newtonl2Res   = Newtonl2Res(nI), &
                                            Source_Fac    = Source_Fac      ,&
                                            gam           = gam             ,&
                                            Mach          = Mach            ,&
                                            u             = u               ,&
                                            p             = p               ,&
                                            c             = c               ,&
                                            gm1           = gm1             ,&
                                            rho           = rho             ,&
                                            rho_Static    = rho_Static      ,&
                                            P_Static      = P_Static        ,&
                                            p0In          = p0In            ,&
                                            rho0In        = rho0In          ,&
                                            Q_Exact       = Q_Exact         ,&
                                            Jac_Curv      = Jac_Curv        ,&
                                            Inv_Jac_Curv  = Inv_Jac_Curv    ,&
                                            dZidX         = dZidX           ,&
                                            dZidt         = dZidt           ,&
                                            nFlow         = nFlow           ,&
                                            Area          = Area            ,&
                                            A1_Inflow     = A1_Inflow       ,&
                                            A2_Inflow     = A2_Inflow       ,&
                                            A3_Inflow     = A3_Inflow       ,&
                                            A1_Outflow    = A1_Outflow      ,&
                                            A2_Outflow    = A2_Outflow      ,&
                                            nDis            = nDis          ,&
                                            A3_Outflow    = A3_Outflow)
                            ELSE IF (nL == 2) THEN 
                                CALL TriDi( iMin            = iMin              ,&
                                            iMax            = iMax              ,&
                                            bMin            = bMin              ,&
                                            bMax            = bMax              ,&
                                            Delta_Q         = Delta_Q           ,&
                                            A_Jac           = A_Jac             ,&
                                            S_Jac           = S_Jac             ,&
                                            Q_n             = Q_n               ,&
                                            Q_np1_k         = Q_np1_k           ,&
                                            Q_np1_kp1       = Q_np1_kp1         ,&
                                            Q_np1           = Q_np1             ,&
                                            Q_nm1           = Q_nm1             ,&
                                            Q_nm2           = Q_nm2             ,&
                                            RHS             = RHS               ,&
                                            time            = time              ,&
                                            delta_x         = delta_x           ,&
                                            delta_Zi        = delta_Zi          ,&
                                            delta_tau       = delta_tau         ,&
                                            delta_t         = delta_t           ,&
                                            DS              = DS                ,&
                                            nD              = nD                ,&
                                            Dis             = Dis               ,&
                                            bet             = bet               ,&
                                            nTau            = nTau              ,&
                                            nStages         = nStages           ,&
                                            Scaling_Fac     = Scaling_Fac       ,&
                                            E               = E                 ,&
                                            E_GPT           = E_GPT             ,&
                                            Newtonl2Res     = Newtonl2Res(nI)   ,&
                                            dQdTau          = dQdTau            ,&
                                            Source_Fac      = Source_Fac        ,&
                                            gam             = gam               ,&
                                            Mach            = Mach              ,&
                                            u               = u                 ,&
                                            p               = p                 ,&
                                            c               = c                 ,&
                                            gm1             = gm1               ,&
                                            rho             = rho               ,&
                                            rho_Static      = rho_Static        ,&
                                            P_Static        = P_Static          ,&
                                            p0In            = p0In              ,&
                                            rho0In          = rho0In            ,&
                                            Q_Exact         = Q_Exact           ,&
                                            Jac_Curv        = Jac_Curv          ,&
                                            Inv_Jac_Curv    = Inv_Jac_Curv      ,&
                                            dZidX           = dZidX             ,&
                                            dZidt           = dZidt             ,&
                                            nFlow           = nFlow             ,&
                                            nDis            = nDis              ,&
                                            Area            = Area              ,&
                                            A1_Inflow       = A1_Inflow         ,&
                                            A2_Inflow       = A2_Inflow         ,&
                                            A3_Inflow       = A3_Inflow         ,&
                                            A1_Outflow      = A1_Outflow        ,&
                                            A2_Outflow      = A2_Outflow        ,&
                                            A3_Outflow      = A3_Outflow)
                            ELSE
                                CYCLE
                            END IF

                            !INCLUDE 'Inc_Unsteady.f90'
                        END DO  !nI

                        DO i = iMin, iMax
                            DO j = bMin, bMax
                                Q_n(i, j)     = Q_np1_kp1(i, j)
                                Q_np1_k(i, j) = Q_np1_kp1(i, j)
                                Q_Mean(i, j)  = Q_n(i, j)
                            END DO
                        END DO
                        
                        CALL CPU_TIME(end_time)
                        total_time = ABS(start_time - end_time)
                        work = work + total_time

                        Residual(nT) = Newtonl2Res(1)

                        WRITE(56, *) nT, Residual(nT), work
                        !WRITE(0, *) DifferStencil, LeftHandSide, nT, Residual(nT), work
                        
                        IF (MOD(nT, 1) == 1) THEN
                            WRITE(6, *) 'Mean ', LeftHandSide, ' ',             &
                                        DifferStencil, ' ',                     &
                                        LHS_Diss, ' ',                          &
                                        Dissipation, nT, nI2, Residual(nT), & 
                                        l2Res(nT)
                        END IF

                        !INCLUDE 'Inc_Unsteady.f90'

                        IF (Residual(nT) < tol) THEN
                            WRITE(6, *) 'Mean '                     ,&
                                        LeftHandSide                ,&
                                        ' '                         ,&
                                        LHS_Diss                    ,&
                                        ' '                         ,&
                                        Dissipation                 ,&
                                        ' '                         ,&
                                        DifferStencil               ,&
                                        "  Converged After: "       ,&
                                        nT                          ,&
                                        "TimeStepsi, Time (s) "     ,&
                                        work

                            !WRITE(56, *)    nT, nI, nI2, Residual(nT), & 
                            !            nu_physical, delta_t, nu_physical, delta_tau                  
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            !WRITE(filename7, '(a15, a8, a5)') 'Initialize_CAA_',&
                            !                                    DifferStencil,  &
                            !                                    '_.dat'
                            !OPEN(68, FILE = filename7, FORM = 'FORMATTED')
                            !WRITE(68, *)    p_Exit   ,&
                            !                iMax    ,&
                            !                iMin    ,&
                            !                bMin    ,&
                            !                bMax    ,&
                            !                gam     ,&
                            !                nD      ,&
                            !                nL
                            !CLOSE(68)
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            !WRITE(filename4, '(a14, a8, a5)') 'Mean_Flow_CAA_',&
                            !                                    DifferStencil,  &
                            !                                    '_.dat'
                            !OPEN(66, FILE = filename4, FORM = 'FORMATTED')
                            !DO i = iMin, iMax
                            !    WRITE(66, *)    Q_n(i, 1)       ,&
                            !                    Q_n(i, 2)       ,&
                            !                    Q_n(i, 3)       ,&
                            !                    Q_Exact(i, 1)   ,&
                            !                    Q_Exact(i, 2)   ,&
                            !                    Q_Exact(i, 3)   ,&
                            !                    x(i)            ,&
                            !                    Jac_Curv(i)
                            !END DO
                            !CLOSE(66)
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            WRITE(filename8, '(a5, a1, a8, a1, a4, a1, i4.4, a12, a8)')       &
                                    LeftHandSide, '_', DifferStencil, '_',      &
                                    Dissipation, '_', nScal, LHS_Diss, '_ROC.dat' 
                            OPEN(57, FILE = TRIM(Current_Directory)//&
                                        &'/TimeStepSol/'&
                                        &//filename8, FORM = 'FORMATTED')
                            DO i = iMin,   iMax
                                WRITE(57, *)    x(i)                            ,&
                                                !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                                                Q_n(i, 1)*Jac_Curv(i)           ,&
                                                !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                                                Q_n(i, 2)/Q_n(i, 1)             ,&
                                                !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                                                ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                                - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                                /Q_n(i, 1)))*Jac_Curv(i)        ,&
                                                !!!!!!!!!!!!Mach Number!!!!!!!!!!!
                                                (Q_n(i, 2)/Q_n(i, 1))/(SQRT(gam*(&
                                                ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                                - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                                /Q_n(i, 1))))/(Q_n(i, 1))))     ,& 
                                                !!!!!!!!!!!!Exact Density!!!!!!!!!
                                                Q_Exact(i, 1)*Jac_Curv(i)       ,&
                                                !!!!!!!!!!!!Exact Velocity!!!!!!!!
                                                Q_Exact(i, 2)/Q_Exact(i, 1)     ,&
                                                !!!!!!!!!!!!Exact Pressure!!!!!!!!
                                                ((gam-1.0_rDef)*(Q_Exact(i, 3) &
                                                - 0.5_rDef*Q_Exact(i, 2)*&
                                                Q_Exact(i, 2)/Q_Exact(i, 1)))*&
                                                Jac_Curv(i)                     ,&
                                                !!!!!!!!!!!!Exact Mach Number!!!!!
                                                (Q_Exact(i, 2)/Q_Exact(i, 1))/&
                                                (SQRT(gam*(((gam-1.0_rDef)*&
                                                (Q_Exact(i, 3) - 0.5_rDef*&
                                                Q_Exact(i, 2)*Q_Exact(i, 2)&
                                                /Q_Exact(i, 1))))/(Q_Exact(i, 1)))), &
                                                i, delta_xNew(i), Jac_Curv(i)
                            END DO
                            CLOSE(57)
                            EXIT
                        ELSE IF (Residual(nT) > 100.0_rDef) THEN
                            WRITE(6, *) 'Mean ', DifferStencil, "Blew Up: ", &
                                        nT, "TimeSteps, Press Any TO Continue"&
                                        ,Residual(nT)
                            EXIT
                        END IF
                    END DO  ! nT
                    CLOSE(56)
                END DO ! nD
            END DO ! nDis
        END DO ! nL
    END DO !nScal
     
    END SUBROUTINE MeanValues

END MODULE GetMeanValues
