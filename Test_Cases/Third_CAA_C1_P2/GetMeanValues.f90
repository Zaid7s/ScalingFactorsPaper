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

    CONTAINS

    SUBROUTINE MeanValues(  Q_Mean              ,&
                            Q_Exact             ,&
                            nL                  ,&
                            nD                  ,&
                            DS                  ,&
                            nDis                ,&
                            kDeltaX             ,&
                            Scaling_Fac         ,&
                            nInitalCondition    ,&
                            nG)
    USE GetGlobalParameter

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT) ::        &
                                                        Q_Mean              ,&
                                                        Q_Exact

    INTEGER, INTENT(IN)::   nL                  ,&
                            nD                  ,&
                            DS                  ,&
                            nDis                ,&
                            nG                  ,&
                            nInitalCondition

    REAL(KIND = rDef), INTENT(IN) ::    kDeltaX    ,&
                                        Scaling_Fac

    INTEGER     ::      iMax3           ,&
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
                        xC

    REAL(KIND = rDef) ::    delta_x     ,&
                            delta_x1    ,&
                            delta_x2    ,&
                            delta_x3    ,&
                            delta_t     ,&
                            delta_tau   ,&
                            gm1         ,&
                            omega1      ,&
                            omega2      ,&
                            time        ,&
                            L2Fac       ,&
                            NewtonL2Fac, &
                            epsi_A      ,&
                            epsi_S      ,&
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
                            P_Static    ,&
                            rho_Static  ,&
                            Factor2     ,&
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
    CHARACTER(LEN = 47):: filename8
    CHARACTER(LEN = 255) :: Current_Directory 

    CALL GET_ENVIRONMENT_VARIABLE('PWD',Current_Directory)
    
    iMin            = 1
    bMin            = 1
    bMax            = 3
    gm1             = 1.0_rDef/(gam-1.0_rDef)
    numTimeSteps    = 2000000
    MaxIteration    = 1
    omega1          = pi/10.0_rDef
    omega2          = pi/50.0_rDef
    time            = 0.0_rDef
    nStages         = 4
    delta_t         = 0.0_rDef

    nT              = 1

    !call get_environment_variable('PWD',cwd)
    !print *, "The current working directory is: ",trim(cwd)
    
    CALL iMaxValue( x0          = xC            ,&
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
    
    INCLUDE 'Inc_CharacterNames.f90'

    WRITE(filename5, '(a12)') 'Jacobian.dat'
    
    OPEN (64, FILE = filename5, FORM = 'FORMATTED')

    CALL DX(  iMin          = iMin          ,&
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
    
    WRITE(filename2, '(a5, a1, a8, a1, a4, a12, a9, a4)') &
            LeftHandSide, '_', DifferStencil, '_', &
            Dissipation, LHS_Diss, '_Converge', '.dat'
    
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
        rho         = Exact(i, 2)*REAL(nInitalCondition, rDef)      + &
                        1.0_rDef*REAL(1 - nInitalCondition, rDef)
        u           = Exact(i, 3)*REAL(nInitalCondition, rDef)      + &
                        0.2_rDef*REAL(1 - nInitalCondition, rDef)
        p           = Exact(i, 1)*REAL(nInitalCondition, rDef)      + &
                        p_Exit*REAL(1 - nInitalCondition, rDef)

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

    nu_physical = 4.0_rDef

    nI2 = 0
    
    DO nT = 1, numTimeSteps

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

            sigma_Inv   = epsi_A*kDeltaX
            delta_tau   = 10000000000.0_rDef 
            delta_t     = nu_physical/SQRT(sigma_Inv&
                                              *sigma_Inv)

            IF (nL == 1) THEN
                ! Call LUSGS_Mean
                CALL LUSGS( iMin            = iMin            ,&
                            iMax            = iMax            ,&
                            bMin            = bMin            ,&
                            bMax            = bMax            ,&
                            Q_n             = Q_n             ,&
                            Q_np1_k         = Q_np1_k         ,&
                            Q_np1           = Q_np1           ,&
                            Q_np1_kp1       = Q_np1_kp1       ,&
                            RHS             = RHS             ,&
                            time            = time            ,&
                            dQdTau          = dQdTau          ,&
                            delta_x         = delta_x         ,&
                            delta_tau       = delta_tau       ,&
                            delta_t         = delta_t         ,&
                            DS              = DS              ,&
                            nD              = nD              ,&
                            Dis             = Dis             ,&
                            L_LowerDiag     = L_LowerDiag     ,&
                            L_Diag          = L_Diag          ,&
                            U_UpperDiag     = U_UpperDiag     ,&
                            U_Diag          = U_Diag          ,&
                            D_Diag          = D_Diag          ,&
                            A_Plus          = A_Plus          ,&
                            A_Minus         = A_Minus         ,&
                            S_Plus          = S_Plus          ,&
                            S_Minus         = S_Minus         ,&
                            S_Jac           = S_Jac           ,&
                            Delta_Q_star2   = Delta_Q_star2   ,&
                            AA              = AA              ,&
                            BB              = BB              ,&
                            EE              = EE              ,&
                            Delta_Q_star    = Delta_Q_star    ,&
                            Delta_Q         = Delta_Q         ,&
                            bet             = bet             ,&
                            nTau            = nTau            ,&
                            nStages         = nStages         ,&
                            nI              = nI              ,&
                            Newtonl2Res     = Newtonl2Res(nI), &
                            Source_Fac      = Source_Fac      ,&
                            Mach            = Mach            ,&
                            u               = u               ,&
                            p               = p               ,&
                            c               = c               ,&
                            gm1             = gm1             ,&
                            rho             = rho             ,&
                            rho_Static      = rho_Static      ,&
                            P_Static        = P_Static        ,&
                            Q_Exact         = Q_Exact         ,&
                            Jac_Curv        = Jac_Curv        ,&
                            Inv_Jac_Curv    = Inv_Jac_Curv    ,&
                            dZidX           = dZidX           ,&
                            dZidt           = dZidt           ,&
                            Area            = Area            ,&
                            nDis            = nDis)
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
                            Mach            = Mach              ,&
                            u               = u                 ,&
                            p               = p                 ,&
                            c               = c                 ,&
                            gm1             = gm1               ,&
                            rho             = rho               ,&
                            rho_Static      = rho_Static        ,&
                            P_Static        = P_Static          ,&
                            Q_Exact         = Q_Exact           ,&
                            Jac_Curv        = Jac_Curv          ,&
                            Inv_Jac_Curv    = Inv_Jac_Curv      ,&
                            dZidX           = dZidX             ,&
                            dZidt           = dZidt             ,&
                            nDis            = nDis              ,&
                            Area            = Area)
            ELSE
                CYCLE
            END IF
            
        END DO  !nI

        DO i = iMin, iMax
            DO j = bMin, bMax
                Q_n(i, j)     = Q_np1_kp1(i, j)
                Q_np1_k(i, j) = Q_np1_kp1(i, j)
                Q_Mean(i, j)  = Q_n(i, j)
            END DO
        END DO
        
        Residual(nT) = Newtonl2Res(1)

        WRITE(56, *) nT, Residual(nT)
        !WRITE(6, *) DifferStencil, LeftHandSide, nT, Residual(nT)
        
        !CALL CreateObject(  object      = thisCFLObject     ,& 
        !                    Residual    = Residual(nT)      ,&
        !                    Fac2        = Factor2)

        !IF (nT == 1) THEN
        !    CALL CFL_Ramping(   object      = thisCFLObject, &
        !                        Factor2     = Factor2)
        !ELSE
        !    CALL CFL_Ramping(   object  = thisCFLObject, &
        !                        CFL     = nu_physical   ,&
        !                        Factor2 = Factor2)
        !END IF


        IF (MOD(nT, 1) == 1) THEN
            WRITE(filename8, '(a5, a1, a8, a1, a4, a1, a12, a8)')       &
                    LeftHandSide, '_', DifferStencil, '_',      &
                    Dissipation, '_', LHS_Diss, '_ROC.dat' 
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
        END IF
        
        IF (MOD(nT, 1) == 1) THEN
            WRITE(6, *) 'Mean ', LeftHandSide, ' ',             &
                        DifferStencil, ' ',                     &
                        LHS_Diss, ' ',                          &
                        Dissipation, nT, nI2, Residual(nT), & 
                        l2Res(nT)
        END IF

        IF (Residual(nT) < tol) THEN
            WRITE(6, *) 'Mean ', LeftHandSide, ' ',  Dissipation, &
                        ' ', DifferStencil, &
                        "  Converged After: ", nT, &
                        "TimeStepsi, And ", nI2, "Computational &
                        Work", LHS_Diss
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
            IF (MOD(nT, 1) == 1) THEN
                WRITE(filename8, '(a5, a1, a8, a1, a4, a1, a12, a8)')       &
                        LeftHandSide, '_', DifferStencil, '_',              &
                        Dissipation, '_', LHS_Diss, '_ROC.dat' 
                OPEN(57, FILE = TRIM(Current_Directory)//&
                            &'/TimeStepSol/'&
                            &//filename8, FORM = 'FORMATTED')
                DO i = iMin,   iMax
                    WRITE(57, *)    x(i),                                   &
                                    !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                                    Q_n(i, 1)*Jac_Curv(i),                  &
                                    !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                                    Q_n(i, 2)/Q_n(i, 1),                    &
                                    !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                                    ((gam-1.0_rDef)*(Q_n(i, 3)              &
                                    - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)          &
                                    /Q_n(i, 1)))*Jac_Curv(i),               &
                                    !!!!!!!!!!!!Mach Number!!!!!!!!!!!
                                    (Q_n(i, 2)/Q_n(i, 1))/(SQRT(gam*(       &
                                    ((gam-1.0_rDef)*(Q_n(i, 3)              &
                                    - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)          &
                                    /Q_n(i, 1))))/(Q_n(i, 1)))),            & 
                                    !!!!!!!!!!!!Exact Density!!!!!!!!!
                                    Q_Exact(i, 1)*Jac_Curv(i),              &
                                    !!!!!!!!!!!!Exact Velocity!!!!!!!!
                                    Q_Exact(i, 2)/Q_Exact(i, 1) ,           &
                                    !!!!!!!!!!!!Exact Pressure!!!!!!!!
                                    ((gam-1.0_rDef)*(Q_Exact(i, 3)          &
                                    - 0.5_rDef*Q_Exact(i, 2)*               &
                                    Q_Exact(i, 2)/Q_Exact(i, 1)))*          &
                                    Jac_Curv(i)                     ,       &
                                    !!!!!!!!!!!!Exact Mach Number!!!!!
                                    (Q_Exact(i, 2)/Q_Exact(i, 1))/          &
                                    (SQRT(gam*(((gam-1.0_rDef)*             &
                                    (Q_Exact(i, 3) - 0.5_rDef*              &
                                    Q_Exact(i, 2)*Q_Exact(i, 2)             &
                                    /Q_Exact(i, 1))))/(Q_Exact(i, 1)))),    &
                                    i, delta_xNew(i), Jac_Curv(i)
                END DO
                CLOSE(57)
            END IF

            EXIT
            DO i = iMin, iMax
                DO j = bMin, bMax
                    Q_Mean(i, j)  = Q_n(i, j)
                END DO
            END DO
        ELSE IF (Residual(nT) > 100.0_rDef) THEN
            WRITE(6, *) 'Mean ', DifferStencil, "Blew Up: ", &
                        nT, "TimeSteps, Press Any TO Continue"&
                        ,Residual(nT)
            EXIT
        END IF
    END DO  ! nT
    CLOSE(56)
     
    END SUBROUTINE MeanValues

END MODULE GetMeanValues
