MODULE GetNormalizeMean

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetLUSGS
    USE GetJacobians
    USE GetTriDi
    USE GetAreaDerivative
    USE GetDX
    USE GetExactSol
    USE GetiMaxValue
    USE Precision_Def
    USE FFT1D

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: NormalizeMean


    CONTAINS

    SUBROUTINE NormalizeMean(   Q_Mean_And_Pet      ,&
                                Q_Mean              ,&
                                Q_Exact             ,&
                                Mean_FFT1           ,&
                                Steps_Per_Cycle     ,&
                                Scaling_Fac         ,&
                                kDeltaX             ,&
                                nL                  ,&
                                nD                  ,&
                                DS                  ,&
                                nDis                ,&
                                nG                  ,&
                                nBound)
    USE GetGlobalParameter

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::    Q_Mean

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    Q_Exact

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)::  &
                                                                Q_Mean_And_Pet

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: Mean_FFT1

    INTEGER, INTENT(IN)::   nL                  ,&
                            nD                  ,&
                            DS                  ,&
                            Steps_Per_Cycle     ,&
                            nG                  ,&
                            nDis                ,&
                            nBound

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
                        nnI             ,&
                        nC              ,&
                        nC3             ,&
                        numOfCycles     ,&
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
                        xC              ,&
                        Step_Number     ,&
                        nnSt            ,&
                        SubIterPerCyc

    REAL(KIND = rDef) ::    delta_x         ,&
                            delta_x1        ,&
                            delta_x2        ,&
                            delta_x3        ,&
                            delta_t         ,&
                            delta_tau       ,&
                            gm1             ,&
                            omega1          ,&
                            omega2          ,&
                            time            ,&
                            Periodic_L2Fac  ,&
                            L2Fac           ,&
                            NewtonL2Fac     ,&
                            epsi_A          ,&
                            epsi_S          ,&
                            sigma_Inv       ,&
                            nu_physical     ,&
                            nu              ,&
                            xx              ,&
                            rho             ,&
                            u               ,&
                            p               ,&
                            Mach            ,&
                            c               ,&
                            eig1_A          ,&
                            eig2_A          ,&
                            eig3_A          ,&
                            eig1_S          ,&
                            eig2_S          ,&
                            eig3_S          ,&
                            P_Static        ,&
                            rho_Static      ,&
                            FinalTime       ,&
                            Cycle_Time      ,&
                            A1_Inflow       ,&
                            A2_Inflow       ,&
                            A3_Inflow       ,&
                            A1_Outflow      ,&
                            A2_Outflow      ,&
                            A3_Outflow

                    
    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE ::     x                   ,&
                                                        x1                  ,&
                                                        delta_xNew          ,&
                                                        Area1               ,&
                                                        Jac_Curv            ,&
                                                        Inv_Jac_Curv        ,&
                                                        NewtonResidual      ,&
                                                        Newtonl2Res         ,&
                                                        bet                 ,&
                                                        Area                ,&
                                                        dAdZi               ,&
                                                        Source_Fac          ,&
                                                        dXdt                ,&
                                                        dZidX               ,&
                                                        dZidt               ,&
                                                        Vel                 ,&
                                                        Vel_Mean            ,&
                                                        Vel_Pet             ,&
                                                        Den                 ,&
                                                        Den_Mean            ,&
                                                        Den_Pet             ,&
                                                        Pres                ,&
                                                        Pres_Mean           ,&
                                                        Pres_Pet            ,&
                                                        Time_FFT            ,&
                                                        Periodic_Res        ,&
                                                        Periodic_SS

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_n                 ,&
                                                        Q_np1               ,&
                                                        Q_nm1               ,&
                                                        Q_nm2               ,&
                                                        dQdTau              ,&
                                                        RHS                 ,&
                                                        E_GPT               ,&
                                                        E                   ,&
                                                        Q_np1_EGPT          ,&
                                                        Q_EGPT              ,&
                                                        Delta_Q_EGPT        ,&
                                                        Q_np1_k             ,&
                                                        Q_np1_kp1           ,&
                                                        Delta_Q_star2       ,&
                                                        Delta_Q_star        ,&
                                                        Delta_Q             ,&
                                                        AA                  ,&
                                                        BB                  ,&
                                                        EE                  ,&
                                                        Dis                 ,&
                                                        Exact               ,&
                                                        Min_Pressure        ,&
                                                        Max_Pressure

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
                                                            
    CHARACTER(LEN = 8)   :: DifferStencil
    CHARACTER(LEN = 4)   :: Dissipation
    CHARACTER(LEN = 5)   :: LeftHandSide
    CHARACTER(LEN = 12)  :: LHS_Diss
    CHARACTER(LEN = 255) :: Current_Directory 
    CHARACTER(LEN = 51):: filename5 ! File 58, Pressure Petrubation on the 
                                    ! start of each cycle AND at every single 
                                    ! step of the final cycle
    
    LOGICAL:: debug

    debug = .TRUE.

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
    nu_physical     = 0.0_rDef

    nnSt            = 4
    nnI             = 0

    ! Where does this come from and Why is it needed? The Inflow/Outflow 
    ! boundary is unsteady now that I'm solving CAA problem. The
    ! Unsteady boundary claims epsilon*SIN(omega*B + omega*time)) 
    ! or epsilon*COS(omega*V +omega*time)). To find the period of a SIN/COS wave
    ! (or as i'm calling it Cycle_Time) we divide 2pi by omega. Note that 
    ! we dont care about 'B' since that only shifts the
    ! functions so no big deal. Nor is epislon important for the period (cycle)
    ! since that deals with the amplitude of the wave. Omega = 0.6*pi 
    
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

    ALLOCATE(   Area(iMax)                                  ,&
                Area1(iMax)                                 ,&
                x1(iMax)                                    ,&
                Source_Fac(iMax)                            ,&
                Q_n(iMax, bMax)                             ,&
                Q_np1(iMax, bMax)                           ,&
                Q_nm1(iMax, bMax)                           ,&
                Q_nm2(iMax, bMax)                           ,&
                dQdTau(iMax, bMax)                          ,&
                dAdZi(iMax)                                 ,&
                RHS(iMax, bMax)                             ,&
                Newtonl2Res(MaxIteration)                   ,&
                NewtonResidual(MaxIteration)                ,&
                E_GPT(iMax, bMax)                           ,&
                E(iMax, bMax)                               ,&
                Dis(iMax, bMax)                             ,&
                Q_np1_k(iMax, bMax)                         ,&
                Mean_FFT1(iMax)                             ,&
                Q_np1_kp1(iMax, bMax)                       ,&
                L_LowerDiag(bMax, bMax, iMax)               ,&
                L_Diag(bMax, bMax, iMax)                    ,&
                U_UpperDiag(bMax, bMax, iMax)               ,&
                U_Diag(bMax, bMax, iMax)                    ,&
                D_Diag(bMax, bMax, iMax)                    ,&
                Exact(iMax, bMax)                           ,&
                Delta_Q_star2(iMax, bMax)                   ,&
                Delta_Q_star(iMax, bMax)                    ,&
                Delta_Q(iMax, bMax)                         ,&
                AA(bMax, bMax)                              ,&
                BB(bMax, bMin)                              ,&
                EE(bMax, bMin)                              ,&
                A_Jac(bMax, bMax, iMax)                     ,&
                A_Plus(bMax, bMax, iMax)                    ,&
                A_Minus(bMax, bMax, iMax)                   ,&
                S_Jac(bMax, bMax, iMax)                     ,&
                S_Plus(bMax, bMax, iMax)                    ,&
                S_Minus(bMax, bMax, iMax)                   ,&
                Jac_Curv(iMax)                              ,&
                Inv_Jac_Curv(iMax)                          ,&
                Vel(iMax)                                   ,&
                Vel_Mean(iMax)                              ,&
                Vel_Pet(iMax)                               ,&
                Den(iMax)                                   ,&
                Den_Mean(iMax)                              ,&
                Den_Pet(iMax)                               ,&
                Pres(iMax)                                  ,&
                Pres_Mean(iMax)                             ,&
                Pres_Pet(iMax)                              ,&
                dXdt(iMax)                                  ,&
                dZidX(iMax)                                 ,&
                dZidt(iMax))

    INCLUDE 'Initialize_Variables.f90'
    INCLUDE 'Inc_CharacterNames.f90'

!! Define Area and Domain
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

    DO i = iMin, iMax
        Source_Fac(i)   = -Jac_Curv(i)*(dAdZi(i))
        Mean_FFT1(i)    = 0.0_rDef
    END DO

    nC              = 0
    Step_Number     = 0

    numOfCycles     = 60
    numTimeSteps    = INT(Steps_Per_Cycle*INT(numOfCycles))

    ALLOCATE(   Q_Store(numTimeSteps, iMax, bMax)               ,&
                Time_FFT(numTimeSteps)                          ,&
                Min_Pressure(Steps_Per_Cycle, iMax)             ,&
                Max_Pressure(Steps_Per_Cycle, iMax)             ,&
                Q_Mean_And_Pet(iMax, bMax)                      ,&
                Periodic_Res(numOfCycles)                       ,&
                Periodic_SS(numOfCycles))

!! Q_Mean_And_Pet**AT THIS POINT**is Q_mean. We use Q_Mean as an 
!! inital condition. 
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_n(i, j) = Q_Mean(i, j)
        END DO
    END DO
    time    = 0.0_rDef
    delta_t = 0.0_rDef

!! Needed For Outflow Boundary Calculation
    P_Static    = p_Exit*Inv_Jac_Curv(1)
    rho_Static  = Q_n(1, 1)

    nI2 = 0

    DO nT = 1, numTimeSteps

        nStages = 4
        ALLOCATE(bet(nStages))
        bet(1)  = (1.0_rDef)/(4.0_rDef) 
        bet(2)  = (1.0_rDef)/(3.0_rDef) 
        bet(3)  = (1.0_rDef)/(2.0_rDef) 
        bet(4)  = 1.0_rDef 
        nu      = 2.0_rDef

        Cycle_Time  = (2.0_rDef)/(0.6_rDef)

        delta_t = Cycle_Time/REAL(Steps_Per_Cycle, rDef)

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

            IF (nL == 1) THEN
                ! Call LUSGS_Mean_Pet
                CALL LUSGS( iMin            = iMin              ,&
                            iMax            = iMax              ,&
                            bMin            = bMin              ,&
                            bMax            = bMax              ,&
                            Q_n             = Q_n               ,&
                            Q_np1_k         = Q_np1_k           ,&
                            Q_np1           = Q_np1             ,&
                            Q_nm1           = Q_nm1             ,&
                            Q_nm2           = Q_nm2             ,&
                            Q_np1_kp1       = Q_np1_kp1         ,&
                            RHS             = RHS               ,&
                            time            = time              ,&
                            dQdTau          = dQdTau            ,&
                            delta_x         = delta_Zi           ,&
                            delta_tau       = delta_tau         ,&
                            delta_t         = delta_t           ,&
                            DS              = DS                ,&
                            nD              = nD                ,&
                            Dis             = Dis               ,&
                            L_LowerDiag     = L_LowerDiag       ,&
                            L_Diag          = L_Diag            ,&
                            U_UpperDiag     = U_UpperDiag       ,&
                            U_Diag          = U_Diag            ,&
                            D_Diag          = D_Diag            ,&
                            A_Plus          = A_Plus            ,&
                            A_Minus         = A_Minus           ,&
                            S_Plus          = S_Plus            ,&
                            S_Minus         = S_Minus           ,&
                            S_Jac           = S_Jac             ,&
                            Delta_Q_star2   = Delta_Q_star2     ,&
                            AA              = AA                ,&
                            BB              = BB                ,&
                            EE              = EE                ,&
                            Delta_Q_star    = Delta_Q_star      ,&
                            Delta_Q         = Delta_Q           ,&
                            bet             = bet               ,&
                            nTau            = nTau              ,&
                            nStages         = nStages           ,&
                            nI              = nI                ,&
                            Newtonl2Res     = Newtonl2Res(nI)   ,&
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
                            Area            = Area              ,&
                            A1_Inflow       = A1_Inflow         ,&
                            A2_Inflow       = A2_Inflow         ,&
                            A3_Inflow       = A3_Inflow         ,&
                            A1_Outflow      = A1_Outflow        ,&
                            A2_Outflow      = A2_Outflow        ,&
                            A3_Outflow      = A3_Outflow        ,&
                            nDis            = nDis              ,&
                            nBound          = nBound)
            ELSE IF (nL == 2) THEN 
                ! Call TriDi_Mean_Pet
                CALL TriDi(  iMin           = iMin              ,&
                             iMax           = iMax              ,&
                             bMin           = bMin              ,&
                             bMax           = bMax              ,&
                             Delta_Q        = Delta_Q           ,&
                             A_Jac          = A_Jac             ,&
                             S_Jac          = S_Jac             ,&
                             Q_n            = Q_n               ,&
                             Q_np1_k        = Q_np1_k           ,&
                             Q_np1_kp1      = Q_np1_kp1         ,&
                             Q_np1          = Q_np1             ,&
                             Q_nm1          = Q_nm1             ,&
                             Q_nm2          = Q_nm2             ,&
                             RHS            = RHS               ,&
                             time           = time              ,&
                             delta_x        = delta_Zi           ,&
                             delta_tau      = delta_tau         ,&
                             delta_t        = delta_t           ,&
                             DS             = DS                ,&
                             nD             = nD                ,&
                             Dis            = Dis               ,&
                             bet            = bet               ,&
                             nTau           = nTau              ,&
                             nStages        = nStages           ,&
                             Scaling_Fac    = Scaling_Fac       ,&
                             E              = E                 ,&
                             E_GPT          = E_GPT             ,&
                             Newtonl2Res    = Newtonl2Res(nI)   ,&
                             dQdTau         = dQdTau            ,&
                             Source_Fac     = Source_Fac        ,&
                             Mach           = Mach              ,&
                             u              = u                 ,&
                             p              = p                 ,&
                             c              = c                 ,&
                             gm1            = gm1               ,&
                             rho            = rho               ,&
                             rho_Static     = rho_Static        ,&
                             P_Static       = P_Static          ,&
                             Q_Exact        = Q_Exact           ,&
                             Jac_Curv       = Jac_Curv          ,&
                             Inv_Jac_Curv   = Inv_Jac_Curv      ,&
                             dZidX          = dZidX             ,&
                             dZidt          = dZidt             ,&
                             nDis           = nDis              ,&
                             Area           = Area              ,&
                             A1_Inflow      = A1_Inflow         ,&
                             A2_Inflow      = A2_Inflow         ,&
                             A3_Inflow      = A3_Inflow         ,&
                             A1_Outflow     = A1_Outflow        ,&
                             A2_Outflow     = A2_Outflow        ,&
                             A3_Outflow     = A3_Outflow        ,&
                             nBound         = nBound)
            ELSE
                CYCLE
            END IF
            
        END DO  !nI

        DO i = iMin, iMax
            DO j = bMin, bMax
                Q_n(i, j)     = Q_np1_kp1(i, j)
                Q_np1_k(i, j) = Q_np1_kp1(i, j)
            END DO
        END DO

        IF (MOD(nT, INT(Steps_Per_Cycle)) == 1) THEN
            Step_Number = Step_Number+1 
        END IF

        time            = time     +   delta_t
        Time_FFT(nT)    = time

        DO i = iMin, iMax
        !!!!!!!!!!!!!!!!!! Pressure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Pres(i)             = ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                    - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                    /Q_n(i, 1)))*Jac_Curv(i)
            Pres_Mean(i)        = ((gam-1.0_rDef)*(Q_Mean(i, 3)     &
                                    - 0.5_rDef*Q_Mean(i, 2)*Q_Mean(i, 2)&
                                    /Q_Mean(i, 1)))*Jac_Curv(i)
            Pres_Pet(i)         = (Pres(i) - Pres_Mean(i))
        END DO

        !!! Write Out Solution on the Final Cycle at every Time Step
        IF (MOD(nT, 5) == 1) THEN
            WRITE(filename5, '(a5, a1, a12, a1, a8, a1, a4, a1, i8.8, a2, a8)')&
                    LeftHandSide,  '_',  LHS_Diss,'_', DifferStencil, '_',           &
                    Dissipation, '_', nT, 'nT', 'Petu.dat' 
            OPEN(58, FILE = TRIM(Current_Directory)//&
                        &'/TimeStepSol/Normalize/'&
                        &//filename5, FORM = 'FORMATTED')
            DO i = iMin,   iMax
                WRITE(58, *) x(i)                               ,&
                            !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                            Den_Pet(i)                          ,&
                            !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                            Vel_Pet(i)                          ,&
                            !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                            Pres_Pet(i)                         ,&
                            !!!!!!!!!!!!Exact Density!!!!!!!!!
                            time                                ,&
                            Step_Number
            END DO
            CLOSE(58)
        END IF

        !!!!!!! Calculate The Periodic Steady State !!!!!!!!!!!!!!!!!!!
        IF (MOD(nT, INT(Steps_Per_Cycle)) == 1) THEN
            nnI                 = Step_Number 
            Periodic_L2Fac = 0.0_rDef
            DO i = iMin, iMax
                Periodic_L2Fac = Periodic_L2Fac+Pres_Pet(i)*Pres_Pet(i)
            END DO
            Periodic_L2Fac = SQRT(Periodic_L2Fac)/REAL(i, rDef)

            Periodic_Res(nnI)   = Periodic_L2Fac
            IF (nnI == 1) THEN
                Periodic_SS(nnI)    = ABS(Periodic_Res(nnI))
            ELSE
                Periodic_SS(nnI)    = ABS(  Periodic_Res(nnI)       -   &
                                            Periodic_Res(nnI-1))
            END IF
            IF ((Periodic_SS(nnI) < tol).OR.(nnI == numOfCycles)) THEN
                numOfCycles     = nnI
                numTimeSteps    = INT(Steps_Per_Cycle*INT(numOfCycles))
            END IF
        END IF
        !!!!!!! End Calculate The Periodic Steady State !!!!!!!!!!!!!!!!!!!

        IF (MOD(nT, 1) == 0) THEN
            WRITE(6, *) 'Normalize Flow Calculation '   ,&
                        LeftHandSide, ' '               ,&
                        LHS_Diss, ' '           ,&
                        DifferStencil, ' '              ,&
                        Dissipation                     ,&
                        nT                              ,&
                        nI                              ,&
                        time                            ,&
                        Step_Number                     ,&
                        nu_physical                     ,&
                        Steps_Per_Cycle                 ,&
                        Periodic_SS(nnI)
        END IF

        !!!!!!! FFT Calculation Starts Here !!!!!!!!!!!!!!!!!
        IF (MOD(nT, Steps_Per_Cycle) == 1) THEN
            SubIterPerCyc   = nI
            nC              = 1
        ELSE
            SubIterPerCyc   = SubIterPerCyc+nI
            nC = nC+1
        END IF

        DO i = iMin, iMax
            Max_Pressure(nC, i) = Pres_Pet(i)
            Min_Pressure(nC, i) = Pres_Pet(i)
        END DO

        DEALLOCATE(bet)

        IF (nT >= numOfCycles*Steps_Per_Cycle) THEN
            WRITE(6, *) "Periodic Steady-State Reached"
            DO i = iMin, iMax
                Mean_FFT1(i) = (MAXVAL(Max_Pressure(:, i))  +          & 
                                MINVAL(Min_Pressure(:, i)))/2.0_rDef
            END DO
            EXIT
        END IF

    END DO ! nT

    DEALLOCATE( Q_Store             ,&
                Time_FFT            ,&
                Min_Pressure        ,&
                Max_Pressure        ,&
                Periodic_Res        ,&
                Periodic_SS)

    DO i = iMin, iMax
        DO j = bMin, bMax
            Q_Mean_And_Pet(i, j)  = Q_n(i, j)
        END DO
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END SUBROUTINE NormalizeMean

END MODULE GetNormalizeMean
