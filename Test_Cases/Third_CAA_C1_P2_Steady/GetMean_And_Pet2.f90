MODULE GetMean_And_Pet
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
    
    PUBLIC:: Mean_And_Pet

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE Mean_And_Pet(Q_Mean_And_Pet  ,&
                            Q_Mean          ,&
                            Q_Peturbation   ,&
                            Q_Exact         ,&
                            nFlow           ,&
                            p_Exit          ,&
                            nL              ,&
                            nD              ,&
                            DS)

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: Q_Mean

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: Q_Mean_And_Pet      ,&
                                                        Q_Exact             ,&
                                                        Q_Peturbation

    INTEGER, INTENT(IN)::   nFlow   ,&
                            nL      ,&
                            nD      ,&
                            DS

    REAL(KIND = rDef), INTENT(IN):: p_Exit

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
                        nnI             ,&
                        nC              ,&
                        nC2             ,&
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
                        nSPC            ,&
                        nStages         ,&
                        nLHS            ,&
                        iQuarter        ,&
                        xC              ,&
                        Steps_Per_Cycle, &
                        Step_Number     ,&
                        nnSt

    REAL(KIND = rDef) ::    delta_x         ,&
                            delta_x1        ,&
                            delta_x2        ,&
                            delta_x3        ,&
                            pi              ,&
                            delta_t         ,&
                            delta_Zi        ,&
                            delta_tau       ,&
                            tol             ,&
                            tol2            ,&
                            gam             ,&
                            gm1             ,&
                            omega1          ,&
                            omega2          ,&
                            time            ,&
                            Periodic_L2Fac  ,&
                            L2Fac           ,&
                            NewtonL2Fac     ,&
                            epsi_A          ,&
                            epsi_S          ,&
                            Scaling_Fac     ,&
                            kDeltaX         ,&
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
                            T0In            ,&
                            p0In            ,&
                            rho0In          ,&
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
                                                        l2Res               ,&
                                                        Residual            ,&
                                                        MaxWaveSpeed        ,&
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
                                                        Pressure_FFT        ,&
                                                        Pressure_FFT_2P     ,&
                                                        Frequency_FFT       ,&
                                                        Time_FFT            ,&
                                                        R_FFT               ,&
                                                        I_FFT               ,&
                                                        M_FFT               ,&
                                                        Max_Pressure_Exact  ,&
                                                        x_Exact             ,&
                                                        Periodic_Res        ,&
                                                        Periodic_SS

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_n                 ,&
                                                        Q_np1               ,&
                                                        Q_nm1               ,&
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
                                                        Q_Peturb_Out        ,&
                                                        Min_Pressure        ,&
                                                        Max_Pressure        ,&
                                                        Mean_Flow_FFT       ,&
                                                        Mean_Flow_FFT_2P    ,&
                                                        R_FFT_2P            ,&
                                                        I_FFT_2P            ,&
                                                        M_FFT_2P            ,& 
                                                        All_FFT

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

    CHARACTER(LEN = 44):: filename   ! File 74, FFT for Exit Pressure 2 Periods
    CHARACTER(LEN = 40):: filename2  ! File 73, FFT for Exit Pressure 2 Periods
    CHARACTER(LEN = 36):: filename3  ! File 70, FFT Everywhere (All time steps)
    CHARACTER(LEN = 40):: filename4  ! File 71, FFT Everywhere (every Cycle)
    CHARACTER(LEN = 38):: filename5  ! File 58, Petrubation
    CHARACTER(LEN = 28):: filename6  ! File 61, P vs Time
    CHARACTER(LEN = 37):: filename7  ! File 72, FFT for Exit Pressure
    CHARACTER(LEN = 36):: filename8  ! File 75, Pressure at all TimeLevels
    CHARACTER(LEN = 40):: filename9  ! File 76, Max/Min Pressure Envelope
    CHARACTER(LEN = 42):: filename10  ! File 77, Exact Max
    CHARACTER(LEN = 49):: filename11  ! File 78, Periodic Steady State
    CHARACTER(LEN = 39):: filename12  ! File 79, FFT Everywhere (every Cycle)

    pi              = 4.0_rDef*ATAN(1.0_rDef)
    xMin            = -10                 
    xMax            =  10               
    delta_Zi        = 1.0_rDef
    iMin            = 1
    bMin            = 1
    bMax            = 3
    tol             = 1E-12_rDef
    tol2            = 1E-20_rDef
    gam             = 1.40_rDef
    gm1             = 1.0_rDef/(gam-1.0_rDef)
    MaxIteration    = 1000
    omega1          = pi/10.0_rDef
    omega2          = pi/50.0_rDef
    time            = 0.0_rDef
    delta_t         = 0.0_rDef
    kDeltaX         = 0.0_rDef  ! Just Needed to Iitialize
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

    p0In            = 0.7975371184_rDef 
    rho0In          = 1.081930199_rDef 
    t0In            = 1.032_rDef 
    
    CALL iMaxValue( xMin        = xMin      ,&
                    xMax        = xMax      ,&
                    x0          = xC        ,&
                    iMin        = iMin      ,&
                    iMax        = iMax      ,&
                    dX_Left     = delta_x1  ,&
                    dX_Mid      = delta_x2  ,&
                    dX_Right    = delta_x3  ,&
                    xNew        = x         ,&
                    deltaX_New  = delta_xNew)

    nPts            = iMax

    ALLOCATE(   Area(iMax)                                  ,&
                Area1(iMax)                                 ,&
                x1(iMax)                                    ,&
                Source_Fac(iMax)                            ,&
                Q_n(iMax, bMax)                             ,&
                Q_np1(iMax, bMax)                           ,&
                Q_nm1(iMax, bMax)                           ,&
                dQdTau(iMax, bMax)                          ,&
                dAdZi(iMax)                                 ,&
                RHS(iMax, bMax)                             ,&
                Newtonl2Res(MaxIteration)                   ,&
                NewtonResidual(MaxIteration)                ,&
                E_GPT(iMax, bMax)                           ,&
                E(iMax, bMax)                               ,&
                Dis(iMax, bMax)                             ,&
                Q_np1_k(iMax, bMax)                         ,&
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

!! Different LHS. 1--> Uses Lower upper Symmetric Gauss Seidel Method by 
!!                      Yoon and Jameson
!!                2--> Uses ADI on the LHS
        IF (nL == 1) THEN
            LeftHandSide = "LUSGS"
        ELSE IF (nL == 2) THEN
            LeftHandSide = "TriDi"
        END IF
!! Different RHS Dissipation.   1--> Second Order Dissipation
!!                              2--> Fourth Order
!!                              3--> Sixth Order
!!                              4--> 8th Order
!!                                  5--> Tenth Order
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
        END IF
!! Different RHS Differncing Stencil.   1--> Second Order
!!                                      2--> Fourth Order
!!                                      3--> Sixth Order
!!                                      4--> RDRP
        IF (DS == 1) THEN
            Scaling_Fac     = 1.8_rDef 
            DifferStencil   = 'SecOrder'
            kDeltaX         = 1.0_rDef
        ELSE IF (DS == 2) THEN
            Scaling_Fac     = 3.791411_rDef
            DifferStencil   = 'FouOrder'
            kDeltaX         = 1.400_rDef
        ELSE IF (DS == 3) THEN
            Scaling_Fac     = 10.567013_rDef
            DifferStencil   = 'SixOrder'
            kDeltaX         = 1.586_rDef
        ELSE IF (DS == 4) THEN
            Scaling_Fac     = 10.6597009_rDef
            DifferStencil   = 'RDRPSten'
            kDeltaX         = 1.664_rDef
        END IF
        
        CALL DX(  xMin          = xMin          ,&
                  xMax          = xMax          ,&
                  iMin          = iMin          ,&
                  iMax          = iMax          ,&
                  delta_xNew    = delta_xNew    ,&
                  xNew          = x             ,&
                  x0            = xC            ,&
                  dX_Left       = delta_x1      ,&
                  dX_Mid        = delta_x2      ,&
                  dX_Right      = delta_x3)
        
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

!WRITE(59, *) nPts
        DO i = iMin, iMax
            Source_Fac(i)  = -Jac_Curv(i)*(dAdZi(i))
        END DO

        WRITE(filename6, '(a5, a1, a8, a1, a4, a9)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_Time.dat'
        OPEN(61, FILE = filename6, FORM = 'FORMATTED')


!! Define Area and Domain


    DO nSPC = 2, 2
        nC2             = 0
        Steps_Per_Cycle = 2**(4+nSPC)
        numOfCycles     = 200
        numTimeSteps    = INT(Steps_Per_Cycle*INT(numOfCycles))

        ALLOCATE(   Q_Store(numTimeSteps, iMax, bMax)               ,&
                    TIme_FFT(numTimeSteps)                          ,&
                    Frequency_FFT(numTimeSteps)                     ,&
                    Q_Peturb_Out(numTimeSteps, bMax)                ,&
                    l2Res(numTimeSteps)                             ,&
                    Residual(numTimeSteps)                          ,&
                    MaxWaveSpeed(numTimeSteps)                      ,&
                    Pressure_FFT(Steps_Per_Cycle)                   ,&
                    Pressure_FFT_2P(Steps_Per_Cycle*nnSt)           ,&
                    All_FFT(numTimeSteps, iMax)                     ,&
                    Mean_Flow_FFT(Steps_Per_Cycle, iMax)            ,&
                    Mean_Flow_FFT_2P(Steps_Per_Cycle*nnSt, iMax)    ,&
                    Min_Pressure(Steps_Per_Cycle, iMax)             ,&
                    Max_Pressure(Steps_Per_Cycle, iMax)             ,&
                    Max_Pressure_Exact(761)                         ,&
                    x_Exact(761)                                    ,&
                    R_FFT(Steps_Per_Cycle*nnSt)                     ,&
                    I_FFT(Steps_Per_Cycle*nnSt)                     ,&
                    M_FFT(Steps_Per_Cycle)                          ,& 
                    R_FFT_2P(Steps_Per_Cycle*nnSt, iMax)            ,&
                    I_FFT_2P(Steps_Per_Cycle*nnSt, iMax)            ,&
                    M_FFT_2P(Steps_Per_Cycle, iMax)                 ,& 
                    Periodic_Res(numOfCycles)                       ,&
                    Periodic_SS(numOfCycles))

        WRITE(filename8, '(a5, a1, a8, a1, a4, a1, i3.3, a13)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_', Steps_Per_Cycle,  '_Pres_All.dat'

        OPEN(75, FILE = '/Users/zaidhs/Documents/PhD/&
                    &Third_CAA_C1_P1/TimeStepSol/'&
                        //filename8, FORM = 'FORMATTED')

        WRITE(filename9, '(a5, a1, a8, a1, a4, a1, i3.3, a17)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_', Steps_Per_Cycle, '_Max_Min_Pres.dat'

        OPEN(76, FILE = '/Users/zaidhs/Documents/PhD/&
                    &Third_CAA_C1_P1/TimeStepSol/'&
                        //filename9, FORM = 'FORMATTED')

        WRITE(filename10, '(a5, a1, a8, a1, a4, a23)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_Max_Min_Pres_Exact.dat'

        OPEN(77, FILE = '/Users/zaidhs/Documents/PhD/&
                    &Third_CAA_C1_P1/TimeStepSol/'&
                        //filename10, FORM = 'FORMATTED')

        WRITE(filename11, '(a5, a1, a8, a1, a4, a1, i3.3, a26)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_', Steps_Per_Cycle, '_Periodic_Steady_State.dat'

        OPEN(78, FILE = '/Users/zaidhs/Documents/PhD/&
                    &Third_CAA_C1_P1/TimeStepSol/'&
                        //filename11, FORM = 'FORMATTED')

        WRITE(filename3, '(a5, a1, a8, a1, a4,  a17)')&
                LeftHandSide, '_', DifferStencil, '_',      &
                Dissipation, '_FFT_TimeStep.dat'
        OPEN(70, FILE = '/Users/zaidhs/Documents/PhD/&
                    &Third_CAA_C1_P1/TimeStepSol/'&
                        //filename3, FORM = 'FORMATTED')

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
        rho_Static  = Q_n(iMin, 1) 

        
        DO nT = 1, numTimeSteps
            nStages = 4
            ALLOCATE(bet(nStages))
            bet(1)  = (1.0_rDef)/(4.0_rDef) 
            bet(2)  = (1.0_rDef)/(3.0_rDef) 
            bet(3)  = (1.0_rDef)/(2.0_rDef) 
            bet(4)  = 1.0_rDef 
            nu      = 2.5_rDef

            Cycle_Time  = 1.0_rDef*(2.0_rDef*pi)/(0.6_rDef*pi)
            Step_Number = INT(time/Cycle_Time+1.0_rDef)

            delta_t = Cycle_Time/REAL(Steps_Per_Cycle, rDef)

            WRITE(filename5, '(a5, a1, a8, a1, a4, a1, i8.8, a2, a8)')&
                    LeftHandSide, '_', DifferStencil, '_',           &
                    Dissipation, '_', nT, 'nT', 'Petu.dat' 

            DO j = bMin, bMax
                DO i = iMin, iMax
                    Q_Store(nT, i, j) = Q_n(i, j)
                END DO
            END DO

            DO j = bMin, bMax
                DO i = iMin, iMax
                    Q_np1_k(i, j)   = Q_Store(nT, i, j)
                END DO
            END DO

            IF (nT == 1) THEN
                DO j = bMin, bMax
                    DO i = iMin, iMax
                        Q_nm1(i, j)   = Q_Store(nT, i, j)
                    END DO
                END DO
            ELSE
                DO j = bMin, bMax
                    DO i = iMin, iMax
                        Q_nm1(i, j)   = Q_Store(nT-1, i, j)
                    END DO
                END DO
            END IF
            
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

                sigma_Inv       = epsi_A*kDeltaX
                
                IF (MOD(nt, 1) == 0) THEN
                    nu_physical     = delta_t*SQRT(sigma_Inv*sigma_Inv)
                    delta_tau       = (nu)/SQRT(sigma_Inv*sigma_Inv &
                                        + (1.0_rDef/delta_t)&
                                        *(1.0_rDef/delta_t))
                END IF

                IF (nL == 1) THEN
                    CALL LUSGS( iMin          = iMin            ,&
                                iMax          = iMax            ,&
                                bMin          = bMin            ,&
                                bMax          = bMax            ,&
                                Q_n           = Q_n             ,&
                                Q_np1_k       = Q_np1_k         ,&
                                Q_np1         = Q_np1           ,&
                                Q_nm1         = Q_nm1           ,&
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
                                L_Diag        = L_Diag          ,&
                                U_UpperDiag   = U_UpperDiag     ,&
                                U_Diag        = U_Diag          ,&
                                D_Diag        = D_Diag          ,&
                                A_Plus        = A_Plus          ,&
                                A_Minus       = A_Minus         ,&
                                S_Plus        = S_Plus          ,&
                                S_Minus       = S_Minus         ,&
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
                
                IF (nI >= 2) THEN
                    NewtonResidual(nI) = ABS(Newtonl2Res(nI) - &
                                    1.0_rDef*Newtonl2Res(nI-1))
                    ELSE
                    NewtonResidual(nI) = Newtonl2Res(nI)
                END IF

                !WRITE(0, *) LeftHandSide, DifferStencil,        &
                !            nT, nI, NewtonResidual(nI),         &
                !            FINDLOC(dQdTau, MAXVAL(dQdTau)),    &
                !            MAXVAL(dQdTau)
                 
                IF (NewtonResidual(nI) > tol2) THEN
                    Q_np1_k = Q_np1_kp1
                ELSE
                    EXIT
                END IF
            END DO

            time = time     +   delta_t
            Time_FFT(nT)        = time

            DO i = iMin, iMax
                DO j = bMin, bMax 
                    Q_n(i, j) = Q_np1_kp1(i, j)
                END DO
            END DO

            DO i = iMin, iMax
            !!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Vel(i)              = Q_n(i, 2)/Q_n(i, 1) 
                Vel_Mean(i)         = Q_Mean(i, 2)/Q_Mean(i, 1)
                Vel_Pet(i)          = Vel(i) - Vel_Mean(i)
            !!!!!!!!!!!!!!!!!! Density   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Den(i)              = Q_n(i, 1)*Jac_Curv(i)
                Den_Mean(i)         = Q_Mean(i, 1)*Jac_Curv(i)
                Den_Pet(i)          = Den(i) - Den_Mean(i)
            !!!!!!!!!!!!!!!!!! Pressure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Pres(i)             = ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                        - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                        /Q_n(i, 1)))*Jac_Curv(i)
                Pres_Mean(i)        = ((gam-1.0_rDef)*(Q_Mean(i, 3)     &
                                        - 0.5_rDef*Q_Mean(i, 2)*Q_Mean(i, 2)&
                                        /Q_Mean(i, 1)))*Jac_Curv(i)
                Pres_Pet(i)         = Pres(i) - Pres_Mean(i)
            END DO

            IF (nT == 1) THEN
                Frequency_FFT(nT)   = 0.0_rDef
            END IF


            IF (MOD(nT, INT(Steps_Per_Cycle)) == 1) THEN
                nnI                 = INT(time/Cycle_Time+1.0_rDef)

                Periodic_L2Fac = 0.0_rDef
                DO i = iMin, iMax
                    Periodic_L2Fac = Periodic_L2Fac+Pres_Pet(i)*Pres_Pet(i)
                END DO
                Periodic_L2Fac = SQRT(Periodic_L2Fac)/REAL(iMax, rDef)

                Periodic_Res(nnI)   = Periodic_L2Fac
                IF (nnI == 1) THEN
                    Periodic_SS(nnI)    = ABS(Periodic_Res(nnI))
                ELSE
                    Periodic_SS(nnI)    = ABS(  Periodic_Res(nnI)       -   &
                                                Periodic_Res(nnI-1))
                END IF


                IF (Periodic_SS(nnI) < tol) THEN
                    numOfCycles     = nnI
                    numTimeSteps    = INT(Steps_Per_Cycle*INT(numOfCycles))
                    DO i = 1, nnI
                        WRITE(78, *)    i               ,&
                                        Periodic_SS(i)
                    END DO
                END IF
            END IF

            IF (MOD(nT, 1) == 0) THEN
                WRITE(6, *) 'Mean and Petrubation ' ,&
                            LeftHandSide, ' '       ,&
                            DifferStencil, ' '      ,&
                            Dissipation             ,&
                            nT                      ,&
                            nI                      ,&
                            time                    ,&
                            nnI                     ,&
                            nu_physical             ,&
                            Steps_Per_Cycle         ,&
                            Periodic_SS(nnI)
            END IF

            
            IF (nSPC == 2) THEN
                IF (MOD(nT, INT(Steps_Per_Cycle/8)) == 1) THEN
                    WRITE(6, *) "Writing Out Solution Here. Check File"
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    OPEN(58, FILE = '/Users/zaidhs/Documents/PhD/&
                                &Third_CAA_C1_P1/TimeStepSol/'&
                                    //filename5, FORM = 'FORMATTED')
                    DO i = iMin,   iMax
                        WRITE(58, *) x(i)                               ,&
                                    !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                                    Den_Pet(i)                          ,&
                                    !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                                    Vel_Pet(i)                          ,&
                                    !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                                    Pres_Pet(i)                         ,&
                                    !!!!!!!!!!!!Exact Density!!!!!!!!!
                                    time, Step_Number
                    END DO
                    CLOSE(58)
                END IF
                !Solution is written out at every time step 
                WRITE(61, *)    nT                                      ,& 
                                !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                                Pres_Pet(iMin)                          ,&
                                Pres_Pet((iMax-1)/3)                  ,&
                                Pres_Pet((iMax-1)/2)                  ,&
                                Pres_Pet((iMax-1)/4)                  ,&
                                Pres_Pet(iMax)                          ,&
                                time                                    ,&
                                INT(time/Cycle_Time+1.0_rDef)           ,&
                                A1_Inflow                               ,&
                                A2_Inflow                               ,&
                                A3_Inflow                               ,&
                                A1_Outflow                              ,&
                                A2_Outflow                              ,&
                                A3_Outflow                              ,&
                                Periodic_L2Fac
                ! Write Out Solution on the Final Cycle at every Time Step
                IF (nT >= Steps_Per_Cycle*numOfCycles-Steps_Per_Cycle) THEN
                    OPEN(58, FILE = '/Users/zaidhs/Documents/PhD/&
                                &Third_CAA_C1_P1/TimeStepSol/'&
                                    //filename5, FORM = 'FORMATTED')
                    DO i = iMin,   iMax
                        WRITE(58, *) x(i)                               ,&
                                    !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                                    Den_Pet(i)                          ,&
                                    !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                                    Vel_Pet(i)                          ,&
                                    !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                                    Pres_Pet(i)                         ,&
                                    time, Step_Number
                    END DO
                    CLOSE(58)
                END IF
            END IF
            

            !! FFT Calculation
            IF (MOD(nT, Steps_Per_Cycle) == 1) THEN
                nC = 1
            ELSE
                nC = nC+1
            END IF

            
            IF (MOD(nT, Steps_Per_Cycle*nnSt) == 1) THEN
                nC2 = 1
            ELSE
                nC2 = nC2+1
            END IF

            DO i = iMin, iMax
                Max_Pressure(nC, i)         = Pres_Pet(i)
                Min_Pressure(nC, i)         = Pres_Pet(i)
            END DO
            
            Pressure_FFT(nC)        = Pres_Pet(iMax) 
            Pressure_FFT_2P(nC2)    = Pres_Pet(iMax)

            DO i = iMin, iMax
                Mean_Flow_FFT(nC, i)        = Pres_Pet(i)
                Mean_Flow_FFT_2P(nC2, i)    = Pres_Pet(i)
            END DO

            IF (nSPC == 2) THEN
                IF (MOD(nT, Steps_Per_Cycle*nnSt) == 0) THEN
                    WRITE(filename2, '(a5, a1, a8, a1, a4, a6, i3.3, a11)')&
                            LeftHandSide, '_', DifferStencil, '_',      &
                            Dissipation, '_Step_', Step_Number, '_FFT_2P.dat'
                    OPEN(73, FILE = '/Users/zaidhs/Documents/PhD/&
                                &Third_CAA_C1_P1/FFT/'&
                                    //filename2, FORM = 'FORMATTED')

                    CALL Real_FFT(  data    = Pressure_FFT_2P       ,&
                                    n       = Steps_Per_Cycle*nnSt  ,&
                                    iSign   = 1_i_def)
                    DO nC2 = 1, Steps_Per_Cycle*nnSt
                        WRITE(73, *)    nC2                     ,&
                                        Pressure_FFT_2P(nC2)
                        WRITE(6, *)    nC2                     ,&
                                        Pressure_FFT_2P(nC2)
                    END DO
                END IF
                CLOSE(73)

                IF (MOD(nT, Steps_Per_Cycle) == 0) THEN
                    WRITE(filename7, '(a5, a1, a8, a1, a4, a6, i3.3, a8)')&
                            LeftHandSide, '_', DifferStencil, '_',      &
                            Dissipation, '_Step_', Step_Number, '_FFT.dat'
                    OPEN(72, FILE = '/Users/zaidhs/Documents/PhD/&
                                &Third_CAA_C1_P1/FFT/'&
                                    //filename7, FORM = 'FORMATTED')
                    WRITE(6, *) "Input"
                    DO nC = 1, Steps_Per_Cycle
                        WRITE(6, *) nC                  ,&
                                    Time_FFT(nC)        ,&
                                    Pressure_FFT(nC)
                    END DO

                    CALL Real_FFT(  data    = Pressure_FFT          ,&
                                    n       = Steps_Per_Cycle       ,&
                                    iSign   = 1_i_def)
                    DO nC = 1, Steps_Per_Cycle
                        WRITE(72, *)    nC                  ,&
                                        Pressure_FFT(nC)
                    END DO
                    WRITE(6, *) "Output"
                    DO nC = 1, Steps_Per_Cycle
                        WRITE(6, *) nC                  ,&
                                    Time_FFT(nC)        ,&
                                    Pressure_FFT(nC)
                    END DO
                END IF
                CLOSE(72)

                IF (MOD(nT, Steps_Per_Cycle*nnSt) == 0) THEN
                    WRITE(filename, '(a5, a1, a8, a1, a4, a6, i3.3, a15)')&
                            LeftHandSide, '_', DifferStencil, '_',      &
                            Dissipation, '_Step_', Step_Number, '_FFT_2P_ALL.dat'
                    OPEN(74, FILE = '/Users/zaidhs/Documents/PhD/&
                                &Third_CAA_C1_P1/FFT/'&
                                    //filename, FORM = 'FORMATTED')

                    DO i = iMin, iMax
                        CALL Real_FFT(  data    = Mean_Flow_FFT_2P(:, i)    ,&
                                        n       = Steps_Per_Cycle*nnSt         ,&
                                        iSign   = 1_i_def)
                    END DO

                    DO nC2 = 1, Steps_Per_Cycle*nnSt
                        WRITE(74, *) (Mean_Flow_FFT_2P(nC2, i), i = iMin, iMax)
                    END DO
                END IF
                CLOSE(74)

                IF (MOD(nT, Steps_Per_Cycle) == 0) THEN
                    WRITE(filename4, '(a5, a1, a8, a1, a4, a6, i3.3, a12)')&
                            LeftHandSide, '_', DifferStencil, '_',      &
                            Dissipation, '_Step_', Step_Number,         &
                            '_FFT_All.dat'
                    OPEN(71, FILE = '/Users/zaidhs/Documents/PhD/&
                                    &Third_CAA_C1_P1/FFT/'&
                                    //filename4, FORM = 'FORMATTED')
                                
                    DO i = iMin, iMax
                        CALL Real_FFT(  data    = Mean_Flow_FFT(:, i)   ,&
                                        n       = Steps_Per_Cycle       ,&
                                        iSign   = 1_i_def)
                    END DO

                    DO nC = 1, Steps_Per_Cycle
                        WRITE(71, *) (Mean_Flow_FFT(nC, i), i = iMin, iMax)
                    END DO

                    DO nC = 1, Steps_Per_Cycle
                        DO i = iMin, iMax
                            All_FFT(nT, i) = Mean_Flow_FFT(nC, i)
                        END DO
                        WRITE(70, *) (All_FFT(nT, i), i = iMin, iMax)
                    END DO

                END IF
                CLOSE(71)
            END IF
            !! End FFT Calculation for Exit Pressure
            
            ! Write The Presssure at each Time Step (75)
            ! Get Max and Min Pressure Envelope along with Mean Flow
            IF (nT >= Steps_Per_Cycle*numOfCycles-Steps_Per_Cycle) THEN
                !WRITE(75, '(401f17.12)') (Max_Pressure(nC, i), i = iMin, iMax)
                WRITE(75, *) (Max_Pressure(nC, i), i = iMin, iMax)
                IF (nT == numTimeSteps) THEN
                    DO i = iMin, iMax
                        WRITE(76, *)    x(i)                                ,&
                                        !!!!!! Max Pressure Envelopee !!!!!!!!!!!
                                        MAXVAL(Max_Pressure(:, i))          ,&     
                                        !!!!!! Min Pressure Envelopee !!!!!!!!!!!
                                        MINVAL(Min_Pressure(:, i))          ,&     
                                        !!!!!! Mean Pressure Envelopee !!!!!!!!!!!
                                        (MAXVAL(Max_Pressure(:, i))  +       &     
                                        MINVAL(Min_Pressure(:, i)))/2.0_rDef, &
                                        !!!!!! Max Pressure Disturbance !!!!!!!!!!!
                                        MAXVAL(Max_Pressure(:, i))   -       &     
                                        (MAXVAL(Max_Pressure(:, i))  +       &     
                                        MINVAL(Min_Pressure(:, i)))/2.0_rDef
                    END DO
                    OPEN(15, FILE = '/Users/zaidhs/Documents/PhD/&
                                &Third_CAA_C1_P1/Exact/CAT1P1.FIG1.DATA'&
                                , FORM = 'FORMATTED', STATUS = 'OLD')
                        DO i = iMin, 761
                            READ(15, *) x_Exact(i)              ,&
                                        Max_Pressure_Exact(i)
                        END DO
                    CLOSE(15)
                    DO i = iMin, 761
                        WRITE(77, *)    x_Exact(i)              ,&
                                        Max_Pressure_Exact(i)
                    END DO

                    WRITE(filename12, '(a5, a1, a8, a1, a4, a4, i3.3, a12)')&
                            LeftHandSide, '_', DifferStencil, '_',      &
                            Dissipation, '_SPC_', Steps_Per_Cycle,         &
                            '_FFT_All.dat'
                    OPEN(79, FILE = '/Users/zaidhs/Documents/PhD/&
                                    &Third_CAA_C1_P1/FFT/'&
                                    //filename12, FORM = 'FORMATTED')
                                
                    DO i = iMin, iMax
                        CALL Real_FFT(  data    = Mean_Flow_FFT(:, i)   ,&
                                        n       = Steps_Per_Cycle       ,&
                                        iSign   = 1_i_def)
                    END DO

                    DO nC = 1, Steps_Per_Cycle
                        WRITE(79, *) (Mean_Flow_FFT(nC, i), i = iMin, iMax)
                    END DO

                    DO nC = 1, Steps_Per_Cycle
                        DO i = iMin, iMax
                            All_FFT(nT, i) = Mean_Flow_FFT(nC, i)
                        END DO
                        WRITE(70, *) (All_FFT(nT, i), i = iMin, iMax)
                    END DO

                CLOSE(71)
                CLOSE(79)
                END IF
            END IF
            DEALLOCATE(bet)
            IF (nT >= numOfCycles*Steps_Per_Cycle) THEN
                WRITE(6, *) "Periodic Steady-State Reached"
                EXIT
            END IF

        END DO

        CLOSE(61)
        CLOSE(70)
        CLOSE(75)
        CLOSE(76)
        CLOSE(77)
        CLOSE(78)

        DEALLOCATE( Q_Store             ,&
                    TIme_FFT            ,&
                    Frequency_FFT       ,&
                    Q_Peturb_Out        ,&
                    l2Res               ,&
                    Residual            ,&
                    MaxWaveSpeed        ,&
                    Pressure_FFT        ,&
                    Pressure_FFT_2P     ,&
                    Min_Pressure        ,&
                    Max_Pressure        ,&
                    Max_Pressure_Exact  ,&
                    x_Exact             ,&
                    All_FFT             ,&
                    Mean_Flow_FFT       ,&
                    Mean_Flow_FFT_2P    ,&
                    R_FFT               ,&
                    I_FFT               ,&
                    M_FFT               ,&
                    R_FFT_2P            ,&
                    I_FFT_2P            ,&
                    M_FFT_2P            ,&
                    Periodic_Res        ,&
                    Periodic_SS)

    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END SUBROUTINE Mean_And_Pet

END MODULE GetMean_And_Pet
