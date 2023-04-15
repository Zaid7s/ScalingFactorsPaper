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

    SUBROUTINE Mean_And_Pet(Q_Peturbation       ,&
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

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: Q_Mean

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: Q_Exact 

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)::  &
                                                                Q_Peturbation

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Mean_FFT1

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
                        nStages         ,&
                        nLHS            ,&
                        iQuarter        ,&
                        xC              ,&
                        Step_Number     ,&
                        nnSt

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
                                                        Time_FFT            ,&
                                                        R_FFT               ,&
                                                        I_FFT               ,&
                                                        M_FFT               ,&
                                                        Check_A             ,&
                                                        Max_Pressure_Exact  ,&
                                                        x_Exact             ,&
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
    CHARACTER(LEN = 12)  :: LHS_Diss

    CHARACTER(LEN = 53):: filename2 ! File 72, P vs Time (every Cycle Seperate File)
    CHARACTER(LEN = 51):: filename5 ! File 58, Pressure Petrubation on the 
                                    ! start of each cycle AND at every single 
                                    ! step of the final cycle
    CHARACTER(LEN = 46):: filename6  ! File 61, P vs Time (All in 1 File)
    CHARACTER(LEN = 54):: filename9  ! File 76, Max/Min Pressure Envelope
    CHARACTER(LEN = 62):: filename11  ! File 78, Periodic Steady State

    CHARACTER(LEN = 255) :: Current_Directory 
    
    LOGICAL:: debug

    debug = .TRUE.

    CALL GET_ENVIRONMENT_VARIABLE('PWD',Current_Directory)

    iMin            = 1
    bMin            = 1
    bMax            = 3
    gm1             = 1.0_rDef/(gam-1.0_rDef)
    MaxIteration    = 400
    omega1          = pi/10.0_rDef
    omega2          = pi/50.0_rDef
    time            = 0.0_rDef
    nStages         = 4
    delta_t         = 0.0_rDef

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
                Q_Peturbation(iMax, bMax)                   ,&
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

        INCLUDE 'Inc_CharacterNames.f90'
        
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

!WRITE(59, *) nPts
        DO i = iMin, iMax
            Source_Fac(i)  = -Jac_Curv(i)*(dAdZi(i))
        END DO

        WRITE(filename6, '(a5, a1, a12, a1, a8, a1, a4, a9)')&
                LeftHandSide, '_', LHS_Diss,'_', DifferStencil, '_',           &
                Dissipation, '_Time.dat'
        OPEN(61, FILE = filename6, FORM = 'FORMATTED')


!! Define Area and Domain

    nC2             = 0
    nC              = 0
    Step_Number     = 0

    numOfCycles     = 400
    numTimeSteps    = INT(Steps_Per_Cycle*INT(numOfCycles))

    ALLOCATE(   Q_Store(numTimeSteps, iMax, bMax)               ,&
                TIme_FFT(numTimeSteps)                          ,&
                Q_Peturb_Out(numTimeSteps, bMax)                ,&
                l2Res(numTimeSteps)                             ,&
                Residual(numTimeSteps)                          ,&
                MaxWaveSpeed(numTimeSteps)                      ,&
                All_FFT(numTimeSteps, iMax)                     ,&
                Mean_Flow_FFT(Steps_Per_Cycle, iMax)            ,&
                Check_A(Steps_Per_Cycle)                        ,&
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    WRITE(filename11, '(a5, a1, a12, a1, a8, a1, a4, a1, i3.3, a26)')&
            LeftHandSide,  '_', LHS_Diss, '_', DifferStencil, '_',           &
            Dissipation, '_', Steps_Per_Cycle, '_Periodic_Steady_State.dat'
    OPEN(78, FILE = TRIM(Current_Directory)//&
                &'/ROC/'&
                &//filename11, FORM = 'FORMATTED')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

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

    nI2 = 0
    
    DO nT = 1, numTimeSteps
        nStages = 4
        ALLOCATE(bet(nStages))
        bet(1)  = (1.0_rDef)/(4.0_rDef) 
        bet(2)  = (1.0_rDef)/(3.0_rDef) 
        bet(3)  = (1.0_rDef)/(2.0_rDef) 
        bet(4)  = 1.0_rDef 
        nu      = 2.0_rDef

        Cycle_Time  = (2.0_rDef*pi)/(0.6_rDef*pi)

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

            sigma_Inv       = epsi_A*kDeltaX
            
            IF (MOD(nt, 1) == 0) THEN
                nu_physical     = delta_t*SQRT(sigma_Inv*sigma_Inv)
                delta_tau       = (nu)/SQRT(sigma_Inv*sigma_Inv &
                                    + (1.0_rDef/delta_t)&
                                    *(1.0_rDef/delta_t))
            END IF

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
                            delta_x         = delta_x           ,&
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
                             delta_x        = delta_x           ,&
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

            NewtonResidual(nI) = Newtonl2Res(nI)

            IF (nT == 1) THEN
                nI2 = nI
            ELSE
                nI2 = nI2 + 1
            END IF
             
            IF (NewtonResidual(nI) > tol) THEN
                DO i = iMin, iMax
                    DO j = bMin, bMax
                        Q_np1_k(i, j) = Q_np1_kp1(i, j)
                    END DO
                END DO
            ELSE
                EXIT
            END IF
        END DO  !nI

        DO i = iMin, iMax
            DO j = bMin, bMax
                Q_n(i, j)     = Q_np1_kp1(i, j)
            END DO
        END DO

        IF (MOD(nT, INT(Steps_Per_Cycle)) == 1) THEN
            Step_Number = Step_Number+1 
        END IF

        time            = time     +   delta_t
        Time_FFT(nT)    = time

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
            Pres_Pet(i)         = (Pres(i) - Pres_Mean(i)) - Mean_FFT1(i)
        END DO

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
                DO i = 1, nnI
                    WRITE(78, *)    i               ,&
                                    Periodic_SS(i)
                END DO
            END IF
        END IF
        !!!!!!! END Calculate The Periodic Steady State !!!!!!!!!!!!!!!!!!!

        IF (MOD(nT, 1) == 1) THEN
            WRITE(6, *) "Writing Out Solution Here. Check File"
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            WRITE(filename5, '(a5, a1, a12, a1, a8, a1, a4, a1, i8.8, a2, a8)')&
                    LeftHandSide,  '_',  LHS_Diss,'_', DifferStencil, '_',           &
                    Dissipation, '_', nT, 'nT', 'Petu.dat' 
            OPEN(58, FILE = TRIM(Current_Directory)//&
                        &'/TimeStepSol/'&
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

        IF (MOD(nT, 1) == 0) THEN
            WRITE(6, *) 'Mean and Petrubation ' ,&
                        LeftHandSide, ' '       ,&
                        LHS_Diss, ' '           ,&
                        DifferStencil, ' '      ,&
                        Dissipation             ,&
                        nT                      ,&
                        nI                      ,&
                        time                    ,&
                        Step_Number             ,&
                        nu_physical             ,&
                        Steps_Per_Cycle         ,&
                        Periodic_SS(nnI)
        END IF

        !Solution is written out at every time step, P vs T 
        WRITE(61, *)    nT                                      ,& 
                        nI                                      ,&
                        nI2                                     ,&
                        Pres_Pet(iMin)                          ,&
                        Pres_Pet(iMax)                          ,&
                        time                                    ,&
                        INT(time/Cycle_Time+1.0_rDef)


        !!! Write Out Solution on the Final Cycle at every Time Step
        !!IF (nT >= Steps_Per_Cycle*numOfCycles-Steps_Per_Cycle) THEN
        !!    WRITE(filename5, '(a5, a1, a12, a1, a8, a1, a4, a1, i8.8, a2, a8)')&
        !!            LeftHandSide,   '_', LHS_Diss,'_', DifferStencil, '_',           &
        !!            Dissipation, '_', nT, 'nT', 'Petu.dat' 
        !!    OPEN(58, FILE = TRIM(Current_Directory)//&
        !!                &'/TimeStepSol/'&
        !!                &//filename5, FORM = 'FORMATTED')
        !!    DO i = iMin,   iMax
        !!        WRITE(58, *) x(i)                               ,&
        !!                    !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
        !!                    Den_Pet(i)                          ,&
        !!                    !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
        !!                    Vel_Pet(i)                          ,&
        !!                    !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
        !!                    Pres_Pet(i)                         ,&
        !!                    time                                ,&
        !!                    Step_Number
        !!    END DO
        !!    CLOSE(58)
        !!END IF
        

        !!!!!!! FFT Calculation Starts Here !!!!!!!!!!!!!!!!!
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
            Max_Pressure(nC, i) = Pres_Pet(i)
            Min_Pressure(nC, i) = Pres_Pet(i)
        END DO

        DO i = iMin, iMax
            Mean_Flow_FFT(nC, i)        = Pres_Pet(i)
            Mean_Flow_FFT_2P(nC2, i)    = Pres_Pet(i)
        END DO
        
        Check_A(nC)         = A2_Inflow


        IF (MOD(nT, Steps_Per_Cycle) == 0) THEN
            WRITE(filename2, '(a5, a1, a12, a1, a8, a1, a4, a6, i3.3, a12)')&
                    LeftHandSide,  '_', LHS_Diss, '_', DifferStencil, '_',      &
                    Dissipation, '_Step_', Step_Number,         &
                    '_PvsTime.dat'
            OPEN(72, FILE = TRIM(Current_Directory)//&
                        &'/FFT/'&
                        &//filename2, FORM = 'FORMATTED')
            DO nC = 1, Steps_Per_Cycle
                WRITE(72, *)    REAL(nC-1, rDef)*&
                                Cycle_Time/&
                                REAL(Steps_Per_Cycle-1, rDef)       ,&
                                Mean_Flow_FFT(nC, iMax)     !Exit Pressure 
            END DO
            CLOSE(72)
        END IF

        ! Update The Mean Flow of the peturbation to 'filter-out' the solution
        IF (MOD(nT, INT(Steps_Per_Cycle)) == 1) THEN
            IF (Step_Number >= 1) THEN
                DO i = iMin, iMax
                    Mean_FFT1(i) = Mean_FFT1(i) +                   &
                                    (MAXVAL(Max_Pressure(:, i))     +& 
                                    MINVAL(Min_Pressure(:, i)))&
                                    /2.0_rDef
                END DO
            ELSE
                DO i = iMin, iMax
                    Mean_FFT1(i) = Mean_FFT1(i)
                END DO
            END IF
        END IF
        
        !!! Write Out Solution for the Final Cycle !!!
        IF (nT >= Steps_Per_Cycle*numOfCycles-Steps_Per_Cycle) THEN
            IF (nT == numTimeSteps) THEN
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                WRITE(filename9, '(a5, a1, a12, a1, a8, a1, a4, a1, i3.3, a17)')&
                        LeftHandSide,  '_', LHS_Diss, '_', DifferStencil, '_',           &
                        Dissipation, '_', Steps_Per_Cycle, '_Max_Min_Pres.dat'
                OPEN(76, FILE = TRIM(Current_Directory)//&
                            &'/TimeStepSol/'&
                            &//filename9, FORM = 'FORMATTED')
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                DO i = iMin, iMax
                    ! Get Max and Min Pressure Envelope along with Mean Flow
                    WRITE(76, *)    x(i)                                ,&
                                    !!!!!! Max Pressure Envelopee !!!!!!!!!
                                    MAXVAL(Max_Pressure(:, i))          ,& 
                                    !!!!!! Min Pressure Envelopee !!!!!!!!!
                                    MINVAL(Min_Pressure(:, i))          ,& 
                                    !!!!!! Mean Pressure Envelopee !!!!!!!!
                                    (MAXVAL(Max_Pressure(:, i))  +       & 
                                    MINVAL(Min_Pressure(:, i)))/2.0_rDef, &
                                    !!!!!! Max Pressure Disturbance !!!!!!!
                                    MAXVAL(Max_Pressure(:, i))   -       & 
                                    (MAXVAL(Max_Pressure(:, i))  +       & 
                                    MINVAL(Min_Pressure(:, i)))/2.0_rDef
                END DO
                CLOSE(76)
                WRITE(filename2, '(a5, a1, a12, a1, a8, a1, a4,a12)')&
                        LeftHandSide,  '_', LHS_Diss, '_', DifferStencil, '_',      &
                        Dissipation, '_PvsTime.dat'
                OPEN(72, FILE = TRIM(Current_Directory)//&
                            &'/TimeStepSol/'&
                            &//filename2, FORM = 'FORMATTED')
                DO nC = 1, Steps_Per_Cycle
                    WRITE(72, *)    REAL(nC-1, rDef)*&
                                    Cycle_Time/&
                                    REAL(Steps_Per_Cycle-1, rDef)       ,&
                                    Mean_Flow_FFT(nC, iMax)     !Exit Pressure 
                END DO
                CLOSE(72)
                END IF
        END IF
        DEALLOCATE(bet)

        IF (nT >= numOfCycles*Steps_Per_Cycle) THEN
            WRITE(6, *) "Periodic Steady-State Reached"
            CLOSE(61)
            CLOSE(78)
            EXIT
        END IF


    END DO

    DEALLOCATE( Q_Store             ,&
                TIme_FFT            ,&
                Q_Peturb_Out        ,&
                l2Res               ,&
                Residual            ,&
                MaxWaveSpeed        ,&
                Min_Pressure        ,&
                Max_Pressure        ,&
                Max_Pressure_Exact  ,&
                x_Exact             ,&
                All_FFT             ,&
                Mean_Flow_FFT       ,&
                Mean_FFT1           ,&
                Mean_Flow_FFT_2P    ,&
                R_FFT               ,&
                I_FFT               ,&
                M_FFT               ,&
                R_FFT_2P            ,&
                I_FFT_2P            ,&
                M_FFT_2P            ,&
                Periodic_Res        ,&
                Periodic_SS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END SUBROUTINE Mean_And_Pet

END MODULE GetMean_And_Pet
