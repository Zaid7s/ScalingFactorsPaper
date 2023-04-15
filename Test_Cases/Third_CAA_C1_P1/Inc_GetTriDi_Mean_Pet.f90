
    INTEGER, INTENT(IN)::   iMin        ,&
                            iMax        ,&
                            bMax        ,&
                            bMin        ,&
                            DS          ,&
                            nD          ,&
                            nStages     ,&
                            nDis        ,&
                            nBound

    INTEGER, INTENT(INOUT):: nTau

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  bet             ,&
                                                    Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt           ,&
                                                    Area


    REAL(KIND = rDef), INTENT(OUT) ::   A1_Inflow       ,&
                                        A2_Inflow       ,&
                                        A3_Inflow       ,&
                                        A1_Outflow      ,&
                                        A2_Outflow      ,&
                                        A3_Outflow
    
    REAL(KIND = rDef), INTENT(IN):: time       ,&
                                    delta_x    ,&
                                    delta_t    ,&
                                    delta_tau  ,&
                                    Scaling_Fac
                                        
    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    S_Jac       ,&
                                                            A_Jac

    REAL(KIND = rDef), INTENT(INOUT):: Newtonl2Res

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS         ,&
                                                            Q_np1_k     ,&
                                                            Q_np1_kp1   ,&
                                                            Q_Exact     ,&
                                                            Dis         ,&
                                                            Q_np1       ,&
                                                            Q_nm1       ,&
                                                            Q_nm2       ,&
                                                            Delta_Q     ,&
                                                            Q_n         ,&
                                                            E           ,&
                                                            E_GPT       ,&
                                                            dQdTau

    REAL(KIND = rDef), INTENT(INOUT)::  Mach        ,&
                                        u           ,&
                                        p           ,&
                                        c           ,&
                                        gm1         ,&
                                        rho         ,&
                                        rho_Static  ,&
                                        P_Static

    INTEGER ::   i, &
                 k, &
                 j, &
                 l, &
                 m, &
                 n

    REAL(KIND = rDef) :: l2Fac

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  RHS2            ,&
                                                        LowerDiag_Bound

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   LowerDiag       ,&
                                                            Diag            ,&
                                                            UpperDiag

    INTEGER, DIMENSION(1):: iiMin, iiMax

    iiMin = iMin
    iiMax = iMax

    ALLOCATE(   LowerDiag(iMax, bMax, bMax), &
                Diag(iMax, bMax, bMax)      ,&
                UpperDiag(iMax, bMax, bMax), &
                RHS2(iMax, bMax)            ,&
                LowerDiag_Bound(bMax, bMax))

    CALL DiagLHS( iMin              = iMin              ,& 
                  iMax              = iMax              ,& 
                  bMin              = bMin              ,& 
                  bMax              = bMax              ,& 
                  delta_t           = delta_t           ,&
                  LowerDiag         = LowerDiag         ,& 
                  Diag              = Diag              ,& 
                  UpperDiag         = UpperDiag         ,&
                  Scaling_Fac       = Scaling_Fac       ,&
                  A_Jac             = A_Jac             ,&
                  S_Jac             = S_Jac             ,&
                  nDis              = nDis              ,&     
                  Source_Fac        = Source_Fac        ,&
                  Jac_Curv          = Jac_Curv          ,&
                  LowerDiag_Bound   = LowerDiag_Bound)


    CALL RungeKutta( iMin            =  iMin            ,&
                     iMax            =  iMax            ,&
                     bMin            =  bMin            ,&
                     bMax            =  bMax            ,&
                     Q_n             =  Q_n             ,&
                     Q_np1_k         =  Q_np1_k         ,&
                     Q_nm1           =  Q_nm1           ,&
                     Q_nm2           =  Q_nm2           ,&
                     dQdTau          =  dQdTau          ,&
                     RHS             =  RHS             ,&
                     time            =  time            ,&
                     delta_x         =  delta_x         ,&
                     delta_tau       =  delta_tau       ,&
                     delta_t         =  delta_t         ,&
                     DS              =  DS              ,&
                     nD              =  nD              ,&
                     Dis             =  Dis             ,&
                     bet             =  bet             ,&
                     nTau            =  nTau            ,&
                     nStages         =  nStages         ,&
                     Newtonl2Res     =  Newtonl2Res     ,&
                     Source_Fac      =  Source_Fac      ,&
                     Mach            =  Mach            ,&
                     u               =  u               ,&
                     p               =  p               ,&
                     c               =  c               ,&
                     gm1             =  gm1             ,&
                     rho             =  rho             ,&
                     rho_Static      =  rho_Static      ,&
                     P_Static        =  P_Static        ,&
                     Q_Exact         =  Q_Exact         ,&
                     Jac_Curv        =  Jac_Curv        ,&
                     Inv_Jac_Curv    =  Inv_Jac_Curv    ,&
                     dZidX           =  dZidX           ,&
                     dZidt           =  dZidt           ,&
                     Area            =  Area            ,&
                     A1_Inflow       =  A1_Inflow       ,&
                     A2_Inflow       =  A2_Inflow       ,&
                     A3_Inflow       =  A3_Inflow       ,&
                     A1_Outflow      =  A1_Outflow      ,&
                     A2_Outflow      =  A2_Outflow      ,&
                     nDis            =  nDis            ,&
                     A3_Outflow      =  A3_Outflow      ,&
                     nBound          =  nBound)


! Get the RHS. This includes the boundaries and dissipation     
    ! CALL RHSVector_Mean_Pet
    CALL RHSVector( iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    Q_n             = Q_n           ,&
                    Q_np1_k         = Q_np1_k       ,&
                    Q_nm1           = Q_nm1         ,&
                    Q_nm2           = Q_nm2         ,&
                    dQdTau          = dQdTau        ,&
                    RHS             = RHS           ,&
                    time            = time          ,&
                    delta_x         = delta_x       ,&
                    delta_tau       = delta_tau     ,&
                    delta_t         = delta_t       ,&
                    DS              = DS            ,&
                    nD              = nD            ,&
                    Dis             = Dis           ,&
                    bet             = bet           ,&
                    nTau            = nTau          ,&
                    nStages         = nStages       ,&
                    Newtonl2Res     = Newtonl2Res   ,&
                    Source_Fac      = Source_Fac    ,&
                    Mach            = Mach          ,&
                    u               = u             ,&
                    p               = p             ,&
                    c               = c             ,&
                    gm1             = gm1           ,&
                    rho             = rho           ,&
                    rho_Static      = rho_Static    ,&
                    P_Static        = P_Static      ,&
                    Q_Exact         = Q_Exact       ,&
                    Jac_Curv        = Jac_Curv      ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv  ,&
                    dZidX           = dZidX         ,&
                    dZidt           = dZidt         ,&
                    Area            = Area          ,&
                    A1_Inflow       = A1_Inflow     ,&
                    A2_Inflow       = A2_Inflow     ,&
                    A3_Inflow       = A3_Inflow     ,&
                    A1_Outflow      = A1_Outflow    ,&
                    A2_Outflow      = A2_Outflow    ,&
                    nDis            = nDis          ,&
                    A3_Outflow      = A3_Outflow    ,&
                    nBound          = nBound)

    CALL FixBoundaries( iMin            =  iMin             ,& 
                        iMax            =  iMax             ,&
                        bMin            =  bMin             ,&
                        bMax            =  bMax             ,&
                        LowerDiag       =  LowerDiag        ,&
                        Diag            =  Diag             ,&
                        UpperDiag       =  UpperDiag        ,&
                        RHS             =  RHS              ,&
                        LowerDiag_Bound =  LowerDiag_Bound)      

!Initialize taQ and dQdTau Vectors
    
    CALL BlockTriMatrix(  numVar       = bMax,      & 
                          iSolveStart  = iiMin,     &
                          iSolveEnd    = iiMax,     &
                          aStart       = iiMin,     &
                          A            = LowerDiag, &
                          bStart       = iiMin,     &
                          B            = Diag,      &
                          cStart       = iiMin,     &
                          C            = UpperDiag, &
                          xStart       = iiMin,     &
                          X            = Delta_Q,   &
                          rStart       = iiMin,     &
                          R            = RHS)
    
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_np1_kp1(i, j) = Delta_Q(i, j) + Q_np1_k(i, j)
        END DO
    END DO

