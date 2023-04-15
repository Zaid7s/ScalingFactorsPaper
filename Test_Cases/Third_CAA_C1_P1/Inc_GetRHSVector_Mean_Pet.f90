
    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS       ,&
                           nD       ,&
                           nStages  ,&
                           nDis     ,&
                           nBound
    
    INTEGER, INTENT(INOUT):: nTau

    REAL(KIND = rDef), INTENT(IN):: delta_x     ,&
                                     delta_t    ,&
                                     delta_tau  ,&
                                     time

    REAL(KIND = rDef), INTENT(INOUT):: Mach        ,&
                                        u           ,&
                                        p           ,&
                                        c           ,&
                                        gm1         ,&
                                        rho         ,&
                                        Newtonl2Res, &
                                        rho_Static  ,&
                                        P_Static
    
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

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS         ,&
                                                            dQdTau      ,&
                                                            Dis         ,&
                                                            Q_n         ,&
                                                            Q_Exact     ,&
                                                            Q_np1_k     ,&
                                                            Q_nm2       ,&
                                                            Q_nm1

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  dEdX            ,&
                                                        dQdT            ,&
                                                        Q_n0            ,&
                                                        E               ,&
                                                        S

    REAL(KIND = rDef) ::    l2Fac       ,&
                            DissiFac    ,&
                            Mean_p      ,&
                            Mean_rho    ,&
                            Mean_u      ,&
                            Mean_c      ,&
                            Mean_Mach   ,&
                            Contra_U    ,&
                            E_tot       ,&
                            AreaFac


    INTEGER:: i, j
    
    ALLOCATE(   E(iMax, bMax)       ,&
                dEdX(iMax, bMax)    ,&
                dQdT(iMax, bMax)    ,&
                Q_n0(iMax, bMax)    ,&
                S(iMax, bMax))

    DissiFac    = 1.0_rDef/delta_t
    
! Initialize E such that it is easier to calculate dEdX
    INCLUDE 'Inc_InitializeFluxes.f90'

    CALL Dissipation_RHS(   iMin        = iMin          ,&
                            iMax        = iMax          ,&
                            bMin        = bMin          ,&
                            bMax        = bMax          ,&
                            delta_X     = delta_X       ,&
                            delta_tau   = delta_t       ,&
                            nD          = nD            ,&
                            Dis         = Dis           ,&
                            Q_n         = Q_np1_k       ,&
                            Jac_Curv    = Jac_Curv) 
    
    CALL dEdXVector(    iMin            = iMin              ,&
                        iMax            = iMax              ,&
                        bMin            = bMin              ,&
                        bMax            = bMax              ,&
                        E               = E                 ,&
                        time            = time              ,&
                        dEdX            = dEdX              ,&
                        delta_X         = delta_X           ,&
                        delta_tau       = delta_tau         ,&
                        delta_t         = delta_t           ,&
                        DS              = DS                ,&
                        Dis             = Dis               ,&
                        Mean_Mach       = Mean_Mach         ,&
                        Mean_u          = Mean_u            ,&
                        Mean_p          = Mean_p            ,&
                        Mean_c          = Mean_c            ,&
                        gm1             = gm1               ,&
                        Mean_rho        = Mean_rho          ,&
                        rho_Static      = rho_Static        ,&
                        P_Static        = P_Static          ,&
                        Q_np1_k         = Q_np1_k           ,&
                        Q_Exact         = Q_Exact           ,&
                        Jac_Curv        = Jac_Curv          ,&
                        Inv_Jac_Curv    = Inv_Jac_Curv      ,&
                        dZidX           = dZidX             ,&
                        dZidt           = dZidt             ,&
                        A1_Inflow       = A1_Inflow         ,&
                        A2_Inflow       = A2_Inflow         ,&
                        A3_Inflow       = A3_Inflow         ,&
                        A1_Outflow      = A1_Outflow        ,&
                        A2_Outflow      = A2_Outflow        ,&
                        A3_Outflow      = A3_Outflow        ,&
                        nBound          = nBound)

    ! CALL dQdT_BDF_Mean_Pet
    CALL dQdT_BDF(  iMin        =  iMin         ,&
                    iMax        =  iMax         ,&
                    bMin        =  bMin         ,&
                    bMax        =  bMax         ,&
                    Q_np1_k     =  Q_np1_k      ,&
                    Q_n         =  Q_n          ,&
                    Q_nm1       =  Q_nm1        ,&
                    Q_nm2       =  Q_nm2        ,&
                    delta_t     =  delta_t      ,&
                    dQdT        =  dQdT         ,&
                    ndQdT       =  ndQdT)

    DO j = bMin, bMax
        DO i = iMin, iMax
            RHS(i, j) = (-(dQdT(i, j)*REAL(nBound, rDef) + dEdX(i, j) - S(i, j) &
                            - Dis(i, j)))
        END DO
    END DO

    l2Fac = 0.0_rDef
    DO j = bMin, bMax
        DO i = iMin, iMax
            l2Fac = l2Fac+RHS(i, j)*RHS(i, j)
        END DO
    END DO
    Newtonl2Res = SQRT(l2Fac)/REAL((iMax)*3, rDef)
