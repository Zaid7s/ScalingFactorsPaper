
    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS       ,&
                           nBound

    REAL(KIND = rDef), INTENT(INOUT):: Mean_Mach   ,&
                                        Mean_u      ,&
                                        Mean_p      ,&
                                        Mean_c      ,&
                                        gm1         ,&
                                        Mean_rho    ,&
                                        rho_Static  ,&
                                        P_Static

    REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                     delta_tau  ,&
                                     delta_t    ,&
                                     time

    REAL(KIND = rDef), INTENT(OUT) ::   A1_Inflow       ,&
                                        A2_Inflow       ,&
                                        A3_Inflow       ,&
                                        A1_Outflow      ,&
                                        A2_Outflow      ,&
                                        A3_Outflow
                                        
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::   E       ,&
                                                        Dis     ,&
                                                        Q_np1_k, &
                                                        Q_Exact

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: dEdX

    INTEGER:: i, j
    

!Inflow Boundary Using Thmpson Style BC and GPT
    DO j = bMin, bMax
        DO i = iMin, iMax 
            dEdX(i, j) = 0.0_rDef 
        END DO
    END DO

    IF ((DS == 6).OR.(DS == 7).OR.(DS == 8)) THEN
        CALL CompactScheme( iMin       =  iMin      ,&
                            iMax       =  iMax      ,&
                            bMin       =  bMin      ,&
                            bMax       =  bMax      ,&
                            E          =  E         ,&
                            dEdX       =  dEdX      ,&
                            delta_x    =  delta_Zi  ,&
                            delta_tau  =  delta_tau, &
                            DS         =  DS        ,&
                            dZidX      =  dZidX) 
    ELSE
        DO j = bMin, bMax
            DO i = iMin, iMax 
                dEdX(i, j) = 0.0_rDef 
            END DO
        END DO
    END IF

    ! Call InflowBC_Mean_Pet
    CALL InflowBC(  iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    E               = E             ,&
                    time            = time          ,&
                    dEdX            = dEdX          ,&
                    delta_x         = delta_Zi      ,&
                    delta_t         = delta_t       ,&
                    delta_tau       = delta_tau     ,&
                    DS              = DS            ,&
                    Mean_Mach       = Mean_Mach     ,&
                    Mean_u          = Mean_u        ,&
                    Mean_p          = Mean_p        ,&
                    Mean_c          = Mean_c        ,&
                    gm1             = gm1           ,&
                    Mean_rho        = Mean_rho      ,&
                    rho_Static      = rho_Static    ,&
                    P_Static        = P_Static      ,&
                    Q_np1_k         = Q_np1_k       ,&
                    Q_Exact         = Q_Exact       ,&
                    Jac_Curv        = Jac_Curv      ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv  ,&
                    dZidX           = dZidX         ,&
                    A_S_BC          = A1_Inflow     ,&
                    A_Plus_BC       = A2_Inflow     ,&
                    A_Minus_BC      = A3_Inflow     ,&
                    nBound          = nBound)

    ! Call OutflowBC_Mean_Pet
!Outflow Boundary Using Thmpson Style BC and GPT
    CALL OutflowBC( iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    E               = E             ,&
                    time            = time          ,&
                    dEdX            = dEdX          ,&
                    delta_x         = delta_Zi      ,&
                    delta_t         = delta_t       ,&
                    delta_tau       = delta_tau     ,&
                    DS              = DS            ,&
                    Mean_Mach       = Mean_Mach     ,&
                    Mean_u          = Mean_u        ,&
                    Mean_p          = Mean_p        ,&
                    Mean_c          = Mean_c        ,&
                    gm1             = gm1           ,&
                    Mean_rho        = Mean_rho      ,&
                    rho_Static      = rho_Static    ,&
                    P_Static        = P_Static      ,&
                    Q_np1_k         = Q_np1_k       ,&
                    Q_Exact         = Q_Exact       ,&
                    Jac_Curv        = Jac_Curv      ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv  ,&
                    dZidX           = dZidX         ,&    
                    A_S_BC          = A1_Outflow    ,&
                    A_Plus_BC       = A2_Outflow    ,&
                    A_Minus_BC      = A3_Outflow    ,&
                    nBound          = nBound)

! dEdX for Interior Domain
    IF ((DS == 6).OR.(DS == 7).OR.(DS == 8)) THEN
        CALL CompactScheme( iMin       =  iMin      ,&
                            iMax       =  iMax      ,&
                            bMin       =  bMin      ,&
                            bMax       =  bMax      ,&
                            E          =  E         ,&
                            time       =  time      ,&
                            dEdX       =  dEdX      ,&
                            delta_x    =  delta_Zi  ,&
                            delta_tau  =  delta_tau, &
                            DS         =  DS        ,&
                            dZidX      =  dZidX) 
    ELSE
        CALL SpatialDerivative( iMin       =  iMin      ,&
                                iMax       =  iMax      ,&
                                bMin       =  bMin      ,&
                                bMax       =  bMax      ,&
                                E          =  E         ,&
                                dEdX       =  dEdX      ,&
                                delta_x    =  delta_Zi  ,&
                                delta_tau  =  delta_tau, &
                                DS         =  DS        ,&
                                dZidX      =  dZidX) 
    END IF
