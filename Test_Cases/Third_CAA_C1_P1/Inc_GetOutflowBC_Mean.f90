    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS

    REAL(KIND = rDef), INTENT(INOUT):: Mean_Mach   ,&
                                        Mean_u      ,&
                                        Mean_p      ,&
                                        Mean_c      ,&
                                        gm1         ,&
                                        Mean_rho    ,&
                                        rho_Static  ,&
                                        P_Static
    
    REAL(KIND = rDef), INTENT(IN)::  delta_x    ,&
                                     delta_t    ,&
                                     delta_tau  ,&
                                     time
                                        
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX
                                        
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::   E       ,&
                                                        Q_np1_k, &
                                                        Q_Exact

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: dEdX

    INTEGER:: i, j, k

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE:: dEdXnoBC1       ,&
                                                    dQdTnoBC        ,&
                                                    dQdT            ,&
                                                    dEdX1           ,&
                                                    dEdX2           ,&
                                                    dEdX3           ,&
                                                    dEdXnoBC2       ,&
                                                    dEdXnoBC3       ,&
                                                    delta_dEdX      ,&
                                                    dEdX_Thompson   ,&
                                                    p               ,&
                                                    rho             ,&
                                                    u
    
    REAL(KIND = rDef) ::    A_S_BC          ,&
                            A_Plus_BC       ,&
                            A_Minus_BC      ,&
                            A_S_noBC        ,&
                            A_Plus_noBC     ,&
                            A_Minus_noBC    ,&
                            p_T             ,&
                            rho_T           ,&
                            u_T             ,&
                            p_TnoBC         ,&
                            rho_TnoBC       ,&
                            u_TnoBC         ,&
                            p_Zi            ,&
                            rho_Zi          ,&
                            u_Zi            ,&
                            Delta_A_S       ,&
                            Delta_A_Plus    ,&
                            P_Outflow       ,&
                            omega           ,&
                            epsi            ,&
                            Bar_u           ,&
                            Bar_p           ,&
                            Bar_rho         ,&
                            tau

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   LowerDiag       ,&
                                                            UpperDiag       ,&
                                                            MainDiag        ,&
                                                            LowerDiag_var   ,&
                                                            UpperDiag_var   ,&
                                                            MainDiag_var

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  RHS     ,&
                                                        TriDi   ,&
                                                        RHS_var

    INTEGER, DIMENSION(1):: iiMin, iiMax

    iiMin = iMin
    iiMax = iMax

    ALLOCATE(   dEdXnoBC1(bMax)     ,&
                dQdTnoBC(bMax)      ,&
                dQdT(bMax)          ,&
                dEdX1(bMax)         ,&
                dEdXnoBC2(bMax)     ,&
                dEdXnoBC3(bMax)     ,&
                delta_dEdX(bMax)    ,&
                dEdX_Thompson(bMax), &
                dEdX2(bMax)         ,&
                dEdX3(bMax)         ,&
                p(iMax)             ,&
                rho(iMax)           ,&
                u(iMax))

! dEdX for Interior Domain
    omega       = 0.6_rDef*pi
    epsi        = 1E-08_rDef
    
    gm1         = 1.0_rDef/(gam-1.0_rDef)

    tau         = delta_t
    Bar_rho    = Q_np1_k(iMax, 1)*Jac_Curv(iMax)
    Bar_u      = (Q_np1_k(iMax, 2)/Q_np1_k(iMax, 1))
    Bar_p      = (gam-1.0_rDef)*(Q_np1_k(iMax, 3) - &
                   0.5_rDef*Q_np1_k(iMax, 2)*Q_np1_k(iMax, 2)/Q_np1_k(iMax, 1))&
                   *Jac_Curv(iMax)

    Mean_rho    = (1.0_rDef/tau)*(Bar_rho*(time+tau) -  &
                                    Bar_rho*time)
    Mean_u      = (1.0_rDef/tau)*(Bar_u*(time+tau)  -   &
                                    Bar_u*time)
    Mean_p      = (1.0_rDef/tau)*(Bar_p*(time+tau)  -   &
                                    Bar_p*time)
    Mean_c      = SQRT(gam*Mean_p/Mean_rho)
    Mean_Mach   = Mean_u/Mean_c

    rho_TnoBC   = 0.0_rDef  
    u_TnoBC     = 0.0_rDef  
    p_TnoBC     = 0.0_rDef  

    DO i = iMin, iMax
        rho(i)     = Q_np1_k(i, 1)*Jac_Curv(i)
        u(i)       = Q_np1_k(i, 2)/Q_np1_k(i, 1)
        p(i)       = (gam-1.0_rDef)*(Q_np1_k(i, 3) - &
                       0.5_rDef*Q_np1_k(i, 2)*Q_np1_k(i, 2)/&
                           Q_np1_k(i, 1))*Jac_Curv(i)              
    END DO

    P_Outflow = p(iMax)/dZidX(iMax)
     
    IF ((DS == 1)) THEN
        DO j = bMin, bMax
            i = iMax
            dEdXnoBC1(j) = -(-3.0_rDef*E(i, j)       &
                            +4.0_rDef*E(i-1, j)   &
                            -         E(i-2, j))  &
                                /(2.0_rDef*delta_x)
        END DO
        i       = iMax
        rho_Zi   = -(-3.0_rDef*rho(i-0)    &
                    +4.0_rDef*rho(i-1)    &
                    -         rho(i-2))   &
                        /(2.0_rDef*delta_x)
        
        u_Zi     = -(-3.0_rDef*u(i-0)      &
                    +4.0_rDef*u(i-1)      &
                    -         u(i-2))     &
                        /(2.0_rDef*delta_x)

        p_Zi     = -(-3.0_rDef*p(i-0)      &
                    +4.0_rDef*p(i-1)      &
                    -         p(i-2))     &
                        /(2.0_rDef*delta_x)

        rho_TnoBC   = -1.0_rDef*(dZidX(i)*u(i)*rho_Zi               + &
                                rho(i)*dZidX(i)*u_Zi)
        u_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*u_Zi                 + &
                                (1.0_rDef/rho(i))*dZidX(i)*p_Zi)
        p_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*p_Zi                 + &
                                    gam*p(i)*dZidX(i)*u_Zi)

        INCLUDE 'IncApplyOutflow.f90'

        DO j = bMin, bMax
            dEdX_Thompson(j)    = -dQdT(j)
            delta_dEdX(j)       = dEdX_Thompson(j) - dEdXnoBC1(j)
            dEdX1(j)            = dEdXnoBC1(j) + delta_dEdX(j)
            dEdX(iMax, j)       = dEdX1(j)
        END DO
    ELSE IF (DS == 2) THEN
        DO j = bMin, bMax
            i = iMax
            dEdXnoBC1(j) = -(-25.0_rDef*E(i-0, j)   &
                             +48.0_rDef*E(i-1, j)    &
                             -36.0_rDef*E(i-2, j)   &
                             +16.0_rDef*E(i-3, j)   &
                              -3.0_rDef*E(i-4, j))  &
                                /(12.0_rDef*delta_x)

            rho_Zi = -(-25.0_rDef*rho(i-0)  &
                      +48.0_rDef*rho(i-1)  &
                      -36.0_rDef*rho(i-2)  &
                      +16.0_rDef*rho(i-3)  &
                       -3.0_rDef*rho(i-4)) &
                        /(12.0_rDef*delta_x)

            u_Zi   = -(-25.0_rDef*u(i-0)  &
                      +48.0_rDef*u(i-1)  &
                      -36.0_rDef*u(i-2)  &
                      +16.0_rDef*u(i-3)  &
                       -3.0_rDef*u(i-4)) &
                         /(12.0_rDef*delta_x)

            p_Zi   = -(-25.0_rDef*p(i-0)  &
                      +48.0_rDef*p(i-1)  &
                      -36.0_rDef*p(i-2)  &
                      +16.0_rDef*p(i-3)  &
                       -3.0_rDef*p(i-4)) &
                         /(12.0_rDef*delta_x)

            rho_TnoBC   = -1.0_rDef*(dZidX(i)*u(i)*rho_Zi               + &
                                    rho(i)*dZidX(i)*u_Zi)
            u_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*u_Zi                 + &
                                    (1.0_rDef/rho(i))*dZidX(i)*p_Zi)
            p_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*p_Zi                 + &
                                    gam*p(i)*dZidX(i)*u_Zi)
        END DO

        INCLUDE 'IncApplyOutflow.f90'

        DO j = bMin, bMax
            i = iMax-1
            dEdXnoBC2(j)    =  -(-3.0_rDef*E(i  + 1, j)    &
                               -10.0_rDef*E(i  - 0, j)    &
                                +18.0_rDef*E(i-1, j)   &
                                 -6.0_rDef*E(i-2, j)   &
                                 +1.0_rDef*E(i-3, j))  &
                                   /(12.0_rDef*delta_x)
        END DO

        DO j = bMin, bMax
            dEdX_Thompson(j)     = -dQdT(j)  
            delta_dEdX(j)        = dEdX_Thompson(j) - dEdXnoBC1(j)
            dEdX1(j)             = dEdXnoBC1(j) + delta_dEdX(j)
            dEdX2(j)             = dEdXnoBC2(j) - &
                                    (1.0_rDef/3.0_rDef)*delta_dEdX(j)
            dEdX(iMax, j)        = dEdX1(j)
            dEdX(iMax-1, j)    = dEdX2(j)
        END DO
    ELSE IF (DS == 3) THEN
        DO j = bMin, bMax
            i = iMax
            dEdXnoBC1(j) = -(-147.0_rDef*E(i+0, j)  &
                           +360.0_rDef*E(i-1, j)   &
                           -450.0_rDef*E(i-2, j)   &
                           +400.0_rDef*E(i-3, j)   &
                           -225.0_rDef*E(i-4, j)   &
                            +72.0_rDef*E(i-5, j)   &
                            -10.0_rDef*E(i-6, j))  &
                                /(60.0_rDef*delta_x)

            rho_Zi        = -(-147.0_rDef*rho(i+0)  &
                           +360.0_rDef*rho(i-1)   &
                           -450.0_rDef*rho(i-2)   &
                           +400.0_rDef*rho(i-3)   &
                           -225.0_rDef*rho(i-4)   &
                            +72.0_rDef*rho(i-5)   &
                            -10.0_rDef*rho(i-6))  &
                                /(60.0_rDef*delta_x)

            u_Zi          = -(-147.0_rDef*u(i+0)  &
                           +360.0_rDef*u(i-1)   &
                           -450.0_rDef*u(i-2)   &
                           +400.0_rDef*u(i-3)   &
                           -225.0_rDef*u(i-4)   &
                            +72.0_rDef*u(i-5)   &
                            -10.0_rDef*u(i-6))  &
                                /(60.0_rDef*delta_x)

            p_Zi          = -(-147.0_rDef*p(i+0)  &
                           +360.0_rDef*p(i-1)   &
                           -450.0_rDef*p(i-2)   &
                           +400.0_rDef*p(i-3)   &
                           -225.0_rDef*p(i-4)   &
                            +72.0_rDef*p(i-5)   &
                            -10.0_rDef*p(i-6))  &
                                /(60.0_rDef*delta_x)

            rho_TnoBC   = -1.0_rDef*(dZidX(i)*u(i)*rho_Zi               + &
                                    rho(i)*dZidX(i)*u_Zi)
            u_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*u_Zi                 + &
                                    (1.0_rDef/rho(i))*dZidX(i)*p_Zi)
            p_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*p_Zi                 + &
                                    gam*p(i)*dZidX(i)*u_Zi)
        END DO

        INCLUDE 'IncApplyOutflow.f90'
        
        DO j = bMin, bMax
            i = iMax-1
            dEdXnoBC2(j) = -(-10.0_rDef*E(i+1, j)   &
                            -77.0_rDef*E(i+0, j)   &
                           +150.0_rDef*E(i-1, j)   &
                           -100.0_rDef*E(i-2, j)   &
                            +50.0_rDef*E(i-3, j)   &
                            -15.0_rDef*E(i-4, j)    &
                             +2.0_rDef*E(i-5, j))  &
                                /(60.0_rDef*delta_x)
        END DO

        DO j = bMin, bMax
            i = iMax-2
            dEdXnoBC3(j) = -( 2.0_rDef*E(i+2, j)   &
                            -24.0_rDef*E(i+1, j)   &
                            -35.0_rDef*E(i+0, j)   &
                            +80.0_rDef*E(i-1, j)   &
                            -30.0_rDef*E(i-2, j)   &
                             +8.0_rDef*E(i-3, j)   &
                             -1.0_rDef*E(i-4, j))  &
                                /(60.0_rDef*delta_x)
        END DO

        DO j = bMin, bMax
            dEdX_Thompson(j)    = -dQdT(j)  
            delta_dEdX(j)       = dEdX_Thompson(j) - dEdXnoBC1(j)
            dEdX1(j)             = dEdXnoBC1(j) + delta_dEdX(j)
            dEdX2(j)            = dEdXnoBC2(j) - &
                                    (1.0_rDef/5.0_rDef)*delta_dEdX(j)
            dEdX3(j)            = dEdXnoBC3(j) + &
                                    (1.0_rDef/10.0_rDef)*delta_dEdX(j)
            dEdX(iMax, j)        = dEdX1(j)
            dEdX(iMax-1, j)    = dEdX2(j)
            dEdX(iMax-2, j)    = dEdX3(j)
        END DO
    ELSE IF (DS == 4) THEN
        DO j = bMin, bMax
            i = iMax
            dEdXnoBC1(j) =-(-119.0_rDef*E(i+0, j)   &
                           +296.0_rDef*E(i-1, j)   &
                           -379.0_rDef*E(i-2, j)   &
                           +344.0_rDef*E(i-3, j)   &
                           -197.0_rDef*E(i-4, j)   &
                            +64.0_rDef*E(i-5, j)    &
                            -9.0_rDef*E(i  - 6, j))  &
                                /(48.0_rDef*delta_x)

            rho_Zi    =     -(-119.0_rDef*rho(i+0)    &
                             +296.0_rDef*rho(i-1)    &
                             -379.0_rDef*rho(i-2)    &
                             +344.0_rDef*rho(i-3)    &
                             -197.0_rDef*rho(i-4)    &
                              +64.0_rDef*rho(i-5)    &
                               -9.0_rDef*rho(i-6))   &
                                /(48.0_rDef*delta_x)

            u_Zi      =     -(-119.0_rDef*u(i+0)    &
                             +296.0_rDef*u(i-1)    &
                             -379.0_rDef*u(i-2)    &
                             +344.0_rDef*u(i-3)    &
                             -197.0_rDef*u(i-4)    &
                              +64.0_rDef*u(i-5)    &
                               -9.0_rDef*u(i-6))   &
                                /(48.0_rDef*delta_x)

            p_Zi      =     -(-119.0_rDef*p(i+0)    &
                             +296.0_rDef*p(i-1)    &
                             -379.0_rDef*p(i-2)    &
                             +344.0_rDef*p(i-3)    &
                             -197.0_rDef*p(i-4)    &
                              +64.0_rDef*p(i-5)    &
                               -9.0_rDef*p(i-6))   &
                                /(48.0_rDef*delta_x)
            
            rho_TnoBC   = -1.0_rDef*(dZidX(i)*u(i)*rho_Zi               + &
                                    rho(i)*dZidX(i)*u_Zi)
            u_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*u_Zi                 + &
                                    (1.0_rDef/rho(i))*dZidX(i)*p_Zi)
            p_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*p_Zi                 + &
                                    gam*p(i)*dZidX(i)*u_Zi)
        END DO

        INCLUDE 'IncApplyOutflow.f90'

        DO j = bMin, bMax
            i = iMax-1
            dEdXnoBC2(j) = -( -9.0_rDef*E(i+1, j)   &
                             -56.0_rDef*E(i+0, j)   &
                            +107.0_rDef*E(i-1, j)   &
                             -64.0_rDef*E(i-2, j)   &
                             +29.0_rDef*E(i-3, j)   &
                              -8.0_rDef*E(i-4, j)   &
                              +1.0_rDef*E(i-5, j))  &
                                /(48.0_rDef*delta_x)
        END DO

        DO j = bMin, bMax
            i = iMax-2
            dEdXnoBC3(j) = -(  1.0_rDef*E(i+2, j)    &
                             -16.0_rDef*E(i+1, j)    &
                             -35.0_rDef*E(i+0, j)    &
                             +72.0_rDef*E(i-1, j)    &
                             -29.0_rDef*E(i-2, j)    &
                              +8.0_rDef*E(i-3, j)    &
                              -1.0_rDef*E(i-4, j))   &
                                /(48.0_rDef*delta_x)
        END DO

        DO j = bMin, bMax
            dEdX_Thompson(j)    = -dQdT(j)  
            delta_dEdX(j)       = dEdX_Thompson(j) - dEdXnoBC1(j)
            dEdX1(j)             = dEdXnoBC1(j) + delta_dEdX(j)
            dEdX2(j)            = dEdXnoBC2(j) - &
                                    (1.0_rDef/9.0_rDef)*delta_dEdX(j)
            dEdX3(j)            = dEdXnoBC3(j) + &
                                    (1.0_rDef/9.0_rDef)*delta_dEdX(j)
            dEdX(iMax, j)        = dEdX1(j)
            dEdX(iMax-1, j)    = dEdX2(j)
            dEdX(iMax-2, j)    = dEdX3(j)
        END DO
    ELSE IF (DS == 5) THEN
        DO j = bMin, bMax
            i = iMax
            dEdXnoBC1(j) =-(-2.192280339_rDef*E(i-0, j)   &
                            +4.748611401_rDef*E(i-1, j)   &
                            -5.108851915_rDef*E(i-2, j)   &
                            +4.461567104_rDef*E(i-3, j)   &
                            -2.833498741_rDef*E(i-4, j)   &
                            +1.128328861_rDef*E(i-5, j)   &
                            -0.203876371_rDef*E(i-6, j))  &
                                /(delta_x)

            rho_Zi      = -(-2.192280339_rDef*rho(i-0)    &
                            +4.748611401_rDef*rho(i-1)    &
                            -5.108851915_rDef*rho(i-2)    &
                            +4.461567104_rDef*rho(i-3)    &
                            -2.833498741_rDef*rho(i-4)    &
                            +1.128328861_rDef*rho(i-5)    &
                            -0.203876371_rDef*rho(i-6))   &
                               /(delta_x)

            u_Zi      = -(  -2.192280339_rDef*u(i-0)    &
                            +4.748611401_rDef*u(i-1)    &
                            -5.108851915_rDef*u(i-2)    &
                            +4.461567104_rDef*u(i-3)    &
                            -2.833498741_rDef*u(i-4)    &
                            +1.128328861_rDef*u(i-5)    &
                            -0.203876371_rDef*u(i-6))   &
                               /(delta_x)

            p_Zi      = -(  -2.192280339_rDef*p(i-0)    &
                            +4.748611401_rDef*p(i-1)    &
                            -5.108851915_rDef*p(i-2)    &
                            +4.461567104_rDef*p(i-3)    &
                            -2.833498741_rDef*p(i-4)    &
                            +1.128328861_rDef*p(i-5)    &
                            -0.203876371_rDef*p(i-6))   &
                               /(delta_x)
            
            rho_TnoBC   = -1.0_rDef*(dZidX(i)*u(i)*rho_Zi               + &
                                    rho(i)*dZidX(i)*u_Zi)
            u_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*u_Zi                 + &
                                    (1.0_rDef/rho(i))*dZidX(i)*p_Zi)
            p_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*p_Zi                 + &
                                    gam*p(i)*dZidX(i)*u_Zi)
        END DO

        INCLUDE 'IncApplyOutflow.f90'

        DO j = bMin, bMax
            i = iMax-1
            dEdXnoBC2(j) =-(-0.20940637445552380_rDef*E(i+1, j)   &
                            -1.08406399497536900_rDef*E(i+0, j)   &
                            +2.14474893675524100_rDef*E(i-1, j)   &
                            -1.38356163287698400_rDef*E(i-2, j)   &
                            +0.76392684656023610_rDef*E(i-3, j)   &
                            -0.27940632043194790_rDef*E(i-4, j)   &
                            +0.04776253942434876_rDef*E(i-5, j))  &
                               /(delta_x)

        END DO

        DO j = bMin, bMax
            i = iMax-2
            dEdXnoBC3(j) = -(0.05644481863938650_rDef*E(i+2, j)    &
                            -0.51192471192015800_rDef*E(i+1, j)    &
                            -0.37038205460987970_rDef*E(i+0, j)    &
                            +1.13854563013693000_rDef*E(i-1, j)    &
                            -0.42076972392956480_rDef*E(i-2, j)    &
                            +0.12838542345539630_rDef*E(i-3, j)    &
                            -0.02029938177210999_rDef*E(i-4, j))   &
                               /(delta_x)
        END DO

        DO j = bMin, bMax
            dEdX_Thompson(j)    = -dQdT(j)
            delta_dEdX(j)       = dEdX_Thompson(j) - dEdXnoBC1(j)
            dEdX1(j)            = dEdXnoBC1(j) + delta_dEdX(j)
            dEdX2(j)            = dEdXnoBC2(j) - &
                                    (0.2342720698336776_rDef)*delta_dEdX(j)
            dEdX3(j)            = dEdXnoBC3(j) + &
                                    (0.09956711350384978_rDef)*delta_dEdX(j)
            dEdX(iMax, j)       = dEdX1(j) 
            dEdX(iMax-1, j)     = dEdX2(j)
            dEdX(iMax-2, j)     = dEdX3(j)
        END DO
    ELSE IF ((DS == 6).OR.(DS == 7).OR.(DS == 8)) THEN  ! Compact Stencil
        DO j = bMin, bMax
            i = iMax
            dEdXnoBC1(j) = dEdX(i, j) 
        END DO
        i       = iMax
        rho_Zi   = -(-3.0_rDef*rho(i-0)    &
                    +4.0_rDef*rho(i-1)    &
                    -         rho(i-2))   &
                        /(2.0_rDef*delta_x)
        
        u_Zi     = -(-3.0_rDef*u(i-0)      &
                    +4.0_rDef*u(i-1)      &
                    -         u(i-2))     &
                        /(2.0_rDef*delta_x)

        p_Zi     = -(-3.0_rDef*p(i-0)      &
                    +4.0_rDef*p(i-1)      &
                    -         p(i-2))     &
                        /(2.0_rDef*delta_x)

        rho_TnoBC   = -1.0_rDef*(dZidX(i)*u(i)*rho_Zi               + &
                                rho(i)*dZidX(i)*u_Zi)
        u_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*u_Zi                 + &
                                (1.0_rDef/rho(i))*dZidX(i)*p_Zi)
        p_TnoBC     = -1.0_rDef*(u(i)*dZidX(i)*p_Zi                 + &
                                    gam*p(i)*dZidX(i)*u_Zi)

        INCLUDE 'IncApplyOutflow.f90'

        DO j = bMin, bMax
            dEdX_Thompson(j)    = -dQdT(j)
            delta_dEdX(j)       = dEdX_Thompson(j) - dEdXnoBC1(j)
            dEdX1(j)            = dEdXnoBC1(j) + delta_dEdX(j)
            dEdX(iMax, j)       = dEdX1(j)
        END DO
    END IF
