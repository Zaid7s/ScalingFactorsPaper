DO i = iMin, iMax
    rho     = Q_np1_k(i, 1)*Jac_Curv(i)
    u       = Q_np1_k(i, 2)/Q_np1_k(i, 1)
    p       = (gam-1.0_rDef)*(Q_np1_k(i, 3) - &
                0.5_rDef*Q_np1_k(i, 2)*Q_np1_k(i, 2)/&
                    Q_np1_k(i, 1))*Jac_Curv(i)              
    E_tot   = Q_np1_k(i, 3)*Jac_Curv(i)
    c       = SQRT(gam*p/rho)
    Mach    = u/c
    AreaFac = Source_Fac(i)/Area(i)

    Contra_U = dZidt(i) + dZidX(i)*u

    E(i, 1) = Q_np1_k(i, 2)*dZidX(i)

    E(i, 2) = (Q_np1_k(i, 2)*Q_np1_k(i, 2)/Q_np1_k(i, 1) + (gam-1.0_rDef)*&
                (Q_np1_k(i, 3) - Q_np1_k(i, 2)*Q_np1_k(i, 2)*0.5_rDef/&
                Q_np1_k(i, 1)))*dZidX(i)

    E(i, 3) = (Q_np1_k(i, 2)/Q_np1_k(i, 1))*&
                (Q_np1_k(i, 3) + (gam-1.0_rDef)*&
                (Q_np1_k(i, 3) - Q_np1_k(i, 2)*Q_np1_k(i, 2)*0.5_rDef/&
                Q_np1_k(i, 1)))*dZidX(i)

    S(i, 1) = AreaFac*Q_np1_k(i, 2)
    S(i, 2) = AreaFac*Q_np1_k(i, 2)*Q_np1_k(i, 2)/Q_np1_k(i, 1)
    S(i, 3) = AreaFac*(Q_np1_k(i, 2)/Q_np1_k(i, 1))*&
                (Q_np1_k(i, 3) + (gam-1.0_rDef)*&
                (Q_np1_k(i, 3) - Q_np1_k(i, 2)*Q_np1_k(i, 2)*0.5_rDef/&
                Q_np1_k(i, 1)))
END DO
