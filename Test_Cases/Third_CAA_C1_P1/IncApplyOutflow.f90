A_S_noBC        = p_TnoBC -   Mean_c*Mean_c*rho_TnoBC
A_Plus_noBC     = p_TnoBC+Mean_rho*Mean_c*u_TnoBC
A_Minus_noBC    = p_TnoBC-Mean_rho*Mean_c*u_TnoBC

A_S_BC      = A_S_noBC 
A_Plus_BC   = A_Plus_noBC
A_Minus_BC  = 2.0_rDef*((P_Static-P_Outflow)/(delta_t)) - A_Plus_BC    !LU-SGS

p_T     = 0.5_rDef*(A_Plus_BC+A_Minus_BC)*Inv_Jac_Curv(iMax)
u_T     = ((A_Plus_BC-A_Minus_BC)/(2.0_rDef*Mean_rho*Mean_c))&
                    *Inv_Jac_Curv(iMax)
rho_T   = (((A_Plus_BC+A_Minus_BC) - 2.0_rDef*A_S_BC)/&
                        (2.0_rDef*Mean_c*Mean_c))*Inv_Jac_Curv(iMax)

dQdT(1) = rho_T 
dQdT(2) = rho(iMax)*u_t+u(iMax)*rho_t 
dQdT(3) = gm1*p_T+0.5_rDef*(rho(iMax)*u(iMax)*u_T + &
                             rho(iMax)*u(iMax)*u_T + &
                               u(iMax)*u(iMax)*rho_t)
