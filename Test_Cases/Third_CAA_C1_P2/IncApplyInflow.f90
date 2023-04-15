A_S_noBC        = p_TnoBC -   Mean_c*Mean_c*rho_TnoBC
A_Plus_noBC     = p_TnoBC+Mean_rho*Mean_c*u_TnoBC
A_Minus_noBC    = p_TnoBC-Mean_rho*Mean_c*u_TnoBC

A_S_BC      = A_S_noBC      + Delta_A_S
A_Plus_BC   = A_Plus_noBC   + Delta_A_Plus
A_Minus_BC  = A_Minus_noBC 

p_T     = 0.5_rDef*(A_Plus_BC+A_Minus_BC)*Inv_Jac_Curv(iMin)
u_T     = ((A_Plus_BC-A_Minus_BC)/(2.0_rDef*Mean_rho*Mean_c))&
                    *Inv_Jac_Curv(iMin)
rho_T   = (((A_Plus_BC+A_Minus_BC) - 2.0_rDef*A_S_BC)/&
                        (2.0_rDef*Mean_c*Mean_c))*Inv_Jac_Curv(iMin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dQdT(1) = rho_T 
dQdT(2) = rho(iMin)*u_t+u(iMin)*rho_t 
dQdT(3) = gm1*p_T+0.5_rDef*(2.0_rDef*rho(iMin)*u(iMin)*u_T + &
                              u(iMin)*u(iMin)*rho_t)

IF (debug) WRITE(0, *) -dQdT(:)
IF (debug) WRITE(0, *) dEdXnoBC1(:)
