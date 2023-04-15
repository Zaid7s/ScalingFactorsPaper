A_S_BC         = 0.0_rDef 
A_Plus_BC      = 0.0_rDef 
A_Minus_BC     = 0.0_rDef 
A_S_noBC       = 0.0_rDef 
A_Plus_noBC    = 0.0_rDef 
A_Minus_noBC   = 0.0_rDef 
p_T            = 0.0_rDef 
rho_T          = 0.0_rDef 
omega          = 0.0_rDef 
epsi           = 0.0_rDef 
u_T            = 0.0_rDef 
p_TnoBC        = 0.0_rDef 
rho_TnoBC      = 0.0_rDef 
u_TnoBC        = 0.0_rDef 
E_tot_TnoBC    = 0.0_rDef 
p_Zi           = 0.0_rDef 
E_tot_Zi       = 0.0_rDef 
rho_Zi         = 0.0_rDef 
u_Zi           = 0.0_rDef 
Delta_A_S      = 0.0_rDef 
Delta_A_Plus   = 0.0_rDef 
Mach           = 0.0_rDef 
P_Inflow       = 0.0_rDef 
rho_Inflow     = 0.0_rDef 
Bar_u          = 0.0_rDef 
Bar_p          = 0.0_rDef 
Bar_rho        = 0.0_rDef 

DO j = bMin, bMax
    dEdXnoBC1(j)        =   0.0_rDef 
    dQdTnoBC(j)         =   0.0_rDef
    dQdT(j)             =   0.0_rDef
    dEdX1(j)            =   0.0_rDef
    dEdX2(j)            =   0.0_rDef
    dEdX3(j)            =   0.0_rDef
    dEdXnoBC2(j)        =   0.0_rDef
    dEdXnoBC3(j)        =   0.0_rDef
    delta_dEdX(j)       =   0.0_rDef
    dEdX_Thompson(j)    =   0.0_rDef
END DO

DO i = iMin, iMax
    E_tot(i)            =   0.0_rDef
    p(i)                =   0.0_rDef
    rho(i)              =   0.0_rDef
    u(i)                =   0.0_rDef
END DO
