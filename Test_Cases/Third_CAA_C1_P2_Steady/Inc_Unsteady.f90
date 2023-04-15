IF (MOD(nT, 500) == 1) THEN
    WRITE(filename8, '(a9, a5, a1, a8, a1, a4, a1, a12, i7.7, a4)')       &
            'Unsteady_', LeftHandSide, '_', DifferStencil, '_',      &
            Dissipation, '_', LHS_Diss, nT, '.dat' 
    OPEN(57, FILE = TRIM(Current_Directory)//&
                &'/TimeStepSol/'&
                &//filename8, FORM = 'FORMATTED')
    DO i = iMin,   iMax
        WRITE(57, *)    x(i)                            ,&
                        !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                        Q_np1_kp1(i, 1)*Jac_Curv(i)           ,&
                        !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                        Q_np1_kp1(i, 2)/Q_np1_kp1(i, 1)             ,&
                        !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                        ((gam-1.0_rDef)*(Q_np1_kp1(i, 3)     &
                        - 0.5_rDef*Q_np1_kp1(i, 2)*Q_np1_kp1(i, 2)&
                        /Q_np1_kp1(i, 1)))*Jac_Curv(i)        ,&
                        !!!!!!!!!!!!Mach Number!!!!!!!!!!!
                        (Q_np1_kp1(i, 2)/Q_np1_kp1(i, 1))/(SQRT(gam*(&
                        ((gam-1.0_rDef)*(Q_np1_kp1(i, 3)     &
                        - 0.5_rDef*Q_np1_kp1(i, 2)*Q_np1_kp1(i, 2)&
                        /Q_np1_kp1(i, 1))))/(Q_np1_kp1(i, 1))))     ,& 
                        !!!!!!!!!!!!Exact Density!!!!!!!!!
                        Q_Exact(i, 1)*Jac_Curv(i)       ,&
                        !!!!!!!!!!!!Exact Velocity!!!!!!!!
                        Q_Exact(i, 2)/Q_Exact(i, 1)     ,&
                        !!!!!!!!!!!!Exact Pressure!!!!!!!!
                        ((gam-1.0_rDef)*(Q_Exact(i, 3) &
                        - 0.5_rDef*Q_Exact(i, 2)*&
                        Q_Exact(i, 2)/Q_Exact(i, 1)))*&
                        Jac_Curv(i)                     ,&
                        !!!!!!!!!!!!Exact Mach Number!!!!!
                        (Q_Exact(i, 2)/Q_Exact(i, 1))/&
                        (SQRT(gam*(((gam-1.0_rDef)*&
                        (Q_Exact(i, 3) - 0.5_rDef*&
                        Q_Exact(i, 2)*Q_Exact(i, 2)&
                        /Q_Exact(i, 1))))/(Q_Exact(i, 1)))), &
                        i, delta_xNew(i), Jac_Curv(i)
    END DO
    CLOSE(57)
END IF
