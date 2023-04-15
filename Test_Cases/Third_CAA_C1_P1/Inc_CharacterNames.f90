IF (DS == 1) THEN
    DifferStencil   = 'SecOrder'
ELSE IF (DS == 2) THEN
    DifferStencil   = 'FouOrder'
ELSE IF (DS == 3) THEN
    DifferStencil   = 'SixOrder'
ELSE IF (DS == 4) THEN
    DifferStencil   = 'RDRPSten'
ELSE IF (DS == 5) THEN
    DifferStencil   = '_DRPSten'
ELSE IF (DS == 6) THEN
    DifferStencil   = 'Compact4'
ELSE IF (DS == 7) THEN
    DifferStencil   = 'PreFac_4'
ELSE IF (DS == 8) THEN
    DifferStencil   = 'PreFac_6'
END IF

IF (nD == 1) THEN
    Dissipation = 'D_02'
ELSE IF (nD == 2) THEN
    Dissipation = 'D_04'
ELSE IF (nD == 3) THEN
    Dissipation = 'D_06'
ELSE IF (nD == 4) THEN
    Dissipation = 'D_08'
ELSE IF (nD == 5) THEN
    Dissipation = 'D_10'
END IF

IF (nDis == 1) THEN            ! First Order Backward Difference
    LHS_Diss = 'LHS_Diss_Off'
ELSEIF (nDis == 2) THEN        ! Second Order Backward Difference
    LHS_Diss = 'LHS_Diss_On_'
END IF

IF (nL == 1) THEN
    LeftHandSide = "LUSGS"
ELSE IF (nL == 2) THEN
    LeftHandSide = "TriDi"
END IF
