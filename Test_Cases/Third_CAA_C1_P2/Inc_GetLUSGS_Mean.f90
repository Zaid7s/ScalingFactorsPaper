INTEGER, INTENT(IN)::   iMin     ,&
                        iMax     ,&
                        bMax     ,&
                        bMin     ,&
                        DS       ,&
                        nD       ,&
                        nStages  ,&
                        nI       ,&
                        nDis

INTEGER, INTENT(INOUT):: nTau

REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                 delta_t    ,&
                                 delta_tau  ,&
                                 time

REAL(KIND = rDef), INTENT(INOUT):: Mach        ,&
                                    u           ,&
                                    p           ,&
                                    c           ,&
                                    gm1         ,&
                                    rho         ,&
                                    rho_Static  ,&
                                    P_Static

REAL(KIND = rDef), INTENT(INOUT):: Newtonl2Res

REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  bet             ,&
                                                Source_Fac      ,&
                                                Jac_Curv        ,&
                                                Inv_Jac_Curv    ,&
                                                dZidX           ,&
                                                dZidt           ,&
                                                Area
REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS             ,&
                                                        Q_np1_k         ,&
                                                        Q_n             ,&
                                                        Q_Exact         ,&
                                                        dQdTau          ,&
                                                        Q_np1           ,&
                                                        Q_np1_kp1       ,&
                                                        Dis             ,&
                                                        Delta_Q_star2   ,&
                                                        Delta_Q_star    ,&
                                                        Delta_Q         ,&
                                                        AA              ,&
                                                        BB              ,&
                                                        EE

REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(INOUT)::  L_LowerDiag     ,&
                                                        L_Diag          ,&
                                                        U_UpperDiag     ,&
                                                        U_Diag          ,&
                                                        D_Diag          ,&
                                                        A_Minus         ,&
                                                        A_Plus          ,&
                                                        S_Jac           ,&
                                                        S_Minus         ,&
                                                        S_Plus

REAL(KIND = rDef) ::    DissiFac        ,&
                        l2Fac           ,&
                        A1_Inflow       ,&
                        A2_Inflow       ,&
                        A3_Inflow       ,&
                        A1_Outflow      ,&
                        A2_Outflow      ,&
                        A3_Outflow

REAL(KIND = rDef), DIMENSION(bMax, bMax, iMax) :: Tau


REAL(KIND = rDef), DIMENSION(iMax, bMax):: dQdT        ,&
                                            dEdX        ,&
                                            E           ,&
                                            RHS2

INTEGER            ::   i               ,&
                        k               ,&
                        j               ,&
                        l               ,&
                        m               ,&
                        n
LOGICAL:: debug

debug = .TRUE.
debug = .FALSE.

! The system of equations solved here is in the form of Ax = B. Using the LUSGS 
! we break down the A Matrix to three matrices:
!                                       - L      Matrix 
!                                       - (D^-1) Matrix
!                                       - U      Matrix
! The B Vector is the RHS. The X Vector is Delta_Q. The solution for Delta_Q
! is obtained over three steps:
!               (1)  L{Delta_Q^**}    = RHS
!               (2)  D^(-1){Delta_Q*} = {Delta_Q^**}
!               (3)  U{Delta_Q}       = {Delta_Q^*}

! Step (1) is easily solved by a Forward Sweep. 
! Step (2) is solved by a direct multiplication since :
!                                       {{D^-1}^-1}{Delta_Q^**} == D{Delta_Q^**}
! Step (3) is solved by a Backward Sweep

! This subroutine creates the L Matrix 
    CALL LDU_Vectors(   iMin           = iMin           ,&
                        iMax           = iMax           ,&
                        bMin           = bMin           ,&
                        bMax           = bMax           ,&
                        delta_x        = delta_x        ,&
                        delta_t        = delta_tau      ,&
                        A_Plus         = A_Plus         ,&
                        A_Minus        = A_Minus        ,&
                        S_Jac          = S_Jac          ,&
                        Source_Fac     = Source_Fac     ,&
                        L_LowerDiag    = L_LowerDiag    ,&
                        U_UpperDiag    = U_UpperDiag    ,&
                        Tau            = Tau            ,&
                        Jac_Curv       = Jac_Curv       ,&
                        Inv_Jac_Curv   = Inv_Jac_Curv   ,&
                        dZidX          = dZidX          ,&
                        dZidt          = dZidt)
    
! Get the RHS. This includes the boundaries and dissipation     
    CALL RHSVector( iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    Q_n             = Q_n           ,&
                    Q_np1_k         = Q_np1_k       ,&
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
                    nDis            = nDis           ,&
                    A3_Outflow      = A3_Outflow)

    CALL SolveLUSGS(L_LowerDiag     =  L_LowerDiag      ,&    
                    U_UpperDiag     =  U_UpperDiag      ,& 
                    Tau             =  Tau              ,&
                    iMin            =  iMin             ,& 
                    iMax            =  iMax             ,& 
                    bMin            =  bMin             ,& 
                    bMax            =  bMax             ,& 
                    RHS             =  RHS              ,& 
                    Delta_Q_star    =  Delta_Q_star     ,&
                    Delta_Q_star2   =  Delta_Q_star2    ,&
                    Delta_Q         =  Delta_Q) 

    
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_np1_kp1(i, j) = Delta_Q(i, j) + Q_np1_k(i, j)
        END DO
    END DO
