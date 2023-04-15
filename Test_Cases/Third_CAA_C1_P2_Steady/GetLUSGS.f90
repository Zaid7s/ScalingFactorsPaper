MODULE GetLUSGS
!!! Working Test Case
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetRHSVector
    USE GetRungeKutta
    USE GetLMatrix
    USE GetUMatrix
    USE GetDMatrix
    USE GetSolveLUSGS
    USE GetLDU_Vectors
    USE ForwardSweep
    USE BackwardSweep

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: LUSGS

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE LUSGS( iMin          ,&
                      iMax          ,&
                      bMin          ,&
                      bMax          ,&
                      Q_n           ,&
                      Q_np1_k       ,&
                      Q_np1         ,&
                      Q_nm1         ,&
                      Q_nm2         ,&
                      Q_np1_kp1     ,&
                      RHS           ,&
                      time          ,&
                      dQdTau        ,&
                      delta_x       ,&
                      delta_Zi      ,&
                      delta_tau     ,&
                      delta_t       ,&
                      DS            ,&
                      nD            ,&
                      Dis           ,&
                      L_LowerDiag   ,&
                      U_UpperDiag   ,&
                      A_Plus        ,&
                      A_Minus       ,&
                      S_Jac         ,&
                      Delta_Q_star2, &
                      AA            ,&
                      BB            ,&
                      EE            ,&
                      Delta_Q_star  ,&
                      Delta_Q       ,&
                      bet           ,&
                      nTau          ,&
                      nStages       ,&
                      nI            ,&
                      Newtonl2Res   ,&
                      Source_Fac    ,&
                      gam           ,&
                      Mach          ,&
                      u             ,&
                      p             ,&
                      c             ,&
                      gm1           ,&
                      rho           ,&
                      rho_Static    ,&
                      P_Static      ,&
                      p0In          ,&
                      rho0In        ,&
                      Q_Exact       ,&
                      Jac_Curv      ,&
                      Inv_Jac_Curv  ,&
                      dZidX         ,&
                      dZidt         ,&
                      nFlow         ,&
                      Area          ,&    
                      A1_Inflow     ,&
                      A2_Inflow     ,&
                      A3_Inflow     ,&
                      A1_Outflow    ,&
                      A2_Outflow    ,&
                      nDis          ,&
                      A3_Outflow) 


    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS       ,&
                           nD       ,&
                           nStages  ,&
                           nI       ,&
                           nFlow    ,&
                           nDis
    
    INTEGER, INTENT(INOUT):: nTau

    REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                     delta_t    ,&
                                     delta_Zi   ,&
                                     delta_tau  ,&
                                     time

    REAL(KIND = rDef), INTENT(INOUT):: gam          ,&
                                        Mach        ,&
                                        u           ,&
                                        p           ,&
                                        c           ,&
                                        gm1         ,&
                                        rho         ,&
                                        rho_Static  ,&
                                        P_Static    ,&
                                        p0In        ,&
                                        rho0In

    REAL(KIND = rDef), INTENT(INOUT):: Newtonl2Res
    
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
    
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS             ,&
                                                            Q_np1_k         ,&
                                                            Q_n             ,&
                                                            Q_Exact         ,&
                                                            dQdTau          ,&
                                                            Q_np1           ,&
                                                            Q_nm1           ,&
                                                            Q_nm2           ,&
                                                            Q_np1_kp1       ,&
                                                            Dis             ,&
                                                            Delta_Q_star2   ,&
                                                            Delta_Q_star    ,&
                                                            Delta_Q         ,&
                                                            AA              ,&
                                                            BB              ,&
                                                            EE

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(INOUT)::  L_LowerDiag     ,&
                                                            U_UpperDiag     ,&
                                                            A_Minus         ,&
                                                            A_Plus          ,&
                                                            S_Jac
    
    REAL(KIND = rDef):: DissiFac        ,&
                        l2Fac           ,&
                        comp_time1      ,&
                        comp_time2      ,&
                        comp_time3      ,&
                        start_time1     ,&
                        start_time2     ,&
                        start_time3     ,&
                        start_time4     ,&
                        start_time5     ,&
                        start_time6     ,&
                        start_time7     ,&
                        end_time1       ,&
                        end_time2       ,&
                        end_time3       ,&
                        end_time4       ,&
                        end_time5       ,&
                        end_time6       ,&
                        end_time7

    REAL(KIND = rDef), DIMENSION(bMax, bMax, iMax) :: Tau

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
    CALL CPU_TIME(start_time1)
    !!!!!!!INCLUDE 'LUSGS_Diagonals.f90'

    CALL LDU_Vectors(   iMin           = iMin           ,&
                        iMax           = iMax           ,&
                        bMin           = bMin           ,&
                        bMax           = bMax           ,&
                        delta_x        = delta_x        ,&
                        delta_Zi       = delta_Zi       ,&
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

    CALL CPU_TIME(end_time1)

! Get the RHS. This includes the boundaries and dissipation     
    CALL CPU_TIME(start_time2)
    CALL RHSVector( iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    Q_n             = Q_n           ,&
                    Q_np1_k         = Q_np1_k       ,&
                    Q_nm1           = Q_nm1         ,&
                    Q_nm2           = Q_nm2         ,&
                    dQdTau          = dQdTau        ,&
                    RHS             = RHS           ,&
                    time            = time          ,&
                    delta_x         = delta_x       ,&
                    delta_Zi        = delta_Zi      ,&
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
                    gam             = gam           ,&
                    Mach            = Mach          ,&
                    u               = u             ,&
                    p               = p             ,&
                    c               = c             ,&
                    gm1             = gm1           ,&
                    rho             = rho           ,&
                    rho_Static      = rho_Static    ,&
                    P_Static        = P_Static      ,&
                    p0In            = p0In          ,&
                    rho0In          = rho0In        ,&
                    Q_Exact         = Q_Exact       ,&
                    Jac_Curv        = Jac_Curv      ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv  ,&
                    dZidX           = dZidX         ,&
                    dZidt           = dZidt         ,&
                    nFlow           = nFlow         ,&
                    Area            = Area          ,&
                    A1_Inflow       = A1_Inflow     ,&
                    A2_Inflow       = A2_Inflow     ,&
                    A3_Inflow       = A3_Inflow     ,&
                    A1_Outflow      = A1_Outflow    ,&
                    A2_Outflow      = A2_Outflow    ,&
                    nDis            = nDis           ,&
                    A3_Outflow      = A3_Outflow)
    CALL CPU_TIME(end_time2)

    CALL CPU_TIME(start_time3)
    !!! INCLUDE 'Sweeping_3Steps.f90'
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

    CALL CPU_TIME(end_time3)
    
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_np1_kp1(i, j) = Delta_Q(i, j) + Q_np1_k(i, j)
        END DO
    END DO
        
    comp_time1 = ABS(end_time1 - start_time1)
    comp_time2 = ABS(end_time2 - start_time2)
    comp_time3 = ABS(end_time3 - start_time3)


    END SUBROUTINE LUSGS

END MODULE GetLUSGS
