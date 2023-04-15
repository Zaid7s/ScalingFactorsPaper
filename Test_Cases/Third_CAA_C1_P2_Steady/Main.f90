PROGRAM Main
    !! Hi
    USE, INTRINSIC:: ISO_FORTRAN_ENV

    USE GetMeanValues

    IMPLICIT NONE
    INTEGER, PARAMETER:: rDef = REAL64

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_Mean              ,&
                                                        Q_Exact             ,&
                                                        Q_Mean_Store        ,&
                                                        Q_Mean_And_Pet_Store, &
                                                        Q_Mean_And_Pet      ,&
                                                        Q_Peturbation

    REAL(KIND = rDef):: p_Exit  ,&
                        gam

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE:: MaxDisl2Res_All
    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE:: MaxDisl2Res_Table

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE::  x           ,&
                                                    Jac_Curv

    CHARACTER(LEN = 8)  :: DifferStencil

    CHARACTER(LEN = 28) ::  filename1         ,&
                            filename2

    INTEGER::   nFlow           ,&
                i               ,&
                j               ,&
                iMin            ,&
                iMax            ,&
                bMin            ,&
                bMax            ,&
                nL              ,&
                nD              ,&
                DS              ,&
                Steps_Per_Cycle, &
                nSPC            ,&
                nG              ,&
                MaxGrid         ,&
                MaxStepPerCycle

    INTEGER, DIMENSION(:), ALLOCATABLE:: SubIterPerCyc_All
    INTEGER, DIMENSION(:, :), ALLOCATABLE:: SubIterPerCyc_Table

    LOGICAL:: debug

    MaxGrid         = 7
    MaxStepPerCycle = 5

    ALLOCATE(   MaxDisl2Res_All(MaxGrid)                        ,&
                MaxDisl2Res_Table(MaxStepPerCycle, MaxGrid)     ,&
                SubIterPerCyc_All(MaxGrid)                      ,&
                SubIterPerCyc_Table(MaxStepPerCycle, MaxGrid))

    debug = .TRUE.
    debug = .FALSE.

    DO nG = 9, 9  ! Grid 
        DO DS = 1, 8  ! 4 Max, Differencing Stencil 
            nFlow = 1
            CALL MeanValues(Q_Mean  = Q_Mean    ,&
                            Q_Exact = Q_Exact   ,&
                            nFlow   = nFlow     ,&
                            p_Exit  = p_Exit    ,&
                            nL      = nL        ,&
                            nD      = nD        ,&
                            DS      = DS        ,&
                            nG      = nG)
            DEALLOCATE(Q_Mean, Q_Exact)
        END DO
    END DO
END PROGRAM Main
