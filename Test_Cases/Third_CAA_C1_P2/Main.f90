PROGRAM Main

    USE, INTRINSIC:: ISO_FORTRAN_ENV

    USE GetMeanValues
    USE GetNormalizeMean
    USE GetMean_And_Pet
    USE GetGlobalParameter  !   Defines floating point precision 
                            !   and other parameter
    IMPLICIT NONE

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_Mean              ,&
                                                        Q_Exact             ,&
                                                        Q_Mean_And_Pet      ,&
                                                        Q_Peturbation

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE :: Mean_FFT1

    CHARACTER(LEN = 28) ::  filename1         ,&
                            filename2

    INTEGER::   Steps_Per_Cycle     ,&
                nG                  ,&
                nL                  ,&
                DS                  ,&
                nDis                ,&
                nD                  ,&
                nI                  ,&
                nBound

    LOGICAL:: debug

    REAL(KIND = rDef) ::    kDeltaX        ,&
                            Scaling_Fac

    CHARACTER(LEN = 8)  :: DifferStencil
    CHARACTER(LEN = 4)  :: Dissipation
    CHARACTER(LEN = 5)  :: LeftHandSide
    CHARACTER(LEN = 12) :: LHS_Diss

    nI    = 1 !If 0, we calculate the MeanFlow from Scratch
              !If 1, caluclate MeanFlow from Exact solution
              !1 is much more efficient

    Steps_Per_Cycle     = 128

    debug = .TRUE.
    debug = .FALSE.

    DO DS = 1, 6
        !! Different RHS Differncing Stencil.   1--> Second Order
        !!                                      2--> Fourth Order
        !!                                      3--> Sixth Order
        !!                                      4--> RDRP
        !!                                      5--> DRP
        DO nL = 1, 1
        !! Different LHS. 1--> Uses Lower upper Symmetric Gauss Seidel Method by 
        !!                      Yoon and Jameson
        !!                2--> Uses ADI on the LHS
            DO nDis = 1, 1
                !! Different Implicit Dissipation.  1 --> No Dissipation
                !!                                  2 --> Second Order
                DO nD = 5, 5
                    IF (DS == 1) THEN
                        Scaling_Fac     = 1.0_rDef 
                        kDeltaX         = 1.0_rDef
                    ELSE IF (DS == 2) THEN
                        Scaling_Fac     = 2.12179602298_rDef
                        kDeltaX         = 1.400_rDef
                    ELSE IF (DS == 3) THEN
                        Scaling_Fac     = 6.352796985_rDef
                        kDeltaX         = 1.586_rDef
                    ELSE IF (DS == 4) THEN
                        Scaling_Fac     = 6.36206870932_rDef
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 5) THEN
                        Scaling_Fac     = 4.97487147161_rDef
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 6) THEN
                        Scaling_Fac     = 3.800180_rDef
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 7) THEN
                        CYCLE
                        Scaling_Fac     = 3.800180_rDef
                        kDeltaX         = 1.664_rDef
                    ELSE IF (DS == 8) THEN
                        Scaling_Fac     = 10.04379_rDef
                        kDeltaX         = 1.664_rDef
                    ELSE
                        EXIT
                    END IF

                    INCLUDE 'Inc_CharacterNames.f90'

                    IF ((nL == 1).AND.(nDis == 2)) THEN
                        CYCLE
                    END IF

                    DO nG = 9, 9  ! Grid 
                        CALL MeanValues(Q_Mean              = Q_Mean        ,&
                                        Q_Exact             = Q_Exact       ,&
                                        nL                  = nL            ,&
                                        nD                  = nD            ,&
                                        DS                  = DS            ,&
                                        nDis                = nDis          ,&
                                        kDeltaX             = kDeltaX       ,&
                                        Scaling_Fac         = Scaling_Fac   ,&
                                        nInitalCondition    = nI            ,&
                                        nG                  = nG)

                        nBound = 0
                        CALL NormalizeMean( Q_Mean_And_Pet      = Q_Mean_And_Pet    ,&
                                            Q_Mean              = Q_Mean            ,&
                                            Q_Exact             = Q_Exact           ,&
                                            Mean_FFT1           = Mean_FFT1         ,&
                                            Steps_Per_Cycle     = Steps_Per_Cycle   ,&
                                            Scaling_Fac         = Scaling_Fac       ,&
                                            kDeltaX             = kDeltaX           ,&
                                            nL                  = nL                ,&
                                            nD                  = nD                ,&
                                            DS                  = DS                ,&
                                            nDis                = nDis              ,&
                                            nG                  = nG                ,&
                                            nBound              = nBound)

                        nBound = 1
                        CALL Mean_And_Pet(  Q_Peturbation       = Q_Peturbation     ,&
                                            Q_Mean              = Q_Mean            ,&
                                            Q_Exact             = Q_Exact           ,&
                                            Mean_FFT1           = Mean_FFT1         ,&
                                            Steps_Per_Cycle     = Steps_Per_Cycle   ,&
                                            Scaling_Fac         = Scaling_Fac       ,&
                                            kDeltaX             = kDeltaX           ,&
                                            nL                  = nL                ,&
                                            nD                  = nD                ,&
                                            DS                  = DS                ,&
                                            nDis                = nDis              ,&
                                            nG                  = nG                ,&
                                            nBound              = nBound)
                    END DO !nG Grid
                    DEALLOCATE(Q_Mean           ,& 
                               Q_Exact          ,&
                               Q_Mean_And_Pet   ,&
                               Q_Peturbation)
                END DO !nD  Explicit Dissipation
            END DO !nDis Implicit Dissipation
        END DO !nL LHS
    END DO !DS Difference
END PROGRAM Main
