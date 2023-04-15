WRITE(filename4, '(a5, a1, a8, a1, a4, a6, i3.3, a12)')&
        LeftHandSide, '_', DifferStencil, '_',      &
        Dissipation, '_Step_', Step_Number,         &
        '_FFT_All.dat'
OPEN(71, FILE = TRIM(Current_Directory)//&
            &'/FFT/'&
            &//filename4, FORM = 'FORMATTED')
            
DO i = iMin, iMax
    CALL Real_FFT(  data    = Mean_Flow_FFT(:, i)   ,&
                    n       = Steps_Per_Cycle       ,&
                    iSign   = 1_i_def)
END DO

DO nC = 1, Steps_Per_Cycle
    WRITE(71, *) (Mean_Flow_FFT(nC, i), i = iMin, iMax)
END DO

DO nC = 1, Steps_Per_Cycle
    DO i = iMin, iMax
        All_FFT(nT, i) = Mean_Flow_FFT(nC, i)
    END DO
    WRITE(70, *) (All_FFT(nT, i), i = iMin, iMax)
END DO
