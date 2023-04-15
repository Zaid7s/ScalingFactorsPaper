DO j = bMin, bMax
    DO i = iMin, iMax
        Q_Store(nT, i, j) = Q_n(i, j)
    END DO
END DO

DO j = bMin, bMax
    DO i = iMin, iMax
        Q_np1_k(i, j)   = Q_Store(nT, i, j)
    END DO
END DO

IF (nT == 1) THEN
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_nm1(i, j)   = Q_Store(nT, i, j)
            Q_nm2(i, j)   = Q_Store(nT, i, j)
        END DO
    END DO
ELSE IF (nT == 2) THEN
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_nm1(i, j)   = Q_Store(nT-1, i, j)
            Q_nm2(i, j)   = Q_Store(nT-1, i, j)
        END DO
    END DO
ELSE
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_nm1(i, j)   = Q_Store(nT-1, i, j)
            Q_nm2(i, j)   = Q_Store(nT-2, i, j)
        END DO
    END DO
END IF
