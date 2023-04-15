MODULE GetdQdT_BDF

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetGlobalParameter
    
    IMPLICIT NONE
    PRIVATE

    PUBLIC:: dQdT_BDF

    INTERFACE dQdT_BDF
        MODULE PROCEDURE dQdT_BDF_Mean
        MODULE PROCEDURE dQdT_BDF_Mean_Pet
    END INTERFACE dQdT_BDF

    CONTAINS
         
    SUBROUTINE dQdT_BDF_Mean(   iMin            ,&
                                iMax            ,&
                                bMin            ,&
                                bMax            ,&
                                Q_np1_k         ,&
                                Q_n             ,&
                                delta_t         ,&
                                dQdT            ,&
                                ndQdT)

    INTEGER, INTENT(IN)::   iMin        ,&
                            iMax        ,&
                            bMin        ,&
                            bMax        ,&
                            ndQdT

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::   Q_np1_k     ,&
                                                        Q_n

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::   dQdT 

    REAL(KIND = rDef), INTENT(IN):: delta_t

    INTEGER:: i, j

    DO j = bMin, bMax
        DO i = iMin, iMax 
            dQdT(i, j) = 0.0_rDef 
        END DO
    END DO

    DO j = bMin, bMax
        DO i = iMin, iMax
            dQdT(i, j) = (Q_np1_k(i, j) - Q_n(i, j))/(delta_t)
        END DO
    END DO

    END SUBROUTINE dQdT_BDF_Mean
         
    SUBROUTINE dQdT_BDF_Mean_Pet(   iMin            ,&
                                    iMax            ,&
                                    bMin            ,&
                                    bMax            ,&
                                    Q_np1_k         ,&
                                    Q_n             ,&
                                    Q_nm1           ,&
                                    Q_nm2           ,&
                                    delta_t         ,&
                                    dQdT            ,&
                                    ndQdT)

    INTEGER, INTENT(IN)::   iMin        ,&
                            iMax        ,&
                            bMin        ,&
                            bMax        ,&
                            ndQdT

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::   Q_np1_k     ,&
                                                        Q_n         ,&
                                                        Q_nm1       ,&
                                                        Q_nm2

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::   dQdT 

    REAL(KIND = rDef), INTENT(IN):: delta_t

    INTEGER:: i, j

    DO j = bMin, bMax
        DO i = iMin, iMax 
            dQdT(i, j) = 0.0_rDef 
        END DO
    END DO

    IF (ndQdT == 1) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax
                dQdT(i, j) = (Q_np1_k(i, j) - Q_n(i, j))/(delta_t)
            END DO
        END DO
    ELSEIF (ndQdT == 2) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax
                dQdT(i, j) = (  3.0_rDef*Q_np1_k(i, j)      -&
                                4.0_rDef*Q_n(i, j)          +&
                                1.0_rDef*Q_nm1(i, j))&
                                /(2.0_rDef*delta_t)
            END DO
        END DO
    ELSEIF (ndQdT == 3) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax
                dQdT(i, j) = (  10.0_rDef*Q_np1_k(i, j)      -&
                                15.0_rDef*Q_n(i, j)          +&
                                 6.0_rDef*Q_nm1(i, j)        -&
                                 1.0_rDef*Q_nm2(i, j))&
                                /(6.0_rDef*delta_t)
            END DO
        END DO
    END IF

    END SUBROUTINE dQdT_BDF_Mean_Pet
         
END MODULE GetdQdT_BDF
