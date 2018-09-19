      SUBROUTINE MKDSTRB(NPTS, R0, R1, IDSTR, NP, P, X, Y)
C     Создаем вектор отсчетов и значений функции распределения на 
C     интервале R0..R1
C     тип распределения задается IDSTR
C     NPTS      - число отсчетов
C     P(NP)      - вектор с параметрами распределения
C     R0,R1     - левая и правая границы интервалов размеров частиц
C     IDSTR     - тип распределения:
C                 1 - степенное P(1) - gamma
C                 2 - логнормальное P(1) - r0, P(2) - sigma
      IMPLICIT NONE
      INTEGER NPTS, IDSTR, NP
      REAL*8  R0, R1, X(NPTS), Y(NPTS), P(NP)
CF2PY INTENT(IN)  NPTS, R0, R1, IDSTR, P
CF2PY INTENT(OUT) X, Y
      REAL*8  DR
      INTEGER I

      DR = (R1 - R0) / DBLE(NPTS - 1)
      
      DO I=1, NPTS
        X(I) = DBLE(I-1)*DR+R0
      END DO

      IF (IDSTR .EQ. 1) THEN
C Степенное распределение
        DO I=1, NPTS
          Y(I) = (X(I)/X(1))**P(1)
        END DO
     
      ELSE IF (IDSTR .EQ. 2) THEN
C логнормальное распределение
      
      ELSE
        PRINT *, 'asa'
        STOP
      END IF

      RETURN
      
      END