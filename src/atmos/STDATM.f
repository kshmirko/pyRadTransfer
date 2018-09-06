      SUBROUTINE STDATM1(Z, T0, P0, T, P, RHO)
      IMPLICIT NONE

      REAL*8  Z, T, P, RHO, T0, P0
CF2PY INTENT(IN)  Z, T0, P0
CF2PY INTENT(OUT) T, P, RHO

      IF (Z .LT.  11000.D0) THEN
        T = 15.04-0.00649*Z+273.15D0
        P = 101.29*(T/288.08D0)**5.256
      ELSE IF (Z .LT. 25000.D0) THEN
        T = -56.46+273.15D0
        P = 22.65 * EXP(1.73D0-0.000157*Z)
      ELSE 
        T = -131.21+0.00299*Z + 273.15
        P = 2.488 * (T/216.6)**(-11.388D0)
      END IF
      
      
      RHO = P/(0.2869*T)
      P = P *1000.08D0
C      RHO = (P/T)*(T0/P0)

      RETURN
      END

      SUBROUTINE STDATM2(Z, T0, P0, T, P, RHO, N)
      IMPLICIT NONE
CF2PY INTENT(IN)  Z, T0, P0, N
CF2PY INTENT(OUT) T, P, RHO
      INTEGER N, I
      REAL*8  T(N), P(N), RHO(N), Z(N), T0, P0

      EXTERNAL  STDATM1

      DO I=1, N

        CALL STDATM1(Z(I), T0, P0, T(I), P(I), RHO(I))

      END DO

      RETURN

      END