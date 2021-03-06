C     EVALUATE LEGENDRE POLYNOMIAL AT SPECIFIC X
C     COEF - LEGENDRE COEFFICIENTS, I.E. (2*i+1)Pm
      FUNCTION LEGEVAL(X, N, COEF) 
      IMPLICIT NONE
      INTEGER N
      REAL*8  X, COEF(0:N), LEGEVAL
      
      INTEGER J
      REAL*8  P0, P1, P2, SUM
      
      P0 = 1.0D0
      P1 = X
      LEGEVAL = COEF(0)*P0+COEF(1)*P1
      
      DO J=2, N
         P2 = DBLE(2*J-1) / DBLE(J) * X * P1 - DBLE(J-1) / DBLE(J) * P0     
         LEGEVAL = LEGEVAL + COEF(J) * P2
         P0 = P1
         P1 = P2
      END DO
      
      RETURN
      END

C     ------------------------------------------------------------------
C     SAME AS ABOVE, BUT X IS A VECTOR
C     EVALUATE LEGENDTE POLYNOMIAL FOR GIVEN POINTS X(N) AND PUT THE
C     RESULT INTO Y(N)
C     ------------------------------------------------------------------
      SUBROUTINE LEGVAL_A(M, COEF, N, X, Y)
      IMPLICIT NONE
      
      INTEGER N, M
      REAL*8  X(N), Y(N), COEF(0:M), LEGEVAL
      EXTERNAL LEGEVAL
      
      INTEGER I
      
      DO 10 I=1, N
         Y(I) = LEGEVAL(X(I), M, COEF)
 10   CONTINUE
      
      RETURN
      
      END
      
