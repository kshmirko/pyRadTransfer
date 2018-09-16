      SUBROUTINE RAYEV(EV, M)
C матрица эванса в коэффициентах лежандра
      IMPLICIT NONE
      
CF2PY INTENT(IN)  M
CF2PY INTENT(OUT) EV
      INTEGER M
      REAL*8  EV(0:M,6)
      
      EV = 0.0D0
      EV(0:2,1) = (/1.0, 0.0, 0.5/)
      EV(0:2,2) = (/-0.5, 0.0, 0.5/)
      EV(0:2,3) = (/0.0, 1.5  , 0./)
      EV(0:2,5) = EV(0:2,1)
      EV(0:2,6) = EV(0:2,3)
      END