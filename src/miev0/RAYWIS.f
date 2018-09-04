      SUBROUTINE RAYWIS(WIS, M)
      IMPLICIT NONE
      
CF2PY INTENT(IN)  M
CF2PY INTENT(OUT) WIS
      INTEGER M
      REAL*8  WIS(0:M,4)
      
      WIS = 0.0D0
      WIS(0:2,1) = (/0.75, 0., 0./)
      WIS(0:2,2) = (/0.25, 0., 0.1/)
      WIS(0:2,3) = (/0.0, 0.25  , 0./)
      END
      