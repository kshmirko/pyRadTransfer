      SUBROUTINE READ_MEAS_FILE(FNAME, MXMEAS, ANG, IV, QV, N)
      IMPLICIT NONE
CF2PY INTETN(IN)    MXMEAS, FNAME
CF2PY INTENT(OUT)   ANG, IV, QV, N
      INTEGER       MXMEAS, N, I, IERR
      CHARACTER*64  FNAME
      REAL*8        ANG(MXMEAS), IV(MXMEAS), QV(MXMEAS)

      PRINT *, FNAME, MXMEAS
      OPEN(UNIT=150, FILE=FNAME, STATUS='OLD', IOSTAT=IERR)

      IF (IERR .NE. 0) THEN
        PRINT *, "Error"
        STOP
      ENDIF

      I=1
      DO WHILE ((IERR .EQ. 0).AND.(I.LT.MXMEAS))

        READ(150, *, IOSTAT=IERR) ANG(I), IV(I), QV(I)
C        PRINT *, ANG(I), IV(I), QV(I), IERR
        I=I+1        
      
      ENDDO

      N = I-2

      CLOSE(150)

      RETURN

      END

