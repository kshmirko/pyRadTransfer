      SUBROUTINE MAKESCATFILE(FNAME, MIDX, R0, R1, GAMMA, NPTS, WL,
     .                        TAUA, TAUM, NMOMS, TAUE_T, TAUS_T, SSA_T, 
     .                        EV_T)
CF2PY INTENT(IN) FNAME, MIDX, R0, R1, GAMMA, TAUA, TAUM, NMOMS
CF2PY INTENT(IN) NPTS
CF2PY INTENT(OUT) TAUE_T, TAUS_T, SSA_T, EV_T
      CHARACTER*64  FNAME
      COMPLEX*16    MIDX
      REAL*8        R0, R1, GAMMA, WL, TAUM, TAUA, DR
      REAL*8        TAUE_T, TAUS_T, SSA_T, EV_T(0:NMOMS, 6), VOL
      REAL*8        EV_M(0:NMOMS, 6), WIS_A(0:NMOMS, 4), EXT, SCA, ASY
      REAL*8        X(NPTS), Y(NPTS), EV_A(0:NMOMS, 6)
      INTEGER       NMOMS, NPTS, IERR, I
      REAL*8        OMEGA_A, TAUS_A

      EXTERNAL RAYEV, DISTR, WISCOMBE2EVANS

      
      TAUE_T = TAUA+TAUM

      DR = (R1-R0)/DBLE(NPTS-1)
      
      DO I=1, NPTS
        X(I) = DBLE(I-1)*DR+R0
        Y(I) = (X(I)/R0)**GAMMA
      END DO
      
      
      CALL DISTR(NPTS, X, Y, MIDX, WL, NMOMS, WIS_A, EXT, SCA, ASY, 
     .                VOL, IERR)

      
      IF (IERR /= 0) THEN
        STOP
      END IF

      CALL WISCOMBE2EVANS(NMOMS, WIS_A, EV_A)
      CALL RAYEV(EV_M, NMOMS)

      OMEGA_A = SCA/EXT
      TAUS_A = TAUA * OMEGA_A

      TAUS_T = TAUS_A + TAUM
      SSA_T = TAUS_T/TAUE_T

      
      EV_T = (EV_A*TAUS_A+EV_M*TAUM)/TAUS_T
      
      OPEN(UNIT=100, FILE=FNAME, STATUS='UNKNOWN')

      WRITE(100, '(E12.4)') TAUE_T
      WRITE(100, '(E12.4)') TAUS_T
      WRITE(100, '(E12.4)') SSA_T
      WRITE(100, '(I4)') NMOMS

      DO I=0, NMOMS
        WRITE(100, '(I4, 6E12.4)') I, EV_T(I, :)
      END DO
      
      CLOSE(100)
      
      END
      

      SUBROUTINE MAKESCATFILE1(FNAME, MIDX, R0, R1, GAMMA, NPTS, WL,
     .                        TAUA, TAUM, SSA_A, NMOMS, TAUE_T, TAUS_T, 
     .                        SSA_T, EV_T)
CF2PY INTENT(IN) FNAME, MIDX, R0, R1, GAMMA, TAUA, TAUM, NMOMS, SSA_A
CF2PY INTENT(IN) NPTS
CF2PY INTENT(OUT) TAUE_T, TAUS_T, SSA_T, EV_T
      CHARACTER*64  FNAME
      COMPLEX*16    MIDX
      REAL*8        R0, R1, GAMMA, WL, TAUM, TAUA, DR, SSA_A
      REAL*8        TAUE_T, TAUS_T, SSA_T, EV_T(0:NMOMS, 6), VOL
      REAL*8        EV_M(0:NMOMS, 6), WIS_A(0:NMOMS, 4), EXT, SCA, ASY
      REAL*8        X(NPTS), Y(NPTS), EV_A(0:NMOMS, 6)
      INTEGER       NMOMS, NPTS, IERR, I
      REAL*8        OMEGA_A, TAUS_A

      EXTERNAL RAYEV, DISTR, WISCOMBE2EVANS

      
      TAUE_T = TAUA+TAUM

      DR = (R1-R0)/DBLE(NPTS-1)
      
      DO I=1, NPTS
        X(I) = DBLE(I-1)*DR+R0
        Y(I) = (X(I)/R0)**GAMMA
      END DO
      
      
      CALL DISTR(NPTS, X, Y, MIDX, WL, NMOMS, WIS_A, EXT, SCA, ASY, 
     .                VOL, IERR)

      
      IF (IERR /= 0) THEN
        STOP
      END IF

      CALL WISCOMBE2EVANS(NMOMS, WIS_A, EV_A)
      CALL RAYEV(EV_M, NMOMS)

      OMEGA_A = SCA/EXT

C     Если SSA_A на отрезке [0..1], то примем это значение за альбедо 
C     однократного рассеяния
      IF ((SSA_A .LE. 1.0D0).AND.(SSA_A .GE. 0.0D0)) OMEGA_A = SSA_A 

      TAUS_A = TAUA * OMEGA_A

      TAUS_T = TAUS_A + TAUM
      SSA_T = TAUS_T/TAUE_T

      
      EV_T = (EV_A*TAUS_A+EV_M*TAUM)/TAUS_T
      
      OPEN(UNIT=100, FILE=FNAME, STATUS='UNKNOWN')

      WRITE(100, '(E12.4)') TAUE_T
      WRITE(100, '(E12.4)') TAUS_T
      WRITE(100, '(E12.4)') SSA_T
      WRITE(100, '(I4)') NMOMS

      DO I=0, NMOMS
        WRITE(100, '(I4, 6E12.4)') I, EV_T(I, :)
      END DO
      
      CLOSE(100)
      
      END
      