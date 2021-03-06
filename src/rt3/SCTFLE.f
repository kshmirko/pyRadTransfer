      SUBROUTINE MAKE_SCAT_FILE(FNAME, NMOMS, PWIS, TAUAE, TAUME, SSA,
     .                          SSAT, EV_T)
C Создаем файл рассеяния на основе входных параметров
      IMPLICIT NONE
CF2PY INTENT(IN)  FNAME, NMOMS, PWIS, TAUAE, TAUME, SSA
CF2PY INTENT(OUT) SSAT, EV_T
      CHARACTER*64  FNAME
      INTEGER       NMOMS, I
      REAL*8        PWIS(0:NMOMS, 4), TAUAE, SSA, TAUME, SSAT
      EXTERNAL RAYEV, WISCOMBE2EVANS
      REAL*8        EV_M(0:NMOMS, 6), EV_A(0:NMOMS, 6), EV_T(0:NMOMS, 6)
      REAL*8        TAUAS, TAUTS, TAUTE


      CALL WISCOMBE2EVANS(NMOMS, PWIS, EV_A)
      CALL RAYEV(EV_M, NMOMS)

      TAUAS = TAUAE*SSA
      TAUTS = TAUME+TAUAS

      TAUTE = TAUAE+TAUME
      EV_T = (TAUAS*EV_A+TAUME*EV_M)/TAUTS

      SSAT = TAUTS / TAUTE

      OPEN(UNIT=100, FILE=FNAME, STATUS='UNKNOWN')

      WRITE(100, '(E12.4)') TAUTE
      WRITE(100, '(E12.4)') TAUTS
      WRITE(100, '(E12.4)') SSAT
      WRITE(100, '(I4)') NMOMS

      DO I=0, NMOMS
        WRITE(100, '(I4, 6E12.4)') I, EV_T(I, :)
      END DO
      
      CLOSE(100)

      END