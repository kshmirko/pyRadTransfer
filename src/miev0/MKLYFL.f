      SUBROUTINE MAKELAYFL(LAYF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .                      TAUM, NMOMS)

      COMPLEX*16  MIDX
      REAL*8      R0, R1, GAMMA, WL, TAUA, TAUM, TAUE_T, TAUS_T, SSA_T
      REAL*8      EV_T(0:NMOMS, 6)
      INTEGER     NPTS, NMOMS
      CHARACTER*64  LAYF, SCATF

      WRITE(SCATF,'(A)') 'scat_file'


      CALL MAKESCATFILE(SCATF, MIDX, R0, R1, GAMMA, NPTS, WL,
     .                        TAUA, TAUM, NMOMS, TAUE_T, TAUS_T, SSA_T, 
     .                        EV_T)


      OPEN(UNIT=200, FILE=LAYF, STATUS='UNKNOWN')
      WRITE(200, "(2F7.2,F7.3,A,A)") 1.0, 0.0, 0.0, '   ', 
     .                                "'"//TRIM(SCATF)//"'"
      WRITE(200, "(2F7.2,F7.3,A,A)") 0.0, 0.0, 0.0, '   ', "'         '"  
      CLOSE(200)
      END