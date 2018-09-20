      SUBROUTINE RUN(LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .               TAUM, NMOMS, DIRECT_FLUX, DIRECT_MU, 
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU, 
     .               MU, IV, QV)
      

CF2PY INTENT(IN)  LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA
CF2PY INTENT(IN)  DIRECT_FLUX, DIRECT_MU, NUMMU
CF2PY INTENT(IN)  TAUM, NMOMS, QUAD_TYPE, DELTAM, GROUND_ALBEDO
CF2PY INTENT(OUT) MU, IV, QV

      CHARACTER*64  LAYF, OUTF
      COMPLEX*16    MIDX, GROUND_INDEX
      REAL*8        R0, R1, GAMMA, WL, TAUA, TAUM
      INTEGER       NPTS, NMOMS

      EXTERNAL RT3, MAKELAYFL, POSTPROCESS

      INTEGER NSTOKES, NUMMU, AZIORDER, SRC_CODE, NUMAZIMUTHS
C      PARAMETER(NUMMU=32)   
      REAL*8  MU_VALUES(2*NUMMU), DIRECT_FLUX, DIRECT_MU, GROUND_TEMP,
     .        GROUND_ALBEDO, SKY_TEMP, MU(2*NUMMU), IV(2*NUMMU), 
     .        QV(2*NUMMU)
      
      CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1

      


      CALL MAKELAYFL(LAYF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, TAUM, 
     .            NMOMS)

      
      NSTOKES = 2
      AZIORDER = 6
      SRC_CODE = 1
      GROUND_TYPE = 'L'
      GROUND_INDEX = (1.0, 0.0)
      UNITS = 'W'
      OUTPOL = 'IQ'
      NUMAZIMUTHS = 2

      CALL RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES,
     .                    SRC_CODE, LAYF, OUTF,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WL, UNITS, OUTPOL,
     .                    NUMAZIMUTHS)


      CALL POSTPROCESS(OUTF, NUMMU, MU, IV, QV)

      END


      SUBROUTINE RUN0(LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .               TAUM, SSA_A, NMOMS, DIRECT_FLUX, DIRECT_MU, 
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU, 
     .               MU, IV, QV)
      

CF2PY INTENT(IN)  LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA
CF2PY INTENT(IN)  DIRECT_FLUX, DIRECT_MU, NUMMU, SSA_A
CF2PY INTENT(IN)  TAUM, NMOMS, QUAD_TYPE, DELTAM, GROUND_ALBEDO
CF2PY INTENT(OUT) MU, IV, QV

      CHARACTER*64  LAYF, OUTF
      COMPLEX*16    MIDX, GROUND_INDEX
      REAL*8        R0, R1, GAMMA, WL, TAUA, TAUM, SSA_A
      INTEGER       NPTS, NMOMS

      EXTERNAL RT3, MAKELAYFL1, POSTPROCESS

      INTEGER NSTOKES, NUMMU, AZIORDER, SRC_CODE, NUMAZIMUTHS
C      PARAMETER(NUMMU=32)   
      REAL*8  MU_VALUES(2*NUMMU), DIRECT_FLUX, DIRECT_MU, GROUND_TEMP,
     .        GROUND_ALBEDO, SKY_TEMP, MU(2*NUMMU), IV(2*NUMMU), 
     .        QV(2*NUMMU)
      
      CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1

      


      CALL MAKELAYFL1(LAYF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, TAUM, 
     .            SSA_A, NMOMS)

      
      NSTOKES = 2
      AZIORDER = 6
      SRC_CODE = 1
      GROUND_TYPE = 'L'
      GROUND_INDEX = (1.0, 0.0)
      UNITS = 'W'
      OUTPOL = 'IQ'
      NUMAZIMUTHS = 2

      CALL RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES,
     .                    SRC_CODE, LAYF, OUTF,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WL, UNITS, OUTPOL,
     .                    NUMAZIMUTHS)


      CALL POSTPROCESS(OUTF, NUMMU, MU, IV, QV)

      END


      SUBROUTINE RUN1(LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .               NMOMS, DIRECT_MU,
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, 
     .               NUMMU, MU, IV, QV)
      

CF2PY INTENT(IN)  LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA
CF2PY INTENT(IN)  DIRECT_MU, NUMMU
CF2PY INTENT(IN)  NMOMS, QUAD_TYPE, DELTAM, GROUND_ALBEDO
CF2PY INTENT(OUT) MU, IV, QV

      CHARACTER*64  LAYF, OUTF
      COMPLEX*16    MIDX, GROUND_INDEX
      REAL*8        R0, R1, GAMMA, WL, TAUA, TAUM, GROUND_ALBEDO
      INTEGER       NPTS, NMOMS, IERR, NUMMU
      CHARACTER     QUAD_TYPE*1, DELTAM*1
      REAL*8        DIRECT_FLUX, DIRECT_MU, MU(2*NUMMU), IV(2*NUMMU), 
     .              QV(2*NUMMU)

      REAL*8   GETTAUM
      EXTERNAL RUN, GETTAUM


C       PRINT '(F12.4)', WL
      TAUM = GETTAUM(WL)
C       PRINT *, TAUM

      DIRECT_FLUX = 3.141592653589793D0*DIRECT_MU
      CALL RUN(LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .               TAUM, NMOMS, DIRECT_FLUX, DIRECT_MU,
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU,
     .               MU, IV, QV)


      END


      SUBROUTINE RUN2(LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .               SSA_A, NMOMS, DIRECT_MU,
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, 
     .               NUMMU, MU, IV, QV)
      

CF2PY INTENT(IN)  LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA
CF2PY INTENT(IN)  DIRECT_MU, NUMMU, SSA_A
CF2PY INTENT(IN)  NMOMS, QUAD_TYPE, DELTAM, GROUND_ALBEDO
CF2PY INTENT(OUT) MU, IV, QV

      CHARACTER*64  LAYF, OUTF
      COMPLEX*16    MIDX, GROUND_INDEX
      REAL*8        R0, R1, GAMMA, WL, TAUA, TAUM, GROUND_ALBEDO, SSA_A
      INTEGER       NPTS, NMOMS, IERR, NUMMU
      CHARACTER     QUAD_TYPE*1, DELTAM*1
      REAL*8        DIRECT_FLUX, DIRECT_MU, MU(2*NUMMU), IV(2*NUMMU), 
     .              QV(2*NUMMU)

      REAL*8   GETTAUM
      EXTERNAL RUN, GETTAUM


C       PRINT '(F12.4)', WL
      TAUM = GETTAUM(WL)
C       PRINT *, TAUM

      DIRECT_FLUX = 3.141592653589793D0*DIRECT_MU
      CALL RUN0(LAYF, OUTF, MIDX, R0, R1, GAMMA, NPTS, WL, TAUA, 
     .               TAUM, SSA_A, NMOMS, DIRECT_FLUX, DIRECT_MU,
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU,
     .               MU, IV, QV)


      END


      SUBROUTINE POSTPROCESS(FNAME, NUMMU, MU, Iv, Qv)
C     Загружаем результаты расчетов из файла и передаем их в вызывающую
C     подпрограмму
C     
CF2PY INTENT(IN)  FNAME, NUMMU
CF2PY INTENT(OUT) MU, Iv, Qv
      CHARACTER*64  FNAME
      INTEGER NUMMU
      REAL*8  MU(2*NUMMU), Iv(2*NUMMU), Qv(2*NUMMU), TMP(2:NUMMU)
      REAL*8  z_i, phi_i, mu_i, I_i, Q_i, Temp

      CHARACTER*72 BUFFER
      INTEGER I,J

      OPEN(UNIT=100, FILE=FNAME, STATUS='OLD')

      DO 5 I=1, 11
        READ(100, '(A)') BUFFER
 5    CONTINUE

      I=1
      DO 10
        READ(100, *, IOSTAT=IERR) z_i, phi_i, mu_i, I_i, Q_i

        IF(IERR/=0) EXIT

        if ((mu_i>0) .and. (mu_i /= 2.0)) then
          MU(I) = ACOS(mu_i)*(COS(phi_i*0.017453292519943))/
     .              0.017453292519943

C          MU(I) = mu_i*(COS(phi_i*0.017453292519943))
          Iv(I) = I_i
          Qv(I) = Q_i
          I=I+1
        end if
      
 10   CONTINUE

      MU(NUMMU+1:) = MU(2*NUMMU:NUMMU:-1)
      Iv(NUMMU+1:) = Iv(2*NUMMU:NUMMU:-1)
      Qv(NUMMU+1:) = Qv(2*NUMMU:NUMMU:-1)

C Simply reverse arrays
      MU = MU(2*NUMMU:1:-1)
      Iv = Iv(2*NUMMU:1:-1)
      Qv = Qv(2*NUMMU:1:-1)
      
C Reverse arrays 
C       DO i = 1, NUMMU       ! Head moves forward
C         j = 2*NUMMU -i + 1     ! J = Tail moves backward
C         Temp  = MU(i)     ! swap elements
C         MU(i)  = MU(j)
C         MU(j)  = Temp

C         Temp  = Iv(i)     ! swap elements
C         Iv(i)  = Iv(j)
C         Iv(j)  = Temp

C         Temp  = Qv(i)     ! swap elements
C         Qv(i)  = Qv(j)
C         Qv(j)  = Temp

C       END DO

C      MU = MU(:-1)
C      Iv = Iv(:-1)
C      Qv = Qv(:-1)
      
      CLOSE(100)

      END