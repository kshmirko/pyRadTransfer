      SUBROUTINE RUNINV(NMOMS, NPTS, MAXMEAS, WL, MIDX, R0, R1, NP, P,
     .                  IDSTR, GROUND_ALBEDO, NWW, WW, TAUE_A, SZA, 
     .                  DIRECT_FLUX, MEAS_FILE, XO, PA1, PA, I2, L0, L1)
CF2PY INTENT(IN)  NMOMS, NPTS, MAXMEAS, WL, MIDX, R0, R1, NP, P, IDSTR
CF2PY INTENT(IN)  GROUND_ALBEDO, NWW, WW, TAUE_A, SZA, DIRECT_FLUX
CF2PY INTENT(IN)  MEAS_FILE
CF2PY INTENT(OUT) XO, PA1, PA, I2, L0, L1
      IMPLICIT NONE
      INTEGER NMOMS, NPTS, NUMMU, MAXMEAS
      REAL*8 DEG2RAD, RAD2DEG, DELTA
      PARAMETER(NUMMU=32, DEG2RAD=0.017453292519943D0, 
     .          RAD2DEG=1.0D0/DEG2RAD, DELTA=0.1D0)

      REAL*8  X(NPTS), Y(NPTS), WL, PMOM(0:NMOMS, 4), EXT, SCA, ASY, VOL
      INTEGER IERR, IDSTR, I, NM, K, NP, NWW
      COMPLEX*16  MIDX
      REAL*8  SSA_A, R0, R1, P(NP), TAUE_A, TAUE_M, SSAT, 
     .        EV_T(0:NMOMS, 6), DIRECT_FLUX, DIRECT_MU,
     .        GROUND_ALBEDO, SZA, EV_A(0:NMOMS,6), EV_R(0:NMOMS,6),
     .        XO(2*NUMMU)
      CHARACTER*64  SCAT_FILE, LAYER_FILE, OUT_FILE, MEAS_FILE
      CHARACTER     QUAD_TYPE*1, DELTAM*1

      REAL*8  MU0(2*NUMMU), I0(2*NUMMU), Q0(2*NUMMU)
      REAL*8  MU1(2*NUMMU), I1(2*NUMMU), Q1(2*NUMMU)
      REAL*8  I2(2*NUMMU), Q2(2*NUMMU), DL(2*NUMMU)
      REAL*8  AMEAS(MAXMEAS), IMEAS(MAXMEAS), QMEAS(MAXMEAS)
      REAL*8  L0(2*NUMMU), L1(2*NUMMU), L(2*NUMMU), PA(2*NUMMU), 
     .        PM(2*NUMMU), PA1(2*NUMMU)
      INCLUDE 'GAUSCOEF.f'
      EXTERNAL DISTR, MKDSTRB, GETTAUM, MAKE_SCAT_FILE, MAKE_LAYER_FILE,
     .         RUNRT3, PWL_VALUE_1D, READ_MEAS_FILE, LEGVAL_A
      EXTERNAL RAYEV, WISCOMBE2EVANS

      REAL*8  GETTAUM
      REAL*8  WW(NWW), SSA_RET


C     Fill X and Y with distribution function
      CALL MKDSTRB(NPTS, R0, R1, IDSTR, NP, P, X, Y)

C     Calculates scattering matrix and coefficients
      CALL DISTR(NPTS, X, Y, MIDX, WL, NMOMS, PMOM, EXT, SCA, ASY, 
     .                VOL, IERR)

C     SINGLE SCATTERING ALBEDO FOR AEROSOLS
C      SSA_A = SCA / EXT
      IF (IERR .NE. 0) THEN
        PRINT *, 'ERROR CALCULATING OPTICAL PROPERTIES'
        STOP
      ENDIF


      CALL WISCOMBE2EVANS(NMOMS, PMOM, EV_A)
      CALL RAYEV(EV_R, NMOMS)

C     MAKE AEROSOL AND RAYLEIGHT PHASE FUNCTION
      CALL LEGVAL_A(NMOMS, EV_A(:,1), 2*NUMMU, XI, PA)
      CALL LEGVAL_A(NMOMS, EV_R(:,1), 2*NUMMU, XI, PM)


C     SSA_A = 0.8

      TAUE_M = GETTAUM(WL)
C      TAUE_A = 0.01
      SCAT_FILE = 'scat_file'
      LAYER_FILE = 'atmos.lay'

C     Save layer file
      CALL MAKE_LAYER_FILE(LAYER_FILE, SCAT_FILE)


C     ******************************************************************
C      MEAS_FILE = 'meas.dat'
C     Приведем измерения к единой сетке      
      CALL READ_MEAS_FILE(MEAS_FILE, MAXMEAS, AMEAS, IMEAS, QMEAS, NM)
      
      AMEAS = COS(AMEAS*0.017453292519943)
 
      CALL PWL_VALUE_1D ( NM, AMEAS(NM:1:-1), IMEAS(NM:1:-1), 2*NUMMU, 
     .                    XI, I2, 0.0D0, 0.0D0 )


C     Подготовка к расчету
C      SZA=74.0D0

      OUT_FILE = 'rt3.out'
      QUAD_TYPE = 'G'
      DELTAM = 'N'
C      DIRECT_FLUX = 27.0D0
      DIRECT_MU =COS(SZA*DEG2RAD)

      DO 100 K=1, NWW
      
      SSA_A = WW(K)
C     Create scattering file
      CALL MAKE_SCAT_FILE(SCAT_FILE, NMOMS, PMOM, TAUE_A, TAUE_M, SSA_A,
     .                          SSAT, EV_T)


      CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU, 
     .               QUAD_TYPE, DELTAM, 0.0D0, NUMMU, 
     .               MU0, I0, Q0)


      CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU, 
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU, 
     .               MU1, I1, Q1)


      MU0 = SZA - MU0

      DO I=1, 2*NUMMU
        IF ((MU0(I) .LT. 0.0D0) .OR. (MU0(I).GT.180.0D0)) THEN
          I0(I) = 0.0D0
          Q0(I) = 0.0D0
          I1(I) = 0.0D0
          Q1(I) = 0.0D0
          IF (MU0(I) .LT. 0.0D0) THEN
            MU0(I) = 0.0D0
          ELSE IF (MU0(I) .GT. 180.0D0) THEN
            MU0(I) = 180.0D0
          ENDIF
        ENDIF
      ENDDO

C     Convert degrees to cosines
      MU0 = COS(MU0*0.017453292519943)

C     Интерполируем на узлы гауссовой сетки наши измерения и результаты расчета


      CALL PWL_VALUE_1D ( 2*NUMMU, MU0, I0, 2*NUMMU, 
     .                    XI, L0, 0.0D0, 0.0D0 )

      CALL PWL_VALUE_1D ( 2*NUMMU, MU0, I1, 2*NUMMU, 
     .                    XI, L1, 0.0D0, 0.0D0 )


C     DL - вклад подстилающей поверхности
      DL = L1 - L0
      PA1(:) = PA(:)

C       PRINT '(10F9.4)', L(5:14), L0(5:14), L1(5:14)
      DO I=1, 2*NUMMU
        IF (I2(I) .NE. 0.0D0) THEN
          PA1(I) = ((I2(I) - DL(I))/L0(I))*PA(I)+((I2(I)-L1(I))/L0(I))*
     .              TAUE_M*PM(I)/(SSA_A*TAUE_A)
          
        ENDIF
        
      ENDDO
C      print '(10F9.4)', PA1
C      PRINT "('W0(',I2,')',F10.6)", K, SUM(PA*WEIGHT)
C     Интегрирование фазовой функци
      SSA_RET = SUM(PA1*WEIGHT) 
      PRINT "('W1(',I2,')',F10.6)", K, SSA_RET

      IF (ABS(SSA_RET - 2.0).LT.DELTA) THEN
        DO I=1, 2*NUMMU
          PRINT '(5F10.3)', XI(I), PA1(I), PA(I), L1(I), I2(I)
        ENDDO
        EXIT
      ENDIF
 100  CONTINUE
      XO = XI
      
      RETURN
      END
