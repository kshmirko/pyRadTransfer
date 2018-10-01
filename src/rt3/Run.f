      SUBROUTINE RUNDRV()
C NMOMS - количество коэффициентов разложения по полиномам Лежандра
C NPTS  - количество точек табулирования функции разложения
C NUMMU - количество углов для радиационного расчета [0-90]
C DEG2RAD, RAD2DEG  - коэффициенты перевода градусов в радианы и обратно
C X(NPTS), Y(NPTS)  - абсцисса и ордината функции распределения
C WL, EXT, SCA, ASY, VOL  - длина волны и параметры распределения в
C     приведении к единичной концентрации
C     PMOM(0:NMOMS,4)     - коэффициенты лежвндра матрицы рассеяния
C     MIDX                - показатель проеломления частицы
C
      IMPLICIT NONE
      INTEGER NMOMS, NPTS, NUMMU, MAXMEAS
      REAL*8 DEG2RAD, RAD2DEG
      PARAMETER(NMOMS=40, NPTS=101, NUMMU=32,
     .        DEG2RAD=0.017453292519943D0, RAD2DEG=1.0D0/DEG2RAD,
     .        MAXMEAS=20)

      REAL*8  X(NPTS), Y(NPTS), WL, PMOM(0:NMOMS, 4), EXT, SCA, ASY, VOL
      INTEGER IERR, IDSTR, I, NM, K
      COMPLEX*16  MIDX
      REAL*8  SSA_A, R0, R1, P(1), TAUE_A, TAUE_M, SSAT,
     .        EV_T(0:NMOMS, 6), DIRECT_FLUX, DIRECT_MU,
     .        GROUND_ALBEDO, SZA, EV_A(0:NMOMS,6), EV_R(0:NMOMS,6)
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
      REAL*8  WW(3)
      DATA WW /0.8, 1.0, 0.9/

      R0 = 0.1D0
      R1 = 1.0D0
      IDSTR = 1
      P(1)  = -3.5D0

C     Fill X and Y with distribution function
      CALL MKDSTRB(NPTS, R0, R1, IDSTR, 1, P, X, Y)

      MIDX = (1.6, 0.00)
      WL  = 0.750

C     Calculates scattering matrix and coefficients
      CALL DISTR(NPTS, X, Y, MIDX, WL, NMOMS, PMOM, EXT, SCA, ASY,
     .                VOL, IERR)

      IF (IERR .NE. 0) THEN
        PRINT *, 'ERROR CALCULATING OPTICAL PROPERTIES'
        STOP
      ENDIF


      CALL WISCOMBE2EVANS(NMOMS, PMOM, EV_A)
      CALL RAYEV(EV_R, NMOMS)

C     MAKE AEROSOL AND RAYLEIGHT PHASE FUNCTION
      CALL LEGVAL_A(NMOMS, EV_A(:,1), 2*NUMMU, XI, PA)
      CALL LEGVAL_A(NMOMS, EV_R(:,1), 2*NUMMU, XI, PM)

C     SINGLE SCATTERING ALBEDO FOR AEROSOLS
      SSA_A = SCA / EXT
C     SSA_A = 0.8

      TAUE_M = GETTAUM(WL)
      TAUE_A = 0.01
      SCAT_FILE = 'scat_file'
      LAYER_FILE = 'atmos.lay'

C     Save layer file
      CALL MAKE_LAYER_FILE(LAYER_FILE, SCAT_FILE)


C     ******************************************************************
      MEAS_FILE = 'meas.dat'
C     Приведем измерения к единой сетке
      CALL READ_MEAS_FILE(MEAS_FILE, MAXMEAS, AMEAS, IMEAS, QMEAS, NM)
      AMEAS = COS(AMEAS*0.017453292519943)

      CALL PWL_VALUE_1D ( NM, AMEAS(NM:1:-1), IMEAS(NM:1:-1), 2*NUMMU,
     .                    XI, I2, 0.0D0, 0.0D0 )

C     Подготовка к расчету
      SZA=74.0D0

      OUT_FILE = 'rt3.out'
      QUAD_TYPE = 'G'
      DELTAM = 'N'
      DIRECT_FLUX = 27.0D0
      DIRECT_MU =COS(SZA*DEG2RAD)

      DO K=1, 3
      SSA_A = WW(K)
C     Create scattering file
      CALL MAKE_SCAT_FILE(SCAT_FILE, NMOMS, PMOM, TAUE_A, TAUE_M, SSA_A,
     .                          SSAT, EV_T)


      GROUND_ALBEDO = 0.0
      CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU,
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU,
     .               MU0, I0, Q0)


      GROUND_ALBEDO = 0.3
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

      DO I=1, 2*NUMMU
        IF (I2(I) .NE. 0.0D0) THEN
          PA1(I) = ((I2(I) - DL(I))/L0(I))*PA(I)+((I2(I)-L1(I))/L0(I))*
     .              TAUE_M*PM(I)/(SSA_A*TAUE_A)
          PRINT *, L0(I)
        ENDIF

      ENDDO

      PRINT "('W(',I2,')',F10.6)", K, SUM(PA1*WEIGHT)

      ENDDO
      RETURN
      END
      
