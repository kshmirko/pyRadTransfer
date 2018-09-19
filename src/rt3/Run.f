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
      INTEGER IERR, IDSTR, I, NM
      COMPLEX*16  MIDX
      REAL*8  SSA_A, R0, R1, P(1), TAUE_A, TAUE_M, SSAT, 
     .        EV_T(0:NMOMS, 6), DIRECT_FLUX, DIRECT_MU,
     .        GROUND_ALBEDO, SZA
      CHARACTER*64  SCAT_FILE, LAYER_FILE, OUT_FILE, MEAS_FILE
      CHARACTER     QUAD_TYPE*1, DELTAM*1

      REAL*8  MU0(2*NUMMU), I0(2*NUMMU), Q0(2*NUMMU)
      REAL*8  MU1(2*NUMMU), I1(2*NUMMU), Q1(2*NUMMU)
      REAL*8  I2(2*NUMMU), Q2(2*NUMMU)
      REAL*8  AMEAS(MAXMEAS), IMEAS(MAXMEAS), QMEAS(MAXMEAS)

      EXTERNAL DISTR, MKDSTRB, GETTAUM, MAKE_SCAT_FILE, MAKE_LAYER_FILE,
     .         RUNRT3, PWL_VALUE_1D, READ_MEAS_FILE
      REAL*8  GETTAUM

      R0 = 0.1D0
      R1 = 1.0D0
      IDSTR = 1
      P(1)  = -3.5D0

C     Fill X and Y with distribution function
      CALL MKDSTRB(NPTS, R0, R1, IDSTR, 1, P, X, Y)

      MIDX = (1.4, 0.01)
      WL  = 0.750

C     Calculates scattering matrix and coefficients
      CALL DISTR(NPTS, X, Y, MIDX, WL, NMOMS, PMOM, EXT, SCA, ASY, 
     .                VOL, IERR)

      IF (IERR .NE. 0) THEN
        PRINT *, 'ERROR CALCULATING OPTICAL PROPERTIES'
        STOP
      ENDIF

C     SINGLE SCATTERING ALBEDO FOR AEROSOLS
      SSA_A = SCA / EXT
      SSA_A = 0.8

      TAUE_M = GETTAUM(WL)
      TAUE_A = 0.1
      SCAT_FILE = 'scat_file'
      LAYER_FILE = 'atmos.lay'

C     Save layer file
      CALL MAKE_LAYER_FILE(LAYER_FILE, SCAT_FILE)


C     ******************************************************************

C     Create scattering file

      CALL MAKE_SCAT_FILE(SCAT_FILE, NMOMS, PMOM, TAUE_A, TAUE_M, SSA_A,
     .                          SSAT, EV_T)

C     Подготовка к расчету
      SZA=74.0D0

      OUT_FILE = 'rt3.out'
      QUAD_TYPE = 'G'
      DELTAM = 'N'
      DIRECT_FLUX = 27.0D0
      DIRECT_MU =COS(SZA*DEG2RAD)
      GROUND_ALBEDO = 0.0

      CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU, 
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU, 
     .               MU0, I0, Q0)


      GROUND_ALBEDO = 0.1
      CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU, 
     .               QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU, 
     .               MU1, I1, Q1)

      MEAS_FILE = 'meas.dat'
C     Приведем измерения к единой сетке      
      CALL READ_MEAS_FILE(MEAS_FILE, MAXMEAS, AMEAS, IMEAS, QMEAS, NM)

      AMEAS(:NM) = SZA-AMEAS(:NM)
      print *, AMEAS(:NM)

      CALL PWL_VALUE_1D ( NM, AMEAS(NM:1:-1), IMEAS(NM:1:-1), 2*NUMMU, 
     .                    MU0, I2 )
      CALL PWL_VALUE_1D ( NM, AMEAS(NM:1:-1), QMEAS(NM:1:-1), 2*NUMMU, 
     .                    MU0, Q2 )

      DO I=1, NM
        PRINT '(2F15.3)', AMEAS(I), IMEAS(I)
      ENDDO


      DO I=1, NUMMU*2
        PRINT '(I5, 3F12.5)', I, MU0(I), I0(I), I2(I)
      ENDDO

      RETURN
      END


      SUBROUTINE CHANGE_DEG(N, ANG, SZA)
      IMPLICIT NONE
      INTEGER N, I
      REAL*8  ANG(N), SZA

      DO I=1, N
        ANG(I) = SZA-ANG(I)
      ENDDO

      END