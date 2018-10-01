      SUBROUTINE RUNINV(NMOMS, NPTS, MAXMEAS, WL, MIDX, R0, R1, NP, P,
     $     IDSTR, GROUND_ALBEDO, NWW, WW, TAUE_A, SZA,
     $     DIRECT_FLUX, MEAS_FILE, XO, PA1, PA, QA1, QA, I2, L0, L1, Q2,
     $     QL0, QL1)
C     ------------------------------------------------------------------
C     подпрограмма-драйвер для запуска процедуры восстановления яркости
CF2PY INTENT(IN)  NMOMS, NPTS, MAXMEAS, WL, MIDX, R0, R1, NP, P, IDSTR
CF2PY INTENT(IN)  GROUND_ALBEDO, NWW, WW, TAUE_A, SZA, DIRECT_FLUX
CF2PY INTENT(IN)  MEAS_FILE
CF2PY INTENT(OUT) XO, PA1, PA, QA1, QA, I2, L0, L1, Q2, QL0, QL1
      IMPLICIT NONE
      INTEGER NMOMS, NPTS, NUMMU, MAXMEAS
      REAL*8 DEG2RAD, RAD2DEG, DELTA
      PARAMETER(NUMMU=32, DEG2RAD=0.017453292519943D0,
     $     RAD2DEG=1.0D0/DEG2RAD, DELTA=0.1D0)

      REAL*8  X(NPTS), Y(NPTS), WL, PMOM(0:NMOMS, 4), EXT, SCA, ASY, VOL
      INTEGER IERR, IDSTR, I, NM, K, NP, NWW
      COMPLEX*16  MIDX
      REAL*8  SSA_A, R0, R1, P(NP), TAUE_A, TAUE_M, SSAT,
     $     EV_T(0:NMOMS, 6), DIRECT_FLUX, DIRECT_MU,
     $     GROUND_ALBEDO, SZA, EV_A(0:NMOMS,6), EV_R(0:NMOMS,6),
     $     XO(2*NUMMU)
      CHARACTER*64  SCAT_FILE, LAYER_FILE, OUT_FILE, MEAS_FILE
      CHARACTER     QUAD_TYPE*1, DELTAM*1

      REAL*8  MU0(2*NUMMU), I0(2*NUMMU), Q0(2*NUMMU)
      REAL*8  MU1(2*NUMMU), I1(2*NUMMU), Q1(2*NUMMU)
      REAL*8  I2(2*NUMMU), Q2(2*NUMMU), DL(2*NUMMU)
      REAL*8  AMEAS(MAXMEAS), IMEAS(MAXMEAS), QMEAS(MAXMEAS)
      REAL*8  L0(2*NUMMU), L1(2*NUMMU), L(2*NUMMU), PA(2*NUMMU),
     $     PM(2*NUMMU), PA1(2*NUMMU), QA1(2*NUMMU), QA(2*NUMMU),
     $     QM(2*NUMMU), QL0(2*NUMMU), QL1(2*NUMMU)
      INCLUDE 'GAUSCOEF.f'
      EXTERNAL DISTR, MKDSTRB, GETTAUM, MAKE_SCAT_FILE, MAKE_LAYER_FILE,
     $     RUNRT3, PWL_VALUE_1D, READ_MEAS_FILE, LEGVAL_A
      EXTERNAL RAYEV, WISCOMBE2EVANS

      REAL*8  GETTAUM

C     WW(NWW)  - передаваемый массив реперных значений альбедо
C     однократного рассеяния  для рассчета коррекции фазовой функции
C
C     PINT(NWW)- значение интеграла фазовой функции
      REAL*8  WW(NWW), PINT(NWW)
      REAL*8  SSA_RET


C     Запрлняем переменные X и Y степенным распределением с заданными
C     параметрами
      CALL MKDSTRB(NPTS, R0, R1, IDSTR, NP, P, X, Y)

C     Вычисляем по заданному распределению матрицу рассеяния и
C     оптические коэффициенты
      CALL DISTR(NPTS, X, Y, MIDX, WL, NMOMS, PMOM, EXT, SCA, ASY,
     $     VOL, IERR)

C     проверка на корректность завершения предыдущей команды
      IF (IERR .NE. 0) THEN
         PRINT *, 'INCREASE <NMOMS> VALUE'
         STOP
      ENDIF

C     Преобразуем матрицу рассеяния в формате Эванса в формат Вискомба.
C     Последняя матрица нормирована и первые два столбца соответствую
C     фазовым функциям I и Q
      CALL WISCOMBE2EVANS(NMOMS, PMOM, EV_A)

C     Заполняем матрицу эванса для молекулярного (Рэлеевского)
C     рассеяния. По размеру она идентична матрице аэрозольного рассеяния
C     EV_A, что дает возможность выполнить смешивание рассеивающих
C     свойств разных частиц
      CALL RAYEV(EV_R, NMOMS)

C     Вычисляем фазовые функции интенсивности для аэрозольного и
C     молекулярного рассеяния. Вычисляем в точках интегрирования по
C     методу Гаусса.
      CALL LEGVAL_A(NMOMS, EV_A(:,1), 2*NUMMU, XI, PA)
      CALL LEGVAL_A(NMOMS, EV_A(:,2), 2*NUMMU, XI, QA)
      CALL LEGVAL_A(NMOMS, EV_R(:,1), 2*NUMMU, XI, PM)
      CALL LEGVAL_A(NMOMS, EV_R(:,2), 2*NUMMU, XI, QM)


C     Вычисляем Рэлеевскую оптическую толщу
      TAUE_M = GETTAUM(WL)

      SCAT_FILE = 'scat_file'
      LAYER_FILE = 'atmos.lay'

C     Save layer file
      CALL MAKE_LAYER_FILE(LAYER_FILE, SCAT_FILE)

C     ------------------------------------------------------------------
C     Приведем измерения к единой сетке
      CALL READ_MEAS_FILE(MEAS_FILE, MAXMEAS, AMEAS, IMEAS, QMEAS, NM)

      AMEAS = COS(AMEAS*0.017453292519943)

      CALL PWL_VALUE_1D ( NM, AMEAS(NM:1:-1), IMEAS(NM:1:-1), 2*NUMMU,
     $     XI, I2, 0.0D0, 0.0D0 )

      CALL PWL_VALUE_1D ( NM, AMEAS(NM:1:-1), -QMEAS(NM:1:-1), 2*NUMMU,
     $     XI, Q2, 0.0D0, 0.0D0 )

      OUT_FILE = 'rt3.out'
      QUAD_TYPE = 'G'
      DELTAM = 'N'
      DIRECT_MU =COS(SZA*DEG2RAD)

      DO 100 K=1, NWW

         SSA_A = WW(K)

C     Создаем файл с коэффициентами рассеяния
         CALL MAKE_SCAT_FILE(SCAT_FILE, NMOMS, PMOM, TAUE_A, TAUE_M,
     $        SSA_A,SSAT, EV_T)

C     Запускаем расчет в отсутствие отражения от подстилающей
C     поверхности
         CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU,
     .        QUAD_TYPE, DELTAM, 0.0D0, NUMMU,
     .        MU0, I0, Q0)

C     Запускаем расчет, полагая альбедо подстилающей поверхности равным
C     GROUND_ALBEDO
         CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU,
     .        QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU,
     .        MU1, I1, Q1)


C     Теперь MU0 - угол рассеяния
         MU0 = SZA - MU0

C     зануляем элементы фазовой функции рассеяния для тех элементов,
C     которые лежат вне диапазона [0..180] градусов
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

C     переводим градусы в радианы
         MU0 = COS(MU0*0.017453292519943)

C     Интерполируем на узлы гауссовой сетки наши измерения и результаты
C     расчета
         CALL PWL_VALUE_1D ( 2*NUMMU, MU0, I0, 2*NUMMU,
     $        XI, L0, 0.0D0, 0.0D0 )

         CALL PWL_VALUE_1D ( 2*NUMMU, MU0, I1, 2*NUMMU,
     $        XI, L1, 0.0D0, 0.0D0 )


C     DL - вклад подстилающей поверхности
         DL = L1 - L0
         PA1(:) = PA(:)

C     Вычисляем коррекцию фазовой функции
         DO I=1, 2*NUMMU

            IF (I2(I) .NE. 0.0D0) THEN
               PA1(I) = ((I2(I) - DL(I))/L0(I))*PA(I)+((I2(I)-L1(I))
     $              /L0(I))*TAUE_M*PM(I)/(SSA_A*TAUE_A)

            ENDIF

         ENDDO

C     Интегрирование фазовой функци
         SSA_RET = SUM(PA1*WEIGHT)
C     Сохраняем значение интеграла
         PINT(K) = SSA_RET

         PRINT 1000, K, SSA_RET



 100  CONTINUE

      XO = XI
C     Проверим, есть ли пересечение графика WW-PINT с осью PINT = 2
      IF ((PINT(1) .GT. 2.0D0).AND.(PINT(NWW).LT.2.0D0)) THEN
C     Найдем альбедо, где PINT = 2.0D0
         SSA_RET = 2.0D0
         CALL PWL_VALUE_1D ( NWW, PINT(NWW:1:-1), WW(NWW:1:-1), 1,
     $        SSA_RET, SSA_A, 0.0D0, 0.0D0 )
         PRINT *, SSA_A
C     Нашли альбедо, повторим вычисления
C     Создаем файл с коэффициентами рассеяния
         CALL MAKE_SCAT_FILE(SCAT_FILE, NMOMS, PMOM, TAUE_A, TAUE_M,
     $        SSA_A,SSAT, EV_T)

C     Запускаем расчет в отсутствие отражения от подстилающей
C     поверхности
         CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU,
     .        QUAD_TYPE, DELTAM, 0.0D0, NUMMU,
     .        MU0, I0, Q0)

C     Запускаем расчет, полагая альбедо подстилающей поверхности равным
C     GROUND_ALBEDO
         CALL RUNRT3(LAYER_FILE, OUT_FILE, WL, DIRECT_FLUX, DIRECT_MU,
     .        QUAD_TYPE, DELTAM, GROUND_ALBEDO, NUMMU,
     .        MU1, I1, Q1)


C     Теперь MU0 - угол рассеяния
         MU0 = SZA - MU0

C     зануляем элементы фазовой функции рассеяния для тех элементов,
C     которые лежат вне диапазона [0..180] градусов
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

C     переводим градусы в радианы
         MU0 = COS(MU0*DEG2RAD)

C     Интерполируем на узлы гауссовой сетки наши измерения и результаты
C     расчета
         CALL PWL_VALUE_1D ( 2*NUMMU, MU0, I0, 2*NUMMU,
     $        XI, L0, 0.0D0, 0.0D0 )

         CALL PWL_VALUE_1D ( 2*NUMMU, MU0, I1, 2*NUMMU,
     $        XI, L1, 0.0D0, 0.0D0 )

         CALL PWL_VALUE_1D ( 2*NUMMU, MU0, Q0, 2*NUMMU,
     $        XI, QL0, 0.0D0, 0.0D0 )

         CALL PWL_VALUE_1D ( 2*NUMMU, MU0, Q1, 2*NUMMU,
     $        XI, QL1, 0.0D0, 0.0D0 )


C     DL - вклад подстилающей поверхности
         DL = L1 - L0
         PA1(:) = PA(:)
         QA1(:) = QA(:)
C     Вычисляем коррекцию фазовой функции
         DO I=1, 2*NUMMU

            IF (I2(I) .NE. 0.0D0) THEN
               PA1(I) = ((I2(I) - DL(I))/L0(I))*PA(I)+((I2(I)-L1(I))
     $              /L0(I))*TAUE_M*PM(I)/(SSA_A*TAUE_A)
            ENDIF
            IF (Q2(I) .NE. 0.0D0) THEN
               QA1(I) = Q2(I)/QL0(I)*QA(I)+((Q2(I)-QL0(I))/QL0(I))
     $              *TAUE_M*QM(I)/(SSA_A*TAUE_A)

            ENDIF

         ENDDO

C     Интегрирование фазовой функци
         SSA_RET = SUM(PA1*WEIGHT)


         PRINT 1000, 0, SSA_RET

C     Проверяем близость нормы полученной фазовой функции к 2.0 и если
C     это так, печатаем её, исходное приближение, яркость с учетом
C     подстилающей поверхности, измеренную яркость
C

         PRINT 1002, 'MU', 'PA1', 'PA', 'QA1', 'QA', 'L1', 'I2',
     $           'QL0', 'QL1', 'Q2'
         DO I=1, 2*NUMMU
            
            IF (I2(I) .NE. 0.0D0) THEN
               PRINT 1001, XI(I), PA1(I), PA(I), QA1(I), QA(I), L1(I),
     $              I2(I),  QL0(I), QL1(I), Q2(I)
            ENDIF

         ENDDO

       ENDIF


 1000 FORMAT(I3, 2X, 'NORM OF CORRECTED PHASE = ', F10.4)

 1001 FORMAT(10F10.3)
 1002 FORMAT(10A10)
      RETURN
      END
