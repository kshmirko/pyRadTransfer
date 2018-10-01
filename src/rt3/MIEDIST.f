C     Можно уменьшить размер требуемой памяти для программы
      SUBROUTINE DISTR(M, X, Y, MIDX, WL, MXMDM, PMOM, EXT, SCA, ASY, 
     .                VOL, IERR)
C     Расчет коэффициентов лежандра для фазовых функций некоторого
C     распределения частиц по размерам для заданной длины волны и 
C     показателя преломления материала чатиц
C     Помимо фазовой функции вычисляюются коэффициенты экстинкции, рассеяния
C     и ассиметрии, а также объем частиц заданного распределения в 
C     пересчете на единичную концентрацию (1 см3)
C     
C     Выходные параметры
C     ------------------
C     Ext - коэффициент экстинкции [m^-1]
C     Sca - коэффициент рассеяния [m^-1]
C     Asy - коэффициент ассиметрии [m^-1]
C     Vol - объем [m^-1]
C     pmom - фазовые функции в виде коэффициентов лежандра 
C     IERR = 0  NO ERRORS
C     IERR = -1 MAXIMUM DIMENSION LESS THAN REQUIRED AMOUNT OF LEGENDRE 
C     COEFFICIENTS, INCREASE MXMDM that 2*INT(MAX(XX))<=MXMDM
     
      IMPLICIT NONE
CF2PY INTENT(IN) X, Y, MIDX, WL, MXMDM, M
CF2PY INTENT(OUT) PMOM, EXT, SCA, ASY, VOL, IERR

      INTEGER M, MXMDM, IERR
      REAL*8  EXT, SCA, ASY, VOL, X(M), Y(M), XX, WL
      REAL*8  PMOM(0:MXMDM, 4)
      COMPLEX*16  MIDX
      
      REAL*8  PI, QEXT, QSCA, GQSC
      PARAMETER(PI=3.141592653589793D0)
      
      REAL*8  PMOM_TMP(M, 0:MXMDM, 4), PMOM_I(0:MXMDM, 4)
      REAL*8  EXTA(M), SCATA(M), ASYA(M), VOLA(M)
      REAL*8  K, XSQUARED, XCUBED
      
      INTEGER I, J, L
      REAL*8 TRAPZ, NORM, GFACTOR
      EXTERNAL MIEDRV, TRAPZ
      
      
      
      PMOM = 0.0D0
      EXT = 0.0D0
      SCA = 0.0D0
      ASY = 0.0D0
      VOL = 0.0D0
      
      K = 2.0D0 * PI / WL
      
      DO I=1, M
        XX = K * X(I)
        
C       Расчет параметров Ми для каждого конкретного размера частиц
        CALL MIEDRV(XX, MIDX, QEXT, QSCA, GQSC, PMOM_I, MXMDM, IERR)
        IF (IERR /= 0) THEN
          RETURN
        END IF
        
        XSQUARED = X(I) * X(I) * 1.0D-12
        XCUBED = XSQUARED * X(I) * 1.0D-6
        
        EXTA(I) = QEXT*XSQUARED*Y(I)
        SCATA(I) = QSCA*XSQUARED*Y(I)
        ASYA(I) = GQSC*SCATA(I)
        VOLA(I) = XCUBED*Y(I)
        
C       Формируем матрицу для дальнейшего усреднения
        DO J=1, 4
          
          DO L=0, MXMDM
            
            PMOM_TMP(I, L, J) = PMOM_I(L,J)*Y(I)
          
          END DO
        
        END DO
      
      END DO
      
      NORM = TRAPZ(X, Y, M)
      EXT = TRAPZ(X, EXTA, M) / NORM
      SCA = TRAPZ(X, SCATA, M) / NORM
      ASY = TRAPZ(X, ASYA, M) / NORM
      VOL = TRAPZ(X, VOLA, M) / NORM

      
      
C     Integrate over first dimension
      
      DO I=1, 4
        DO J=0, MXMDM
          
          PMOM(J,I) = TRAPZ(X, PMOM_TMP(:, J, I), M) / NORM
  
        END DO
      END DO

     
      IERR = 0
      RETURN
      END 
      
