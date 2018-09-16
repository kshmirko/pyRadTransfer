      SUBROUTINE MIEDRV(XX, MIDX, QEXT, QSCA, GQSC, PMOM, MXMDM, IERR)
C     Подпрограмма-драйвер для вычисления оптических свойств сферических
C     частиц
C     использует  
C     IERR = 0  NO ERRORS
C     IERR = -1 MAXIMUM DIMENSION LESS THAN REQUIRED AMOUNT OF LEGENDRE C               COEFFICIENTS
      IMPLICIT NONE  
      ! DUMMY PARAMETERS
CF2PY INTENT(IN)  XX, MIDX, MXMDM
CF2PY INTENT(OUT) QEXT, QSCA, GQSC, PMOM, IERR
      REAL*8                        XX
      COMPLEX*16                    MIDX
      INTEGER                       MXMDM
      REAL*8                        QEXT, QSCA, GQSC
      REAL*8                        PMOM(0:MXMDM,4)
      REAL                          PMOM_TMP(0:MXMDM, 4)
  
      ! actually used parametrs
      LOGICAL  ANYANG, PERFCT, PRNT(2)
      INTEGER  IPOLZN, NUMANG, NMOM
      REAL     GQSC_S, MIMCUT, QEXT_S, QSCA_S, SPIKE, XMU
      COMPLEX  CREFIN, SFORW, SBACK, S1(1), S2(1), TFORW(2), TBACK(2)
  
      EXTERNAL MIEV0
      INTEGER I, J, IERR
      
      ANYANG = .TRUE.
      PERFCT = .FALSE.
      PRNT   = .FALSE.
      NUMANG = 1
      NMOM   = 2*INT(XX)
      
      IF (NMOM.GT.MXMDM) THEN
        IERR = -1
        RETURN
      END IF
      
      IF (NMOM .LT.3) NMOM=3
      
      CREFIN = CMPLX(midx)
      MIMCUT = 1.0E-8
      IPOLZN = +1234
  
      XMU = 0.0
      
      PMOM_TMP = 0.0D0
  
      call MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     .            NMOM, IPOLZN, MXMDM, PRNT, QEXT_S, QSCA_S, GQSC_S,
     .            PMOM_TMP, SFORW, SBACK, S1, S2, TFORW, TBACK,
     .            SPIKE )
  
      ! COPY VALUES
      QEXT = DBLE(QEXT_S)
      QSCA = DBLE(QSCA_S)
      GQSC = DBLE(GQSC_S)
      PMOM = 0.0D0
      
      ! COPY MOMENTS
      PMOM = 0.0D0
      
      DO I=1, 4
        DO J=0, NMOM
          PMOM(J,I) = DBLE(PMOM_TMP(J,I))
        END DO
      END DO
      
      IERR = 0
      RETURN
  
      END
      
      
