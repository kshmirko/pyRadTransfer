      SUBROUTINE PWL_BASIS_1D ( ND, XD, NI, XI, BK )

C*********************************************************************72
C
CC PWL_BASIS_1D EVALUATES A 1D PIECEWISE LINEAR BASIS FUNCTION.
C
C  LICENSING:
C
C    THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.
C
C  MODIFIED:
C
C    04 OCTOBER 2017
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, INTEGER ND, THE NUMBER OF DATA POINTS.
C
C    INPUT, DOUBLE PRECISION XD(ND), THE DATA POINTS.
C
C    INPUT, INTEGER NI, THE NUMBER OF INTERPOLATION POINTS.
C
C    INPUT, DOUBLE PRECISION XI(NI), THE INTERPOLATION POINTS.
C
C    OUTPUT, DOUBLE PRECISION BK(NI), THE BASIS FUNCTION AT THE 
C    INTERPOLATION POINTS.
C
      IMPLICIT NONE

      INTEGER ND
      INTEGER NI

      DOUBLE PRECISION BK(NI,ND)
      INTEGER I
      INTEGER J
      DOUBLE PRECISION T
      DOUBLE PRECISION XD(ND)
      DOUBLE PRECISION XI(NI)

      DO J = 1, ND
        DO I = 1, NI
          BK(I,J) = 0.0D+00
        END DO
      END DO

      IF ( ND .EQ. 1 ) THEN
        DO I = 1, NI
          BK(I,1) = 1.0D+00
        END DO
        RETURN
      END IF

      DO I = 1, NI

        DO J = 1, ND

          IF ( J .EQ. 1 .AND. XI(I) .LE. XD(J) ) THEN

            T = ( XI(I) - XD(J) ) / ( XD(J+1) - XD(J) )
            BK(I,J) = 1.0D+00 - T

          ELSE IF ( J .EQ. ND .AND. XD(J) .LE. XI(I) ) THEN

            T = ( XI(I) - XD(J-1) ) / ( XD(J) - XD(J-1) )
            BK(I,J) = T

          ELSE IF ( XD(J-1) .LE. XI(I) .AND. XI(I) .LE. XD(J) ) THEN

            T = ( XI(I) - XD(J-1) ) / ( XD(J) - XD(J-1) )
            BK(I,J) = T

          ELSE IF ( XD(J) .LE. XI(I) .AND. XI(I) .LE. XD(J+1) ) THEN

            T = ( XI(I) - XD(J) ) / ( XD(J+1) - XD(J) )
            BK(I,J) = 1.0D+00 - T

          END IF

        END DO

      END DO

      RETURN
      END
      SUBROUTINE PWL_INTERP_1D ( ND, XD, YD, NI, XI, YI )

C*********************************************************************72
C
CC PWL_INTERP_1D EVALUATES THE PIECEWISE LINEAR INTERPOLANT.
C
C  DISCUSSION:
C
C    THE PIECEWISE LINEAR INTERPOLANT L(ND,XD,YD)(X) IS THE PIECEWISE
C    LINEAR FUNCTION WHICH INTERPOLATES THE DATA (XD(I),YD(I)) FOR I = 1
C    TO ND.
C
C  LICENSING:
C
C    THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.
C
C  MODIFIED:
C
C    22 SEPTEMBER 2012
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, INTEGER ND, THE NUMBER OF DATA POINTS.
C    ND MUST BE AT LEAST 1.
C
C    INPUT, DOUBLE PRECISION XD(ND), THE DATA POINTS.
C
C    INPUT, DOUBLE PRECISION YD(ND), THE DATA VALUES.
C
C    INPUT, INTEGER NI, THE NUMBER OF INTERPOLATION POINTS.
C
C    INPUT, DOUBLE PRECISION XI(NI), THE INTERPOLATION POINTS.
C
C    OUTPUT, DOUBLE PRECISION YI(NI), THE INTERPOLATED VALUES.
C
      IMPLICIT NONE

      INTEGER ND
      INTEGER NI

      INTEGER I
      INTEGER K
      DOUBLE PRECISION T
      DOUBLE PRECISION XD(ND)
      DOUBLE PRECISION YD(ND)
      DOUBLE PRECISION XI(NI)
      DOUBLE PRECISION YI(NI)

      DO I = 1, NI
        YI(I) = 0.0D+00
      END DO

      IF ( ND .EQ. 1 ) THEN
        DO I = 1, NI
          YI(I) = YD(1)
        END DO
        RETURN
      END IF

      DO I = 1, NI

        IF ( XI(I) .LE. XD(1) ) THEN

          T = ( XI(I) - XD(1) ) / ( XD(2) - XD(1) )
          YI(I) = ( 1.0D+00 - T ) * YD(1) + T * YD(2)

        ELSE IF ( XD(ND) .LE. XI(I) ) THEN

          T = ( XI(I) - XD(ND-1) ) / ( XD(ND) - XD(ND-1) )
          YI(I) = ( 1.0D+00 - T ) * YD(ND-1) + T * YD(ND)

        ELSE

          DO K = 2, ND

            IF ( XD(K-1) .LE. XI(I) .AND. XI(I) .LE. XD(K) ) THEN

              T = ( XI(I) - XD(K-1) ) / ( XD(K) - XD(K-1) )
              YI(I) = ( 1.0D+00 - T ) * YD(K-1) + T * YD(K)
              EXIT

            END IF

          END DO

        END IF

      END DO

      RETURN
      END
      SUBROUTINE PWL_VALUE_1D ( ND, XD, YD, NI, XI, YI, LEFT, RIGHT )

C*********************************************************************72
C
CC PWL_VALUE_1D EVALUATES THE PIECEWISE LINEAR INTERPOLANT.
C
C  DISCUSSION:
C
C    THE PIECEWISE LINEAR INTERPOLANT L(ND,XD,YD)(X) IS THE PIECEWISE
C    LINEAR FUNCTION WHICH INTERPOLATES THE DATA (XD(I),YD(I)) FOR I = 1
C    TO ND.
C
C  LICENSING:
C
C    THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.
C
C  MODIFIED:
C
C    04 OCTOBER 2017
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, INTEGER ND, THE NUMBER OF DATA POINTS.
C    ND MUST BE AT LEAST 1.
C
C    INPUT, DOUBLE PRECISION XD(ND), THE DATA POINTS.
C
C    INPUT, DOUBLE PRECISION YD(ND), THE DATA VALUES.
C
C    INPUT, INTEGER NI, THE NUMBER OF INTERPOLATION POINTS.
C
C    INPUT, DOUBLE PRECISION XI(NI), THE INTERPOLATION POINTS.
C
C    OUTPUT, DOUBLE PRECISION YI(NI), THE INTERPOLATED VALUES.
C
        IMPLICIT NONE

        INTEGER ND
        INTEGER NI

        INTEGER I
        INTEGER K
        DOUBLE PRECISION T
        DOUBLE PRECISION XD(ND)
        DOUBLE PRECISION YD(ND)
        DOUBLE PRECISION XI(NI)
        DOUBLE PRECISION YI(NI)
        DOUBLE PRECISION LEFT, RIGHT

        DO I = 1, NI
          YI(I) = 0.0D+00
        END DO

        IF ( ND .EQ. 1 ) THEN
          DO I = 1, NI
            YI(I) = YD(1)
          END DO
          RETURN
        END IF

        DO I = 1, NI

          IF ( XI(I) .EQ. XD(1) ) THEN

            T = ( XI(I) - XD(1) ) / ( XD(2) - XD(1) )
            YI(I) = ( 1.0D+00 - T ) * YD(1) + T * YD(2)

          ELSE IF ( XD(ND) .LE. XI(I))  THEN

C             T = ( XI(I) - XD(ND-1) ) / ( XD(ND) - XD(ND-1) )
C             YI(I) = ( 1.0D+00 - T ) * YD(ND-1) + T * YD(ND)
            YI(I) = RIGHT
          ELSE IF ( XI(I) .LE. XD(1))  THEN

C             T = ( XI(I) - XD(ND-1) ) / ( XD(ND) - XD(ND-1) )
C             YI(I) = ( 1.0D+00 - T ) * YD(ND-1) + T * YD(ND)
            YI(I) = LEFT
          ELSE

            DO K = 2, ND

              IF ( XD(K-1) .LE. XI(I) .AND. XI(I) .LE. XD(K) ) THEN

                T = ( XI(I) - XD(K-1) ) / ( XD(K) - XD(K-1) )
                YI(I) = ( 1.0D+00 - T ) * YD(K-1) + T * YD(K)
                GO TO 10

              END IF

            END DO

10          CONTINUE

          END IF

        END DO
  
        RETURN
      END