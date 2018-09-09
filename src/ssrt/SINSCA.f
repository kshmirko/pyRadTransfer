      REAL*8 FUNCTION  SINSCA( DITHER, LAYRU, NLYR, PHASE, OMEGA, TAU,
     &                       UMU, UMU0, UTAU, FBEAM, PI )

c        Calculates single-scattered intensity from EQS. STWL (65b,d,e)

c                I N P U T   V A R I A B L E S

c        DITHER   10 times machine precision
c
c        LAYRU    index of UTAU in multi-layered system
c
c        NLYR     number of sublayers
c
c        PHASE    phase functions of sublayers
c
c        OMEGA    single scattering albedos of sublayers
c
c        TAU      optical thicknesses of sublayers
c
c        UMU      cosine of emergent angle
c
c        UMU0     cosine of incident zenith angle
c
c        UTAU     user defined optical depth for output intensity
c
c        FBEAM   incident beam radiation at top
c
c        PI       3.1415...
c
c   Called by- INTCOR
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..
CF2PY INTENT(HIDE) NLYR
CF2PY INTENT(IN) LAYRU, NLYR, DITHER, FBEAM, PI, UMU0, UTAU
CF2PY INTENT(IN) OMEGA, PHASE, TAU, UMU

      INTEGER   LAYRU, NLYR
      REAL*8      DITHER, FBEAM, PI, UMU, UMU0, UTAU
c     ..
c     .. Array Arguments ..

      REAL*8      OMEGA( NLYR ), PHASE( NLYR ), TAU( 0:NLYR-1 )
c     ..
c     .. Local Scalars ..

      INTEGER   LYR
      REAL*8      EXP0, EXP1
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, EXP
c     ..


      SINSCA = 0.
      EXP0 = EXP( -UTAU/UMU0 )

      IF( ABS( UMU+UMU0 ).LE.DITHER ) THEN

c                                 ** Calculate downward intensity when
c                                 ** UMU=UMU0, Eq. STWL (65e)

         DO 10 LYR = 1, LAYRU - 1
            SINSCA = SINSCA + OMEGA( LYR ) * PHASE( LYR ) *
     &               ( TAU( LYR ) - TAU( LYR-1 ) )
   10    CONTINUE

         SINSCA = FBEAM / ( 4.*PI * UMU0 ) * EXP0 * ( SINSCA +
     &            OMEGA( LAYRU )*PHASE( LAYRU )*( UTAU-TAU(LAYRU-1) ) )

         RETURN

      END IF


      IF( UMU.GT.0. ) THEN
c                                 ** Upward intensity, Eq. STWL (65b)
         DO 20 LYR = LAYRU, NLYR

            EXP1 = EXP( -( ( TAU( LYR )-UTAU )/UMU + TAU( LYR )/UMU0 ) )
            SINSCA = SINSCA + OMEGA( LYR )*PHASE( LYR )*( EXP0 - EXP1 )
            EXP0 = EXP1

   20    CONTINUE

      ELSE
c                                 ** Downward intensity, Eq. STWL (65d)
         DO 30 LYR = LAYRU, 1, -1

            EXP1 = EXP( -( ( TAU(LYR-1)-UTAU )/UMU + TAU(LYR-1)/UMU0 ) )
            SINSCA = SINSCA + OMEGA( LYR )*PHASE( LYR )*( EXP0 - EXP1 )
            EXP0 = EXP1

   30    CONTINUE

      END IF

      SINSCA = FBEAM / ( 4.*PI * ( 1. + UMU/UMU0 ) ) * SINSCA


      RETURN
      END
