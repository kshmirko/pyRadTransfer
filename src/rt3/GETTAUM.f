      FUNCTION GETTAUM(WL)
      REAL*8 WL, GETTAUM
CF2PY INTENT(IN)  WL
      GETTAUM = 0.008569D0/(WL**4.0)*
     .          (1.0D0+0.0113D0/(WL**2.0D0)+0.00013D0/(WL**4.0D0))
      RETURN
      END