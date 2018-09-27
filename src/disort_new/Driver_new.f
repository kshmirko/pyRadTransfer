c     -------------------------------------------------------------------
c     Python wrapper to the DISORT radiative transfer solver
c
c     Author: Sebastian Gimeno Garcia
c
c
c     License:
c
c     Do whatever you want with this piece of code. Enjoy it. If you
c     find it helpful, think about the authors of DISORT and drink to
c     their health, and why not, also to mine.
c
c     If you find any bug, please let me now.
c      
c     Ref:
c     
c     K. Stamnes, SC. Tsay, W. Wiscombe and K. Jayaweera, Numerically
c     stable algorithm for discrete-ordinate-method radiative
c     transfer in multiple scattering and emitting layered media,
c     Appl Opt 27 (1988) (12), pp. 2502–2509.
c     -------------------------------------------------------------------
      
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Driver based on:
c      
c     RCS version control information:
c     $Header: DISOTEST.f,v 2.1 2000/04/03 21:21:55 laszlo Exp $
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE  RUN_1(
     I     MAXCLY, MAXMOM, MAXCMU,
     I     MAXUMU, MAXPHI, MAXULV,
     I     USRANG, USRTAU, IBCND, ONLYFL, PRNT,
     I     PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,
     I     DTAUC, SSALB, IPHAS, GG, TEMPER, WVNMLO, WVNMHI,
     I     UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
     I     FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
     I     EARTH_RADIUS, H_LYR,
     O     RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,
     O     ACCUR, HEADER,
     O     RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     O     ALBMED, TRNMED)
      IMPLICIT NONE
      INTEGER MAXCLY, MAXMOM, MAXCMU, MAXUMU, MAXPHI, MAXULV, IBCND
      LOGICAL USRANG, USRTAU, ONLYFL, PRNT(5), PLANK, LAMBER, DELTAMPLUS
     $     , DO_PSEUDO_SPHERE
      REAL    DTAUC( MAXCLY ), SSALB( MAXCLY ), PMOM( 0:MAXMOM, MAXCLY )
     $     , TEMPER( 0:MAXCLY ), WVNMLO, WVNMHI, UTAU( MAXULV ), UMU0,
     $     PHI0, UMU( MAXUMU ), PHI( MAXPHI ), FBEAM, FISOT, ALBEDO,
     $     BTEMP, TTEMP, TEMIS, EARTH_RADIUS, H_LYR( 0:MAXCLY ),
     $     RHOQ(MAXCMU/2, 0:MAXCMU/2, 0:(MAXCMU-1)),
     $     RHOU(MAXUMU, 0:MAXCMU/2, 0:(MAXCMU-1)),
     $     RHO_ACCURATE(MAXUMU, MAXPHI),
     $     EMUST(MAXUMU), BEMST(MAXCMU/2),
     $     ACCUR
      CHARACTER HEADER*127
      REAL RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     $     DFDT( MAXULV ), UAVG( MAXULV ), UU( MAXUMU, MAXULV, MAXPHI ),
     $     ALBMED( MAXUMU ), TRNMED( MAXUMU )
      
CF2PY INTENT(IN) MAXCLY, MAXMOM, MAXCMU, MAXUMU, MAXPHI, MAXULV
CF2PY INTENT(IN) USRANG, USRTAU, IBCND, ONLYFL, PRIN, PLANK,
CF2PY INTENT(IN) LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE, DTAUC,
CF2PY INTENT(IN) SSALB, PMOM, TEMPER, PHI, FBEAM, FISOT,
CF2PY INTENT(IN) WVNMLO, WVNMHI, UTAU, UMU0, PHI0, UMU, 
CF2PY INTENT(IN) ALBEDO, BTEMP, TTEMP, TEMIS, EARTH_RADIUS, H_LYR 
CF2PY INTENT(OUT) RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,
CF2PY INTENT(OUT) ACCUR, HEADER, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
CF2PY INTENT(OUT) UU, ALBMED, TRNMED
     
c    Runs test problems for DISORT and checks answers. These
c    problems test almost all logical branches in DISORT.

c    It is HIGHLY recommended that you use the code below as a template
c    for creating your own CALLs to DISORT, rather than starting from
c    scratch.  This will prevent a lot of mistakes and ensure that every
c    input argument gets a value.  Note in particular how GETMOM is
c    sometimes called to fill an array section of PMOM (for one layer);
c    several people have done this incorrectly in attempting to write it
c    ab initio (passing array sections for arrays that do not start at
c    element 1 is tricky).

c    Note that the ratio to the 'correct answer' may occasionally be
c    significantly different from unity -- even so different that
c    the ratio just prints as ****** rather than a number.  However,
c    this mostly occurs for values of flux or intensity that are very
c    small compared to the forcing functions (that is, small compared
c    to internal thermal emission and/or radiation incident at the
c    boundaries).  The printed number 'SERIOUSLY NON-UNIT RATIOS'
c    attempts to count just the cases where there is a real disagreement
c    and not those where quantitites are down at their noise level
c    (defined as 10^(-6) times their maximum value).

c    Further documentation can be found in the file DISOTEST.doc.


c  Routines called :

c    DISORT:   The discrete ordinates radiative transfer program

c    BDREF:    Sets bidirectional reflectance of lower boundary

c    GETMOM:   Sets phase function Legendre coefficients

c    PRTFIN:   Prints fluxes and intensities and their ratios to
c              the correct values

c    CHEKDO:   Data block containing correct fluxes and intensities

c    RATIO :   Ratio of calculated to correct value with underflow
c              and overflow protection (kept in file DISORT.f)

c       INPUT: IPHAS   Phase function options
c                      1 : Isotropic
c                      2 : Rayleigh
c                      3 : Henyey-Greenstein with asymmetry factor GG
c                      4 : Haze L as specified by Garcia/Siewert
c                      5 : Cloud C.1 as specified by Garcia/Siewert

c              GG      Asymmetry factor for Henyey-Greenstein case

c              NMOM    Index of highest Legendre coefficient needed
c                        ( = number of streams 'NSTR'  chosen
c                         for the discrete ordinate method)
c+---------------------------------------------------------------------+
c
c                 ** DISORT I/O specifications **

      INTEGER  IPHAS( MAXCLY )
      REAL     GG( MAXCLY )
      INTEGER  LC
c+---------------------------------------------------------------------+


c+---------------------------------------------------------------------+
c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX
c     ..

      DO LC = 1, MAXCLY
         CALL  GETMOM( IPHAS( LC ), GG( LC ), MAXMOM, PMOM(0,LC) )
      END DO 
      HEADER = 'Python wrapper to the DISORT radiative transfer solver'

      CALL DISORT( MAXCLY, MAXMOM, MAXCMU, 
     &     MAXUMU, MAXPHI, MAXULV,
     &     USRANG, USRTAU, IBCND, ONLYFL, PRNT,
     &     PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,
     &     DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,
     &     UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
     &     FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
     &     EARTH_RADIUS, H_LYR, 
     &     RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,
     &     ACCUR,  HEADER,
     &     RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &     ALBMED, TRNMED )
      
      RETURN
      END   
