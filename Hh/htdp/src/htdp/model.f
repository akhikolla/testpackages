*****************************************************************************
      SUBROUTINE MODEL

*** Obtain parameters defining crustal motion model

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF

      A = 6.378137D06
      F = 1.D0 / 298.257222101D0
      E2 = 0.6694380022903146D-2
      AF = A / (1.D0 -F)
      EPS = F*(2.D0 - F) / ((1.D0 -F)**2)
      PI = 4.D0 * DATAN(1.D0)
      RHOSEC = (180.D0 * 3600.D0) / PI
      TWOPI = PI + PI

C*** Set default reference epoch to Jan. 1, 2010
      IYRREF = 2010
      IMOREF = 1
      IDYREF = 1
      CALL IYMDMJ (IYRREF, IMOREF, IDYREF, MJD)
      ITREF = MJD * 24 * 60

      CALL GETBDY

      RETURN
      END
