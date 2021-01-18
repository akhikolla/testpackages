
C***************************************************

      INTEGER FUNCTION IUNGRD(IREGN, I, J, IVEC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NUMGRD = 10)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),
     1          GRDLY(NUMGRD), GRDUY(NUMGRD),
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)

      IUNGRD = NBASE(IREGN) +
     1      3 * ((J - 1) * (ICNTX(IREGN) + 1) +  (I - 1)) + IVEC

      RETURN
      END
