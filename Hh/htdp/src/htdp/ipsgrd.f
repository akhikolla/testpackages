
C***************************************************
      INTEGER FUNCTION IPSGRD(IGRID, I, J, IVEC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NUMPSG = 1)
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)

      IPSGRD = NBASEP(IGRID) +
     1      3 * ((J - 1) * (ICNTPX(IGRID) + 1) +  (I - 1)) + IVEC

      RETURN
      END
