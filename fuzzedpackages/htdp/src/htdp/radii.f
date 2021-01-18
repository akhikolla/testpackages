*********************************************************************
      SUBROUTINE RADII(YLAT,RADMER,RADPAR)
C
C  Computes the radius of curvature in the meridian
C  and the radius of curvature in a parallel of latitude
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COSLAT = DCOS(YLAT)
      DENOM = DSQRT(1.D0 + EPS*COSLAT*COSLAT)
      RADMER = AF/(DENOM**3)
      RADPAR = AF*COSLAT/DENOM
      RETURN
      END
