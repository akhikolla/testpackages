*******************************************************************
      SUBROUTINE GETREG(X0,YKEEP,JREGN)

*** Determine the region containing a given point.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NMREGN = 17)
      COMMON /BNDRY/ X(4000), Y(4000), NPOINT(30)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC

      Y0 = TWOPI - YKEEP
      IF (Y0 .lt. 0.d0) Y0 = Y0 + TWOPI
      IR = 0     
    1 IR = IR + 1
      IF(IR .GT. NMREGN) THEN
         JREGN = 0
         RETURN
      ENDIF
      IBEGIN = NPOINT(IR)
      NUMVER = NPOINT(IR + 1) - IBEGIN
      CALL POLYIN(X0,Y0,X(IBEGIN),Y(IBEGIN), NUMVER, NTEST)
      IF(NTEST .EQ. 0) GO TO 1
      JREGN = IR       

      RETURN
      END
