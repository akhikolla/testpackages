****************************************************************************
      SUBROUTINE XTO08 (X,Y,Z,RLAT,WLON,EHT08,DATE,IOPT)

*** Converts X,Y,Z in specified datum to latitude and
*** longitude (in radians) and height (meters) in ITRF2008
*** datum with longitude positive west.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      LOGICAL FRMXYZ
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

*** Convert to cartesian coordinates in ITRF2008   
c     if (iopt .eq. 0 .or. iopt .eq. 1) then
      if (iopt .eq. 15) then
         x2 = x
         y2 = y
         z2 = z
      elseif (iopt .eq. 1) then
         call toit94(x,y,z,x1,y1,z1,date,iopt)
         call frit94(x1,y1,z1,x2,y2,z2,date,15)
      else 
         call toit94_iers(x,y,z,x1,y1,z1,date,iopt)
         call frit94_iers(x1,y1,z1,x2,y2,z2,date,15)
      endif

***Convert to geodetic coordinates
C     IF(.NOT.FRMXYZ(X2,Y2,Z2,RLAT,ELON,EHT08))STOP 666
      IF(.NOT.FRMXYZ(X2,Y2,Z2,RLAT,ELON,EHT08)) THEN
         CALL REXIT('Failed to converge in FRMXYZ')
      ENDIF
      WLON = -ELON
 100  IF(WLON .LT. 0.D0) THEN
          WLON = WLON + TWOPI
          GO TO 100
      ENDIF
      RETURN
      END
