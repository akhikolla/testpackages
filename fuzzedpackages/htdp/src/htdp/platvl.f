****************************************
      SUBROUTINE PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
*** Compute the ITRF2008 velocity at point on plate = IPLATE
***    with coordinates X, Y, Z (in meters)
***    The resulting velocities--VX, VY, and VZ--will be in meters/yr
***    References 
***     Altamimi et al. 2012 = JGR (Paper on ITRF2008 plate motion)
***     DeMets et al. 2010 = Geophysical Journal Int'l, vol 181, 
***     Snay 2003 = SALIS, Vol 63, No 1 (Paper on Frames for Pacific)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER(4) (I-N)
      DIMENSION WX(7), WY(7), WZ(7)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

*** IPLATE = 1 --> North America (from Altamimi et al. 2012)
***          2 --> Caribbean (from Altamimi et al. 2012)
***          3 --> Pacific (from Altamimi et al. 2012)
***          4 --> Juan de Fuca (from DeMets et al. 2010)
***          5 --> Cocos (from DeMets et al. 2010)
***          6 --> Mariana (from Snay, 2003)
***          7 --> Philippine Sea (from DeMets et al. 2010)
      DATA WX /0.170D-9, 0.238D-9, -1.993D-9, 
     1         6.626D-9, -10.390D-9, -.097D-9, -0.841D-9/
      DATA WY /-3.209D-9, -5.275D-9, 5.023D-9,
     1         11.708D-9, -14.954D-9,  .509D-9, 3.989D-9/
      DATA WZ /-0.485D-9, 3.219D-9,-10.501D-9, 
     1        -10.615D-9,  9.148D-9,-1.682D-9, -10.626D-9/

      IF (IPLATE .LE. 0 .OR. IPLATE .GT. 7) THEN
C         WRITE (LUOUT, 1) IPLATE
C   1     FORMAT(' Improper plate ID in PLATVL = ', I6)
C         STOP
          CALL REXIT('Improper plate ID in routine PLATVL')
      ENDIF

      VX = -WZ(IPLATE) * Y + WY(IPLATE) * Z
      VY =  WZ(IPLATE) * X - WX(IPLATE) * Z
      VZ = -WY(IPLATE) * X + WX(IPLATE) * Y

*** The parameters--WX, WY, and WZ--refer to ITRF2000
*** for the Mariana Plate (Snay, 2003). Hence,
*** for this plate, VX, VY, and VZ, correspond to ITRF2000.
*** The following code converts these to ITRF2008 velocities for
*** this plate.
      IF (IPLATE .EQ. 6) THEN
         VX = VX*1000.d0
         VY = VY*1000.d0
         VZ = VZ*1000.d0
         CALL VTRANF(X, Y, Z, VX, VY, VZ, 11, 15)
         VX = VX/1000.d0
         VY = VY/1000.d0
         VZ = VZ/1000.d0
*** The following translations rates are added per Altamimi et al. (2012)
*** for the other six plates
      ELSE
         VX = 0.00041d0 + VX
         VY = 0.00022d0 + VY
         VZ = 0.00041d0 + VZ
      ENDIF

      RETURN
      END
