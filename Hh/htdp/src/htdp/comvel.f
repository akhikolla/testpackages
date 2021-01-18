*****************************************************************
      SUBROUTINE COMVEL(YLAT,YLON,JREGN,VN,VE,VU)
C
C Compute the NAD_83(CORS96) velocity at a point in mm/yr    !Not anymore since 09/12/2014
                                                             !Now the velocity refer to ITRF2008
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NUMGRD = 10)
      parameter (NMREGN = 17)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1,I2,I3,I4,I5,I6
      COMMON /VGRID/ B(210000)
      DIMENSION WEI(2,2), VEL(2,2,3)

c     WRITE (6, 1001) JREGN
c1001 FORMAT( 'JREGN = ', I6)

      IF(JREGN .GT. NUMGRD .AND. JREGN .LE. NMREGN) THEN
*** Use tectonic plate model to compute velocity relative
***    to ITRF2008
        IPLATE = JREGN - NUMGRD
        ELON = - YLON
        HT = 0.D0
        CALL TOXYZ(YLAT, ELON, HT, X, Y, Z)
        CALL PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
        VX = VX * 1000.D0
        VY = VY * 1000.D0
        VZ = VZ * 1000.D0
*** Convert ITRF2008 velocity to NAD_83(CORS96) velocity
c       CALL VTRANF(X, Y, Z, VX, VY, VZ, 15, 1)          !No Do not. Leave the velocity in ITRF2008
        CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)

      ELSEIF(JREGN .GE. 1 .AND. JREGN .LE. NUMGRD) THEN

C*** Get indices for the lower left hand corner of the grid
C*** and get the weights for the four corners
        CALL GRDWEI (YLON, YLAT, JREGN, I, J, WEI)

C*** Get the velocity vectors at the four corners
        CALL GRDVEC (JREGN, I, J, VEL, B)

        VN = WEI(1,1) * VEL(1,1,1) + WEI(1,2) * VEL(1,2,1)
     *     + WEI(2,1) * VEL(2,1,1) + WEI(2,2) * VEL(2,2,1)

        VE = WEI(1,1) * VEL(1,1,2) + WEI(1,2) * VEL(1,2,2)
     *     + WEI(2,1) * VEL(2,1,2) + WEI(2,2) * VEL(2,2,2)
  
        VU = WEI(1,1) * VEL(1,1,3) + WEI(1,2) * VEL(1,2,3)
     *     + WEI(2,1) * VEL(2,1,3) + WEI(2,2) * VEL(2,2,3)

C*** If the point in one of the four Alaskan regions,
C*** then set its vertical velocity to 0.0
        IF(JREGN .GE. 7 .AND. JREGN .LE. 10) THEN
           VU = 0.D0
        ENDIF

C*** If the point is in one of the first ten regions, then
c*** the velocity grids contain the ITRF2008 velocity.
c*** Hence, the following code transforms this ITRF2008 velocity
c*** to the corresponding NAD 83 (CORS96) velocity. 
C
c       IF(JREGN .LE. NUMGRD) THEN                            !Starting09/12/2014, the velocities
c          ELON = - YLON                                      !coming out of this routine are in ITRF2008
c          HT = 0.D0
c          CALL TOXYZ(YLAT, ELON, HT, X, Y, Z)
c          CALL TOVXYZ(YLAT, ELON, VN, VE, VU, VX, VY, VZ)
c          CALL VTRANF( X, Y, Z, VX, VY, VZ, 15, 1)
c          CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)
c       ENDIF

      ELSE
C       WRITE(LUOUT,100) JREGN
C 100   FORMAT(' Improper region identifier ',I4,'in COMVEL.')
C       STOP
        CALL REXIT('Improper region identifier in routine COMVEL')
      ENDIF
      RETURN
      END
