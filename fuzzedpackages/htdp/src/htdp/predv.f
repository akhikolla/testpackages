******************************************************************
      subroutine PREDV(ylat, ylon, eht, date, iopt,
     1   jregn, vn, ve, vu) 

** Predict velocity in iopt reference frame       

** ylat       input - north latitude (radians)
** ylon       input - west longitude (radians)
** eht        input - ellipsoid height (meters)
** date       input - date (decimal years)
** iopt       input - reference frame
** jregn      output - deformation region
** vn         output - northward velocity in mm/yr
** ve         output - eastward velocity in mm/yr
** vu         output - upward velocity in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      logical  Is_iopt_NAD83

** Get reference latitude (RLAT) and reference longitude (RLON)

C  The following 2 lines were added on 07/22/2015 after Rich found this bug
         elon = -ylon
         call TOXYZ(ylat, elon, eht, x, y, z)

c        IF(IOPT .EQ. 0 .OR. IOPT .EQ. 1) THEN   !Velocity grids are No longer in NAD83
         IF(IOPT .EQ. 15) THEN                   !They are in ITRF2008
             RLAT = YLAT
             RLON = YLON
         ELSE
c            elon = -ylon                         !Was a bug, commented out on 07222015
c            call TOXYZ(ylat, elon, eht, x, y, z) !Was a bug, commented out on 07222015
c            CALL XTONAD(X,Y,Z,RLAT,RLON,EHTNAD,DATE,IOPT)   !No longer NAD83
             CALL XTO08 (X,Y,Z,RLAT,RLON,EHTNAD,DATE,IOPT)   !Positions should be in ITRF2008
         ENDIF

** Get deformation region

         CALL GETREG(RLAT,RLON,JREGN)
         IF (JREGN .EQ. 0) THEN
            VN = 0.D0
            VE = 0.D0
            VU = 0.D0
            RETURN
         ENDIF
         CALL COMVEL( RLAT, RLON, JREGN, VN, VE, VU)       !Those velocities are in ITRF2008

** Convert  velocity to reference of iopt, if iopt != NAD83   !No, was in ITRF2008, not in NAD83  (since 09/12/2014)

         Is_iopt_NAD83 = (iopt == 1)
c        IF (IOPT .NE. 0 .AND. IOPT .NE. 1) THEN
         IF (IOPT .NE. 15) THEN
            CALL TOVXYZ( YLAT, ELON, VN, VE, VU, VX, VY, VZ)
            if (Is_iopt_NAD83) then
c             CALL VTRANF( X, Y, Z, VX, VY, VZ, 1, IOPT)
              CALL VTRANF( X, Y, Z, VX, VY, VZ, 15, IOPT)
            else
              CALL VTRANF_IERS( X, Y, Z, VX, VY, VZ, 15, IOPT)
            endif
            CALL TOVNEU( YLAT, ELON, VX, VY, VZ, VN, VE, VU)
         ENDIF

         RETURN
         END
