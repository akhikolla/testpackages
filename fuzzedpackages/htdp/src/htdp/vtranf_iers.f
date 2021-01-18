******************************************************
      SUBROUTINE VTRANF_IERS(X,Y,Z,VX,VY,VZ, IOPT1, IOPT2)

*** Convert velocity from reference frame of IOPT1 to 
*** reference frame of IOPT2.
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (numref = 15)
      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)

      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      IF(IOPT1 .le. numref .and. IOPT2. le. numref
     &   .and. IOPT1 .gt. 0 .and. IOPT2 .gt. 0 ) THEN

*** Convert from mm/yr to m/yr
         VX = VX /1000.d0
         VY = VY / 1000.d0
         VZ = VZ / 1000.d0

*** From IOPT1 to ITRF94 
*** (following equations use approximations assuming
*** that rotations and scale change are small)
         WX = -drx1(iopt1)           
         WY = -dry1(iopt1)      
         WZ = -drz1(iopt1)      
         DS = -dscale1(iopt1)
         VX = VX - dtx1(iopt1) + DS*X + WZ*Y - WY*Z
         VY = VY - dty1(iopt1) - WZ*X  +DS*Y + WX*Z
         VZ = VZ - dtz1(iopt1) + WY*X - WX*Y + DS*Z

*** From ITRF94 to IOPT2 reference frame
*** (following equations use approximations assuming
***  that rotations and scale change are small)
         WX = drx1(iopt2)
         WY = dry1(iopt2)
         WZ = drz1(iopt2)
         DS = dscale1(iopt2)
         VX = VX + dtx1(iopt2) + DS*X + WZ*Y - WY*Z
         VY = VY + dty1(iopt2) - WZ*X + DS*Y + WX*Z 
         VZ = VZ + dtz1(iopt2) + WY*X - WX*Y + DS*Z

*** FROM m/yr to mm/yr
         VX = VX * 1000.d0
         VY = VY * 1000.d0
         VZ = VZ * 1000.d0

      ELSE
C        write(luout,*) ' Improper reference frame in routine vtranf'
C        stop
         CALL REXIT('Improper reference frame in routine VTRANF_IERS')
      ENDIF

      RETURN
      END
