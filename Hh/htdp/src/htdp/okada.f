****************************************************
      SUBROUTINE OKADA(X1,X2,XL,DU,W,DIP,
     1     U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)

************************************************************
*  This subroutine computes displacements at the point X1,X2
*  on the Earth's surface due to 1.0 meter of right-lateral 
*  strike slip (SS) and 1.0 meter of normal dip slip (DS) 
*  along a rectangular fault.
*
*  The rectangular fault dips in the direction of the positive
*  X2-axis.  The rectangle's strike parallels the X1-axis.
*  With the X3-axis directed upward out of the Earth, the X1-,
*  X2-, and X3-axes form a right-handed system.
*
*  The equations of dislocation theory are employed whereby
*  Earth is represented an a homogeneous, isotropic half-space
*  with a Poisson ratio of PNU.
*
*  REFERENCE: Okada, Y., Surface deformation due to shear and
*    tensile faults in a half-space, Bulletin of the 
*    Seismological Society of America, vol. 75, pp. 1135-1154 (1985)
*
*  The X3 = 0 plane corresponds to the Earth's surface. The plane's
*  origin is located directly above the midpoint of the rectangle's
*  upper edge.
*
*  INPUT:
*    X1,X2 - Location in meters
*    XL    - Rectangle's half-length in meters
*    DU    - Vertical depth to rectangle's upper edge in meters
*            (always positive or zero)
*    W     - Rectangle's width in meters
*    DIP   - Rectangle's dip in radians (always between 0 and PI/2
*
*  OUTPUT
*    U1SS  - Displacement in X1-direction due to 1.0 meters
*            of right-lateral strike slip
*    U2SS  - Displacement in X2-direction due to 1.0 meters
*            of right-lateral strike slip
*    U3SS  - Displacement in X3-direction due to 1.0 meters
*            of right-lateral strike slip
*    U1DS  - Displacement in X1-direction due to 1.0 meters
*            of normal dip slip
*    U2DS  - Displacement in X2-direction due to 1.0 meters
*            of normal dip slip
*    U3DS  - Displacement in X3-direction due to 1.0 meters
*            of normal dip slip
*******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      LOGICAL VERT   
      PI = 3.141593D0
      TWOPI = PI + PI
      PNU = 0.25D0
      RATIO = 1.D0 - 2.D0*PNU

      IF(DABS(PI/2.D0 - DIP) .LT. .01D0)THEN
               DIPK = -PI/2.D0
               VERT = .TRUE.
      ELSE
               DIPK = -DIP
               VERT = .FALSE.
      ENDIF

      SDIP = DSIN(DIPK)
      CDIP = DCOS(DIPK)
      P = X2*CDIP + DU*SDIP
      Q = X2*SDIP - DU*CDIP

      PSI = X1 + XL
      ETA = P
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1    VERT,U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)

      PSI = X1 + XL
      ETA = P - W
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS - C1SS
      U2SS = U2SS - C2SS
      U3SS = U3SS - C3SS
      U1DS = U1DS - C1DS 
      U2DS = U2DS - C2DS
      U3DS = U3DS - C3DS

      PSI = X1 - XL
      ETA = P
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS - C1SS
      U2SS = U2SS - C2SS
      U3SS = U3SS - C3SS
      U1DS = U1DS - C1DS
      U2DS = U2DS - C2DS
      U3DS = U3DS - C3DS

      PSI = X1 - XL
      ETA = P - W
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS + C1SS
      U2SS = U2SS + C2SS
      U3SS = U3SS + C3SS
      U1DS = U1DS + C1DS
      U2DS = U2DS + C2DS
      U3DS = U3DS + C3DS
      RETURN
      END
