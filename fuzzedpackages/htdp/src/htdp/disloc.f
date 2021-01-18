*************************************************************
      SUBROUTINE DISLOC (YLAT,YLON,STRIKE,HL,EQLAT,EQLON,
     &          SS,DS,DIP,DEPTH,WIDTH,DNORTH,DWEST,DUP)

*** Compute 3-dimensional earthquake displacement at point
*** using dislocation theory
*
*   INPUT:
*         YLAT = Latitude in radians (positive north)
*         YLON = Longitude in radians (positive west)
*         STRIKE = strike in radians clockwise from north such
*                  that the direction of dip is pi/2 radians
*                  counterclockwise from the direction of strike
*         HL   = Half-length in meters
*         EQLAT = Latitude in radians of midpoint of the
*                 rectangle's upper edge (positive north)
*         EQLON = Longitude in radians of midpoint of the
*                 rectangle's upper edge (positive west)
*         SS = strike slip in meters (positive = right lateral)
*         DS = dip slip  in meters (positive = normal faulting)
*         DIP = dip in radians
*         DEPTH = Vertical depth of rectangle's upper edge
*                 in meters
*         WIDTH = width of rectangle in meters
*
*   OUTPUT:
*         DNORTH = northward displacement in radians
*         DWEST = westward displacement in radians
*         DUP = upward displacement in meters
************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)

*** Compute radii of curvature at fault center
      CALL RADII (EQLAT, RMER, RPAR)

*** Compute planar coordinates in meters
      DLAT = (YLAT - EQLAT) * RMER
      DLON = (YLON - EQLON) * RPAR
      COSSTR = DCOS(STRIKE)
      SINSTR = DSIN(STRIKE)
      X1 = COSSTR*DLAT - SINSTR*DLON
      X2 = SINSTR*DLAT + COSSTR*DLON

*** Compute displacements in fault-oriented coordinates
      CALL OKADA(X1,X2,HL,DEPTH,WIDTH,DIP,U1SS,U2SS,
     &      U3SS,U1DS,U2DS,U3DS)
      U1 = U1SS*SS + U1DS*DS
      U2 = U2SS*SS + U2DS*DS
      DUP = U3SS*SS + U3DS*DS

*** Convert horizontal displacements to radians
*** in north-west coordinate system
      DNORTH = ( COSSTR*U1 + SINSTR*U2) / RMER
      DWEST  = (-SINSTR*U1 + COSSTR*U2) / RPAR

      RETURN
      END
