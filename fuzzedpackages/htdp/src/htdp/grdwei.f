*******************************************************************

C     SUBROUTINE GRDCHK (POSX, POSY, INSIDE)

C
C ROUTINE CHECKS IF THE POINT HAVING COORDINATES (POSX, POSY)
C IS WITHIN THE REGION SPANNED BY THE GRID
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     IMPLICIT INTEGER(4) (I-N)
C     LOGICAL INSIDE
C     COMMON /CDGRID/ GRDLX, GRDUX, GRDLY, GRDUY, ICNTX, ICNTY

C     INSIDE = .TRUE.

C     IF (POSX .LT. GRDLX .OR. POSX .GT. GRDUX) THEN
C        INSIDE = .FALSE.
C     ENDIF
C     IF (POSY .LT. GRDLY .OR. POSY .GT. GRDUY) THEN
C        INSIDE = .FALSE.
C     ENDIF

C     RETURN
C     END

*******************************************************************

      SUBROUTINE GRDWEI (YLON, YLAT, JREGN, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDWEI
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
C              AND COMPUTES NORMALIZED WEIGHTS FOR 
C              BI-LINEAR INTERPOLATION OVER A PLANE
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C YLON         LONGITUDE OF POINT IN RADIANS, POSITIVE WEST
C YLAT         LATITUDE OF POINT IN RADIANS, POSITIVE NORTH
C JREGN        ID OF GEOGRAPHIC REGION CONTAINING POINT
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C WEI          A TWO BY TWO ARRAY CONTAINING THE NORMALIZED WEIGHTS
C              FOR THE CORNER VECTORS
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   COMVEL
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      /CDGRID/, /CONST/
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9511.09, RAS, MODIFIED FOR HTDP
C::9712.05, RAS, MODIFIED TO ACCOUNT FOR MULTIPLE GRIDS
C********1*********2*********3*********4*********5*********6*********7**
    
C**** COMPUTES THE WEIGHTS FOR AN ELEMENT IN A GRID

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NUMGRD = 10)
      DIMENSION WEI(2,2)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD), 
     1          GRDLY(NUMGRD), GRDUY(NUMGRD),
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C*** Convert input coordinates to degrees
      POSX = (TWOPI - YLON) * 180.D0 / PI
      POSY = YLAT * 180.D0 / PI

C*** Obtain indices for the lower-left corner of the cell
C*** containing the point
      STEPX = (GRDUX(JREGN) - GRDLX(JREGN)) / ICNTX(JREGN)
      STEPY = (GRDUY(JREGN) - GRDLY(JREGN)) / ICNTY(JREGN)
      I = IDINT((POSX - GRDLX(JREGN))/STEPX) + 1
      J = IDINT((POSY - GRDLY(JREGN))/STEPY) + 1
c     write(6,1001) JREGN, I, J
c1001 format(1x, 'jregn = ', I5 /
c    1       1x, ' i = ', I5 /
c    1       1x, ' j = ', I5)

C*** Compute the limits of the grid cell 
      GRLX = GRDLX(JREGN) + (I - 1) * STEPX
      GRUX = GRLX + STEPX                    
      GRLY = GRDLY(JREGN) + (J - 1) * STEPY                
      GRUY = GRLY + STEPY                     

C*** Compute the normalized weights for the point               
      DENOM = (GRUX - GRLX) * (GRUY - GRLY)
      WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
      WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
      WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
      WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM

      RETURN
      END
