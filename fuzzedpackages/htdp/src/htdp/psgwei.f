*******************************************************************
      SUBROUTINE PSGWEI (POSX, POSY, K, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
C              AND COMPUTES NORMALIZED WEIGHTS FOR 
C              BI-LINEAR INTERPOLATION OVER A PLANE
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C POSX         LONGITUDE OF POINT IN DEGREES, POSITIVE EAST
C POSY         LATITUDE OF POINT IN DEGREES, POSITIVE NORTH
C K            ID OF EARTHQUAKE GRID                   
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
C    THIS MODULE CALLED BY:   PSDISP
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      /PSGRID/, /CONST/
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
      parameter (NUMPSG = 1)
      DIMENSION WEI(2,2)
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG), 
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C*** Obtain indices for the lower-left corner of the cell
C*** containing the point
      STEPX = (PSGUX(K) - PSGLX(K)) / ICNTPX(K)
      STEPY = (PSGUY(K) - PSGLY(K)) / ICNTPY(K)
      I = IDINT((POSX - PSGLX(K))/STEPX) + 1
      J = IDINT((POSY - PSGLY(K))/STEPY) + 1
c     write(6,1001) K, I, J
c1001 format(1x, 'quake = ', I5 /
c    1       1x, ' i = ', I5 /
c    1       1x, ' j = ', I5)

C*** Compute the limits of the grid cell 
      GRLX = PSGLX(K) + (I - 1) * STEPX
      GRUX = GRLX + STEPX                    
      GRLY = PSGLY(K) + (J - 1) * STEPY                
      GRUY = GRLY + STEPY                     

C*** Compute the normalized weights for the point               
      DENOM = (GRUX - GRLX) * (GRUY - GRLY)
      WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
      WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
      WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
      WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM

      RETURN
      END
