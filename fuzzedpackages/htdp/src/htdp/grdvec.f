
C*********************************************************************
C
      SUBROUTINE GRDVEC (JREGN, I, J, VEL, B)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDVEC
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     THIS SUBROUTINE RETRIEVES THE APPROXIMATE VALUES OF THE
C              GRID NODE VELOCITIES FOR GRID (I,J) 
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C JREGN        ID OF GEOGRAPHIC REGION CORRESPONDING TO GRID
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C B            THE ARRAY CONTAINING ALL THE APPROXIMATE VALUES
C              FOR THE ADJUSTMENT
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C VEL          A TWO BY TWO ARRAY CONTAINING THE VELOCITY VECTORS
C              FOR THE CORNERS OF THE GRID
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
C    COMMON BLOCKS USED:      NONE     
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9712.05, RAS, MODIFIED FOR HTDP (version 2.2)
C********1*********2*********3*********4*********5*********6*********7**


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      DIMENSION VEL(2,2,3), B(*)

      DO 30 II = 0,1
         DO 20 IJ = 0,1
            DO 10 IVEC = 1, 3
               INDEX = IUNGRD(JREGN, I + II, J + IJ, IVEC)
               VEL(II + 1, IJ + 1, IVEC) = B(INDEX)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE   

      RETURN
      END
