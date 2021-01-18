
C*********************************************************************
      SUBROUTINE GRDAMP (K, I, J, AMP, PS)
C********1*********2*********3*********4*********5*********6*********7**
C
C PURPOSE:     THIS SUBROUTINE RETRIEVES THE AMPLITUDES OF THE FOUR
C              GRID NODES OFGRID K WHERE I,J ARE THE INDICES OF
C              THE LOWER LEFT HAND CORNER
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C
C K            ID OF EARTHQUAKE CORRESPONDING TO GRID
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C PS           THE ARRAY CONTAINING ALL THE GRIDDED AMPLITUDES
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C AMP          A TWO BY TWO ARRAY CONTAINING THE 3D AMPLITUDES
C              FOR THE CORNERS OF THE GRID
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
C    COMMON BLOCKS USED:      NONE     
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::2011.08.17, RAS, ORIGINAL CREATION FOR TRANS4D
C********1*********2*********3*********4*********5*********6*********7**


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      DIMENSION AMP(2,2,3), PS(*)

      DO 30 II = 0,1
         DO 20 IJ = 0,1
            DO 10 IVEC = 1, 3
               INDEX = IPSGRD(K, I + II, J + IJ, IVEC)
               AMP(II + 1, IJ + 1, IVEC) = PS(INDEX)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE   

      RETURN
      END
