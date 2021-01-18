************************************
      SUBROUTINE IYMDMJ( IYR, IMON, IDAY, MJD )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:       IYMDMJ
C VERSION:    Sep. 17, 2010
C WRITTEN BY: R. SNAY (after M. SCHENEWERK)
C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE 
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IDAY              DAY
C IMON              MONTH
C IYR               YEAR
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C MJD               MODIFIED JULIAN DATE 
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C A                 TEMPORARY STORAGE
C B                 TEMPORARY STORAGE
C C                 TEMPORARY STORAGE
C D                 TEMPORARY STORAGE
C IMOP              TEMPORARY STORAGE
C IYRP              TEMPORARY STORAGE
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY: GENERAL USE
C
C       THIS MODULE CALLS:     DINT
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:       
C
C       REFERENCES:            DUFFETT-SMITH, PETER  1982, 'PRACTICAL
C                              ASTRONOMY WITH YOUR CALCULATOR', 2ND
C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                              NEW YORK, P.9
C
C       COMMENTS:              THIS SUBROUTINE REQUIRES THE FULL YEAR,
C                              I.E. 1992 RATHER THAN 92.  
C
C********1*********2*********3*********4*********5*********6*********7**
C::LAST MODIFICATION
C::8909.06, MSS, DOC STANDARD IMPLIMENTED
C::9004.17, MSS, CHANGE ORDER YY MM DD
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(4) (I-N)
C
      INTEGER(4)     A, B, C, D

      IYRP = IYR
C
C........  0.0  EXPLICIT INITIALIZATION
C
      IF( IMON .LT. 3 ) THEN
        IYRP= IYRP - 1
        IMOP= IMON + 12
      ELSE
        IMOP= IMON
      END IF
C
C........  1.0  CALCULATION
C
      A=  IYRP*0.01D0
      B=  2 - A + DINT( A*0.25D0 )
      C=  365.25D0*IYRP
      D=  30.6001D0*(IMOP + 1)
      MJD =  (B + C + D + IDAY - 679006) 
C      
      RETURN
      END
