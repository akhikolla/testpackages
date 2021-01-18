
*******************************************************************
      SUBROUTINE GETBDY

*** Obtain coordinates for vertices that form the polygons
*** that correspond to the boundaries for the regions.       
*** Region 1 is the San Andreas fault in central California         
*** Region 2 is southern California
*** Region 3 is Northern California
*** Region 4 is the Pacific Noerthwest
*** Region 5 is western CONUS
*** Region 6 is CONUS
*** Region 7 is St. Elias, Alaska
*** Region 8 is south-central Alaska
*** Region 9 is southeast Alaska
*** Region 10 is All Mainland Alaska
*** Region 11 is the North American plate 
*** Region 12 is the Caribbean plate
*** Region 13 is the Pacific plate
*** Region 14 is the Juan de Fuca plate
*** Region 15 is the Cocos plate
*** Region 16 is the Mariana plate
*** REGION 17 is the Philippine Sea plate

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NMREGN = 17)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /BNDRY/ X(4000), Y(4000), NPOINT(30)
 
      IEND = NPOINT(NMREGN + 1) - 1  
         DO 10 J = 1, IEND              
           X(J) = (X(J) * 3600.D0)/RHOSEC
           Y(J) = (Y(J) * 3600.D0)/RHOSEC
   10    CONTINUE
      RETURN
      END
