***************************************************************
      SUBROUTINE NEWCOR(YLAT,YLON,HTOLD,MIN1,MIN2,
     1         YLAT3,YLON3,HTNEW,DN,DE,DU,VN, VE, VU)

*** Predict coordinates at time MIN2 given coordinates at time MIN1.
*** Predict displacements from time MIN1 to time MIN2.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)

      HT = HTOLD

      CALL COMPSN(YLAT1,YLON1,HT1,YLAT,YLON,HT,
     1            MIN1,VN, VE, VU)

      CALL COMPSN(YLAT2,YLON2,HT2,YLAT,YLON,HT,
     1            MIN2, VN, VE, VU)

      YLAT3 = YLAT + YLAT2 - YLAT1
      YLON3 = YLON + YLON2 - YLON1
      HTNEW   = HT   + HT2   - HT1

      CALL RADII(YLAT,RADMER,RADPAR)

      DN =  RADMER * (YLAT2 - YLAT1)
      DE = -RADPAR * (YLON2 - YLON1)
      DU =  HT2 - HT1

      RETURN
      END
