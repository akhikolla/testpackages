*****************************************************************
      SUBROUTINE RADR8T (YLAT,VN,VE,VNR,VER)

C Convert horizontal velocities from mm/yr to rad/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)

      CALL RADII (YLAT,RADMER,RADPAR)
      VNR = VN / (1000.D0 * RADMER)
      VER = VE / (1000.D0 * RADPAR)

      RETURN
      END
