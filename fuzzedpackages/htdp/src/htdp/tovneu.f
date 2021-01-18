*****************************************************
      SUBROUTINE TOVNEU(GLAT,GLON,VX,VY,VZ,VN,VE,VU)

*** Convert velocities from vx,vy,vz to vn,ve,vu

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      VN = -SLAT*CLON*VX - SLAT*SLON*VY + CLAT*VZ
      VE = -SLON*VX + CLON*VY
      VU = CLAT*CLON*VX + CLAT*SLON*VY + SLAT*VZ

      RETURN
      END
