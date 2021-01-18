***************************************************
      SUBROUTINE TOVXYZ(GLAT,GLON,VN,VE,VU,VX,VY,VZ)
*** Convert velocities from vn,ve,vu to vx,vy,vz
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      VX = -SLAT*CLON*VN - SLON*VE + CLAT*CLON*VU
      VY = -SLAT*SLON*VN + CLON*VE + CLAT*SLON*VU
      VZ =  CLAT*VN + SLAT*VU

      RETURN
      END
