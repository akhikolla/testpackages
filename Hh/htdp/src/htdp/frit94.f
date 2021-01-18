*************************************************************
      subroutine frit94(x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts ITRF94 cartesian coordinates to cartesian
*** coordinates in the specified reference frame for the
*** given date
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (x1, y1, z1) --> input ITRF94 coordiates (meters)
*** (x2, y2, z2) --> output coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of output reference frame

      implicit double precision (a-h, o-z)
      implicit INTEGER(4) (i-n)
      parameter (numref = 15)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc(iopt)
      tranx = tx(iopt) + dtx(iopt)*dtime
      trany = ty(iopt) + dty(iopt)*dtime
      tranz = tz(iopt) + dtz(iopt)*dtime
      rotnx  = rx(iopt) + drx(iopt)*dtime
      rotny  = ry(iopt) + dry(iopt)*dtime
      rotnz  = rz(iopt) + drz(iopt)*dtime
      ds     = 1.d0 + scale(iopt) + dscale(iopt)*dtime

      x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1

      return
      end
