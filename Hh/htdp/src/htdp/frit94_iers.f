*************************************************************************************
      subroutine frit94_IERS (x1, y1, z1, x2, y2, z2, date, jopt)

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

      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)


      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc1(iopt)
      tranx = tx1(iopt) + dtx1(iopt)*dtime
      trany = ty1(iopt) + dty1(iopt)*dtime
      tranz = tz1(iopt) + dtz1(iopt)*dtime
      rotnx  = rx1(iopt) + drx1(iopt)*dtime
      rotny  = ry1(iopt) + dry1(iopt)*dtime
      rotnz  = rz1(iopt) + drz1(iopt)*dtime
      ds     = 1.d0 + scale1(iopt) + dscale1(iopt)*dtime

      x2 = tranx + ds*x1    + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1    + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1

      return
      end
