*****************************************************
      SUBROUTINE PSDISP(YLAT, YLON, MIN, DNORTH, DEAST, DUP)
********
*   Compute total postseismic displacement for all earthquakes
*
* INPUT
*   YLAT       latitude of point in radians, positive north
*   YLON       longitude of point in radians, positive west
*   MIN        modified julian date of reference epoch for new coordinates
*               in minutes
*
*   DNORTH     Total northward postseismic displacement at point during
*              period from ITREF to MIN in meters
*   DEAST      TOTAL eastward postseismic displacement
*   DUP        Total upward postseismic displacement 
*******

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER(4) (I-N)
      LOGICAL INSIDE

      parameter (NUMPSG = 1)
      COMMON /CONST/ A, F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)
      COMMON /PGRID/ PS(18000)
      DIMENSION ITEQ(NUMPSG)      
      DIMENSION TAU(NUMPSG)   
      DIMENSION WEI(2,2)
      DIMENSION AMP(2,2,3)

*** Relaxation constant (in years) for 2002 Denali earthquake
      TAU(1) = 5.0D0

*** Modofied Julian Date (in minutes) for the 2002 Denali earthquake
      IYEAR = 2002
      IMO = 11
      IDAY = 3
      CALL IYMDMJ(IYEAR,IMO,IDAY, MJD)
      ITEQ(1) = MJD*60*24

      DNORTH = 0.0D0
      DEAST = 0.0D0
      DUP = 0.0D0

      DO K = 1, NUMPSG
*** Check if the point is inside the grid
         POSX = YLON*180.d0/PI
         POSX = 360.d0 - POSX
         IF (POSX .GT. 360.D0) POSX = POSX - 360.D0
         POSY = YLAT*180.D0/PI
         CALL GRDCHK(POSX, POSY, PSGLX(K), PSGUX(K),
     1            PSGLY(K), PSGUY(K), INSIDE)

         IF (INSIDE ) THEN
*** Get the indices for the lower left-hand corner of the grid
         CALL PSGWEI(POSX,POSY,K,I,J,WEI)
*** Get the displacement amplitude at the four corners
         CALL GRDAMP(K,I,J,AMP,PS)

         ANORTH = WEI(1,1)*AMP(1,1,1) + WEI(1,2)*AMP(1,2,1)
     1          + WEI(2,1)*AMP(2,1,1) + WEI(2,2)*AMP(2,2,1)     
         AEAST  = WEI(1,1)*AMP(1,1,2) + WEI(1,2)*AMP(1,2,2)
     1          + WEI(2,1)*AMP(2,1,2) + WEI(2,2)*AMP(2,2,2)
         AUP    = WEI(1,1)*AMP(1,1,3) + WEI(1,2)*AMP(1,2,3)
     1          + WEI(2,1)*AMP(2,1,3) + WEI(2,2)*AMP(2,2,3)

c        write (6, 30) anorth, aeast, aup
c  30    format (1x, 3f15.5)

*** Convert amplitudes from mm to meters
         ANORTH = ANORTH / 1000.D0
         AEAST =  AEAST  / 1000.D0
         AUP =    AUP    / 1000.D0

         IF (MIN .GT. ITEQ(K)) THEN
            DTIME = DBLE(MIN - ITEQ(K))/(60.D0*24.D0*365.D0)
            FACTOR = 1.D0 - DEXP(-DTIME/TAU(K))
            DNORTH = DNORTH + ANORTH*FACTOR
            DEAST = DEAST + AEAST*FACTOR
            DUP = DUP + AUP*FACTOR
         ENDIF
         IF (ITREF .GT. ITEQ(K)) THEN
            DTIME = DBLE(ITREF - ITEQ(K))/(60.D0*24.D0*365.D0)
            FACTOR = 1.D0 - DEXP(-DTIME/TAU(K))
            DNORTH = DNORTH - ANORTH*FACTOR
            DEAST = DEAST - AEAST*FACTOR
            DUP = DUP - AUP*FACTOR
         ENDIF
         ENDIF
      ENDDO
      RETURN
      END
