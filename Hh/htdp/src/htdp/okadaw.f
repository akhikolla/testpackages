*************************************************************
      SUBROUTINE OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1           VERT,U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      LOGICAL VERT   

      YBAR = ETA*CDIP + Q*SDIP
      DBAR = ETA*SDIP - Q*CDIP
      R = DSQRT(PSI*PSI + ETA*ETA + Q*Q)
      X = DSQRT(PSI*PSI + Q*Q)
      IF(DABS(Q) .LE. 0.1d0) THEN
         TERM = 0.D0
      ELSE
         TERM = DATAN(PSI*ETA/(Q*R))
      ENDIF

      IF(VERT) THEN
         F5 = -RATIO*PSI*SDIP/(R + DBAR)
         F4 = -RATIO*Q/(R + DBAR)
         F3 = 0.5D0*RATIO*(ETA/(R + DBAR)
     1          + YBAR*Q/((R + DBAR)*(R + DBAR))
     2          - DLOG(R + ETA))
         F1 = -0.5D0*RATIO*PSI*Q/
     1        ((R + DBAR)*(R + DBAR))
      ELSE
         IF(DABS(PSI) .LE. 0.1D0) then
            F5 = 0.d0
         ELSE
            F5 = 2.D0*RATIO*
     1      DATAN((ETA*(X+Q*CDIP)+X*(R+X)*SDIP)/(PSI*(R+X)*CDIP))
     2          /CDIP
         ENDIF
         F4 = RATIO*(DLOG(R+DBAR)-SDIP*DLOG(R+ETA))/CDIP
         F3 = RATIO*(YBAR/(CDIP*(R+DBAR)) - DLOG(R+ETA))
     1         + SDIP*F4/CDIP
         F1 = -RATIO*(PSI/(CDIP*(R+DBAR))) - SDIP*F5/CDIP
      ENDIF
         F2 = -RATIO*DLOG(R+ETA) - F3

      U1SS = -(PSI*Q/(R*(R+ETA))
     1         + TERM + F1*SDIP)/TWOPI
      U2SS = -(YBAR*Q/(R*(R+ETA))
     1         + Q*CDIP/(R+ETA)
     2         + F2*SDIP)/TWOPI
      U3SS = -(DBAR*Q/(R*(R+ETA))
     1         + Q*SDIP/(R+ETA)
     2         + F4*SDIP)/TWOPI
      U1DS = -(Q/R - F3*SDIP*CDIP)/TWOPI
      U2DS = -(YBAR*Q/(R*(R+PSI))
     1         + CDIP*TERM - F1*SDIP*CDIP)/TWOPI
      U3DS = -(DBAR*Q/(R*(R+PSI))
     1         + SDIP*TERM - F5*SDIP*CDIP)/TWOPI
      RETURN
      END
