
************************************************************************
      SUBROUTINE COMPSN(YLATT,YLONT,HTT,YLAT,YLON,HT,
     1                  MIN,VN, VE, VU)

*** Compute the position of a point at specified time
*** Upon input VN, VE, and VU are in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(4) (I-N)
      parameter (NDLOC = 2195)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /QPARM/ STRIKE(NDLOC), HL(NDLOC), EQLAT(NDLOC),
     1          EQLON(NDLOC), SSLIP(NDLOC), DSLIP(NDLOC),
     1          DIP(NDLOC), DEPTH(NDLOC), WIDTH(NDLOC),
     1               EQLATR(50),EQLONR(50),EQRAD(50),
     1               ITEQK(50),NLOC(50),NFP(50),NUMEQ
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

** Compute the contribution due to constant velocity
         DTIME = DBLE(MIN - ITREF) / 525960.D0
         CALL RADR8T(YLAT,VN,VE,VNR,VER)
         YLATT = YLAT + VNR*DTIME
         YLONT = YLON - VER*DTIME
         HTT   = HT + ((VU * DTIME) /1000.D0)
       
** Compute the contribution due to earthquakes.
** It is assumed that the components of displacement,
** DNORTH,DWEST,DUP, do not vary from one reference
** frame to another given the accuracy of dislocation
** models.
      DO 10 I = 1, NUMEQ
          IF(ITEQK(I) .GT. ITREF) THEN 
               NTIME = 1
          ELSE
               NTIME = 0
          ENDIF
          IF(MIN .LT. ITEQK(I)) NTIME = NTIME - 1
          IF(NTIME .NE. 0) THEN 
             CALL RADII(EQLATR(I),RADMER,RADPAR)
             DDLAT = (YLAT - EQLATR(I))*RADMER
             DDLON = (YLON - EQLONR(I))*RADPAR
             DIST = DSQRT(DDLAT*DDLAT + DDLON*DDLON)
             IF(DIST .LE. EQRAD(I)) THEN
               ISTART = NLOC(I)
               IEND = NLOC(I) + NFP(I) - 1
               DO 5 JREC = ISTART,IEND
                 CALL DISLOC(YLAT,YLON,STRIKE(JREC),HL(JREC),
     &              EQLAT(JREC),EQLON(JREC),SSLIP(JREC),     
     &              DSLIP(JREC),DIP(JREC),DEPTH(JREC),
     &              WIDTH(JREC),DNORTH,DWEST,DUP)     
                 YLATT = YLATT + NTIME*DNORTH
                 YLONT = YLONT + NTIME*DWEST     
                 HTT   = HTT   + NTIME*DUP
    5          CONTINUE
             ENDIF
          ENDIF  
   10 CONTINUE

*** Compute contribution due to postseismic deformation
      CALL PSDISP(YLAT, YLON, MIN, DNORTH, DEAST, DUP)
      CALL RADII (YLAT, RMER, RPAR)
      YLATT = YLATT + DNORTH/RMER
      YLONT = YLONT - DEAST/RPAR
      HTT  = HTT + DUP
      RETURN
      END
