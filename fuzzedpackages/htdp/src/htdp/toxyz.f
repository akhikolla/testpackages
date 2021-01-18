***************************************************************
C     SUBROUTINE GETTIM(MONTH, IDAY, IYEAR, DATE, MINS, TEST)
C
*** Read month-day-year and convert to decimal years
*** and Julian time in minutes      
***    MONTH      input - number from 1 to 12
***    IDAY       input - number from 1 to 31
***    IYEAR      input - must be after 1906
***    DATE       output - corresponding time in decimal years
***    MINS       output - corresponding julian time in minutes
***    TEST       output - if (true) then there is an error
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     IMPLICIT INTEGER(4) (I-N)
C     LOGICAL TEST
C     COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

C     READ(LUIN,*) MONTH,IDAY,IYEAR
    
C     IF(IYEAR .le. 1906) THEN
C         WRITE(LUOUT,10)
C  10     FORMAT(' The model is not valid for dates prior ',
C    1           'to 1906.'/)
C         TEST = .TRUE.
C         RETURN
C     ENDIF

C     IF(MONTH .le. 0 .or. MONTH .gt. 12) THEN
C         WRITE(LUOUT,20)
C  20     FORMAT(' Improper month specified.'/)
C         TEST = .TRUE.
C         RETURN
C     ENDIF

C     IF(IDAY .le. 0 .or. IDAY .gt. 31) THEN
C         WRITE(LUOUT,30)
C  30     FORMAT(' Improper day specified.'/)
C         TEST = .TRUE.
C         RETURN
C     ENDIF

C     CALL TOTIME(IYEAR, MONTH, IDAY, MINS)
C     CALL TOTIME(IYEAR, 1, 1, MIN00)
C     DATE = DBLE(IYEAR) + DBLE(MINS - MIN00)/525600.D0
C     TEST = .FALSE.
C     RETURN
C     END

C************************************************************************
      SUBROUTINE TOXYZ(glat,glon,eht,x,y,z)
 
*** compute x,y,z
*** ref p.17 geometric geodesy notes vol 1, osu, rapp
 
      implicit double precision(a-h,o-z)
      common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
 
      slat=dsin(glat)
      clat=dcos(glat)
      w=dsqrt(1.d0-e2*slat*slat)
      en=a/w
 
      x=(en+eht)*clat*dcos(glon)
      y=(en+eht)*clat*dsin(glon)
      z=(en*(1.d0-e2)+eht)*slat
 
      return
      end
