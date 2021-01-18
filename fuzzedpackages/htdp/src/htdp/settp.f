
**********************************************************
      subroutine SETTP

*** Specify transformation parameters from ITRF94
*** to other reference frames
**************************************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97
C  The latter parameters were added to HTDP in 09/2014. They will be used to transform between
C  ITRF systems. They will not be used if the transformation involves NAD83 or WGS84 (transit).
C  They will be used for the Pacific branches of NAD83.


      implicit double precision (a-h, o-z)
      implicit INTEGER(4) (i-n)
      parameter (numref = 15)

      common /const/ a, f, e2, eps, af, pi, twopi, rhosec
      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)
      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)

C  Parameters computed with the IGS values of ITRF96==>ITRF97

*** From ITRF94 to NAD 83
      tx(1) = 0.9910d0        
      ty(1) = -1.9072d0      
      tz(1) = -.5129d0      
      dtx(1) = 0.d0        
      dty(1) = 0.d0       
      dtz(1) = 0.d0      
      rx(1) = 1.25033d-7
      ry(1) = 0.46785d-7
      rz(1) = 0.56529d-7
      drx(1) = 0.00258d-7
      dry(1) = -.03599d-7
      drz(1) = -.00153d-7
      scale(1) = 0.d0
      dscale(1) = 0.0d0
      refepc(1) = 1997.0d0

*** From ITRF94 to ITRF88
      tx(2) = 0.018d0
      ty(2) = 0.000d0
      tz(2) = -.092d0
      dtx(2) = 0.0d0    
      dty(2) = 0.0d0
      dtz(2) = 0.0d0
      rx(2) = -.0001d0 / rhosec
      ry(2) = 0.0d0                 
      rz(2) = 0.0d0 
      drx(2) = 0.0d0                  
      dry(2) = 0.0d0              
      drz(2) = 0.0d0                   
      scale(2) = 0.74d-8
      dscale(2) = 0.0d0
      refepc(2) = 1988.0d0

*** From ITRF94 to ITRF89
      tx(3) = 0.023d0
      ty(3) = 0.036d0
      tz(3) = -.068d0
      dtx(3) = 0.0d0    
      dty(3) = 0.0d0
      dtz(3) = 0.0d0
      rx(3) = 0.0d0
      ry(3) = 0.0d0                 
      rz(3) = 0.0d0 
      drx(3) = 0.0d0                  
      dry(3) = 0.0d0              
      drz(3) = 0.0d0                   
      scale(3) = 0.43d-8
      dscale(3) = 0.0d0
      refepc(3) = 1988.0d0

*** From ITRF94 to ITRF90
      tx(4) = 0.018d0
      ty(4) = 0.012d0
      tz(4) = -.030d0
      dtx(4) = 0.0d0    
      dty(4) = 0.0d0
      dtz(4) = 0.0d0
      rx(4) = 0.0d0
      ry(4) = 0.0d0                 
      rz(4) = 0.0d0 
      drx(4) = 0.0d0                  
      dry(4) = 0.0d0              
      drz(4) = 0.0d0                   
      scale(4) = 0.09d-8
      dscale(4) = 0.0d0
      refepc(4) = 1988.0d0

*** From ITRF94 to ITRF91
      tx(5) = 0.020d0
      ty(5) = 0.016d0
      tz(5) = -.014d0
      dtx(5) = 0.0d0    
      dty(5) = 0.0d0
      dtz(5) = 0.0d0
      rx(5) = 0.0d0
      ry(5) = 0.0d0                 
      rz(5) = 0.0d0 
      drx(5) = 0.0d0                  
      dry(5) = 0.0d0              
      drz(5) = 0.0d0                   
      scale(5) = 0.06d-8
      dscale(5) = 0.0d0
      refepc(5) = 1988.0d0

*** From ITRF94 to ITRF92
      tx(6) = 0.008d0
      ty(6) = 0.002d0
      tz(6) = -.008d0
      dtx(6) = 0.0d0    
      dty(6) = 0.0d0
      dtz(6) = 0.0d0
      rx(6) = 0.0d0
      ry(6) = 0.0d0                 
      rz(6) = 0.0d0 
      drx(6) = 0.0d0                  
      dry(6) = 0.0d0              
      drz(6) = 0.0d0                   
      scale(6) = -.08d-8
      dscale(6) = 0.0d0
      refepc(6) = 1988.0d0

*** From ITRF94 to ITRF93
      tx(7) = 0.006d0
      ty(7) = -.005d0
      tz(7) = -.015d0
      dtx(7) = -.0029d0
      dty(7) = 0.0004d0
      dtz(7) = 0.0008d0
      rx(7) = 0.00039d0 / rhosec
      ry(7) = -.00080d0 / rhosec
      rz(7) = 0.00096d0 / rhosec
      drx(7) = .00011d0 / rhosec
      dry(7) = .00019d0 / rhosec
      drz(7) =-.00005d0 / rhosec
      scale(7) = 0.04d-8
      dscale(7) = 0.0d0
      refepc(7) = 1988.0d0

*** From ITRF94 to ITRF96
      tx(8) = 0.d0
      ty(8) = 0.d0
      tz(8) = 0.d0
      dtx(8) = 0.d0
      dty(8) = 0.d0
      dtz(8) = 0.d0
      rx(8) = 0.d0
      ry(8) = 0.d0
      rz(8) = 0.d0
      drx(8) = 0.d0
      dry(8) = 0.d0
      drz(8) = 0.d0
      scale(8) = 0.d0
      dscale(8) = 0.0d0
      refepc(8) = 1996.0d0

*** From ITRF94 to ITRF97 (based on IGS adopted values)
*** According to IERS:  ITRF97 = ITRF96 = ITRF94
      tx(9) = 0.00207d0
      ty(9) = 0.00021d0
      tz(9) = -0.00995d0
      dtx(9) = -0.00069d0
      dty(9) = 0.00010d0
      dtz(9) = -0.00186d0
      rx(9) = -0.00012467d0 / rhosec
      ry(9) = 0.00022355d0 / rhosec
      rz(9) = 0.00006065d0 / rhosec
      drx(9) = -0.00001347d0 / rhosec
      dry(9) = 0.00001514d0 / rhosec
      drz(9) = -0.00000027d0 / rhosec
      scale(9) = 0.93496d-9
      dscale(9) = 0.19201d-9
      refepc(9) = 1997.0d0

*** From ITRF94 to WGS 72 (composition of ITRF94 -> NAD_83 -> WGS_72)
      tx(10) = 0.9910d0
      ty(10) = -1.9072d0
      tz(10) = -.5129d0 - 4.5d0
      dtx(10) = 0.d0
      dty(10) = 0.d0
      dtz(10) = 0.d0
      rx(10) = 1.25033d-7
      ry(10) = 0.46785d-7
      rz(10) = 0.56529d-7 + 26.85868d-7
      drx(10) = 0.00258d-7
      dry(10) = -.03599d-7
      drz(10) = -.00153d-7
      scale(10) = 0.d0 - 0.2263d-6
      dscale(10) = 0.0d0
      refepc(10) = 1997.0d0

*** From ITRF94 to ITRF00
*** assumes that ITRF94 = ITRF96 and
*** uses IGS values for ITRF96 -> ITRF97
*** and IERS values for ITRF97 -> ITRF00
      tx(11) = -.00463d0
      ty(11) = -.00589d0
      tz(11) = +.00855d0
      dtx(11) = -0.00069d0
      dty(11) = 0.00070d0
      dtz(11) = -0.00046d0
      rx(11) = -.00012467d0 / rhosec
      ry(11) = 0.00022355d0 / rhosec
      rz(11) = 0.00006065d0 / rhosec
      drx(11) = -0.00001347d0 / rhosec
      dry(11) = 0.00001514d0 / rhosec
      drz(11) = 0.00001973d0 / rhosec
      scale(11) = -0.61504d-9
      dscale(11) = 0.18201d-9
      refepc(11) = 1997.0d0

*** From ITRF94 to PACP00
*** use PA/ITRF00 rotation rates from Beavan et al., (2002)
      tx(12) = 0.9056d0
      ty(12) = -2.0200d0
      tz(12) = -0.5516d0
      dtx(12) = -.00069d0
      dty(12) = 0.00070d0
      dtz(12) = -0.00046d0
      rx(12) = 0.027616d0 / rhosec
      ry(12) = 0.013692d0 / rhosec
      rz(12) = 0.002773d0 / rhosec
      drx(12) = -.000397d0 / rhosec
      dry(12) = 0.001022d0 / rhosec
      drz(12) = -.002166d0 / rhosec
      scale(12) = -0.61504d-9
      dscale(12) = 0.18201d-9
      refepc(12) = 1997.0d0

*** From ITRF94 to MARP00
*** Use velocity of GUAM
      tx(13) = 0.9056d0
      ty(13) = -2.0200d0
      tz(13) = -0.5516d0
      dtx(13) = -0.00069d0
      dty(13) = 0.00070d0
      dtz(13) = -0.00046d0
      rx(13) = .028847d0 / rhosec
      ry(13) = .010644d0 / rhosec
      rz(13) = 0.008989d0 / rhosec
      drx(13) = -0.000033d0 / rhosec
      dry(13) =  0.000120d0 / rhosec
      drz(13) = -0.000327d0 / rhosec
      scale(13) = -0.61504d-9
      dscale(13) = 0.18201d-9
      refepc(13) = 1997.00d0

*** From ITRF94 to ITRF2005
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
      tx(14) = -0.00533d0
      ty(14) = -0.00479d0
      tz(14) =  0.00895d0
      dtx(14) = -0.00049d0
      dty(14) =  0.00060d0
      dtz(14) =  0.00134d0
      rx(14) = -.00012467d0 / rhosec
      ry(14) = 0.00022355d0 / rhosec
      rz(14) = 0.00006065d0 / rhosec
      drx(14) = -0.00001347d0 / rhosec
      dry(14) = 0.00001514d0 / rhosec
      drz(14) = 0.00001973d0 / rhosec
      scale(14) = -0.77504d-9
      dscale(14) = 0.10201d-9
      refepc(14) = 1997.0d0

*** From ITRF94 to ITRF2008 (also IGS08 and IGB08)
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008 (and IGS08 and IGB08)
      tx(15) = -0.00243d0
      ty(15) = -0.00389d0
      tz(15) =  0.01365d0
      dtx(15) = -0.00079d0
      dty(15) =  0.00060d0
      dtz(15) =  0.00134d0
      rx(15) = -0.00012467d0 / rhosec
      ry(15) =  0.00022355d0 / rhosec
      rz(15) =  0.00006065d0 / rhosec
      drx(15) = -0.00001347d0 / rhosec
      dry(15) =  0.00001514d0 / rhosec
      drz(15) =  0.00001973d0 / rhosec
      scale(15) = -1.71504d-9
      dscale(15) = 0.10201d-9
      refepc(15) = 1997.0d0

C*************************************************************************************************************************
C  Parameters computed with the IERS values of ITRF96==>ITRF97

*** From ITRF94 to NAD 83
      tx1(1) = 0.9910d0        
      ty1(1) = -1.9072d0      
      tz1(1) = -.5129d0      
      dtx1(1) = 0.d0        
      dty1(1) = 0.d0       
      dtz1(1) = 0.d0      
      rx1(1) = 1.25033d-7
      ry1(1) = 0.46785d-7
      rz1(1) = 0.56529d-7
      drx1(1) = 0.00258d-7
      dry1(1) = -.03599d-7
      drz(1) = -.00153d-7
      scale(1) = 0.d0
      dscale1(1) = 0.0d0
      refepc1(1) = 1997.0d0

*** From ITRF94 to ITRF88
      tx1(2) = 0.018d0
      ty1(2) = 0.000d0
      tz1(2) = -.092d0
      dtx1(2) = 0.0d0    
      dty1(2) = 0.0d0
      dtz1(2) = 0.0d0
      rx1(2) = -.0001d0 / rhosec
      ry1(2) = 0.0d0                 
      rz1(2) = 0.0d0 
      drx1(2) = 0.0d0                  
      dry1(2) = 0.0d0              
      drz1(2) = 0.0d0                   
      scale1(2) = 0.74d-8
      dscale1(2) = 0.0d0
      refepc1(2) = 1988.0d0

*** From ITRF94 to ITRF89
      tx1(3) = 0.023d0
      ty1(3) = 0.036d0
      tz1(3) = -.068d0
      dtx1(3) = 0.0d0    
      dty1(3) = 0.0d0
      dtz1(3) = 0.0d0
      rx1(3) = 0.0d0
      ry1(3) = 0.0d0                 
      rz1(3) = 0.0d0 
      drx1(3) = 0.0d0                  
      dry1(3) = 0.0d0              
      drz1(3) = 0.0d0                   
      scale1(3) = 0.43d-8
      dscale1(3) = 0.0d0
      refepc1(3) = 1988.0d0

*** From ITRF94 to ITRF90
      tx1(4) = 0.018d0
      ty1(4) = 0.012d0
      tz1(4) = -.030d0
      dtx1(4) = 0.0d0    
      dty1(4) = 0.0d0
      dtz1(4) = 0.0d0
      rx1(4) = 0.0d0
      ry1(4) = 0.0d0                 
      rz1(4) = 0.0d0 
      drx1(4) = 0.0d0                  
      dry1(4) = 0.0d0              
      drz1(4) = 0.0d0                   
      scale1(4) = 0.09d-8
      dscale1(4) = 0.0d0
      refepc1(4) = 1988.0d0

*** From ITRF94 to ITRF91
      tx1(5) = 0.020d0
      ty1(5) = 0.016d0
      tz1(5) = -.014d0
      dtx1(5) = 0.0d0    
      dty1(5) = 0.0d0
      dtz1(5) = 0.0d0
      rx1(5) = 0.0d0
      ry1(5) = 0.0d0                 
      rz1(5) = 0.0d0 
      drx1(5) = 0.0d0                  
      dry1(5) = 0.0d0              
      drz1(5) = 0.0d0                   
      scale1(5) = 0.06d-8
      dscale1(5) = 0.0d0
      refepc1(5) = 1988.0d0

*** From ITRF94 to ITRF92
      tx1(6) = 0.008d0
      ty1(6) = 0.002d0
      tz1(6) = -.008d0
      dtx1(6) = 0.0d0    
      dty1(6) = 0.0d0
      dtz1(6) = 0.0d0
      rx1(6) = 0.0d0
      ry1(6) = 0.0d0                 
      rz1(6) = 0.0d0 
      drx1(6) = 0.0d0                  
      dry1(6) = 0.0d0              
      drz1(6) = 0.0d0                   
      scale1(6) = -.08d-8
      dscale1(6) = 0.0d0
      refepc1(6) = 1988.0d0

*** From ITRF94 to ITRF93
      tx1(7) = 0.006d0
      ty1(7) = -.005d0
      tz1(7) = -.015d0
      dtx1(7) = -.0029d0
      dty1(7) = 0.0004d0
      dtz1(7) = 0.0008d0
      rx1(7) = 0.00039d0 / rhosec
      ry1(7) = -.00080d0 / rhosec
      rz1(7) = 0.00096d0 / rhosec
      drx1(7) = .00011d0 / rhosec
      dry1(7) = .00019d0 / rhosec
      drz1(7) =-.00005d0 / rhosec
      scale1(7) = 0.04d-8
      dscale1(7) = 0.0d0
      refepc1(7) = 1988.0d0

*** From ITRF94 to ITRF96
      tx1(8) = 0.d0
      ty1(8) = 0.d0
      tz1(8) = 0.d0
      dtx1(8) = 0.d0
      dty1(8) = 0.d0
      dtz1(8) = 0.d0
      rx1(8) = 0.d0
      ry1(8) = 0.d0
      rz1(8) = 0.d0
      drx1(8) = 0.d0
      dry1(8) = 0.d0
      drz1(8) = 0.d0
      scale1(8) = 0.d0
      dscale1(8) = 0.0d0
      refepc1(8) = 1996.0d0

*** From ITRF94 to ITRF97 (based on IERS adopted values)
*** According to IERS:  ITRF97 = ITRF96 = ITRF94
      tx1(9)     = 0.00000d0
      ty1(9)     = 0.00000d0
      tz1(9)     = 0.00000d0
      dtx1(9)    = 0.00000d0
      dty1(9)    = 0.00000d0
      dtz1(9)    = 0.00000d0
      rx1(9)     = 0.00000000d0 
      ry1(9)     = 0.00000000d0
      rz1(9)     = 0.00000000d0 
      drx1(9)    = 0.00000000d0 
      dry1(9)    = 0.00000000d0 
      drz1(9)    = 0.00000000d0 
      scale1(9)  = 0.d0          
      dscale1(9) = 0.d0          
      refepc1(9) = 2000.0d0

*** From ITRF94 to WGS 72 (composition of ITRF94 -> NAD_83 -> WGS_72)
      tx1(10) = 0.9910d0
      ty1(10) = -1.9072d0
      tz1(10) = -.5129d0 - 4.5d0
      dtx1(10) = 0.d0
      dty1(10) = 0.d0
      dtz1(10) = 0.d0
      rx1(10) = 1.25033d-7
      ry1(10) = 0.46785d-7
      rz1(10) = 0.56529d-7 + 26.85868d-7
      drx1(10) = 0.00258d-7
      dry1(10) = -.03599d-7
      drz1(10) = -.00153d-7
      scale1(10) = 0.d0 - 0.2263d-6
      dscale1(10) = 0.0d0
      refepc1(10) = 1997.0d0

*** From ITRF94 to ITRF00
*** assumes that         ITRF94 = ITRF96 and
*** uses IERS values for ITRF96 -> ITRF97
*** and  IERS values for ITRF97 -> ITRF00
      tx1(11)     = -.00670d0
      ty1(11)     = -.00430d0
      tz1(11)     = +.02270d0
      dtx1(11)    = 0.00000d0
      dty1(11)    = 0.00060d0
      dtz1(11)    = 0.00140d0
      rx1(11)     = 0.00000000d0 / rhosec
      ry1(11)     = 0.00000000d0 / rhosec
      rz1(11)     = +0.00006000d0 / rhosec
      drx1(11)    = 0.00000000d0 / rhosec
      dry1(11)    = 0.00000000d0 / rhosec
      drz1(11)    = +0.00002000d0 / rhosec
      scale1(11)  = -1.58000d-9
      dscale1(11) = -0.01000d-9
      refepc1(11) = 2000.0d0

*** From ITRF94 to PACP00
*** use PA/ITRF00 rotation rates from Beavan et al., (2002)
      tx1(12)     = 0.9035d0
      ty1(12)     = -2.0202d0
      tz1(12)     = -0.5417d0
      dtx1(12)    = 0.00000d0
      dty1(12)    = 0.00060d0
      dtz1(12)    = 0.0014d0
      rx1(12)     = 0.027741d0 / rhosec
      ry1(12)     = 0.013469d0 / rhosec
      rz1(12)     = 0.002712d0 / rhosec
      drx1(12)    = -.000384d0 / rhosec
      dry1(12)    = 0.001007d0 / rhosec
      drz1(12)    = -.002166d0 / rhosec
      scale1(12)  = -1.55000d-9
      dscale1(12) = -0.010000d-9
      refepc1(12) = 1997.0d0

*** From ITRF94 to MARP00
*** Use velocity of GUAM
      tx1(13) = 0.9035d0
      ty1(13) = -2.0202d0
      tz1(13) = -0.5417d0
      dtx1(13) = -0.00000d0
      dty1(13) = 0.00060d0
      dtz1(13) = 0.00140d0
      rx1(13) = .028971d0 / rhosec
      ry1(13) = .01042d0 / rhosec
      rz1(13) = 0.008928d0 / rhosec
      drx1(13) = -0.00002d0 / rhosec
      dry1(13) =  0.000105d0 / rhosec
      drz1(13) = -0.000327d0 / rhosec
      scale1(13) = -1.55000d-9
      dscale1(13) = -0.01000d-9
      refepc1(13) = 1997.00d0

*** From ITRF94 to ITRF2005
*** assumes that         ITRF94   = ITRF96
*** uses IERS values for ITRF96   -> ITRF97
*** uses IERS values for ITRF97   -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
      tx1(14)     = -0.00680d0
      ty1(14)     = -0.00350d0
      tz1(14)     =  0.0285d0
      dtx1(14)    =  0.00020d0
      dty1(14)    =  0.00050d0
      dtz1(14)    =  0.00320d0
      rx1(14)     = 0.00000000d0 / rhosec
      ry1(14)     = 0.00000000d0 / rhosec
      rz1(14)     = +0.00006000d0 / rhosec
      drx1(14)    = 0.00000000d0 / rhosec
      dry1(14)    = 0.00000000d0 / rhosec
      drz1(14)    = +0.00002000d0 / rhosec
      scale1(14)  = -1.98000d-9
      dscale1(14) = -0.09000d-9
      refepc1(14) = 2000.0d0

*** From ITRF94 to ITRF2008 (also IGS08 and IGB08)
*** assumes that ITRF94 = ITRF96
*** uses IERS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008 (and IGS08 and IGB08)
      tx1(15)     = -0.00480d0
      ty1(15)     = -0.00260d0
      tz1(15)     =  0.03320d0
      dtx1(15)    = -0.00010d0
      dty1(15)    =  0.00050d0
      dtz1(15)    =  0.00320d0
      rx1(15)     =  0.00000000d0 / rhosec
      ry1(15)     =  0.00000000d0 / rhosec
      rz1(15)     = +0.00006000d0 / rhosec
      drx1(15)    = -0.00000000d0 / rhosec
      dry1(15)    =  0.00000000d0 / rhosec
      drz1(15)    = +0.00002000d0 / rhosec
      scale1(15)  = -2.92d-9
      dscale1(15) = -0.09d-9
      refepc1(15) = 2000.0d0

      return
      end
