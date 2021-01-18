*********************************************************
      SUBROUTINE SETRF

*** Specify arrays that may be used to convert
*** a reference frame identifier in the blue book
*** to a reference frame identifier in HTDP and back

      IMPLICIT INTEGER(4) (I-N)
      parameter ( numref = 15 )
      COMMON /REFCON/ IRFCON(28), JRFCON(numref)

*** From blue book identifier to HTDP indentifier
*** WGS 72 Precise
      IRFCON(1) = 10

*** WGS 84 (orig) Precise (set  equal to NAD 83)
      IRFCON(2) = 1

*** WGS 72 Broadcast
      IRFCON(3) = 10

*** WGS 84 (orig) Broadcast (set equal to NAD 83)
      IRFCON(4) = 1

*** ITRF89
      IRFCON(5) = 3

*** PNEOS 90 or NEOS 91.25 (set equal to ITRF90)
      IRFCON(6) = 4

*** NEOS 90 (set equal to ITRF90)
      IRFCON(7) = 4

*** ITRF91
      IRFCON(8) = 5

*** SIO/MIT 92.57 (set equal to ITRF91)
      IRFCON(9) = 5

*** ITRF91
      IRFCON(10) = 5

*** ITRF92
      IRFCON(11) = 6

*** ITRF93
      IRFCON(12) = 7

*** WGS 84 (G730) Precise (set equal to ITRF91)
      IRFCON(13) = 5

*** WGS 84 (G730) Broadcast (set equal to ITRF91)
      IRFCON(14) = 5

*** ITRF94
      IRFCON(15) = 8

*** WGS 84 (G873) Precise  (set equal to ITRF94)
      IRFCON(16) = 8

*** WGS 84 (G873) Broadcast (set equal to ITRF94)
      IRFCON(17) = 8

*** ITRF96
      IRFCON(18) = 8

*** ITRF97
      IRFCON(19) = 9

*** IGS97
      IRFCON(20) = 9

*** ITRF00
      IRFCON(21) = 11

*** IGS00
      IRFCON(22) = 11

*** WGS 84 (G1150)
      IRFCON(23) = 11

*** IGb00
      IRFCON(24) = 11

*** ITRF2005
      IRFCON(25) = 14

*** IGS05
      IRFCON(26) = 14

*** ITRF2008 or IGS08
      IRFCON(27) = 15

*** IGB08
      IRFCON(28) = 15

*** From HTDP identifier to blue book identifier
*** NAD 83 (set equal to WGS 84 (transit))
      JRFCON(1) = 2

*** ITRF88 (set equal to ITRF89)
      JRFCON(2) = 5

*** ITRF89
      JRFCON(3) = 5

*** ITRF90 (set equal to NEOS 90)
      JRFCON(4) = 7

*** ITRF91
      JRFCON(5) = 8

*** ITRF92
      JRFCON(6) = 11

*** ITRF93
      JRFCON(7) = 12

*** ITRF96 (= ITRF94)
      JRFCON(8) = 18

*** ITRF97
      JRFCON(9) = 19

*** WGS 72
      JRFCON(10) = 1

*** ITRF00
      JRFCON(11) = 21

*** NAD 83 (PACP00)
      JRFCON(12) = 2

*** NAD 83 (MARP00)
      JRFCON(13) = 2

*** ITRF2005 or IGS05
      JRFCON(14) = 26

*** ITRF2008 or IGS08
      JRFCON(15) = 27

*** IGB08
      JRFCON(15) = 28

      RETURN
      END
