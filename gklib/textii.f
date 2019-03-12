C
C
C
      SUBROUTINE textii(ITI,IER,STRING)
C
C     ASCII computer
C     , or   leave TEXT unchanged
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) STRING
C
 1000 FORMAT(A)
C
Code ...
C
      IER=NTXER1
      LINE=' '
      READ(ITI,1000,END=20,ERR=30,IOSTAT=IER) LINE
      IF (LENGTH(LINE).GT.LNBFIO) GOTO 40
   20 IF (LENGTH(LINE).GT.0 .AND.
     +    LINE(1:1).NE.',') THEN
            STRING = LINE
      ENDIF
      IER=NTXER0
   30 RETURN
C
   40 IER=NTXER2
C
      RETURN
      END
