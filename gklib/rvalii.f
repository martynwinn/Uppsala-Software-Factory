C
C
C
      SUBROUTINE rvalii(ITI,IER,N,RVALUE)
C
      INCLUDE 'ioin.incl'
C
      REAL RVALUE(*)
C
 1010 FORMAT(A)
C
Code ...
C
      IER=NTXER1
      LINE=' '
      READ (ITI,1010,END=20,ERR=30,IOSTAT=IER) LINE
      IF (LENGTH(LINE).GT.LNBFIO) GOTO 40
   20 IF (LENGTH(LINE).GT.0)
     +  READ(LINE,*,END=25,ERR=30,IOSTAT=IER)(RVALUE(I),I=1,N)
   25 IER=NTXER0
   30 RETURN
C
   40 IER=NTXER2
C
      RETURN
      END
