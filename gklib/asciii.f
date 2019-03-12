C
C
C
      SUBROUTINE asciii(ITI,IER,N,ASCI)
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) ASCI(*)
      CHARACTER*(LNBFIO+1) B
C
 1010 FORMAT(A)
C
Code ...
C
      IER=NTXER1
      LINE=' '
      READ(ITI,1010,END=20,ERR=30,IOSTAT=IER) LINE
      IF (LENGTH(LINE).GT.LNBFIO) GOTO 40
   20 K=LENGTH(LINE)
      IF (K.GT.0) THEN
        J=1
        K=LEN(LINE)
        DO 35 I=1,N
          IF (J.GT.K) GOTO 45
          READ(LINE(J:K),1010,END=45,ERR=30,IOSTAT=IER) B
          IF(LENGTH(B).GT.0 .AND.
CCC     +       B(1:1).NE.'/' .AND.
     +       B(1:1).NE.',') THEN
               ASCI(I)=B(1:LEN(ASCI(I)))
          ENDIF
          J=1+J+LEN(ASCI(I))
   35   CONTINUE
      ENDIF
   45 IER=NTXER0
   30 RETURN
C
   40 IER=NTXER2
C
      RETURN
      END
