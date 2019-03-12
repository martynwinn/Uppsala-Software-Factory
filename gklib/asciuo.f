C
C
C
      SUBROUTINE asciuo(ITO,IER,TEXT,N,ASCI)
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) ASCI(*)
      CHARACTER*(*) TEXT
C
 1000 FORMAT(A,' ')
C
Code ...
C
      LENTE = MAX (1,LENGTH(TEXT))
      LENMAX = LNBFIO - LENTE - 4
      IER=NTXER1
      IF (N.LE.0) THEN
        IER=NTXER3
      ELSE
        LINE= TEXT(1:LENTE)//' ('
        DO 15 I=1,N
          I1=LENGTH(LINE)+2
          IF ((I1+LEN(ASCI(I))).GT.LENMAX) THEN
            IER=NTXER2
            RETURN
          ENDIF
          WRITE(LINE(I1:),1000,ERR=999) ASCI(I)
   15   CONTINUE
        IER=NTXER0
        LINE = LINE (1:LENGTH(LINE))//') '
        CALL IOLINE(ITO,IER,LINE)
      ENDIF
C
      RETURN
C
  999 CONTINUE
      IER=NTXER4
      RETURN
C
      END
