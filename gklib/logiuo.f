C
C
C
      SUBROUTINE logiuo(ITO,IER,TEXT,N,LOGI)
C
      INCLUDE 'ioin.incl'
C
      PARAMETER(LENFMT=2)
C
      LOGICAL LOGI(*)
      CHARACTER*(*) TEXT
C
 1000 FORMAT(L1,' ')
C
Code ...
C
      LENTE = MAX (1,LENGTH(TEXT))
      LENMAX = LNBFIO - LENTE - 4
      IER=NTXER1
      IF (N.LE.0) THEN
        IER=NTXER3
      ELSE
        IF (LENFMT*N.LE.LENMAX) THEN
          LINE= TEXT(1:LENTE)//' ('
          I1=LENTE+3
          DO 15 I=1,N
            WRITE(LINE(I1:),1000,ERR=999) LOGI(I)
            I1=I1+LENFMT
   15     CONTINUE
          IER=NTXER0
          LINE = LINE (1:LENGTH(LINE))//') '
          CALL IOLINE(ITO,IER,LINE)
        ELSE
          IER=NTXER2
        ENDIF
      ENDIF
C
      RETURN
C
  999 CONTINUE
      IER=NTXER4
      RETURN
C
      END
