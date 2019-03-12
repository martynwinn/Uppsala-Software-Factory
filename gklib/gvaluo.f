C
C
C
      SUBROUTINE gvaluo(ITO,IER,TEXT,N,RVALUE)
C
      INCLUDE 'ioin.incl'
C
      PARAMETER (LENFMT=12)
      PARAMETER (BIGGY=1.0E6, SMALLY=0.1)
C
      REAL RVALUE(*)
      CHARACTER*(*) TEXT
C
 1000 FORMAT(1P,E11.3,' ')
 1001 FORMAT(0P,F11.3,' ')
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
            RABS = ABS (RVALUE(I))
            IF (RABS.GT.BIGGY.OR.RABS.LT.SMALLY) THEN
              WRITE(LINE(I1:),1000,ERR=999) RVALUE(I)
            ELSE
              WRITE(LINE(I1:),1001,ERR=999) RVALUE(I)
            ENDIF
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
