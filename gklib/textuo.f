C
C
C
      SUBROUTINE textuo(ITO,IER,TEXT,STRING)
C
      INCLUDE 'ioin.incl'
C
C     ASCII computer
C     / or , or   leave TEXT unchanged
C
      integer lenesv
c
      CHARACTER*(*) STRING
      CHARACTER*(*) TEXT
C
 1000 FORMAT(A)
C
Code ...
C
      LENTE = MAX (1,LENGTH(TEXT))
      lenesv = max(1,length(string))
      LENMAX = LNBFIO - LENTE - 4
      IER=NTXER1
      IF (LENGTH(STRING).LE.LENMAX) THEN
        LINE= TEXT(1:LENTE)//' ('//STRING(1:lenesv)//') '
        IER=NTXER0
        CALL IOLINE(ITO,IER,LINE)
      ELSE
        IER=NTXER2
      ENDIF
C
      RETURN
C
  999 CONTINUE
      IER=NTXER4
      RETURN
C
      END
