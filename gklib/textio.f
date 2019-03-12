C
C
C
      SUBROUTINE textio(ITI,ITO,IER,TEXT,STRING)
C
C     ASCII computer
C     / or , or   leave TEXT unchanged
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) STRING
      CHARACTER*(*) TEXT
C
      CALL TEXTUO(ITO,IER,TEXT,STRING)
      IF (IER.EQ.0) CALL TEXTII(ITI,IER,STRING)
C
      RETURN
      END
