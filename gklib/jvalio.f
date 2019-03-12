C
C
C
      SUBROUTINE jvalio(ITI,ITO,IER,TEXT,N,JVALUE)
C
      INCLUDE 'ioin.incl'
C
      integer jvalue(*)
      CHARACTER*(*) TEXT
C
      CALL JVALUO(ITO,IER,TEXT,N,JVALUE)
      IF (IER.EQ.0) CALL JVALII(ITI,IER,N,JVALUE)
C
      RETURN
      END
