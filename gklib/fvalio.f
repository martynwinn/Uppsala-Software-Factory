C
C
      SUBROUTINE fvalio(ITI,ITO,IER,TEXT,N,RVALUE)
C
      INCLUDE 'ioin.incl'
C
      real rvalue(*)
      CHARACTER*(*) TEXT
C
      CALL FVALUO(ITO,IER,TEXT,N,RVALUE)
      IF (IER.EQ.0) CALL FVALII(ITI,IER,N,RVALUE)
C
      RETURN
      END
