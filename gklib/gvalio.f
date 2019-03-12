C
C
      SUBROUTINE gvalio(ITI,ITO,IER,TEXT,N,RVALUE)
C
      INCLUDE 'ioin.incl'
C
      real rvalue(*)
      CHARACTER*(*) TEXT
C
      CALL GVALUO(ITO,IER,TEXT,N,RVALUE)
      IF (IER.EQ.0) CALL GVALII(ITI,IER,N,RVALUE)
C
      RETURN
      END
