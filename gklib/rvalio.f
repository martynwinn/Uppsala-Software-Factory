C
C
C
      SUBROUTINE rvalio(ITI,ITO,IER,TEXT,N,RVALUE)
C
      INCLUDE 'ioin.incl'
C
      real rvalue(*)
      CHARACTER*(*) TEXT
C
      CALL RVALUO(ITO,IER,TEXT,N,RVALUE)
      IF (IER.EQ.0) CALL RVALII(ITI,IER,N,RVALUE)
C
      RETURN
      END
