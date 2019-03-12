C
C..........................................................................
C
      SUBROUTINE ivalio(ITI,ITO,IER,TEXT,N,IVALUE)
C
      INCLUDE 'ioin.incl'
C
      integer ivalue(*)
      CHARACTER*(*) TEXT
C
      CALL IVALUO(ITO,IER,TEXT,N,IVALUE)
      IF (IER.EQ.0) CALL IVALII(ITI,IER,N,IVALUE)
C
      RETURN
      END
