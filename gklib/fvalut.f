C
C
C
      SUBROUTINE fvalut(TEXT,N,RVALUE)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      rval   .... array of real        variables (must fit in F8.3)
C
      INCLUDE 'ioin.incl'
C
      REAL RVALUE(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = .false.
      CALL FVALUO(ITXOUT,IER,TEXT,N,RVALUE)
C
C ... in case of format overflow, resort to GVALUO
C
      IF (IER.EQ.NTXER4) THEN
        CALL GVALUO(ITXOUT,IER,TEXT,N,RVALUE)
      ENDIF
C
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' fvaluo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
