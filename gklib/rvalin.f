C
C
C
      SUBROUTINE rvalin(TEXT,N,RVALUE)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      rval   .... array of real        variables
C
      INCLUDE 'ioin.incl'
C
      REAL RVALUE(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = xinter()
      CALL RVALUO(ITXOUT,IER,TEXT,N,RVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' rvaluo -- i/o error',IER)
        RETURN
      ENDIF
C
      CALL RVALII(ITXINP,IER,N,RVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' rvalii -- i/o error',IER)
      ENDIF
C
      RETURN
      END
