C
C
C
      SUBROUTINE rvalut(TEXT,N,RVALUE)
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
      dollar = .false.
      CALL RVALUO(ITXOUT,IER,TEXT,N,RVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' rvaluo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
