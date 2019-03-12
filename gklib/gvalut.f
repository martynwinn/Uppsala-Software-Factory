C
C
C
      SUBROUTINE gvalut(TEXT,N,RVALUE)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      rval   .... array of real        variables (F or E format)
C
      INCLUDE 'ioin.incl'
C
      REAL RVALUE(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = .false.
      CALL GVALUO(ITXOUT,IER,TEXT,N,RVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' gvaluo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
