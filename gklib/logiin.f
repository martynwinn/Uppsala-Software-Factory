C
C
C
      SUBROUTINE logiin(TEXT,N,LVAL)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      lval   .... array of logical       variables
C
      INCLUDE 'ioin.incl'
C
      LOGICAL LVAL(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = xinter()
      CALL LOGIUO(ITXOUT,IER,TEXT,N,LVAL)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' logiuo -- i/o error',IER)
        RETURN
      ENDIF
C
      CALL LOGIII(ITXINP,IER,N,LVAL)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' logiii -- i/o error',IER)
      ENDIF
C
      RETURN
      END
