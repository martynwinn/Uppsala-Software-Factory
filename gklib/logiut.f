C
C
C
      SUBROUTINE logiut(TEXT,N,LVAL)
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
      dollar = .false.
      CALL LOGIUO(ITXOUT,IER,TEXT,N,LVAL)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' logiuo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
