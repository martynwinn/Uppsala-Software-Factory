C
C
C
      SUBROUTINE jvalin(TEXT,N,JVALUE)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      jval   .... array of integer     variables
C
      INCLUDE 'ioin.incl'
C
      INTEGER JVALUE(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = xinter()
      CALL JVALUO(ITXOUT,IER,TEXT,N,JVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' jvaluo -- i/o error',IER)
        RETURN
      ENDIF
C
      CALL JVALII(ITXINP,IER,N,JVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' jvalii -- i/o error',IER)
      ENDIF
C
      RETURN
      END
