C
C..........................................................................
C
      SUBROUTINE ivalut(TEXT,N,IVALUE)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      ival   .... array of integer     variables
C
      INCLUDE 'ioin.incl'
C
      INTEGER IVALUE(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = .false.
      CALL IVALUO(ITXOUT,IER,TEXT,N,IVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' ivaluo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
