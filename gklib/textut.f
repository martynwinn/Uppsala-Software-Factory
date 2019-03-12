C
C
C
      SUBROUTINE textut(TEXT,STRING)
C
C      text   .... output   character*(*) string
C      string ....          character*(*) variable
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) STRING
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = .false.
      CALL TEXTUO(ITXOUT,IER,TEXT,STRING)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' textuo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
