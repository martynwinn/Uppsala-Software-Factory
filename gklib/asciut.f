C
C
C
      SUBROUTINE asciut(TEXT,N,ACHR)
C
C      text   .... output   character*(*) string
C      n      .... # input  variables
C      achr   .... array of character*(*) variables
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) ACHR(*)
      CHARACTER*(*) TEXT
C
Code ...
C
      dollar = .false.
      CALL ASCIUO(ITXOUT,IER,TEXT,N,ACHR)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' asciuo -- i/o error',IER)
      ENDIF
C
      RETURN
      END
