C
C..........................................................................
C
      SUBROUTINE chkrio(ITO,TEXT,IER)
C
      INCLUDE 'ioin.incl'
C
      CHARACTER*(*) TEXT
C
Code ...
C
      IF (IER.EQ.NTXER0) THEN
         WRITE(ITO,1000) TEXT
      ELSE IF (IER.EQ.NTXER1) THEN
         WRITE(ITO,1010) TEXT
      ELSE IF (IER.EQ.NTXER2) THEN
         WRITE(ITO,1020) TEXT
      ELSE IF (IER.EQ.NTXER3) THEN
         WRITE(ITO,1030) TEXT
      ELSE IF (IER.EQ.NTXER4) THEN
         WRITE (ITO,1040) TEXT
      ELSE IF (IER.EQ.NTXER5) THEN
         WRITE (ITO,1050) TEXT
      ELSE
         WRITE(ITO,1999) TEXT,IER
      ENDIF
c
      call flusho(ito)
C
      RETURN
C
 1000 FORMAT(A,' no error')
 1010 FORMAT(A,' zero trip')
 1020 FORMAT(A,' buffer overflow')
 1030 FORMAT(A,' no variables')
 1040 FORMAT(A,' format overflow')
 1050 FORMAT(A,' output buffer overflow')
C
 1999 FORMAT(A,' = ',I8,' (unknown)')
C
      END
