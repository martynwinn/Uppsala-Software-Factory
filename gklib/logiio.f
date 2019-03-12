C
C
C
      SUBROUTINE logiio(ITI,ITO,IER,TEXT,N,LOGI)
C
      INCLUDE 'ioin.incl'
C
      logical logi(*)
      CHARACTER*(*) TEXT
C
      CALL LOGIUO(ITO,IER,TEXT,N,LOGI)
      IF (IER.EQ.0) CALL LOGIII(ITI,IER,N,LOGI)
C
      RETURN
      END
