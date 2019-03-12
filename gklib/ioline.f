C
C..........................................................................
C
      SUBROUTINE ioline(ITO,IER,MYTEXT)
C
C ... gjk@900731 - output
c     gjk@951105 - more "line-break markers" (], }, etc.)
C
      INCLUDE 'ioin.incl'
C
      PARAMETER (MAXLEN=ISYSNM-5,MINLEN=40)
C
      CHARACTER*(*) MYTEXT
      CHARACTER*(MAXLEN) MYLINE
C
Code ...
C
      LENTOT = LENGTH (MYTEXT) + 1
c
C ... mytext always ends with a space
C
      IF (LENTOT.LE.MAXLEN) THEN
        if (dollar) then
          WRITE (ITO,'(A,$)') MYTEXT(1:LENTOT)
        else
          WRITE (ITO,'(A)') MYTEXT(1:LENTOT)
        end if
      ELSE
        ISTART = 1
  100   CONTINUE
        ISTOP  = MIN( (ISTART + MAXLEN - 1), LENTOT)
        IF (ISTOP.LT.ISTART) GOTO 130
        DO 110 I=ISTOP,ISTART+MINLEN,-1
          IF (MYTEXT(I:I).EQ.' ' .OR.
     +        MYTEXT(I:I).EQ.')' .OR.
     +        MYTEXT(I:I).EQ.'-' .OR.
     +        MYTEXT(I:I).EQ.']' .OR.
     +        MYTEXT(I:I).EQ.'}' .OR.
     +        MYTEXT(I:I).EQ.';' .OR.
     +        MYTEXT(I:I).EQ.':' .OR.
     +        MYTEXT(I:I).EQ.',') THEN
            NOWSTP = I
            GOTO 120
          ENDIF
  110   CONTINUE
        NOWSTP = ISTOP
C
  120   CONTINUE
CCCD       WRITE (ITXOUT,*) ' start, stop ', ISTART, NOWSTP
        IF (ISTART.EQ.1) THEN
          WRITE (ITO,'(A)',ERR=888) MYTEXT(ISTART:NOWSTP)
        ELSE
          if (dollar .and. nowstp.ge.lentot) then
            WRITE (ITO,'(2X,A,$)',ERR=888) MYTEXT(ISTART:NOWSTP)
          else
            WRITE (ITO,'(2X,A)',ERR=888) MYTEXT(ISTART:NOWSTP)
          end if
        ENDIF
C
        ISTART = NOWSTP + 1
        IF (ISTART.LE.LENTOT) GOTO 100
  130   CONTINUE
C
      ENDIF
c
      call flusho (ito)
C
      RETURN
C
  888 CONTINUE
      IER = NTXER5
      RETURN
C
      END
