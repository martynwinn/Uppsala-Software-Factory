C
C IO SUBROUTINES
C
C Written by and courtesy of:
C
C Rolf Boelens & Gerard Kleywegt
C Department of NMR Spectroscopy
C University of Utrecht
C Padualaan 8
C NL-358 CH Utrecht
C The Netherlands
C
C..........................................................................
C
C   standardised interactive I/O routines
C
C   881104 /rb
C   890915 /gk    (FVAL../GVAL..)
C   900731 /gk    (ERRor traps)
C   900731 /gk/rb (buffer via IOLINE)
C   900731 /gk    (JVAL.. routines)
c   930316 /gk    (if interactive and input, use $ in format)
c   991220 /gk    (lower-case subroutine names for linux)
C
C..........................................................................
C   standard input/output :
C
C    IVALIN (text,n,ival)
C    JVALIN (text,n,ival)
C    RVALIN (text,n,rval)
C    FVALIN (text,n,rval)
C    GVALIN (text,n,rval)
C    LOGIIN (text,n,lval)
C    ASCIIN (text,n,achr)
C    TEXTIN (text,line)            ..... one      line input
C
C      input  .... unit ITXINP
C      output .... unit ITXOUT
C
C..........................................................................
C   standard output only:
C
C    IVALUT (text,n,ival)
C    JVALUT (text,n,ival)
C    RVALUT (text,n,rval)
C    FVALUT (text,n,rval)
C    GVALUT (text,n,rval)
C    LOGIUT (text,n,lval)
C    ASCIUT (text,n,achr)
C    TEXTUT (text,line)
C
C      output .... unit ITXOUT
C
C..........................................................................
C   input/output :
C
C    IVALIO (iti,ito,ier,text,n,ival)
C    JVALIO (iti,ito,ier,text,n,ival)
C    RVALIO (iti,ito,ier,text,n,rval)
C    FVALIO (iti,ito,ier,text,n,rval)
C    GVALIO (iti,ito,ier,text,n,rval)
C    LOGIIO (iti,ito,ier,text,n,lval)
C    ASCIIO (iti,ito,ier,text,n,achr)
C    TEXTIO (iti,ito,ier,text,line)    ..... one      line input
C
C      input  .... unit iti
C      output .... unit ito
C      error  .... ier.ne.0
C..........................................................................
C   input only:
C
C    IVALII (iti,ier,n,ival)
C    JVALII (iti,ier,n,ival)
C    RVALII (iti,ier,n,rval)
C    FVALII (iti,ier,n,rval)
C    GVALII (iti,ier,n,rval)
C    LOGIII (iti,ier,n,lval)
C    ASCIII (iti,ier,n,achr)
C    TEXTII (iti,ier,line)    ..... one      line input
C
C      input .... unit iti
C      error  .... ier.ne.0
C..........................................................................
C   output only:
C
C    IVALUO (ito,ier,text,n,ival)
C    JVALUO (ito,ier,text,n,ival)
C    RVALUO (ito,ier,text,n,rval)
C    FVALUO (ito,ier,text,n,rval)
C    GVALUO (ito,ier,text,n,rval)
C    LOGIUO (ito,ier,text,n,lval)
C    ASCIUO (ito,ier,text,n,achr)
C    TEXTUO (ito,ier,text,line)    ..... one      line output
C
C      output .... unit ito
C      error  .... ier.ne.0
C
C..........................................................................
C   error:
C
C    CHKRIO(ito,text,ier)
C
C      output .... unit ito
C      ier    .... error to be interpreted
C
C..........................................................................
C   output:
C
C    IOLINE (ITO,IER,TEXT)
C
C      output .... unit ito
C      ier    .... returned error number (or unchanged)
C      text   .... string to be printed
C
c    FLUSHO (IUNIT)
c
c      iunit  .... output channel to be flushed
c
C..........................................................................
C   general:
C
C    output
C      text   .... output   character*(*) string
C
C    input
C      n      .... # input  variables
C      ival   .... array of integer     variables
C      rval   .... array of real        variables
C      lval   .... array of logical       variables
C      achr   .... array of character*(*) variables
C      line   ....          character*(*) variable
C
C
C   restriction : only one line (LNBFIO characters) input/output
C   900731 ...... NO MORE ! Implemented routine IOLINE
C
C..........................................................................
C
      SUBROUTINE ivalin(TEXT,N,IVALUE)
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
      dollar = xinter()
      CALL IVALUO(ITXOUT,IER,TEXT,N,IVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' ivaluo -- i/o error',IER)
        RETURN
      ENDIF
C
      CALL IVALII(ITXINP,IER,N,IVALUE)
      IF (IER.NE.NTXER0) THEN
        CALL CHKRIO(ITXERR,' ivalii -- i/o error',IER)
      ENDIF
C
      RETURN
      END
