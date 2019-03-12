      SUBROUTINE xopen (IFUNIT,FNAME,FSTAT,FFORM,LRETRY,IERR)
C
C -----------------------------------------------------------------------------
C
C --- XOPEN (IFUNIT,FNAME,FSTAT,FFORM,LRETRY,IERR)
C
C     general subroutine for opening files; should even work with Unix !
C
C ... Gerard J Kleywegt @ 900504/900907/920403/920820/930303/981021
c
c ... 991220 - don't support CARRIAGECONTROL, READONLY any longer (linux)
C
C --- ALTERNATIVE ENTRIES -----------------------------------------------------
C
C     XOPEN  (IFUNIT,FNAME,FSTAT,FFORM,LRETRY,IERR)
C     ...... the general one; will suit most purposes
C
C     XOLASC (IFUNIT,FNAME,LRETRY,IERR)
C     ...... quicky to open an OLD ASCII ('O', 'F') file
C
C     XEXIST (IFUNIT,FNAME,LRETRY,IERR)
C     ...... same as XOLASC, but no error messages printed
C            (value of LRETRY is ignored !)
C
C     XOPXNA (IFUNIT,FNAME,LRETRY,IERR)
C     ...... quicky to open a NEW ASCII ('N', 'F') file
C
C     XOPXUA (IFUNIT,FNAME,LRETRY,IERR)
C     ...... quicky to open an UNKNOWN ASCII ('U', 'F') file
C
C     XOPXOB (IFUNIT,FNAME,LRETRY,IERR)
C     ...... quicky to open an OLD BINARY ('O', 'U') file
C
C     XOPXNB (IFUNIT,FNAME,LRETRY,IERR)
C     ...... quicky to open a NEW BINARY ('N', 'U') file
C
C     XOPXUB (IFUNIT,FNAME,LRETRY,IERR)
C     ...... quicky to open an UNKNOWN BINARY ('U', 'U') file
C
C     XOPENL (IFUNIT,FNAME,FSTAT,FFORM,MRECSI,LRETRY,IERR)
C     ...... as XOPEN; specify RECORDSIZE
C
C     XOPENA (IFUNIT,FNAME,FSTAT,FFORM,FACCS,LRETRY,IERR)
C     ...... as XOPEN; specify ACCESS
C
C     XOPENC (IFUNIT,FNAME,FSTAT,FFORM,FCARC,LRETRY,IERR)
C     ...... as XOPEN; specify CARRIAGECONTROL
C
C     XOPENR (IFUNIT,FNAME,FSTAT,FFORM,MRONLY,LRETRY,IERR)
C     ...... as XOPEN; specify READONLY
C
C     XOPENX (IFUNIT,FNAME,FSTAT,FFORM,FACCS,FCARC,
C    +        MRECSI,MRONLY,LRETRY,MIOS,IERR)
C     ...... as XOPEN; specify (virtually) everything
C
C
C --- PARAMETERS --------------------------------------------------------------
C
C   * IFUNIT = unit number              IN RANGE 1-99
C
C   * FNAME  = filename                 OPTIONAL
C
C   * FSTAT  = status                   may be 'O' (old), 'N' (new),
C                                       'U' (unknown) or 'S' (scratch)
C                                       DEFAULT='UNKNOWN'
C
C   * FFORM  = (un)formatted file       may be 'U' or 'F'
C                                       DEFAULT='FORMATTED'
C
C   * LRETRY = retry flag
C
C   * IERR   = returned error code      will be 0 if okay
C
C
C --- ALTERNATIVE PARAMETERS --------------------------------------------------
C
C   * MRECSI = record length            0 = DON'T USE = DEFAULT
C
C   * FACCS  = access                   may be 'D' (direct), 'K' (keyed),
C                                       'A' (append) or 'S' (sequential);
C                                       'D' & 'K' only with UNFORM files
C                                       DEFAULT='SEQUENTIAL'
C
C   * FCARC  = carriage control         may 'F' (fortran), 'L' (list) or
C                                       'N' (none)
C                                       DEFAULT='LIST'(ascii)/'NONE'(binary)
C
C   * MRONLY = readonly specifier       DEFAULT=.FALSE.
C                                       only for OLD/UNKNOWN files !
C
C   * MIOS   = value of IOSTAT          RETURNED
C
C
C --- USE OF LRETRY -----------------------------------------------------------
C
C     If LRETRY is TRUE (should only be done in interactive programs !),
C     then, if the opening fails, the following happens:
C        OLD/NAME => prompt user for other file name
C        OLD/UNIT => nothing
C        NEW/NAME => ask if okay to open file as OLD (Unix !!!)
C                    if not ask if other file name
C        NEW/UNIT => nothing
C        UNKNOWN  => treated the same as OLD
C        SCRATCH  => treated the same as NEW
C
C
C --- ERROR CODES -------------------------------------------------------------
C
C     Returned error codes :
C     IERR =  0 => opening was successful !
C     IERR =  1 => unit number less than or equal to zero
C     IERR =  2 => error while opening OLD or UNKNOWN by NAME
C     IERR =  3 => error while opening OLD or UNKNOWN by UNIT
C     IERR =  4 => error while opening NEW by NAME
C     IERR =  5 => error while opening NEW by UNIT
C     IERR = 10 => invalid record length
C     IERR = 11 => KEYED or DIRECT access attempted with formatted file
C     IERR = 12 => FORTRAN or LIST carriage control attempted with
C                  unformatted file
C
C
C -----------------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER NILL
      PARAMETER (NILL=0)
C
      INTEGER MIOS,IOS,NRECSI,IERR,IFUNIT,MRECSI,MYLEN,MYSLEN
      INTEGER MYFLEN,MYALEN,MYCLEN,IOGENV,li,leng1
C
      LOGICAL LRETRY,LRONLY,MRONLY
      LOGICAL LUSURE,lprint,lchangename
C
      CHARACTER FNAME*(*),FSTAT*(*),FFORM*(*),FACCS*(*),FCARC*(*)
      CHARACTER MYNAME*120,MYSTAT*30,MYFORM*30,MYMESS*250
      CHARACTER MYACCS*30,MYCARC*30,MENTRY*6,info*250
      character sysmsg*132
C
Code ...
C
      MYNAME = FNAME
      MYSTAT = FSTAT
      MYFORM = FFORM
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = ' XOPEN'
      lprint = .true.
      GOTO 69
C
C ... ALTERNATIVE ENTRIES ...
C
C ... OLD ASCII
C
      ENTRY xolasc (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'OLD'
      MYFORM = 'F'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XOLASC'
      lprint = .true.
      GOTO 69
C
C ... XEXIST (OLD ASCII)
C
      ENTRY xexist (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'OLD'
      MYFORM = 'F'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XEXIST'
      lprint = .false.
      GOTO 69
C
C ... NEW ASCII
C
      ENTRY xopxna (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'NEW'
      MYFORM = 'F'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XOPXNA'
      lprint = .true.
      GOTO 69
C
C ... UNKNOWN ASCII
C
      ENTRY xopxua (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'UNKNOWN'
      MYFORM = 'F'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XOPXUA'
      lprint = .true.
      GOTO 69
C
C ... OLD BINARY
C
      ENTRY xopxob (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'OLD'
      MYFORM = 'U'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XOPXOB'
      lprint = .true.
      GOTO 69
C
C ... NEW BINARY
C
      ENTRY xopxnb (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'NEW'
      MYFORM = 'U'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XOPXNB'
      lprint = .true.
      GOTO 69
C
C ... UNKNOWN BINARY
C
      ENTRY xopxub (IFUNIT,FNAME,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = 'UNKNOWN'
      MYFORM = 'U'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      MENTRY = 'XOPXUB'
      lprint = .true.
      GOTO 69
C
C ... specify record length
C
      ENTRY xopenl (IFUNIT,FNAME,FSTAT,FFORM,MRECSI,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = FSTAT
      MYFORM = FFORM
      MENTRY = 'XOPENL'
      NRECSI = MRECSI
      MYACCS = ' '
      MYCARC = ' '
      LRONLY = .FALSE.
      lprint = .true.
      GOTO 69
C
C ... specify access
C
      ENTRY xopena (IFUNIT,FNAME,FSTAT,FFORM,FACCS,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = FSTAT
      MYFORM = FFORM
      MENTRY = 'XOPENA'
      MYACCS = FACCS
      NRECSI = NILL
      MYCARC = ' '
      LRONLY = .FALSE.
      lprint = .true.
      GOTO 69
C
C ... specify carriage control
C
      ENTRY xopenc (IFUNIT,FNAME,FSTAT,FFORM,FCARC,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = FSTAT
      MYFORM = FFORM
      MENTRY = 'XOPENC'
      MYCARC = FCARC
      NRECSI = NILL
      MYACCS = ' '
      LRONLY = .FALSE.
      lprint = .true.
      call errcon ('CARRIAGECONTROL no longer supported (XOPENC)')
      GOTO 69
C
C ... specify readonly
C
      ENTRY xopenr (IFUNIT,FNAME,FSTAT,FFORM,MRONLY,LRETRY,IERR)
      MYNAME = FNAME
      MYSTAT = FSTAT
      MYFORM = FFORM
      MENTRY = 'XOPENR'
      NRECSI = NILL
      MYACCS = ' '
      MYCARC = ' '
ccc      LRONLY = MRONLY
      LRONLY = .false.
      call errcon ('READONLY no longer supported (XOPENR)')
      lprint = .true.
      GOTO 69
C
C ... specify all attributes
C
      ENTRY xopenx (IFUNIT,FNAME,FSTAT,FFORM,FACCS,FCARC,
     +              MRECSI,MRONLY,LRETRY,MIOS,IERR)
      MYNAME = FNAME
      MYSTAT = FSTAT
      MYFORM = FFORM
      MENTRY = 'XOPENX'
      MYACCS = FACCS
      MYCARC = FCARC
      NRECSI = MRECSI
ccc      LRONLY = MRONLY
      LRONLY = .false.
      call errcon ('READONLY no longer supported (XOPENR)')
      lprint = .true.
      call errcon ('CARRIAGECONTROL no longer supported (XOPENX)')
      GOTO 69
C
Checks ...
C
   69 CONTINUE
      IERR = NILL
      lchangename = .false.
C
C ... unit number
C
      IF (IFUNIT.LE.NILL) THEN
        IERR = 1
        WRITE (MYMESS,6000) MENTRY, IFUNIT
        CALL ERRCON (MYMESS)
        RETURN
      ENDIF
C
C ... record length
C
      IF (NRECSI.LT.NILL) THEN
        IERR = 10
        WRITE (MYMESS,6030) MENTRY, NRECSI
        CALL ERRCON (MYMESS)
        RETURN
      ENDIF
C
C ... status
C
      CALL UPCASE (MYSTAT)
      IF (MYSTAT(1:1).EQ.'O') THEN
        MYSTAT = 'OLD'
      ELSE IF (MYSTAT(1:1).EQ.'N') THEN
        MYSTAT = 'NEW'
      ELSE IF (MYSTAT(1:1).EQ.'S') THEN
        MYSTAT = 'SCRATCH'
      ELSE
        MYSTAT = 'UNKNOWN'
      ENDIF
C
C ... form
C
      CALL UPCASE (MYFORM)
      IF (MYFORM(1:1).EQ.'U') THEN
        MYFORM = 'UNFORMATTED'
      ELSE
        MYFORM = 'FORMATTED'
      ENDIF
C
C ... access
C
      CALL UPCASE (MYACCS)
      IF (MYACCS(1:1).EQ.'D') THEN
        MYACCS = 'DIRECT'
      ELSE IF (MYACCS(1:1).EQ.'K') THEN
        MYACCS = 'KEYED'
      ELSE IF (MYACCS(1:1).EQ.'A') THEN
        MYACCS = 'APPEND'
      ELSE
        MYACCS = 'SEQUENTIAL'
      ENDIF
C
      IF (MYFORM(1:1).NE.'U' .AND.
     +    (MYACCS(1:1).EQ.'D' .OR. MYACCS(1:1).EQ.'K')) THEN
        IERR = 11
        WRITE (MYMESS,6040) MENTRY, MYACCS(1:LENG1(MYACCS))
        CALL ERRCON (MYMESS)
        RETURN
      ENDIF
C
C ... carriage control
C
c      CALL UPCASE (MYCARC)
c      IF (MYCARC(1:1).EQ.'L') THEN
c        MYCARC = 'LIST'
c      ELSE IF (MYCARC(1:1).EQ.'N') THEN
c        MYCARC = 'NONE'
c      ELSE IF (MYCARC(1:1).EQ.'F') THEN
c        MYCARC = 'FORTRAN'
c      ELSE
c        IF (MYFORM(1:1).EQ.'F') THEN
c          MYCARC = 'LIST'
c        ELSE
c          MYCARC = 'NONE'
c        ENDIF
c      ENDIF
C
c      IF (MYFORM(1:1).EQ.'U' .AND.
c     +    (MYCARC(1:1).EQ.'L' .OR. MYCARC(1:1).EQ.'F')) THEN
c        IERR = 12
c        WRITE (MYMESS,6050) MENTRY, MYCARC(1:LENG1(MYCARC))
c        CALL ERRCON (MYMESS)
c        RETURN
c      ENDIF
c
c ... first close the unit (just to be sure)
c
      close (ifunit)
C
C ... the actual opening bit ..................................................
C
    5 CONTINUE
      IERR = NILL
C
      MYLEN  = LENG1(MYNAME)
      MYSLEN = LENG1(MYSTAT)
      MYFLEN = LENG1(MYFORM)
      MYALEN = LENG1(MYACCS)
      MYCLEN = LENG1(MYCARC)
c
      write (info,6900) ifunit,mystat(1:myslen),
     +  mycarc(1:myclen),myform(1:myflen),
     +  myaccs(1:myalen)
      call pretty (info)
      li = leng1(info)
      if (lronly) then
        info(li+1:) = ' READONLY '
      end if
      li = leng1(info)
      if (nrecsi .gt. 0) then
        write (info(li+1:),'(a,i10)') ' RECL=',nrecsi
      end if
      call pretty (info)
      li = leng1 (info)
c
 6900 format (' UNIT=',i6,' STATUS=',a,
     +  ' CAR_CONTROL=',a,' FORM=',a,
     +  ' ACCESS=',a)
C
CCCD     CALL TEXTUT ('  entry :',MENTRY)
CCCD     CALL IVALUT ('  unit  :',1,IFUNIT)
CCCD     IF (MYLEN.GT.0) CALL TEXTUT ('  file  :',MYNAME)
CCCD     CALL TEXTUT ('  stat  :',MYSTAT)
CCCD     CALL TEXTUT ('  carc  :',MYCARC)
CCCD     CALL TEXTUT ('  accs  :',MYACCS)
CCCD     CALL TEXTUT ('  form  :',MYFORM)
CCCD     CALL IVALUT ('  recl  :',1,NRECSI)
CCCD     CALL LOGIUT ('  reado :',1,LRONLY)
CCCD     CALL LOGIUT ('  retry :',1,LRETRY)
C
      IF (MYSTAT(1:3).EQ.'OLD' .OR. MYSTAT(1:3).EQ.'UNK') THEN
C
C --- OLD or UNK/NAME
C
        IF (MYLEN.GT.NILL) THEN
          IF (NRECSI.GT.NILL) THEN
            IF (LRONLY) THEN
              OPEN (IFUNIT,
     +          FILE=MYNAME,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          RECL=NRECSI,
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
ccc     +          READONLY,
     +          IOSTAT=IOS,
     +          ERR=10)
            ELSE
              OPEN (IFUNIT,
     +          FILE=MYNAME,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          RECL=NRECSI,
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +          IOSTAT=IOS,
     +          ERR=10)
            ENDIF
          ELSE
            IF (LRONLY) THEN
              OPEN (IFUNIT,
     +          FILE=MYNAME,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
ccc     +          READONLY,
     +          IOSTAT=IOS,
     +          ERR=10)
            ELSE
              OPEN (IFUNIT,
     +          FILE=MYNAME,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +          IOSTAT=IOS,
     +          ERR=10)
            ENDIF
          ENDIF
          IERR = NILL
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          fname = myname
          RETURN
C ... ERROR
   10     CONTINUE
          IERR = 2
          if (lprint) then
            WRITE (MYMESS,6010) MENTRY, IOS, MYSTAT(1:MYSLEN),
     +                          MYNAME(1:MYLEN)
            CALL ERRCON (MYMESS)
            call textut (' OPEN :',info)
            call gkerr (sysmsg)
            call textut (' Error :',sysmsg)
            IF (LRETRY) THEN
              IF (LUSURE('Open with other name')) THEN
                CALL textin (' File ?',MYNAME)
                lchangename = .true.
                GOTO 5
              ENDIF
            ENDIF
          end if
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          fname = myname
          RETURN
C
C --- OLD or UNK/UNIT
C
        ELSE
          IF (IOGENV (IFUNIT,MYNAME) .GT. NILL) THEN
            IF (NRECSI.GT.NILL) THEN
              IF (LRONLY) THEN
                OPEN (IFUNIT,
     +            FILE=MYNAME,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            RECL=NRECSI,
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
ccc     +            READONLY,
     +            IOSTAT=IOS,
     +            ERR=20)
              ELSE
                OPEN (IFUNIT,
     +            FILE=MYNAME,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            RECL=NRECSI,
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +            IOSTAT=IOS,
     +            ERR=20)
              ENDIF
            ELSE
              IF (LRONLY) THEN
                OPEN (IFUNIT,
     +            FILE=MYNAME,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
ccc     +            READONLY,
     +            IOSTAT=IOS,
     +            ERR=20)
              ELSE
                OPEN (IFUNIT,
     +            FILE=MYNAME,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +            IOSTAT=IOS,
     +            ERR=20)
              ENDIF
            ENDIF
          ELSE
            IF (NRECSI.GT.NILL) THEN
              IF (LRONLY) THEN
                OPEN (IFUNIT,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            RECL=NRECSI,
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
ccc     +            READONLY,
     +            IOSTAT=IOS,
     +            ERR=20)
              ELSE
                OPEN (IFUNIT,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            RECL=NRECSI,
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +            IOSTAT=IOS,
     +            ERR=20)
              ENDIF
            ELSE
              IF (LRONLY) THEN
                OPEN (IFUNIT,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
ccc     +            READONLY,
     +            IOSTAT=IOS,
     +            ERR=20)
              ELSE
                OPEN (IFUNIT,
     +            STATUS=MYSTAT(1:MYSLEN),
     +            FORM=MYFORM(1:MYFLEN),
     +            ACCESS=MYACCS(1:MYALEN),
ccc     +            CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +            IOSTAT=IOS,
     +            ERR=20)
              ENDIF
            ENDIF
          ENDIF
          IERR = NILL
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          RETURN
C ... ERROR
   20     CONTINUE
          IERR = 3
          if (lprint) then
            WRITE (MYMESS,6020) MENTRY, IOS, MYSTAT(1:MYSLEN),IFUNIT
            CALL ERRCON (MYMESS)
            call textut (' OPEN :',info)
            call gkerr (sysmsg)
            call textut (' Error :',sysmsg)
          end if
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          RETURN
        ENDIF
C
      ELSE IF (MYSTAT(1:3).EQ.'NEW' .OR. MYSTAT(1:3).EQ.'SCR') THEN
C
C --- NEW or SCR/NAME
C
        IF (MYLEN.GT.NILL) THEN
          IF (NRECSI.GT.NILL) THEN
            OPEN (IFUNIT,
     +        FILE=MYNAME,
     +        STATUS=MYSTAT(1:MYSLEN),
     +        FORM=MYFORM(1:MYFLEN),
     +        RECL=NRECSI,
     +        ACCESS=MYACCS(1:MYALEN),
ccc     +        CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +        IOSTAT=IOS,
     +        ERR=60)
          ELSE
            OPEN (IFUNIT,
     +        FILE=MYNAME,
     +        STATUS=MYSTAT(1:MYSLEN),
     +        FORM=MYFORM(1:MYFLEN),
     +        ACCESS=MYACCS(1:MYALEN),
ccc     +        CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +        IOSTAT=IOS,
     +        ERR=60)
          ENDIF
          IERR = NILL
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          fname = myname
          RETURN
C ... ERROR
   60     CONTINUE
          IERR = 4
          if (lprint) then
            WRITE (MYMESS,6010) MENTRY, IOS, MYSTAT(1:MYSLEN),
     +                          MYNAME(1:MYLEN)
            CALL ERRCON (MYMESS)
            call textut (' OPEN :',info)
            call gkerr (sysmsg)
            call textut (' Error :',sysmsg)
            IF (LRETRY) THEN
              IF (LUSURE('Open file as OLD')) THEN
                MYSTAT = 'OLD'
                GOTO 5
              ENDIF
              IF (LUSURE('Open with other name')) THEN
                CALL textin (' File ?',MYNAME)
                lchangename = .true.
                GOTO 5
              ENDIF
            ENDIF
          end if
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          fname = myname
          RETURN
C
C --- NEW or SCR/UNIT
C
        ELSE
          IF (IOGENV (IFUNIT,MYNAME) .GT. NILL) THEN
            IF (NRECSI.GT.NILL) THEN
              OPEN (IFUNIT,
     +          FILE=MYNAME,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          RECL=NRECSI,
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +          IOSTAT=IOS,
     +          ERR=70)
            ELSE
              OPEN (IFUNIT,
     +          FILE=MYNAME,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +          IOSTAT=IOS,
     +          ERR=70)
            ENDIF
          ELSE
            IF (NRECSI.GT.NILL) THEN
              OPEN (IFUNIT,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          RECL=NRECSI,
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +          IOSTAT=IOS,
     +          ERR=70)
            ELSE
              OPEN (IFUNIT,
     +          STATUS=MYSTAT(1:MYSLEN),
     +          FORM=MYFORM(1:MYFLEN),
     +          ACCESS=MYACCS(1:MYALEN),
ccc     +          CARRIAGECONTROL=MYCARC(1:MYCLEN),
     +          IOSTAT=IOS,
     +          ERR=70)
            ENDIF
          ENDIF
          IERR = NILL
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          RETURN
C ... ERROR
   70     CONTINUE
          IERR = 5
          if (lprint) then
            WRITE (MYMESS,6020) MENTRY, IOS, MYSTAT(1:MYSLEN),IFUNIT
            CALL ERRCON (MYMESS)
            call textut (' OPEN :',info)
            call gkerr (sysmsg)
            call textut (' Error :',sysmsg)
            IF (LRETRY) THEN
              IF (LUSURE('Open file as OLD')) THEN
                MYSTAT = 'OLD'
                GOTO 5
              ENDIF
            ENDIF
          end if
          IF (MENTRY.EQ.'XOPENX') MIOS = IOS
          RETURN
        ENDIF
      ENDIF
C
      IF (MENTRY.EQ.'XOPENX') MIOS = IOS
c
      if (lchangename) FNAME = MYNAME
c
      RETURN
C
C --- formats
C
 6000 FORMAT(A,' - invalid unit number : ',I3)
 6010 FORMAT(A,' - error # ',I3,' while opening ',A,' file : ',A)
 6020 FORMAT(A,' - error # ',I3,' while opening ',A,' unit : ',I3)
 6030 FORMAT(A,' - invalid record length : ',I3)
 6040 FORMAT(A,' - ',A,' access only with unformatted files')
 6050 FORMAT(A,' - ',A,' carriage control only with formatted files')
C
      END
