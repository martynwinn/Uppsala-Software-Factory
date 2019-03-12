c
c
c
      subroutine getcom (pro,line,munit,lecho)
c
c ... handle command I/O and macros
c
      implicit NONE
c
      integer macunt
      parameter (macunt=61)
c
      integer length,munit,i,j,ierr,nempty
c
      logical xinter,lmacro,ldone,lecho,lretry
c
      character pro*(*),line*(*),myline*256
c
code ...
c
      nempty = 0
c
      lmacro = (munit .ne. 5)
c
c ... initialisation macro ???
c
      if (pro .eq. '*INIT*') goto 130
c
    1 continue
c
      line = ' '
c
ccc      print *,' *1* munit, pro ',munit,' ',pro
c
      lmacro = (munit .ne. 5)
c
c ... read one line; handle end-of-macro
c
      if (lmacro) then
        read (munit,'(a)',err=110,end=120) line
        goto 100
c
  110   call errcon ('While reading line from macro file')
        goto 1
c
  120   call prompt (' ... End of macro file')
        close (munit)
        munit = munit - 1
        if (munit .lt. macunt) then
          call prompt (' ... Control returned to terminal')
          munit = 5
        end if
        goto 1
c
  100   continue
      else
        call prompt (pro)
        read (*,'(a)') line
      end if
c
ccc      print *,' *2* line ',line(1:20)
c
      if (length(line) .lt. 1) then
        nempty = nempty + 1
        if (nempty .ge. 50) then
          call errstp ('GETCOM - Fifty consecutive empty input lines')
        end if
        goto 1
      end if
c
      nempty = 0
c
  130 continue
      myline = line
      call remspa (myline)
c
ccc      print *,' *3* line ',line(1:20)
c
      call dohist (line,ldone)
c
ccc      print *,' *4* line ',line(1:20)
c
      if (ldone) goto 1
      myline = line
      call remspa (myline)
      call upcase (myline)
c
c ... if comment or macro or not interactive, echo command line
c
      if (myline(1:1) .eq. '!' .or. lmacro .or.
     +    (.not.xinter()) .or. lecho) then
        call textut (' Command > ',line)
      end if
c
ccc      print *,' *5* line ',line(1:20)
c
c ... comment line ?
c
      if (myline(1:1) .eq. '!') then
        goto 1
      end if
c
c ... system command
c
ccc      print *,' *6* line ',line(1:20)
c
      if (myline(1:1) .eq. '$') then
        i = index (line,'$')
        line (i:i) = ' '
        if (length(line) .lt. 1) then
          call errcon ('Empty system command ignored')
          goto 1
        end if
        call textut (' Spawn system command :',line)
ccc        myline = line
        call gksys (line(1:length(line)))
        goto 1
      end if
c
ccc      print *,' *7* line ',line(1:20)
c
c
      if (myline(1:1) .ne. '@') return
c
c ... new macro
c
      j=length(line)
      do i=1,j-1
        if (line(i:i) .eq. '@') then
          myline = line(i+1:)
          goto 200
        end if
      end do
      if (lmacro) then
        call errcon ('@ without filename inside macro')
        goto 1
      else
        myline = ' '
        call textin (' Name of macro file ?',myline)
      end if
c
  200 continue
      call remspa (myline)
      if (lmacro) then
        j = munit + 1
      else
        j = macunt
      end if
      lretry = (xinter() .and. (.not.lmacro) .and.
     +          pro .ne. '*INIT*')
      call xopxoa (j,myline,lretry,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening macro file')
        close (j)
c
c ... if startup macro doesn't exist, bail out
c
        if (pro .eq. '*INIT*') return
        goto 1
      end if
c
c ... macro opened successfully
c
      munit = j
      call textut (' ... Opened macro file :',myline)
      call ivalut (' ... On unit :',1,munit)
c
      goto 1
c
      end
