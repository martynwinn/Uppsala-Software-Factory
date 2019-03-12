c
c
c
      subroutine dohist (line,ldone)
c
c ... DOHIST - manage command history
c
c ... # ? = list history
c     # on = history on
c     # of = history off
c     # <num> = if <num> positive, repeat command nr <num>
c               if <num> = 0, repeat previous command
c               if <num> = -1, repeat penultimate command, etc.
c     NOTE: rest of line is appended
c
      implicit none
c
      integer i,ierr,idum,iptr,imin,length,ll
c
      logical ldone,lex
c
      character line*(*),myline*256,junk*20
c
      integer maxhis,hislen
      parameter (maxhis=500,hislen=256)
      integer numhis,hisctr(maxhis)
      logical ldohis
      character histor(maxhis)*(hislen)
      common /dohis1/ numhis,hisctr
      common /dohis2/ ldohis
      common /dohis3/ histor
c
code ...
c
 6000 format (1x,i8,2x,a)
 6010 format (1x,a)
c
      ldone = .false.
c
      if (length(line) .lt. 1) return
c
ccc      if (line(1:1) .eq. '!') return
c
      myline = line
      call pretty (myline)
      call upcase (myline)
c
      ldone = .true.
c
      if (myline(1:6) .eq. '*INIT*') then
        numhis = 0
        ldohis = .true.
        do i=1,maxhis
          hisctr (i) = -1
          histor (i) = ' '
        end do
        return
      end if
c
      if (myline(1:1) .eq. '#') then
        if (index(myline(2:),'?') .gt. 0 .or.
     +      index(myline(2:),'!') .gt. 0) then
          lex = (index(myline(2:),'!') .gt. 0)
          imin = numhis + 1
          iptr = -1
          do i=1,maxhis
            if (hisctr(i).lt.imin .and. hisctr(i).gt.0) then
              imin = hisctr(i)
              iptr = i
            end if
          end do
          if (iptr .gt. 0) then
            do i=iptr,maxhis
              if (hisctr(i).gt.0 .and.
     +            length(histor(i)).gt.0) then
                if (lex) then
                  write (*,6010)
     +              histor(i)(1:length(histor(i)))
                else
                  write (*,6000) hisctr(i),
     +              histor(i)(1:length(histor(i)))
                end if
              end if
            end do
            if (iptr .gt. 1) then
              do i=1,iptr-1
                if (hisctr(i).gt.0 .and.
     +              length(histor(i)).gt.0) then
                  if (lex) then
                    write (*,6010)
     +                histor(i)(1:length(histor(i)))
                  else
                    write (*,6000) hisctr(i),
     +                histor(i)(1:length(histor(i)))
                  end if
                end if
              end do
            end if
          else
            call errcon ('No commands in history list yet')
          end if
        else if (index(myline(2:5),'ON') .gt. 0) then
          ldohis = .true.
          call prompt (' Command history ON')
        else if (index(myline(2:5),'OF') .gt. 0) then
          ldohis = .false.
          call prompt (' Command history OFF')
        else
c
          if (index(myline(2:5),'#') .gt. 0) then
            idum = 0
            ierr = 0
            junk = '#'
            ll = 1
            iptr = 1 + index(myline(2:5),'#')
          else
            call str2i (myline(2:),idum,ierr)
            if (ierr .eq. 0) then
              write (junk,*) idum
              call remspa (junk)
              ll = length(junk)
              iptr = index(line,junk(1:ll))
            end if
          end if
c
          if (ierr .eq. 0) then
            if (idum.le.0) idum = numhis+idum
c
            do i=1,maxhis
              if (idum .eq. hisctr(i)) then
                if (iptr .le. 0) then
                  line = histor(i)
                else
                  line = histor(i)(1:length(histor(i))) //
     +                   line(iptr+ll:)
                end if
                goto 100
              end if
            end do
            call errcon ('History command number not found')
          else
            call errcon ('Invalid history command number')
          end if
        end if
        return
      end if
c
  100 continue
      ldone = .false.
      if (.not. ldohis) return
      numhis = numhis + 1
      imin = hisctr(1)
      iptr = 1
      do i=2,maxhis
        if (hisctr(i).lt.imin) then
          imin = hisctr(i)
          iptr = i
        end if
      end do
      histor (iptr) = line
      hisctr (iptr) = numhis
c
      return
      end
