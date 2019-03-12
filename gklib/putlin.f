c
c=====================================================
c
      subroutine putlin (iunit,line,lstop,ierr)
c
c ... PUTLIN - write LINE (1:leng1) to stream IUNIT
c
      implicit none
c
      integer iunit,leng1,ierr
c
      logical lstop
c
      character line*(*)
c
code ...
c
      ierr = 0
      write (iunit,'(a)',err=999) line(1:leng1(line))
      return
c
  999 continue
      ierr = -1
      if (lstop) then
        call errcon ('While writing line')
        call jvalut (' Unit nr >',1,iunit)
        call textut (' Line >',line)
        call errstp ('PUTLIN - While writing line')
      else
        call errcon ('While writing line')
      end if
c
      return
      end
