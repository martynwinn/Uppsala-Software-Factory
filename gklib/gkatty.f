c
c ===========================================================================
c
      subroutine gkatty (iunit,ltty,ntty)
c
      implicit none
c
c --- GKATTY (...) => returns info about tty
c
c --- G J Kleywegt @ 920403
c
      integer iunit,length
c
      logical ltty,isatty
c
      character ntty*(*),ttynam*(99)
c
code ...
c
      ntty = ' '
      ltty = isatty (iunit)
      if (ltty) ntty = ttynam (iunit)
c
      if (length(ntty) .lt. 1) then
        ntty = 'Not a tty'
        ltty = .false.
      end if
c
      return
      end
