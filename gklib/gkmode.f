c
c ==========================================================================
c
      subroutine gkmode (mymode)
c
c --- GKMODE (MYMODE) - gets mode (interactive, batch or socket)
c
c --- G J Kleywegt @ 920915
c
      implicit none
c
      integer maxarg,maxlen
      parameter (maxarg = 64, maxlen = 128)
c
      integer nargs,i
c
      logical linter,lbatch,lsocket
c
      character clargs(maxarg)*(maxlen),mymode*(*),mode*20
      character line*(maxlen)
c
      common /coliar/ nargs
      common /colia2/ clargs
      common /promod/ linter,lbatch,lsocket
      common /promo2/ mode
c
code ...
c
      linter = .true.
      lbatch = .false.
      lsocket = .false.
c
      if (nargs .le. 0) goto 99
c
      do i=1,nargs
        line = clargs (i)
        call remspa (line)
        call upcase (line)
        if (line(1:2) .eq. '-I') then
          linter = .true.
          lbatch = .false.
          lsocket = .false.
          goto 99
        else if (line(1:2) .eq. '-B') then
          linter = .false.
          lbatch = .true.
          lsocket = .false.
          goto 99
        else if (line(1:2) .eq. '-S') then
          linter = .false.
          lbatch = .false.
          lsocket = .true.
          goto 99
        end if
      end do
c
   99 continue
      if (linter) then
        mode = 'interactive'
      else if (lbatch) then
        mode = 'batch'
      else
        mode = 'socket'
      end if
c
      mymode = mode
c
cc      print *,linter,lbatch,lsocket
c
      return
      end
