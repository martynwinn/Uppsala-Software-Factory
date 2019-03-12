c
c ==========================================================================
c
      subroutine gkgarg (iarg,text,ierr)
c
c --- GKGARG (IARG,TEXT) - get command line argument nr IARG
c
c --- G J Kleywegt @ 950117
c
      implicit none
c
      integer maxarg,maxlen
      parameter (maxarg = 64, maxlen = 128)
c
      integer nargs,iarg,ierr
c
      character clargs(maxarg)*(maxlen),text*(*)
c
      common /coliar/ nargs
      common /colia2/ clargs
c
code ...
c
      if (iarg .gt. nargs) then
        ierr = -1
        return
      end if
c
      if (nargs .le. 0) then
        ierr = -2
        return
      end if
c
      text = clargs (iarg)
      ierr = 0
c
      return
      end
