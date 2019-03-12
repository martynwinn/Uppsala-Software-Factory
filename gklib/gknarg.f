c
c ==========================================================================
c
      subroutine gknarg (narg)
c
c --- GKNARG (NARG) - get nr of command line arguments
c
c --- G J Kleywegt @ 950117
c
      implicit none
c
      integer maxarg,maxlen
      parameter (maxarg = 64, maxlen = 128)
c
      integer nargs,narg
c
      common /coliar/ nargs
c
code ...
c
      narg = nargs
c
      return
      end
