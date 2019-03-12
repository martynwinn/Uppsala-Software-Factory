c
c ==========================================================================
c
      subroutine gkargs ()
c
c --- GKARGS () - get command line arguments (if any) and store in common block
c
c --- G J Kleywegt @ 920915
c
      implicit none
c
      integer maxarg,maxlen
      parameter (maxarg = 64, maxlen = 128)
c
      integer nargs,iargc,i
c
      character clargs(maxarg)*(maxlen)
c
      common /coliar/ nargs
      common /colia2/ clargs
c
      data nargs /0/
c
code ...
c
      nargs = iargc ()
c
      if (nargs .gt. 0) then
        if (nargs .gt. maxarg) then
          call errcon ('GKARGS - too many command line arguments')
          nargs = maxarg
        end if
        do i=1,nargs
          call getarg (i,clargs(i))
        end do
      end if
c
      return
      end
