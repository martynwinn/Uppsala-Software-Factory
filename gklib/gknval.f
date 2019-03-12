c
c ===========================================================================
c
      subroutine gknval (nam,val,ierr)
c
      implicit none
c
c --- GKNVAL (...) => returns value of environment variable
c
c --- G J Kleywegt @ 920403
c
      integer ierr,length
c
      character nam*(*),val*(*)
c
code ...
c
      if (length(nam) .lt. 1) then
        ierr = -1
        return
      end if
c
      call getenv (nam,val)
c
      if (length(val) .lt. 1) then
        ierr = -2
        return
      end if
c
      ierr = 0
c
      return
      end
