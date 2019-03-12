c
c ===========================================================================
c
      subroutine errint (why)
c
c ... error is fatal if socket or batch, but not if interactive
c
      implicit none
c
      logical xinter
c
      character why*(*)
c
code ...
c
      if (xinter()) then
        call errcon (why)
      else
        call errstp (why)
      end if
c
      return
      end
