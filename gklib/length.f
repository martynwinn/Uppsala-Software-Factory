c
c ... functions.f - Fortran functions
c
c ... gerard j kleywegt @ 960409
c
c
c ===========================================================================
c
      integer function length (string)
c
      implicit none
c
      integer ls,ml,i
c
      character string*(*)
c
code ...
c
      ls = len (string)
      ml = 0
      if (ls .le. 0) goto 999
c
      do i=ls,1,-1
        if (string(i:i) .ne. ' ') then
          ml = i
          goto 999
        end if
      end do
c
  999 continue
      length = ml
c
      return
      end
