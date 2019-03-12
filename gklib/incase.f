c
c ========================================================================
c
      subroutine incase (str)
c
c ... INCASE (STRING) - set initial caps in names etc.
c     e.g. "a b john von dunder-klump ph.d." will become
c          "A B John Von Dunder-Klump Ph.D."
c
      implicit none
c
      integer i,ll,length
c
      character str*(*)
c
code ...
c
      ll = length(str)
c
      call pretty (str)
      call locase (str)
      call upcase (str(1:1))
c
      do i=1,ll-1
        if (index ('- .,',str(i:i)) .gt. 0) then
          call upcase (str(i+1:i+1))
        end if
      end do
c
      return
      end
