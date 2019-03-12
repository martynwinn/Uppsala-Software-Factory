c
c ========================================================================
c
      subroutine str2i (str,ival,ierr)
c
c ... STR2I - read integer from string
c
      implicit none
c
      integer ival,ierr
c
      character str*(*)
c
code ...
c
      read (str,*,err=99) ival
      ierr = 0
      return
c
   99 continue
      ierr = -1
      call errcon ('While reading integer from input string')
c
      return
      end
