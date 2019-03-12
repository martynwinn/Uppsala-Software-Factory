c
c ========================================================================
c
      subroutine str2r (str,xval,ierr)
c
c ... STR2R - read real from string
c
      implicit none
c
      real xval
c
      integer ierr
c
      character str*(*)
c
code ...
c
      read (str,*,err=99) xval
      ierr = 0
      return
c
   99 continue
      ierr = -1
      call errcon ('While reading real from input string')
c
      return
      end
