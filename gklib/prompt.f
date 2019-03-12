c
c ========================================================================
c
      subroutine prompt (string)
c
      implicit none
c
      integer ls,length
c
      character*(*) string
c
code ...
c
      ls = length (string)
      if (ls .le. 0) return
c
      if (string(1:1) .eq. '0') then
        write (*,'(/1x,a)') string(2:ls)
      else if (string(1:1) .eq. '1') then
        write (*,'(//1x,a)') string(2:ls)
      else if (string(1:1) .eq. '$') then
        write (*,'(1x,a,a1,$)') string(2:ls),' '
      else
        write (*,'(a)') string(1:ls)
      end if
c
      return
      end
