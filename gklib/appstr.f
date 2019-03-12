c
c ========================================================================
c
      subroutine appstr (str1,str2)
c
      implicit none
c
c === APPSTR (STR1,STR2) = append STR2 to STR1
c
c === G J Kleywegt @ 910606/920403
c
      integer length
c
      character*(*) str1,str2
c
code ...
c
      if (length(str1) .lt. 1) then
        str1 = str2
      else
        str1 = str1(1:length(str1))//str2
      end if
c
      return
      end
