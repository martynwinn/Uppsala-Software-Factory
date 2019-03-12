c
c ========================================================================
c
      subroutine locase (string)
c
c === LOCASE (STRING) = convert STRING into lowercase
c
c === G J Kleywegt @ 920916
c
      implicit none
c
      integer lenstr,iup,i,length
c
      character string*(*)
c
code ...
c
      lenstr = length(string)
      if (lenstr .le. 0) return
c
      iup = ichar ('a') - ichar ('A')      
c
      do i=1,lenstr
        if (string(i:i).ge.'A' .and. string(i:i).le.'Z')
     +    string (i:i) = char (iup + ichar(string(i:i)))
      end do
c
      return
      end
