c
c ... strings.f - string manipulation routines
c
c ... gerard j kleywegt @ 960409
c
c
c ========================================================================
c
      subroutine upcase (string)
c
c === UPCASE (STRING) = convert STRING into uppercase
c
c === G J Kleywegt @ 920221
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
      iup = ichar ('A') - ichar ('a')      
c
      do i=1,lenstr
        if (string(i:i).ge.'a' .and. string(i:i).le.'z')
     +    string (i:i) = char (iup + ichar(string(i:i)))
      end do
c
      return
      end
