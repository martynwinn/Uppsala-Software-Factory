c
c ========================================================================
c
      logical function isitok (string)
c
c === G J Kleywegt @ 920306,19
c
      implicit none
c
      character string*(*),answer*1
c
code ...
c
      answer = 'Y'
      call asciin (string,1,answer)
      isitok = (.not.(answer.eq.'n' .or. answer.eq.'N'))
c
      return
      end
