c
c ===========================================================================
c
      subroutine remspa (string)
c
c === REMSPA (STRING) = remove spaces from STRING
c
c === G J Kleywegt @ 910606,11
c
      implicit none
c
      integer lenstr,icnt,length
c
      character string*(*),space*1
c
      data space/' '/
c
code ...
c
      lenstr = length(string)
      icnt = 0
c
   20 continue
        icnt = icnt + 1
        if (icnt.gt.lenstr) goto 999
c
        if (string(icnt:icnt) .eq. space) then
          string(icnt:lenstr-1) = string(icnt+1:lenstr)
          string (lenstr:lenstr) = space
          icnt   = icnt - 1
          lenstr = lenstr - 1
        end if
c
      goto 20
c
c ... exit point of this routine
c
  999 continue
c
      return
      end
