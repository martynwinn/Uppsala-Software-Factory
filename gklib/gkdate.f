c
c ===========================================================================
c
      subroutine gkdate (str24)
c
      implicit none
c
c --- GKDATE (STR24) => returns date and time
c
c --- G J Kleywegt @ 920403
c
      character str24*24
c
code ...
c
      call fdate (str24)
c
      return
      end
