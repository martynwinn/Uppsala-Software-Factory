c
c ==========================================================================
c
      subroutine gkerr (string)
c
c --- GKERR (string) - get system error message
c
c --- G J Kleywegt @ 930303
c
      implicit none
c
      character string*(*)
c
code ...
c
      string = ' '
      call gerror (string)
c
      return
      end
