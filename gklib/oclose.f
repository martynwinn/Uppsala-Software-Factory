c
c ==========================================================================
c
      subroutine oclose (unit)
c
c --- OCLOSE (UNIT) - send CLOSE CONNECTION signal to O
c
c --- G J Kleywegt @ 920916
c
      implicit none
c
      integer unit
c
      logical xsocket
c
code ...
c
      if (xsocket()) then
        write (unit,'(a)',err=999) 'CLS'
        call flusho (unit)
      end if
c
      return
c
  999 continue
      call errint ('OCLOSE - could not send CLOSE signal to O')
c
      return
      end
