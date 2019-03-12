c
c ==========================================================================
c
      subroutine oready (unit)
c
c --- OREADY (UNIT) - send READY signal to O
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
        write (unit,'(a)',err=999) 'RDY'
        call flusho (unit)
      end if
c
      return
c
  999 continue
      call errint ('OREADY - could not send READY signal to O')
c
      return
      end
