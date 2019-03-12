c
c ==========================================================================
c
      subroutine oack (unit)
c
c --- OACK (UNIT) - send ACKNOWLEDGE signal to O
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
        write (unit,'(a)',err=999) 'ACK'
        call flusho (unit)
      end if
c
      return
c
  999 continue
      call errint (
     +  'OACK - could not send ACKNOWLEDGE signal to O')
c
      return
      end
