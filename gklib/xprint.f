c
c ==========================================================================
c
      subroutine xprint (unit,text)
c
c --- XPRINT (UNIT,TEXT) - write TEXT to UNIT and flush
c                          if socket, send to O instead
c
c --- G J Kleywegt @ 920916
c
      implicit none
c
      integer unit,ltx,length
c
      logical xsocket
c
      character text*(*)
c
code ...
c
      if (xsocket()) then
        call oprint (unit,text)
      else
        ltx = length(text)
        if (ltx .lt. 1) then
          write (unit,*)
        else
          write (unit,'(1x,a)') text(1:ltx)
        end if
        call flusho(unit)
      end if
c
      return
      end
