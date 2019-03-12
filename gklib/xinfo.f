c
c ==========================================================================
c
      subroutine xinfo (unit,text)
c
c --- XINFO (UNIT,TEXT) - write TEXT to UNIT and flush
c                         unless socket
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
      if (.not. xsocket()) then
        ltx = length(text)
        if (ltx .lt. 1) then
          write (unit,*)
        else
          write (unit,'(1x,a)') text(1:ltx)
        end if
      end if
c
      return
      end
