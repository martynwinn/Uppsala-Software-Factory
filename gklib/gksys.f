c
c ==========================================================================
c
      subroutine gksys (str)
c
c --- GKSYS (str) - execute system command
c
c --- G J Kleywegt @ 930608
c
      implicit none
c
      integer length
c
      character str*(*)
c
code ...
c
      if (length(str) .gt. 0) then
        call system (str)
      end if
c
      return
      end
