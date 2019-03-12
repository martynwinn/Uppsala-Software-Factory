c
c ===========================================================================
c
      subroutine spawno (command)
c
      implicit none
c
c === SPAWNO (COMMAND) => pass string COMMAND to the operating system
c
c === G J Kleywegt @ 920519
c
      character command*(*)
c
code ...
c
      call gksys (command)
c
      return
      end
