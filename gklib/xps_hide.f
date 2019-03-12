c
c ================================================================
c
      subroutine xps_hide (log)
c
c ... switch white-out of areas outside plot on or off
c
      include 'xps.incl'
c
      logical log
c
code ...
c
      pshide = log
c
      return
      end
