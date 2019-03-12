c
c ================================================================
c
      subroutine xps_solid ()
c
c ... reset solid line type
c
      include 'xps.incl'
c
code ...
c
      psdash = 0
      call xps_dash ()
c
      return
      end
