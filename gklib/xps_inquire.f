c
c ================================================================
c
      subroutine xps_inquire (rxmin,rxmax,rymin,rymax)
c
c ... return actual PostScript coordinate limits for plot area
c
      include 'xps.incl'
c
      real rxmin,rxmax,rymin,rymax
c
code ...
c
      rxmin = psxmin
      rymin = psymin
      rxmax = psxmax
      rymax = psymax
c
      return
      end
