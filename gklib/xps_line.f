c
c ================================================================
c
      subroutine xps_line (x1,y1,x2,y2)
c
      implicit none
c
c ... draw a line
c
      real x1,x2,y1,y2
c
code ...
c
      call xps_move (x1,y1)
      call xps_draw (x2,y2)
c
      return
      end
