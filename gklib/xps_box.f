c
c ================================================================
c
      subroutine xps_box (x1,x2,y1,y2)
c
      implicit none
c
c ... draw a box
c
      real x1,x2,y1,y2
c
code ...
c
      call xps_move (x1,y1)
      call xps_draw (x2,y1)
      call xps_draw (x2,y2)
      call xps_draw (x1,y2)
      call xps_draw (x1,y1)
c
      return
      end
