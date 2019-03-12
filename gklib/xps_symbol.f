c
c ================================================================
c
      subroutine xps_symbol (isym,x1,x2,y1,y2)
c
c ... plot a symbol
c
c ... 0 = default = box
c     1 = plus
c     2 = cross
c     3 = Z
c     4 = diamond
c     5 = triangle point up
c     6 = triangle point down
c     7 = triangle point left
c     8 = triangle point right
c
      include 'xps.incl'
c
      real x1,x2,y1,y2,xm,ym
c
      integer isym
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      xm = 0.5*(x1+x2)
      ym = 0.5*(y1+y2)
c
c ... 1 = PLUS
c
      if (isym .eq. 1) then
        call xps_move (x1,ym)
        call xps_draw (x2,ym)
        call xps_move (xm,y1)
        call xps_draw (xm,y2)
        call xps_stroke ()
c
c ... 2 = CROSS
c
      else if (isym .eq. 2) then
        call xps_move (x1,y1)
        call xps_draw (x2,y2)
        call xps_move (x2,y1)
        call xps_draw (x1,y2)
        call xps_stroke ()
c
c ... 3 = Z
c
      else if (isym .eq. 3) then
        call xps_move (x1,y2)
        call xps_draw (x2,y2)
        call xps_draw (x1,y1)
        call xps_draw (x2,y1)
        call xps_stroke ()
c
c ... 4 = diamond
c
      else if (isym .eq. 4) then
        call xps_move (x1,ym)
        call xps_draw (xm,y2)
        call xps_draw (x2,ym)
        call xps_draw (xm,y1)
        call xps_draw (x1,ym)
        call xps_stroke ()
c
c ... 5 = triangle, point up
c
      else if (isym .eq. 5) then
        call xps_move (x1,y1)
        call xps_draw (xm,y2)
        call xps_draw (x2,y1)
        call xps_draw (x1,y1)
        call xps_stroke ()
c
c ... 6 = triangle, point down
c
      else if (isym .eq. 6) then
        call xps_move (x1,y2)
        call xps_draw (x2,y2)
        call xps_draw (xm,y1)
        call xps_draw (x1,y2)
        call xps_stroke ()
c
c ... 7 = triangle, point left
c
      else if (isym .eq. 7) then
        call xps_move (x2,y1)
        call xps_draw (x1,ym)
        call xps_draw (x2,y2)
        call xps_draw (x2,y1)
        call xps_stroke ()
c
c ... 8 = triangle, point right
c
      else if (isym .eq. 8) then
        call xps_move (x1,y2)
        call xps_draw (x2,ym)
        call xps_draw (x1,y1)
        call xps_draw (x1,y2)
        call xps_stroke ()
c
c ... room for more
c
      else
c
c ... 0 = DEFAULT = BOX
c
        call xps_move (x1,y1)
        call xps_draw (x1,y2)
        call xps_draw (x2,y2)
        call xps_draw (x2,y1)
        call xps_draw (x1,y1)
        call xps_stroke ()
      end if
c
      return
      end
