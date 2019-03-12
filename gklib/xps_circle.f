c
c ================================================================
c
      subroutine xps_circle (cx,cy,crad)
c
c ... draw circle centred at (CX,CY) with radius CRAD
c
      include 'xps.incl'
c
      real cx,cy,crad,phi
c
      integer i
c
code ...
c
      call xps_move (cx+crad,cy)
      do i=0,360,1
        phi = degtor * float(i)
        call xps_draw (cx+crad*cos(phi),cy+crad*sin(phi))
      end do
c
      return
      end
