c
c ================================================================
c
      subroutine xps_circle_part (cx,cy,crad,i1,i2,i3)
c
c ... draw partial circle centred at (CX,CY) with radius CRAD
c
      include 'xps.incl'
c
      real cx,cy,crad,phi
c
      integer i,i1,i2,i3
c
code ...
c
      phi = degtor * float(i1)
      call xps_move (cx+crad*cos(phi),cy+crad*sin(phi))
      do i=i1,i2,i3
        phi = degtor * float(i)
        call xps_draw (cx+crad*cos(phi),cy+crad*sin(phi))
      end do
c
      return
      end
