c
c ================================================================
c
      subroutine xps_ellipse_part (xc,yc,xr,yr,i1,i2,i3)
c
c ... draw a partial ellipse
c
      include 'xps.incl'
c
      real xc,yc,xr,yr,phi
c
      integer i,i1,i2,i3
c
code ...
c
      phi = degtor * float(i1)
      call xps_move (xc+xr*cos(phi),yc+yr*sin(phi))
      do i=i1,i2,i3
        phi = degtor * float(i)
        call xps_draw (xc+xr*cos(phi),yc+yr*sin(phi))
      end do
c
      return
      end
