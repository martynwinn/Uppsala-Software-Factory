c
c ================================================================
c
      subroutine xps_ellipse (xc,yc,xr,yr)
c
c ... draw an ellipse
c
      include 'xps.incl'
c
      real xc,yc,xr,yr,phi
c
      integer i
c
code ...
c
      call xps_move (xc+xr,yc)
      do i=0,360,1
        phi = degtor * float(i)
        call xps_draw (xc+xr*cos(phi),yc+yr*sin(phi))
      end do
c
      return
      end
