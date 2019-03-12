c
c ================================================================
c
      subroutine xps_polar (xmin,xmax)
c
c ... set up for use of polar coordinates
c
      include 'xps.incl'
c
      real xmin,xmax,dx
c
      integer i
c
code ...
c
      pspola = .true.
      poxmin = min (xmin,xmax)
      poxmax = max (xmin,xmax)
      dx = poxmax-poxmin
      call xps_scale  (-dx,dx,-dx,dx)
c
      call xps_move (poxmax,0.0)
      do i=0,360
        call xps_draw (poxmax,float(i))
      end do
c
      do i=0,360,30
        call xps_move (poxmax,float(i))
        call xps_draw (poxmin+0.98*dx,float(i))
      end do
c
      return
      end
