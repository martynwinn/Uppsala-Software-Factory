c
c ================================================================
c
      subroutine xps_polyline (n,xd,yd)
c
      implicit none
c
c ... draw a polygon
c
      integer n,i
c
      real xd(*),yd(*)
c
code ...
c
      if (n .lt. 2) return
c
      call xps_move (xd(1),yd(1))
      do i=1,n
        call xps_draw (xd(i),yd(i))
      end do
c
      return
      end
