c
c ===========================================================================
c
      subroutine fratra (tra)
c
c ... set fractional translations to machine precision
c     implemented: 1/3, 2/3, 1/6, 5/6
c
      implicit none
c
      double precision x1o3,x2o3,x1o6,x5o6
      parameter (x1o3=1.000000D0/3.000000D0)
      parameter (x2o3=2.000000D0/3.000000D0)
      parameter (x1o6=1.000000D0/6.000000D0)
      parameter (x5o6=5.000000D0/6.000000D0)
c
      real tra(3)
c
      integer i
c
code ...
c
      do i=1,3
c
        if (int(10000.0*tra(i)) .eq. 3333) then
          tra(i) = x1o3
        else if (nint(10000.0*tra(i)) .eq. 6667) then
          tra(i) = x2o3
        else if (nint(10000.0*tra(i)) .eq. 1667) then
          tra(i) = x1o6
        else if (nint(10000.0*tra(i)) .eq. 8333) then
          tra(i) = x5o6
        end if
c
      end do
c
      return
      end
