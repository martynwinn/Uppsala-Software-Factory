c
c
c
      subroutine rgb2o (r,g,b,icol)
c
c ... RGB2O - convert RGB values [0,1] into an O integer colour code
c
      implicit none
c
      real r,g,b
c
      integer icol,i,j,k
c
code ...
c
      i = max (0, min (nint (255.0 * r), 255))
      j = max (0, min (nint (255.0 * g), 255))
      k = max (0, min (nint (255.0 * b), 255))
c
      icol = i*256*256 + j*256 + k
c
      return
      end
