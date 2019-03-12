c
c
c
      subroutine nn8 (a, sizx, sizy, sizz, x, value, errcod)
c
c --- "ZERO-ORDER" interpolation: take nearest grid point
c
      implicit none
c
      integer errcod, sizx, sizy, sizz
c
      real a(sizx, sizy, sizz), x(3), value
c
      integer i,j,k
c
code ...
c
      errcod = 1
c
      i = nint(x(1))
      j = nint(x(2))
      k = nint(x(3))
c
      if (i .le. 0) goto 10
      if (j .le. 0) goto 10
      if (k .le. 0) goto 10
      if (i .gt. sizx) goto 10
      if (j .gt. sizy) goto 10
      if (k .gt. sizz) goto 10
c
      value = a(i,j,k)
      errcod = 0
c
   10 continue
c
      return
      end
