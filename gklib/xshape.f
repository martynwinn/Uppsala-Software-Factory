c
c ========================================================================
c
      subroutine xshape (xdata,ydata,ndata,match)
c
c ... XSHAPE - calculate how well the SHAPES of two arrays of
c              values match
c
c     match = cos(alpha) = cos (angle (xdata, ydata)) =
c           = dot-product(xdata,ydata) / SQRT(|xdata|)*SQRT(|ydata|)
c     where - dot-product(xdata,ydata) = SUM (xdata(i)*ydata(i))
c             |xdata| = SUM (xdata(i)**2)
c             |ydata| = SUM (ydata(i)**2)
c
c     match will always be in [-1.0,+1.0]
c     match = -1.0 means same shape but opposite signs
c     match = +1.0 means same shape
c
c     note that the value of MATCH is INDEPENDENT of the
c     MAGNITUDES of XDATA compared to those of YDATA !!!
c
c     note: if both averages are ZERO, then MATCH is equal to
c           the correlation coefficient between XDATA and YDATA
c           (I think)
c
c ... Gerard Kleywegt @ 920923
c
      implicit none
c
      real xdata(*),ydata(*),match
      double precision sumxx,sumyy,sumxy,xmatch,dzero
c
      integer i,ndata
c
code ...
c
      dzero = 0.0D0
c
      sumxx = dzero
      sumyy = dzero
      sumxy = dzero
      match = -1.0
c
      if (ndata .lt. 1) then
        call errcon ('XSHAPE - Not enough data points')
        call jvalut (' Nr of data points :',1,ndata)
        return
      end if
c
      do i=1,ndata
        sumxx = sumxx + (xdata(i)*xdata(i))
        sumxy = sumxy + (xdata(i)*ydata(i))
        sumyy = sumyy + (ydata(i)*ydata(i))
      end do
c
      if (sumxx .eq. dzero .and. sumyy .eq. dzero) then
        match = 1.0
      else if (sumxx .eq. dzero .or. sumyy .eq. dzero) then
        match = -1.0
      else
        xmatch = sumxy / (dsqrt(sumxx) * dsqrt(sumyy))
        match = xmatch
      end if
c
      return
      end
