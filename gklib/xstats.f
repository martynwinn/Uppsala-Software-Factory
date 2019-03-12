c
c ========================================================================
c
      subroutine xstats (xdata,ndata,save,ssdv,sxmin,sxmax,sxtot)
c
c ... XSTATS - calc average etc. of an array
c
c ... Gerard Kleywegt @ 920715
c
      implicit none
c
      real xdata(*),save,ssdv,sxmin,sxmax,sxtot,qqq
      double precision ave,sdv,xnd,sqr,xmin,xmax,xtot
c
      integer i,ndata
c
code ...
c
      if (ndata .lt. 1) then
        call errcon ('XSTATS - Not enough data points')
        call jvalut (' Nr of data points :',1,ndata)
        return
      end if
c
      save  = xdata(1)
      ssdv  = 0.0
      sxmin = xdata(1)
      sxmax = xdata(1)
      sxtot = xdata(1)
c
      if (ndata .eq. 1) return
c
      ave = xdata(1)
      sdv = xdata(1)*xdata(1)
      xmin = xdata(1)
      xmax = xdata(1)
      xtot = xdata(1)
      do i=2,ndata
        ave = ave + xdata(i)
        sdv = sdv + xdata(i)*xdata(i)
        if (xdata(i) .lt. xmin) xmin = xdata(i)
        if (xdata(i) .gt. xmax) xmax = xdata(i)
ccc        xmin = min (xmin,xdata(i))
ccc        xmax = max (xmax,xdata(i))
      end do
      xtot = ave
c
      xnd = float (ndata)
      ave = ave / xnd
c
ccc      print *,' N, tot, ave ',ndata,xtot,ave,xtot/float(ndata)
ccc      print *,' min max ',xmin,xmax
c
      if (ndata .eq. 2) then
        sdv = dabs (0.5D0*(xdata(2)-xdata(1)))
      else
        sqr = (sdv/xnd) - (ave*ave)
        if (sqr .lt. 0.0D0) then
          sqr = 0.0D0
          do i=1,ndata
            sqr = sqr + (xdata(i)-ave)**2
          end do
          sqr = sqr / xnd
          if (sqr .lt. 0.0D0) then
            call errcon ('XSTATS - Root of negative number')
            qqq = sqr
            call rvalut (' Value :',1,qqq)
            sdv = 0.0D0
          else
            sdv = dsqrt (sqr)
          end if
        else
          sdv = dsqrt (sqr)
        end if
      end if
c
      save = ave
      ssdv = sdv
      sxmin = xmin
      sxmax = xmax
      sxtot = xtot
c
      return
      end
