c
c ... numeric.f - numeric subroutines
c
c ... gerard j kleywegt @ 960409
c
c ========================================================================
c
      subroutine avesdv (xdata,ndata,save,ssdv)
c
c ... AVESDV - calc average and standard deviation of an array
c
c ... Gerard Kleywegt @ 920708
c
      implicit none
c
      real xdata(*),save,ssdv,qqq
      double precision ave,sdv,xnd,sqr
c
      integer i,ndata
c
code ...
c
      save = 0.0
      ssdv = 0.0
c
      if (ndata .le. 1) then
        call errcon ('AVESDV - Not enough data points')
        call jvalut (' Nr of data points :',1,ndata)
        return
      end if
c
      ave = 0.0D0
      sdv = 0.0D0
      do i=1,ndata
        ave = ave + xdata(i)
        sdv = sdv + xdata(i)*xdata(i)
      end do
c
      xnd = float (ndata)
      ave = ave / xnd
      if (ndata .eq. 2) then
        sdv = abs (0.5D0*(xdata(2)-xdata(1)))
      else
        sqr = (sdv/xnd) - (ave*ave)
        if (sqr .lt. 0.0D0) then
          call errcon ('AVESDV - Root of negative number')
          qqq = sqr
          call rvalut (' Value :',1,qqq)
          sdv = 0.0D0
        else
          sdv = dsqrt (sqr)
        end if
      end if
c
      save = ave
      ssdv = sdv
c
      return
      end
