c
c ========================================================================
c
      subroutine xstat2 (xdata,ndata,srms,shav)
c
c ... XSTAT2 - calc RMS and harmonic average of an array
c
c ... Gerard Kleywegt @ 990823/010906
c
      implicit none
c
      real xdata(*),srms,shav
      double precision ave,sdv
c
      integer i,ndata,notnil
c
code ...
c
      srms  = xdata(1)
      shav  = xdata(1)
c
      if (ndata .le. 1) then
        call errcon ('XSTAT2 - Not enough data points')
        call jvalut (' Nr of data points :',1,ndata)
        return
      end if
c
      ave = 0.0D0
      sdv = 0.0D0
      notnil = 0
      do i=1,ndata
        ave = ave + xdata(i)*xdata(i)
        if (xdata(i) .ne. 0.0) then
          notnil = notnil + 1
          sdv = sdv + 1.0D0/xdata(i)
        end if
      end do
c
      ave = dsqrt (ave / float (ndata))
      if (notnil .gt. 0) then
        sdv = 1.0D0 / (sdv / float (notnil))
        if (notnil .ne. ndata) then
          call errcon ('XSTAT2 - Zeroes ignored for harmonic average')
          call jvalut (' Nr of zeroes :',1,(ndata-notnil))
        end if
      else
        call errcon (
     +    'XSTAT2 - Array of zeros - harmonic average undefined')
      end if
c
      srms = ave
      shav = sdv
c
      return
      end
