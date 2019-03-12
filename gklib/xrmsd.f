c
c ========================================================================
c
      subroutine xrmsd (xdata,ydata,ndata,srmsd)
c
c ... XRMSD - calc Root-Mean-Square Distance between two arrays
c
c ... Gerard Kleywegt @ 920922
c
      implicit none
c
      real xdata(*),ydata(*),srmsd
      double precision sumsqd,rmsd
c
      integer i,ndata
c
code ...
c
      sumsqd = 0.0D0
      rmsd = 0.0D0
      srmsd = 0.0
c
      if (ndata .lt. 1) then
        call errcon ('XRMSD - Not enough data points')
        call jvalut (' Nr of data points :',1,ndata)
        return
      end if
c
      do i=1,ndata
        sumsqd = sumsqd + (xdata(i)-ydata(i))**2
      end do
c
      rmsd = dsqrt (sumsqd / float (ndata))
      srmsd = rmsd
c
      return
      end
