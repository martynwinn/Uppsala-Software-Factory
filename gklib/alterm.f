c
c
c
      subroutine alterm (mask,na,nb,nc,from,to)
c
c ... set all mask points which are 'FROM' to 'TO'
c
      implicit none
c
      integer from,to,i1,i2,i3,na,nb,nc
      integer mask(*)
c
code ...
c
      i1=na*nb*nc
      do i2=1,i1
        if (mask(i2) .eq. from) mask(i2) = to
      end do
c
      return
      end
