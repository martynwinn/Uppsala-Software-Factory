c
c
c
      subroutine inimsk (mask,na,nb,nc,ns)
c
c ... initialise a mask with a value of NS for all points
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,ns
      integer mask(*)
c
code ...
c
      i1 = na*nb*nc
      do i2=1,i1
        mask (i2) = ns
      end do
c
      return
      end
