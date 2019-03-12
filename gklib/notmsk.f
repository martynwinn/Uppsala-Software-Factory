c
c
c
      subroutine notmsk (mask,na,nb,nc)
c
c ... logical NOT of a mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      integer mask(*)
c
code ...
c
      i1 = na*nb*nc
      do i2=1,i1
        mask(i2) = 1 - mask(i2)
      end do
c
      return
      end
