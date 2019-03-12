c
c
c
      subroutine copmij (mask1,mask2,na,nb,nc)
c
c ... copy mask2 to mask1
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      integer mask1(*),mask2(*)
c
code ...
c
      i1 = na*nb*nc
      do i2=1,i1
        mask1 (i2) = mask2 (i2)
      end do
c
      return
      end
