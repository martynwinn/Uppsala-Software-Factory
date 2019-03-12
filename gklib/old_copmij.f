c
c
c
      subroutine old_copmij (mask1,mask2,na,nb,nc)
c
c ... copy mask2 to mask1
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      integer mask1(na,nb,nc),mask2(na,nb,nc)
c
code ...
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            mask1 (i1,i2,i3) = mask2 (i1,i2,i3)
          end do
        end do
      end do
c
      return
      end
