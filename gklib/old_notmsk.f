c
c
c
      subroutine old_notmsk (mask,na,nb,nc)
c
c ... logical NOT of a mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      integer mask(na,nb,nc)
c
code ...
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            mask(i1,i2,i3) = 1 - mask(i1,i2,i3)
          end do
        end do
      end do
c
      return
      end
