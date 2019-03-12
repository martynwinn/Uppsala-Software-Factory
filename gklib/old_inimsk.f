c
c
c
      subroutine old_inimsk (mask,na,nb,nc,ns)
c
c ... initialise a mask with a value of NS for all points
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,ns
      integer mask(na,nb,nc)
c
code ...
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            mask (i1,i2,i3) = ns
          end do
        end do
      end do
c
      return
      end
