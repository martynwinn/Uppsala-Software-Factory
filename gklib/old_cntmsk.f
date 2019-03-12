c
c
c
      subroutine old_cntmsk (mask,na,nb,nc,ns)
c
c ... count mask points
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,ns
      integer mask(na,nb,nc)
c
code ...
c
      ns = 0
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. 1) ns = ns + 1
          end do
        end do
      end do
c
      return
      end
