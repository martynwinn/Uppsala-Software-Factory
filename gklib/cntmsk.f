c
c
c
      subroutine cntmsk (mask,na,nb,nc,ns)
c
c ... count mask points
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,ns
      integer mask(*)
c
code ...
c
      ns = 0
c
      i1 = na*nb*nc
      do i2=1,i1
        if (mask(i2) .eq. 1) ns = ns + 1
      end do
c
      return
      end
