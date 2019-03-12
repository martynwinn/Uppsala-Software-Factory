c
c
c
      subroutine smooth (mask,shadow,na,nb,nc,fac)
c
c ... smooth mask 
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,fac
      integer mask(*),shadow(*)
c
code ...
c
      call cntnbr (mask,shadow,na,nb,nc,0,1)
c
      i1 = na*nb*nc
      do i2=1,i1
        if (shadow(i2) .ge. fac) mask(i2) = 1
      end do
c
      return
      end
