c
c
c
      subroutine cutmsk (mask,shadow,na,nb,nc,fac)
c
c ... slice bits off mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,fac
      integer mask(*),shadow(*)
c
code ...
c
      call cntnbr (mask,shadow,na,nb,nc,1,0)
c
      i1 = na*nb*nc
      do i2=1,i1
        if (shadow(i2) .ge. fac) mask(i2) = 0
      end do
c
      return
      end
