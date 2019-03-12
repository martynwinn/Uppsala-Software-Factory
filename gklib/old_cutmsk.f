c
c
c
      subroutine old_cutmsk (mask,shadow,na,nb,nc,fac)
c
c ... slice bits off mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,fac
      integer mask(na,nb,nc),shadow(na,nb,nc)
c
code ...
c
      call cntnbr (mask,shadow,na,nb,nc,1,0)
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            if (shadow(i1,i2,i3) .ge. fac) mask(i1,i2,i3) = 0
          end do
        end do
      end do
c
      return
      end
