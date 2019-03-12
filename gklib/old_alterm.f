c
c
c
      subroutine old_alterm (mask,na,nb,nc,from,to)
c
c ... set all mask points which are 'FROM' to 'TO'
c
      implicit none
c
      integer from,to,i1,i2,i3,na,nb,nc
      integer mask(na,nb,nc)
c
code ...
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. from) mask(i1,i2,i3) = to
          end do
        end do
      end do
c
      return
      end
