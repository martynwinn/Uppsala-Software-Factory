c
c
c
      subroutine copmap (map1,map2,na,nb,nc)
c
c ... copy map2 to map1
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      real map1(*),map2(*)
c
code ...
c
c      do i1=1,na
c        do i2=1,nb
c          do i3=1,nc
c            map1 (i1,i2,i3) = map2 (i1,i2,i3)
c          end do
c        end do
c      end do
c
      i1 = na*nb*nc
      do i2=1,i1
        map1(i2) = map2(i2)
      end do
c
      return
      end
