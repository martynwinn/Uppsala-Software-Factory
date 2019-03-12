c
c
c
      subroutine inimap (map,na,nb,nc,ns)
c
c ... initialise a map with a value of NS for all points
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      real map(*),ns
c
code ...
c
c      do i1=1,na
c        do i2=1,nb
c          do i3=1,nc
c            map (i1,i2,i3) = ns
c          end do
c        end do
c      end do
c
      i1 = na*nb*nc
      do i2=1,i1
        map(i2) = ns
      end do
c
      return
      end
