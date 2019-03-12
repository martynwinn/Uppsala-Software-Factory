c
c
c
      real function dist (i1,i2,xyz)
c
c ... calculate distance between atoms I1 and I2; XYZ(3,*) holds coords
c
      implicit none
c
      real xyz(3,*)
c
      integer i1,i2,i
c
code ...
c
      dist = 0.0
c
      do i=1,3
        dist = dist + (xyz(i,i1) - xyz(i,i2))**2
      end do
c
      dist = sqrt (dist)
c
      return
      end
