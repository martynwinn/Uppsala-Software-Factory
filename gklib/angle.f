c
c
c
      real function angle (i1,i2,i3,xyz)
c
c ... calculate angle of atoms I1-I2-I3; XYZ(3,*) holds coords
c
      implicit none
c
      real rtodeg
      parameter (rtodeg = 360.0/6.2831853071796)
c
      real xyz(3,*),v1(3),v2(3),xx,yy,xy,qqq
c
      integer i1,i2,i3,j
c
code ...
c
c --- form the two interatomic vectors
c
      do j=1,3
        v1(j) = xyz(j,i2) - xyz(j,i1)
        v2(j) = xyz(j,i2) - xyz(j,i3)
      end do
c
c --- cos (bond_angle) = (V1.dot.V2) / (|V1|*|V2|)
c
      xx = 0.0
      yy = 0.0
      xy = 0.0
      do j=1,3
        xx = xx + v1(j)**2
        yy = yy + v2(j)**2
        xy = xy + v1(j)*v2(j)
      end do
c
      if (abs(xy) .le. 1.0E-10) then
        angle = 0.0
      else
        qqq = max (-1.0, min (1.0, xy/(sqrt(xx)*sqrt(yy))))
        angle = acos( qqq ) * rtodeg
      end if
c
      return
      end
