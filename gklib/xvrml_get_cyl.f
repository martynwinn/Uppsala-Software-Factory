c
c
c
      subroutine xvrml_get_cyl (x1,x2,a,b,c,d,e,xc)
c
c ... figure out parameters for a cylinder ("stick")
c
      implicit none
c
      real x1(3),x2(3),xc(3),a,b,c,d,e
      real y2(3)
c
      integer i
c
code ...
c
c ... XC = centroid of X1 and X2
c
      do i=1,3
        xc (i) = 0.5 * (x2(i) + x1(i))
        y2 (i) = x2(i) - x1(i)
      end do
c
c ... D is length of vector X1->X2
c
      d = sqrt (y2(1)**2 + y2(2)**2 + y2(3)**2)
      do i=1,3
        y2(i) = y2(i) / d
      end do
      c = sqrt(y2(1)**2 + y2(3)**2)
c
c ... E = rotation angle
c
      e = - acos ( y2(2) )
c
c ... A,B,C = direction of rotation axis (undetermined if rotation = 0)
c
      if (abs(e) .le. 0.01) then
        a = 1.0
        b = 0.0
        c = 1.0
      else
        a = -y2(3)/c
        b = 0.0
        c = y2(1)/c
      end if
c
      return
      end
