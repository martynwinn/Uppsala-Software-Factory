c
c ==========================================================================
c
      subroutine vecang (vec1,vec2,alpha,ierr)
c
c ... calculate the angle between two vectors
c
      implicit none
c
      real twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      real vec1(3),vec2(3),alpha,ip,l1,l2,cosa
c
      integer ierr,i
c
code ...
c
      ierr = -1
      ip = 0.0
      l1 = 0.0
      l2 = 0.0
c
      do i=1,3
        ip = ip + vec1(i)*vec2(i)
        l1 = l1 + vec1(i)*vec1(i)
        l2 = l2 + vec2(i)*vec2(i)
      end do
c
      if (l1 .le. 0.0 .or. l2 .le. 0.0) then
        alpha = -9999.999
        ierr = -1
        return
      end if
c
      cosa = ip / sqrt(l1*l2)
      cosa = max (-1.0, min (1.0, cosa))
      alpha = rtodeg * acos ( cosa )
      call fixang (alpha)
      ierr = 0
c
      return
      end
