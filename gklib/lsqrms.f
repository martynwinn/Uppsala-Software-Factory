c
c
c
      subroutine lsqrms (x,y,n,r)
c
c ... LSQRMS = calculate RMSD for N equivalenced atoms
c
      implicit none
c
      integer n,i
c
      real x(3,n),y(3,n),r
c
code ...
c
      r = 0.0
      if (n .le. 1) return
c
      do i=1,n
        r = r + (x(1,i)-y(1,i))**2 + (x(2,i)-y(2,i))**2 +
     +          (x(3,i)-y(3,i))**2
      end do
c
      r = sqrt (r/float(n))
c
      return
      end
