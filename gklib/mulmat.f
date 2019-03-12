c
c ===========================================================================
c
      subroutine mulmat (a,b,c)
c
c ... MULMAT - multiply two 3*3 matrices
c
      real a(3,*),b(3,*),c(3,*)
c
      integer i,j
c
code ...
c
      do i=1,3
        do j=1,3
          c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j) + a(i,3)*b(3,j)
        end do
      end do
c
      return
      end
