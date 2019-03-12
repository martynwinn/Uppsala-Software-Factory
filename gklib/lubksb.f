c
c ========================================================================
c
      subroutine lubksb (a,n,np,indx,b)
c
c ... solve Ax=b with N*N matrix A, physical dim NP*NP
c     B input = RHS vector b; B output = solution vector x
c     INDX integer *N
c
c ... Numerical Recipes, section 2.3, pp. 36
c
      implicit none
c
      integer n,np
c
      real a(np,np),b(n)
      real sum
c
      integer indx(n)
      integer i,j,ii,ll
c
code ...
c
      ii = 0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii .ne. 0) then
          do j=ii,i-1
            sum = sum - a(i,j)*b(j)
          end do
        else if (sum .ne. 0.0) then
          ii = i
        end if
        b(i) = sum
      end do
c
      do i=n,1,-1
        sum = b(i)
        if (i .lt. n) then
          do j=i+1,n
            sum = sum - a(i,j)*b(j)
          end do
        end if
        b(i) = sum/a(i,i)
      end do
c
      return
      end
