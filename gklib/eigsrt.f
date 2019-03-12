c
c ========================================================================
c
      subroutine eigsrt (d,v,n,np)
c
c ... simple sort of eigenvectors by decreasing eigenvalue
c
c ... Numerical Recipes, pp. 348
c
      implicit none
c
      integer n,np
c
      real d(np),v(np,np)
c
      real p
c
      integer i,j,k
c
code ...
c
      do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
          if (d(j) .ge. p) then
            k=j
            p=d(j)
          end if
        end do
        if (k.ne.i) then
          d(k)=d(i)
          d(i)=p
          do j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
          end do
        end if
      end do
c
      return
      end
