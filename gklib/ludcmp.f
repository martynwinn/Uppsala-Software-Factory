c
c ========================================================================
c
      subroutine ludcmp (a,n,np,indx,d,vv,ierr)
c
c ... LU-decomposition of a N*N matrix A, physical dim NP*NP
c     INDX integer *N, D real scalar, VV real *N
c     IERR = 0 means success
c
c ... Numerical Recipes, section 2.3, pp. 35
c
      implicit none
c
      real tiny
      parameter (tiny=1.0e-20)
c
      integer n,np
c
      real a(np,np),vv(n)
      real d,sum,aamax,dum
c
      integer indx(n)
      integer i,j,k,ierr,imax
c
code ...
c
      d = 1.0
      ierr = -1
c
      do i=1,n
        aamax = 0.0
        do j=1,n
          if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
        end do
        if (aamax .le. 0) then
          call errcon ('Singular matrix in LUDCMP - ignore results')
          return
        end if
        vv(i) = 1.0/aamax
      end do
c
      do j=1,n
        if (j .gt. 1) then
          do i=1,j-1
            sum = a(i,j)
            if (i .gt. 1) then
              do k=1,i-1
                sum = sum - a(i,k)*a(k,j)
              end do
              a(i,j) = sum
            end if
          end do
        end if
c
        aamax = 0.0
        do i=j,n
          sum = a(i,j)
          if (j .gt. 1) then
            do k=1,j-1
              sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j) = sum
          end if
          dum = vv(i)*abs(sum)
          if (dum .ge. aamax) then
            imax = i
            aamax = dum
          end if
        end do
c
        if (j .ne. imax) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          end do
          d = - d
          vv(imax) = vv(j)
        end if
c
        indx(j) = imax
        if (j .ne. n) then
          if (a(j,j) .eq. 0.0) a(j,j) = tiny
          dum = 1.0/a(j,j)
          do i=j+1,n
            a(i,j) = a(i,j)*dum
          end do
        end if
      end do
c
      if (a(n,n) .eq. 0.0) a(n,n) = tiny
      ierr = 0
c
      return
      end
