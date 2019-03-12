c
c ===========================================================================
c
      subroutine invncs (numncs,rtncs,rtinv)
c
      implicit none
c
      integer numncs,i,j,ierr
c
      real rtncs(12,numncs),rtinv(12,numncs)
      real x(3),x1(3),x2(3)
c
code ...
c
      ierr = 0
c
      do i=1,numncs
        do j=1,9
          rtinv(j,i) = rtncs(j,i)
        end do
        call matinv (rtinv(1,i), 3, x, x1, x2)
        call mulmtx (rtinv(1,i), rtncs(10,i), rtinv(10,i), 3, 3, 1)
        do j=10,12
          rtinv(j,i) = -rtinv(j,i)
        end do
        write (*,6000) i
        call anancs (1,rtncs(1,i),.true.,ierr)
        write (*,6010) i
        call anancs (1,rtinv(1,i),.true.,ierr)
      end do
c
 6000 format (/' ***** Operator nr ',i3,' *****'/)
 6010 format (/' ***** INVERSE Operator nr ',i3,' *****'/)
c
      return
      end
