c
c
c
      subroutine coprim (rmask,imask,na,nb,nc,mode)
c
c ... copy REAL to INTEGER mask/map or vice versa
c
      implicit none
c
      integer na,nb,nc,mode,i,j,k
      real rmask(*)
      integer imask(*)
c
code ...
c
      if (mode .ge. 0) then
c
c ... copy from REAL to INTEGER
c
        i = na*nb*nc
        do j=1,i
          imask (j) = nint (rmask(j))
        end do
      else
c
c ... copy from INTEGER to REAL
c
        i = na*nb*nc
        do j=1,i
          rmask (j) = float (imask(j))
        end do
      end if
c
      return
      end
