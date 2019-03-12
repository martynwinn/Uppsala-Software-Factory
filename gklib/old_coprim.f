c
c
c
      subroutine old_coprim (rmask,imask,na,nb,nc,mode)
c
c ... copy REAL to INTEGER mask/map or vice versa
c
      implicit none
c
      integer na,nb,nc,mode,i,j,k
      real rmask(na,nb,nc)
      integer imask(na,nb,nc)
c
code ...
c
      if (mode .ge. 0) then
c
c ... copy from REAL to INTEGER
c
        do i=1,na
          do j=1,nb
            do k=1,nc
              imask (i,j,k) = nint (rmask(i,j,k))
            end do
          end do
        end do
      else
c
c ... copy from INTEGER to REAL
c
        do i=1,na
          do j=1,nb
            do k=1,nc
              rmask (i,j,k) = float (imask(i,j,k))
            end do
          end do
        end do
      end if
c
      return
      end
