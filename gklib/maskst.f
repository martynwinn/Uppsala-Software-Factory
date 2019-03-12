c
c
c
      subroutine maskst (mask,ext1,ext2,ext3)
c
c ---	Get statistics on the mask.
c ---	This version just says how many points are in it
c
      implicit none
c
      integer ext1, ext2, ext3
      integer mask (ext1, ext2, ext3)
      integer i, j, k, ct
c
code ...
c
      ct = 0
      do k=1,ext3
        do j=1,ext2
          do i=1,ext1
            if (mask(i,j,k) .eq. 1) ct = ct+1
          end do
        end do
      end do
c
      call jvalut (' Number of points in mask :',1,ct)
c
      return
      end
