c
c
c
      subroutine voxvol (cell,grid,cvol,vxv)
c
c ... VOXVOL - calculate voxel volume
c
      implicit NONE
c
      real cell(6),cvol,vxv
c
      integer grid(3),i
c
code ...
c
      call celvol (cell,cvol)
c
      vxv = cvol
      do i=1,3
        vxv = vxv / float(grid(i))
      end do
c
      return
      end
