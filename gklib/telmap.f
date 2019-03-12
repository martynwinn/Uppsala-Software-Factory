c
c=====================================================
c
      subroutine telmap (grid,orgn,ext,cell)
c
      implicit none
c
      real cell(6)
c
      integer grid(3),orgn(3),ext(3)
c
code ...
c
      call fvalut (' Cell axes   (A) :',3,cell(1))
      call fvalut (' Cell angles (d) :',3,cell(4))
      call ivalut (' Grid axes (pts) :',3,grid)
      call ivalut (' Origin    (pts) :',3,orgn)
      call ivalut (' Extent    (pts) :',3,ext)
c
      return
      end
