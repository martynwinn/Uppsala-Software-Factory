c
c
c
      subroutine wamore (file,iunit,map,nx,ny,nz,
     +                   origin,grid,cell,ierr)
c
c ... WAMORE - write binary AMORE map
c
c ... GJ Kleywegt @ 980113
c
      implicit none
c
      integer extent(3),iunit,origin(3),grid(3),ierr,nx,ny,nz
      integer j,ix,iy,iz
c
      real map(nx,ny,nz)
      real cell(6),a(3,3),xlo(3),xhi(3),flo(3),fhi(3)
c
      logical xinter
c
      character file*(*)
c
code ...
c
      ierr = 0
      call xopxub (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) return
c
      extent(1) = nx
      extent(2) = ny
      extent(3) = nz
c
      call orthog (cell, a, 0)
c
c ... calc origin and top in fractionals
c
      do j=1,3
        flo (j) = float(origin(j)) * 
     +             (cell(j)/float(grid(j))) / cell(j)
        fhi (j) = float(extent(j)+origin(j)-1) *
     +             (cell(j)/float(grid(j))) / cell(j)
      end do
c
      call mulmtx (a, flo, xlo, 3, 3, 1)
      call mulmtx (a, fhi, xhi, 3, 3, 1)
c
      call fvalut (' XYZ min (A):',3,xlo)
      call fvalut (' XYZ max (A):',3,xhi)
c
      write (iunit,err=8000) extent(1),extent(2),extent(3),
     +  xlo(1),xlo(2),xlo(3),xhi(1),xhi(2),xhi(3),
     +  cell(4),cell(5),cell(6)
c
      do iz=1,nz
        write (iunit,err=8000) ((map(ix,iy,iz),ix=1,nx),iy=1,ny)
      end do
c
      ierr = 0
      close (iunit)
      return
c
 8000 continue
      call errcon ('While writing binary AMORE map')
      ierr = -1
      close (iunit)
c
      return
      end
