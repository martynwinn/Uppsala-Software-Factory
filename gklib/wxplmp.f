c
c
c
      subroutine wxplmp (file,iunit,map,nx,ny,nz,
     +                   origin,grid,cell,ave,sdv,ierr)
c
c ... WXPLMAP - write ASCII X-PLOR (4) map
c
      implicit none
c
      integer nx,ny,nz,iunit,origin(3),grid(3),ierr
      integer length,i,j,k
c
      real map(nx,ny,nz)
      real cell(6),ave,sdv
c
      logical xinter
c
      character file*(*),line*100
c
code ...
c
      ierr = 0
      call xopxua (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
c
      write (iunit,'(/i8)',err=8000) 2
      write (iunit,'(3a)',err=8000) 'REMARKS FILENAME="',
     +  file(1:length(file)),'"'
      write (iunit,'(9a)',err=8000) 'REMARKS ',
     +  line(1:length(line))
c
      write (iunit,'(9i8)',err=8000)
     +  grid(1),origin(1),origin(1)+nx-1,
     +  grid(2),origin(2),origin(2)+ny-1,
     +  grid(3),origin(3),origin(3)+nz-1
c
 6000 format (1p,6e12.5)
      write (iunit,6000,err=8000) (cell(i),i=1,6)
c
      write (iunit,'(a3)',err=8000) 'ZYX'
c
      do k=1,nz
        write (iunit,'(i8)',err=8000) k-1
        write (iunit,6000,err=8000)
     +    ((map(i,j,k),i=1,nx),j=1,ny)
      end do
c
      write (iunit,'(i8)',err=8000) -9999
      write (iunit,'(1x,1p,6e12.5)',err=8000) ave,sdv
c
      ierr = 0
      close (iunit)
      return
c
 8000 continue
      call errcon ('While writing ASCII X-PLOR map')
      ierr = -1
      close (iunit)
c
      return
      end
