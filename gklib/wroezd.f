c
c
c
      subroutine wroezd (map,exta,extb,extc,origin,grid,cell,
     +                   iunit,scale,form,ierr)
c
c ... write OLDEZD file
c
      implicit none
c
      real cell(6),scale
c
      integer exta,extb,extc,origin(3),grid(3),iunit
      integer i,j,k,ierr
      real map (exta,extb,extc)
c
      character form*(*)
c
code ...
c
      ierr = -1
c
      write (iunit,6000,err=900) (origin(i),i=1,3)
      write (iunit,6000,err=900) exta,extb,extc
      write (iunit,6000,err=900) (grid(i),i=1,3)
      write (iunit,6010,err=900) (cell(i),i=1,6)
      write (iunit,6020,err=900) form,scale
c
 6000 format (3i5)
 6010 format (6f10.3)
 6020 format ('MAP ',a,' ',1pe14.6)
c
      write (iunit,form,err=900)
     +  (((map(i,j,k),i=1,exta), j=1,extb), k=1,extc)
c
      ierr = 0
      goto 910
c
  900 continue
      call errcon ('While writing OLDEZD file')
c
  910 continue
      close (iunit)
c
      return
      end
