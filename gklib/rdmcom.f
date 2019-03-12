c
c
c
      subroutine rdmcom (mask,size,space,lun,ierr)
c
c ... read mask in compressed format
c
      implicit none
c
      integer lun,size,ierr,i1,i2,i,ns,np,space
      integer mask(space)
c
      character line*80
c
code ...
c
      ierr = -1
c
      if (size .gt. space) then
        call errcon (' Space exceeded in RDMCOM')
        call jvalut (' Required :',1,size)
        call jvalut (' Maximum  :',1,space)
        close (lun)
        return
      end if
c
      ns = 0
      np = 0
c
      call jvalut (' Grid points :',1,size)
c
      do i=1,size
        mask (i) = 0
      end do
c
   10 continue
      read (lun,*,err=9000,end=100) i1,i2
      do i=i1,i2
        mask (i) = 1
      end do
      ns = ns + 1
      np = np + (i2-i1+1)
      goto 10
c
  100 continue
c
      call jvalut (' Stretches   :',1,ns)
      call jvalut (' Mask points :',1,np)
c
      ierr = 0
c
      return
c
 9000 continue
      backspace (lun)
      read (lun,'(a)') line
      call upcase (line)
      if (line(1:3) .eq. 'END') goto 100
c
      return
c
      end
