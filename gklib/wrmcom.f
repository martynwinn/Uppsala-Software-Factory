c
c
c
      subroutine wrmcom (mask,size,lun,ierr)
c
c ... write mask in compressed format
c
      implicit none
c
      integer lun,size,ierr,i1st,i,ns,np,length
      integer mask(size)
c
      logical l1
c
      character line*80
c
code ...
c
      ierr = -1
c
      call jvalut (' Grid points :',1,size)
c
      l1 = .false.
      i1st = 0
      ns = 0
      np = 0
c
      do i=1,size
        if (l1) then
          if (mask(i) .ne. 1) then
            write (line,*) i1st,i-1
            call pretty (line)
            write (lun,*,err=9000) (line(1:length(line)))
            l1 = .false.
            ns = ns + 1
            np = np + (i-i1st)
          end if
        else
          if (mask(i) .eq. 1) then
            l1 = .true.
            i1st = i
          end if
        end if
      end do
c
      if (l1) then
        write (lun,*,err=9000) i1st,size
        ns = ns + 1
        np = np + (i-i1st)
      end if
c
      call jvalut (' Stretches   :',1,ns)
      call jvalut (' Mask points :',1,np)
c
      ierr = 0
c
      return
c
 9000 continue
c
      return
c
      end
