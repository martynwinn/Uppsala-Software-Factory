c
c
c
      subroutine wroncs (iunit,file,ncs,maxncs,rt,ierr)
c
c ... WRONCS - write one or many O NCS operators to a file
c
      implicit none
c
      integer iunit,ncs,maxncs,ierr,i,j,length
c
      real rt(12,maxncs)
c
      logical xinter
c
      character file*(*),line*128
c
code ...
c
      ierr = 0
      call xopxua (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening O NCS operator file')
        return
      end if
c
c ... write NCS operator(s)
c
      call stamp (line)
      write (iunit,'(9a)',err=900) '! ',line(1:length(line))
c
      do i=1,ncs
        write (line,'(a,i6)') '.LSQ_RT_',i
        call remspa (line)
        call appstr (line,'  R  12  (3f15.8)')
        write (iunit,'(a)',err=900) '!'
        write (iunit,'(a)',err=900) line(1:length(line))
        write (iunit,'(3f15.8)',err=900) (rt(j,i),j=1,12)
      end do
      goto 999
c
  900 continue
      call errcon ('While writing NCS operator file')
      ierr = -1
      close (iunit)
      return
c
  999 continue
      close (iunit)
      write (*,*)
      call jvalut (' Nr of NCS operators written :',1,ncs)
c
      return
      end
