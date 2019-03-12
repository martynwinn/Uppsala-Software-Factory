c
c
c
      subroutine rdoncs (iunit,file,ncs,maxncs,rt,ierr)
c
c ... RDONCS - read one or many O NCS operators from one file
c
      implicit none
c
      integer iunit,ncs,maxncs,ierr,maxopt,nopt,i
      parameter (maxopt = 4)
c
      real rt(12,maxncs)
c
      logical xinter
c
      character file*(*),line*128,optpar(maxopt)*80
c
code ...
c
      ierr = 0
      call xopxoa (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening O NCS operator file')
        return
      end if
c
c ... look for next NCS operator
c
   10 continue
      if (ncs .ge. maxncs) then
        call errcon ('Maximum number of NCS operators read')
        call jvalut (' Maximum :',1,maxncs)
        ierr = 0
        return
      end if
c
c ... scan file
c
   20 continue
      read (iunit,'(a)',end=900,err=990) line
c
c ... ignore comments
c
      if (line (1:1) .eq. '!') goto 20
c
      call pretty (line)
      call upcase (line)
c
c ... parse datablock header
c
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (ierr .ne. 0) goto 20
      if (nopt .ne. 4) goto 20
      if (optpar(2)(1:1) .ne. 'R') goto 20
      if (optpar(3)(1:2) .ne. '12') goto 20
c
c ... read operator
c
c ... 970506 - jump to 990 (instead of 900) on unexpected EOF !!!
c
      ncs = ncs + 1
      write (*,*)
      call textut (' ==> Read NCS operator :',optpar(1))
      read (iunit,optpar(4),end=990,err=990) (rt(i,ncs),i=1,12)
c
      call anancs (1,rt(1,ncs),.true.,ierr)
      if (ierr .ne. 0) then
        call errstp ('RDONCS - In NCS operator; removed')
        ncs = ncs - 1
      end if
c
      goto 10
c
  990 continue
      call errcon ('While reading NCS operator file')
      call prompt (' Check the file format and make sure the file')
      call prompt (' ends with a <NEWLINE> !')
      ierr = -1
      goto 999
c
  900 continue
      ierr = 0
c
  999 continue
      close (iunit)
      write (*,*)
      call jvalut (' Nr of NCS operators now :',1,ncs)
c
      return
      end
