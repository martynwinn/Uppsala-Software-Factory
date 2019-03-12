c
c
c
      subroutine vrmlca (iunit,filnam,vrdist,ivrml,cavrml,ierr)
c
      implicit none
c
      integer iunit,ivrml,ierr,nca,ntot,i,j,imax,leng1
c
      real x(3),xo(3)
      real vrdist,xdum,distce
c
      logical xinter
c
      character filnam*(*),cavrml*4,line*256
c
code ...
c
      call xopxoa (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file')
        return
      end if
c
      nca = 0
      ntot = 0
c
   10 continue
      read (iunit,'(a)',err=10,end=20) line
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') goto 10
      if (line(13:16) .ne. cavrml) goto 10
c
cATOM     29  C   SER     4      24.627  20.568  41.091  1.00 10.39      1CBS 251
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      read (line(31:54),'(3f8.3)',err=10,end=10) x(1),x(2),x(3)
      if (nca .gt. 0) then
        xdum = distce (x,xo)
        if (xdum .gt. vrdist) then
          write (ivrml,6020,err=9999)
          write (ivrml,6030,err=9999)
          do i=1,nca,20
            imax = min(nca,i+19)
            write (line,6040) (j-1,j=i,imax)
            call remspa (line)
            write (ivrml,'(a)',err=9999) line(1:leng1(line))
          end do
          write (ivrml,6050,err=9999)
          nca = 0
        end if
      end if
c
      nca = nca + 1
      ntot = ntot + 1
      if (nca .eq. 1) then
        write (ivrml,6000,err=9999)
      end if
c
      write (line,6010) x(1),x(2),x(3)
      call pretty (line)
      write (ivrml,'(a)',err=9999) line(1:leng1(line))
      xo(1) = x(1)
      xo(2) = x(2)
      xo(3) = x(3)
c
      goto 10
c
   20 continue
c
      if (nca .gt. 0) then
        write (ivrml,6020,err=9999)
        write (ivrml,6030,err=9999)
        do i=1,nca,20
          imax = min(nca,i+19)
          write (line,6040) (j-1,j=i,imax)
          call remspa (line)
          write (ivrml,'(a)',err=9999) line(1:leng1(line))
        end do
        write (ivrml,6050,err=9999)
      end if
c
      call jvalut (' Nr of atoms written :',1,ntot)
      ierr = 0
c
      return
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  IndexedLineSet { coordIndex [ ')
 6040 format (20(i8,','))
 6050 format ('-1 ] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      ierr = -1
c
      return
      end
