c
c
c
      subroutine xvrml_box8 (xyz,lline)
c
c ... 8 points are: 000, 001, 010, 011, 100, 101, 110, 111
c     (so they are easy to generate with 3 nested DO loops)
c
      include 'xvrml.incl'
c
      real xyz(3,8)
c
      integer i,j,leng1
c
      logical lline
c
      character line*256
c
code ...
c
      if (lxvopn) then
c
        write (ixvrml,6000,err=9999)
        do i=1,8
          write (line,6010) (xyz(j,i),j=1,3)
          call pretty (line)
          write (ixvrml,'(a)',err=9999) line(1:leng1(line))
        end do
        write (ixvrml,6020,err=9999)
        if (lline) then
          write (ixvrml,6030,err=9999)
        else
          write (ixvrml,6031,err=9999)
        end if
c
        write (line,6040) 0,2,6,4,0,-1
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
        write (line,6040) 1,5,7,3,1,-1
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
        write (line,6040) 0,1,3,2,0,-1
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
        write (line,6040) 2,3,7,6,2,-1
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
        write (line,6040) 4,6,7,5,4,-1
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
        write (line,6040) 0,4,5,1,0,-1
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
        write (ixvrml,6050,err=9999)
      end if
c
      return
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  IndexedLineSet { coordIndex [ ')
 6031 format ('  IndexedFaceSet { coordIndex [ ')
 6040 format (20(i8,','))
 6050 format ('-1 ] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
