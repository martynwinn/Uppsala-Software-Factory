c
c
c
      subroutine xvrml_pointset (nat,xyz)
c
      include 'xvrml.incl'
c
      integer nat
c
      real xyz(3,nat)
c
      integer i,leng1
c
      character line*256
c
code ...
c
      if (lxvopn) then
        write (ixvrml,6000,err=9999)
        do i=1,nat
          write (line,6010) xyz(1,i),xyz(2,i),xyz(3,i)
          call pretty (line)
          write (ixvrml,'(a)',err=9999) line(1:leng1(line))
        end do
        write (ixvrml,6020,err=9999)
        write (ixvrml,6030,err=9999) nat
      end if
c
      return
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  PointSet { startIndex 0 '/
     +  '    numPoints ',i8,' } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
