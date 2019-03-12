c
c
c
      subroutine xvrml_polyline (nat,xyz)
c
      include 'xvrml.incl'
c
      integer nat
c
      real xyz(3,nat)
c
      integer i,j,leng1,imax
c
      character line*256
c
code ...
c
      if (lxvopn) then
        if (nat .eq. 1) then
          call xvrml_encode_rgb (rxrgb1,rxrgb2,rxrgb3,j)
          call xvrml_plus (1,xyz(1,1),j,0.3)
        else
          write (ixvrml,6000,err=9999)
          do i=1,nat
            write (line,6010) xyz(1,i),xyz(2,i),xyz(3,i)
            call pretty (line)
            write (ixvrml,'(a)',err=9999) line(1:leng1(line))
          end do
          write (ixvrml,6020,err=9999)
          write (ixvrml,6030,err=9999)
          do i=1,nat,20
            imax = min(nat,i+19)
            write (line,6040) (j-1,j=i,imax)
            call remspa (line)
            write (ixvrml,'(a)',err=9999) line(1:leng1(line))
          end do
          write (ixvrml,6050,err=9999)
        end if
      end if
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
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
