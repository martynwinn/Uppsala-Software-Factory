c
c
c
      subroutine xvrml_face_surf (nx,ny,z,x1,x2,y1,y2)
c
      include 'xvrml.incl'
c
      integer nx,ny
c
      real z(nx,ny)
      real x1,x2,y1,y2,z1,z2,dx,dy,sx,sy,sz
c
      integer i,j,leng1,i1,i2,j1,j2
c
      character line*256
c
code ...
c
      if (lxvopn) then
c
        z1 = z(1,1)
        z2 = z(1,1)
        do i=1,nx
          do j=1,ny
            if (z(i,j) .lt. z1) z1 = z(i,j)
            if (z(i,j) .gt. z2) z2 = z(i,j)
          end do
        end do
        sx = 100.0/(x2-x1)
        sy = 100.0/(y2-y1)
        sz = 100.0/(z2-z1)
        call xvrml_scale (sx,sy,sz)
c
        dx = (x2-x1)/float(nx-1)
        dy = (y2-y1)/float(ny-1)
        write (ixvrml,6000,err=9999)
        do i=1,nx
          do j=1,ny
            write (line,6010) x1+float(i)*dx,y1+float(j)*dy,z(i,j)
            call pretty (line)
            write (ixvrml,'(a)',err=9999) line(1:leng1(line))
          end do
        end do
        write (ixvrml,6020,err=9999)
        write (ixvrml,6030,err=9999)
c
        do j=1,ny-1
          do i=1,nx-1
c
            i1 = (j-1)*nx+i-1
            i2 = i1 + 1
            j1 = j*nx+i-1
            j2 = j1+1
c
            write (line,6040) i1,j1,i2,i1,-1
            call pretty (line)
            write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
            write (line,6040) i2,j1,j2,i2,-1
            call pretty (line)
            write (ixvrml,'(a)',err=9999) line(1:leng1(line))
c
          end do
        end do
c
        write (ixvrml,6050,err=9999)
      end if
c
      return
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
c 6030 format ('  IndexedLineSet { coordIndex [ ')
 6030 format ('  IndexedFaceSet { coordIndex [ ')
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
