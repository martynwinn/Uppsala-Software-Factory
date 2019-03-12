c
c
c
      subroutine xvrml_plus (nat,xyz,icol,size)
c
      include 'xvrml.incl'
c
      integer nat
c
      real xyz(3,nat)
      real size,r,g,b
c
      integer icol(nat)
      integer i
c
code ...
c
      if (lxvopn) then
        do i=1,nat
          call xvrml_decode_rgb (r,g,b,icol(i))
          write (ixvrml,6000,err=9999) r,g,b,
     +      xyz(1,i)-size,xyz(2,i),xyz(3,i),
     +      xyz(1,i)+size,xyz(2,i),xyz(3,i),
     +      xyz(1,i),xyz(2,i)-size,xyz(3,i),
     +      xyz(1,i),xyz(2,i)+size,xyz(3,i),
     +      xyz(1,i),xyz(2,i),xyz(3,i)-size,
     +      xyz(1,i),xyz(2,i),xyz(3,i)+size
        end do
      end if
c
      return
c
 6000 format ('Separator { '/
     +  '  Material { diffuseColor ',3(f6.3,1x),'}'/
     +  '  Coordinate3 { point [ ',2(3(f8.3,1x),',')/
     +  '  ',3(3(f8.3,1x),','),3(f8.3,1x),' ] }'/
     +  '  IndexedLineSet { coordIndex [ ',
     +  ' 0, 1, -1, 2, 3, -1, 4, 5, -1 ] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
