c
c
c
      subroutine xvrml_wire (nat,xyz,maxbnd,nbonds,bndptr)
c
      include 'xvrml.incl'
c
      integer nat,maxbnd
c
      real xyz(3,nat)
c
      integer bndptr(maxbnd,nat)
      integer nbonds(nat)
      integer i,j,leng1
c
      character line*80
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
        write (ixvrml,6030,err=9999)
        do i=1,nat
          if (nbonds(i) .gt. 0) then
            do j=1,nbonds(i)
              write (line,6040) i-1,bndptr(j,i)-1
              call remspa (line)
              write (ixvrml,'(a)',err=9999) line(1:leng1(line))
            end do
          end if
        end do
        write (ixvrml,6050,err=9999)
        do i=1,nat
          if (nbonds(i) .eq. 0) then
            call xvrml_encode_rgb (rxrgb1,rxrgb2,rxrgb3,j)
            call xvrml_plus (1,xyz(1,i),j,0.3)
          end if
        end do
      end if
c
      return
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  IndexedLineSet { coordIndex [ ')
 6040 format (1x,i8,',',i8,',-1,')
 6050 format ('-1 ] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
