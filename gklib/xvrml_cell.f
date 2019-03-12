c
c
c
      subroutine xvrml_cell (cell,iox,ioy,ioz,lline)
c
      include 'xvrml.incl'
c
      real xyz(3,8),cell(6),a(9),dum(3),dum2(3)
c
      integer ix,iy,iz,iox,ioy,ioz,leng1,l
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
        call orthog (cell,a,0)
        l = 0
        do ix=0,1
          dum(1) = float (ix+iox)
          do iy=0,1
            dum(2) = float (iy+ioy)
            do iz=0,1
              dum(3) = float (iz+ioz)
              call mulmtx (a,dum,dum2,3,3,1)
              write (line,6010) dum2(1),dum2(2),dum2(3)
              call pretty (line)
              write (ixvrml,'(a)',err=9999) line(1:leng1(line))
              l = l + 1
              xyz(1,l) = dum2(1)
              xyz(2,l) = dum2(2)
              xyz(3,l) = dum2(3)
            end do
          end do
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
c
        if (iox.eq.0.and.ioy.eq.0.and.ioz.eq.0) then
          call xvrml_text (xyz(1,1),xyz(2,1),xyz(3,1),'O')
          call xvrml_text (xyz(1,2),xyz(2,2),xyz(3,2),'c')
          call xvrml_text (xyz(1,3),xyz(2,3),xyz(3,3),'b')
          call xvrml_text (xyz(1,5),xyz(2,5),xyz(3,5),'a')
        end if
c
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
