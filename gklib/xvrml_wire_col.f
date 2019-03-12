c
c
c
      subroutine xvrml_wire_col (nat,xyz,maxbnd,nbonds,bndptr,icol)
c
      include 'xvrml.incl'
c
      integer nat,maxbnd
c
      real xyz(3,nat)
      real r,g,b,rr,gg,bb
c
      integer bndptr(maxbnd,nat)
      integer nbonds(nat),icol(nat)
      integer i,j,k
c
code ...
c
      if (lxvopn) then
        do i=1,nat
          if (nbonds(i) .gt. 0) then
            do j=1,nbonds(i)
              k = bndptr(j,i)
c
ccc              print *,i,k,icol(i),icol(k)
c
              if (icol(i) .eq. icol(k)) then
                call xvrml_decode_rgb (r,g,b,icol(i))
                write (ixvrml,6000,err=9999) r,g,b,
     +            xyz(1,i),xyz(2,i),xyz(3,i),
     +            xyz(1,k),xyz(2,k),xyz(3,k)
              else
                call xvrml_decode_rgb (r,g,b,icol(i))
                call xvrml_decode_rgb (rr,gg,bb,icol(k))
                write (ixvrml,6010,err=9999)
     +            xyz(1,i),xyz(2,i),xyz(3,i),
     +            0.5*(xyz(1,i)+xyz(1,k)),
     +            0.5*(xyz(2,i)+xyz(2,k)),
     +            0.5*(xyz(3,i)+xyz(3,k)),
     +            xyz(1,k),xyz(2,k),xyz(3,k),
     +            r,g,b,rr,gg,bb
              end if
            end do
          else if (nbonds(i) .eq. 0) then
            call xvrml_plus (1,xyz(1,i),icol(i),0.3)
          end if
        end do
      end if
c
      return
c
 6000 format ('Separator {'/
     +  '  Material { diffuseColor ',3(f6.3,1x),'}'/
     +  '  Coordinate3 { point [ ',3(f8.3,1x),','/
     +  '  ',3(f8.3,1x),' ] }'/
     +  '  IndexedLineSet { coordIndex [ 0, 1, -1 ] } }')
c
 6010 format ('Separator { Coordinate3 { point [ ',3(f8.3,1x),','/
     +  '  ',3(f8.3,1x),',',3(f8.3,1x),' ] }'/
     +  '  Material { diffuseColor ',3(f6.3,1x),'}'/
     +  '  IndexedLineSet { coordIndex [ 0, 1, -1 ] }'/
     +  '  Material { diffuseColor ',3(f6.3,1x),'}'/
     +  '  IndexedLineSet { coordIndex [ 1, 2, -1 ] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
