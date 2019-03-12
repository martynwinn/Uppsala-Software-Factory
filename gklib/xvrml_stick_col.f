c
c
c
      subroutine xvrml_stick_col (nat,xyz,maxbnd,nbonds,bndptr,
     +                            strad,icol)
c
      include 'xvrml.incl'
c
      integer nat,maxbnd
c
      real xyz(3,nat)
      real xc(3),xd(3)
      real strad,a,b,c,d,e,r,g,bb
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
              if (icol(i) .eq. icol(k)) then
                call xvrml_decode_rgb (r,g,bb,icol(i))
                call xvrml_get_cyl (xyz(1,i),xyz(1,k),a,b,c,d,e,xc)
                write (ixvrml,6000,err=9999) r,g,bb,
     +            xc(1),xc(2),xc(3),a,b,c,e,strad,d
              else
                call xvrml_decode_rgb (r,g,bb,icol(i))
                call xvrml_get_cyl (xyz(1,i),xyz(1,k),a,b,c,d,e,xc)
                call xvrml_get_cyl (xyz(1,i),xc,a,b,c,d,e,xd)
                write (ixvrml,6000,err=9999) r,g,bb,
     +            xd(1),xd(2),xd(3),a,b,c,e,strad,d
                call xvrml_decode_rgb (r,g,bb,icol(k))
                call xvrml_get_cyl (xyz(1,i),xyz(1,k),a,b,c,d,e,xc)
                call xvrml_get_cyl (xc,xyz(1,k),a,b,c,d,e,xd)
                write (ixvrml,6000,err=9999) r,g,bb,
     +            xd(1),xd(2),xd(3),a,b,c,e,strad,d
              end if
            end do
          end if
        end do
      end if
c
      return
c
 6000 format ('Separator { '/
     +  '  Material { diffuseColor ',3(f6.3,1x)/
     +  '            emissiveColor 0.0 0.0 0.0 }'/
     +  '  Translation { translation ',3f8.3,' }'/
     +  '  Rotation { rotation ',3f9.5,1x,f10.5,' }'/
     +  '  Cylinder { radius ',f8.3,' height ',f8.3,' } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
