c
c
c
      subroutine xvrml_stick (nat,xyz,maxbnd,nbonds,bndptr,strad)
c
      include 'xvrml.incl'
c
      integer nat,maxbnd
c
      real xyz(3,nat)
      real xc(3)
      real strad,a,b,c,d,e
c
      integer bndptr(maxbnd,nat)
      integer nbonds(nat)
      integer i,j,k
c
code ...
c
      if (lxvopn) then
        do i=1,nat
          if (nbonds(i) .gt. 0) then
            do j=1,nbonds(i)
              k = bndptr(j,i)
              call xvrml_get_cyl (xyz(1,i),xyz(1,k),a,b,c,d,e,xc)
              write (ixvrml,6000,err=9999)
     +          rxrgb1,rxrgb2,rxrgb3,
     +          xc(1),xc(2),xc(3),a,b,c,e,strad,d
            end do
          end if
        end do
      end if
c
      return
c
 6000 format ('Separator { '/
     +  '  Material {  diffuseColor ',3f8.3/
     +  '             emissiveColor 0.0 0.0 0.0 }'/
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
