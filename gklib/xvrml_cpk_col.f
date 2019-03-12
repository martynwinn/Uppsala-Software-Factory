c
c
c
      subroutine xvrml_cpk_col (nat,xyz,rad,icol)
c
      include 'xvrml.incl'
c
      integer nat
c
      real xyz(3,nat),rad(nat)
      real r,g,b
c
      integer icol(nat)
      integer i
c
code ...
c
      if (lxvopn) then
        do i=1,nat
          call xvrml_decode_rgb (r,g,b,icol(i))
          write (ixvrml,6000,err=9999) r,g,b,xyz(1,i),
     +      xyz(2,i),xyz(3,i),rad(i)
        end do
      end if
c
 6000 format ('Separator { '/
     +  '  Material { diffuseColor ',3(f6.3,1x)/
     +  '            emissiveColor 0.0 0.0 0.0 }'/
     +  '  Translation { translation ',3(f8.3,1x),'}'/
     +  '  Sphere { radius ',f8.3,' } }')
c
      return
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
