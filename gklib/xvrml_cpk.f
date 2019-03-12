c
c
c
      subroutine xvrml_cpk (nat,xyz,rad)
c
      include 'xvrml.incl'
c
      integer nat
c
      real xyz(3,nat),rad(nat)
c
      integer i
c
code ...
c
      if (lxvopn) then
        do i=1,nat
          write (ixvrml,6000,err=9999) rxrgb1,rxrgb2,rxrgb3,
     +      xyz(1,i),xyz(2,i),xyz(3,i),rad(i)
        end do
      end if
c
 6000 format ('Separator { '/
     +  '  Material {  diffuseColor ',3f8.3/
     +  '             emissiveColor 0.0 0.0 0.0 }'/
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
