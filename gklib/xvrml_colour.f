c
c
c
      subroutine xvrml_colour (rgb1,rgb2,rgb3)
c
      include 'xvrml.incl'
c
      real rgb1,rgb2,rgb3
c
code ...
c
      rxrgb1 = rgb1
      rxrgb2 = rgb2
      rxrgb3 = rgb3
c
      if (lxvopn) then
        write (ixvrml,6000,err=9999) rxrgb1,rxrgb2,rxrgb3,
     +    rxrgb1,rxrgb2,rxrgb3
      end if
c
 6000 format ('Material { diffuseColor ',3f8.3/
     +        '          emissiveColor ',3f8.3,' }')
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
