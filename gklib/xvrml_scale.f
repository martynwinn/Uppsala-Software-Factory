c
c
c
      subroutine xvrml_scale (sx,sy,sz)
c
      include 'xvrml.incl'
c
      real sx,sy,sz
c
      integer leng1
c
      character line*256
c
code ...
c
      if (lxvopn) then
        write (line,6000) sx,sy,sz
        call pretty (line)
        write (ixvrml,'(a)',err=9999) line(1:leng1(line))
      end if
c
      return
c
 6000 format ('Scale { scaleFactor ',3f20.10,' }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
