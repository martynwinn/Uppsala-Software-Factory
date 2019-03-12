c
c
c
      subroutine xvrml_text (x,y,z,text)
c
      include 'xvrml.incl'
c
      real x,y,z
c
      integer leng1
c
      character line*256,text*(*)
c
code ...
c
      if (lxvopn) then
c
        if (leng1(text) .gt. 0) then
          write (line,6000,err=9999) x,y,z
          call pretty (line)
          write (ixvrml,'(a)',err=9999) line(1:leng1(line))
          write (ixvrml,6010,err=9999) text(1:leng1(text))
        end if
c
      end if
c
      return
c
 6000 format ('Separator { Translation { translation ',3f8.3,' }')
 6010 format (' AsciiText { string ["',a,'"] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
