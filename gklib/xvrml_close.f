c
c
c
      subroutine xvrml_close ()
c
      include 'xvrml.incl'
c
code ...
c
      if (lxvopn) then
        close (ixvrml)
        lxvopn = .false.
      end if
      call prompt (' Closed VRML file')
c
      return
      end
