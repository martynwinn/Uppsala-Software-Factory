c
c ================================================================
c
      subroutine xps_stroke ()
c
c ... stroke
c
      include 'xps.incl'
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      write (psunit,'(a1)') 'S'
c
      return
      end
