c
c ================================================================
c
      subroutine xps_offset (xoff,yoff)
c
c ... define offset
c
      include 'xps.incl'
c
      real xoff,yoff
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      psxoff = xoff
      psyoff = yoff
c
      return
      end
