c
c ================================================================
c
      subroutine xps_label (xlabel,ylabel)
c
c ... set label for X and Y
c
      include 'xps.incl'
c
      character*(*) xlabel,ylabel
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      psxlab = xlabel
      psylab = ylabel
c
      return
      end
