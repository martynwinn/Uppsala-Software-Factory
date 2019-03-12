c
c ================================================================
c
      subroutine xps_legend (gkline)
c
c ... add a line to the figure legend
c
      include 'xps.incl'
c
      character gkline*(*)
c
code ...
c
      if (gkline .eq. ' ') return
      psline = gkline
      call pretty (psline)
c
      if (psntxt .ge. maxtxt) then
        write (*,*) 'ERROR - Too many text lines for legend'
        return
      end if
c
      psntxt = psntxt + 1
      pstext (psntxt) = psline
c
      return
      end
