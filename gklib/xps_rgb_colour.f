c
c ================================================================
c
      subroutine xps_rgb_colour (r,g,b)
c
c ... set colour
c
      include 'xps.incl'
c
      real r,g,b
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      if (r.gt.1.0 .and. g.gt.1.0 .and. b.gt.1.0) then
        r=r/256.
        g=g/256.
        b=b/256.
      end if
c
      if (r .lt. 0.0 .or. r .gt. 1.0 .or.
     +    g .lt. 0.0 .or. g .gt. 1.0 .or.
     +    b .lt. 0.0 .or. b .gt. 1.0) then
        write (*,*) 'ERROR - RGB value(s) outside 0-1 range'
        return
      end if
c
      write (psline,61) r,g,b,'C'
      call pretty (psline)
      write (psunit,60) psline(1:length(psline))
c
   60 format (A)
   61 format (3(f10.4,1x),a)
c
      return
      end
