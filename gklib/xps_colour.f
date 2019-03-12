c
c ================================================================
c
      subroutine xps_colour (icol)
c
c ... set colour
c
      include 'xps.incl'
c
c ... nr of pre-defined colours
c
      integer maxcol
      parameter (maxcol=8)
c
      character cols(0:maxcol-1)*(10)
c
      real r,g,b
c
      integer icol
c
      data cols /'black','red','green','yellow','blue',
     +           'magenta','cyan','white'/
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      pscolr = icol
c
      if (pscolr .ge. 0 .and. pscolr .le. (maxcol-1)) then
        write (psunit,60) cols(pscolr)
      else if (pscolr .ge. maxcol) then
ccc        call errcon ('setrgbcolor not yet implemented')
ccc        write (psline,61) 0.0,0.0,0.0,'C'
        call xvrml_col_rgb (pscolr,r,g,b)
        write (psline,61) r,g,b,'C'
        call pretty (psline)
        write (psunit,60) psline(1:length(psline))
      else if (pscolr .lt. 0) then
        call errcon ('sethsvcolor not yet implemented')
        write (psline,61) 0.0,0.0,0.0,'K'
        call pretty (psline)
        write (psunit,60) psline(1:length(psline))
      end if
c
   60 format (A)
   61 format (3(f10.2,1x),a)
c
      return
      end
