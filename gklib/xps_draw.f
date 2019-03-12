c
c ================================================================
c
      subroutine xps_draw (rx,ry)
c
c ... draw line segment
c
      include 'xps.incl'
c
      real rx,ry,zx,zy
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      call xps_mapol (rx,ry,zx,zy)
c
      psx = psxmin + psxsca*(zx+psxoff-psrxin)
      psy = psymin + psysca*(zy+psyoff-psryin)
c
      write (psline,61) psx,psy,'L'
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
c
      psnuml = psnuml + 1
c
c ... if 100 subsequent line draws, flush buffer (STROKE)
c     and move back to this point & reset counter
c
      if (psnuml .ge. 100) then
        write (psunit,'(a1)') 'S'
        write (psline,61) psx,psy,'M'
        call pretty (psline)
        write (psunit,'(a)') psline(1:length(psline))
        psnuml = 0
      end if
c
   61 format (f10.2,1x,f10.2,1x,a)
c
      return
      end
