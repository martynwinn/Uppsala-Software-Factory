c
c ================================================================
c
      subroutine xps_move (rx,ry)
c
c ... move
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
      write (psline,61) 'S',psx,psy,'M'
      psnuml = 0
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
c
   61 format (a,1x,f10.2,1x,f10.2,1x,a)
c
      return
      end
