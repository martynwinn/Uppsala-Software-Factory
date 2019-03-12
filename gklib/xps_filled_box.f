c
c ================================================================
c
      subroutine xps_filled_box (x1,x2,y1,y2,r,g,b)
c
c ... draw a coloured box
c
      include 'xps.incl'
c
      real x1,x2,y1,y2,px1,px2,py1,py2,zx,zy,r,g,b
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      call xps_mapol (x1,y1,zx,zy)
      px1 = psxmin + psxsca*(zx+psxoff-psrxin)
      py1 = psymin + psysca*(zy+psyoff-psryin)
c
      call xps_mapol (x2,y2,zx,zy)
      px2 = psxmin + psxsca*(zx+psxoff-psrxin)
      py2 = psymin + psysca*(zy+psyoff-psryin)
c
 6000 format (3(1x,f6.3),1x,a,8(1x,f8.2),1x,a)
      write (psline,6000) r,g,b,'C',px1,py1,px1,py2,
     +  px2,py2,px2,py1,'FilledBox'
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
c
      return
      end
