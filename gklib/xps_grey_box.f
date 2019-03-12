c
c ================================================================
c
      subroutine xps_grey_box (x1,x2,y1,y2)
c
c ... draw a grey box
c
      include 'xps.incl'
c
      real x1,x2,y1,y2,px1,px2,py1,py2,zx,zy
c
      character type*10
c
code ...
c
      type = 'GreyBox'
      goto 10
c
      entry xps_dark_box (x1,x2,y1,y2)
      type = 'DarkBox'
      goto 10
c
      entry xps_light_box (x1,x2,y1,y2)
      type = 'LightBox'
      goto 10
c
c ... common for all entries
c
   10 continue
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
 6000 format (8(1x,f8.2),1x,a)
      write (psline,6000) px1,py1,px1,py2,px2,py2,px2,py1,type
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
c
      return
      end
