c
c ================================================================
c
      subroutine xps_axes ()
c
c ... figure out divisions along the two axes
c
      include 'xps.incl'
c
      real xperti,yperti,xnow,axnow,xgk,ynow,aynow
c
      integer nxtick,nytick,ierror,it,icol
c
      character gktext*20
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
c ... labels in black
c
      write (psunit,'(a)') 'black'
c
c ... label X-axis
c
      call xtxgk1 (psrxin,psrxax,maxtix,xperti,nxtick,ierror)
      do it=0,nxtick+1
        xnow = float(int((psrxin + float(it)*xperti)/xperti)) * xperti
        if (xnow .ge. psrxin .and. xnow .le. psrxax) then
          axnow = xnow
          psx = psxmin + psxsca*(xnow-psrxin)
          write (psline,61) psx,psymin,'M'
          call pretty (psline)
          write (psunit,'(a)') psline(1:length(psline))
          write (psline,61) psx,psymin-psytsi,'L'
          call pretty (psline)
          write (psunit,'(a)') psline(1:length(psline))
          call xtxgk2 (axnow,gktext)
          xgk = fudge1*float(length(gktext))
          psx = psxmin + psxsca*(xnow-psrxin) - xgk
          write (psline,61) psx,psymin-2*psytsi-fudge3,'M'
          call pretty (psline)
          write (psunit,'(a)') psline(1:length(psline))
          write (psunit,60)
     +      '('//gktext(1:length(gktext))//') 12.00 Print'
        end if
      end do
c
c ... label X-axis
c
      call xtxgk1 (psryin,psryax,maxtiy,yperti,nytick,ierror)
      do it=0,nytick+1
        ynow = float( int( (psryin + float(it)*yperti) / yperti) )
     +         * yperti
        if (ynow .ge. psryin .and. ynow .le. psryax) then
          aynow = ynow
          psy = psymin + psysca*(ynow-psryin)
          write (psline,61) psxmax,psy,'M'
          call pretty (psline)
          write (psunit,'(a)') psline(1:length(psline))
          write (psline,61) psxmax+psxtsi,psy,'L'
          call pretty (psline)
          write (psunit,'(a)') psline(1:length(psline))
          call xtxgk2 (aynow,gktext)
          write (psline,61) psxmax+psxtsi+fudge2,psy-fudge3,'M'
          call pretty (psline)
          write (psunit,'(a)') psline(1:length(psline))
          write (psunit,60)
     +      '('//gktext(1:length(gktext))//') 12.00 Print'
        end if
      end do
c
c ... restore colour
c
      icol = pscolr
      call xps_colour (icol)
c
   60 format (A)
   61 format (f10.2,1x,f10.2,1x,a)
c
      return
      end
