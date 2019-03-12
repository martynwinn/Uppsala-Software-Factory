c
c ================================================================
c
      subroutine xps_close ()
c
c ... write last bits & close the file
c
      include 'xps.incl'
c
      integer i,icol
c
      real ty,px1,px2,py1,py2
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
   60 format (A)
   61 format (f10.2,1x,f10.2,1x,a)
c
      write (psunit,60)
     +  '% start of XPS_GRAF footer stuff'
c
c ... stroke & restore solid line type
c
      write (psunit,60) 'S [1 0] 0 setdash'
c
c ... white out anything that ended up outside the plot area
c     (yeah, I know, poor and lazy way of hidden line removal ...)
c
      if (pshide) then
c
        icol = pscolr
        call xps_colour (7)
c
        px1 = -100.0
        px2 = 1000.0
c
        write (psunit,60) '% white out everything below plot box'
        write (psline,6000) px1,px1,px2,px1,
     +    px2,psymin,px1,psymin,'FilledBox'
        call pretty (psline)
        write (psunit,'(a)') psline(1:length(psline))
c
        write (psunit,60) '% white out everything above plot box'
        write (psline,6000) px1,psymax,px2,psymax,
     +    px2,px2,px1,px2,'FilledBox'
        call pretty (psline)
        write (psunit,'(a)') psline(1:length(psline))
c
        write (psunit,60) '% white out everything left of box'
        write (psline,6000) px1,px1,psxmin,px1,
     +    psxmin,px2,px1,px2,'FilledBox'
        call pretty (psline)
        write (psunit,'(a)') psline(1:length(psline))
c
        write (psunit,60) '% white out everything right of box'
        write (psline,6000) psxmax,px1,px2,px1,
     +    px2,px2,psxmax,px2,'FilledBox'
        call pretty (psline)
        write (psunit,'(a)') psline(1:length(psline))
c
        call xps_colour (icol)
c
      end if
c
 6000 format (8(1x,f8.2),1x,a)
c
c ... label axes
c
      write (psunit,60) '% label X axis and Y axis'
      call xps_axes ()
c
c ... stroke & restore black colour
c
      write (psunit,60) 'black'
c
c ... draw corner at (25,725) and (525,25)
c
      write (psunit,60) '% little corners top-left and bottom-right'
      write (psunit,61) 25.,715.,'M'
      write (psunit,61) 25.,725.,'L'
      write (psunit,61) 35.,725.,'L S'
      write (psunit,61) 515.,25.,'M'
      write (psunit,61) 525.,25.,'L'
      write (psunit,61) 525.,35.,'L S'
c
c ... draw axis labels & file name
c
      write (psunit,60) '% axis labels'
      write (psunit,61) psxmin,psymin-40.0,'M'
cc      write (psunit,60)
cc     +  '(X = '//psxlab(1:length(psxlab))//') 15.00 Center'
      write (psunit,60)
     +  '(X = '//psxlab(1:length(psxlab))//') 15.00 Print S'
c
      write (psunit,61) psxmin,psymin-60.0,'M'
cc      write (psunit,60)
cc     +  '(Y = '//psylab(1:length(psylab))//') 15.00 CenterRot90'
      write (psunit,60)
     +  '(Y = '//psylab(1:length(psylab))//') 15.00 Print S'
c
      write (psunit,60) '% plot header'
      write (psunit,61) psxmin,psymax+20.0,'M'
      write (psunit,60)
     +  '(File '//psfile(1:length(psfile))//') 15.00 Bold S'
c
      write (psunit,61) psxmin,psymax+40.0,'M'
      write (psunit,60)
     +  '(D.) 15.00 Greek S'
c
      write (psunit,61) psxmin+20,psymax+40.0,'M'
      write (psunit,60)
     +  '(Created by '//psprog(1:length(psprog))//
     +  ' at '//psdate(1:length(psdate))//' for '//
     +  psuser(1:length(psuser))//' *** XPS_GRAF '//
     +  psvers//' *** (c) GJ Kleywegt, 1992-2009 ***) 9.00 Print S'
c
c ... legend
c
      if (psntxt .gt. 0) then
        write (psunit,60) '% remarks as legend'
        do i=1,psntxt
          ty = psymin - 80.0 - float ( (i-1)*15 )
          write (psunit,61) psxmin,ty,'M'
          write (psunit,60)
     +      '(' // pstext(i)(1:length(pstext(i))) //
     +      ') 12.00 Print S'
        end do
      end if
c
c ... last bits for the PostScript file
c
      write (psunit,60) '% almost done now'
      write (psunit,60) 'grestore stroke'
      write (psunit,60) 'showpage'
      write (psunit,60) '%%Trailer'
c
      close (psunit)
c
c ... reset legend counter
c
      psntxt = 0
c
      psdash = 0
      pscolr = 0
      psopen = .false.
      pshide = .false.
c
      return
      end
