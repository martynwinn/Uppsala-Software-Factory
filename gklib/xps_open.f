c
c ================================================================
c
      subroutine xps_open (iunit,fname,pname)
c
c ... open the file and write the initial bits
c
      include 'xps.incl'
c
      integer iunit,ll
c
      logical xinter
c
      character fname*(*),pname*(*)
c
code ...
c
      if (iunit .gt. 0) psunit = iunit
c
      if (length(fname) .lt. 1) then
        psfile = 'out.ps'
        call textin (' PostScript file ?',psfile)
      else
        psfile = fname
      end if
c
      if (length(pname) .lt. 1) then
        psprog = 'UNKNOWN PROGRAM'
      else
        psprog = pname
      end if
c
c ... open PS file as unknown
c
      call xopxua (psunit,psfile,xinter(),pserror)
      if (pserror .ne. no_error) then
        write (*,*) 'ERROR - Could not open output PostScript file'
        return
      end if
c
c      write (*,*) 'PostScript file opened; writing header'
c
      psopen = .true.
      call gkdate (psdate)
c
      call gkuser (psuser)
      ll = length (psuser)
      if (ll .lt. 1) psuser = 'an UNKNOWN user'
c
      call textut (psrout,psvers)
      call textut (' Opened PostScript file :',psfile)
      call textut (' Date    :',psdate)
      call textut (' User    :',psuser)
      call textut (' Program :',psprog)
c
   60 format (a)
   61 format (f10.2,1x,f10.2,1x,a)
c
      write (psunit,60)
     +  '%!PS-Adobe-3.0 EPSF-3.0'
      write (psunit,60)
     +  '%%Title: '//psfile(1:length(psfile))
      write (psunit,60)
     +  '%%Creator: '//psprog(1:length(psprog))//
     +  ' (c) 1992-2009 Gerard J. Kleywegt'
      write (psunit,60)
     +  '%%CreationDate: '//psdate(1:length(psdate))
      write (psunit,60)
     +  '%%DocumentFonts: Times-Roman Times-Bold Symbol'
      write (psunit,60)
     +  '%%Pages: 1'
      write (psunit,60)
     +  '%%BoundingBox: 10.0 10.0 535.0 785.0'
      write (psunit,60)
     +  '%%EndComments'
      write (psunit,60)
     +  '/Line { moveto lineto stroke } bind def'
      write (psunit,60)
     +  '/M { moveto } bind def'
      write (psunit,60)
     +  '/L { lineto } bind def'
      write (psunit,60)
     +  '/S { stroke } bind def'
      write (psunit,60)
     +  '/C { stroke setrgbcolor } bind def'
      write (psunit,60)
     +  '/K { stroke sethsvcolor } bind def'
      write (psunit,60)
     +  '/black {0 0 0 C} bind def'
      write (psunit,60)
     +  '/red {1 0 0 C} bind def'
      write (psunit,60)
     +  '/green {0 1 0 C} bind def'
      write (psunit,60)
     +  '/yellow {1 1 0 C} bind def'
      write (psunit,60)
     +  '/blue {0 0 1 C} bind def'
      write (psunit,60)
     +  '/magenta {1 0 1 C} bind def'
      write (psunit,60)
     +  '/cyan {0 1 1 C} bind def'
      write (psunit,60)
     +  '/white {1 1 1 C} bind def'
      write (psunit,60)
     +  '/Center {'
      write (psunit,60)
     +  'dup /Times-Roman findfont exch scalefont setfont'
      write (psunit,60)
     +  'exch stringwidth pop -2 div exch -3 div rmoveto } bind def'
      write (psunit,60)
     +  '/Print { /Times-Roman findfont exch scalefont '//
     +  'setfont show } bind def'
      write (psunit,60)
     +  '/Bold { /Times-Bold findfont exch scalefont '//
     +  'setfont show } bind def'
      write (psunit,60)
     +  '/Greek { /Symbol findfont exch scalefont '//
     +  'setfont show } bind def'
      write (psunit,60)
     +  '/CenterRot90 {'
      write (psunit,60)
     +  'dup /Times-Roman findfont exch scalefont setfont'
      write (psunit,60)
     +  'exch stringwidth pop -2 div exch 3 div exch rmoveto'
      write (psunit,60)
     +  '} bind def'
      write (psunit,60)
     +  '/FilledBox { gsave M L L L fill S '//
     +  'grestore } bind def'
      write (psunit,60)
     +  '/GreyBox { gsave 0.5 setgray M L L L fill S '//
     +  'grestore } bind def'
      write (psunit,60)
     +  '/DarkBox { gsave 0.2 setgray M L L L fill S '//
     +  'grestore } bind def'
      write (psunit,60)
     +  '/LightBox { gsave 0.8 setgray M L L L fill S '//
     +  'grestore } bind def'
      write (psunit,60)
     +  '%%EndProlog'
      write (psunit,60)
     +  '% change next line to move the plot on the page'
      write (psunit,60)
     +  '     0.0 0.0 translate'
      write (psunit,60)
     +  '% change next line to scale the plot'
      write (psunit,60)
     +  '     1.0 1.0 scale'
      write (psunit,60)
     +  '     1 setlinecap 1 setlinejoin'
      write (psunit,60)
     +  '     1 setlinewidth 0 setgray'
      write (psunit,60)
     +  '%%Page: 1 1'
      write (psunit,60)
     +  '% begin plot box'
      write (psunit,60)
     +  'black'
      write (psunit,61) psxmin,psymin,'M'
      write (psunit,61) psxmax,psymin,'L'
      write (psunit,61) psxmax,psymax,'L'
      write (psunit,61) psxmin,psymax,'L'
      write (psunit,60)
     +  'closepath gsave'
      write (psunit,60)
     +  '% end plot box'
      write (psunit,60)
     +  'gsave 1.0000 setgray fill grestore S'
      write (psunit,60)
     +  '% end of XPS_GRAF header stuff'
c
cc      write (psunit,60)
cc     +  'clip newpath'
c
      return
      end
