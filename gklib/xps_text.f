c
c ================================================================
c
      subroutine xps_text (rx,ry,ptsize,text)
c
c ... print text with Times-Roman font
c
      include 'xps.incl'
c
      real rx,ry,ptsize
c
      character text*(*)
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      if (length(text) .le. 0) return
c
      call xps_move (rx,ry)
      write (psline,'(a,a,a,f10.2,a)')
     +  '(',text(1:length(text)),') ',max(1.0,ptsize),' Print'
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
      call xps_stroke ()
c
      return
      end
