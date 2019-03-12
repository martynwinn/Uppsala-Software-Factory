c
c ================================================================
c
      subroutine xps_ps_comment (text)
c
c ... add a comment to the PostScript file
c
      include 'xps.incl'
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
      write (psline,'(2a)') '% PROGAM COMMENT: ',text(1:length(text))
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
c
      return
      end
