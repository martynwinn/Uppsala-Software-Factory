c
c ================================================================
c
      subroutine xps_delete ()
c
c ... close and DELETE the file
c
      include 'xps.incl'
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      close (psunit,status='DELETE')
c
      return
      end
