c
c ================================================================
c
      subroutine xps_dash ()
c
c ... change line dash type
c
      include 'xps.incl'
c
      integer ib,iw,is
c
      save ib,iw,is
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      psdash = psdash + 1
c
c ... first call sets line type to solid
c
      if (psdash .eq. 1) then
        ib = 1
        iw = 0
        is = 1
c
c ... subsequent calls give Black/White 2 2, 3 2, 4 2, 5 3,
c     6 3, 7 4, 8 4, 9 4, 10 4, ...
c
      else
        ib = ib + 1
        if (ib .gt. (2*iw)) iw = iw + 1
        iw = max (4, iw)
        is = max (1, ib/2)
      end if
c
      call xps_stroke ()
      write (psline,'(a,i3,a,i3,a,i3,a)')
     +  '[ ',ib,' ',iw,' ] ',is,' setdash'
      call pretty (psline)
      write (psunit,'(a)') psline(1:length(psline))
c
ccc      print *,psline(1:length(psline))
c
      return
      end
