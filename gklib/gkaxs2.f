c
c=====================================================
c
      subroutine gkaxs2 (xnow,mytext)
c
c ... XTXGK2 (...) - generate axis-label string from float
c
      implicit none
c
      real zero
      parameter (zero = 0.0)
c
      real xnow
c
      character mytext*(*)
c
      entry xtxgk2 (xnow,mytext)
c
code ...
c
      if ((xnow .ge. 0.001 .and. xnow .le. 9999.999) .or.
     +    (xnow .le. -0.001 .and. xnow .ge. -9999.999) .or.
     +    (xnow .eq. zero)) then
c
        xnow = 0.0001*float (nint(10000.0*xnow))
c
        write (mytext,'(f10.3)') xnow
        if (mytext(10:10) .eq. '0') then
          mytext(10:10) = ' '
          if (mytext(9:9) .eq. '0') then
            mytext(9:9) = ' '
            if (mytext(8:8) .eq. '0') then
              mytext(7:8) = '  '
            end if
          end if
        end if
      else
        write (mytext,'(1pe8.1)') xnow
      end if
      call remspa (mytext)
c
      return
      end
