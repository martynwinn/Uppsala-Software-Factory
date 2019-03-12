c
c
c
      subroutine mapclo (iunit)
c
c ... close a map which has been opened for reading
c
      integer iunit
c
c ... for proper map closure
c
      integer mgtype(100)
      common /mgunit/ mgtype
      save /mgunit/
c
code ...
c
      if (iunit .lt. 0 .or. iunit .gt. 100) then
        call errcon ('MAPCLO - invalid unit number')
        call jvalut (' I/O unit :',1,iunit)
        return
      end if
c
      if (mgtype(iunit) .le. 0) then
        call errcon ('MAPCLO - unit not open')
        call ivalut (' I/O unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 1) then
        call ivalut (' Closing PROTEIN map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 2) then
        call ivalut (' Closing FFT-Y map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 3) then
        call ivalut (' Closing TENEYCK2 map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 4) then
        call ivalut (' Closing BINARY CCP4 map on unit :',1,iunit)
        call mrclos (iunit)
      else if (mgtype(iunit) .eq. -4) then
        call ivalut (' Closing ASCII CCP4 map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 5) then
        call ivalut (' Closing XPLOR map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 6) then
        call ivalut (' Closing OLDEZD map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 7) then
        call ivalut (' Closing MASK on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 8) then
        call ivalut (' Closing (NEW)EZD map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 9) then
        call ivalut (' Closing BINXPLOR map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 13) then
        call ivalut (' Closing TNT map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 14) then
        call ivalut (' Closing PHASES map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 15) then
        call ivalut (' Closing FSMASK map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 20) then
        call ivalut (' Closing EM08 map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 22) then
        call ivalut (' Closing MPI (EM) map on unit :',1,iunit)
        close (iunit)
      else if (mgtype(iunit) .eq. 23) then
        call ivalut (
     +    ' Closing AMBER-style X-PLOR map on unit :',1,iunit)
        close (iunit)
      else
        call errcon ('MAPCLO - unknown map type')
        call ivalut (' Map type :',1,mgtype(iunit))
        close (iunit)
      end if
c
      mgtype (iunit) = 0
c
      return
      end
