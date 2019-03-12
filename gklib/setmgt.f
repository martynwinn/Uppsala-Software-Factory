c
c
c
      subroutine setmgt (iunit,maptyp)
c
      implicit none
c
      integer iunit,maptyp
c
c ... for proper map closure
c
      integer mgtype(100)
      common /mgunit/ mgtype
      save /mgunit/
c
code ...
c
      if (iunit .gt. 0 .and. iunit .le. 100) then
        mgtype (iunit) = maptyp
      else
        call errcon ('SETMGT - Invalid unit number')
        call ivalut (' Unit :',1,iunit)
      end if
c
      return
      end
