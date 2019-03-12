c
c=====================================================
c
      subroutine gkaxs1 (xlo,xhi,nxmax,xperti,nxtick,ierr)
c
c ... XTXGK1 (...) - calculate spacing for tick-marks
c
      implicit none
c
      real zero
      parameter (zero = 0.0)
c
      integer xtxdbg
      parameter (xtxdbg = 0)
c
      real xlo,xhi,xperti,xrange,xmult
c
      integer nxmax,nxtick,ierr,id1
c
      entry xtxgk1 (xlo,xhi,nxmax,xperti,nxtick,ierr)
c
code ...
c
      ierr = 0
c
      xrange = xhi-xlo
      xperti = xrange/float(nxmax)
c
      if (xperti .gt. zero) then
        id1= nint(alog10(xperti))
        xperti = 10.0**id1
      else if (xperti .lt. zero) then
        id1= nint(alog10(-xperti))
        xperti = -1.0*(10.0**id1)
      else
        ierr = -1
        return
      end if
c
      xmult = 5.0
  123 continue
      nxtick = int (xrange/xperti)
      if (nxtick .gt. nxmax) then
        xperti = xmult*xperti
        xmult = 7.0 - xmult
        goto 123
      end if
c
      if (xtxdbg .gt. 0) write (*,*)
     +  ' x per tick, # ticks : ',xperti,nxtick
c
      return
      end
