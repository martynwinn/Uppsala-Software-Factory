c
c ========================================================================
c
      subroutine rlohi (flt1,flt2)
c
      implicit none
c
c === RLOHI (R1,R2) = make R1 the smaller of {R1,R2}
c
c === G J Kleywegt @ 910611
c
      real flt1,flt2,fdum
c
code ...
c
      if (flt1 .gt. flt2) then
        fdum = flt1
        flt1 = flt2
        flt2 = fdum
      end if
c
      return
      end
