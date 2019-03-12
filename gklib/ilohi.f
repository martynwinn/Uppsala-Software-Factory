c
c ========================================================================
c
      subroutine ilohi (int1,int2)
c
      implicit none
c
c === ILOHI (I1,I2) = make I1 the smaller of {I1,I2}
c
c === G J Kleywegt @ 910611
c
      integer idum,int1,int2
c
code ...
c
      if (int1 .gt. int2) then
        idum = int1
        int1 = int2
        int2 = idum
      end if
c
      return
      end
