c
c ========================================================================
c
      subroutine iswap (int1,int2)
c
      implicit none
c
c === ISWAP (I1,I2) = swap two integers
c
c === G J Kleywegt @ 910911
c
      integer int1,int2,idum
c
code ...
c
      idum = int1
      int1 = int2
      int2 = idum
c
      return
      end
