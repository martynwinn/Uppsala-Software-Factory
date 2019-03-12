c
c ========================================================================
c
      subroutine rswap (r1,r2)
c
      implicit none
c
c === RSWAP (R1,R2) = swap two reals
c
c === G J Kleywegt @ 910911
c
      real r1,r2,rdum
c
code ...
c
      rdum = r1
      r1 = r2
      r2 = rdum
c
      return
      end
