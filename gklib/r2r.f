c
c
c
      subroutine r2r (intx,realx)
c
c ... "intx" is a real into which an integer has been put
c     convert it back to a real "realx"
c     note: intx and realx may be the same real variable
c
c ... this routine is needed since the mask read routine
c     thinks it gets an integer array to store the mask
c     values in; in fact, it is a real array which becomes
c     incomprehensible when used a real array after the
c     call to the mask read routine.  this routine also
c     assumes that the first argument is integer, whereas
c     in fact it is a real again.
c     this makes Fortran interpret the variable as an integer;
c     we convert it to real and store it in the same variable
c     again
c
      implicit none
c
      integer intx,idum
c
      real realx
c
code ...
c
      idum = intx
      realx = float(idum)
c
      return
      end
