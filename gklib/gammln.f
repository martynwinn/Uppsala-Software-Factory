c
c
c
      real function gammln (xx)
c
c ... GAMMLN - Numerical recipes, p. 157
c
      implicit none
c
      real xx
c
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
c
      integer j
c

      data cof,stp/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      data half,one,fpf /0.5D0,1.0D0,5.5D0/
c
code ...
c

      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
c
      return
      end
