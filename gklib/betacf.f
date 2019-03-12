c
c
c
      real function betacf (a,b,x)
c
c ... BETACF - Numerical Recipes, p. 168
c
      implicit none
c
      integer itmax
      parameter (itmax=2000)
c
      real eps
      parameter (eps=3.0E-7)
c
      real a,b,x,am,bm,az,qab,qap,qam,bz,em,tem,d,ap,bp
      real app,bpp,aold
c
      integer m
c
code ...
c
      am = 1.0
      bm = 1.0
      az = 1.0
      qab = a+b
      qap = a+1.0
      qam = a-1.0
      bz = 1.0 - qab*x/qap
      do m=1,itmax
        em = m
        tem = em + em
        d = em*(b-m)*x/((qam+tem)*(a+tem))
        ap = az + d*am
        bp = bz + d*bm
        d = - (a+em)*(qab+em)*x/((a+tem)*(qap+tem))
        app = ap + d*az
        bpp = bp + d*bz
        aold = az
        am = ap/bpp
        bm = bp/bpp
        az = app/bpp
        bz = 1.0
        if (abs(az-aold) .lt. eps*abs(az)) goto 1
      end do
      call errcon (
     +  'BETACF - argument too big or ITMAX too small')
      betacf = 0.0
      return
c
    1 continue
      betacf = az
c
      return
      end
