c
c
c
      subroutine nonpar (d1,d2,n,w1,w2,d,zd,probd,rs,probrs)
c
c ... NONPAR - non-parametric correlation
c
c     Numerical Recipes, pp. 490
c
      implicit none
c
      integer n
      real d1(n),d2(n),w1(n),w2(n)
c
      real erfcc,betai
      real d,zd,probd,rs,probrs,sf,sg,en,en3n,aved,fac,vard,t,df
c
      integer j
c
code ...
c
      do j=1,n
        w1(j)=d1(j)
        w2(j)=d2(j)
      end do
c
      call sort2 (n,w1,w2)
      call crank (n,w1,sf)
      call sort2 (n,w2,w1)
      call crank (n,w2,sg)
c
      d=0.0
      do j=1,n
        d = d + (w1(j)-w2(j))**2
      end do
      en = n
      en3n = en**3 - en
      aved = en3n/6.0 - (sf+sg)/12.0
      fac = (1.0-sf/en3n)*(1.0-sg/en3n)
      vard = ((en-1.0)*en**2*(en+1.0)**2/36.0)*fac
      zd = (d-aved)/sqrt(vard)
      probd = erfcc (abs(zd)/1.4142136)
      rs = (1.0 - (6.0/en3n)*(d+0.5*(sf+sg)))/fac
      t = rs*sqrt((en-2.0)/((1.0+rs)*(1.0-rs)))
      df = en - 2.0
      probrs = betai (0.5*df,0.5,df/(df+t**2))
c
      return
      end
