c
c
c
      subroutine linint (a, sizx, sizy, sizz, x, value, errcod)
c
c --- linear interpolation using 8 nearest grid points
c
      implicit none
c
      integer errcod, sizx, sizy, sizz
c
      real a(sizx, sizy, sizz), x(3), value
      real t,u,v,mt,mu,mv
c
      integer i0,j0,k0
c
code ...
c
      errcod = 1
c
      i0 = int(x(1))
      j0 = int(x(2))
      k0 = int(x(3))
c
      if (i0 .le. 0) goto 10
      if (j0 .le. 0) goto 10
      if (k0 .le. 0) goto 10
      if (i0 .ge. sizx) goto 10
      if (j0 .ge. sizy) goto 10
      if (k0 .ge. sizz) goto 10
c
      t = x(1)-float(i0)
      u = x(2)-float(j0)
      v = x(3)-float(k0)
      mt = 1.0-t
      mu = 1.0-u
      mv = 1.0-v
c
      value = mu*mt*mv*a(i0,j0,k0)   + t*u*v*a(i0+1,j0+1,k0+1) +
     +        t*mu*mv*a(i0+1,j0,k0)  + mt*u*mv*a(i0,j0+1,k0) +
     +        mt*mu*v*a(i0,j0,k0+1)  + t*mu*v*a(i0+1,j0,k0+1) +
     +        t*u*mv*a(i0+1,j0+1,k0) + mt*u*v*a(i0,j0+1,k0+1)
c
      errcod = 0
c
   10 continue
c
      return
      end
