c
c ===========================================================================
c
      subroutine xplpol (ang,mat)
c
c ... convert X-PLOR (Rossmann 1962) polar angles PHI, PSI, KAPPA
c     into a rotation matrix
c
      implicit none
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real ang(3),mat(3,3)
c
      double precision phi,psi,kap,ck,sk,cf,sf,cp,sp,cck
c
code ...
c
      phi = ang(1)*degtor
      psi = ang(2)*degtor
      kap = ang(3)*degtor
c
      ck = cos(kap)
      sk = sin(kap)
      cf = cos(phi)
      sf = sin(phi)
      cp = cos(psi)
      sp = sin(psi)
c
      cck = 1.0D0 - ck
c
      mat (1,1) = ck+cck*sp*sp*cf*cf
      mat (2,1) = sp*(cck*cp*cf  + sk*sf)
      mat (3,1) = -cck*sp*sp*cf*sf + sk*cp
c
      mat (1,2) = sp*(cck*cp*cf - sk*sf)
      mat (2,2) = ck + cck*cp*cp
      mat (3,2) = -sp*(cck*cp*sf + sk*cf)
c
      mat (1,3) = -cck*sp*sp*cf*sf - sk*cp
      mat (2,3) = -sp*(cck*cp*sf - sk*cf)
      mat (3,3) = ck+cck*sp*sp*sf*sf
c
      return
      end
