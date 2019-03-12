c
c ===========================================================================
c
      subroutine ccppol(ang,amat)
c
c ... Rams' routine
c
c... 
c   In the spherical polar convention angles Omega and Phi
c    define the lattitude (relative to polar axis z) and longi-
c    tude (from x) respectively of the rotation axis. The angle
c    of rotation about this axis is Chi.
c
c    The rotation matrix is defined as follows:
c
c    | l**2+(m**2+n**2)cosChi   lm(1-cosChi)-nsinChi     nl(1-cosChi)+msinChi |
c    | lm(1-cosChi)+nsinChi     m**2+(l**2+n**2)cosChi   mn(1-cosChi)-lsinChi |
c    | nl(1-cosChi)-msinChi     mn(1-cosChi)+lsinChi    n**2+(l**2+m**2)cosChi |
c
c    where l m n are the direction cosines of the axis about which the rotation
c    through Chi takes place.
c
c                   | l |        | sinOmega*cosPhi |
c                   | m |   =    | sinOmega*sinPhi |
c                   | n |        | cosOmega        |
c.    that part is pinched again from the ccp4 overview .doc
c...
c
      implicit none
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real ang(3),amat(3,3)
      double precision coschi,sinchi,al,am,an,omega,phi,chi,done
c
code ...
c
      done = 1.0d0
      omega = ang(1)*degtor
      phi = ang(2)*degtor
      chi = ang(3)*degtor
c
      al = sin(Omega)*cos(Phi)
      am = sin(Omega)*sin(Phi)
      an = cos(Omega)
c
      coschi=cos(chi)
      sinchi=sin(chi)
c
      amat(1,1)=al**2+(am**2+an**2)*cosChi
      amat(1,2)=al*am*(done-cosChi)-an*sinChi 
      amat(1,3)=an*al*(done-cosChi)+am*sinChi
      amat(2,1)=al*am*(done-cosChi)+an*sinChi
      amat(2,2)=am**2+(al**2+an**2)*cosChi
      amat(2,3)=am*an*(done-cosChi)-al*sinChi
      amat(3,1)=an*al*(done-cosChi)-am*sinChi
      amat(3,2)=am*an*(done-cosChi)+al*sinChi
      amat(3,3)=an**2+(al**2+am**2)*cosChi
c
      return
      end
