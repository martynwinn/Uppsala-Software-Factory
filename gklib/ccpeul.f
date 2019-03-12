c
c ===========================================================================
c
      subroutine ccpeul(ang,amat)
c
c ... Rams' routine
c
c.+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c.    calculate rotaion matrix from alpha beta gamma output of ccp4.
c.    input alpha, beta and gamma
c.     output  matrix mat(3,3)
c..     alpha beta and gamma are the eulerian angles of ccp4
c.    note: the matrix has been taken from ccp4 documnetation
c.  to read alpha beta and gamma of the eulierian angle output from almn and rotate them to
c.. generate a new pdb file
c..  the folloiwng is the matrix used for this purpose (Taken from the ccp4 manual) - overview.doc
c.  A= alpha, B= beta and  G=gamma
c        |cosAcosBcosG-sinAsinG  -cosAcosBsinG-sinAcosG  cosAsinB|
c        |sinAcosBcosG+cosAsinG  -sinAcosBsinG+cosAcosG  sinAsinB|
c        |      -sinBcosG               sinBsinG           cosB  |
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real ang(3),amat(3,3)
      double precision al,be,ga,cosa,cosb,cosg,sina,sinb,sing
c
code ...
c	
      al = ang(1)*degtor
      be = ang(2)*degtor
      ga = ang(3)*degtor
c
      cosa=cos(al)
      sina=sin(al)
      cosb=cos(be)
      sinb=sin(be)
      cosg=cos(ga)
      sing=sin(ga)
c
      amat(1,1)= cosa*cosb*cosg - sina*sing
      amat(1,2)= -(cosa*cosb*sing) - sina*cosg
      amat(1,3)= cosa*sinb
      amat(2,1)= sina*cosb*cosg + cosa*sing
      amat(2,2)= -(sina*cosb*sing) + cosa*cosg
      amat(2,3)= sina*sinb
      amat(3,1)= -sinb*cosg
      amat(3,2)= sinb*sing
      amat(3,3)= cosb		
c
      return
      end
