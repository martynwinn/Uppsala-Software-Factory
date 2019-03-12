c
c
c
      real function spline (DX,Y0,Y1,Z0,Z1)
c
      implicit none
c
      real r6
      parameter (r6 = 1.000/6.000)
c
      real DX,Y0,Y1,Z0,Z1,c1,c2,c3,c4,hx
c
code ...
c
      C1 = R6*Z0
      C2 = R6*Z1
      C3 = Y0-C1
      C4 = Y1-C2
      HX = 1.-DX
c
      SPLINE = (C1*HX**2+C3)*HX + (C2*DX**2+C4)*DX
c
      RETURN
      END
