c
c
c
      real function tangle (i1,i2,i3,i4,xyz)
c
c ... TANGLE - calculate dihedral angle (adapted from PROLSQ)
c
c ... calculate torsion of atoms I1-I2-I3-I4; XYZ(3,*) holds coords
c
      implicit none
c
      real rtodeg
      parameter (rtodeg = 360.0/6.2831853071796)
c
      real xyz(3,*),q(3,3),u(3),v(3),w(3),a(3),b(3)
      real wmag,t,s,det3
c
      integer i1,i2,i3,i4,i,j
c
code ...
c
C----FORM THE TRIAD OF INTERATOMIC VECTORS
c
      DO I=1,3
        U(I) = xyz(I,i1) - xyz(I,i2)
        V(I) = xyz(I,i4) - xyz(I,i3)
        W(I) = xyz(I,i3) - xyz(I,i2)
      end do
C
C----FIND A = U.CROSS.W
c
      A(1) = U(2)*W(3) - U(3)*W(2)
      A(2) = U(3)*W(1) - U(1)*W(3)
      A(3) = U(1)*W(2) - U(2)*W(1)
C
C----FIND B = V.CROSS.W
c
      B(1) = V(2)*W(3) - V(3)*W(2)
      B(2) = V(3)*W(1) - V(1)*W(3)
      B(3) = V(1)*W(2) - V(2)*W(1)
C
C----FIND !W!
c
      WMAG = SQRT(W(1)**2+W(2)**2+W(3)**2)
C
C---- (A.CROSS.B).DOT.W = !A!*!B!*!W!*SIN(ANGLE)
c
      DO J=1,3
        Q(1,J) = A(J)
        Q(2,J) = B(J)
        Q(3,J) = W(J)
      end do
      S = DET3(Q)
C
C---- (A.DOT.B)*!W! = !A!*!B!*!W!*COS(ANGLE)
c
      T = (A(1)*B(1) + A(2)*B(2) + A(3)*B(3))*WMAG
C
C----TORSION ANGLE RESULT IN DEGREES
c
      tangle = ATAN2(S,T) * rtodeg
c
      return
      end
