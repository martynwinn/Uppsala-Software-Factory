c
c
c
      subroutine c2cart (nat,xyz,int,ref)
c
c ... C2CART - convert internal to Cartesian coordinates
c
c ... NAT = number of atoms
c     XYZ (3,*) = cartesian coords of atoms
c     INT (3,*) = internal coords of atoms (distance,angle,torsion)
c     REF (3,*) = pointers to atoms with respect to which the
c                 internal coords are defined
c
      implicit none
c
      real twopi,pi,degtor,rtodeg,tiny
c
      parameter (twopi=6.2831853071796)
      parameter (pi=0.50*twopi)
      parameter (degtor=twopi/360.0)
      parameter (rtodeg=360.0/twopi)
      parameter (tiny=1.0E-9)
c
      real xyz(3,*),int(3,*)
c
      integer ref(3,*),nat
c
code ...
c
      if (nat .lt. 1) return
c
      xyz (1,1) = 0.0
      xyz (2,1) = 0.0
      xyz (3,1) = 0.0
      if (nat .lt. 2) return
c
      xyz (1,2) = int(1,2)
      xyz (2,2) = 0.0
      xyz (3,2) = 0.0
c
      if (abs(xyz(1,2)) .lt. tiny) then
        call prompt (
     +    ' C2CART - WARNING - Invalid X-coordinate for atom #2')
        call rvalut (' X #2 :',1,xyz(1,2))
      end if
      if (nat .lt. 3) return
c
      xyz (1,3) = xyz(1,2) - int(1,3)*cos(degtor*int(2,3))
      xyz (2,3) = int(1,3)*sin(degtor*int(2,3))
      xyz (3,3) = 0.0
c
      if (abs(xyz(1,3)) .lt. tiny) then
        call prompt (
     +    ' C2CART - WARNING - Invalid X-coordinate for atom #3')
        call rvalut (' X #3 :',1,xyz(1,3))
      end if
      if (abs(xyz(2,3)) .lt. tiny) then
        call prompt (
     +    ' C2CART - WARNING - Invalid Y-coordinate for atom #3')
        call rvalut (' Y #3 :',1,xyz(2,3))
      end if
      if (abs(xyz(1,3)-xyz(1,2)) .lt. tiny) then
        call prompt (
     +    ' C2CART - WARNING - Identical X-coords atom #2 and #3')
        call rvalut (' X #2 :',1,xyz(1,2))
        call rvalut (' X #3 :',1,xyz(1,3))
      end if
      if (nat .lt. 4) return
c
      call c2car4 (nat,xyz,int,ref)
c
      return
      end
