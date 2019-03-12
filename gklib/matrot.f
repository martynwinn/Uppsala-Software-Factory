c
c=======================================================================
c
	subroutine matrot (axis, deg, s)
c
c  Matrix_Rotation  --  Given the axis (1,2,3 ~ x,y,z), and the angle 
c  of rotation around this axis, this routine returns the 3x3 
c  rotation matrix. Has code from Alwyns Frodo RRX, RRY and RRZ routines.
c
	implicit none
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
	integer axis
	real deg, s(9)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double precision phi, sa, ca
c
	phi = deg*degtor
	sa = sin(phi)
	ca = cos(phi)
c
	if (axis .eq. 3) then
	   s(1) = ca
	   s(2) = sa
	   s(3) = 0.0
	   s(4) = -sa
	   s(5) = ca
	   s(6) = 0.0
	   s(7) = 0.0
	   s(8) = 0.0
	   s(9) = 1.0
	elseif (axis .eq. 2) then
	   s(1) = ca
	   s(2) = 0.0
	   s(3) = -sa
	   s(4) = 0.0
	   s(5) = 1.0
	   s(6) = 0.0
	   s(7) = sa
	   s(8) = 0.0
	   s(9) = ca
	elseif (axis .eq. 1) then
	   s(1) = 1.0
	   s(2) = 0.0
	   s(3) = 0.0
	   s(4) = 0.0
	   s(5) = ca
	   s(6) = sa
	   s(7) = 0.0
	   s(8) = -sa
	   s(9) = ca
	else
	   call prompt (' matrot: error, strange axis number.')
	endif
c
	return
	end
