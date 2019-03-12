c
c
c
	subroutine swpord (a, b, ix, iy)
c
	implicit none
c
	integer ix, iy
	real a(ix, iy), b(iy, ix)
	integer i, j
c
	do 100 i=1,ix
	do 100 j=1,iy
100	  b(j,i) = a(i,j)
	return
	end
