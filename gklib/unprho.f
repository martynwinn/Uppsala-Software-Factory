c
c
c
	subroutine unprho (sav, ix, iy, iz, level, rho, jx, jy, uvw)
c
c ---	Unpack from in-core array to a 213 type level
c
	implicit none
c
	integer ix, iy, iz, level, jx, jy, uvw(3)
	real sav(ix,iy,iz), rho(jx, jy)
	integer i, j
c
code ...
c
        if (uvw(1) .eq. 1 .and. uvw(2) .eq. 2 .and. uvw(3) .eq. 3) then
          do 100 i=1,iy
          do 100 j=1,ix
100         rho(j,i) = sav(j,i,level)
        else if 
     $	(uvw(1) .eq. 1 .and. uvw(2) .eq. 3 .and. uvw(3) .eq. 2) then
          do 110 i=1,iz
          do 110 j=1,ix
110         rho(j,i) = sav(j,level,i)
        else if 
     $	(uvw(1) .eq. 2 .and. uvw(2) .eq. 1 .and. uvw(3) .eq. 3) then
          do 120 i=1,ix
          do 120 j=1,iy
120         rho(j,i) = sav(i,j,level)
        else if 
     $	(uvw(1) .eq. 2 .and. uvw(2) .eq. 3 .and. uvw(3) .eq. 1) then
          do 130 i=1,iz
          do 130 j=1,iy
130         rho(j,i) = sav(level,j,i)
        else if 
     $	(uvw(1) .eq. 3 .and. uvw(2) .eq. 1 .and. uvw(3) .eq. 2) then
          do 140 i=1,ix
          do 140 j=1,iz
140         rho(j,i) = sav(i,level,j)
        else if 
     $	(uvw(1) .eq. 3 .and. uvw(2) .eq. 2 .and. uvw(3) .eq. 1) then
          do 150 i=1,iy
          do 150 j=1,iz
150         rho(j,i) = sav(level,i,j)
        end if
c
	return
	end
