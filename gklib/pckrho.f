c
c
c
	subroutine pckrho (sav, ix, iy, iz, level, rho, jx, jy, uvw)
c
c          call pckrho (savrho, imxyz(1), imxyz(2), imxyz(3), i,
c     $          rho, imxyz(iuvw(1)), imxyz(iuvw(2)), iuvw)
c
c ---	Pack from a general level into in-core array, orderered (x,y,z)
c
	implicit none
c
	integer ix, iy, iz, level, jx, jy, uvw(3)
	real sav(ix,iy,iz), rho(jx, jy)
	integer i, j
c
code ...
c
cxyz
c        print *,'IX,IY,IZ    - ',ix, iy, iz
c        print *,'JX,JY,LEVEL - ',jx,jy,level
c        print *,'UVW         - ',uvw
c
	if (uvw(1) .eq. 1 .and. uvw(2) .eq. 2 .and. uvw(3) .eq. 3) then
	  do 100 i=1,iy
	  do 100 j=1,ix
100	    sav(j,i,level) = rho(j,i)
	else if (uvw(1) .eq. 1 .and. uvw(2) .eq. 3 .and. uvw(3) .eq. 2) then
	  do 110 i=1,iz
	  do 110 j=1,ix
110	    sav(j,level,i) = rho(j,i)
	else if (uvw(1) .eq. 2 .and. uvw(2) .eq. 1 .and. uvw(3) .eq. 3) then
	  do 120 j=1,iy
	  do 120 i=1,ix
120	    sav(i,j,level) = rho(j,i)
	else if (uvw(1) .eq. 2 .and. uvw(2) .eq. 3 .and. uvw(3) .eq. 1) then
	  do 130 i=1,iz
	  do 130 j=1,iy
130	    sav(level,j,i) = rho(j,i)
	else if (uvw(1) .eq. 3 .and. uvw(2) .eq. 1 .and. uvw(3) .eq. 2) then
	  do 140 i=1,ix
cxyz
c          call ivalut (' I(x) now :',1,i)
	  do 140 j=1,iz
cxyz
c            call ivalut (' J(z) now :',1,j)
c            call rvalut (' RHO(j,i) :',1,rho(j,i))
c            call rvalut (' SAV(I,LEVEL,J) :',1,sav(i,level,j))
140	    sav(i,level,j) = rho(j,i)
	else if (uvw(1) .eq. 3 .and. uvw(2) .eq. 2 .and. uvw(3) .eq. 1) then
	  do 150 i=1,iy
	  do 150 j=1,iz
150	    sav(level,i,j) = rho(j,i)
	end if
c
	return
	end
