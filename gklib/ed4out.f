c
c
c
      subroutine ed4out 
     $  (file, lun, ed, ext1, ext2, ext3, orgn, grid, uvw, cell, spgrp)
c
c ---	Write out a CCP4 map
c ---	This version writes out Y-sections, Z-fastest in each section.
c ---	The space group variable is wrong
c
      implicit none
c
      integer maxrho
      parameter (maxrho=1000000)
c
      character file*(*)
      integer lun, orgn(3), ext1, ext2, ext3, grid(3), uvw(3), spgrp
      real ed(ext1, ext2, ext3), cell(6)
c ---	Work variables
      character myline*128
      integer e(3), i, j, k, mode
      real rho(maxrho), rhoav, rhomin, rhomax, rhosq
c
code ...
c
c ... for MAP, mode = 2
c
      mode = 2
      goto 6206
c
c ... for MASK, mode = 0
c
      entry edxout
     $  (file, lun, ed, ext1, ext2, ext3, orgn, grid, uvw, cell, spgrp)
c
      mode = 0
c
 6206 continue
c
      e(1) = ext1
      e(2) = ext2
      e(3) = ext3
c
      call stamp (myline)
      call textut (' Stamp :',myline)
c
c ---	Header
c
      call mwrhdl (lun,file,myline,e(uvw(3)),uvw,grid,orgn(uvw(3)),
     $  orgn(uvw(1)),(orgn(uvw(1))+e(uvw(1))-1), 
     $  orgn(uvw(2)),(orgn(uvw(2))+e(uvw(2))-1), cell, spgrp, mode)
c
      call mttcpy (myline)
c
c ---	Write out the y-sections, z fastest
c
      do 100 j=1,e(uvw(3))
        call unprho 
     $    (ed, ext1, ext2, ext3, j, rho, e(uvw(1)), e(uvw(2)), uvw)
        call mspew (lun, rho)
100   continue
c
      rhomin = 99999999.
      rhomax = -99999999.
      rhosq = 0.
      rhoav = 0.
      do 200 k=1,ext3
      do 200 j=1,ext2
      do 200 i=1,ext1
        if (rhomin .gt. ed(i,j,k)) rhomin = ed(i,j,k)
        if (rhomax .lt. ed(i,j,k)) rhomax = ed(i,j,k)
        rhosq = rhosq+ ed(i,j,k)*ed(i,j,k)
        rhoav = rhoav+ ed(i,j,k)
200   continue
c
c ---	Close the file
c
      call mclose (lun, rhomin, rhomax, rhoav, rhosq)
      call prompt (' Map written out')
c
      return
      end
