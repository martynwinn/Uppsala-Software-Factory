c
c
c
      subroutine edin 
     $  (file, lun, ed, orgn, ext, grid, uvw, cell, spgrp, space)
c
c ---	Read in a CCP4 map
c
      implicit none
c
      integer maxrho
      parameter (maxrho=1000000)
c
      character file*(*)
      integer lun, orgn(3), ext(3), grid(3), uvw(3), space, spgrp
      real ed(space), cell(6)
c ---	Work variables
      integer i, j, iccp4
      logical incore
      real rho(maxrho)
c
code ...
c
      incore = .false.
      call prompt (' Read header')
      iccp4=4
c
      call maphdr (file, lun, iccp4, orgn, ext, grid, uvw,
     +  cell, spgrp, rho, maxrho, ed, space, incore)
      call prompt (' Header done')
      if (.not. incore) then
        call prompt (' Not in core')
        j = 0
        call prompt (' Reading levels')
        do 100 i=1, ext(uvw(3))
          if (i .eq. i/10*10) 
     $      call ivalut (' Level number :', 1, i)
c     $	    call prmpti ('('' Level number'',i5)', 1, i)
          call mapin (file, lun, 4, orgn, ext, grid, uvw, cell, 
     +      rho, maxrho, rho, 1, incore, j)
          call pckrho (ed, ext(1), ext(2), ext(3), i, 
     $      rho, ext(uvw(1)), ext(uvw(2)), uvw)
100     continue
      end if
c
      call prompt (' Map read OK')
c
c ... close map
c
      call mapclo (lun)
c
      return
      end
