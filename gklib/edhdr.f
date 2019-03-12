c
c
c
      subroutine edhdr 
     $  (file, lun, orgn, ext, grid, uvw, cell, spgrp, ed, space)
c
c ---	Read header of a CCP4 map
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
      logical incore
      real rho(maxrho)
c
code ...
c
      incore = .false.
      call prompt (' Read header')
      call maphdr (file, lun, 4, orgn, ext, grid, uvw, cell, spgrp,
     $  rho, maxrho, ed, space, incore)
      call prompt (' Header done')
c
      return
      end
