c
c
c
      subroutine edout  (file, lun, ed, orgn, ext, grid, uvw, 
     $                   cell, spgrp)
c
c ---	Write out a CCP4 map
c
      implicit none
c
      character file*(*)
      integer lun, orgn(3), ext(3), grid(3), uvw(3), spgrp
      real ed(*), cell(6)
c
code ...
c
      call ed4out 
     $ (file, lun, ed, ext(1), ext(2), ext(3), orgn, grid, uvw, 
     $  cell, spgrp)
c
      return
      end
