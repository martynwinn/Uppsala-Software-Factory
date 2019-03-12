c
c
c
      subroutine edinx
     $  (file, lun, ed, orgn, ext, grid, uvw, cell, spgrp, space)
c
c ---	Read in an X-plor map
c
      implicit none
c
      integer maxrho
      parameter (maxrho=1000000)
c
      character file*(*)
      integer lun, orgn(3), ext(3), grid(3), uvw(3), space, spgrp
      real ed(*), cell(6)
c ---	Work variables
      integer i, j
      logical incore
      real rho(maxrho),savrho(maxrho),ave,sdv
c
code ...
c
      incore = .false.
      call prompt (' Reading Header')
      call maphdr (file, lun, 5, orgn, ext, grid, uvw, cell, spgrp,
     $  rho, maxrho, ed, space, incore)
      call prompt (' Header done')
c
      if (.not. incore) then
        call prompt (' Not in core')
        j = 0
        call prompt (' Reading a level')
        do 100 i=1, ext(uvw(3))
          call mapin (file, lun, 5, orgn, ext, grid, uvw, cell, 
     $      rho, maxrho, savrho, 1, incore, j)
          call pckrho (ed, ext(1), ext(2), ext(3), i, 
     $      rho, ext(uvw(1)), ext(uvw(2)), uvw)
100     continue
      end if
c
c ... 951113 - check if it's an X-PLOR4 map
c
      ave = 0.0
      sdv = -1.0
c
      read (lun,'(i8)',err=200,end=200) i
      if (i .ne. -9999) goto 200
      read (lun,'(1x,2e12.5)',err=200,end=200) ave,sdv
      if (ave .eq. 0.0 .or. sdv .le. 0.0) goto 200
c
      call prompt (' This is an X-PLOR version 4 map !')
      call rvalut (' Average density :',1,ave)
      call rvalut (' Sigma level     :',1,sdv)
      call prompt (' Restoring original density values (e-/A3)')
      do i=1,ext(1)*ext(2)*ext(3)
        ed(i) = (sdv * ed(i)) + ave
      end do
c
  200 continue
      call prompt (' Map read OK')
c
c ... close map
c
      call mapclo (lun)
c
      return
      end
