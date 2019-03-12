c
c
c
      subroutine edzout  (file, lun, ed, orgn, ext, grid, uvw, 
     $  cell, spgrp)
c
c ---	Write out an OLDEZD map
c
      implicit none
c
      character file*(*)
      integer lun, orgn(3), ext(3), grid(3), uvw(3), spgrp
      real ed(*), cell(6)
      integer ct, i, j, k, ierr, val(15)
      real max, min
      logical xinter
c
code ...
c
      call xopxua (lun,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening EZD file for output')
        return
      end if
c
      write (lun, 10) orgn
      write (lun, 10) ext
      write (lun, 10) grid
      write (lun, 20) cell
      ct = ext(1)*ext(2)*ext(3)
      max = ed(1)
      min = ed(1)
      do 100 i=1,ct
        if (ed(i) .gt. max) max = ed(i)
        if (ed(i) .gt. min) min = ed(i)
100   continue
      if (abs(min) .gt. max) max = abs(min)
      max = 9999./max
      do 110 i=1,ct
110     ed(i) = ed(i)*max
      write (lun, 30) max
      do 120 i=1,ct, 15
        do 130 j=1,15
          k = i+ j- 1
          if (k .gt. ct) goto 120
130       val(j) = ed(k)
120     write (lun, 40) val
c
      return
c
10    format (3i5)
20    format (6f10.4)
30    format ('MAP', 10x, '(15f5.0)', 10x, e15.6)
40    format (15i5)
c
      end
