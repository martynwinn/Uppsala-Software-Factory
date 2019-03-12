c
c
c
      subroutine maskin (lun,mask,orgn,ext,grid,cell,space,ierr)
c
c --- read a mask
c
      implicit none
c
      integer lun, orgn(3), ext(3), grid(3), space, ierr, length
      integer mask (1), idum
c
      real cell(6)
c
      character line*80
c
code ...
c
      ierr = -1
c
c ... OLD FORMAT
c
      read (lun, 10,err=100,end=9000) orgn
      call prompt (' Reading mask (old format)')
      read (lun, 10,err=9000,end=9000) ext
      read (lun, 10,err=9000,end=9000) grid
      read (lun, 20,err=9000,end=9000) cell
c
      call maskqq (lun,mask,ext(1),ext(2),ext(3),space,'OLD',ierr)
      call maskst (mask,ext(1),ext(2),ext(3))
      goto 900
c
  100 continue
      rewind (lun)
  105 read (lun,6000,err=9000,end=9000) line
      if (length(line).lt.1) goto 105
      if (line(1:1) .eq. '!') goto 105
c
      if (line(1:15) .eq. 'COMPRESSED_MASK') goto 200
      rewind (lun)
      if (line(1:8)  .eq. 'NEW_MASK') goto 110
      if (line(1:11) .eq. '.MASK_INPUT') goto 300
c
      call errcon ('Unknown mask format !!!')
      call textut (' Header :',line)
      call prompt (' Assuming O mask format')
      goto 300
c
c ... NEW FORMAT
c
  110 continue
      read (lun,6000,err=9000,end=9000) line
      if (length(line) .lt. 1) goto 110
      call upcase (line)
      if (line(1:1) .eq. '!') goto 110
c
      if (line(1:6) .eq. 'ORIGIN') then
        read (line(7:),*,err=9000,end=9000) orgn
      else if (line(1:6) .eq. 'EXTENT') then
        read (line(7:),*,err=9000,end=9000) ext
      else if (line(1:4) .eq. 'GRID') then
        read (line(5:),*,err=9000,end=9000) grid
      else if (line(1:4) .eq. 'CELL') then
        read (line(5:),*,err=9000,end=9000) cell
      else if (line(1:3) .eq. 'MAP') then
        goto 120
      else if (line(1:8) .eq. 'NEW_MASK') then
        call prompt (' Reading mask (new format)')
      end if
c
      goto 110
c
  120 continue
c
      call maskqq (lun,mask,ext(1),ext(2),ext(3),space,'NEW',ierr)
      call maskst (mask,ext(1),ext(2),ext(3))
      goto 900
c
c ... COMPRESSED FORMAT
c
  200 continue
      call prompt (' Reading mask (compressed format)')
      read (lun, 10,err=9000,end=9000) orgn
      read (lun, 10,err=9000,end=9000) ext
      read (lun, 10,err=9000,end=9000) grid
      read (lun, 20,err=9000,end=9000) cell
c
      call rdmcom (mask,(ext(1)*ext(2)*ext(3)),space,lun,ierr)
      if (ierr .ne. 0) goto 9000
      call maskst (mask,ext(1),ext(2),ext(3))
      goto 900
c
c ... O_MASK FORMAT
c
  300 continue
      read (lun,6000,err=9000,end=9000) line
      if (length(line) .lt. 1) goto 300
      call upcase (line)
      call pretty (line)
      if (line(1:1) .eq. '!') goto 300
c
      if (line(1:6) .eq. 'ORIGIN') then
        read (line(7:),*,err=9000,end=9000) orgn
      else if (line(1:6) .eq. 'EXTENT') then
        read (line(7:),*,err=9000,end=9000) ext
      else if (line(1:4) .eq. 'GRID') then
        read (line(5:),*,err=9000,end=9000) grid
      else if (line(1:4) .eq. 'CELL') then
        read (line(5:),*,err=9000,end=9000) cell
      else if (line(1:8) .eq. 'EXPLICIT') then
        call prompt (' Format : explicit')
        goto 330
      else if (line(1:10) .eq. 'COMPRESSED') then
        call prompt (' Format : compressed')
        goto 320
      else if (line(1:11) .eq. '.MASK_INPUT') then
        call prompt (' Reading mask (O format)')
      else if (line(1:5) .eq. 'CHECK') then
        read (line(6:),*,err=9000,end=9000) idum
        call jvalut (' Nr of mask points when written :',1,idum)
      else
        call errcon ('Unrecognised line in mask file')
        call textut (' Line :',line)
      end if
c
      goto 300
c
c ... COMPRESSED
c
  320 continue
c
      call rdmcom (mask,(ext(1)*ext(2)*ext(3)),space,lun,ierr)
      if (ierr .ne. 0) goto 9000
      call maskst (mask,ext(1),ext(2),ext(3))
      goto 900
c
c ... EXPLICIT
c
  330 continue
c
      call maskqq (lun,mask,ext(1),ext(2),ext(3),space,'NEW',ierr)
      if (ierr .ne. 0) goto 9000
      call maskst (mask,ext(1),ext(2),ext(3))
      goto 900
c
c ... ROOM FOR MORE FORMATS
c
  900 continue
      ierr = 0
c
      return
c
 9000 continue
      call errcon ('While reading mask')
      return
c
10    format (3i5)
20    format (6f10.4)
c
6000  format (a)
c
      end
