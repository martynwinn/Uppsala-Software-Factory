c
c
c
      subroutine wrmask (lun,mask,x,y,z,origin,grid,cell,how,ierr)
c
c ---   Write out the mask in modified AVGSYS format
c ---   Alwyn Jonews 11-Nov-90
c
      implicit none
c
      real cell(6)
c
      integer lun, x, y, z, grid(3), origin(3)
      integer mask (x,y,z)
      integer i, j, k, length, ierr
c
      character how*(*),line*128
c
code ...
c
      ierr = -1
      call stamp (line)
c
      if (how(1:3) .eq. 'NEW') then
c
        call prompt (' Writing mask (new format)')
        write (lun,6000,err=9000) 'NEW_MASK'
        write (lun,6000,err=9000)
     +    ('! '//line(1:length(line)))
c
        write (line,110,err=9000) 'ORIGIN',origin
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,110,err=9000) 'EXTENT',x,y,z
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,110,err=9000) 'GRID',grid
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,120,err=9000) 'CELL',cell
        write (lun,6000,err=9000) (line(1:length(line)))
c
        write (lun,6000,err=9000) 'MAP'
        write (lun,130,err=9000)
     +    (((mask(i,j,k),i=1,x), j=1,y), k=1,z)
c
      else if (how(1:3) .eq. 'OMA') then
c
        call prompt (' Writing mask (compressed O format)')
        write (lun,6000,err=9000) '.MASK_INPUT'
        write (lun,6000,err=9000)
     +    ('! '//line(1:length(line)))
c
        write (line,110,err=9000) 'ORIGIN',origin
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,110,err=9000) 'EXTENT',x,y,z
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,110,err=9000) 'GRID',grid
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,120,err=9000) 'CELL',cell
        write (lun,6000,err=9000) (line(1:length(line)))
c
        write (lun,6000,err=9000) 'COMPRESSED'
        call wrmcom (mask,(x*y*z),lun,ierr)
        write (lun,6000,err=9000) 'END'
        if (ierr .ne. 0) goto 9000
c
      else if (how(1:3) .eq. 'OEX') then
c
        call prompt (' Writing mask (explicit O format)')
        write (lun,6000,err=9000) '.MASK_INPUT'
        write (lun,6000,err=9000)
     +    ('! '//line(1:length(line)))
c
        write (line,110,err=9000) 'ORIGIN',origin
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,110,err=9000) 'EXTENT',x,y,z
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,110,err=9000) 'GRID',grid
        write (lun,6000,err=9000) (line(1:length(line)))
        write (line,120,err=9000) 'CELL',cell
        write (lun,6000,err=9000) (line(1:length(line)))
c
        write (lun,6000,err=9000) 'EXPLICIT'
        write (lun,130,err=9000)
     +    (((mask(i,j,k),i=1,x), j=1,y), k=1,z)
c
      else if (how(1:3) .eq. 'COM') then
c
        call prompt (' Writing mask (compressed format)')
        write (lun,6000,err=9000) 'COMPRESSED_MASK'
        write (lun, 10,err=9000) origin,(line(1:length(line)))
        write (lun, 10,err=9000) x, y, z
        write (lun, 10,err=9000) grid
        write (lun, 20,err=9000) cell
        call wrmcom (mask,(x*y*z),lun,ierr)
        if (ierr .ne. 0) goto 9000
c
      else
c
c ... default - OLD format
c
        call prompt (' Writing mask (old format)')
        write (lun, 10,err=9000) origin,(line(1:length(line)))
        write (lun, 10,err=9000) x, y, z
        write (lun, 10,err=9000) grid
        write (lun, 20,err=9000) cell
        write (lun, 30,err=9000)
     +    (((mask(i,j,k),i=1,x), j=1,y), k=1,z)
c
      end if
c
      call prompt (' Mask write OK')
      ierr = 0
c
      return
c
c ... WRITE errors
c
 9000 continue
      call errcon ('While writing MASK file')
      return
c
10    format (3i5,1x,a)
20    format (6f10.4)
30    format (40i2)
c
110   format (a,3(1x,i10))
120   format (a,6(1x,f10.4))
130   format (80i1)
c
6000  format (a)
c
      end
