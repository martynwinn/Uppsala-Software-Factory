      program mapfix
c
c ---	alter spacegroup and or title of a CCP4 map
c
c ...   Gerard Kleywegt @ 930507
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'MAPFIX', vers = '040701/3.0.1')
c
      integer maxsiz
      parameter (maxsiz = maxgk1)
c
c      pointer (iaptr,mapa)
c
c      real mapa(1)
c      integer malloc
c
#ifdef ALPHA
      integer*8 iaptr
      integer*8 fmalloc
#else
      integer iaptr
      integer fmalloc
#endif
c
      integer nb, mapsize
c
code ...
c
      call gainit (prognm,vers)
c
      mapsize = maxsiz
      call extint ('MAPSIZE',mapsize)
      mapsize = max (mapsize, minsiz)
      call jvalut (' Allocate maps of size  :',1,mapsize)
c
c ... WRDBYT accounts for 4 or 8 bytes per word
c
      nb = wrdbyt*mapsize
      iaptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0) then
#else
      if (iaptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call domapf (%val(iaptr), mapsize)
c
      call ffree (iaptr)
c
   99 continue
      call gkquit ()
c
      end
c
c
c
      subroutine domapf (mapa, maxsiz)
c
c ---	alter spacegroup and or title of a CCP4 map
c
c ...   Gerard Kleywegt @ 930507
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz
c      parameter (maxsiz = maxgk1)
c
c ---	Map A data structures
      integer orgna(3), exta(3), grida(3), uvwa(3)
      real mapa(maxsiz), cella(6)
c
      character file*80,newmap*80
c
      integer spgrp,length
c
code ...
c
c      call gainit (prognm,vers)
c
c      call jvalut (' Max size of map :',1,maxsiz)
c
c ... 950118 - check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
      file = ' '
      call textin (' Input map ?',file)
      call textut (' Input map :',file)
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
c ---	read in map a
      call edin (file, 1, mapa,  orgna, exta, grida, uvwa,
     $           cella, spgrp, maxsiz)
      call telmap (grida,orgna,exta,cella)
c
      newmap = 'out.E'
      write (*,*)
      call textin (' Output map ?',newmap)
      call textut (' Output map :',newmap)
      if (length(newmap) .gt. 0) then
        call ccpout 
     $    (newmap, 3, mapa, exta(1), exta(2), exta(3),
     $     orgna, grida, uvwa, cella, spgrp)
      else
        call prompt (' No filename; no new map')
      end if
c
      close (1)
      close (3)
c
      return
c
      end
c
c
c
      subroutine ccpout
     $  (file, lun, ed, ext1, ext2, ext3, orgn, grid, uvw, cell, spgrp)
c
c ---   Write out a CCP4 map
c
      implicit none
c
      integer maxrho
      parameter (maxrho=512000)
c
      character file*(*)
      integer lun, orgn(3), ext1, ext2, ext3, grid(3), uvw(3), spgrp
      real ed(ext1, ext2, ext3), cell(6)
c ---   Work variables
      character title*80,file2*80
      integer e(3), i, j, k, ierr
      real rho(maxrho), rhoav, rhomin, rhomax, rhosq
c
code ...
c
      e(1) = ext1
      e(2) = ext2
      e(3) = ext3
      title = 'RAVE/CCP4 map (symm-ops added by MAPFIX)'
      write (*,*)
      call textin (' Title ?',title)
      call textut (' Title :',title)
      call mttcpy (title)
c
      write (*,*)
      call ivalin (' UVW (write-order axes) ?',3,uvw)
      call ivalut (' UVW (write-order axes) :',3,uvw)
c
      write (*,*)
      call ivalin (' Spacegroup ?',1,spgrp)
      call ivalut (' Spacegroup :',1,spgrp)
      spgrp = max (1,spgrp)
c
c ... get CCP4 environment variable CLIBD
c
      call gknval ('CLIBD',file2,ierr)
      if (ierr .ne. 0) then
        file2 = '/nfs/public/packages/ccp4/lib/data/symop.lib'
      else
        call appstr (file2,'/symop.lib')
      end if
      call textin (' Library file with symm-ops ?',file2)
      call textut (' Library file with symm-ops :',file2)
      call xopxoa (12,file2,.true.,j)
      if (j.ne.0) then
        call errcon ('While opening file')
        return
      end if
      write (*,*)
c
c --- Header
c
      call mwrhdl (lun, file, title, e(uvw(3)), uvw, grid, orgn(uvw(3)),
     $  orgn(uvw(1)), (orgn(uvw(1))+e(uvw(1))-1),
     $  orgn(uvw(2)), (orgn(uvw(2))+e(uvw(2))-1), cell, spgrp, 2)
c
      call mttcpy (title)
      call msyput2 (12,spgrp,lun)
c
c ---   Write out the y-sections, z fastest
c
      do 100 j=1,e(uvw(3))
        call unprho
     $    (ed, ext1, ext2, ext3, j, rho, e(uvw(1)), e(uvw(2)), uvw)
        call mspew (lun, rho)
100   continue
c
      call prompt (' Calculating min, max etc.')
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
c ---   Close the file
c
      call mclose (lun, rhomin, rhomax, rhoav, rhosq)
      call prompt (' Map written out')
c
      return
      end
c
c
c
      SUBROUTINE MSYPUT2(IST,LSPGRP,IUNIT)
C     ===================================
C
C---- Read symmetry operator file from stream IST, find entry for
C     space-group LSPGRP. Copy symmetry operators to map stream
C     IUNIT, leaving space at head of file for NBHDR items of
C     header record. Puts number of characters of symmetry
C     information NSYMBT into header record in com  MOHDR.
C
C
C
C---- Map header common
C
C     .. Parameters ..
      INTEGER LUNOUT
      PARAMETER (LUNOUT=6)
C     ..
C     .. Scalar Arguments ..
      INTEGER IST,IUNIT,LSPGRP
C     ..
C     .. Scalars in Common ..
      REAL AMAX,AMEAN,AMIN,ARMS
      INTEGER ISPG,ITMHDR,ITMSC1,LSKFLG,MODE,NC,NC1,NCHITM,NLAB,NR,NR1,
     +        NS,NS1,NSYMBT
C     ..
C     .. Arrays in Common ..
      REAL CEL,SKWMAT,SKWTRN
      INTEGER JUNK,LABELS,LSTRM,MAPCRS,NXYZ
C     ..
C     .. Local Scalars ..
      INTEGER I,IFAIL,ISG,NBHDR,NBLIN,NCLIN,NLIN
C     ..
C     .. Local Arrays ..
      INTEGER JLINE(20)
C     ..
C     .. External Functions ..
      INTEGER NBYTXX
      EXTERNAL NBYTXX
C     ..
C     .. External Subroutines ..
      EXTERNAL CCPOPN,QSEEK,QWRITE
C     ..
C     .. Common blocks ..
      COMMON /MOHDR/NC,NR,NS,MODE,NC1,NR1,NS1,NXYZ(3),CEL(6),MAPCRS(3),
     +       AMIN,AMAX,AMEAN,ISPG,NSYMBT,LSKFLG,SKWMAT(3,3),SKWTRN(3),
     +       JUNK(17),ARMS,NLAB,LABELS(20,10),NCHITM,ITMHDR,ITMSC1
      COMMON /MSTRM/LSTRM(12)
C     ..
C     .. Save statement ..
      SAVE /MSTRM/,/MOHDR/
C     ..
C     .. Data statements ..
      DATA NBHDR/256/
C     ..
C
      NCLIN = NBYTXX(20)
C
C---- Open symmetry file
C
      IFAIL = 1
cc      CALL CCPOPN(IST,'SYMOP',5,1,LDUM,IFAIL)
      IF (IFAIL.LT.0) THEN
C
C---- Error conditions
C
        WRITE (LUNOUT,FMT=6002)
        WRITE (LUNOUT,FMT=6004)
      ELSE
C
C---- Position map file to before symmetry operators
C
        CALL QSEEK(LSTRM(IUNIT),2,1,ITMHDR)
C
C---- Calculate number of items
C     / line (allowing for number of characters /
C
        NBLIN = (NCLIN+NCHITM-1)/NCHITM
   10   CONTINUE
C
C---- Find correct space-group in file.
C     Each space-group has header line of space-group number,
C     number of line of symmetry operations
C
        READ (IST,FMT=*,END=30) ISG,NLIN
        IF (ISG.EQ.LSPGRP) THEN
          GO TO 40
        ELSE
C
C---- Skip NLIN lines
C
          DO 20 I = 1,NLIN
            READ (IST,FMT=*)
   20     CONTINUE
          GO TO 10
        END IF
   30   WRITE (LUNOUT,FMT=6002)
        WRITE (LUNOUT,FMT=6006) LSPGRP
        GO TO 60
C
C---- Space-group found, copy NLIN lines of symmetry
C     operators (NCLIN characters / line) to output file
C
   40   CONTINUE
        DO 50 I = 1,NLIN
          READ (IST,FMT=6000) JLINE
          CALL QWRITE(LSTRM(IUNIT),JLINE,NBLIN)
   50   CONTINUE
C
C---- Number of characters of symmetry information
C
        NSYMBT = NLIN*NCLIN
C
C---- Position of first section
C
        ITMSC1 = NSYMBT/NCHITM + ITMHDR + 1
C
        REWIND IST
        RETURN
      END IF
c
   60 continue
cc      WRITE (LUNOUT,FMT=6008)
      call errcon ('In symmetry file')
cc      CALL CCPERR(1, '**SYMMETRY FILE ERROR** in maplib.for')
C
C---- Format statements
C
 6000 FORMAT (20A4)
 6002 FORMAT (/' **SYMMETRY FILE ERROR**')
 6004 FORMAT (/' **MSYPUT2: ERROR IN OPENING SYMOP FILE**')
 6006 FORMAT (/' **MSYPUT2: NO SYMMETRY INFORMATION FOR SPACE GROUP NUM',
     +       'BER',I4,' IN SYMOP FILE**')
 6008 FORMAT (/' **PROGRAM TERMINATED**')
C
C
      END
