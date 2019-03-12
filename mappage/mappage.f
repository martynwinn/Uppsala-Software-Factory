	program mappage
c
c ---	Author: Alwyn Jones,27-jun-81
c ---	Modifications:
c ---	2/3-Aug-91 unixed by Alwyn Jones
c ---	Too many people have fiddled around with this program, so now
c	I waste some time to clean it up for running on Unix machines.
c ---	Significant (and worthwhile) Modifications:
c	August 1991 Phil Evans, added ccp4 ascii map input, and some
c	little/big endian issues.
c	February 1987 mk (Aarhus) (CCP4 format added)
c	~1983 Bruce Bush added explicit routines to read Protein
c	style maps- they had already been written by W.Steigemann.
c
c --- Spring, 1993 - Gerard Kleywegt - This program is now
c       obsolete !!! Use MAPMAN instead (has more input map options
c       plus facilities to easily set the dynamic range - although
c       I implemented these in this version as well)
c
c ---	NOTE: Some things could not be checked. In particular, Protein
c	style maps do not exist on our Unix machines yet. The FFT-Y
c	style maps will cause the program to halt since they are in the
c	wrong order, and I since I see no use for them I have not
c	bothered doing an in-core save (as done for CCP4). Anyway, I
c	have no example for de-bug purposes.

c ---	If you tinker with this program, you should know the following
c	1.Each density point is stored as a byte in the direct access
c	  file. This means one needs to do a conversion. This is done
c	  in MAKBYT which is computer specific. The version of MAKBYT
c         provided here should work on most (all?) machines
c	2.The random access file OPEN command could depend on your
c	  computer system. See switch section below setting KRECL
c	3.A call is made to a routine called BYTSWP. This is necessary
c	  on Unix machines because I decided to make Unix O read the map
c	  files produced by the Oatly/Cambillau Silicon Graphics Frodo.
c         Byte-swapping of the header record is conditional on the
c         "endedness" of the machine, tested by subroutine LITEND
c	4.If you add a new map format, only MAPHDR and MAPIN need to
c	  be changed. The program must have maps sections where
c	  Y is fastest changing, then X, then Z (i.e. z-sectioned).
c	  If yours is not, check what I do for CCP4 maps.
c	  If you add a new map format, we would appreciate getting a
c	  copy of the changes needed.
c
	implicit none
C-*
c   Switch for changing version: set here
c      krecl, record length for output direct access file
c          = 128 words = 512 bytes
         integer krecl
C         parameter (krecl=128*4)
        parameter (krecl=128)
C-*
	common /page/ ioxyz(3), imxyz(3), iuvw(3), mygrid(3),
     $	plus, nx, ny, nz, cmnd, indxra,  prod, cell(6), rhomax, rhomin,
     $	scale
	integer cmnd, indxra, ioxyz, imxyz, iuvw, mygrid, 
     $	nx, ny, nz, plus
	real cell, prod, rhomax, rhomin, scale

	common /dens/ space, slice(8,90000), rho(90000)
	integer slice, space
	real rho

	common /core/ incore, sizsav, savrho(10 000 000)
	logical incore
	integer sizsav
	real savrho

	character file*80
	integer i, j, level
	integer*2 first(256)
        logical  litend
c
      call gainit ('MAPPAGE','930615/1.0')
c
        call prompt (' ...')
        call prompt (' ... WHY ARE YOU STILL USING MAPPAGE ???')
        call prompt (' ...')
        call prompt (' ... THIS PROGRAM IS TERRIBLY OBSOLETE')
        call prompt (' ...')
        call prompt (' ... USE MAPMAN OR xdlMAPMAN INSTEAD')
        call prompt (' ...')
c
c	write(6,10)
	space = 90000
	sizsav = 10 000 000
	prod = -1.
	call defvol (file, first, krecl)
	nz = 0
	do 100 i=1,imxyz(3)
c ---	Read in one level
	  call mapin (file, 11, cmnd,
     $                ioxyz, imxyz, mygrid, iuvw, cell,
     $                rho, space, savrho, sizsav, incore, level)

	  do 120 j=1,imxyz(iuvw(1))*imxyz(iuvw(2))
120	    rho(j) = rho(j)*scale
	  if (i .eq. 1 .and. prod .le. 0.) then
c ---	    Calculate prod/plus based on first level
	    rhomax = -99999.
	    rhomin = +99999.
	    do 110 j=1,imxyz(iuvw(1))*imxyz(iuvw(2))
	      if (rho(j) .gt. rhomax) rhomax = rho(j)
	      if (rho(j) .lt. rhomin) rhomin = rho(j)
110	    continue
	    prod = 255./(rhomax-rhomin)
c ---	    There is a maximum value of PROD fixed by packing 100*PROD into
c ---	    an integer*2 header
	    prod = min (prod, 32767./100.)
	    plus = -rhomin*prod
	    write (6,20) prod, plus
	    first(16) = 100*prod
	    first(17) = plus
	  end if
	  if (i .eq. 1) then
c ---	  Write out first record
	    indxra = 1
C+ Byte-swap header for little-ended machines only
	    if (litend()) call bytswp (first)
	    write (12, rec=indxra) first
	  end if
100	  call paged (rho, imxyz(1), imxyz(2), i)
	call rest (imxyz(1),imxyz(2),imxyz(3))
c
c ... close map
c
        call mapclo (11)
c
        call gkquit
c
	stop
10	format(/,t10,'*** Mappage program ****')
20	format(/,' Byte = ',F10.3,' * density + ',I10)
	end
c
c
c
	subroutine defvol (file1, first, krecl)
c ---	Define the parameters of our envelope,and write out the first record
	implicit none
	character file1*80
	integer*2 first(256)
        integer krecl
	common /page/ ioxyz(3), imxyz(3), iuvw(3), mygrid(3),
     $	plus, nx, ny, nz, cmnd, indxra,  prod, cell(6), rhomax, rhomin,
     $	scale
	integer cmnd, indxra, ioxyz, imxyz, iuvw, mygrid, 
     $	nx, ny, nz, plus
	real cell, prod, rhomax, rhomin, scale

	common /dens/ space, slice(8,90000), rho(90000)
	integer slice, space
	real rho

	common /core/ incore, sizsav, savrho(10 000 000)
	logical incore
	integer sizsav
	real savrho

        
	character file2*80, option(9)*10, what*10
	integer i, i1, i2, nopt, length
      real ppgk(2),xplus

	data option /'PROTEIN   ', 'FFT-Y     ', 'TENEYCK2  ', 
     $	             'CCP4      ', 'X-PLOR    ', 'EZD       ',
     +               'MASK      ', 'NEWEZD    ', 'BINXPLOR  '/
	data nopt /9/
c
code ...
c
      file1 = ' '
      call textin (' Input map file ?',file1)
      if (length(file1) .lt. 1) call errstp (
     +  'No filename provided')
c
      file2 = ' '
      call textin (' O-style map file ?',file2)
      if (length(file2) .lt. 1) call errstp (
     +  'No filename provided')
c
C+PRE
C  Note that there is a dispute between different computers about how
C  to count record-lengths in unformatted direct access files, either
C  in bytes or words (=4 bytes). The record length here needs to be
C  512 bytes == 128 words. KRECL is set at head of main program
c
c ... 970630 - remove MAXREC qualifier (PDP11 days !)
c
c	open (12, file=file2, status='unknown', form='unformatted',
c     $     access='direct', recl=krecl, maxrec=32767)
	open (12, file=file2, status='unknown', form='unformatted',
     $     access='direct', recl=krecl)
C-PRE
      call asciut (' Supported map types :',nopt,option)
      what = 'CCP4'
      call textin (' Type of your input map ?',what)
      if (length(what) .lt. 1) what = 'CCP4'
      call upcase (what)
c ---	Decide what command was given
	do 100 i=1,nopt
	  if (what .eq. option(i)) then
            cmnd = i
            goto 110
          end if
100	continue
      call errstp ('Unrecognised map type')

110	continue
c ---	Read map header
	call maphdr (file1, 11, cmnd,
     $                     ioxyz, imxyz, mygrid, iuvw, cell, i,
     $                     rho, space, savrho, sizsav, incore)
c
c ... fix for MASKS
c
        if (cmnd .eq. 7) then
          call prompt (' Scale, Prod and Plus fixed for maps')
          scale = 1.0
          prod = 100.0
          plus = 1
          goto 6996
	end if
c
      what = 'N'
      call textin (' Do you want to scale the density ?',what)
      if (length(what) .lt. 1) what = 'N'
      call upcase (what)
c
	if (what(1:1) .eq. 'Y' ) then
c ---	read scale constant
          scale = 1.0
          call fvalin (' Scale factor ?',1,scale)
          if (abs(scale). le. 1.0e-9) scale = 0.01
          call fvalut (' SCALE :',1,scale)
	else
	  scale = 1.
	end if
c
      write (*,'(99(1x,a/))')
     + ' Prod & Plus options:',
     + ' D = don''t care (as it used to be)',
     + ' S = set prod and plus',
     + ' F = figure out from dynamic range',
     + ' NOTE: "PLUS" is an integer number'
      what = 'D'
      call textin (' Option for prod/plus ?',what)
      if (length(what) .lt. 1) what = 'D'
      call upcase (what)
c
	if (what(1:1) .eq. 'D' ) goto 6996
c
        if (what(1:1) .eq. 'F') then
          ppgk(1) = 0.0
          ppgk(2) = 100.0
          call fvalin (' Dynamic density range after scaling ?',2,ppgk)
          call rlohi (ppgk(1),ppgk(2))
          prod = 255.0 / (ppgk(2)-ppgk(1))
          xplus = -ppgk(1)*prod
          ppgk(1) = prod
          ppgk(2) = xplus
        else
          ppgk(1) = 1.0
          ppgk(2) = 0.0
        end if
c
          call fvalin (' Values for PROD and PLUS ?',2,ppgk)
          prod = ppgk(1)
          plus = nint(ppgk(2))
          if (abs(prod). le. 1.0e-9) prod = 0.01
	  if (prod .gt. 327.66) then
	    plus = plus*(327.66/prod)
	    prod = 327.66
	  end if
          call fvalut (' PROD  :',1,prod)
          call ivalut (' PLUS  :',1,plus)
          call fvalut (' SCALE :',1,scale)
          ppgk(1) = -plus/prod
          ppgk (2) = (255.0-plus)/prod
          call fvalut (' Dynamic range (scaled density) :',2,ppgk)
          ppgk(1)=ppgk(1)/scale
          ppgk(2)=ppgk(2)/scale
          call fvalut (' Dynamic range (unscaled)       :',2,ppgk)
c
 6996 continue

c ---	That should define everything for me,
c ---	Check if directions are correct
	if (incore) then
	  iuvw(1) = 2
	  iuvw(2) = 1
	  iuvw(3) = 3
	end if

	if (iuvw(1) .ne. 2) goto 1000
	if (iuvw(2) .ne. 1) goto 1000
	if (iuvw(3) .ne. 3) goto 1000
	

c ---	Fill up our first record
	i1 = 80
	i2 = 100
	do 120 i=1,3
	  first(i  ) = ioxyz(i)
	  first(i+3) = imxyz(i)
	  first(i+6) = mygrid(i)
	  first(i+9) = i1*cell(i)
	  first(i+12)= i1*cell(i+3)
120	continue
	first(16) = i2*prod
	first(17) = plus
	first(18) = i1
	first(19) = i2
	do 130 i=20,256
130	  first(i) = 0
c ---	How many of my funny pages are there?
	nx = imxyz(1)/8
	if(mod(imxyz(1),8) .ge. 1)nx = nx+1
	ny = imxyz(2)/8
	if(mod(imxyz(2),8) .ge. 1)ny = ny+1
	nz = imxyz(3)/8
	if(mod(imxyz(3),8) .ge. 1)nz = nz+1
	write (6,90) nx,ny,nz
90	format(' Number of pages along x , y & z ',3I5)
	return

1000	call errstp ('Input map is NOT in Z-sections')
	stop
	end
c
c
c
	subroutine paged (rholev, ix, iy, ilev)
c ---	A new level of density, store it and if necessary
c	write it out as 3-d non-overlapping boxes of 8*8*8 values
c
	implicit none
	integer ix, iy
	real rholev
	dimension rholev (iy,ix)

	common /page/ ioxyz(3), imxyz(3), iuvw(3), mygrid(3),
     $	plus, nx, ny, nz, cmnd, indxra,  prod, cell(6), rhomax, rhomin,
     $	scale
	integer cmnd, indxra, ioxyz, imxyz, iuvw, mygrid, 
     $	nx, ny, nz, plus
	real cell, prod, rhomax, rhomin, scale

	common /dens/ space, slice(8,90000), rho(90000)
	integer slice, space
	real rho

	integer i, ict, i1, ilev, j, jct, j1, j2, j3, 
     $	k, k1, k2, k3, value
	integer*2 irecrd(256)
	byte record(512)
	equivalence (record, irecrd)

	i1 = mod(ilev,8)
	if(i1 .eq. 0)i1 = 8
	ict = 0
	do 100 i=1,iy
	do 100 j=1,ix
	  ict = ict+1
	  value = rholev(i,j)* prod+ plus
	  if(value .gt. 255) value = 255
	  if(value .lt.   0) value =   0
	  slice(i1,ict) = value
100	 continue
c ---	Now pick out our non-overlapping bricks ?
	if(i1 .ne. 8)return

	value = 0
c ---	Loop over possible y-pages
	do 110 j=1,ny
	  j1 = (j-1)*8+1
	  j2 =   j  *8
c ---     Loop over possible x-pages
	  do 110 k=1,nx
c ---	  Now get our loop parameters
	    k1 = (k-1)*8+1
	    k2 =   k  *8
	    ict = 0
c ---	    Loop over z-levels
	    do 120 i=1,8
c ---	    Loop over y-indeces of current page
	    do 120 j3=j1,j2
	      jct = (j3-1)*ix+k1-1
c ---	      Loop over x-indeces of current page
	      do 120 k3=k1,k2
	        ict = ict+1
	        jct = jct+1
c ---	        If either direction over edge,pack record
	        if(j3 .gt. iy  .or.  k3 .gt. ix)  then
	          call makbyt (record(ict), value)
	        else
	          call makbyt (record(ict), slice(i,jct))
	        end if
120	    continue
	    indxra = indxra+1
	    call bytswp (record)
	    write (12, rec=indxra) irecrd
110	  continue

	return
	end
c
c
c
	subroutine rest (ix, iy, ilev)
c ---	Write out rest of density
c ---	write it out as 3-d elements.
	implicit none
	integer ix, iy, ilev

	common /page/ ioxyz(3), imxyz(3), iuvw(3), mygrid(3),
     $	plus, nx, ny, nz, cmnd, indxra,  prod, cell(6), rhomax, rhomin,
     $	scale
	integer cmnd, indxra, ioxyz, imxyz, iuvw, mygrid, 
     $	nx, ny, nz, plus
	real cell, prod, rhomax, rhomin, scale

	common /dens/ space, slice(8,90000), rho(90000)
	integer slice, space
	real rho

	integer i, ict, i1, j, jct, j1, j2, j3, 
     $	k, k1, k2, k3
	byte record(512)
	integer*2 irecrd(256)
	equivalence (record(1),irecrd(1))
 
	i1 = mod(ilev,8)
	if(i1 .eq. 0)i1= 8
	ict = 0
c ---	Now pick out our overlapping bricks .
c ---	Loop over possible y-pages
	do 100 j=1,ny
	  j1 = (j-1)*8+1
	  j2 =   j  *8
c ---     Loop over possible x-pages
	   do 100 k=1,nx
c ---	   Now get our loop parameters
	     k1 = (k-1)*8+1
	     k2 =   k  *8
	     ict = 0
c ---	     Loop over z-levels
	     do 110 i=1,8
c ---	     Loop over y-indeces of current page
	     do 110 j3 = j1,j2
	       jct = (j3-1)*ix+k1-1
c ---	       Loop over x-indeces of current page
	       do 110 k3 = k1,k2
	         ict = ict+1
	         jct = jct+1
c ---	         If either direction over edge,pack record
c	         or if z - direction over edge,pack record
	         if(j3 .gt. iy  .or.  k3 .gt. ix .or.
     $		                       i .gt. i1 ) then
		   call makbyt (record(ict), plus)
	         else
	           call makbyt (record(ict), slice(i,jct))
	         end if
110	     continue
	     indxra = indxra+1
	     call bytswp (record)
	     write (12, rec=indxra) irecrd
100	continue

	return
	end


c=======================================================================
	subroutine bytswp (rec)
	implicit none
	byte rec(2,256), one
	integer i
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	do 100 i=1,256
	  one = rec (1,i)
	  rec(1,i) = rec(2,i)
100	  rec(2,i) = one
	return
	end

c=======================================================================
	subroutine makbyt (value, ivalue)
c ---	This is a machine specific routine.
c	Depending on whether the machine is big or little endian
c	the last assignment statement will have to be changed.
c ---	On Vax,IBM  PC's need one assignment
c	Other machines use the other
	implicit none
	byte value
	integer ivalue
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	integer i
	byte j(4)
	equivalence (i,j)

C+PRE This should work on all machines
        value = ivalue

c	i = ivalue
c ---	Most unix machines and IBM mainframes
c	value = j(4)
c ---	VAX/IBM PC
c	value = j(1)
	return
	end
c
c
c
      logical function litend ()
C
C Check endedness, return true if little-ended (Vax-like),
C false if big-ended (IBM-like)
C
      integer i
c
      byte b(4)
c
      equivalence (i,b(1))
c
code ...
c
      i=1
      litend = (b(1) .ne. 0)
c
      return
      end
