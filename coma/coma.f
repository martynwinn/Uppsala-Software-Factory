      program coma
c
c ... COrrelation MAsk
c
c ... Gerard Kleywegt @ 971202
c
c ... Version 0.6 @ 980710
c
c ... f77_6d_comp coma.f
c     f77_6d_link -o COMA coma.o ../gklib/6d_kleylib /home/gerard/lib/libccp4_irix62_32.a
c     strip COMA
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'COMA', vers = '080625/1.0.2')
c
      integer maxsiz, maxmsk
      parameter (maxsiz = maxgk1)
      parameter (maxmsk = maxgk1)
c
c      pointer (iaptr,mapin)
c      pointer (ibptr,mapout)
c      pointer (icptr,mapdum)
c      pointer (idptr,mapi)
c      pointer (ieptr,maski)
c      pointer (igptr,mapj)
c
c      real mapin(1),mapout(1),mapdum(1),mapi(1),mapj(1)
c      integer maski(1), malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr,idptr,ieptr,igptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr,idptr,ieptr,igptr
      integer fmalloc
#endif
c
      integer nb, mapsize, masksize, minmsk
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
      masksize = maxmsk
      call extint ('MASKSIZE',masksize)
      masksize = max (masksize, minsiz)
      call jvalut (' Allocate masks of size :',1,masksize)
c
c ... WRDBYT accounts for 4 or 8 bytes per word
c
      nb = wrdbyt*mapsize
      iaptr = fmalloc (nb)
      ibptr = fmalloc (nb)
      icptr = fmalloc (nb)
c
      minmsk = masksize / 10
      nb = wrdbyt*minmsk
      call jvalut (' Allocate mini maps/masks of size :',1,minmsk)
      idptr = fmalloc (nb)
      ieptr = fmalloc (nb)
      igptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0 .or.
     +    idptr .eq. 0 .or. ieptr .eq. 0 .or. igptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0 .or.
     +    idptr .le. 0 .or. ieptr .le. 0 .or. igptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call docoma (%val(iaptr), %val(ibptr), %val(icptr),
     +             %val(idptr), %val(igptr), %val(ieptr),
     +             mapsize, masksize, minmsk)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
      call ffree (idptr)
      call ffree (ieptr)
      call ffree (igptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine docoma (mapin, mapout, mapdum,
     +                   mapi, mapj, maski,
     +                   maxsiz, maxmsk, minmap)
c
c ... COrrelation MAsk
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz,maxmsk,minmap,maxncs,maxsym
      parameter (maxncs = maxgk2)
      parameter (maxsym = maxgk2)
c
      real mapin(maxsiz),mapout(maxsiz),mapj(minmap),mapi(minmap)
      real mapdum(maxsiz)
      real grid(3),cella(6),cellb(6),spacb(3),testar(2,3)
      real rtatob(12,maxncs),rtsym(12,maxsym)
      real radius
c
      integer maski(minmap)
      integer orgna(3),exta(3),grida(3),gridb(3),uvwa(3)
      integer orgnx(3),extx(3)
      integer spgrp,ierr,iunit,length,ndum
      integer ctrt,ctsym,errcod,idum,i,j,itype
c
      character file*256
      character fmt*80,par*25,partyp*1
c
      logical xinter,linter
c
      data testar / 0.0,0.5, 0.0,0.5, 0.0,0.5 /
      data radius /4.5/, itype/1/
c
code ...
c
      iunit = 11
      linter = xinter()
c
c ... check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
c ... input CCP4 map
c
      write (*,*)
      file = ' '
      call textin (' Input CCP4 map file ?',file)
      call textut (' Input CCP4 map file :',file)
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
c ... read in map
c
      call edin (file, iunit, mapin,  orgna, exta, grida, uvwa,
     +           cella, spgrp, maxsiz)
c
      close (iunit)
      call telmap (grida,orgna,exta,cella)
c
      do i=1,3
        grid (i) = cella (i) / float (grida(i))
      end do
c
c ... assume grid = 1/3 * resolution; sphere radius = 1.5 * resolution
c     ==> guess radius as 1.5 * 3.0 * smallest spacing
c         and round to 0.1
c
      call fvalut (' Grid spacing (A):',3,grid)
      radius = 1.5 * 3.0 * min(grid(1),grid(2),grid(3))
      radius = 0.1 * float(nint(10.0*radius))
c
c ... get spacegroup symmetry operators
c
      write (*,*)
      file = 'symop.o'
      call textin (' File with spacegroup symmetry operators ?',file)
      call textut (' File with spacegroup symmetry operators :',file)
      if (file .eq. ' ') then
        call errcon ('No filename provided')
        return
      end if
c
      call osymop (1,file,errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening O datablock file')
        return
      end if
      close (1)
c
      call opoodb (1,file,par,partyp,j,fmt,errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening O datablock file')
        return
      end if
c
      ctsym = j/12
      if (ctsym .gt. maxsym) then
        call errcon ('Too many spacegroup operators')
        return
      end if
      if ((12*ctsym) .ne. j) then
        call errcon ('Invalid nr of elements; must be N*12')
        return
      end if
      do 220 j=1,ctsym
220     read (1, fmt,err=999,end=999) (rtsym(i, j),i=1,12)
c
      do j=1,ctsym
        call fratra (rtsym(10,j))
      end do
      close (1)
c
      call anasgs (ctsym,rtsym,.true.,ierr)
      if (ierr .ne. 0) then
        call errcon ('In spacegroup symmetry operators')
        return
      end if
c
c ... get NCS operators
c
      ctrt = 0
c
210   continue
      file = ' '
      write (*,*)
      call textin (' File with NCS RT-operator(s) (RETURN to quit) ?',
     +  file)
      call textut (' File with NCS RT-operator(s) :',file)
c
      if (file .eq. ' ') goto 200
c
c ... read one or many NCS operators
c
      call rdoncs (1,file,ctrt,maxncs,rtatob,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading NCS operators')
        return
      end if
c
      goto 210
c
200   continue
      write (*,*)
      if (ctrt .lt. 2) then
        call errcon ('Fewer than 2 NCS operators')
        return
      end if
c
c      call anancs (ctrt,rtatob,.true.,ierr)
c      if (ierr .ne. 0)
c     +  call errstp ('In NCS operators')
c
c ... get limits of test area (fractional coords)
c
      write (*,*)
      call prompt (
     +  ' Define the area where you expect your "reference"')
      call prompt (
     +  ' molecule to be in terms of fractional coordinates.')
      call fvalin (' Fractional search range A ?',2,testar(1,1))
      call fvalin (' Fractional search range B ?',2,testar(1,2))
      call fvalin (' Fractional search range C ?',2,testar(1,3))
c
      call rlohi (testar(1,1),testar(2,1))
      call rlohi (testar(1,2),testar(2,2))
      call rlohi (testar(1,3),testar(2,3))
c
      call fvalut (' Fractional search range A :',2,testar(1,1))
      call fvalut (' Fractional search range B :',2,testar(1,2))
      call fvalut (' Fractional search range C :',2,testar(1,3))
c
c ... get sphere radius
c
      write (*,*)
      call prompt (
     +  ' Sphere radius ~ 1.5 * resolution of trusted phases')
      call fvalin (' Sphere radius (A) ?',1,radius)
      if (radius .le. 1.0) then
        radius = 4.0
        call errcon ('Radius too small; reset to 4.0 A !')
      end if
      call fvalut (' Sphere radius (A) :',1,radius)
c
c ... type of correlation map
c
      write (*,*)
      call prompt (
     +  ' Which type of correlation map do you want ?')
      call prompt (
     +  ' (Note: R1 = R2 for two-fold NCS)')
      call prompt (' 1 = R1 map (<CC(1,J); J=2,N>)')
      call prompt (' 2 = R2 map (<CC(I,J); I=1,N-1; J=I+1,N>)')
      call ivalin (' Type of map (1,2) ?',1,itype)
      if (itype .ne. 2) itype = 1
      call ivalut (' Type of map :',1,itype)
c
c ... output map name
c
      write (*,*)
      file = 'coma.E'
      call textin (' Output CCP4 correlation map file ?',file)
      call textut (' Output CCP4 correlation map file :',file)
c
c ... same cell and grid for sphere maps and masks and output map
c
      do i=1,3
        cellb (i) = cella (i)
        cellb (i+3) = cella (i+3)
        gridb (i) = nint (cellb(i)/(radius/3.0))
        spacb (i) = cellb(i) / float(gridb(i))
      end do
c
      write (*,*)
      call ivalut (' Output map grid :',3,gridb)
      call fvalut (' Spacing (A)     :',3,spacb)
c
c ... calc origin and extent of output map (test area)
c
      do i=1,3
        orgnx(i) = -1 + int (testar(1,i)*gridb(i))
        j = 1 + int(testar(2,i)*gridb(i))
        extx(i) = j - orgnx(i) + 1
      end do
c
      write (*,*)
      call jvalut (' Origin :',3,orgnx)
      call jvalut (' Extent :',3,extx)
c
      idum = extx(1)*extx(2)*extx(3)
      call jvalut (' Calculated output map size :',1,idum)
      if (idum .gt. maxsiz) then
        call errcon ('Too big for map')
        return
      else if (idum .gt. maxmsk) then
        call errcon ('Too big for mask')
        return
      end if
c
c ... initialise output map
c
      do i=1,idum
        mapout (i) = 0.0
      end do
c
c      call ivalut (' Nr NCS  :',1,ctrt)
c      call ivalut (' SymOps  :',1,ctsym)
c
c ... now go do everything
c
      if (itype .eq. 1) then
        call calcr1 (
     +    mapin,orgna,exta(1),exta(2),exta(3),grida,cella,grid,
     +    mapout,orgnx,extx(1),extx(2),extx(3),gridb,cellb,spacb,
     +    mapi,mapj,maski,minmap,
     +    ctsym,rtsym,
     +    ctrt,rtatob,
     +    radius)
c
      else
c
c ... calculate how many small maps will fit into the big one
c
        j = 1 + (radius/min(spacb(1),spacb(2),spacb(3)))
        idum = (2*j+7)**3 + 10
        ndum = maxsiz / ( idum + 10 )
        call jvalut (' Size of sphere maps :',1,idum)
        call jvalut (' Number storeable    :',1,ndum)
c
        call calcr2 (
     +    mapin,orgna,exta(1),exta(2),exta(3),grida,cella,grid,
     +    mapout,orgnx,extx(1),extx(2),extx(3),gridb,cellb,spacb,
     +    mapi,mapj,maski,minmap,mapdum,idum,ndum,
     +    ctsym,rtsym,
     +    ctrt,rtatob,
     +    radius)
c
      end if
c
c ... write output map
c
      write (*,*)
      call edout (file,iunit,mapout,orgnx,extx,gridb,uvwa,cellb,spgrp)
c
c ... all done
c
      return
c
  999 continue
      call errcon ('While reading operators')
      return
c
      end
c
c
c
      subroutine calcr1 (
     +  mapin,orgna,exta1,exta2,exta3,grida,cella,grid,
     +  mapout,orgnx,extx1,extx2,extx3,gridb,cellb,spacb,
     +  mapi,mapj,maski,minmap,
     +  ctsym,rtsym,
     +  ctrt,rtatob,
     +  radius)
c
      implicit none
c
      integer maxhis
      parameter (maxhis=41)
c
      integer exta1,exta2,exta3,orgna(3),grida(3)
      real mapin(exta1,exta2,exta3),cella(6),grid(3)
c
      integer extx1,extx2,extx3,orgnx(3),gridb(3)
      real mapout(extx1,extx2,extx3),cellb(6),spacb(3)
c
      integer minmap
      real mapi(minmap),mapj(minmap)
      integer maski(minmap)
c
      integer ctsym
      real rtsym (12,ctsym)
c
      integer ctrt
      real rtatob (12,ctrt)
c
      real radius
c
c ... local variables
c
      real x(3),x2(3),a(3,3),b(3,3),rtunit(12),g(3)
      real avy(2),avysq(2),avxy(2),avxmy(2)
      real sumcc,xdum,cc,r1,r2,rmsd,gmin,value,volume,voxva
c
      integer orgnz(3),extz(3),ihisto(maxhis)
      integer i,i1,i2,i3,j,ct,ierr,ncc,off,bobo1,bobo2,bobo3,ngk
c
      logical lmask
c
      data rtunit / 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0,
     +  0.0,0.0,0.0 /
c
code ...
c
      do i=1,maxhis
        ihisto (i) = 0
      end do
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cella, a, 0)
      call orthog (cella, b, 1)
c
c      call ivalut (' Nr NCS  :',1,ctrt)
c      call ivalut (' SymOps  :',1,ctsym)
c
      do i=1,3
        g(i) = cellb(i)/float(gridb(i))
      end do
      gmin = min (g(1),g(2),g(3))
      off = 1 + radius/gmin
c
      lmask = .false.
      ncc = ctrt - 1
      bobo1 = extx1*extx2*extx3
      bobo2 = bobo1/20
      bobo3 = bobo2
      ngk = 0
c
c      call jvalut (' Nr of points to test :',1,bobo1)
c
c ... loop over test area
c
      do i1 = orgnx(1),orgnx(1)+extx1-1
        do i2 = orgnx(2),orgnx(2)+extx2-1
          do i3 = orgnx(3),orgnx(3)+extx3-1
c
c           print *,'I1,2,3 = ',i1,i2,i3
c
c ... calculate coordinates for this point
c
c            x(1) = (i1-1+ orgnx(1))*spacb(1)/cellb(1)
c            x(2) = (i2-1+ orgnx(2))*spacb(2)/cellb(2)
c            x(3) = (i3-1+ orgnx(3))*spacb(3)/cellb(3)
c
            x(1) = (i1)*spacb(1)/cellb(1)
            x(2) = (i2)*spacb(2)/cellb(2)
            x(3) = (i3)*spacb(3)/cellb(3)
            call mulmtx (a, x, x2, 3, 3, 1)
c
c            print *,' FRACT ',x
c            print *,' CART  ',x2
c
c ... generate spherical mask around this point
c
c            do i=1,3
c              g(i) = cellb(i)/float(spacb(i))
c              orgnz(i) = nint (x(i)- radius/g(i))-1
c              j = nint (x(i)+ radius/g(i))+1
c              extz(i) = j - orgnz(i)+1
c            end do
c
            orgnz (1) = i1 - off - 3
            orgnz (2) = i2 - off - 3
            orgnz (3) = i3 - off - 3
            do i=1,3
              extz (i) = 2*(off+3) + 1
            end do
c
c ... only need to generate the actual spherical mask once !!!
c
            if (.not. lmask) then
c
              ct = 1
              do i=1,3
                ct = ct * extz(i)
              end do
c
              if (ct .gt. minmap) then
                call errstp ('Sphere map/mask too big')
              end if
c
              call grpmsk (cellb,gridb,radius,x2,orgnz,
     +                     extz(1),extz(2),extz(3),maski,ierr)
c
              j = 0
              do i=1,extz(1)*extz(2)*extz(3)
                if (maski(i).eq.1) j=j+1
              end do
c
c              call ivalut (' Nr mask points :',1,j)
c              call ivalut (' Zorigin :',3,orgnz)
c              call ivalut (' Zextent :',3,extz)
c              call ivalut (' Nr NCS  :',1,ctrt)
c              call ivalut (' SymOps  :',1,ctsym)
c
              lmask = .true.
c
            end if
c
c ... calculate R1 coefficient
c
            sumcc = 0.0
c
c ... get density inside sphere in mol 1
c
            call sub02x (
     +        mapin,exta1,exta2,exta3,orgna,grid,cella,
     +        mapi,extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +        maski, rtunit, 1, rtsym, ctsym,
     +        rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
c
c ... loop over the other molecules
c
            do i=2,ctrt
c
c ... get density for sphere in mol "i"
c
              call sub02x (
     +          mapin,exta1,exta2,exta3,orgna,grid,cella,
     +          mapj,extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +          maski, rtatob(1,i), 1, rtsym, ctsym,
     +          rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
c
c ... calc corr coeff with profile map
c
              call simsim (mapi, mapj, maski,
     +          extz(1), extz(2), extz(3),
     +          cc,r1,r2,rmsd)
c
c              print *,i1,i2,i3,cc,r1
c
              sumcc = sumcc + cc
c
            end do
c
c            print *,'Indices = ',i1-orgnx(1)+1,
c     +              i2-orgnx(2)+1,
c     +              i3-orgnx(3)+1
c
            value = sumcc / float(ncc)
            mapout (i1-orgnx(1)+1,
     +              i2-orgnx(2)+1,
     +              i3-orgnx(3)+1) = value
c
            i = 1 + int ((value+1.0)/0.05)
            ihisto (i) = ihisto (i) + 1
c
            ngk = ngk + 1
            if (ngk .eq. bobo3) then
              bobo3 = bobo3 + bobo2
              xdum = 100.0 * float(ngk)/float(bobo1)
              call fvalut (' Progress (%) :',1,xdum)
            end if
c
          end do
        end do
      end do
c
c ... histogram etc.
c
      write (*,*)
      call jvalut (' Nr of points :',1,ngk)
      call voxvol (cellb,gridb,volume,voxva)
      call rvalut (' Cell volume (A**3)  :',1,volume)
      call rvalut (' Voxel volume (A**3) :',1,voxva)
      j = 0
      do i=maxhis,1,-1
        j = j + ihisto(i)
        if (ihisto(i) .gt. 0) then
          write (*,6000) i,0.05*float(i)-1.0,ihisto(i),
     +      j,float(j)*voxva,
     +      nint(100.0*( volume - ctsym*ctrt*j*voxva ) / volume)
        end if
      end do
      write (*,*)
      call prompt (
     +  ' "Cumul Vol" ignores any mask overlap !')
      call prompt (
     +  ' "Solv Cont" ignores overlap, assumes one domain,')
      call prompt (
     +  '   improper NCS, and all NCS-RT and SymmOps used here !')
      write (*,*)
c
 6000 format (' Bin ',i2,' ... R1 >= ',f5.2,' Nr ',i8,
     +  ' Cum ',i8,' Cum Vol ',1pe12.4,' Solv Cont (%) ',i5)
c
      return
      end
c
c
c
      subroutine calcr2 (
     +  mapin,orgna,exta1,exta2,exta3,grida,cella,grid,
     +  mapout,orgnx,extx1,extx2,extx3,gridb,cellb,spacb,
     +  mapi,mapj,maski,minmap,mapdum,istore,nstore,
     +  ctsym,rtsym,
     +  ctrt,rtatob,
     +  radius)
c
      implicit none
c
      integer maxhis
      parameter (maxhis=41)
c
      integer exta1,exta2,exta3,orgna(3),grida(3)
      real mapin(exta1,exta2,exta3),cella(6),grid(3)
c
      integer extx1,extx2,extx3,orgnx(3),gridb(3)
      real mapout(extx1,extx2,extx3),cellb(6),spacb(3)
c
      integer minmap
      real mapi(minmap),mapj(minmap)
      integer maski(minmap)
c
      integer istore,nstore
      real mapdum(istore,nstore)
c
      integer ctsym
      real rtsym (12,ctsym)
c
      integer ctrt
      real rtatob (12,ctrt)
c
      real radius
c
c ... local variables
c
      real x(3),x2(3),a(3,3),b(3,3),rtunit(12),g(3)
      real avy(2),avysq(2),avxy(2),avxmy(2)
      real sumcc,xdum,cc,r1,r2,rmsd,gmin,value,volume,voxva
c
      integer orgnz(3),extz(3),ihisto(maxhis)
      integer i,i1,i2,i3,j,k,ct,ierr,ncc,off,bobo1,bobo2,bobo3,ngk
c
      logical lmask
c
      data rtunit / 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0,
     +  0.0,0.0,0.0 /
c
code ...
c
      do i=1,maxhis
        ihisto (i) = 0
      end do
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cella, a, 0)
      call orthog (cella, b, 1)
c
c      call ivalut (' Nr NCS  :',1,ctrt)
c      call ivalut (' SymOps  :',1,ctsym)
c
      do i=1,3
        g(i) = cellb(i)/float(gridb(i))
      end do
      gmin = min (g(1),g(2),g(3))
      off = 1 + radius/gmin
c
      lmask = .false.
      ncc = (ctrt * (ctrt - 1)) / 2
      bobo1 = extx1*extx2*extx3
      bobo2 = bobo1/20
      bobo3 = bobo2
      ngk = 0
c
c      call jvalut (' Nr of points to test :',1,bobo1)
c
c ... loop over test area
c
      do i1 = orgnx(1),orgnx(1)+extx1-1
        do i2 = orgnx(2),orgnx(2)+extx2-1
          do i3 = orgnx(3),orgnx(3)+extx3-1
c
c           print *,'I1,2,3 = ',i1,i2,i3
c
c ... calculate coordinates for this point
c
c            x(1) = (i1-1+ orgnx(1))*spacb(1)/cellb(1)
c            x(2) = (i2-1+ orgnx(2))*spacb(2)/cellb(2)
c            x(3) = (i3-1+ orgnx(3))*spacb(3)/cellb(3)
c
            x(1) = (i1)*spacb(1)/cellb(1)
            x(2) = (i2)*spacb(2)/cellb(2)
            x(3) = (i3)*spacb(3)/cellb(3)
            call mulmtx (a, x, x2, 3, 3, 1)
c
c            print *,' FRACT ',x
c            print *,' CART  ',x2
c
c ... generate spherical mask around this point
c
c            do i=1,3
c              g(i) = cellb(i)/float(spacb(i))
c              orgnz(i) = nint (x(i)- radius/g(i))-1
c              j = nint (x(i)+ radius/g(i))+1
c              extz(i) = j - orgnz(i)+1
c            end do
c
            orgnz (1) = i1 - off - 3
            orgnz (2) = i2 - off - 3
            orgnz (3) = i3 - off - 3
            do i=1,3
              extz (i) = 2*(off+3) + 1
            end do
c
c ... only need to generate the actual spherical mask once !!!
c
            if (.not. lmask) then
c
              ct = 1
              do i=1,3
                ct = ct * extz(i)
              end do
c
              if (ct .gt. minmap) then
                call errstp ('Sphere map/mask too big')
              end if
c
              call grpmsk (cellb,gridb,radius,x2,orgnz,
     +                     extz(1),extz(2),extz(3),maski,ierr)
c
              j = 0
              do i=1,extz(1)*extz(2)*extz(3)
                if (maski(i).eq.1) j=j+1
              end do
c
c              call ivalut (' Nr mask points :',1,j)
c              call ivalut (' Zorigin :',3,orgnz)
c              call ivalut (' Zextent :',3,extz)
c              call ivalut (' Nr NCS  :',1,ctrt)
c              call ivalut (' SymOps  :',1,ctsym)
c
              lmask = .true.
c
            end if
c
c ... calculate R2 coefficient
c
            sumcc = 0.0
c
            do i=1,min(ctrt,nstore)
              call sub02x (
     +          mapin,exta1,exta2,exta3,orgna,grid,cella,
     +          mapdum(1,i),extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +          maski, rtatob(1,i), 1, rtsym, ctsym,
     +          rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
            end do
c
            do k=1,ctrt-1
c
c ... get density for sphere in mol K
c
              if (k .gt. nstore) then
              call sub02x (
     +          mapin,exta1,exta2,exta3,orgna,grid,cella,
     +          mapi,extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +          maski, rtatob(1,k), 1, rtsym, ctsym,
     +          rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
              end if
c
c ... loop over the other molecules
c
              do i=k+1,ctrt
c
c ... get density for sphere in mol "i"
c
                if (i .gt. nstore) then
                call sub02x (
     +            mapin,exta1,exta2,exta3,orgna,grid,cella,
     +            mapj,extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +            maski, rtatob(1,i), 1, rtsym, ctsym,
     +            rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
                end if
c
c ... calc corr coeff with profile map
c
                if (i .le. nstore .and. k .le. nstore) then
                  call simsim (mapdum(1,i), mapdum(1,k), maski,
     +              extz(1), extz(2), extz(3),
     +              cc,r1,r2,rmsd)
                else if (i .le. nstore) then
                  call simsim (mapdum(1,i), mapj, maski,
     +              extz(1), extz(2), extz(3),
     +              cc,r1,r2,rmsd)
                else if (k .le. nstore) then
                  call simsim (mapi, mapdum(1,k), maski,
     +              extz(1), extz(2), extz(3),
     +              cc,r1,r2,rmsd)
                else
                  call simsim (mapi, mapj, maski,
     +              extz(1), extz(2), extz(3),
     +              cc,r1,r2,rmsd)
                end if
c
c                print *,i1,i2,i3,cc,r1
c
                sumcc = sumcc + cc
c
              end do
c
            end do
c
            value = sumcc / float(ncc)
            mapout (i1-orgnx(1)+1,
     +              i2-orgnx(2)+1,
     +              i3-orgnx(3)+1) = value
c
            i = 1 + int ((value+1.0)/0.05)
            ihisto (i) = ihisto (i) + 1
c
            ngk = ngk + 1
            if (ngk .eq. bobo3) then
              bobo3 = bobo3 + bobo2
              xdum = 100.0 * float(ngk)/float(bobo1)
              call fvalut (' Progress (%) :',1,xdum)
            end if
c
          end do
        end do
      end do
c
c ... histogram etc.
c
      write (*,*)
      call jvalut (' Nr of points :',1,ngk)
      call voxvol (cellb,gridb,volume,voxva)
      call rvalut (' Cell volume (A**3)  :',1,volume)
      call rvalut (' Voxel volume (A**3) :',1,voxva)
      j = 0
      do i=maxhis,1,-1
        j = j + ihisto(i)
        if (ihisto(i) .gt. 0) then
          write (*,6000) i,0.05*float(i)-1.0,ihisto(i),
     +      j,float(j)*voxva,
     +      nint(100.0*( volume - ctsym*ctrt*j*voxva ) / volume)
        end if
      end do
      write (*,*)
      call prompt (
     +  ' "Cumul Vol" ignores any mask overlap !')
      call prompt (
     +  ' "Solv Cont" ignores overlap, assumes one domain,')
      call prompt (
     +  '   improper NCS, and all NCS-RT and SymmOps used here !')
      write (*,*)
c
 6000 format (' Bin ',i2,' ... R2 >= ',f5.2,' Nr ',i8,
     +  ' Cum ',i8,' Cum Vol ',1pe12.4,' Solv Cont (%) ',i5)
c
      return
      end
c
c
c
      subroutine grpmsk (cell,grid,radius,xyzr,origin,
     +                   exta1,exta2,exta3,mask,ierr)
c
      implicit none
c
      integer i,j,k,l,m,extent(3),grid(3),origin(3),ix2(3),ip(3)
      integer ct,i1,i2,i3,off,ierr,size,inx,nxy
      integer exta1,exta2,exta3
      integer mask(exta1,exta2,exta3)
c
      real cell(6),spacing(3),distce,x(3),xp(3),xyzr(3)
      real a(3,3),b(3,3),c(3,3)
      real x1(3),x2(3),gmin
      real radius
c
code ...
c
      ierr = -1
c
      extent(1)=exta1
      extent(2)=exta2
      extent(3)=exta3
c
      do i=1,3
        spacing(i) = cell(i)/float(grid(i))
      end do
c
c --- Blank mask
c
      do i1=1,exta1
        do i2=1,exta2
          do i3=1,exta3
            mask(i1,i2,i3) = 0
          end do
        end do
      end do
c
c ... B = cart to fract
c     C = fract to cart
c
      call orthog (cell, c, 0)
      call orthog (cell, b, 1)
c
c --- Do the actual mask job
c
      gmin = min (spacing(1),spacing(2),spacing(3))
      off = 1 + radius/gmin
c
c ... X1 = fractional coordinates of the point
c
      call mulmtx (b, xyzr, x1, 3, 3, 1)
c
c ... IX2 = coordinates in grid points
c
      do i=1,3
        ix2 (i) = nint(x1(i)*cell(i)/spacing(i))
      end do
c
c      call fvalut (' Cart :',3,xyzr)
c      call fvalut (' Frac :',3,x1)
c      call ivalut (' Grid :',3,ix2)
c      print *,off,radius,gmin
c
      m = 0
      do j=-off,off	
        do k=-off,off	
          do l=-off,off
c
c ... IP = grid point to be tested
c
            ip (1) = ix2 (1) + l
            ip (2) = ix2 (2) + k
            ip (3) = ix2 (3) + j
c
c ... XP = fractional coords
c
            xp(1) = float(ip(1))*spacing(1)/cell(1)
            xp(2) = float(ip(2))*spacing(2)/cell(2)
            xp(3) = float(ip(3))*spacing(3)/cell(3)
c
c ... X2 = cartesian coords
c
            call mulmtx (c, xp, x2, 3, 3, 1)
c      call fvalut (' XP :',3,xp)
c      call fvalut (' X2 :',3,x2)
            if (distce(xyzr,x2) .gt. radius) goto 130
c            i1 = nint(xp(1)/g(1))- origin(1) +1
c            i2 = nint(xp(2)/g(2))- origin(2) +1
c            i3 = nint(xp(3)/g(3))- origin(3) +1
            i1 = ip(1) - origin(1) + 1
            i2 = ip(2) - origin(2) + 1
            i3 = ip(3) - origin(3) + 1
c            print *,' ..... I1,2,3 ',i1,i2,i3
            if (i1 .le. 0) goto 130
            if (i2 .le. 0) goto 130
            if (i3 .le. 0) goto 130
            if (i1 .gt. extent(1)) goto 130
            if (i2 .gt. extent(2)) goto 130
            if (i3 .gt. extent(3)) goto 130
ccc            inx = (i3-1)*nxy + (i2-1)*extent(1) + i1
            mask(i1,i2,i3) = 1
            m = m + 1
 130        continue
          end do
        end do
      end do
c
      call ivalut (' Nr of mask points in sphere :',1,m)
c
c --- DONE
c
      ierr = 0
c
      return
c
      end
c
c
c
      subroutine simsim (mapi, mapj, maski, ext1, ext2, ext3,
     +      ccoef,sumfo,sumfc,rmsd)
c
      implicit none
c
      integer ext1, ext2, ext3, ntot, i1,i2,i3
      integer maski(ext1, ext2, ext3)
c
      real mapi(ext1, ext2, ext3),mapj(ext1, ext2, ext3)
      real avx,avy,avxy,avxsq,avysq,sumfo,sumfc,fofc,rmsd
      real avalue,bvalue,aabs,babs,f,ccoef,shape
c
code ...
c
      ntot = 0
      avx = 0.
      avy = 0.
      avxy = 0.
      avxsq = 0.
      avysq = 0.
      sumfo = 0.
      sumfc = 0.
      fofc = 0.
      rmsd = 0.
c
      do i1=1,ext1
        do i2=1,ext2
          do i3=1,ext3
            if (maski(i1,i2,i3) .eq. 1) then
              ntot = ntot + 1
              avalue = mapi(i1,i2,i3)
              bvalue = mapj(i1,i2,i3)
              aabs = abs(avalue)
              babs = abs(bvalue)
              avx = avx+ avalue
              avy = avy+ bvalue
              avxsq = avxsq+ avalue**2
              avysq = avysq+ bvalue**2
              avxy = avxy+ avalue*bvalue
              fofc = fofc + abs(aabs-babs)
              sumfo = sumfo + aabs
              sumfc = sumfc + babs
              rmsd = rmsd + (avalue-bvalue)**2
            end if
          end do
        end do
      end do
c
c      call jvalut (' Nr of points in common grid :',1,ntot)
      if (ntot .le. 0) return
c
      f = float(ntot)
      ccoef = (avxy/f- avx*avy/(f*f))/ 
     +    (sqrt(avxsq/f- (avx/f)**2)* sqrt(avysq/f- (avy/f)**2))
c      call fvalut (' Correlation coefficient :',1,ccoef)
c
      sumfo = fofc/sumfo
c      call fvalut (' R-factor w.r.t. map 1   :',1,sumfo)
c
      sumfc = fofc/sumfc
c      call fvalut (' R-factor w.r.t. map 2   :',1,sumfc)
c
      rmsd = sqrt (rmsd/f)
c      call rvalut (' RMS difference          :',1,rmsd)
c
      shape = avxy / (sqrt(avxsq)*sqrt(avysq))
c      call fvalut (' Shape similarity index  :',1,shape)
c
c      call prompt (' R-factors based on UNSCALED data !')
c
      return
      end
c
c
c
      subroutine sub02x (
     +    mapin,exta1,exta2,exta3,orgna,spaca,cella,
     +    mapout,extb1,extb2,extb3,orgnb,spacb,cellb,
     +    maski, rtatob, ctrt, rtsym,ctsym, or2or, ierr,
     +    avy,avysq,avxy,avxmy,lposit)
c
c ... average into a (possibly) different spacegroup
c
c ... Gerard Kleywegt, 930317,18
c     based on code by Alwyn Jones, 1991
c
      implicit none
c
      real small
      parameter (small=1.0e-9)
c
      integer orgna(3), exta1, exta2, exta3
      real mapin(exta1, exta2, exta3),spaca(3),cella(6)
c
      integer orgnb(3), extb1, extb2, extb3
      real mapout(extb1, extb2, extb3),spacb(3),cellb(6)
      integer maski(extb1, extb2, extb3)
c
      real rtatob(12,*),rtsym(12,*),or2or(12),val1,f,cc
      real avy(*),avysq(*),avxy(*),avxmy(*),avx,avxsq,avax
      integer ctrt,ctsym,ierr,nbit,ngjk
c
      integer errcod, i, j, k, l, loop, m, ext(3), ijk(3)
c      integer bobo1,bobo2,bobo3
      integer nerr1,nerr2
c
      real forgna(3), fexta(3), mapoutit(4,4,4)
      real value, x(3), gexta(3)
c
      real f2ca(3,3),c2fa(3,3),f2cb(3,3),c2fb(3,3)
      real fmsk(3),cmsk(3),cmska(3),ncmska(3),fcmska(3),xdum
c
      logical lposit
c
      equivalence (ijk(1),i)
      equivalence (ijk(2),j)
      equivalence (ijk(3),k)
c
code ...
c
      ierr = -1
      nbit = 0
c
      call orthog (cella, f2ca, 0)
      call orthog (cella, c2fa, 1)
      call orthog (cellb, f2cb, 0)
      call orthog (cellb, c2fb, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
c
c ... get edges of map A in fractional coordinates
c
      do i=1,3
        forgna(i) = float(orgna(i))         *spaca(i)/cella(i)
        fexta(i)  = float(ext(i)+orgna(i)-1)*spaca(i)/cella(i)
        gexta(i)  = float(ext(i)+orgna(i)-2)*spaca(i)/cella(i)
      end do
c
c      call fvalut (' FORGNA :',3,forgna)
c      call fvalut (' FEXTA  :',3,fexta)
c      call fvalut (' GEXTA  :',3,gexta)
c
      avx   = 0.0
      avax  = 0.0
      avxsq = 0.0
      do loop=1,ctrt
        avy(loop)    = 0.0
        avysq(loop)  = 0.0
        avxy (loop)  = 0.0
        avxmy (loop) = 0.0
      end do
      ngjk = 0
      nerr1 = 0
      nerr2 = 0
c
c      call cntmsk (maski,extb1,extb2,extb3,bobo1)
c      call jvalut (' Points in mask :',1,bobo1)
c      bobo2 = bobo1/10
c      bobo3 = bobo2
c
c ... loop over the mask points
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        mapout (i,j,k) = 0.0
        if (maski(i,j,k) .ne. 1) goto 100
c
        ngjk = ngjk + 1
c
c        if (ngjk .eq. bobo3) then
c          xdum = 100.0 * float(ngjk-1)/float(bobo1)
c          call fvalut (' Progress (% mask) :',1,xdum)
c          bobo3 = bobo3 + bobo2
c        end if
c
c ... FMSK = fractional coords of mask point in B
c
        do l=1,3
          fmsk (l) = float(ijk(l)-1+orgnb(l))*spacb(l)/cellb(l)
        end do
c
c ... CMSK = cartesian coords of mask point in B
c     CMSKA = ditto, but now in A
c
        call mulmtx (f2cb,fmsk,cmsk,3,3,1)
        call vecrtv (cmsk,cmska,1,or2or(1),or2or(10))
c
c ... loop over the NCS operators in B
c
        do 120 loop=1,ctrt
c
c ... NCMSKA = ditto, in A, but in NCS mate
c     FCMSKA = ditto, but in fractional coords
c
          call vecrtv (cmska,ncmska,1,rtatob(1,loop),rtatob(10,loop))
          call mulmtx (c2fa,ncmska,fcmska,3,3,1)
c
c ... force this point into the asymmetric unit of map A
c
          call frcsym (fcmska,forgna,gexta,rtsym,ctsym,errcod)
c
c ... if not possible, it is close to the extent; build a 4x4x4 map
c
          if (errcod .ne. 0) then
c
            nbit = nbit + 1
            call bldbit (
     +        mapin,exta1,exta2,exta3,orgna,spaca,cella,
     +        forgna,fexta,rtsym,ctsym,mapoutit,fcmska,errcod)
c
            if (errcod .ne. 0) then
              if (nerr1 .lt. 10) then
                call errcon ('Serious FRCSYM error')
                call ivalut (' Mask point in reference  :',3,ijk)
                call fvalut (' Fractional               :',3,fmsk)
                call fvalut (' Cartesian                :',3,cmsk)
                call fvalut (' Cartesian in second xtal :',3,cmska)
                call fvalut (' Cartesian of NCS mate    :',3,ncmska)
                call fvalut (' Fractional               :',3,fcmska)
              else if (nerr1 .eq. 10) then
                call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
              end if
              nerr1 = nerr1 + 1
              goto 120
            end if
c
c ... interpolate in 4x4x4 map
c
            do l=1,3
              x(l) = fcmska(l)*cella(l)/spaca(l)- float(orgna(l))+ 1
              m = x(l)
              x(l) = x(l)- float(m-1)
            end do
            call intrpl (mapoutit, 4, 4, 4, x, value, errcod)
c
c ... if first FRCSYM worked okay, do straightforward interpolation
c
          else
c
            do l=1,3
              x(l) = fcmska(l)*cella(l)/spaca(l)- float(orgna(l)) + 1
            end do
c
            call intrpl (mapin, exta1, exta2, exta3, x, value, errcod)
c
          end if
c
          if (errcod .eq. 0) then
            mapout (i,j,k) = mapout (i,j,k) + value
c
              if (loop .eq. 1) then
                val1 = value
                avx = avx + value
                avxsq = avxsq + value*value
                avax  = avax + abs(value)
              end if
c
              avy (loop)   = avy(loop)   + value
              avysq (loop) = avysq(loop) + value*value
              avxy (loop)  = avxy(loop)  + value*val1
              avxmy (loop) = avxmy(loop) + abs(value-val1)
c
          else
            if (nerr2 .lt. 10) then
              call errcon ('Interpolation error')
              call ivalut (' Mask point in reference  :',3,ijk)
              call fvalut (' Fractional               :',3,fmsk)
              call fvalut (' Cartesian                :',3,cmsk)
              call fvalut (' Cartesian in second xtal :',3,cmska)
              call fvalut (' Cartesian of NCS mate    :',3,ncmska)
              call fvalut (' Fractional               :',3,fcmska)
            else if (nerr2 .eq. 10) then
              call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
            end if
            nerr2 = nerr2 + 1
          end if
  120   continue
  100 continue
c
      ierr = 0
c
c      call jvalut (' Calls to BLDBIT      :',1,nbit)
      if (nerr1.gt.0) call jvalut (' Severe FRCSYM errors :',1,nerr1)
      if (nerr1.gt.0) call jvalut (' Interpolation errors :',1,nerr2)
c
c ... print correlation coefficients for the various operators
c
c      call jvalut (' Nr of mask points :',1,ngjk)
      f = float (ngjk)
      do loop=1,ctrt
        cc = (avxy(loop)/f- avx*avy(loop)/(f*f))/ 
     +    (sqrt(avxsq/f- (avx/f)**2) *
     +     sqrt(avysq(loop)/f- (avy(loop)/f)**2))
c        write (*,'(a,i3,a,f8.5)') ' Corr. coeff. for operator ',
c     +    loop,' = ',cc
        cc = avxmy(loop) / avax
c        write (*,'(a,i3,a,f8.5)') ' R-factor for operator ',
c     +    loop,' w.r.t. operator 1 = ',cc
      end do
c
c ... average (if necessary)
c
      if (lposit) then
        call rvalut (' Positivity; set <= 0 to:',1,small)
        xdum = 1.0 / float(ctrt)
        do 300 k=1,extb3
        do 300 j=1,extb2
        do 300 i=1,extb1
c ... changed next line for HP
          if (maski(i,j,k) .eq. 1) then
            if (mapout(i,j,k) .gt. 0.0) then
              mapout(i,j,k) = mapout(i,j,k) * xdum
            else
              mapout(i,j,k) = small
            end if
          end if
  300   continue
      else
c        call prompt (' Averaging without positivity constraint')
        if (ctrt .gt. 1) then
          xdum = 1.0 / float(ctrt)
          do 200 k=1,extb3
          do 200 j=1,extb2
          do 200 i=1,extb1
200         mapout(i,j,k) = mapout(i,j,k) * xdum
        end if
      end if
c
      return
      end
