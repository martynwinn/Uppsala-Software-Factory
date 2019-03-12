      program prof
c
c ... Electron-density PROFiles
c
c ... Gerard Kleywegt @ 961127
c
c ... Version 0.1 @ 961127
c             0.2 @ 961209
c             0.3 @ 961210
c             0.4 @ 961212
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'PROF', vers = '080625/1.1.1')
c
      integer maxsiz, maxmsk
      parameter (maxsiz = maxgk1)
      parameter (maxmsk = maxgk1)
c
c      pointer (iaptr,mapa)
c      pointer (ibptr,mapb)
c      pointer (icptr,maska)
c      pointer (idptr,mapc)
c      pointer (ieptr,maskb)
c      pointer (ifptr,maskc)
c      pointer (igptr,mapd)
c
c      real mapa(1),mapb(1),mapc(1),mapd(1)
c      integer maska(1),maskb(1),maskc(1), malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr,idptr,ieptr,ifptr,igptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr,idptr,ieptr,ifptr,igptr
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
c
      nb = wrdbyt*masksize
      icptr = fmalloc (nb)
c
      minmsk = masksize / 10
      nb = wrdbyt*minmsk
      call jvalut (' Allocate mini maps/masks of size :',1,minmsk)
      idptr = fmalloc (nb)
      ieptr = fmalloc (nb)
      ifptr = fmalloc (nb)
      igptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0 .or.
     +    idptr .eq. 0 .or. ieptr .eq. 0 .or. ifptr .eq. 0 .or.
     +    igptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0 .or.
     +    idptr .le. 0 .or. ieptr .le. 0 .or. ifptr .le. 0 .or.
     +    igptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call doprof (%val(iaptr),%val(ibptr),%val(icptr),
     +             %val(idptr),%val(ieptr),%val(igptr),
     +             %val(ifptr),
     +             mapsize, masksize, minmsk)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
      call ffree (idptr)
      call ffree (ieptr)
      call ffree (ifptr)
      call ffree (igptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine doprof (mapa, mapb, maskx,
     +                   mapy, mapz, maskz, cntr,
     +                   maxsiz, maxmsk, minmap)
c
c ... Electron-density PROFiles
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz,maxatm,maxmsk,maxgrp,maxlib,minmap,maxapg,maxrt
      parameter (maxatm = 10000)
      parameter (maxgrp = 100)
      parameter (maxapg = 100)
      parameter (maxrt  = 1000)
      parameter (maxlib = maxgrp * maxapg)
c
      real mapa(maxsiz),mapb(maxsiz),mapz(minmap),mapy(minmap)
      real xyz(3,maxatm),libxyz(3,maxlib),shadow(3,maxapg)
      real grid(3),cella(6),cellb(6),rt(12,maxrt),rtsym(12)
      real spacb(3),rmsres(maxrt)
      real xdum,radius,rms,avy,avysq,avxy,avxmy
      real cc,r1,r2,rmsd,qdum
c
      integer maskx(maxmsk),maskz(minmap),cntr(minmap),bakptr(maxatm)
      integer grptr1(maxgrp),grptr2(maxgrp),copptr(maxrt),cobptr(maxrt)
      integer orgna(3),exta(3),grida(3),gridb(3),uvwa(3)
      integer orgnx(3),extx(3),orgnz(3),extz(3)
      integer numgrp,numlib,ig,ia,ib,ic,ncopy,io1,io2,io3
      integer spgrp,natoms,ierr,iunit,length,leng1
      integer i,j,k,idum
c
      character atmlin(maxatm)*80,liblin(maxlib)*80,grpnam(maxgrp)*80
      character file*80,line*256,pfile*80
c
      logical lmaska(maxlib)
      logical xinter,linter
c
      data rtsym /1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.0,0.0,0.0/
c
code ...
c
      iunit = 11
c
c ... radius around atoms (A)
c
      radius = 1.6
c
c ... cell and 0.5 A grid for profile map and masks
c
      do i=1,3
        cellb (i) = 100.0
        cellb (i+3) = 90.0
        gridb (i) = 200
        spacb (i) = cellb(i) / float(gridb(i))
      end do
c
      linter = xinter()
c
c ... 950118 - check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
      call jvalut (' Max nr of atoms in model  :',1,maxatm)
      call jvalut (' Max nr of profile groups  :',1,maxgrp)
      call jvalut (' Max nr of atoms per group :',1,maxapg)
      call jvalut (' Max nr of occurences      :',1,maxrt)
      call fvalut (' Atom masking radius (A)   :',1,radius)
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
      call edin (file, iunit, mapa,  orgna, exta, grida, uvwa,
     +           cella, spgrp, maxsiz)
c
      close (iunit)
      call telmap (grida,orgna,exta,cella)
c
      do i=1,3
        grid (i) = cella (i) / float (grida(i))
      end do
      call fvalut (' Grid spacing (A):',3,grid)
c
c ... read structure
c
      write (*,*)
      file = 'm1.pdb'
      call textin (' Name of model PDB file ?',file)
      call textut (' Name of model PDB file :',file)
c
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('Unable to open file')
        return
      end if
c
      call pdboni ('P',iunit,maxatm,natoms,xyz,atmlin,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading file')
        return
      end if
c
      if (natoms .lt. 3) then
        call errcon ('Less than three atoms')
        return
      end if
c
      do i=1,natoms
        write (atmlin(i)(61:66),'(f6.1)') -111.1
        bakptr (i) = -1
      end do
c
      close (iunit)
c
c ... read library (default to $GKLIB/prof.lib)
c
      write (*,*)
      file = 'prof.lib'
      call gklibf (file)
      call textin (' Name of PROF library file ?',file)
      call textut (' Name of PROF library file :',file)
c
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('Unable to open file')
        return
      end if
c
 6000 format (A)
c
      numgrp = 0
      numlib = 0
c
   10 continue
      read (iunit,6000,end=40,err=66) line
      if (line(1:6) .ne. 'PROFIL') goto 10
      numgrp = numgrp + 1
      grpnam(numgrp) = line(8:)
      grptr1 (numgrp) = numlib + 1
      grptr2 (numgrp) = -1
      call textut (' Group :',grpnam(numgrp))
c
c ... read this group
c
   20 continue
      read (iunit,6000,err=66) line
c
      if (line(1:6) .eq. 'ENDPRO') goto 30
c
c ... ATOM   -> use for LSQ and MASK
c     REMARK -> only use for LSQ
c
      if (line(1:6) .eq. 'ATOM  ') then
        numlib = numlib + 1
        liblin (numlib) = line
        read (line(31:54),'(3f8.3)') (libxyz(j,numlib),j=1,3)
        lmaska (numlib) = .true.
      else if (line(1:6) .eq. 'REMARK') then
        numlib = numlib + 1
        liblin (numlib) = line
        read (line(31:54),'(3f8.3)') (libxyz(j,numlib),j=1,3)
        lmaska (numlib) = .false.
      end if
      goto 20
c
c ... check if at least 3 atoms
c
   30 continue
      grptr2 (numgrp) = numlib
      idum = grptr2(numgrp)-grptr1(numgrp)+1
      if (idum .lt. 3) then
        numgrp = numgrp - 1
        call errcon ('Group has fewer than 3 atoms')
      else if (idum .gt. maxapg) then
        numgrp = numgrp - 1
        call errcon ('Group has too many atoms')
      end if
      goto 10
c
   66 continue
      call errcon ('While reading library file')
      return
c
   40 continue
      close (iunit)
c
      if (numgrp .lt. 1) then
        call errcon ('No groups found')
        return
      end if
c
      idum = grptr2(numgrp)-grptr1(numgrp)+1
      if (idum .le. 0) then
        grptr2 (numgrp) = numlib
        if ( idum .lt. 3) then
          numgrp = numgrp - 1
          call errcon ('Group has fewer than 3 atoms')
        else if (idum .gt. maxapg) then
          numgrp = numgrp - 1
          call errcon ('Group has too many atoms')
        end if
      end if
c
      write (*,*)
      file = 'prof.E'
      call textin (' Output CCP4 map file ?',file)
      call textut (' Output CCP4 map file :',file)
c
      write (*,*)
      pfile = 'prof.pdb'
      call textin (' Output PDB file ?',pfile)
      call textut (' Output PDB file :',pfile)
c
c ... calc origin and extent of mask/map around groups 
c
      call sizmsk (cellb,gridb,numlib,maxlib,libxyz,
     +             orgnx,extx,ierr,lmaska)
c
      write (*,*)
      call jvalut (' Origin :',3,orgnx)
      call jvalut (' Extent :',3,extx)
c
      idum = extx(1)*extx(2)*extx(3)
      call jvalut (' Calculated map/masksize :',1,idum)
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
        mapb (i) = 0.0
        maskx (i) = 0
      end do
c
c ... loop over groups
c
      do ig = 1, numgrp
c
        write (*,*)
        write (*,*) '*****************************************'
        write (*,*)
        call textut (' Group :',grpnam(ig))
c
        ncopy = 0
c
c ... set up minimap and minimask
c
        idum = grptr2(ig)-grptr1(ig)+1
        call sizmsk (cellb,gridb,idum,
     +               idum,libxyz(1,grptr1(ig)),orgnz,extz,ierr,lmaska)
c
c        call jvalut (' Origin :',3,orgnz)
c        call jvalut (' Extent :',3,extz)
c
        idum = extz(1)*extz(2)*extz(3)
c        call jvalut (' Calculated map/masksize :',1,idum)
c
        do i=1,idum
          mapz (i) = 0.0
          mapy (i) = 0.0
          maskz (i) = 0
          cntr (i) = 0
        end do
c
c ... generate mask around current group
c
        idum = grptr2(ig)-grptr1(ig)+1
        call grpmsk (cellb,gridb,idum,idum,radius,
     +               libxyz(1,grptr1(ig)),orgnz,extz,
     +               minmap,maskz,ierr,lmaska)
c        call jvalut (' Nr of atoms in mask  :',1,idum)
c
        j = 0
        idum = extz(1)*extz(2)*extz(3)
        do i=1,idum
          if (maskz(i) .eq. 1) j = j + 1
        end do
c        call jvalut (' Nr of points in mask :',1,j)
c
c ... loop over structure to find copies of this group
c
        ia = 0
  100   continue
        ia = ia + 1
        if (ia .gt. natoms) goto 199
c
c ... same residue type ?
c
        if (atmlin(ia)(18:20) .ne. liblin(grptr1(ig))(18:20)) then
          goto 100
        end if
c
c        call textut (' Check :',atmlin(ia)(1:30))
c
c ... find last atom of this residue
c
        do j=ia+1,natoms
          if (atmlin(j)(22:26) .ne. atmlin(ia)(22:26)) then
            ib = j - 1
            goto 110
          end if
        end do
        ib = natoms
  110   continue
        ic = ib
c
c ... get pointers to library atoms
c
        do j=grptr1(ig),grptr2(ig)
          idum = j-grptr1(ig)+1
          do k=ia,ib
            if (atmlin(k)(13:16) .eq. liblin(j)(13:16)) then
              shadow (1,idum) = xyz (1,k)
              shadow (2,idum) = xyz (2,k)
              shadow (3,idum) = xyz (3,k)
              if (lmaska(j)) bakptr (k) = j
              goto 120
            end if
          end do
          call errcon (' Atom not found in structure')
          call textut (' Atom    :',liblin(j)(1:30))
          call textut (' Residue :',atmlin(ia)(1:30))
          ia = ib
          goto 100
c
  120     continue
        end do
c
c ... add it (NCOPY <= MAXRT ?)
c
        ncopy = ncopy + 1
        copptr (ncopy) = ia
        cobptr (ncopy) = ic
c
c ... get RT-operator
c
        idum = grptr2(ig)-grptr1(ig)+1
        call lsqgjk (shadow,libxyz(1,grptr1(ig)),
     +               idum,rms,rt(1,ncopy),ierr)
c
        rmsres (ncopy) = rms
c
c        write (*,6010) atmlin(ia)(18:26),rms
c
c 6010 format (1x,a,' ... RMSD = ',f8.3,' A')
c
c        call fvalut (' RMSD :',1,rms)
c
c ... extract density around this copy
c
        call sub02x (
     +    mapa,exta(1),exta(2),exta(3),orgna,grid,cella,
     +    mapy,extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +    maskz, rtsym, 1, rtsym, 1,
     +    rt(1,ncopy),ierr,avy,avysq,avxy,avxmy,.false.)
c
        idum = extz(1)*extz(2)*extz(3)
        do i=1,idum
          if (maskz(i) .eq. 1) mapz (i) = mapz (i) + mapy (i)
        end do
c
        ia = ib
        goto 100
c
c ... end of loop over structure
c
  199   continue
c
c ... average profile if more than one copy
c
        if (ncopy .le. 0) goto 299
c
        xdum = 1.0 / float(ncopy)
        idum = extz(1)*extz(2)*extz(3)
        qdum = 0.0
        j = 0
        do i=1,idum
          if (maskz(i) .eq. 1) then
            mapz (i) = mapz (i) * xdum
            qdum = qdum + mapz (i)
            j = j + 1
          end if
        end do
c
        call jvalut (' Nr of map points :',1,j)
        qdum = qdum / float(j)
        call rvalut (' Average density  :',1,qdum)
c
c ... add to complete output map of all groups
c
        io1 = orgnz(1)-orgnx(1)
        io2 = orgnz(2)-orgnx(2)
        io3 = orgnz(3)-orgnx(3)
        call fillin (mapb,extx(1),extx(2),extx(3),
     +               mapz,extz(1),extz(2),extz(3),
     +               maskz,maskx,io1,io2,io3)
c
c ... compare to individual maps
c
        do ic=1,ncopy
c
          ia = copptr(ic)
c
c          call textut (' Check :',atmlin(ia)(1:30))
c
c ... extract density around this copy again
c
        call sub02x (
     +    mapa,exta(1),exta(2),exta(3),orgna,grid,cella,
     +    mapy,extz(1),extz(2),extz(3),orgnz,spacb,cellb,
     +    maskz, rtsym, 1, rtsym, 1,
     +    rt(1,ic),ierr,avy,avysq,avxy,avxmy,.false.)
c
c ... calc corr coeff with profile map
c
          call simsim (mapy, mapz, maskz, extz(1), extz(2), extz(3),
     +      cc,r1,r2,rmsd)
c
c          write (*,6020) atmlin(ia)(18:26),cc,r2,rmsd
c
c 6020 format (1x,a,' ... CC = ',f6.3,' ... R = ',f6.3,
c     +        ' ... RMSD = ',1pe12.4)
c
          write (*,6030) atmlin(ia)(18:26),cc,r2,rmsres(ic)
c
 6030 format (1x,a,' ... CC = ',f6.3,' ... R = ',f6.3,
     +  ' ... LSQ-RMSD = ',f6.3,' A')
c
c ... set B-factors to Correlation Coefficient
c
           do i=copptr(ic),cobptr(ic)
             if (bakptr(i) .gt. 0) then
               write (atmlin(i)(61:66),'(f6.2)') 100.0 * cc
               bakptr (i) = -1
             end if
           end do
c
c          call fvalut (' Correlation coefficient    :',1,cc)
c          call fvalut (' R-factor w.r.t. <profile>  :',1,r2)
c          call fvalut (' RMSD(rho) w.r.t. <profile> :',1,rmsd)
c
        end do
c
  299   continue
c
c ... end loop over groups
c
      end do
c
c ... adjust combined map
c
      idum = extx(1)*extx(2)*extx(3)
      do i=1,idum
        if (maskx(i) .gt. 1) then
          mapb (i) = mapb (i) / float(maskx(i))
        end if
      end do
c
c ... write map containing all groups
c
      write (*,*)
      spgrp = 1
      call edout (file,iunit,mapb,orgnx,extx,gridb,uvwa,cellb,spgrp)
c
c ... write new PDB file with B = CC (group dens, profile dens)
c
      close (iunit)
      call xopxua (iunit,pfile,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening output PDB file')
        return
      end if
c
      write (*,*)
      call prompt (
     +  ' Writing PDB file with B = CC (group dens, profile dens)')
      call stamp (line)
      write (iunit,6000,err=900) 'REMARK File: '//pfile(1:leng1(pfile))
      write (iunit,6000,err=900) 'REMARK '//line(1:leng1(line))
      write (iunit,6000,err=900)
     +  'REMARK B = 100.0 * CC (group dens, profile dens)'
      do i=1,natoms
        write (iunit,6000,err=900) atmlin(i)(1:leng1(atmlin(i)))
      end do
      write (iunit,6000,err=900) 'END   '
      close (iunit)
c
      return
c
  900 continue
      call errcon ('While writing output PDB file')
c
      return
      end
c
c
c
      subroutine pdboni (typ,iunit,maxatm,natoms,xyz,
     +  atmlin,ierr)
c
      implicit none
c
      integer iunit,maxatm,natoms,ierr,i,jerr,nopt,j,nh
c
      real xyz(3,maxatm)
c
      logical lhydro
c
      character atmlin(maxatm)*80
      character typ*1,line*80,optpar(4)*80
c
code ...
c
      ierr = -1
      natoms = 0
      nh = 0
      if (typ .eq. 'B') goto 1000
c
200   continue
        read (iunit,'(a)',err=6900,end=210) line
        if (line(1:6) .ne. 'ATOM  ' .and.
     +      line(1:6) .ne. 'HETATM') goto 200
        natoms = natoms+1
        if (natoms .gt. maxatm) then
          call errcon ('Too many atoms - rest skipped')
          call jvalut (' Max nr of atoms :',1,maxatm)
          ierr = 0
          return
        end if
c
        atmlin(natoms) = line
c
c ... strip hydrogen atoms
c
        if (lhydro(atmlin(natoms)(13:16))) then
          nh = nh + 1
          natoms = natoms - 1
          goto 200
        end if
c
c ... okay, read coordinates
c
        read (line,'(30x,3f8.3)')  (xyz(i,natoms),i=1,3)
c
      goto 200
c
cATOM     25  CB  ALA     5      -2.700   2.106  -0.485  1.00 20.00   6
c123456789012345678901234567890123456789012345678901234567890
c
210   continue
      call jvalut (' Number of atoms  :',1,natoms)
      call jvalut (' Nr of H stripped :',1,nh)
      ierr = 0
      return
c
1000  continue
c
      jerr = 0
c
      read (iunit,'(a)',err=6910) line
      call extrop (line,nopt,4,optpar,jerr)
      if (jerr .ne. 0) goto 6910
      if (nopt .lt. 4) goto 6910
c
      call textut (' Data block name :',optpar(1))
      call textut (' Data block type :',optpar(2))
      call textut (' Nr of values    :',optpar(3))
      call textut (' Format          :',optpar(4))
c
      call str2i (optpar(3),natoms,jerr)
      if (jerr .ne. 0) goto 6910
c
      natoms = natoms / 4
      if (natoms .gt. maxatm) then
        call errcon ('Too many atoms - rest skipped')
        call jvalut (' Max nr of atoms :',1,maxatm)
        natoms = maxatm
      end if
c
      read (iunit,optpar(4),err=6910) ((xyz(i,j),i=1,3), j=1,natoms)
c
      call jvalut (' Number of atoms :',1,natoms)
      ierr = 0
      return
c
c --- Read errors
c
 6900 continue
      call errcon ('While reading PDB file - operation aborted')
      return
c
 6910 continue
      call errcon ('While reading BONES file - operation aborted')
      return
c
      end
c
c
c
      subroutine sizmsk (cell,grid,ctmask,maxatm,xyzr,
     +                   origin,extent,ierr,lmaska)
c
      implicit none
c
      integer maxatm
      integer i,j,ctmask,extent(3),grid(3),origin(3),ct,ierr
c
      real cell(6),g(3),xyzr(3,maxatm)
      real b(3,3),ming(3),minr,maxg(3),maxr
      real x1(3),x2(3),minc(3),maxc(3)
c
      logical lmaska(maxatm)
c
code ...
c
      ierr = -1
c
      minr = 5.0
      maxr = 5.0
c
c --- Check the coords to see the required origin and extent
c
      call orthog (cell, b, 1)
c
      do i=1,3
        ming(i) =  9999999.0
        maxg(i) = -9999999.0
        minc(i) =  9999999.0
        maxc(i) = -9999999.0
      end do
c
      do 150 i=1,ctmask
        do 160 j=1,3
160       x1(j) = xyzr(j,i)
        call mulmtx (b, x1, x2, 3, 3, 1)
        do 170 j=1,3
c
          if (lmaska(i)) then
            minc (j) = min (minc(j),xyzr(j,i))
            maxc (j) = max (minc(j),xyzr(j,i))
          end if
c
          x1(j) = x2(j)*float(grid(j))
c
          if (lmaska(i)) then
            ming (j) = min (ming(j),x1(j))
            maxg (j) = max (maxg(j),x1(j))
          end if
c
170     continue
150   continue
c
c      call fvalut (' Lower bounds (coordinates) :',3,minc)
c      call fvalut (' Upper bounds (coordinates) :',3,maxc)
c      call fvalut (' Lower bounds (grid points) :',3,ming)
c      call fvalut (' Upper bounds (grid points) :',3,maxg)
c
      do i=1,3
        g(i) = cell(i)/float(grid(i))
        origin(i) = nint (ming(i)- minr/g(i))-1
        j = nint (maxg(i)+ maxr/g(i))+1
        extent(i) = j- origin(i)+1
      end do
c
c --- Now extend by some grid points around this origin and extent
c
      j = 5
      ct = 1
      do i=1,3
        extent(i) = extent(i) + 2*j
        origin(i) = origin(i) - j
        ct = ct * extent(i)
      end do
c
c --- Report
c
c      call jvalut (' Mask origin :',3,origin)
c      call jvalut (' Mask extent :',3,extent)
c      call jvalut (' Grid points :',1,ct)
c      call jvalut (' Mask grid   :',3,grid)
c      call fvalut (' Mask cell   :',6,cell)
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
      subroutine grpmsk (cell,grid,ctmask,maxatm,radius,
     +                   xyzr,origin,extent,
     +                   size,mask,ierr,lmaska)
c
      implicit none
c
      integer maxatm
c
      integer i,j,k,l,ctmask,extent(3),grid(3),origin(3)
      integer ct,i1,i2,i3,off,ierr,size,inx,nxy
      integer mask(size)
c
      real cell(6),g(3),distce,x(3),xp(3),xyzr(3,maxatm)
      real a(3,3),b(3,3),c(3,3),fake(6)
      real x1(3),x2(3),gmin
      real radius
c
      logical lmaska(maxatm)
c
code ...
c
      ierr = -1
c
      do i=1,3
        g(i) = cell(i)/float(grid(i))
      end do
c
c --- Blank mask
c
      ct = extent(3)*extent(2)*extent(1)
      do k=1,ct
        mask(k) = 0
      end do
c
      do i=1,3
        fake(i) = 1.0
        fake(i+3) = cell(i+3)
      end do
      call orthog (fake, a, 1)
      call orthog (fake, c, 0)
      call orthog (cell, b, 1)
c
c --- Do the actual mask job
c
      nxy = extent(1)*extent(2)
      gmin = min (g(1),g(2),g(3))
c
      do 120 i=1,ctmask
        if (.not. lmaska(i)) goto 125
        do j=1,3
          x(j) = xyzr(j,i)
        end do
        call mulmtx (a, x, x1, 3, 3, 1)
        off = 1 + radius/gmin
        do 130 j=-off,off	
          do 130 k=-off,off	
            do 130 l=-off,off
              xp(1) = x1(1)+ l*g(1)
              xp(2) = x1(2)+ k*g(2)
              xp(3) = x1(3)+ j*g(3)
              call mulmtx (c, xp, x2, 3, 3, 1)
              if (distce(x,x2) .gt. radius) goto 130
              i1 = nint(xp(1)/g(1))- origin(1) +1
              i2 = nint(xp(2)/g(2))- origin(2) +1
              i3 = nint(xp(3)/g(3))- origin(3) +1
              if (i1 .le. 0) goto 130
              if (i2 .le. 0) goto 130
              if (i3 .le. 0) goto 130
              if (i1 .gt. extent(1)) goto 130
              if (i2 .gt. extent(2)) goto 130
              if (i3 .gt. extent(3)) goto 130
              inx = (i3-1)*nxy + (i2-1)*extent(1) + i1
              mask(inx) = 1
 130    continue
 125    continue
 120  continue
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
      subroutine simsim (mapy, mapz, maskz, ext1, ext2, ext3,
     +      ccoef,sumfo,sumfc,rmsd)
c
      implicit none
c
      integer ext1, ext2, ext3, ntot, i1,i2,i3
      integer maskz(ext1, ext2, ext3)
c
      real mapy(ext1, ext2, ext3),mapz(ext1, ext2, ext3)
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
            if (maskz(i1,i2,i3) .eq. 1) then
              ntot = ntot + 1
              avalue = mapy(i1,i2,i3)
              bvalue = mapz(i1,i2,i3)
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
      subroutine fillin (mapb,extx1,extx2,extx3,
     +                   mapz,extz1,extz2,extz3,
     +                   maskz,maskx,io1,io2,io3)
c
      implicit none
c
      integer extx1,extx2,extx3,extz1,extz2,extz3,io1,io2,io3
c
      real mapb(extx1,extx2,extx3),mapz(extz1,extz2,extz3)
c
      integer maskx(extx1,extx2,extx3),maskz(extz1,extz2,extz3)
      integer i1,i2,i3
c
code ...
c
      do i1=1,extz1
        do i2=1,extz2
          do i3=1,extz3
            if (maskz(i1,i2,i3) .eq. 1) then
              mapb (i1+io1,i2+io2,i3+io3) = 
     +          mapb (i1+io1,i2+io2,i3+io3) + mapz (i1,i2,i3)
              maskx (i1+io1,i2+io2,i3+io3) = 
     +          maskx (i1+io1,i2+io2,i3+io3) + 1
            end if
          end do
        end do
      end do
c
      return
      end
c
c
c
      subroutine sub02x (
     +    mapa,exta1,exta2,exta3,orgna,spaca,cella,
     +    mapb,extb1,extb2,extb3,orgnb,spacb,cellb,
     +    maskb, rtbtoa, ctrt, rtsym,ctsym, or2or, ierr,
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
      real mapa(exta1, exta2, exta3),spaca(3),cella(6)
c
      integer orgnb(3), extb1, extb2, extb3
      real mapb(extb1, extb2, extb3),spacb(3),cellb(6)
      integer maskb(extb1, extb2, extb3)
c
      real rtbtoa(12,*),rtsym(12,*),or2or(12),val1,f,cc
      real avy(*),avysq(*),avxy(*),avxmy(*),avx,avxsq,avax
      integer ctrt,ctsym,ierr,nbit,ngjk
c
      integer errcod, i, j, k, l, loop, m, ext(3), ijk(3)
c      integer bobo1,bobo2,bobo3
      integer nerr1,nerr2
c
      real forgna(3), fexta(3), mapbit(4,4,4)
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
c      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
c      call jvalut (' Points in mask :',1,bobo1)
c      bobo2 = bobo1/10
c      bobo3 = bobo2
c
c ... loop over the mask points
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        mapb (i,j,k) = 0.0
        if (maskb(i,j,k) .ne. 1) goto 100
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
          call vecrtv (cmska,ncmska,1,rtbtoa(1,loop),rtbtoa(10,loop))
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
            call bldbit (mapa,exta1,exta2,exta3,orgna,spaca,cella,
     +                   forgna,fexta,rtsym,ctsym,mapbit,fcmska,errcod)
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
            call intrpl (mapbit, 4, 4, 4, x, value, errcod)
c
c ... if first FRCSYM worked okay, do straightforward interpolation
c
          else
c
            do l=1,3
              x(l) = fcmska(l)*cella(l)/spaca(l)- float(orgna(l)) + 1
            end do
c
            call intrpl (mapa, exta1, exta2, exta3, x, value, errcod)
c
          end if
c
          if (errcod .eq. 0) then
            mapb (i,j,k) = mapb (i,j,k) + value
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
      call jvalut (' Calls to BLDBIT      :',1,nbit)
      call jvalut (' Severe FRCSYM errors :',1,nerr1)
      call jvalut (' Interpolation errors :',1,nerr2)
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
          if (maskb(i,j,k) .eq. 1) then
            if (mapb(i,j,k) .gt. 0.0) then
              mapb(i,j,k) = mapb(i,j,k) * xdum
            else
              mapb(i,j,k) = small
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
200         mapb(i,j,k) = mapb(i,j,k) * xdum
        end if
      end if
c
      return
      end
