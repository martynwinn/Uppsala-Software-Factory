      program solex
c
c ... SOLution EXtractor for use with ESSENS
c
c ... Gerard Kleywegt @ 961120
c
c ... f77sgi solex.f
c     f77 -o SOLEX solex.o ~/progs/gklib/kleylib /public/src/ccp4/lib/libccp4.a ; strip SOLEX 
c
c ... f77al solex.f
c     f77 -o SOLEX solex.o ~/progs/gklib/alpha_kleylib /public/src/ccp4/lib/src/libccp4.a ; strip SOLEX 
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'SOLEX', vers = '050622/2.0.2')
c
      integer maxsiz, maxmsk
      parameter (maxsiz = maxgk1)
      parameter (maxmsk = maxgk1)
c
c      pointer (iaptr,mapa)
c      pointer (ibptr,maskb)
c
c      real mapa(1)
c      integer maskb(1), malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr
      integer fmalloc
#endif
c
      integer nb, mapsize, masksize
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
c
      nb = wrdbyt*masksize
      ibptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call dosole (%val(iaptr),%val(ibptr),mapsize,masksize)
c
      call ffree (iaptr)
      call ffree (ibptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine dosole (mapa, maskx, maxsiz, maxmsk)
c
c ... SOLution EXtractor for use with ESSENS
c
c ... Gerard Kleywegt @ 961120
c
c ... f77sgi solex.f
c     f77 -o SOLEX solex.o ~/progs/gklib/kleylib /public/src/ccp4/lib/libccp4.a ; strip SOLEX 
c
c ... f77al solex.f
c     f77 -o SOLEX solex.o ~/progs/gklib/alpha_kleylib /public/src/ccp4/lib/src/libccp4.a ; strip SOLEX 
c
      implicit none
c
      include 'maxdim.incl'
c
      integer    maxsiz,maxatm,maxnew,maxmsk
c      parameter (maxsiz =    maxgk1)
      parameter (maxatm =     10000)
      parameter (maxnew = 10*maxatm)
c
      real mapa(maxsiz)
      real xyz(3,maxatm), rotxyz(3,maxatm)
      real newxyz(3,maxnew),newb(maxnew)
      real cell(6),grid(3),xt(3),cella(6),cellx(6)
      real xdum,q,xcut,xave,xsdv
c
      integer maskx(maxmsk),resid(maxatm),newptr(maxnew)
      integer orgna(3),exta(3),grida(3),uvwa(3)
      integer orgnx(3),extx(3),gridx(3)
      integer spgrp,natoms,ierr,iunit,junit,kunit,length
      integer i,j,mode,isum,jsum,idum,nmax
c
      character atmlin(maxatm)*80
      character file*80,line*80,base*80,dejavu*80
c
      logical xinter,lokay,linter,leuler,lalpha
c
code ...
c
      iunit = 11
      junit = 12
      kunit = 13
c
      linter = xinter()
c
c      call gainit (prognm,vers)
c
c ... 950118 - check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
c      call jvalut (' Max size of map :',1,maxsiz)
      call jvalut (' Max nr of atoms :',1,maxatm)
c
c ... ESSENS score map
c
      write (*,*)
      file = ' '
      call textin (' Input ESSENS score map file ?',file)
      call textut (' Input ESSENS score map file :',file)
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
c ... read in map
c
      call edin (file, iunit, mapa,  orgna, exta, grida, uvwa,
     $           cella, spgrp, maxsiz)
c
      close (iunit)
      call telmap (grida,orgna,exta,cella)
c
      do i=1,3
        grid (i) = cella (i) / float (grida(i))
      end do
      call fvalut (' Grid spacing (A):',3,grid)
c
      do i=1,6
        cell (i) = cella (i)
      end do
c
c ... ESSENS rotation file
c
      write (*,*)
      file = ' '
      call textin (' Input ESSENS rotation file ?',file)
      call textut (' Input ESSENS rotation file :',file)
c
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
c ... read rotation file
c
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening file')
        return
      end if
c
      read (iunit,'(A)',err=999,end=999) line
      if (line (1:19) .ne. 'ROTATION_MAP_ESSENS') then
        call errcon ('Not an ESSENS rotation file')
        return
      end if
c
      read (iunit,'(A)',err=999,end=999) line
      if (line(1:5) .eq. 'EULER') then
        leuler = .true.
      else if (line(1:5) .eq. 'POLAR') then
        leuler = .false.
      else
        call errcon ('No EULER or POLAR keyword found')
        return
      end if
c
      read (iunit,*,err=999,end=999) xave,xsdv
      read (iunit,*,err=999,end=999) (orgnx(i),i=1,3)
      read (iunit,*,err=999,end=999) (extx(i),i=1,3)
      read (iunit,*,err=999,end=999) (gridx(i),i=1,3)
      read (iunit,*,err=999,end=999) (cellx(i),i=1,6)
c
      idum = extx(1)*extx(2)*extx(3)
      if ( idum .gt. maxmsk) then
        call errcon ('Rotation file too big')
        call jvalut (' Available :',1,maxmsk)
        call jvalut (' Requested :',1,idum)
        return
      end if
c
c ... check if map and rotation file have same cell and grid
c
      lokay = .true.
c
      do i=1,3
        lokay = (lokay .and. (orgna(i) .eq. orgnx(i)))
        lokay = (lokay .and. (exta(i)  .eq. extx(i)))
        lokay = (lokay .and. (grida(i) .eq. gridx(i)))
        lokay = (lokay .and. (abs(cella(i)-cellx(i)) .le. 0.01))
        lokay = (lokay .and.
     +           (abs(cella(i+3)-cellx(i+3)) .le. 0.01))
      end do
c
      if (.not. lokay) then
        call errcon (
     +    'Map and rotation file do NOT have the same cell and grid !')
        return
      end if
c
      call prompt (' Reading rotation file ...')
      call slurpm (iunit,maskx,idum,ierr,isum)
      if (ierr .ne. 0) goto 999
c
  997 continue
      read (iunit,'(A)',err=999,end=998) line
      if (line (1:8) .ne. 'CHECKSUM') goto 997
c
      read (line(9:),*,err=999,end=999) jsum
      if (isum .ne. jsum) then
        call jvalut (' Checksum in file    :',1,jsum)
        call jvalut (' Checksum calculated :',1,isum)
        call errcon (' Checksums do not match; file corrupted ?')
        return
      end if
c
      close (iunit)
      call telmap (gridx,orgnx,extx,cellx)
      go to 10
c
  998 continue
      call errcon ('Checksum line not found')
      return
c
  999 continue
      call errcon ('While reading rotation file')
      return
c
   10 continue
c
c ... read template/search probe
c
      write (*,*)
      file = 'alpha.pdb'
      call textin (' Name of template PDB file ?',file)
      call textut (' Name of template PDB file :',file)
c
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('Unable to open file')
        return
      end if
c
      call pdboni ('P',iunit,maxatm,natoms,xyz,resid,atmlin,ierr)
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
      close (iunit)
c
c ... get centre-of-gravity
c
      do i=1,3
        xt (i) = 0.0
      end do
c
      do i=1,natoms
        do j=1,3
          xt(j) = xt(j) + xyz(j,i)
        end do
      end do
      do i=1,3
        xt(i)=xt(i)/float(natoms)
      end do
      call fvalut (' Centre of gravity of atoms :',3,xt)
c
      xdum = 9999.999
      idum = -1
      do i=1,natoms
        if (atmlin(i)(13:16) .eq. ' CA ') then
          q = (xyz(1,i)-xt(1))**2 + (xyz(2,i)-xt(2))**2 +
     +        (xyz(3,i)-xt(3))**2
          if (q .lt. xdum) then
            idum = i
            xdum = q
          end if
        end if
      end do
c
      if (idum .gt. 0) then
        call prompt (' Using central CA atom as pivot point')
        line = atmlin (idum)
        call pretty (line)
        call textut (' Atom :',line)
        xt(1) = xyz(1,idum)
        xt(2) = xyz(2,idum)
        xt(3) = xyz(3,idum)
      else
        call prompt (' Using centre-of-gravity as pivot point')
      end if
c
      call fvalut (' Pivot point :',3,xt)
c
      xdum = 0.0
      do i=1,natoms
        q = 0.0
        do j=1,3
          xyz(j,i) = xyz(j,i) - xt(j)
          q = q + (xyz(j,i))**2
        end do
        xdum = max (xdum, sqrt(q))
      end do
      call fvalut (' Furthest atom (A):',1,xdum)
c
c ... mode ?
c
      write (*,*)
      write (*,*) 'SOLEX can extract top solutions in one of 3 modes:'
      write (*,*) '1 = each to a separate file (e.g., domain)'
      write (*,*) '2 = all to one file (e.g., ligand)'
      write (*,*) '3 = try to trace alpha-helices'
      write (*,*) '4 = try to trace beta-strands'
ccc      call prompt (' ==> Mode 3 not implemented yet !!!')
      mode = 1
      call ivalin (' Solution extraction mode ?',1,mode)
      call ivalut (' Solution extraction mode :',1,mode)
c
c      call gkuser (line)
c      if (mode .eq. 3) then
c        if (line(1:6) .ne. 'gerard') then
c          call errcon ('Not implemented yet; using mode 2')
c          mode = 2
c        end if
c      end if
c
      if (mode .lt. 1 .or. mode .gt. 4) then
        call errcon ('Invalid value for mode')
        return
      end if
c
      lalpha = .false.
      if (mode .eq. 3) then
        lalpha = .true.
      else if (mode .eq. 4) then
        mode = 3
        lalpha = .false.
      end if
c
      write (*,*)
      if (mode .eq. 1) then
        base = 'solex_'
        call textin (' Basename for output PDB files ?',base)
        call textut (' Basename for output PDB files :',base)
      else
        if (mode .eq. 2) then
          base = 'solex.pdb'
        else if (mode .eq. 3) then
          base = 'beta.pdb'
          if (lalpha) base = 'alpha.pdb'
        end if
        call textin (' Name for output PDB file ?',base)
        call textut (' Name for output PDB file :',base)
      end if
c
c ... prepare for trace mode
c
      if (mode .eq. 3) then
        write (*,*)
        dejavu = 'beta.sse'
        if (lalpha) dejavu = 'alpha.sse'
        call textin (' Name for output DEJAVU SSE file ?',dejavu)
        call textut (' Name for output DEJAVU SSE file :',dejavu)
c
        write (*,*)
        call prompt (' Renaming " CA?" atoms to " CA " ...')
        do i=1,natoms
          if (atmlin(i)(13:15) .eq. ' CA') atmlin(i)(13:16) = ' CA '
        end do
        call prompt (' Removing all but " CA " atoms ...')
        j = 0
        do i=1,natoms
          if (atmlin(i)(13:16) .eq. ' CA ') then
            j = j + 1
            if (j .lt. i) then
              atmlin (j) = atmlin (i)
              resid (j) = resid (i)
              xyz(1,j) = xyz(1,i)
              xyz(2,j) = xyz(2,i)
              xyz(3,j) = xyz(3,i)
            end if
          end if
        end do
        call ivalut (' Nr of remaining atoms :',1,j)
        if (j .lt. 2) then
          call errcon ('Less than 2 CA atoms')
          return
        end if
        natoms = j
      end if
c
c ... cut-off ?
c
      write (*,*)
      xcut = xave + 2.0*xsdv
      call rvalin (' Cut-off score for top solutions ?',1,xcut)
      call rvalut (' Cut-off score for top solutions :',1,xcut)
c
c ... max nr of solutions ?
c
      write (*,*)
      nmax = 50
      call ivalin (' Max nr of top solutions ?',1,nmax)
      call ivalut (' Max nr of top solutions :',1,nmax)
      if (nmax .lt. 1) then
        call errcon ('Really funny, mate !')
        return
      end if
c
      if (mode.eq.3 .and. (nmax*natoms) .gt. maxnew) then
        call errcon ('Value too high')
        nmax = (maxnew / natoms) - 1
        call ivalut (' Value reset to :',1,nmax)
      end if
c
      call sub97x (mapa, maskx, exta(1), exta(2), exta(3), 
     +  orgna, cell, grid, natoms, xyz, rotxyz, xcut, nmax, 
     +  iunit, base, mode, maxnew, newxyz, leuler, linter,
     +  xave, xsdv, atmlin, newb, newptr, lalpha, dejavu)
c
      return
c
      end
c
c
c
      subroutine sub97x 
     + (mapa, maskx, exta1, exta2, exta3, orgna, 
     +  cell, grid, natoms, xyz, rxyz, xcut, nmax,
     +  iunit, base, mode, maxnew, newxyz, leuler, linter,
     +  xave, xsdv, atmlin, newb, newptr, lalpha, dejavu)
c
c --- do it
c
      implicit none
c
      integer orgna(3), exta1, exta2, exta3, natoms
      integer maxnew, mode, nmax
c
      real mapa(exta1, exta2, exta3)
      real xyz(3,natoms), rxyz(3,natoms), newxyz(3,maxnew)
      real a(3,3), b(3,3), newb(maxnew)
      real cell(6), grid(3), rt(9), del(3), x1(3), x(3)
c
      real total,user,sys,q,xcut,xave,xsdv,vsol
c
      integer maskx(exta1, exta2, exta3)
      integer newptr(maxnew),ext(3)
c
      integer i,i1,i2,i3,j,l,k,newatm,ipart
      integer msol,iunit,leng1,ierr,isol,jsol,ksol,nsol
c
      logical leuler,linter,lalpha
c
      character line*256
      character base*(*), dejavu*(*)
      character atmlin(natoms)*80
c
code ...
c
      write (*,*)
      call prompt (' Extracting solutions ...')
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell, a, 0)
      call orthog (cell, b, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
c
      nsol = 0
      newatm = 0
      q = xcut - 0.01 * abs(xcut) - 1.0
c
 4323 continue
      nsol = nsol + 1
c
      if (nsol .gt. nmax) then
        call prompt (' Max number of solutions done ...')
        nsol = nsol - 1
        goto 4399
      end if
c
      isol = 0
      jsol = 0
      ksol = 0
      msol = 0
      vsol = q
c
      do k=2,exta3-1
        do j=2,exta2-1
          do i=2,exta1-1
c
c ... is it an evaluated point ?
c
            if (maskx(i,j,k) .eq. -1) goto 4324
c
c ... is it a high value ?
c
            if (mapa(i,j,k) .le. vsol) goto 4324
c
c ... is it a local maximum in the original score map ?
c
            do i1=-1,1
              do i2=-1,1
                do i3=-1,1
                  if (i1.eq.0 .and. i2.eq.0 .and. i3.eq.0)
     +              goto 4326
                  if (mapa(i,j,k) .le. mapa(i+i1,j+i2,k+i3))
     +              goto 4324
 4326             continue
                end do
              end do
            end do
c
c ... all conditions satisfied
c
            vsol = mapa(i,j,k)
            isol = i
            jsol = j
            ksol = k
            msol = maskx(i,j,k)
c
 4324       continue
          end do
        end do
      end do
c
      if (vsol .le. xcut) then
        call prompt (' Next solution scores lower than cut-off')
        write (*,6100) nsol,vsol,isol,jsol,ksol,i1,i2,i3
        nsol = nsol - 1
        goto 4399
      end if
c
      call packut (i1,i2,i3,ipart,msol)
      write (*,6100) nsol,vsol,isol,jsol,ksol,i1,i2,i3
c
 6100 format (' Sol # ',i6,' Value ',1pe12.4,0p,' @ ',3i6,
     +  ' Rot ',3i6)
 6200 format ('REMARK Sol # ',i6,' Value ',1pe12.4,0p,' @ ',3i6,
     +  ' Rot ',3i6)
c
c ... mark this one as done
c
      maskx (isol,jsol,ksol) = -1
c
c ... write to file ?
c
      if (mode .eq. 1) then
        write (line,'(a,i8,a)') base(1:leng1(base)),nsol,'.pdb'
        call remspa (line)
        call xopxua (iunit,line,linter,ierr)
        call stamp (line)
        write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        write (iunit,6200) nsol,vsol,isol,jsol,ksol,i1,i2,i3
      else if (mode .eq. 2) then
        if (nsol .eq. 1) then
          call xopxua (iunit,base,linter,ierr)
          call stamp (line)
          write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        end if
        write (iunit,6200) nsol,vsol,isol,jsol,ksol,i1,i2,i3
      end if
c
      del(1) = float(i1)
      del(2) = float(i2)
      del(3) = float(i3)
c
      if (leuler) then
        call ccpeul (del,rt)
      else
        call ccppol (del,rt)
      end if
c
c ... generate best-fitting orientation/position
c
      do i=1,natoms
c ... rotate
        call mulmtx (rt,xyz(1,i),rxyz(1,i),3,3,1)
c ... fractionalise
        call mulmtx (b,rxyz(1,i),x1,3,3,1)
c ... convert to grid points, add best offset and add origin
        x(1) = isol + orgna(1) - 1 + x1(1)*cell(1)/grid(1)
        x(2) = jsol + orgna(2) - 1 + x1(2)*cell(2)/grid(2)
        x(3) = ksol + orgna(3) - 1 + x1(3)*cell(3)/grid(3)
c ... back to fractional
        do l=1,3
          x1(l) = x(l)*grid(l)/cell(l)
        end do
c ... orthogonalise
        call mulmtx (a,x1,rxyz(1,i),3,3,1)
      end do
c
      if (mode .eq. 1 .or. mode .eq. 2) then
c
c ... TO DO ? If mode = 2, also increase residue numbers ???
c
        do i=1,natoms
          write (atmlin(i)(31:66),'(3f8.3,2f6.2)')
     +      (rxyz(j,i),j=1,3),float(nsol)/100.0,(vsol-xave)/xsdv
          write (iunit,'(a)') atmlin(i)(1:leng1(atmlin(i)))
        end do
c
        if (mode .eq. 1) then
          write (iunit,'(a)') 'END'
          close (iunit)
        end if
c
      else if (mode .eq. 3) then
c
        do i=1,natoms
          newatm = newatm + 1
          newb   (newatm) = (vsol-xave)/xsdv
          newptr (newatm) = i
          newxyz (1,newatm) = rxyz (1,i)
          newxyz (2,newatm) = rxyz (2,i)
          newxyz (3,newatm) = rxyz (3,i)
        end do
c
      end if
c
c ... blank out solution ?
c
      goto 4323
c
c ... done
c
 4399 continue
c
      if (mode .eq. 2 .and. nsol .gt. 0) then
        write (iunit,'(a)') 'END'
        close (iunit)
      end if
c
      if (mode .eq. 1 .or. mode .eq. 1) then
        write (*,*)
        if (nsol .le. 0) then
          call prompt (' No solutions found')
        else if (nsol .eq. 1) then
          call prompt (' Top solution written')
        else
          call prompt (' Top solutions written')
        end if
      end if
c
      if (mode .eq. 3) then
        call jvalut (' Nr of CA atoms generated :',1,newatm)
        if (newatm .gt. natoms) then
          call prompt (' Trying to connect them ...')
          call tracer (newatm,natoms,newb,newptr,newxyz,atmlin,
     +      iunit,base,linter,nsol,maxnew,lalpha,dejavu)
        end if
      end if
c
      call gkdcpu (total,user,sys)
      write (*,'(/a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      return
      end
c
c
c
      subroutine pdboni (typ,iunit,maxatm,natoms,xyz,resid,
     +  atmlin,ierr)
c
      implicit none
c
      integer iunit,maxatm,natoms,ierr,i,jerr,nopt,j,nh
      integer resid(maxatm)
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
        read (line,'(22x,i4)') resid(natoms)
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
      subroutine slurpm (iunit,mask,nsize,ierr,nsum)
c
      implicit none
c
      integer iunit,nsize,ierr,i,nsum
      integer mask(*)
c
code ...
c
      ierr = 0
      nsum = 0
c
      read (iunit,*,err=999,end=999) (mask(i),i=1,nsize)
c
      mask (nsize+1) = 0
      do i=1,nsize,2
        nsum = nsum + mask(i) - mask(i+1)
      end do
c
      return
c
  999 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine tracer (newatm,natoms,newb,newptr,newxyz,atmlin,
     +                   iunit,base,linter,nsol,maxnew,lalpha,dejavu)
c
      implicit none
c
      integer maxsol,maxlen
      parameter (maxsol = 2000)
      parameter (maxlen = 2000)
c
      integer newatm,natoms,maxnew
c
      real newb(newatm),newxyz(3,maxnew)
      real mxcaca,partol
      real q,dist,rmsd,qq
      real vector(3,maxsol),angvec(maxsol,maxsol),closes(maxsol,maxsol)
      real nitocj(maxsol,maxsol)
c
      integer numatm(maxsol),ndummy(maxsol)
      integer pointer(maxlen,maxsol)
      integer newptr(newatm)
      integer iunit,nsol,leng1,junit
      integer i,j,k,l,m,i1,j1,i2,j2,kk,ll,idum
c
      logical ldone(maxsol)
      logical linter,lextend,lalpha
c
      character atmlin(natoms)*80,base*(*),line*128,dejavu*(*)
      character ssenam*6,ssetyp*6,quote*1
c
      data quote /''''/
c
code ...
c
c ... max distance between "identical" CA atoms
c
      mxcaca = 1.5
c
c ... tolerance for parallel vectors
c
      partol = 30.0
c
      write (*,*)
      call prompt (' Starting TRACER ...')
c
      call ivalut (' Nr of solutions :',1,nsol)
      call ivalut (' Max allowed     :',1,maxsol)
c
      if (nsol .gt. maxsol) then
        call errcon ('Too many potential solutions')
        call ivalut (' Only use top :',1,maxsol)
        nsol = maxsol
      end if
c
      do i=1,nsol
        ldone (i) = .false.
        numatm (i) = natoms
        i1 = natoms*(i-1)
        do j=1,natoms
          pointer (j,i) = i1 + j
        end do
      end do
c
c ... main loop
c
   10 continue
c
c ... remove redundant solutions
c
      write (*,*)
      goto 1100
      call prompt (' Remove redundant solutions ...')
c
      do i=1,nsol
        if (ldone(i)) goto 1030
        do j=1,nsol
          if (ldone(j)) goto 1020
          if (numatm(j) .gt. numatm(i)) goto 1020
          if (i .eq. j) goto 1020
c
          do l=0,(numatm(i)-numatm(j))
            rmsd = 0.0
            do k=1,numatm(j)
              q = dist (pointer(1,i)+k+l-1,pointer(1,j)+k-1,newxyz)
              rmsd = rmsd + q*q
              if (q .gt. mxcaca) goto 1015
            end do
            rmsd = sqrt (rmsd/float(numatm(j)))
            write (*,*) 'Solution ',j,' = ',i,' forward, rmsd: ',rmsd
            ldone(j) = .true.
            goto 1020
c
 1015       continue
          end do
c
          do l=0,(numatm(i)-numatm(j))
            rmsd = 0.0
            do k=1,numatm(j)
              q = dist (pointer(1,i)+k+l-1,
     +                  pointer(numatm(j),j)-k+1,newxyz)
              rmsd = rmsd + q*q
              if (q .gt. mxcaca) goto 1025
            end do
            rmsd = sqrt (rmsd/float(numatm(j)))
            write (*,*) 'Solution ',j,' = ',i,' backward, rmsd: ',rmsd
            ldone(j) = .true.
            goto 1020
c
 1025       continue
          end do
c
 1020     continue
        end do
 1030   continue
      end do
c
c ... calculate vectors, angles and distances
c
 1100 continue
      call prompt (' Calculate vectors, angles and distances ...')
c
      do i=1,nsol
        if (.not. ldone (i)) then
          i1 = pointer (1,i)
          i2 = pointer (numatm(i),i)
c          print *,' vector ',i,i1,i2
          vector (1,i) = newxyz(1,i2) - newxyz(1,i1)
          vector (2,i) = newxyz(2,i2) - newxyz(2,i1)
          vector (3,i) = newxyz(3,i2) - newxyz(3,i1)
c          call fvalut (' Vector :',3,vector(1,i))
        end if
      end do
c
      do i=1,nsol
        if (ldone(i)) goto 1130
        i1 = pointer (1,i)
        i2 = pointer (numatm(i),i)
        do j=1,nsol
          if (ldone(j)) goto 1120
          if (i.eq.j) goto 1120
          j1 = pointer (1,j)
          j2 = pointer (numatm(j),j)
c
c ... angle of direction vectors
c
          call vecang (vector(1,i),vector(1,j),q,k)
          angvec (i,j) = q
          angvec (j,i) = q
c          print *,i,j,' vecang'
c
c ... distance of closest approach
c
          q = 999.99
          do k=1,numatm(i)
            do l=1,numatm(j)
              kk = pointer (k,i)
              ll = pointer (l,j)
              q = min (q, dist (kk,ll,newxyz))
            end do
          end do
          closes (i,j) = q
          closes (j,i) = q
c          print *,i,j,' closes'
c
c ... distance N(i) -> C(j)
c
          nitocj (i,j) = dist (i1,j2,newxyz)
          nitocj (j,i) = dist (i2,j1,newxyz)
c          print *,i,j,' nitocj'
c
 1120     continue
        end do
 1130   continue
      end do
c
c ... try to extend
c
      call prompt (' Try to extend bits ...')
c
      lextend = .false.
c
      do i=1,nsol-1
        if (ldone(i)) goto 1200
        i1 = pointer (1,i)
        i2 = pointer (numatm(i),i)
        do j=i+1,nsol
          if (ldone(j)) goto 1210
          j1 = pointer (1,j)
          j2 = pointer (numatm(j),j)
c
 1230     continue
          if (angvec(i,j) .le. partol) then
            if (closes(i,j) .le. mxcaca) then
              if (nitocj(i,j) .gt. nitocj(j,i)) then
                goto 2000
              else
                goto 3000
              end if
            end if
          else if ( (180.0-angvec(i,j)) .le. partol) then
            if (closes(i,j) .le. mxcaca) then
              goto 4000
            end if
          end if
          goto 1210
c
c === APPEND J AT THE C TERMINUS OF I
c
 2000     continue
          print *,'Append ',j,' at C-terminus of ',i
ccc          print *,angvec(i,j),closes(i,j),nitocj(i,j),nitocj(j,i)
c
c ... find atom closest to C-term
c
          qq = mxcaca
          kk = 0
          do k=1,numatm(j)
            q = dist (pointer(numatm(i),i),pointer(k,j),newxyz)
            if (q .lt. qq) then
              qq = q
              kk = k
            end if
          end do
c
          if (kk .gt. 0 .and. kk .lt. numatm(j)) then
            do k=kk+1,numatm(j)
              numatm(i) = numatm(i) + 1
              pointer (numatm(i),i) = pointer (k,j)
            end do
          end if
c
ccc          print *,i,' now has ',numatm(i),' atoms'
          ldone (j) = .true.
          lextend = .true.
          goto 1210
c
c === PREPEND J AT THE N TERMINUS OF I
c
 3000     continue
          print *,'Prepend ',j,' at N-terminus of ',i
ccc          print *,angvec(i,j),closes(i,j),nitocj(i,j),nitocj(j,i)
c
c ... find atom closest to N-term
c
          qq = mxcaca
          kk = 0
          do k=1,numatm(j)
            q = dist (pointer(1,i),pointer(k,j),newxyz)
            if (q .lt. qq) then
              qq = q
              kk = k
            end if
          end do
c
          if (kk .gt. 1) then
            idum = 0
            do k=1,kk-1
              idum = idum + 1
              ndummy (idum) = pointer (k,j)
            end do
c
            do k=1,numatm(i)
              idum = idum + 1
              ndummy (idum) = pointer (k,i)
            end do
c
            numatm (i) = idum
            do k=1,idum
              pointer (k,i) = ndummy (k)
            end do
          end if
c
ccc          print *,i,' now has ',numatm(i),' atoms'
          ldone (j) = .true.
          lextend = .true.
          goto 1210
c
c === REVERT DIRECTION OF J AND TRY AGAIN
c
 4000     continue
          print *,'Invert ',j,' for use with ',i
c
          do k=1,numatm(j)
            ndummy(k) = pointer(k,j)
          end do
          do k=1,numatm(j)
            pointer (k,j) = ndummy (numatm(j)-k+1)
          end do
          j1 = pointer (1,j)
          j2 = pointer (numatm(j),j)
          nitocj (i,j) = dist (i1,j2,newxyz)
          nitocj (j,i) = dist (i2,j1,newxyz)
          angvec (i,j) = 180.0 - angvec(i,j)
c
          goto 1230
c
 1210     continue
          if (lextend) goto 1300
        end do
 1200   continue
      end do
c
c ... any changes ?
c
 1300 continue
      if (lextend) goto 10
c
c ... sort & write traces
c
      write (*,*)
      call prompt (' Could not extend any further')
      idum = 0
      do i=1,nsol
        if (.not.ldone(i)) idum = idum + 1
      end do
      call ivalut (' Nr of stretches left :',1,idum)
c
c ... open PDB file
c
      call xopxua (iunit,base,linter,k)
      call stamp (line)
      write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
c
c ... open DEJAVU SSE file
c
      junit = iunit + 11
      call xopxua (junit,dejavu,linter,k)
      call stamp (line)
      write (junit,'(a1,1x,a)') '!',line(1:leng1(line))
      write (junit,'(a1)') '!'
      write (junit,'(a)') 'MOL   bone'
      write (junit,'(a)') 'NOTE  auto-generated by SOLEX'
      write (junit,'(9a)') 'PDB   ',base(1:leng1(base))
      write (junit,'(a1)') '!'
      ssenam (1:1) = 'B'
      if (lalpha) ssenam (1:1) = 'A'
      ssetyp = 'BETA  '
      if (lalpha) ssetyp = 'ALPHA '
c
      idum = 0
      l = 0
c
 8000 continue
      i1 = 0
      i2 = 0
      do i=1,nsol
        if (.not. ldone(i)) then
          if (numatm(i) .gt. i1) then
            i1 = numatm(i)
            i2 = i
          end if
        end if
      end do
      if (i1 .le. 0) goto 9000
c
c ... strip off loose ends
c
      if (i1 .le. natoms) goto 9000
c
      idum = idum + 1
      print *,'Stretch ',idum,' has ',i1,' atoms (solution ',i2,' )'
      write (iunit,'(a6,1x,a,i6,a,i6)') 'REMARK','Stretch # ',idum,
     +  ' --- Nr of atoms ',i1
c
c ... add to PDB file
c
      do i=1,i1
        j = pointer(i,i2)
        k = newptr (j)
        line = atmlin (k)
        l = l + 1
c
        write (line(31:66),'(3f8.3,2f6.2)')
     +      (newxyz(m,j),m=1,3),float(idum),newb(j)
        write (line(23:26),'(i4)') l
        write (iunit,'(a)') line(1:leng1(line))
c
      end do
c
c ... add to DEJAVU SSE file
c
      write (ssenam(2:6),'(I5)') idum
      call remspa (ssenam)
      j = pointer(1,i2)
      k = pointer(i1,i2)
      write (line,6000) ssetyp,
     +  quote,ssenam(1:leng1(ssenam)),quote,
     +  quote,l-i1+1,quote,   quote,l,quote,   i1,
     +  (newxyz(m,j),m=1,3),(newxyz(m,k),m=1,3)
      call pretty (line(7:))
      write (junit,'(a)') line(1:leng1(line))
c
 6000 format (a6, 1x,a1,a,a1, 1x,a1,i6,a1, 1x,a1,i6,a1, i6,
     +        6f10.2)
c
c ... BETA  'B5'  '41' '47' 7 48.66 66.03 20.31 37.64 70.56 35.31
c ... ALPHA 'A1'  '61' '68' 8 56.92 69.65 44.14 49.58 62.42 45.02
c
      ldone (i2) = .true.
      goto 8000
c
 9000 continue
      write (iunit,'(a)') 'END'
      close (iunit)
      call prompt (' Written PDB file')
c
      write (junit,'(a)') 'ENDMOL'
      close (junit)
      call prompt (' Written DEJAVU SSE file')
c
      return
      end
