      program essens
c
c ... real-space convolution of structural templates with electron density
c
c ... Gerard Kleywegt & Alwyn Jones @ 950710
c
c modifications for parallel computing Kay Diederichs 11/99
c modification for faster sorting Kay Diederichs 11/99
c to use parallelization with SGI compilers, compile with f77 -mp -O3
c to use parallelization with PGI compilers on Linux, compile with f77 -mp -fast
c to use parallelization with Alpha compilers, compile _and_link_everything_ 
c (including the gksub routines) with f77 -omp -O3
c 02/01 KD get rid of mp_block. Change mp_numthreads to omp_get_max_threads
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'ESSENS', vers = '040701/3.0.5')
c
      include 'maxdim.incl'
c
      integer maxsiz, maxmsk
      parameter (maxsiz = maxgk1)
      parameter (maxmsk = maxgk1)
c
c      pointer (iaptr,mapa)
c      pointer (ibptr,mapb)
c      pointer (icptr,maskb)
c      pointer (idptr,maskx)
c
c      real mapa(1), mapb(1)
c      integer maskb(1), maskx(1), malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr,idptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr,idptr
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
      ibptr = fmalloc (nb)
c
      nb = wrdbyt*masksize
      icptr = fmalloc (nb)
      idptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0 .or.
     +    idptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0 .or.
     +    idptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call doesse (%val(iaptr),%val(ibptr),%val(icptr),%val(idptr),
     +             mapsize, masksize)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
      call ffree (idptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine doesse (mapa, mapb, maskb, maskx, maxsiz, maxmsk)
c
c ... real-space convolution of structural templates with electron density
c
c ... Gerard Kleywegt & Alwyn Jones @ 950710
c
      implicit none
c
      include 'maxdim.incl'
c
      integer    maxsiz,maxatm,maxmsk
c      parameter (maxsiz = maxgk1)
      parameter (maxatm = 10000)
c
      real mapa(maxsiz), mapb(maxsiz)
      real xyz(3,maxatm), rotxyz(3,maxatm), dens(maxatm)
      real rot(3,3),cell(6),grid(3),xt(3),cella(6),cellx(6)
      real xdum,q,xcut
c
      integer maskb(maxmsk),maskx(maxmsk)
      integer orgna(3),exta(3),grida(3),uvwa(3)
      integer orgnx(3),extx(3),gridx(3)
      integer delind(3,maxatm)
      integer spgrp,natoms,ierr,iunit,junit,kunit,length
      integer i,j,idum,nx,ny,nz,n,kmin,ntot,neva
c
      character atmnam(maxatm)*4,atmlin(maxatm)*80
      character mapf*80,mapm*80,rotf*80,file*80,line*80,t*1
      character answer*1
c
      logical xinter,leuler,lokay,linter,locmax
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
c ... get input map
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
      call edin   (file, iunit, mapa,  orgna, exta, grida, uvwa,
     $             cella, spgrp, maxsiz)
c
      close (iunit)
      call telmap (grida,orgna,exta,cella)
c
      do i=1,3
        grid (i) = cella (i) / float (grida(i))
      end do
      call fvalut (' Grid spacing (A):',3,grid)
c
      if (exta(1)*exta(2)*exta(3) .gt. maxmsk) then
        call errcon ('Internal mask to small')
        call jvalut (' Allocated :',1,maxmsk)
        call jvalut (' Required  :',1,exta(1)*exta(2)*exta(3))
        return
      end if
c
      do i=1,6
        cell (i) = cella (i)
      end do
c
c ... optionally: read mask
c
      write (*,*)
      write (*,*) 'Mask is *OPTIONAL* !'
      file = ' '
      call textin (' Input O/MAMA mask file ?',file)
      call textut (' Input O/MAMA mask file :',file)
c
      if (length(file) .lt. 1) then
c
c ... no mask
c
        call prompt (' No mask file provided - use all points in map')
c
        do i=1,maxmsk
          maskx (i) = 1
        end do
c
        do i=1,3
          orgnx(i) = orgna(i)
          gridx(i) = grida(i)
          extx(i) = exta(i)
          cellx(i) = cella(i)
          cellx(i+3) = cella(i+3)
        end do
c
      else
c
c ... read mask
c
        call xopxoa (iunit,file,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening mask')
          return
        end if
        call maskin 
     $    (iunit, maskx, orgnx, extx, gridx, cellx, maxmsk, ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading mask')
          return
        end if
c
        close (iunit)
        call telmap (gridx,orgnx,extx,cellx)
c
c ... check if map and mask have same cell and grid
c
        lokay = .true.
        do i=1,3
          lokay = (lokay .and. (grida(i) .eq. gridx(i)))
          lokay = (lokay .and. (abs(cella(i)-cellx(i)) .le. 0.01))
          lokay = (lokay .and.
     +             (abs(cella(i+3)-cellx(i+3)) .le. 0.01))
        end do
        if (.not. lokay) then
          call errcon (
     +      'Map and mask do NOT have the same cell and grid !')
          return
        end if
      end if
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
      call pdboni ('P',iunit,maxatm,natoms,xyz,atmnam,atmlin,ierr)
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
      write (*,*)
      call fvalut (' Centre of gravity of atoms :',3,xt)
c
      xdum = 9999.999
      idum = -1
      do i=1,natoms
        if (atmnam(i) .eq. ' CA ') then
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
      xdum = -1.0
      idum = 1
      do i=1,natoms
        q = 0.0
        do j=1,3
          xyz(j,i) = xyz(j,i) - xt(j)
          q = q + (xyz(j,i))**2
        end do
        q = sqrt (q)
        if (q .gt. xdum) then
          xdum = q
          idum = i
        end if
      end do
      write (*,*)
      call fvalut (' Furthest atom (A):',1,xdum)
      line = atmlin (idum)
      call pretty (line)
      call textut (' Atom :',line)
      nx = 2 + int(xdum/grid(1))
      ny = 2 + int(xdum/grid(2))
      nz = 2 + int(xdum/grid(3))
      write (*,*)
      write (*,'(1x,a,3i3)') 'Safety border NX/NY/NZ = ',nx,ny,nz
      call prompt (' Check if probe can be rotated in map ...')
      write (*,'(1x,a,i6,a,i6,a)') 'X - extent ',exta(1),
     +  ' >= 2*NX+1 = ',2*nx+1,' ?'
      write (*,'(1x,a,i6,a,i6,a)') 'Y - extent ',exta(2),
     +  ' >= 2*NY+1 = ',2*ny+1,' ?'
      write (*,'(1x,a,i6,a,i6,a)') 'Z - extent ',exta(3),
     +  ' >= 2*NZ+1 = ',2*nz+1,' ?'
c
      if (exta(1)-2*nx .le. 0) then
        call errcon ('Not enough points along X')
        return
      else if (exta(2)-2*ny .le. 0) then
        call errcon ('Not enough points along Y')
        return
      else if (exta(3)-2*nz .le. 0) then
        call errcon ('Not enough points along Z')
        return
      end if
      idum = (exta(1)-2*nx) * (exta(2)-2*ny) * (exta(3)-2*nz)
c
c ... get rotation ranges
c
      do i=1,3
        rot (1,i) =    0.0
        rot (2,i) =  359.0
        rot (3,i) =   10.0
      end do
c
      write (*,*)
      call prompt (' Rotations may be in Euler or Polar angles')
      t = 'Y'
      call textin (' Use Euler angles (Y/N) ?',t)
      call upcase (t)
      leuler = (t .ne. 'N')
c
      if (leuler) then
        rot (2,2) = 179.0
        call prompt (' Rotations are Euler angles in degrees')
        call fvalin (' Rotation ALPHA start, end, step ?',3,rot(1,1))
        call fvalin (' Rotation BETA  start, end, step ?',3,rot(1,2))
        call fvalin (' Rotation GAMMA start, end, step ?',3,rot(1,3))
        call fvalut (' Rotation ALPHA start, end, step :',3,rot(1,1))
        call fvalut (' Rotation BETA  start, end, step :',3,rot(1,2))
        call fvalut (' Rotation GAMMA start, end, step :',3,rot(1,3))
      else
        rot (2,3) = 179.0
        call prompt (' Rotations are Polar angles in degrees')
        call fvalin (' Rotation PHI   start, end, step ?',3,rot(1,1))
        call fvalin (' Rotation PSI   start, end, step ?',3,rot(1,2))
        call fvalin (' Rotation KAPPA start, end, step ?',3,rot(1,3))
        call fvalut (' Rotation PHI   start, end, step :',3,rot(1,1))
        call fvalut (' Rotation PSI   start, end, step :',3,rot(1,2))
        call fvalut (' Rotation KAPPA start, end, step :',3,rot(1,3))
      end if
c
      n = 1
      do i=1,3
        call rlohi (rot(1,i),rot(2,i))
        rot (3,i) = max (0.1, rot(3,i))
        n = n * nint ((rot(2,i)-rot(1,i))/rot(3,i) + 1.0)
      end do
      call jvalut (' Number of rotations :',1,n)
c
c ... get value of K
c
      kmin = 6 * natoms
      kmin = kmin / 10
      kmin = max(kmin,3)
      write (*,*)
      call prompt (' K = Nr of atoms to sum in K-minimum sum function')
      call ivalin (' Value of K ?',1,kmin)
      kmin = max (1, min (natoms,kmin))
      call ivalut (' K :',1,kmin)
c
c ... get density cut-off
c
      xcut = mapa(1)
      do i=2,exta(1)*exta(2)*exta(3)
        if (mapa(i).lt.xcut) xcut=mapa(i)
      end do
      write (*,*)
      call rvalut (' Minimum density in input map     :',1,xcut)
      call rvalin (' Lower density cut-off ?',1,xcut)
      call rvalut (' Lower density cut-off :',1,xcut)
c
      write (*,*)
      ntot = exta(1)*exta(2)*exta(3)
      neva = 0
      do i=1,exta(1)*exta(2)*exta(3)
        if (mapa(i).ge.xcut) neva = neva + 1
      end do
      call jvalut (' Grid points in map     :',1,ntot)
      call jvalut (' Grid points >= cut-off :',1,neva)
c
c ... only look at local density maxima ?
c
      write (*,*)
      call prompt (' The calculations can be sped up considerably')
      call prompt (' by only considering local electron density')
      call prompt (' maxima (although some weaker features may')
      call prompt (' be missed !)')
      answer = 'N'
      call textin (' Only consider local maxima (Y/N) ?',answer)
      call textut (' Only consider local maxima (Y/N) :',answer)
      call upcase (answer)
      locmax = (answer .eq. 'Y')
      if (locmax) then
        call prompt (' Only local maxima will be considered')
      else
        call prompt (' All selected map points will be considered')
      end if
c
      write (*,*)
      file = 'score.E'
      call textin (' Output score map file ?',file)
      call textut (' Output score map file :',file)
      mapf = file
c
      write (*,*)
      file = 'display.E'
      call textin (' Output display map file ?',file)
      call textut (' Output display map file :',file)
      mapm = file
c
      write (*,*)
      file = 'temp.rot'
      call textin (' Output rotation file ?',file)
      call textut (' Output rotation file :',file)
      rotf = file
c
c ... Note: MASKX and MASKY are the same masks, but used for
c     two different purposes; MASKY is only used once MASKX
c     is no longer needed, so this is safe; since they have
c     different sizes, they have to be used as separate arrays
c     in the subroutine, but there is now no need to allocate
c     yet another mask
c
      call sub03x (mapa, mapb, maskb, exta(1), exta(2), exta(3), 
     +  orgna, cell, grid, leuler, natoms, xyz, rot, nx, ny, nz,
     +  delind, rotxyz, xcut, kmin, dens, 
     +  maskx, extx(1), extx(2), extx(3), orgnx,
     +  iunit, grida, rotf, maskx, locmax, ierr)
c
      if (ierr .ne. 0) return
c
c ... write the convolution map ("score")
c
      write (*,*)
      call prompt (' Writing convolution (score) map ...')
      close (iunit)
      call edout (mapf,iunit,mapb,orgna,exta,grida,uvwa,cella,spgrp)
c
c ... write the display map ("Morten")
c
      write (*,*)
      call prompt (' Writing display (Morten) map ...')
      close (iunit)
      call edout (mapm,iunit,mapa,orgna,exta,grida,uvwa,cella,spgrp)
      close (iunit)
c
      return
c
      end
c
c
c
      subroutine sub03x 
     + (mapa, mapb, maskb, exta1, exta2, exta3, orgna, 
     +  cell, grid, leuler, natoms, xyz, rot, nx, ny, nz,
     +  delind, rxyz, xcut, kmin, dens,
     +  maskx, mx, my, mz, orgnx, iunit, grida, rotf, masky,
     +  locmax, ierr)
c
c --- do it
c
      implicit none
c
      integer orgna(3), exta1, exta2, exta3, nx, ny, nz, natoms
      integer mx, my, mz, orgnx(3)
c
      real mapa(exta1, exta2, exta3)
      real mapb(exta1, exta2, exta3)
      real xyz(3,natoms), rxyz(3,natoms)
      real a(3,3), b(3,3), rot(3,3)
      real cell(6), grid(3), rt(9), del(3), dens(natoms)
c
      real total,user,sys,q,xcut,oneovf,oneof2,f
      real sumy,sumysq,qinit,perc,ws,qmax,amax,bmax,cmax
      real fkmin,one8th,qq,q3,q4,q5,xdum
c
      integer maskb(exta1, exta2, exta3)
      integer masky(exta1, exta2, exta3)
      integer maskx(mx, my, mz)
      integer delind (3,natoms),ilim(2,3),grida(3),ext(3)
c
      integer n1,n2,n3,i,i1,i2,i3,j,l,ngk,k,imax,jmax,kmax,nn
      integer ndel,ntot,ni,nj,nk,mi,mj,mk,kmin,is,js,ks,ms,ipm
      integer npo,nne,lmin,lmax,io1,io2,io3,rotpck,ii,jj,kk
      integer msol,iunit,leng1,ierr,ipart

      integer ithread, nthread, omp_get_max_threads, lenthr, mk1, nk1
c
      logical leuler,first,xinter,locmax
c
      character line*256
      character rotf*(*)
c
code ...
c
      ierr = 0
c
c ... zero output map
c
      write (*,*)
      call prompt (' Initialising output map ...')
c
      qinit = -9999.99
c
      do k=1,exta3
        do j=1,exta2
          do i=1,exta1
            mapb  (i,j,k) = qinit
            maskb (i,j,k) = 0
          end do
        end do
      end do
c
c ... best fit
c
      qmax = -9.9E29
      imax=0
      jmax=0
      kmax=0
      amax=0.0
      bmax=0.0
      cmax=0.0
c
      f = float(natoms)
      oneovf = 1.0 / f
      oneof2 = oneovf*oneovf
c
      one8th = 1.0 / 8.0
      fkmin = 1.0 / float (kmin)
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
      n1 = nint ( (rot(2,1)-rot(1,1)) / rot(3,1) ) + 1
      n2 = nint ( (rot(2,2)-rot(1,2)) / rot(3,2) ) + 1
      n3 = nint ( (rot(2,3)-rot(1,3)) / rot(3,3) ) + 1
      ntot = n1*n2*n3
      ndel = max(1,nint(0.01*float(ntot)))
      ngk = 0
      first = .true.
c
c ... get safety limits (also if non-orthorhombic cell ;-)
c
      write (*,*)
      call prompt (' Determining safety borders ...')
c
      do i1=1,n1
        del(1) = rot (1,1) + (i1-1)*rot(3,1)
        do i2=1,n2
          del(2) = rot (1,2) + (i2-1)*rot(3,2)
          do i3=1,n3
            del(3) = rot (1,3) + (i3-1)*rot(3,3)
c
            if (leuler) then
              call ccpeul (del,rt)
            else
              call ccppol (del,rt)
            end if
c
            call getdel (natoms,xyz,rxyz,rt,b,cell,grid,delind)
c
            do i=1,natoms
              nx = max (nx, iabs(delind(1,i)))
              ny = max (ny, iabs(delind(2,i)))
              nz = max (nz, iabs(delind(3,i)))
            end do
c
          end do
        end do
      end do
c
      write (*,'(1x,a,3i3)') 'Safety border X/Y/Z = ',nx,ny,nz
      call gkdcpu (total,user,sys)
      write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      mk = nz+1
      nk = exta3-nz-1
      mj = ny+1
      nj = exta2-ny-1
      mi = nx+1
      ni = exta1-nx-1
c
c ... set up mask
c
      write (*,*)
      call prompt (' Setting up mask ...')
c
      lmin = max (orgna(1) + mi, orgnx(1))
      lmax = min (orgna(1) + ni, orgnx(1) + mx - 1)
      if (lmin .ge. lmax) then
        print *,' (X) ... ',orgna(1) + mi,orgnx(1),
     +          orgna(1) + ni,orgnx(1) + mx - 1,lmin,lmax     
        call errcon (
     +    'User mask does not overlap with central map region')
        ierr = -1
        return
      end if
      ilim (1,1) = lmin - orgna(1) + 1
      ilim (2,1) = lmax - orgna(1) + 1
      io1 = (lmin - orgnx(1) + 1) - ilim (1,1)
c
      lmin = max (orgna(2) + mj, orgnx(2))
      lmax = min (orgna(2) + nj, orgnx(2) + my - 1)
      if (lmin .ge. lmax) then
        print *,' (Y) ... ',orgna(2) + mj,orgnx(2),
     +         orgna(2) + nj,orgnx(2) + my - 1,lmin,lmax     
        call errcon (
     +    'User mask does not overlap with central map region')
        ierr = -1
        return
      end if
      ilim (1,2) = lmin - orgna(2) + 1
      ilim (2,2) = lmax - orgna(2) + 1
      io2 = (lmin - orgnx(2) + 1) - ilim (1,2)
c
      lmin = max (orgna(3) + mk, orgnx(3))
      lmax = min (orgna(3) + nk, orgnx(3) + mz - 1)
      if (lmin .ge. lmax) then
        print *,' (Z) ... ',orgna(3) + mk,orgnx(3),
     +         orgna(3) + nk,orgnx(3) + mz - 1,lmin,lmax     
        call errcon (
     +    'User mask does not overlap with central map region')
        ierr = -1
        return
      end if
      ilim (1,3) = lmin - orgna(3) + 1
      ilim (2,3) = lmax - orgna(3) + 1
      io3 = (lmin - orgnx(3) + 1) - ilim (1,3)
c
      nn = 0
      do i1=ilim(1,1),ilim(2,1)
        do i2=ilim(1,2),ilim(2,2)
          do i3=ilim(1,3),ilim(2,3)
            if (maskx(i1+io1,i2+io2,i3+io3) .eq. 1) then
              maskb(i1,i2,i3) = 1
              nn = nn + 1
            else
              maskb(i1,i2,i3) = 0
            end if
          end do
        end do
      end do
c
      call jvalut (' Grid points inside mask and map :',1,nn)
c
      if (nn .le. 0) then
        call errcon ('No points left to try - check mask and input')
        ierr = -1
        return
      end if
c
c ... now apply density cut-off (and ignore safety borders)
c
      write (*,*)
      call prompt (' Applying density cut-off ...')
      nn = 0
      do k=mk,nk
        do j=mj,nj
          do i=mi,ni
            if (maskb(i,j,k) .eq. 1) then
              if (mapa(i,j,k).ge.xcut) then
                maskb(i,j,k) = 2
                nn = nn + 1
              else
                maskb(i,j,k) = 0
              end if
            end if
          end do
        end do
      end do
c
c ... now all valid mask points should have a value of 2
c    reset the points that don't
c    initialise MASKX array for storing best rotations
c
      write (*,*)
      call prompt (' Marking grid points to use ...')
      nn = 0
c
      do k=1,exta3
        do j=1,exta2
          do i=1,exta1
c
c ... NOTE: at this stage MASKX (which is equivalenced to MASKY) has
c           outlived its usefullness, so we can recycle it here
c
            masky (i,j,k) = -1
c
            if (maskb(i,j,k) .eq. 2) then
              maskb (i,j,k) = 1
              nn = nn + 1
            else
              maskb (i,j,k) = 0
            end if
          end do
        end do
      end do
c
      call jvalut (' After density cut-off and borders :',1,nn)
c
      if (nn .le. 0) then
        call errcon ('No points left to try - check mask and input')
        ierr = -1
        return
      end if
c
c ... if requested, select only local maxima
c
      if (locmax) then
c
        write (*,*)
        call prompt (' Picking local maxima ...')
        nn = 0
c
        do k=2,exta3-1
          do j=2,exta2-1
            do i=2,exta1-1
              if (maskb(i,j,k) .ne. 1) goto 6206
              xdum = mapa (i,j,k)
c
		  if (xdum .lt. mapa(i-1,j-1,k-1)) goto 6210
		  if (xdum .lt. mapa(i-1,j-1,k  )) goto 6210
		  if (xdum .lt. mapa(i-1,j-1,k+1)) goto 6210
		  if (xdum .lt. mapa(i-1,j  ,k-1)) goto 6210
		  if (xdum .lt. mapa(i-1,j  ,k  )) goto 6210
		  if (xdum .lt. mapa(i-1,j  ,k+1)) goto 6210
		  if (xdum .lt. mapa(i-1,j+1,k-1)) goto 6210
		  if (xdum .lt. mapa(i-1,j+1,k  )) goto 6210
		  if (xdum .lt. mapa(i-1,j+1,k+1)) goto 6210
c
		  if (xdum .lt. mapa(i  ,j-1,k-1)) goto 6210
		  if (xdum .lt. mapa(i  ,j-1,k  )) goto 6210
		  if (xdum .lt. mapa(i  ,j-1,k+1)) goto 6210
		  if (xdum .lt. mapa(i  ,j  ,k-1)) goto 6210
		  if (xdum .lt. mapa(i  ,j  ,k+1)) goto 6210
		  if (xdum .lt. mapa(i  ,j+1,k-1)) goto 6210
		  if (xdum .lt. mapa(i  ,j+1,k  )) goto 6210
		  if (xdum .lt. mapa(i  ,j+1,k+1)) goto 6210
c
		  if (xdum .lt. mapa(i+1,j-1,k-1)) goto 6210
		  if (xdum .lt. mapa(i+1,j-1,k  )) goto 6210
		  if (xdum .lt. mapa(i+1,j-1,k+1)) goto 6210
		  if (xdum .lt. mapa(i+1,j  ,k-1)) goto 6210
		  if (xdum .lt. mapa(i+1,j  ,k  )) goto 6210
		  if (xdum .lt. mapa(i+1,j  ,k+1)) goto 6210
		  if (xdum .lt. mapa(i+1,j+1,k-1)) goto 6210
		  if (xdum .lt. mapa(i+1,j+1,k  )) goto 6210
		  if (xdum .lt. mapa(i+1,j+1,k+1)) goto 6210
c
              nn = nn + 1
              goto 6206
c
 6210         continue
              maskb(i,j,k) = 0
c
 6206         continue
            end do
          end do
        end do
c
      call jvalut (' After picking local maxima :',1,nn)
      end if
c
      if (nn .le. 0) then
        call errcon ('No points left to try - check mask and input')
        ierr = -1
        return
      end if
c
      call gkdcpu (total,user,sys)
      write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      write (*,*)
      call prompt (' Please be patient for a while ...')
      write (*,*)
      nthread=1
c an OpenMP Fortran compiler considers the next line as conditional compilation 
c (indicated by C$) if using the -mp flag (SGI or PGI/Linux) , or -omp
c (Alpha) on the f77 line
C$    nthread=omp_get_max_threads()
C$    if (nthread.lt.1) nthread=1
C$    print*,'nthread,mk,nk,nk-mk+1=',nthread,mk,nk,nk-mk+1 
      lenthr=(nk-mk+1)/nthread
      if (lenthr*nthread.lt.nk-mk+1) then
        lenthr=lenthr+1
      endif
c
c try to evenly distribute the job among processes. the last process gets 
c least to do.
C$OMP PARALLEL DO PRIVATE(ithread,mk1,nk1) SHARED(mi,mj,n1,n2,n3,ni,nj,
c$omp& ndel,ntot,exta1,exta2,exta3,leuler,natoms,kmin,nthread) 
      do ithread=1,nthread
        mk1=mk+(ithread-1)*lenthr
        if (mk1.le.nk) then
          nk1=mk1+lenthr-1
          if (nk1.gt.nk) nk1=nk
          call do_par(mi,mj,mk1,n1,n2,n3,ni,nj,nk1,b,ithread,
     +      ndel,ntot,
     +      mapa, mapb, maskb, exta1, exta2, exta3,
     +      cell, grid, leuler, natoms, xyz, rot, 
     +      kmin, masky )
        endif
      end do
      write (*,*)
      call prompt (' Done all rotations ... post-processing ...')
c
      call gkdcpu (total,user,sys)
      write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      do i=1,exta1
        do j=1,exta2
          do k=1,exta3
            mapb (i,j,k) = mapb (i,j,k) * one8th * fkmin
          end do
        end do
      end do
c
c ... set all unset points to be equal to the minimum in the map
c
      write (*,*)
      call prompt (' Getting statistics ...')
c
      q = 9.99E29
      qq = -q
      sumysq = 0.0
      sumy = 0.0
      nn = 0
      npo = 0
      nne = 0
c
      do k=mk,nk
        do j=mj,nj
          do i=mi,ni
c
            if (maskb(i,j,k).eq.1) then
              nn = nn + 1
              if (nn .eq. 1) then
                q = mapb(i,j,k)
                qq = mapb(i,j,k)
              end if
              if (mapb(i,j,k).lt.q) q = mapb(i,j,k)
              if (mapb(i,j,k).gt.qq) qq = mapb(i,j,k)
              sumy   = sumy + mapb(i,j,k)
              sumysq = sumysq + mapb(i,j,k)*mapb(i,j,k)
              if (mapb(i,j,k) .gt. 0.0) then
                npo = npo + 1
              else if (mapb(i,j,k) .lt. 0.0) then
                nne = nne + 1
              end if
            end if
c
          end do
        end do
      end do
c
      write (*,*)
      call jvalut (' Nr of points in masked map :',1,nn)
      call jvalut (' Nr of points with positive signal :',1,npo)
      call jvalut (' Nr of points with negative signal :',1,nne)
c
      write (*,*)
      call rvalut (' Minimum signal in masked map :',1,q)
      call rvalut (' Maximum signal in masked map :',1,qq)
c
      oneovf = 1.0 / float(nn)
      oneof2 = oneovf * oneovf
      q3 = sumy * oneovf
      write (*,*)
      call rvalut (' Average signal in masked map :',1,q3)
      q4 = sumysq*oneovf - sumy*sumy*oneof2
      q5 = 0.0
      if (q4 .gt. 0.0) then
        q4 = sqrt(q4)
        write (*,*)
        call rvalut (' St.dev. signal in masked map :',1,q4)
        q5 = (qq-q3) / q4
        call rvalut (' ( Max - Ave ) / St.dev.      :',1,q5)
      else
        write (*,*)
        call errcon ('While computing standard deviation')
      end if
c
      write (*,*)
      call prompt (' Resetting all unset points to minimum ...')
c
      do k=1,exta3
        do j=1,exta2
          do i=1,exta1
            if (maskb(i,j,k).eq.0) mapb(i,j,k) = q
          end do
        end do
      end do
c
      call gkdcpu (total,user,sys)
      write (*,'(/a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
c ... write file with best rotations for every point
c
      write (*,*)
      call prompt (' Writing rotation map ...')
c
      call xopxua (iunit,rotf,xinter(),ierr)
      if (ierr .ne. 0) goto 6969
c
      call stamp (line)
      write (iunit,'(a,1x,a)',err=6969) 'ROTATION_MAP_ESSENS',
     +  line(1:leng1(line))
c
      if (leuler) then
        write (iunit,'(a)',err=6969) 'EULER angles'
      else
        write (iunit,'(a)',err=6969) 'POLAR angles'
      end if
c
      write (iunit,'(1p,2e14.6,1x,a)',err=6969) q3,q4,'AVE SDV'
      write (iunit,'(3i5,1x,a)',err=6969) (orgna(i),i=1,3),'Origin'
      write (iunit,'(3i5,1x,a)',err=6969) exta1,exta2,exta3,'Extent'
      write (iunit,'(3i5,1x,a)',err=6969) (grida(i),i=1,3),'Grid'
      write (iunit,'(6f10.4,1x,a)',err=6969) (cell(i),i=1,6),'Cell'
      call dumpit (iunit,masky,exta1*exta2*exta3)
c
      call gkdcpu (total,user,sys)
      write (*,'(/a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
 6969 continue
      close (iunit)
c
      write (*,*)
      call prompt (' Applying Morten''s trick ...')
c
c ... copy score map back to start map
c
      do k=1,exta3
        do j=1,exta2
          do i=1,exta1
            mapa (i,j,k) = mapb (i,j,k)
          end do
        end do
      end do
c
      do k=mk,nk
        do j=mj,nj
          do i=mi,ni
c
            if (maskb(i,j,k) .eq. 0) goto 2342
c
c ... get the *original* score for this point
c
            q = mapb (i,j,k)
c
c ... get the rotation angles that gave the best score
c
            msol = masky (i,j,k)
            call packut (i1,i2,i3,ipart,msol)
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
            call getdel (natoms,xyz,rxyz,rt,b,cell,grid,delind)
c
c ... make a "print" of the current orientation in the display map
c
            do l=1,natoms
c
              ii = i+delind(1,l)
              jj = j+delind(2,l)
              kk = k+delind(3,l)
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)+1
              jj = j+delind(2,l)
              kk = k+delind(3,l)
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)+1
              jj = j+delind(2,l)+1
              kk = k+delind(3,l)
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)+1
              jj = j+delind(2,l)+1
              kk = k+delind(3,l)+1
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)
              jj = j+delind(2,l)+1
              kk = k+delind(3,l)
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)
              jj = j+delind(2,l)+1
              kk = k+delind(3,l)+1
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)
              jj = j+delind(2,l)
              kk = k+delind(3,l)+1
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
              ii = i+delind(1,l)+1
              jj = j+delind(2,l)
              kk = k+delind(3,l)+1
              if (q .gt. mapa (ii,jj,kk)) then
                if (maskb(ii,jj,kk) .eq. 0) maskb(ii,jj,kk) = 2
                mapa (ii,jj,kk) = q
              end if
c
            end do
c
 2342       continue
c
          end do
        end do
      end do
c
      call gkdcpu (total,user,sys)
      write (*,'(/a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      write (*,*)
      call prompt (' Getting statistics ...')
c
      q = 9.99E29
      qq = -q
      sumysq = 0.0
      sumy = 0.0
      nn = 0
      npo = 0
      nne = 0
c
      do k=mk,nk
        do j=mj,nj
          do i=mi,ni
c
            if (maskb(i,j,k).gt.0) then
c
              maskb(i,j,k) = 1
c
              nn = nn + 1
              if (nn .eq. 1) then
                q = mapa(i,j,k)
                qq = mapa(i,j,k)
              end if
              if (mapa(i,j,k).lt.q) q = mapa(i,j,k)
              if (mapa(i,j,k).gt.qq) qq = mapa(i,j,k)
              sumy   = sumy + mapa(i,j,k)
              sumysq = sumysq + mapa(i,j,k)*mapa(i,j,k)
              if (mapa(i,j,k) .gt. 0.0) then
                npo = npo + 1
              else if (mapa(i,j,k) .lt. 0.0) then
                nne = nne + 1
              end if
            end if
c
          end do
        end do
      end do
c
      write (*,*)
      call jvalut (' Nr of points in display map :',1,nn)
      call jvalut (' Nr of points with positive signal :',1,npo)
      call jvalut (' Nr of points with negative signal :',1,nne)
c
      write (*,*)
      call rvalut (' Minimum signal in display map :',1,q)
      call rvalut (' Maximum signal in display map :',1,qq)
c
      oneovf = 1.0 / float(nn)
      oneof2 = oneovf * oneovf
      q3 = sumy * oneovf
      write (*,*)
      call rvalut (' Average signal in display map :',1,q3)
      q4 = sumysq*oneovf - sumy*sumy*oneof2
      q5 = 0.0
      if (q4 .gt. 0.0) then
        q4 = sqrt(q4)
        write (*,*)
        call rvalut (' St.dev. signal in display map :',1,q4)
        q5 = (qq-q3) / q4
        call rvalut (' ( Max - Ave ) / St.dev.       :',1,q5)
      else
        write (*,*)
        call errcon ('While computing standard deviation')
      end if
c
      write (*,*)
      call prompt (' Resetting all unset points to minimum ...')
c
      do k=1,exta3
        do j=1,exta2
          do i=1,exta1
            if (maskb(i,j,k).eq.0) mapa(i,j,k) = q
          end do
        end do
      end do
c
      call gkdcpu (total,user,sys)
      write (*,'(/a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      return
      end
c
c
c
      subroutine getdel (natoms,xyz,rxyz,rt,b,cell,grid,delind)
c
      implicit none
c
      integer natoms,i,l
      integer delind(3,natoms)
c
      real xyz(3,natoms),rxyz(3,natoms),b(3,3),rt(9),cell(6)
      real grid(3),x1(3),x(3)
c
code ...
c
      do i=1,natoms
c
c ... apply rotation
c
        call mulmtx (rt,xyz(1,i),rxyz(1,i),3,3,1)
c
c ... fractionalise
c
        call mulmtx (b,rxyz(1,i),x1,3,3,1)
c
c ... get offset in grid points
c
        do l=1,3
          x(l) = x1(l)*cell(l)/grid(l)
          delind (l,i) = int (x(l))
        end do
c
      end do
c
      return
      end
c
c
c
      subroutine pdboni (typ,iunit,maxatm,natoms,xyz,atmnam,
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
      character atmnam(maxatm)*4,atmlin(maxatm)*80
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
        atmnam(natoms) = line(13:16)
c
c ... strip hydrogen atoms
c
        if (lhydro(atmnam(natoms))) then
          nh = nh + 1
          natoms = natoms - 1
          goto 200
        end if
c
c ... okay, read coordinates
c
        read (line, '(30x,3f8.3)')  (xyz(i,natoms),i=1,3)
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
      subroutine dumpit (iunit,maskx,size)
c
      implicit none
c
      integer maxlen
      parameter (maxlen = 256)
c
      integer iunit,size,i,j,k,nsum,leng1
      integer maskx(*)
c
      character line*(maxlen)
c
c ...
c
      nsum = 0
      maskx (size+1) = 0
c
      do i=1,size,12
        j = i + 11
        if (j .gt. size) j = size
        write (line,'(12i20)',err=9000) (maskx(k),k=i,j)
        call pretty (line)
        write (iunit,'(A)',err=9000) line(1:leng1(line))
        do k=i,j,2
          nsum = nsum + maskx(k) - maskx(k+1)
        end do
      end do
c
      write (iunit,'(a,1x,i20)',err=9000) 'CHECKSUM',nsum
c
 9000 continue
      close (iunit)
c
      return
      end
c
      subroutine do_par(mi,mj,mk,n1,n2,n3,ni,nj,nk,b,ithread,
     +          ndel,ntot,
     +          mapa, mapb, maskb, exta1, exta2, exta3, 
     +          cell, grid, leuler, natoms, xyz, rot, 
     +          kmin, masky )
c
c --- out-sourced only to allow for parallelization 
c KD 26.11.99
c
      implicit none
c
      integer maxatm
      parameter (maxatm = 10000)
      real rxyz(3,maxatm), dens(maxatm)
      integer delind(3,maxatm)

      integer exta1, exta2, exta3, natoms
c
      real mapa(exta1, exta2, exta3)
      real mapb(exta1, exta2, exta3)
      real xyz(3,natoms) 
      real b(3,3), rot(3,3)
      real cell(6), grid(3), rt(9), del(3)
c
      real q, perc, total, user, sys
c
      integer maskb(exta1, exta2, exta3)
      integer masky(exta1, exta2, exta3)
c
      integer n1,n2,n3,i,i1,i2,i3,j,l,k
      integer ni,nj,nk,mi,mj,mk,kmin
      integer rotpck

      integer ithread,ngk,ndel,ntot
c
      logical leuler,first
      
      if (natoms.gt.maxatm) stop 'natom.gt.maxatm in do_par'
      first=.true.

      ngk=0

      do i1=1,n1
        del(1) = rot (1,1) + (i1-1)*rot(3,1)
        do i2=1,n2
          del(2) = rot (1,2) + (i2-1)*rot(3,2)
          do i3=1,n3
            del(3) = rot (1,3) + (i3-1)*rot(3,3)
c
            if (ithread.eq.1) then
            ngk = ngk + 1
            if (ngk .eq. ngk/ndel*ndel) then
              perc = 100.0*float(ngk)/float(ntot)
              write (*,'(/a,i8,a,f6.1/a,3f8.2)')
     +          ' Rotation # ',ngk,' --- % = ',perc,' Angles  : ',
     +          del(1),del(2),del(3)
              call gkdcpu (total,user,sys)
              write (*,'(a,3f10.1)') ' CPU total/user/sys :',
     +          total,user,sys
              if (perc .gt. 0.9) then
                perc = (100.0 - perc) * total
                write (*,'(a,3f10.1)')
     +            ' Expected CPU time left (s, m, h) : ',perc,perc/60.0,
     +            perc/3600.0
              end if
c
c ... flush output
c
              call flusho (6)
            end if
            end if
c
c ... encode current rotations (integer) into a single integer
c

            call packin (nint(del(1)),nint(del(2)),nint(del(3)),
     +                   0,rotpck)
c
c ... generate rotation matrix
c
            if (leuler) then
              call ccpeul (del(1),rt)
            else
              call ccppol (del(1),rt)
            end if
c
c ... rotate molecule & get offset of each atom's 3 indices
c
            call getdel (natoms,xyz,rxyz,rt,b,cell,grid,delind)
c
            if (ithread.eq.1.and.first) then
              do i=1,min(10,natoms)
                write (*,*)
                call ivalut (' Atom :',1,i)
                call fvalut (' XYZ  :',3,xyz(1,i))
                call fvalut (' ROT  :',3,rxyz(1,i))
                call ivalut (' delI :',3,delind(1,i))
              end do
            end if
c
c ... loop over masked map points (K-minimum; SUM function)
c
            do k=mk,nk
              do j=mj,nj
                do i=mi,ni
c
                  if (maskb(i,j,k).eq.0) goto 1235
c
c ... the innermost of the seven nested loops ...
c
                  do l=1,natoms
                    dens (l) =  (
     +       mapa (i+delind(1,l),j+delind(2,l),k+delind(3,l)) +
     +       mapa (i+delind(1,l)+1,j+delind(2,l),k+delind(3,l)) +
     +       mapa (i+delind(1,l)+1,j+delind(2,l)+1,k+delind(3,l)) +
     +       mapa (i+delind(1,l)+1,j+delind(2,l)+1,k+delind(3,l)+1) +
     +       mapa (i+delind(1,l),j+delind(2,l)+1,k+delind(3,l)) +
     +       mapa (i+delind(1,l),j+delind(2,l)+1,k+delind(3,l)+1) +
     +       mapa (i+delind(1,l),j+delind(2,l),k+delind(3,l)+1) +
     +       mapa (i+delind(1,l)+1,j+delind(2,l),k+delind(3,l)+1)
     +                              )
c
                  end do
c
c ... sort density values (in-lined SHELL-SORT algorithm)
c
c                  is=1
c  10              continue
c                  if (is.gt.natoms) goto 20
c                  ms=2*is-1
c                  is=2*is
c                  goto 10
c  20              continue
c                  ms=ms/2
c                  if (ms.eq.0) goto 60
c                  ks=natoms-ms
c                  do js=1,ks
c                    is=js
c  50                continue
c                    ipm=is+ms
c                    if (dens(ipm).ge.dens(is)) goto 40
cc
cc ... swap densities
cc
c                    ws=dens(is)
c                    dens(is)=dens(ipm)
c                    dens(ipm)=ws
cc
c                    is=is-ms
c                    if (is.ge.1) goto 50
c  40                continue
c                  end do
c                  goto 20
cc
c  60              continue
c
        if (kmin.ne.natoms) call select(kmin,natoms,dens)
c ... calculate score by summing K lowest densities
c
                  q = 0.0
                  do l=1,kmin
                    q = q + dens(l)
                  end do
c
c ... better than current maximum ???
c
                  if (q. gt. mapb (i,j,k)) then
                    mapb (i,j,k) = q
                    masky (i,j,k) = rotpck
                  else if (first) then
                    mapb (i,j,k) = q
                    masky (i,j,k) = rotpck
                  end if
 1235             continue
c
                end do
              end do
            end do
c
            first = .false.
c
c ... end of loop over rotations
c
          end do
        end do
      end do
      return
      end
c
      subroutine select(k,n,arr)
c select the lowest k out of n elements of array arr
c this algorithm takes order(n) CPU time
c see Numerical Recipes chapter 8.5

      implicit none
      INTEGER k,n
      REAL arr(n)
      INTEGER i,ir,j,l,mid
      REAL a,temp
      l=1
      ir=n
1     if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
        endif
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      goto 1
      END
