      program ncs6d
c
c --- Try to find NCS operators through a 6-dimensional search
c
c ... Gerard Kleywegt @ 930331
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'NCS6D', vers = '080625/3.1.2')
c
      include 'maxdim.incl'
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
      call jvalut (' Allocate map of size    :',1,mapsize)
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
      call doncs6 (%val(iaptr), mapsize)
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
      subroutine doncs6 (mapa, maxsiz)
c
c --- Try to find NCS operators through a 6-dimensional search
c
c ... Gerard Kleywegt @ 930331
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz,maxatm,maxsym
c      parameter (maxsiz = maxgk1)
      parameter (maxsym = maxgk2)
      parameter (maxatm = 50000)
c
      integer orgna(3), exta(3), grida(3), uvwa(3)
      real mapa(maxsiz), cella(6), xyz(3,maxatm), dens(maxatm)
      real rotxyz(3,maxatm)
c
      real cell(6), grid(3), rtbtoa(12), rtsym(12,maxsym)
      real x0(3), xt(3)
      real rot(3,3),tra(3,3)
c
      character file*80, fmt*80, par*25, partyp*1
      character t*1
c
      integer ctsym,errcod,i,j,k,spgrp,natoms,ierr,n,length,nkeep
c
      logical xinter,leuler
c
code ...
c
c      call gainit (prognm,vers)
c
c      call jvalut (' Max size of map         :',1,maxsiz)
      call jvalut (' Max nr of (BONES) atoms :',1,maxatm)
      call jvalut (' Max nr of symmetry ops  :',1,maxsym)
c
      write (*,*)
      file = ' '
      call textin (' Map file ?',file)
      call textut (' Map file :',file)
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
c --- read in map
c
      call edin   (file, 1, mapa,  orgna, exta, grida, uvwa,
     $             cella, spgrp, maxsiz)
c
      close (1)
      call telmap (grida,orgna,exta,cella)
c
      do 100 i=1,3
100     grid(i) = cella(i)/float(grida(i))
      do 110 i=1,6
110     cell(i) = cella(i)
c
      t = 'B'
      write (*,*)
      call textin (' Read BONES or PDB file (B/P) ?',t)
      call textut (' Read BONES or PDB file (B/P) :',t)
      call upcase (t)
      if (t .ne. 'P' .and. t .ne. 'B') then
        call errcon ('Invalid selection')
        return
      end if
c
      file = ' '
      if (t .eq. 'P') then
        call textin (' Name of PDB file ?',file)
        call textut (' Name of PDB file :',file)
      else
        call textin (' Name of BONES file ?',file)
        call textut (' Name of BONES file :',file)
      end if
c
      call xopxoa (2,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('Unable to open file')
        return
      end if
c
      call pdboni (t,2,maxatm,natoms,xyz,ierr)
      close (2)
      if (ierr .ne. 0) then
        call errcon ('While reading file')
        return
      end if
c
      nkeep = 10
      write (*,*)
      write (*,*)
     +  'To speed up the calculations, you can opt to use only'
      write (*,*)
     +  'every Nth (bones) atom. This will give a speed-up factor'
      write (*,*)
     +  'equal to N. Suggested range for N is 3-20.'
      call jvalin (' Value for N ?',1,nkeep)
      nkeep = max (1, nkeep)
      call jvalut (' Value for N :',1,nkeep)
c
      if (nkeep .ne. 1) then
        call prompt (' Removing all atoms but every Nth ...')
        j = 0
        k = 0
        do i = 1,natoms
          k = k + 1
          if (k .eq. nkeep) then
            k = 0
            j = j + 1
            xyz (1,j) = xyz (1,i)
            xyz (2,j) = xyz (2,i)
            xyz (3,j) = xyz (3,i)
          end if
        end do
        natoms = j
        call jvalut (' Number of atoms now:',1,natoms)
      end if
c
      if (natoms .le. 3) then
        call errcon (' Too few atoms !')
        return
      end if
c
      file = ' '
      write (*,*)
      call textin (' File with symmetry operators ?',file)
      call textut (' File with symmetry operators :',file)
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
        call errcon ('Too many space group operators.')
        return
      end if
      if ((12*ctsym) .ne. j) then
        call errcon ('Invalid nr of elements; must be N*12')
        return
      end if
      do j=1,ctsym
        read (1, fmt) (rtsym(i,j),i=1,12)
        call fratra (rtsym(10,j))
      end do
      close (unit=1)
      call anasgs (ctsym,rtsym,.true.,ierr)
      if (ierr .ne. 0) then
        call errcon ('In spacegroup symmetry operators')
        return
      end if
c
      do i=1,3
        x0 (i) = 0.0
        xt (i) = 0.0
      end do
      write (*,*)
      call fvalin (' Approximate centre of NCS mate ?',3,x0)
      call fvalut (' Approximate centre of NCS mate :',3,x0)
c
      do i=1,natoms
        do j=1,3
          xt(j) = xt(j) + xyz(j,i)
        end do
      end do
      do i=1,3
        xt(i)=xt(i)/float(natoms)
      end do
      call fvalut (' Centre of gravity of atoms     :',3,xt)
      call fvalut (' Approximate centre of NCS mate :',3,x0)
      do i=1,12
        rtbtoa(i) = 0.0
      end do
      rtbtoa(1) = 1.0
      rtbtoa(5) = 1.0
      rtbtoa(9) = 1.0
      rtbtoa(10) = x0(1) - xt(1)
      rtbtoa(11) = x0(2) - xt(2)
      rtbtoa(12) = x0(3) - xt(3)
      write (*,'(a,4(/1x,3f15.6))')
     +  ' Initial operator :',rtbtoa
c
      do i=1,3
        rot (1,i) =    0.0
        rot (2,i) =  359.0
        rot (3,i) =   10.0
        tra (1,i) =  -10.0
        tra (2,i) =   10.0
        tra (3,i) =    2.0
      end do
c
      write (*,*)
      call prompt (' Rotations may be Euler or Polar')
      t = 'Y'
      call textin (' Use Euler angles (Y/N) ?',t)
      call upcase (t)
      call textut (' Use Euler angles (Y/N) :',t)
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
      write (*,*)
      call prompt (' Translations are in A relative to NCS-centre')
      call fvalin (' Translation X start, end, step  ?',3,tra(1,1))
      call fvalin (' Translation Y start, end, step  ?',3,tra(1,2))
      call fvalin (' Translation Z start, end, step  ?',3,tra(1,3))
      call fvalut (' Translation X start, end, step  :',3,tra(1,1))
      call fvalut (' Translation Y start, end, step  :',3,tra(1,2))
      call fvalut (' Translation Z start, end, step  :',3,tra(1,3))
c
      n = 1
      do i=1,3
        call rlohi (tra(1,i),tra(2,i))
        tra (3,i) = max (0.1, tra(3,i))
        call rlohi (rot(1,i),rot(2,i))
        rot (3,i) = max (0.1, rot(3,i))
        n = n * nint ((tra(2,i)-tra(1,i))/tra(3,i) + 1.0)
        n = n * nint ((rot(2,i)-rot(1,i))/rot(3,i) + 1.0)
      end do
      call jvalut (' Estimated number of trials :',1,n)
c
      write (*,*)
      call prompt (' You have a choice of interpolation methods:')
      call prompt (' N = density at nearest grid point')
      call prompt (' A = average density at 8 nearest grid points')
      call prompt (' L = linear interpolation using 8 points')
      call prompt (' F = full, 64-point spline interpolation')
      t = 'L'
      call textin (' Interpolation type (N/A/L/F) ?',t)
      call upcase (t)
      call textut (' Interpolation type (N/A/L/F) :',t)
      if (index ('NALF',t) .le. 0) t='L'
c
      file = 'rt_best.o'
      write (*,*)
      call textin (' File for best RT-operator ?',file)
      call textut (' File for best RT-operator :',file)
      call xopxua (2,file,xinter(),errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening RT file')
        return
      end if
c
      write (*,*)
      call prompt (' Please be patient for a while ...')
      write (*,*)
c
c ---	Now go do everything. This is to make subscripting easier
c
      call sub03x (mapa, exta(1), exta(2), exta(3), orgna, 
     $  cell, grid, rtsym, ctsym,leuler,
     $  natoms,xyz,rotxyz,rtbtoa,rot,tra,xt,x0,dens,t)
c
      return
c
      end
c
c
c
      subroutine sub03x 
     $ (mapa, exta1, exta2, exta3, orgna, 
     $ cell, grid, rtsym, ctsym, leuler,
     $ natoms,xyz,rxyz,rtbtoa,rot,tra,xmo,xncs,dens,type)
c
c ---	Find the operator by maximising the correlation coefficient
c
      implicit none
c
      integer mxkeep
      parameter (mxkeep=100)
c
      real scores(mxkeep),minsco,rtkeep(12,mxkeep),xdum
c
c ---	Map A data structures
      integer orgna(3), exta1, exta2, exta3
      real mapa(exta1, exta2, exta3)
c
      real cell(6), grid(3), rtbtoa(12), rot(3,3), tra(3,3)
      integer ctsym,natoms,j1,j2,j3,n1,n2,n3,m1,m2,m3,ext(3)
      real xyz(3,natoms), rtsym(12,ctsym), tdel(3)
      real rxyz(3,natoms), xncs(3)
c
      integer ctmask,errcod,i,i1,i2,i3,j,l,m,ngk,igk,lclose,leng1
      integer nerr1,nerr2
c
      real a(3,3), b(3,3),forgn(3),fext(3),rtatob(12)
      real value, x(3), x1(3), x2(3),mapbit(4,4,4)
      real avx,avxsq,avxy,avy,avysq,ccoef,f,peak,peakrt(12),rt(12)
      real del(3),xmo(3),xp(3),xo(3)
      real gext(3),total,user,sys,dens(natoms),xgk(3)
      real gave,gsdv
c
      logical bit,leuler
c
      character type*1,line*80
c
code ...
c
      minsco = 0.0
      gave = 0.0
      gsdv = 0.0
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
c ... initialise best operators
c
      do i=1,mxkeep
        scores (i) = minsco
        do j=1,9,3
          rtkeep (j,i) = 1.0
          rtkeep (j+1,i) = 0.0
          rtkeep (j+2,i) = 0.0
        end do
        do j=10,12
          rtkeep (j,i) = 0.0
        end do
      end do
c
c ---	Generate envelope edges in fractional cell coords
c
      do 130 i=1,3
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
130     gext(i)  = float(ext(i)+orgna(i)-2)*grid(i)/cell(i)
c
      call fvalut (' FORGN :',3,forgn)
      call fvalut (' FEXT  :',3,fext)
      call fvalut (' GEXT  :',3,gext)
c
      do 510 j=1,9
510     rtatob(j) = rtbtoa(j)
      call matinv (rtatob, 3, x, x1, x2)
      call mulmtx (rtatob, rtbtoa(10), rtatob(10), 3, 3, 1)
      do 520 j=10,12
520     xo(j-9) = -rtatob(j)
c
      call fvalut (' Rotation origin :',3,xo)
      call fvalut (' Centre of atoms :',3,xmo)
c
      call mulmtx (rtbtoa, xmo, xp, 3, 3, 1)
      do 610 i=1,3
610     xp(i) = xp(i)+ rtbtoa(i+9)
c
      peak = -999.
      ngk = 0
      nerr1 = 0
      nerr2 = 0
c
c --- Loop over atoms
c
        do 69 i=1,natoms
c
            dens (i) = -1.0E5
            call mulmtx (b, xyz(1,i), x1, 3, 3, 1)
c
c --- Is it in the envelope, if not make it (if possible)
c
            call frcsym (x1, forgn, fext, rtsym, ctsym, errcod)
            if (errcod .ne. 0) then
	      call bldbit (mapa, exta1, exta2, exta3, orgna, grid,
     $	        cell, forgn, fext, rtsym, ctsym, mapbit, x1, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call errcon ('Severe FRCSYM error 1')
                  call fvalut (' Coordinates :',3,xyz(1,i))
                  call fvalut (' Fractional  :',3,x1)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 69
              end if
	      do 225 l=1,3
	        x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
	        m = x(l)
	        x(l) = x(l)- float(m-1)
225           continue
	      call ave8 (mapbit, 4, 4, 4, x, value, errcod)
            else
              do 210 l=1,3
210             x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c --- x is now set of indices to the mapa array
              if (type.eq.'N') then
                call intrp1 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              else if (type.eq.'A') then
                call ave8 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              else if (type.eq.'L') then
                call linint 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              else
                call intrpl 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              end if
            end if
            if (errcod .ne. 0) then
              if (nerr2 .lt. 10) then
                call errcon ('Interpolation error 1')
              else if (nerr2 .eq. 10) then
                  call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
              goto 69
            end if
            dens (i) = value
69        continue
c
      call flusho (6)
c
      n1 = nint ( (rot(2,1)-rot(1,1)) / rot(3,1) ) + 1
      n2 = nint ( (rot(2,2)-rot(1,2)) / rot(3,2) ) + 1
      n3 = nint ( (rot(2,3)-rot(1,3)) / rot(3,3) ) + 1
      m1 = nint ( (tra(2,1)-tra(1,1)) / tra(3,1) ) + 1
      m2 = nint ( (tra(2,2)-tra(1,2)) / tra(3,2) ) + 1
      m3 = nint ( (tra(2,3)-tra(1,3)) / tra(3,3) ) + 1
c
      do 300 i1=1,n1
        del(1) = rot (1,1) + (i1-1)*rot(3,1)
      do 300 i2=1,n2
        del(2) = rot (1,2) + (i2-1)*rot(3,2)
      do 300 i3=1,n3
        del(3) = rot (1,3) + (i3-1)*rot(3,3)
c
c ... generate rotation matrix
c
        if (leuler) then
          call ccpeul (del,rt)
        else
          call ccppol (del,rt)
        end if
c
        write (*,*)
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call fvalut (' Rotations :',3,del)
        write (*,'(a,3(/1x,3f15.6))')
     +    ' Rotation matrix :',(rt(i),i=1,9)
c
ccc        call matana (rt)
c
ccc 950216        call anancs (1,rt,.true.,errcod)
c
c ... rotate molecule
c
        do i=1,natoms
          call mulmtx (rt,xyz(1,i),rxyz(1,i),3,3,1)
        end do
c
        call mulmtx (rt,xmo,x,3,3,1)
        call fvalut (' C-o-G after rotation :',3,x)
        rtbtoa(10) = xncs(1) - x(1)
        rtbtoa(11) = xncs(2) - x(2)
        rtbtoa(12) = xncs(3) - x(3)
        call fvalut (' Initial translation  :',3,rtbtoa(10))
c
        call fvalut (' Min score for the top :',1,minsco)
c
        do 300 j1=1,m1
          tdel(1) = tra (1,1) + (j1-1)*tra(3,1)
          rt(10)=rtbtoa(10)+tdel(1)
        do 300 j2=1,m2
          tdel(2) = tra (1,2) + (j2-1)*tra(3,2)
          rt(11)=rtbtoa(11)+tdel(2)
        do 300 j3=1,m3
c
          ngk = ngk + 1
c
          tdel(3) = tra (1,3) + (j3-1)*tra(3,3)
          rt(12)=rtbtoa(12)+tdel(3)
c
c ---	Set sums to zero
c
        avx = 0.
        avy = 0.
        avxy = 0.
        avxsq = 0.
        avysq = 0.
        ctmask = 0
c
c ---	Loop over atoms
c
        do 100 i=1,natoms
            if (dens(i) .lt. -0.99E5) goto 100
c
            x(1) = rxyz(1,i) + rt(10)
            x(2) = rxyz(2,i) + rt(11)
            x(3) = rxyz(3,i) + rt(12)
c
            call mulmtx (b, x, x1, 3, 3, 1)
c ---	    Is it in the envelope, if not make it (if possible)
            call frcsym (x1, forgn, gext, rtsym, ctsym, errcod)
            if (errcod .ne. 0) then
	      call bldbit (mapa, exta1, exta2, exta3, orgna, grid,
     $	        cell, forgn, fext, rtsym, ctsym, mapbit, x1, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call prompt (' FRCSYM error 2')
                  call jvalut (' Atom nr :',1,i)
                  call fvalut (' Coordinates :',3,x)
                  call fvalut (' Fractional  :',3,x1)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 100
              end if
	      do 125 l=1,3
	        x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
	        m = x(l)
	        x(l) = x(l)- float(m-1)
125           continue
c
              if (type.eq.'N') then
                call intrp1 
     $            (mapbit, 4, 4, 4, x, value, errcod)
              else if (type.eq.'A') then
                call ave8 
     $            (mapbit, 4, 4, 4, x, value, errcod)
              else if (type.eq.'L') then
                call linint 
     $            (mapbit, 4, 4, 4, x, value, errcod)
              else
                call intrpl 
     $            (mapbit, 4, 4, 4, x, value, errcod)
              end if
c
              bit = .true.
            else
              do 110 l=1,3
110             x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c ---	      x is now set of indices to the mapa array
              if (type.eq.'N') then
                call intrp1 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              else if (type.eq.'A') then
                call ave8 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              else if (type.eq.'L') then
                call linint 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              else
                call intrpl 
     $            (mapa, exta1, exta2, exta3, x, value, errcod)
              end if
c
              bit = .false.
            end if
            if (errcod .eq. 0) then
              avx = avx+ dens(i)
              avy = avy+ value
              avxsq = avxsq+ dens(i)**2
              avysq = avysq+ value**2
              avxy = avxy+ dens(i)*value
              ctmask = ctmask+ 1
	    else
              if (nerr2 .lt. 10) then
                call prompt (' Interpolation error 2')
                call ivalut (' Atom :',1,i)
                call fvalut (' XYZ  :',3,xyz(1,i))
                call fvalut (' NCS  :',3,x)
                call fvalut (' X1   :',3,x1)
                call logiut (' BIT  ?',1,bit)
              else if (nerr2 .eq. 10) then
                call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
            end if
100     continue
c
        f = float (ctmask)
        ccoef = (avxy/f- avx*avy/(f*f))/ 
     $    (sqrt(avxsq/f- (avx/f)**2)* sqrt(avysq/f- (avy/f)**2))
c
        gave = gave + ccoef
        gsdv = gsdv + ccoef*ccoef
c
        if (ccoef .gt. minsco) then
          igk = lclose (minsco,scores,mxkeep,0.001,xdum)
          if (igk.ge.1 .and. igk.le.mxkeep) then
            call ivalut (' NEW in the top :',1,mxkeep)
            call fvalut (' Trans vector :',3,rt(10))
            call fvalut (' Correlation coefficient >>>',1,ccoef)
            do l=1,12
              rtkeep(l,igk)=rt(l)
            end do
            scores (igk) = ccoef
            minsco = 2.0
            do l=1,mxkeep
              minsco = min (minsco,scores(l))
            end do
          end if
        end if
c
        if (ccoef .gt. peak) then
          call prompt (' ... NEW MAXIMUM ...')
          call jvalut (' Trial number :',1,ngk)
          call fvalut (' Trans offset :',3,tdel)
          call fvalut (' Trans vector :',3,rt(10))
          call ivalut (' Nr of atoms checked :',1,ctmask)
          call fvalut (' Correlation coefficient >>>',1,ccoef)
          call vecrtv (xmo, xgk, 1, rt(1), rt(10))
          call fvalut (' NCS COG :',3,xgk)
          do 330 l=1,12
330         peakrt(l) = rt(l)
          peak = ccoef
	end if
        call flusho (6)
300   continue
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
        write (*,*)
        call jvalut (' Total number of trials :',1,ngk)
        call fvalut (' Maximum corr coeff     :',1,peak)
        gave = gave / float(ngk)
        call rvalut (' Average corr coeff     :',1,gave)
        gsdv = sqrt ( (gsdv/float(ngk)) - (gave*gave) )
        call rvalut (' Standard deviation     :',1,gsdv)
        write (*,'(/a/)') ' *** TOP SOLUTION ***'
        do 1400 i=1,12
 1400     rtbtoa(i) = peakrt(i)
        write (*,20) peakrt, peak
c
        call anancs (1,peakrt,.true.,errcod)
c
      write (*,'(/a/)') ' *** TOP SCORES (UNSORTED) ***'
      do i=1,mxkeep
        if (scores(i) .gt. 0.0) then
          write (*,29) i,scores(i),(rtkeep(j,i),j=1,12)
        end if
      end do
      write (*,*)
c
c ... 950216 - create operator file
c
      call stamp (line)
      write (2,'(a1,1x,a)') '!',line(1:leng1(line))
c
      write (2,'(a1)') '!'
      write (line,*) '! NCS6D best operator with CC = ',peak
      call pretty (line)
      write (2,'(a)') line(1:leng1(line))
      write (2,'(a)') '.LSQ_RT_NCS6D R 12 (3F15.8)'
      write (2,'(3f15.8)') (rtbtoa(i),i=1,12)
c
      write (2,'(a1)') '!'
      write (line,*) '! Top ',mxkeep,' solutions (UNSORTED)'
      call pretty (line)
      write (2,'(a)') line(1:leng1(line))
c
      write (2,'(a1)') '!'
      do i=1,mxkeep
        write (2,'(a1)') '!'
        write (line,*) '! Solution # ',i,' with CC = ',scores(i)
        call pretty (line)
        write (2,'(a)') line(1:leng1(line))
        write (2,'(a)') '.LSQ_RT_OTHER R 12 (3F15.8)'
        write (2,'(3f15.8)') (rtkeep(j,i),j=1,12)
      end do
c
      write (2,'(a1)') '!'
      close (2)
c
      return
c
10    format (' Shift=',3f12.6, ' | Corr coeff=',f12.6)
20    format (' Best Rotation Matrix'/3(1x,3f12.6/),
     +  ' Best Translation'/1x,3f12.6/
     +  ' Correlation Coefficient = ',f12.6)
25    format (' Best Rotation Matrix'/3(1x,3f12.6/),
     +  ' Best Translation'/1x,3f12.6/
     +  ' Predicted Correlation Coefficient = ',f12.6)
29    format (/' Solution nr ',i3,' Corr coeff = ',f8.5/
     +         ' Rotation Matrix'/3(1x,3f12.6/),
     +         ' Translation'/1x,3f12.6)
c
      end
c
c
c
      subroutine pdboni (typ,iunit,maxatm,natoms,xyz,ierr)
c
      implicit none
c
      integer iunit,maxatm,natoms,ierr,i,jerr,nopt,j
c
      real xyz(3,maxatm)
c
      character typ*1,line*80,optpar(4)*80
c
code ...
c
      ierr = -1
      natoms = 0
      if (typ .eq. 'B') goto 1000
c
200   continue
        read (iunit,'(a)',err=6900,end=210) line
        if (line(1:6) .ne. 'ATOM  ') goto 200
        natoms = natoms+1
        if (natoms .gt. maxatm) then
          call errcon ('Too many atoms - rest skipped')
          call jvalut (' Max nr of atoms :',1,maxatm)
          ierr = 0
          return
        end if
        read (line, '(30x,3f8.3)')  (xyz(i,natoms),i=1,3)
      goto 200
c
210   continue
      call jvalut (' Number of atoms :',1,natoms)
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
      natoms = natoms / 3
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
        subroutine ave8 (a, sizx, sizy, sizz, x, value, errcod)
c ---   average 8 nearest grid points
        implicit none
        integer errcod, sizx, sizy, sizz
        real a(sizx, sizy, sizz), x(3), value,tot
        integer i,j,k,i0,j0,k0,nn
c
code ...
c
        errcod = 1
        i0 = int(x(1))
        j0 = int(x(2))
        k0 = int(x(3))
        tot = 0.0
        nn = 0
        do 10 i=i0,i0+1
        do 10 j=j0,j0+1
        do 10 k=k0,k0+1
        if (i .le. 0) goto 10
        if (j .le. 0) goto 10
        if (k .le. 0) goto 10
        if (i .gt. sizx) goto 10
        if (j .gt. sizy) goto 10
        if (k .gt. sizz) goto 10
        nn = nn + 1
        tot = tot + a(i,j,k)
   10   continue
c
        if (nn.gt.0) then
          value = tot/float(nn)
          errcod = 0
        else
          value = 0
          errcod = -1
        end if
c
        return
        end
c
c
c
      subroutine matana (q)
c
c ... fudged version specifically for CELLO
c
      implicit none
c
      real twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      real q(3,3)
      real r(3,3),det,det3,theta,dtheta,v(3),vlen,p(3),e(3)
      real xsign,qqq
c
      integer i,j
c
code ...
c
c ... transpose alwyn's rotation matrix
c
      do i=1,3
        do j=1,3
          r(i,j)=q(j,i)
ccc          r(i,j)=q(i,j)
        end do
      end do
c
      det = det3(r)
      write (*,'(a,f10.7)') ' Determinant  : ',det
c
c ... find POLAR angles from matrix
c
      qqq=max(-1.0,min(1.0,((r(1,1)+r(2,2)+r(3,3)-1)*0.5)))
      theta = acos (qqq)
      dtheta = theta*rtodeg
      write (*,'(a,f10.4)') ' Net rotation : ',dtheta
c
c ... changed 940316
c
      xsign = 1.0
      if (dtheta .gt. 180.0) then
        dtheta = 360.0 - dtheta
        xsign = -1.0
      end if
c
      if (dtheta .eq. 0.0) then
        write (*,'(a)')
     +    ' ERROR - indeterminate direction cosines'
        goto 200
      else if (abs(dtheta-180.0) .lt. 1.0e-5) then
        v(1) = r(1,1) + 1.0
        v(2) = r(2,1)
        v(3) = r(3,1)
        vlen = sqrt (v(1)**2 + v(2)**2 + v(3)**2)
        do i=1,3
          v(i) = v(i)/vlen
        end do
        write (*,'(a,3f10.4)')
     +    ' Direction cosines of rotation axis : ',
     +    (v(i),i=1,3)
      else
c
c ... changed 940316
c
        vlen = xsign * 0.5 / sin(theta)
        v(1) = r(3,2)-r(2,3)
        v(2) = r(1,3)-r(3,1)
        v(3) = r(2,1)-r(1,2)
        do i=1,3
          v(i) = v(i)*vlen
        end do
        write (*,'(a,3f10.4)')
     +    ' Direction cosines of rotation axis : ',
     +    (v(i),i=1,3)
      end if
c
c ... changed 940316
c
      det = sqrt( abs(1.0-v(3)*v(3)) )
      if (v(3).eq.0.0) then
        p(1) = 9999.999
      else
        p(1) = atan2 (det,v(3)) * rtodeg
      end if
ccc      p(1) = acos (v(3)) * rtodeg
      if (v(1).eq.0.0) then
        p(2) = 9999.999
      else
        p(2) = atan2 (v(2),v(1)) * rtodeg
      end if
      p(3) = dtheta
c
      write (*,'(a,3f10.4)')
     +  ' ALMN Polar angles Omega/Phi/Chi    : ',(p(i),i=1,3)
c
c ... find EULER angles from matrix
c
  200 continue
c
      if (r(1,3).eq.0.0) then
        e(1) = 9999.999
      else
        e(1) = atan2 (r(2,3),r(1,3)) * rtodeg
      end if
cc      print *,'E1 ',e(1)
c
      r(3,3)=max(-1.0,min(1.0,r(3,3)))
      e(2) = acos (r(3,3)) * rtodeg
cc      print *,'E2 ',e(2)
c
      if (r(3,1).eq.0.0) then
        e(3) = 9999.999
      else
        e(3) = atan2 (r(3,2),-r(3,1)) * rtodeg
      end if
cc      print *,'E3 ',e(3)
c
      write (*,'(a,3f10.4)')
     +  ' ALMN Euler angles Alpha/Beta/Gamma : ',(e(i),i=1,3)
c
      return
      end
