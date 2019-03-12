      program spancsi
c
c --- SPANCSI = SPAnish NCS Inquisition
c
c ... Gerard J Kleywegt @ 970522
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'SPANCSI', vers = '080625/1.0.2')
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
c      pointer (idptr,mapc)
c
c      real mapa(1), mapb(1), mapc(1)
c      integer maskb(1), malloc
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
      idptr = fmalloc (nb)
c
      nb = wrdbyt*masksize
      icptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call dospan (%val(iaptr),%val(ibptr),%val(idptr),
     +             %val(icptr), mapsize, masksize)
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
      subroutine dospan (mapa,mapb,mapc,maskb,maxsiz,maxmsk)
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz,maxsym,maxncs,maxmsk
c      parameter (maxsiz = maxgk1)
c      parameter (maxmsk = maxgk1)
      parameter (maxsym = maxgk2)
      parameter (maxncs = maxgk2)
c
c ---      Map A data structures
c
      integer orgna(3), exta(3), grida(3), uvwa(3)
      real mapa(maxsiz), cella(6)
c
c ---      Map B data structure (ditto for Map C)
c
      integer orgnb(3), extb(3), gridb(3), maskb(maxmsk)
      real mapb(maxsiz), mapc(maxsiz), cellb(6)
c
      integer orgnc(3), extc(3), gridc(3)
      real cellc(6)
c
      real cell(6), grid(3), rtatob(12,maxncs), varian(maxncs)
      real rtbtoa(12,maxncs), rtsym(12,maxsym)
      real avy(maxncs),avysq(maxncs),avxy(maxncs),avxmy(maxncs)
      real avea,sdva,xmin,xmax,xtot,x1,x2,x3,x4,x5,x6,x7,x8,x9
      real varmap,factor,avrho,back
c
      character file*80,fmt*80,par*25,partyp*1
      character task*1,fmapa*80,fmapb*80,fmapc*80,fmask*80
c
      integer ctrt, ctsym, errcod, i, j, spgrp, ierr, length
      integer incs, k, l, ioff, jncs, ctin
c
      logical xinter,check,lposit,lback,lprint,lvari
c
code ...
c
c ... 950118 - check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
c      call jvalut (' Max size of maps and masks    :',1,maxsiz)
      call jvalut (' Max nr of spacegroup symm ops :',1,maxsym)
      call jvalut (' Max nr of NCS symm ops        :',1,maxncs)
c
      task = 'B'
      ctrt = 0
      fmapa = ' '
      fmapb = ' '
      fmapc = ' '
      fmask = ' '
c
      lposit = .false.
      lback  = .true.
c
c ... get map
c
      write (*,*)
      call textin (' Map to be used ?',fmapa)
      call textut (' Map to be used :',fmapa)
c
      if (length(fmapa) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
      call edin (fmapa, 1, mapa, orgna, exta, grida, uvwa,
     $           cella, spgrp, maxsiz)
      close (1)
      call telmap (grida,orgna,exta,cella)
c
c ... get mask
c
      write (*,*)
      call textin (' Mask file ?',fmask)
      call textut (' Mask file :',fmask)
      if (length(fmask) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
      call xopxoa (3,fmask,xinter(),errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening mask')
        return
      end if
      call maskin 
     $     (3, maskb, orgnb, extb, gridb, cellb, maxmsk, ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading mask')
        return
      end if
      close (3)
      call telmap (gridb,orgnb,extb,cellb)
c
      check = .true.
      do i=1,3
        check = (check .and. (grida(i) .eq. gridb(i)))
        check = (check .and. (abs(cella(i)-cellb(i)) .le. 0.01))
        check = (check .and.
     +          (abs(cella(i+3)-cellb(i+3)) .le. 0.01))
      end do
      if (.not. check) then
        call errcon ('Map should have same grid etc. as mask')
        return
      end if
c
      write (*,*)
      call textin (' Output CCP4 map (if required) ?',fmapb)
      call textut (' Output CCP4 map :',fmapb)
      if (length(fmapb) .lt. 1) then
        call errcon ('No filename provided')
        return
      end if
c
      do 100 i=1,3
100     grid(i) = cella(i)/float(grida(i))
      do 110 i=1,6
110     cell(i) = cella(i)
c
      file = 'symop.o'
      write (*,*)
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
c --- Get rt b->a
c
210   continue
      file = ' '
      write (*,*)
      call textin (' File with NCS RT-operator(s) (RETURN to quit) ?',
     +  file)
      call textut (' File with NCS RT-operator(s) (RETURN to quit) :',
     +  file)
c
      if (file .eq. ' ') goto 200
c
c ... read one or many NCS operators
c
      call rdoncs (1,file,ctrt,maxncs,rtbtoa,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading NCS operators')
        return
      end if
c
cc      call opoodb (1,file,par,partyp,j,fmt,errcod)
cc      if (errcod .ne. 0) call errstp (
cc     +  'While opening O datablock file')
cc
cc      ctrt = ctrt+1
cc      if (ctrt .gt. maxncs) then
cc        call errstp ('Too many NCS operators')
cc      end if
cc      read (1, fmt,err=999,end=999) (rtbtoa(i, ctrt),i=1,j)
cc      close (1)
c
      goto 210
c
c --- Now go do everything. This is to make subscripting easier
c
200   continue
      write (*,*)
      if (ctrt .lt. 2) then
        call errcon ('No NCS operators')
        return
      end if
c
cc      call anancs (ctrt,rtbtoa,.true.,ierr)
cc      if (ierr .ne. 0)
cc     +  call errstp ('In NCS operators')
c
      task = 'A'
      lposit = .false.
      lprint = .false.
      lback = .false.
      lvari = .false.
      do i=1,ctrt
        varian (i) = -999.999
      end do
c
      l = exta(1)*exta(2)*exta(3)
      call xstats (mapa,l,avea,sdva,xmin,xmax,xtot)
      varmap = sdva * sdva
      call rvalut (' Variance of entire map :',1,varmap)
c
c ... event loop
c
 6969 continue
c
      write (*,'(10(1x,a,/))')
     +  ' ',
     +  'Select one of the following options:',
     +  ' ',
     +  'A = Analyse density for individual NCS molecules',
     +  'C = Correlate density for pairs of NCS molecules',
     +  'B = Both average and expand using variance scaling',
     +  'Q = Quit this program'
c
      call textin (' Option (A/C/S/Q) ?',task)
      call upcase (task)
c
      if (task .eq. 'A') then
c
c ... loop over NCS operators, get density, calc & print statistics
c
        write (*,6100) 'NCS','Aver. dens.','St. dev.','Variance',
     +    'Minimum','Maximum','Sum dens.'
c
 6100 format (/1x,a3,6a13)
 6110 format (1x,i3,1p,6(1x,e12.4))
c
        k = extb(1)*extb(2)*extb(3)
c
        do incs = 1,ctrt
c
          call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $      mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $      cell, grid, rtbtoa(1,incs), 1, rtsym, ctsym, ierr,
     +      avy,avysq,avxy,avxmy,lposit,lprint)
c
c ... copy only masked points to MAPC
c
          l = 0
          do i=1,k
            if (maskb(i) .eq. 1) then
              l = l + 1
              mapc (l) = mapb (i)
            end if
          end do
c
ccc          print *,'#pts, #masked ',k,l
c
          call xstats (mapc,l,avea,sdva,xmin,xmax,xtot)
          write (*,6110) incs,avea,sdva,sdva*sdva,xmin,xmax,xtot
c
          if (.not. lvari) varian (incs) = sdva*sdva
c
        end do
c
        lvari = .true.
c
      else if (task .eq. 'C') then
c
        write (*,6200) 'NCSi','NCSj','RMSD','Corr.coeff.',
     +    'R (wrt NCSi)','R (wrt NCSj)'
c
 6200 format (/2a6,4a13)
 6210 format (2i6,1p,(1x,e12.4),0p,3(1x,f12.3))
c
        k = extb(1)*extb(2)*extb(3)
c
        do incs=2,ctrt
c
ccc          write (*,*)
c
          call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $      mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $      cell, grid, rtbtoa(1,incs), 1, rtsym, ctsym, ierr,
     +      avy,avysq,avxy,avxmy,lposit,lprint)
c
          l = 0
          do i=1,k
            if (maskb(i) .eq. 1) then
              l = l + 1
              mapc (l) = mapb (i)
            end if
          end do
c
          if ((.not. lvari) .and. varian(incs).lt.0.0) then
            call xstats (mapc,l,avea,sdva,xmin,xmax,xtot)
            varian (incs) = sdva*sdva
          end if
c
          ioff = l
          if (2*ioff .gt. maxsiz) then
            call errcon ('Not enough memory for operation !')
            goto 6969
          end if
c
          do jncs=1,incs-1
c
            call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $        mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $        cell, grid, rtbtoa(1,jncs), 1, rtsym, ctsym, ierr,
     +        avy,avysq,avxy,avxmy,lposit,lprint)
c
            l = ioff
            do i=1,k
              if (maskb(i) .eq. 1) then
                l = l + 1
                mapc (l) = mapb (i)
              end if
            end do
c
            if ((.not. lvari) .and. varian(jncs).lt.0.0) then
              call xstats (mapc(ioff+1),ioff,avea,sdva,xmin,xmax,xtot)
              varian (jncs) = sdva*sdva
            end if
c
            call xystat (mapc(1),mapc(ioff+1),ioff,
     +                   x1,x2,x3,x4,x5,x6,x7,x8,x9)
c
            write (*,6210) incs,jncs,x1,x3,x4,x5
c
          end do
c
        end do
c
        lvari = .true.
c
        write (*,*)
c
      else if (task .eq. 'B') then
c
        k = extb(1)*extb(2)*extb(3)
c
c ... if variance not calculated yet, do it now
c
        if (.not. lvari) then
          write (*,6100) 'NCS','Aver. dens.','St. dev.','Variance',
     +      'Minimum','Maximum','Sum dens.'
          do incs = 1,ctrt
            call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $        mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $        cell, grid, rtbtoa(1,incs), 1, rtsym, ctsym, ierr,
     +        avy,avysq,avxy,avxmy,lposit,lprint)
            l = 0
            do i=1,k
              if (maskb(i) .eq. 1) then
                l = l + 1
                mapc (l) = mapb (i)
              end if
            end do
            call xstats (mapc,l,avea,sdva,xmin,xmax,xtot)
            write (*,6110) incs,avea,sdva,sdva*sdva,xmin,xmax,xtot
            varian (incs) = sdva*sdva
          end do
          lvari = .true.
        end if
c
c ... zero output map
c
        write (*,*)
        do i=1,k
          mapc (i) = 0.0
        end do
c
c ... get density for each molecule; scale by variance
c
        do incs = 1,ctrt
c
          call jvalut (' Add density for operator :',1,incs)
c
          call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $      mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $      cell, grid, rtbtoa(1,incs), 1, rtsym, ctsym, ierr,
     +      avy,avysq,avxy,avxmy,lposit,lprint)
c
c ... add masked points to MAPC
c
          factor = varian(1)/varian(incs)
c
          call fvalut ('             Scale factor :',1,factor)
c
          l = 0
          do i=1,k
            mapc (i) = mapc (i) + factor*mapb (i)
          end do
c
        end do
c
c ... average
c
        call prompt (' Average density')
c
        factor = 1.0 / float(ctrt)
c
        do i=1,k
          mapc (i) = factor * mapc (i)
        end do
c
c ... zero output map
c
        call prompt (' Zero output map')
c
        do i=1,exta(1)*exta(2)*exta(3)
          mapa (i) = 0.0
        end do
c
c ... expand each molecule; scale by variance
c
        do incs=1,ctrt
c
          call jvalut (' Expand for operator :',1,incs)
c
          factor = varian(incs)/varian(1)
c
          call fvalut ('        Scale factor :',1,factor)
c
          do i=1,k
            mapb (i) = factor * mapc (i)
          end do
c
c ... call expansion routine
c
          call sub02x (mapa, exta(1), exta(2), exta(3), orgna, grida,
     $      mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $      cell, grid, rtbtoa(1,incs), rtatob(1,incs), 1,
     +      rtsym, ctsym, lback, lprint)
c
        end do
c
c ... flatten solvent
c
        call prompt (' Flatten solvent')
c
        avrho = 0.
        ctin = 0
        k = exta(1)*exta(2)*exta(3)
        do i=1,k
          if (mapa(i) .ne. 0.0) then
            avrho = avrho + mapa(i)
            ctin = ctin + 1
          end if
        end do
        back = -avrho/float(k - ctin)
        avrho = avrho/float(ctin)
c
        call jvalut (' Nr of points in map  :',1,k)
        call jvalut (' Nr of masked points  :',1,ctin)
        call jvalut (' Nr of solvent points :',1,(k-ctin))
        factor = 100.0 * float(k-ctin) / float(k)
        call fvalut (' Solvent fraction (%) :',1,factor)
c
        call rvalut (' Average density inside masks :',1,avrho)
        call rvalut (' Average density in solvent   :',1,back)
        call prompt (' Setting background ...')
c
        avrho = 0.
        do i=1,k
          if (mapa(i) .eq. 0.0) mapa(i) = back
          avrho = avrho + mapa(i)
        end do
        avrho = avrho/float(k)
        call rvalut (' Average density overall      :',1,avrho)
c
c ... write new map
c
        write (*,*)
        call edout 
     $    (fmapb, 4, mapa, orgna, exta, grida, uvwa, cella, spgrp)
c
      else if (task .eq. 'Q') then
c
        goto 6868
c
      else
c
        call errcon ('Invalid option !')
c
      end if
c
      goto 6969
c
 6868 continue
c
c      call gkquit
c
      return
c
c ... error traps
c
  999 continue
      call errcon ('While reading file')
      return
c
      end
c
c
c
      subroutine sub01x 
     $  (mapa, exta1, exta2, exta3, orgna, 
     $   mapb, maskb, extb1, extb2, extb3, orgnb, 
     $   cell, grid, rtbtoa, ctrt, rtsym, ctsym, ierr,
     +   avy, avysq, avxy, avxmy, lposit, lprint)
c
c ---  Build up a map. Given map a, map b and mask b, and rt_b_to_a
c      For every point in mask b, transform to map a, interpolate
c      onto map a . Add mapa value to map b.
c ---  Assumes both maps are on same grid.
c ---  Does not require a unit operator, i.e. uses map b value.
c ---  Alwyn Jones, 21-Feb-91
c
      implicit none
c
      real small
      parameter (small=1.0e-9)
c
c ---      Map A data structures
      integer orgna(3), exta1, exta2, exta3
      real mapa(exta1, exta2, exta3)
c ---      Map B data structure
      integer orgnb(3), extb1, extb2, extb3
      integer maskb(extb1, extb2, extb3)
      real mapb(extb1, extb2, extb3)
      real cell(6), grid(3), rtbtoa(12,*), rtsym(12,*)
      integer ctrt, ctsym,ierr,nbit,ngjk,bobo1,bobo2,bobo3
c
      real avx,avxsq,avy(*),avysq(*),avxy(*),avxmy(*)
      real avax,val1,f,cc
c
      integer errcod, i, j, k, l, loop, m, ext(3),extb(3)
      integer nerr1,nerr2
c
      real a(3,3), b(3,3), forgn(3), fext(3), mapbit(4,4,4),
     $     value, x(3), x1(3), x2(3),gext(3),xdum
c
      logical lposit,lprint
c
code ...
c
      ierr = -1
      nbit = 0
      if (lprint) call prompt (' Average ...')
c
c ---      A = frac to ang, B=Ang to frac
c
      call orthog (cell, a, 0)
      call orthog (cell, b, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
      extb(1) = extb1
      extb(2) = extb2
      extb(3) = extb3
c
c ---      Generate envelope edges in fractional cell coords
c
      do 130 i=1,3
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
130     gext(i)  = float(ext(i)+orgna(i)-2)*grid(i)/cell(i)
c
      if (lprint) then
        call fvalut (' Map  Cell   :',6,cell)
        call fvalut (' Spacing     :',3,grid)
        call ivalut (' Map  origin :',3,orgna)
        call ivalut (' Map  extent :',3,ext)
        call ivalut (' Mask origin :',3,orgnb)
        call ivalut (' Mask extent :',3,extb)
        call fvalut (' FORGN :',3,forgn)
        call fvalut (' FEXT  :',3,fext)
        call fvalut (' GEXT  :',3,gext)
      end if
c
c ---      Loop over mask looking for something
c
      avx   = 0.0
      avax  = 0.0
      avxsq = 0.0
      do loop=1,ctrt
        avy(loop) = 0.0
        avysq(loop) = 0.0
        avxy (loop) = 0.0
        avxmy(loop) = 0.0
      end do
      ngjk = 0
      nerr1 = 0
      nerr2 = 0
c
      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
      if (lprint) call jvalut (' Points in mask :',1,bobo1)
      bobo2 = bobo1/10
      bobo3 = bobo2
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        mapb(i,j,k) = 0.
        if (maskb(i,j,k) .eq. 1) then
          ngjk = ngjk + 1
c
          if (ngjk .eq. bobo3 .and. lprint) then
            xdum = 100.0 * float(ngjk-1)/float(bobo1)
            call fvalut (' Progress (% mask) :',1,xdum)
            bobo3 = bobo3 + bobo2
          end if
c
          x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
          x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
          x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
          call mulmtx (a, x, x2, 3, 3, 1)
c
c ---     Now loop over the number of operators
c
          do 120 loop=1, ctrt
            call vecrtv (x2, x, 1, rtbtoa(1,loop), rtbtoa(10,loop))
            call mulmtx (b, x, x1, 3, 3, 1)
c ---            Is it in the envelope, if not make it (if possible)
            call frcsym (x1, forgn, gext, rtsym, ctsym, errcod)
c
            if (errcod .ne. 0) then
c
c ---       Unable to do it. Probably around the edge so build a 4x4x4 map
              nbit = nbit + 1
              call bldbit (mapa, exta1, exta2, exta3, orgna, grid,
     $          cell, forgn, fext, rtsym, ctsym, mapbit, x1, errcod)
c
              if (errcod .ne. 0) then
c ---      Still unable to do it.
                if (nerr1 .lt. 10) then
                  call errcon (' Severe FRCSYM error')
                  call fvalut (' Mask point :',3,x2)
                  call fvalut (' NCS point  :',3,x1)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 120
              end if
c
              do 125 l=1,3
                x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
                m = x(l)
                x(l) = x(l)- float(m-1)
125           continue
              call intrpl (mapbit, 4, 4, 4, x, value, errcod)
            else
              do 110 l=1,3
110             x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l)) + 1
c
c ---           x is now set of indices to the mapa array
c
              call intrpl 
     $          (mapa, exta1, exta2, exta3, x, value, errcod)
            end if
            if (errcod .eq. 0) then
              mapb(i,j,k) =  mapb(i,j,k) + value
c
              if (loop .eq. 1) then
                val1 = value
                avx   = avx + value
                avax  = avax + abs(value)
                avxsq = avxsq + value*value
              end if
c
ccc              write (*,'(3i4,3f8.2,1pe12.4)') i,j,k,x,value
c
              avy (loop)   = avy(loop) + value
              avysq (loop) = avysq(loop) + value*value
              avxy (loop)  = avxy(loop) + value*val1
              avxmy (loop) = avxmy(loop) + abs(value-val1)
c
            else
              call prompt (' Interpolation error !')
              write (*,'(1x,a,1x,3f8.2,1x,3i4)')
     +          'x & ext',x,exta1,exta2,exta3
              nerr2 = nerr2 + 1
c
            end if
120       continue
        else
        end if
100   continue
c
      ierr = 0
      if (lprint) then
        call jvalut (' Calls to BLDBIT      :',1,nbit)
        call jvalut (' Severe FRCSYM errors :',1,nerr1)
        call jvalut (' Interpolation errors :',1,nerr2)
      end if
c
c ... print correlation coefficients for the various operators
c
      if (lprint) call jvalut (' Nr of mask points :',1,ngjk)
      f = float (ngjk)
      do loop=1,ctrt
        cc = (avxy(loop)/f - avx*avy(loop)/(f*f))/ 
     +    (sqrt(avxsq/f - (avx/f)**2) *
     +     sqrt(avysq(loop)/f - (avy(loop)/f)**2))
        if (lprint) write (*,'(a,i3,a,f8.5)')
     +    ' Corr. coeff. for operator ',loop,' = ',cc
        cc = avxmy(loop) / avax
        if (lprint) write (*,'(a,i3,a,f8.5)') 
     +    ' R-factor for operator ',loop,' w.r.t. operator 1 = ',cc
      end do
c
c --- Normalize by number of additions
c
      if (lposit) then
        if (lprint)
     +    call rvalut (' Positivity; set <= 0 to:',1,small)
        xdum = 1.0 / float(ctrt)
        do 300 k=1,extb3
        do 300 j=1,extb2
        do 300 i=1,extb1
c ... next line changed for HP
          if (maskb(i,j,k) .eq. 1) then
            if (mapb(i,j,k) .gt. 0.0) then
              mapb(i,j,k) = mapb(i,j,k) * xdum
            else
              mapb(i,j,k) = small
            end if
          end if
  300   continue
      else
        if (lprint)
     +    call prompt (' Averaging without positivity constraint')
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
c
c
c
      subroutine sub02x 
     $  (mapa, exta1, exta2, exta3, orgna, grida,
     $   mapb, maskb, extb1, extb2, extb3, orgnb, 
     $   cell, grid, rtbtoa, rtatob, ctrt, rtsym, ctsym,
     +   lback, lprint)
c
c ---  Build up a map. Given map a, map b and mask b, and rt_b_to_a
c      For every point in mask b, transform to map a, interpolate
c      onto map a . Add mapa value to map b.
c ---  Assumes both maps are on same grid.
c ---  Does not require a unit operator, i.e. uses map b value.
c ---  Alwyn Jones, 21-Feb-91
c
      implicit none
c
c ---      Map A data structures
      integer orgna(3), exta1, exta2, exta3, grida(3)
      real mapa(exta1, exta2, exta3)
c ---      Map B data structure
      integer orgnb(3), extb1, extb2, extb3
      integer maskb(extb1, extb2, extb3)
      real mapb(extb1, extb2, extb3)
      real cell(6), grid(3), rtatob(12,*), rtbtoa(12,*), rtsym(12,*)
      integer ctrt, ctsym,bobo1,bobo2,bobo3
c
      integer ctin, errcod, i, ina(3), j, k, l, loop,nasu
      integer i1, j1, k1, ijk1(3),ngk1,ngk2,nerr1,nerr2,nerr3
      integer ext(3), i2, j2, k2, ijk2(3)
c
      integer kk1,kk2,kk3,nmask,nout,ntwice,ijk3(3)
c
      real a(3,3), avrho, b(3,3), back, forgn(3), fext(3), 
     $     value, x(3), x1(3), x2(3), xdum
c
      logical lback, lprint
c
      equivalence (ijk1(1), i1), (ijk1(2), j1), (ijk1(3), k1)
      equivalence (ijk2(1), i2), (ijk2(2), j2), (ijk2(3), k2)
c
code ...
c
      if (lprint) call prompt (' Expand ...')
c
c ---      A = frac to ang, B=Ang to frac
      call orthog (cell, a, 0)
      call orthog (cell, b, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
c ---      Generate envelope edges in fractional cell coords
      do 340 i=1,3
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
340     fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
c
c ---      Get reverse transformations
      do 300 i=1,ctrt
        do 310 j=1,9
310       rtatob(j,i) = rtbtoa(j,i)
        call matinv (rtatob(1,i), 3, x, x1, x2)
        call mulmtx (rtatob(1,i), rtbtoa(10,i), rtatob(10,i), 3, 3, 1)
        do 320 j=10,12
320       rtatob(j,i) = -rtatob(j,i)
300   continue
c
c .. find out how many asymmetric units map A comprises
c
      i = grida(1)*grida(2)*grida(3)
      k = i/ctsym
      j = exta1*exta2*exta3
      xdum = float(j)/float(k)
      nasu = nint (xdum)
      if (abs(xdum-float(nasu)) .gt. 0.2) then
        call jvalut (' Nr of points in unit cell  :',1,i)
        call jvalut (' Nr of symmetry operators   :',1,ctsym)
        call jvalut (' Nr of points in asymm unit :',1,k)
        call jvalut (' Nr of points in map        :',1,j)
        call fvalut (' Nr of asymm units in map   :',1,xdum)
        call jvalut (' Nearest integer            :',1,nasu)
        call prompt (' WARNING - Not an integer number !!')
c        return
      end if
      if (lprint) then
        call fvalut (' Nr of asymm units in map   :',1,xdum)
        call ivalut (' Nr of asymm units in map   :',1,nasu)
      end if
      if (xdum .lt. 1.0) then
        call errcon ('Less than one asymmetric unit')
        return
      end if
c
      ngk1 = 0
      ngk2 = 0
      nerr1 = 0
      nerr2 = 0
      nerr3 = 0
      nout = 0
      ntwice = 0
c
c --- Zero the A array
c
ccc      do 330 k=1,exta3
ccc      do 330 j=1,exta2
ccc      do 330 i=1,exta1
ccc330      mapa(i,j,k) = 0.
c
      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
      if (lprint) call jvalut (' Points in mask :',1,bobo1)
      bobo2 = bobo1/10
      bobo3 = bobo2
c
c --- Loop over mask looking for something
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        if (maskb(i,j,k) .eq. 1) then
          ngk1 = ngk1 + 1
c
          if (ngk1 .eq. bobo3 .and. lprint) then
            xdum = 100.0 * float(ngk1-1)/float(bobo1)
            call fvalut (' Progress (% mask) :',1,xdum)
            bobo3 = bobo3 + bobo2
          end if
c
          x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
          x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
          x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
          call mulmtx (a, x, x2, 3, 3, 1)
c ---          Now loop over the number of operators
          do 120 loop=1, ctrt
            call vecrtv (x2, x, 1, rtbtoa(1,loop), rtbtoa(10,loop))
            call mulmtx (b, x, x1, 3, 3, 1)
            do 110 l=1,3
110           x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c
c ---       x is now set of indices to the mapa array
c ---       Now back interpolate for each of the 8 integral grid points
c           in map a nearest this
c
            do 130 l=1,3
130              ina(l) = x(l)
            do 150 k1=ina(3),ina(3)+1
            do 150 j1=ina(2),ina(2)+1
            do 150 i1=ina(1),ina(1)+1
c
              call  frctrn 
     $              (ijk1, ijk2, orgna, grida, forgn, fext, 
     $              rtsym, ctsym, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call errcon ('FRCTRN error')
                  call ivalut (' Mask indices :',3,ijk1)
                  call ivalut (' Map  indices :',3,ijk2)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCTRN errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 150
              end if
c
c ... check if point already set
c
              if (mapa(i2,j2,k2) .ne. 0.) goto 150
c
              x(1) = (i1-1+ orgna(1))*grid(1)/cell(1)
              x(2) = (j1-1+ orgna(2))*grid(2)/cell(2)
              x(3) = (k1-1+ orgna(3))*grid(3)/cell(3)
c
              call mulmtx (a, x, x1, 3, 3, 1)
              call vecrtv (x1, x, 1, rtatob(1,loop), rtatob(10,loop))
              call mulmtx (b, x, x1, 3, 3, 1)
c
              do 140 l=1,3
140             x(l) = x1(l)*cell(l)/grid(l)- float(orgnb(l))+ 1
c
c ---         x is now set of indices to the mapb array
c
c ... check whether inside or on edge of mask
c
              do l=1,3
                ijk3(l) = x(l)
              end do
              nmask = 0
c
              do kk1=ijk3(1),ijk3(1)+1
                if (kk1 .gt. 0 .and. kk1 .le. extb1) then
                  do kk2=ijk3(2),ijk3(2)+1
                    if (kk2 .gt. 0 .and. kk2 .le. extb2) then
                      do kk3=ijk3(3),ijk3(3)+1
                        if (kk3 .gt. 0 .and. kk3 .le. extb3) then
                          if (maskb(kk1,kk2,kk3) .eq. 1)
     $                      nmask = nmask + 1
                        end if
                      end do
                    end if
                  end do
                end if
              end do
c
c ... at least one neighbour point must be part of the mask
c
              if (nmask .lt. 1) then
c                call ivalut (' Nr of mask nbrs:',1,nmask)
c                call ivalut (' Mask point    :',3,ijk1)
c                call ivalut (' NCS point     :',3,ijk2)
c                call ivalut (' Mask indices  :',3,ijk3)
                nout = nout + 1
                goto 150
              end if
c
              call intrpl 
     $              (mapb, extb1, extb2, extb3, x, value, errcod)
              if (errcod .eq. 0) then
                call  frcval 
     $                (ijk1, mapa, value, exta1, exta2, exta3, 
     $                orgna, grida, forgn, fext, 
     $                rtsym, ctsym, nasu, ntwice, errcod)
                if (errcod .eq. 0) then
                  ngk2 = ngk2 + 1
                else
                  if (nerr3 .lt. 10) then
                    call errcon ('Serious FRCVAL error')
                    call ivalut (' Mask point :',3,ijk1)
                    call ivalut (' NCS point  :',3,ijk2)
                  else if (nerr3 .eq. 10) then
                    call prompt (
     +              ' NOTE: further FRCVAL errors but not listed!!!')
                  end if
                  nerr3 = nerr3 + 1
                end if
              else
                if (nerr2 .lt. 10) then
                  call errcon ('Interpolation error')
                  call ivalut (' Mask point :',3,ijk1)
                  call ivalut (' NCS point  :',3,ijk2)
                else if (nerr2 .eq. 10) then
                  call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
                end if
                nerr2 = nerr2 + 1
              end if
150         continue
120       continue
        end if
100   continue
c
      if (lprint) then
c
        call jvalut (' Points in mask       :',1,ngk1)
        call jvalut (' Set in asymm unit    :',1,ngk2)
        call jvalut (' Total points set     :',1,(ngk2*nasu))
c
        call jvalut (' FRCTRN errors        :',1,nerr1)
        call jvalut (' Interpolation errors :',1,nerr2)
        call jvalut (' FRCVAL errors        :',1,nerr3)
c
        call jvalut (' Points outside mask  :',1,nout)
        call jvalut (' Points set > 1 *     :',1,ntwice)
c
      end if
c
c --- Now fix up the background so that the average density is 0.0
c
      if (.not. lback) then
        if (lprint) call prompt (' Keeping background at zero')
      else
        if (lprint) call prompt (' Calculating background ...')
        avrho = 0.
        ctin = 0
        do 500 k=1,exta3
        do 500 j=1,exta2
        do 500 i=1,exta1
          if (mapa(i,j,k) .ne. 0.0) then
            avrho = avrho + mapa(i,j,k)
            ctin = ctin + 1
          end if
500     continue
        back = -avrho/float(exta1*exta2*exta3 - ctin)
        avrho = avrho/float(ctin)
c
        if (lprint) then
          call rvalut (' Average density inside masks :',1,avrho)
          call rvalut (' Average density in solvent   :',1,back)
          call prompt (' Setting background ...')
        end if
c
        avrho = 0.
        do 510 k=1,exta3
        do 510 j=1,exta2
        do 510 i=1,exta1
          if (mapa(i,j,k) .eq. 0.0) then
            mapa(i,j,k) = back
          end if
          avrho = avrho + mapa(i,j,k)
510     continue
        avrho = avrho/float(exta1*exta2*exta3)
        if (lprint) call rvalut (
     +    ' Average density overall      :',1,avrho)
c
      end if
c
      return
      end
c
c
c
      subroutine frcval (i1, map, value, e1, e2, e3, 
     +                   orgn, grid, min, max, 
     +                   rt, ct, nasu, ntwice, errcod)
c
c ---	FoRCe_VALue into the asymmetric unit
c ---	Given the grid index, force via symops
c	into the envelope. Errcod .ne. 0 if impossible
c ---	Alwyn Jones
c
      implicit none
c
      integer ct, e1, e2, e3, errcod, i1(3), grid(3), orgn(3),nasu
      real min(3), max(3), map(e1, e2, e3), rt(12,*), value
c
      integer i, i2(3), j, nok, ntwice
      real x(3), x1(3), dummy
c
code ...
c
      nok = 0
      errcod = 0
c
      do 100 i=1,3
100     x(i) = float(i1(i)+ orgn(i)- 1)/float(grid(i))
      do 200 i=1,ct
        errcod = 0
        call vecrtv (x, x1, 1, rt(1,i), rt(10,i))
        do 210 j=1,3
220       if (x1(j) .lt. min(j)) then
            x1(j) = x1(j)+ 1.
            goto 220
          end if
c
230       if (x1(j) .gt. max(j)) then
            x1(j) = x1(j)- 1.
            goto 230
          end if
c
          if (x1(j) .lt. min(j)) then
            errcod = 1
          else
            nok = nok + 1
          end if
210     continue
c
        if (errcod .eq. 0) then
          do 250 j=1,3
250         i2(j) = nint(x1(j)*float(grid(j))) - orgn(j) + 1
          if (i2(1) .le. 0 .or. i2(2) .le. 0 .or. i2(3) .le. 0) then
            goto 200
          end if
          if (i2(1) .gt. e1 .or. i2(2) .gt. e2 .or. i2(3) .gt. e3) then
            goto 200
          end if
c
          dummy = map(i2(1),i2(2),i2(3))
          if (dummy .ne. 0.0) then
            if (dummy .ne. value) then
cc              write (*,'(a,3i6,2f10.4)') ' Point set before : ',
cc     +          i2,dummy,value
              if (value.gt.dummy) map(i2(1),i2(2),i2(3)) = value
              ntwice = ntwice + 1
            end if
          else
            map(i2(1),i2(2),i2(3)) = value
          end if
c
        end if
200   continue
c
      errcod = 0
      if (nok .lt. nasu) errcod = 1
c
      return
      end
