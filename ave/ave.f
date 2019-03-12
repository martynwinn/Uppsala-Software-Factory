      program average
c
c ---  Build up a map. Given map a, map b and mask b, and a number of 
c      rt_b_to_a operators.
c      For every point in mask b, transform to map a, interpolate
c      onto map a . Add map a value to map b.
c ---  Assumes both maps are on same grid.
c ---  Alwyn Jones, 21-Feb-91
c ---  17-Sep-91, convert to ccp4
c mok 950120: Dynamic allocation of memory
c mok 961121: Average now called as subroutine.
c gjk 961122: Dyn mem alloc to work on SGI, ALPHA and ESV; support
c             MAPSIZE and MASKSIZE as environment variables and/or
c             as command-line arguments (former MUST be uppercase)
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'AVE', vers = '080625/5.0.2')
c
      include '../maxdim.incl'
c
      integer maxsiz, maxmsk
      parameter (maxsiz = maxgk1)
      parameter (maxmsk = maxgk1)
c
c      pointer (iaptr,mapa)
c      pointer (ibptr,mapb)
c      pointer (icptr,maskb)
c
c      real mapa(1), mapb(1)
c      integer maskb(1), malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr
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
      call doaver (%val(iaptr),%val(ibptr),%val(icptr),
     +             mapsize, masksize)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine doaver (mapa,mapb,maskb,maxsiz,maxmsk)
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
c ---      Map B data structure
c
      integer orgnb(3), extb(3), gridb(3), maskb(maxmsk)
      real mapb(maxsiz), cellb(6)
c
      integer orgnc(3), extc(3), gridc(3)
      real cellc(6)
c
      real cell(6), grid(3), rtatob(12,maxncs)
      real rtbtoa(12,maxncs), rtsym(12,maxsym)
      real avy(maxncs),avysq(maxncs),avxy(maxncs),avxmy(maxncs)
c
      character file*80,fmt*80,par*25,partyp*1
      character task*1,fmapa*80,fmapb*80,fmapc*80,fmask*80
c
      integer ctrt, ctsym, errcod, i, j, spgrp, ierr, length
c
      logical xinter,check,lposit,lback
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
      write (*,*)
      write (*,*) 'Which task ?'
      write (*,*) ' A = Average (normal)'
      write (*,*) ' E = Expand (normal)'
      write (*,*) ' B = Both average and expand (normal)'
      write (*,*) ' M = Mask-less average and expand (one step)'
      write (*,*) ' Q = Quit right now'
      write (*,*) ' T = Test interpolation'
      write (*,*) ' P = Average; enforce positivity'
      write (*,*) ' R = Average and expand; enforce positivity'
      write (*,*) ' X = Expand; keep zero background'
      write (*,*) ' Y = Average and expand; keep zero background'
      write (*,*) ' Z = Average and expand; zero back; positivity'
      write (*,*)
c
      call textin (' Task ?',task)
      call upcase (task)
      call textut (' Task :',task)
c
      if (task .eq. 'P') then
        task = 'A'
        lposit = .true.
      else if (task .eq. 'R') then
        task = 'B'
        lposit = .true.
      else if (task .eq. 'X') then
        task = 'E'
        lback = .false.
      else if (task .eq. 'Y') then
        task = 'B'
        lback = .false.
      else if (task .eq. 'Z') then
        task = 'B'
        lback = .false.
        lposit = .true.
      end if
c
c ... QUIT
c
c      if (task .eq. 'Q') call gkquit
      if (task .eq. 'Q') return
c
c ... Gerard's hidden interpolation test
c
      if (task .eq. 'T') then
        do i=1,3
          orgna (i)  = -10
          grida (i)  = 100
          exta (i)   =  21
          cell (i)   = 100.0
          cell (i+3) =  90.0
        end do
        call intest (mapa, exta(1), exta(2), exta(3), orgna, 
     +    cell, grida, ierr)
c        call gkquit
        return
      end if
c
      if (task .ne. 'A' .and.
     +    task .ne. 'E' .and.
     +    task .ne. 'M' .and.
     +    task .ne. 'B') then
       call errcon ('Invalid task')
       return
      end if
c
      if (task .eq. 'A') then
        call prompt (' ==> Average density')
      else if (task .eq. 'E') then
        call prompt (' ==> Expand density')
      else if (task .eq. 'B') then
        call prompt (' ==> Average and expand density')
      else if (task .eq. 'M') then
        call prompt (' ==> Average and expand without a mask')
      end if
c
      if (task .eq. 'A' .or. task .eq. 'B') then
        if (lposit) then
          call prompt (' ==> Enforce positivity')
        else
          call prompt (' ==> No positivity constraint')
        end if
      end if
c
      if (task .eq. 'E' .or. task .eq. 'B') then
        if (lback) then
          call prompt (' ==> Set background level')
        else
          call prompt (' ==> Keep zero background')
        end if
      end if
c
c ... get map
c
      write (*,*)
      if (task .eq. 'A' .or. task .eq. 'B' .or.
     +    task .eq. 'M') then
        call textin (' Map to be averaged  ?',fmapa)
        call textut (' Map to be averaged  :',fmapa)
      else
        call textin (
     +    ' Example map of your ASU or unit cell ?',fmapa)
        call textut (
     +    ' Example map of your ASU or unit cell :',fmapa)
      end if
c
      if (length(fmapa) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
      if (task .eq. 'A' .or. task .eq. 'B' .or.
     +    task .eq. 'M') then
        call edin (fmapa, 1, mapa, orgna, exta, grida, uvwa,
     $             cella, spgrp, maxsiz)
      else
        call edhdr (fmapa, 1, orgna, exta, grida, uvwa,
     $             cella, spgrp, mapa, maxsiz)
      end if
      close (1)
      call telmap (grida,orgna,exta,cella)
c
c ... get mask
c
      if (task .ne. 'M') then
        write (*,*)
        call textin (' Mask file ?',fmask)
        call textut (' Input mask :',fmask)
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
     $    (3, maskb, orgnb, extb, gridb, cellb, maxmsk, ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading mask')
          return
        end if
        close (3)
        call telmap (gridb,orgnb,extb,cellb)
      end if
c
c ... check
c
      if (task .eq. 'E') then
c
        write (*,*)
        call textin (' Previously averaged map ?',fmapc)
        call textut (' Previously averaged map :',fmapc)
c
        if (length(fmapc) .lt. 1) then
          call errcon ('No file name provided')
          return
        end if
c
        call edin (fmapc, 1, mapb, orgnc, extc, gridc, uvwa,
     $             cellc, spgrp, maxsiz)
        close (1)
        call telmap (gridc,orgnc,extc,cellc)
c
        check = .true.
        do i=1,3
          check = (check .and. (orgnc(i) .eq. orgnb(i)))
          check = (check .and. (gridc(i) .eq. gridb(i)))
          check = (check .and. (extc(i)  .eq. extb(i)))
          check = (check .and. (abs(cellc(i)-cellb(i)) .le. 0.01))
          check = (check .and.
     +             (abs(cellc(i+3)-cellb(i+3)) .le. 0.01))
        end do
        if (.not. check) then
          call errcon ('Map should have same grid etc. as mask')
          return
        end if
      end if
c
      if (task .ne. 'M') then
        check = .true.
        do i=1,3
          check = (check .and. (grida(i) .eq. gridb(i)))
          check = (check .and. (abs(cella(i)-cellb(i)) .le. 0.01))
          check = (check .and.
     +            (abs(cella(i+3)-cellb(i+3)) .le. 0.01))
        end do
        if (.not. check) then
          call errcon ('Map should have same grid etc. as mask')
          return
        end if
      end if
c
      write (*,*)
      if (task.eq.'A') then
        call textin (' Output averaged CCP4 map ?',fmapb)
        call textut (' Output averaged CCP4 map :',fmapb)
      else
        call textin (' Output expanded CCP4 map ?',fmapb)
        call textut (' Output expanded CCP4 map :',fmapb)
      end if
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
      if (ctrt .lt. 1) then
        call errcon ('No NCS operators')
        return
      end if
c
cc      call anancs (ctrt,rtbtoa,.true.,ierr)
cc      if (ierr .ne. 0)
cc     +  call errstp ('In NCS operators')
c
      if (task .eq. 'B') then
c
        call flusho (6)
        call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $    mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $    cell, grid, rtbtoa, ctrt, rtsym, ctsym, ierr,
     +    avy,avysq,avxy,avxmy,lposit)
        write (*,*)
        call flusho (6)
        call sub02x (mapa, exta(1), exta(2), exta(3), orgna, grida,
     $    mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $    cell, grid, rtbtoa, rtatob, ctrt, rtsym, ctsym,lback)
        write (*,*)
        call flusho (6)
        call edout 
     $    (fmapb, 4, mapa, orgna, exta, grida, uvwa, cella, spgrp)
c
      else if (task .eq. 'A') then
c
        call flusho (6)
        call sub01x (mapa, exta(1), exta(2), exta(3), orgna, 
     $    mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $    cell, grid, rtbtoa, ctrt, rtsym, ctsym, ierr,
     +    avy,avysq,avxy,avxmy,lposit)
        write (*,*)
        call flusho (6)
        call edout 
     $    (fmapb, 4, mapb, orgnb, extb, gridb, uvwa, cellb, spgrp)
c
      else if (task .eq. 'E') then
c
        call flusho (6)
        call sub02x (mapa, exta(1), exta(2), exta(3), orgna, grida,
     $    mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $    cell, grid, rtbtoa, rtatob, ctrt, rtsym, ctsym,lback)
        write (*,*)
        call flusho (6)
        call edout 
     $    (fmapb, 4, mapa, orgna, exta, grida, uvwa, cella, spgrp)
c
      else if (task .eq. 'M') then
c
        call flusho (6)
        call sub03x (mapa, exta(1), exta(2), exta(3), orgna, 
     $    mapb, cell, grid, rtbtoa, ctrt, rtsym, ctsym, ierr,
     +    avy,avysq,avxy,avxmy,lposit)
        write (*,*)
        call flusho (6)
        call edout 
     $    (fmapb, 4, mapb, orgna, exta, grida, uvwa, cella, spgrp)
c
      end if
c
c      call gkquit
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
     +   avy, avysq, avxy, avxmy, lposit)
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
      logical lposit
c
code ...
c
      ierr = -1
      nbit = 0
      call prompt (' Average ...')
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
      call fvalut (' Map  Cell   :',6,cell)
      call fvalut (' Spacing     :',3,grid)
      call ivalut (' Map  origin :',3,orgna)
      call ivalut (' Map  extent :',3,ext)
      call ivalut (' Mask origin :',3,orgnb)
      call ivalut (' Mask extent :',3,extb)
      call fvalut (' FORGN :',3,forgn)
      call fvalut (' FEXT  :',3,fext)
      call fvalut (' GEXT  :',3,gext)
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
      call jvalut (' Points in mask :',1,bobo1)
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
          if (ngjk .eq. bobo3) then
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
              if (nerr2 .lt. 10) then
                call prompt (' Interpolation error !')
                write (*,'(1x,a,1x,3f8.2,1x,3i4)')
     +            'x & ext',x,exta1,exta2,exta3
              else if (nerr2 .eq. 10) then
                call prompt (
     +        ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
c
            end if
120       continue
        else
        end if
100   continue
c
      ierr = 0
      call jvalut (' Calls to BLDBIT      :',1,nbit)
      call jvalut (' Severe FRCSYM errors :',1,nerr1)
      call jvalut (' Interpolation errors :',1,nerr2)
c
c ... print correlation coefficients for the various operators
c
      call jvalut (' Nr of mask points :',1,ngjk)
      f = float (ngjk)
      do loop=1,ctrt
        cc = (avxy(loop)/f - avx*avy(loop)/(f*f))/ 
     +    (sqrt(avxsq/f - (avx/f)**2) *
     +     sqrt(avysq(loop)/f - (avy(loop)/f)**2))
        write (*,'(a,i3,a,f8.5)') ' Corr. coeff. for operator ',
     +    loop,' = ',cc
        cc = avxmy(loop) / avax
        write (*,'(a,i3,a,f8.5)') ' R-factor for operator ',
     +    loop,' w.r.t. operator 1 = ',cc
      end do
c
c --- Normalize by number of additions
c
      if (lposit) then
        call rvalut (' Positivity; set <= 0 to:',1,small)
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
        call prompt (' Averaging without positivity constraint')
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
     $   cell, grid, rtbtoa, rtatob, ctrt, rtsym, ctsym, lback)
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
     $     value, x(3), x1(3), x2(3), xdum, factor
c
      logical lback
c
      equivalence (ijk1(1), i1), (ijk1(2), j1), (ijk1(3), k1)
      equivalence (ijk2(1), i2), (ijk2(2), j2), (ijk2(3), k2)
c
code ...
c
      call prompt (' Expand ...')
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
      call fvalut (' Nr of asymm units in map   :',1,xdum)
      call ivalut (' Nr of asymm units in map   :',1,nasu)
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
      do 330 k=1,exta3
      do 330 j=1,exta2
      do 330 i=1,exta1
330      mapa(i,j,k) = 0.
c
      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
      call jvalut (' Points in mask :',1,bobo1)
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
          if (ngk1 .eq. bobo3) then
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
c --- Now fix up the background so that the average density is 0.0
c
      if (.not. lback) then
        call prompt (' Keeping background at zero')
      else
        call prompt (' Calculating background ...')
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
        k = exta1*exta2*exta3
        call jvalut (' Nr of points in map  :',1,k)
        call jvalut (' Nr of masked points  :',1,ctin)
        call jvalut (' Nr of solvent points :',1,(k-ctin))
        factor = 100.0 * float(k-ctin) / float(k)
        call fvalut (' Solvent fraction (%) :',1,factor)
c
        call rvalut (' Average density inside masks :',1,avrho)
        call rvalut (' Average density in solvent   :',1,back)
c
        call prompt (' Setting background ...')
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
        call rvalut (' Average density overall      :',1,avrho)
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
c
c
c
      subroutine intest (mapa, exta, extb, extc, origin, 
     +                   cell, grid, ierr)
c
c ... interpolation test
c
      implicit none
c
      integer maxpnt
      parameter (maxpnt=100000)
c
      integer exta,extb,extc
      real mapa (exta,extb,extc)
c
      integer origin(3),grid(3),ext(3)
      integer ierr,ifunc,iseed,npnt,i,j,k,jlist,ngrid
c
      real values (maxpnt,4),fralim(2,3),a(3,3),b(3,3)
      real cell(6),ggg(3),x(3),y(3),x2(3)
      real f1a,f1b,xdum,rmsd,shap,corr,rf1,rf2
      real avea,sdva,xmin,xmax,xtot,x6,x7,x8,x9
c
code ...
c
      ext (1) = exta
      ext (2) = extb
      ext (3) = extc
      ngrid = exta * extb * extc
c
      write (*,*)
      call prompt (' === Interpolation Test Routine ===')
      write (*,*)
      call jvalut (' Max nr of random test points :',1,maxpnt)
      write (*,*)
c
      iseed = -1
      call ivalin (' Seed for random-number generator ?',1,iseed)
      call ivalut (' Seed for random-number generator :',1,iseed)
      call gkrand (xdum,0.0,0.0,iseed)
      write (*,*)
c
      call ivalut (' Map origin :',3,origin)
      call ivalut (' Map extent :',3,ext)
      call ivalut (' Map grid   :',3,grid)
c
      ifunc = 1
      npnt = maxpnt
      jlist = 10000
c
      f1a = 10.0
      f1b = 500
c
c ... major test loop
c
   10 continue
c
      write (*,'(99(1x,a/))')
     +  ' ',
     +  'Select function to test (or "0" to QUIT):',
     +  ' 1 ==> R = A * EXP (-B*(x^2+y^2+z^2))',
     +  ' 2 ==> R = A * (COS (x) + COS (y) + COS (z))'
c
      call ivalin (' Function ?',1,ifunc)
      call ivalut (' Function :',1,ifunc)
      if (ifunc .le. 0) return
      if (ifunc .gt. 2) then
        call errcon ('Invalid function selected')
        return
      end if
c
      write (*,*)
      call jvalin (' Nr of points to test ?',1,npnt)
      npnt = max (10, min (npnt,maxpnt))
      call jvalut (' Nr of test points :',1,npnt)
c
      write (*,*)
      call jvalin (' List every Nth point; N ?',1,jlist)
      jlist = max (0, min (npnt,jlist))
      call jvalut (' N :',1,jlist)
c
      write (*,*)
      call fvalin (' Map cell   ?',3,cell)
      call fvalut (' Map cell   :',6,cell)
c
      do i=1,3
        ggg(i) = cell(i)/float(grid(i))
        fralim (1,i) = (origin(i)+3)*ggg(i)/cell(i)
        fralim (2,i) = (origin(i)+ext(i)-3)*ggg(i)/cell(i)
      end do
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell, a, 0)
      call orthog (cell, b, 1)
c
      write (*,*)
      call fvalut (' Spacings :',3,ggg)
      call fvalut (' A limits (fractional) :',2,fralim(1,1))
      call fvalut (' B limits (fractional) :',2,fralim(1,2))
      call fvalut (' C limits (fractional) :',2,fralim(1,3))
c
      if (ifunc .eq. 1) then
c
        write (*,*)
        call rvalin (' Value for constant A ?',1,f1a)
        call rvalut (' Value for constant A :',1,f1a)
        call rvalin (' Value for constant B ?',1,f1b)
        call rvalut (' Value for constant B :',1,f1b)
c
        call prompt (' Calculating "map" ...')
        call flusho (6)
c
        do i=1,ext(1)
          do j=1,ext(2)
            do k=1,ext(3)
c
              x(1) = (i-1+ origin(1))*ggg(1)/cell(1)
              x(2) = (j-1+ origin(2))*ggg(2)/cell(2)
              x(3) = (k-1+ origin(3))*ggg(3)/cell(3)
c
              call mulmtx (a, x, x2, 3, 3, 1)
              xdum = x2(1)*x2(1) + x2(2)*x2(2) + x2(3)*x2(3)
              mapa (i,j,k) = f1a * exp(-f1b*xdum)
c
            end do
          end do
        end do
c
      else if (ifunc .eq. 2) then
c
        write (*,*)
        call rvalin (' Value for constant A ?',1,f1a)
        call rvalut (' Value for constant A :',1,f1a)
c
        call prompt (' Calculating "map" ...')
        call flusho (6)
c
        do i=1,ext(1)
          do j=1,ext(2)
            do k=1,ext(3)
c
              x(1) = (i-1+ origin(1))*ggg(1)/cell(1)
              x(2) = (j-1+ origin(2))*ggg(2)/cell(2)
              x(3) = (k-1+ origin(3))*ggg(3)/cell(3)
c
              call mulmtx (a, x, x2, 3, 3, 1)
              xdum = cos(x2(1)) + cos(x2(2)) + cos(x2(3))
              mapa (i,j,k) = f1a * xdum
c
            end do
          end do
        end do
c
      end if
c
      call xstats (mapa,ngrid,avea,sdva,xmin,xmax,xtot)
      write (*,*)
      call jvalut (' Nr of points in map:',1,ngrid)
      write (*,*)
      call prompt (' Testing interpolation routines ...')
      call flusho (6)
c
c ... testing interpolation
c
      j = 0
c
      do i=1,npnt
c
c ... generate random fractional coordinates
c
        call gkrand (x(1),fralim(1,1),fralim(2,1),0)
        call gkrand (x(2),fralim(1,2),fralim(2,2),0)
        call gkrand (x(3),fralim(1,3),fralim(2,3),0)
c
c ... convert to Angstroms
c
        call mulmtx (a, x, x2, 3, 3, 1)
c
c ... evaluate actual function value
c
        if (ifunc .eq. 1) then
          xdum = x2(1)*x2(1) + x2(2)*x2(2) + x2(3)*x2(3)
          values (i,1) = f1a * exp(-f1b*xdum)
        else if (ifunc .eq. 2) then
          xdum = cos(x2(1)) + cos(x2(2)) + cos(x2(3))
          values (i,1) = f1a * xdum
        end if
c
c ... map fractional coordinates to indices in map array
c
        y(1) = x(1)*cell(1)/ggg(1) - float(origin(1)) + 1
        y(2) = x(2)*cell(2)/ggg(2) - float(origin(2)) + 1
        y(3) = x(3)*cell(3)/ggg(3) - float(origin(3)) + 1
c
c ... try the three interpolation methods
c
        call i n t r p l (mapa,exta,extb,extc,y,values(i,2),ierr)
        call linint      (mapa,exta,extb,extc,y,values(i,3),ierr)
        call nn8         (mapa,exta,extb,extc,y,values(i,4),ierr)
c
        j = j + 1
c
        if (j .eq. jlist .or. ierr .ne. 0) then
          if (ierr .ne. 0) call errcon ('Interpolation')
          write (*,6000) i,x,y,(values(i,k),k=1,4)
          if (j .eq. jlist) j = 0
        end if
c
      end do
c
      call flusho (6)
      write (*,*)
      write (*,6110) 'Type  ','     Minimum',
     +  '     Maximum','     Sum','     Average',
     +  '     St.Dev.'
c
      write (*,6200) '"MAP" ',xmin,xmax,xtot,avea,sdva
c
      call xstats (values(1,1),npnt,avea,sdva,xmin,xmax,xtot)
      write (*,6200) 'EXACT ',xmin,xmax,xtot,avea,sdva
c
      call xstats (values(1,1),npnt,avea,sdva,xmin,xmax,xtot)
      write (*,6200) 'INTRPL',xmin,xmax,xtot,avea,sdva
c
      call xstats (values(1,1),npnt,avea,sdva,xmin,xmax,xtot)
      write (*,6200) 'LININT',xmin,xmax,xtot,avea,sdva
c
      call xstats (values(1,1),npnt,avea,sdva,xmin,xmax,xtot)
      write (*,6200) 'NN8   ',xmin,xmax,xtot,avea,sdva
c
      write (*,*)
      write (*,6110) 'Type  ','    RMSD    ',
     +  ' Simil.Index',' Corr.Coeff.','  R-factor/1',
     +  '  R-factor/2'
c
      call xystat (values(1,1),values(1,2),npnt,
     +             rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
      write (*,6100) 'INTRPL',rmsd,shap,corr,rf1,rf2
c
      call xystat (values(1,1),values(1,3),npnt,
     +             rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
      write (*,6100) 'LININT',rmsd,shap,corr,rf1,rf2
c
      call xystat (values(1,1),values(1,4),npnt,
     +             rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
      write (*,6100) 'NN8   ',rmsd,shap,corr,rf1,rf2
c
      write (*,*)
      call flusho (6)
c
      goto 10
c
 6000 format (' # ',i8,' X = ',3f8.4,' I = ',3f8.2/
     +  12x,'V = ',1p,4e12.4)
 6110 format (' ',a6,5a12)
 6100 format (' ',a6,1pe12.4,0p,4f12.5)
 6200 format (' ',a6,1p,5e12.4)
c
      end
c
c
c
      subroutine sub03x 
     $  (mapa, exta1, exta2, exta3, orgna, 
     $   mapb, cell, grid, rtbtoa, ctrt, rtsym, ctsym, ierr,
     +   avy, avysq, avxy, avxmy, lposit)
c
c ---  do mask-less averaging
c
c ... NOTE: as it is now, it only works for proper symmetry !
c
c     need a separate option for mask-less and improper symm
c
c     generate ALL operators A->B, B->A, A->C, B->C etc.
c     use set with most similar density values ???
c
c     e.g., if NCS=2, then try A->B and B->A; use
c     the value most similar to the NULL operator ???
c
c     if more molecules, try all possibilities:
c       suppose it's in mol A, then A-B, A-C, etc
c       suppose it's in mol B, then B-A, B-C, etc
c     afterwards, select the one with smallest st.devn. ???
c     (==> can assign molecule, namely A or B => output
c          mask if A ???)
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
      real mapb(exta1, exta2, exta3)
      real cell(6), grid(3), rtbtoa(12,*), rtsym(12,*)
      integer ctrt, ctsym,ierr,nbit,ngjk,bobo1,bobo2,bobo3
c
      real avx,avxsq,avy(*),avysq(*),avxy(*),avxmy(*)
      real avax,val1,f,cc
c
      integer errcod, i, j, k, l, loop, m, ext(3)
      integer nerr1,nerr2
c
      real a(3,3), b(3,3), forgn(3), fext(3), mapbit(4,4,4),
     $     value, x(3), x1(3), x2(3),gext(3),xdum
c
      logical lposit
c
code ...
c
      ierr = -1
      nbit = 0
      call prompt (' Average & Expand ...')
c
c ---      A = frac to ang, B=Ang to frac
c
      call orthog (cell, a, 0)
      call orthog (cell, b, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
c
c ---      Generate map edges in fractional cell coords
c
      do 130 i=1,3
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
130     gext(i)  = float(ext(i)+orgna(i)-2)*grid(i)/cell(i)
c
      call fvalut (' Map  Cell  :',6,cell)
      call fvalut (' Spacing    :',3,grid)
      call ivalut (' Map origin :',3,orgna)
      call ivalut (' Map extent :',3,ext)
      call fvalut (' FORGN :',3,forgn)
      call fvalut (' FEXT  :',3,fext)
      call fvalut (' GEXT  :',3,gext)
c
c ---      Loop over map
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
      bobo1 = exta1*exta2*exta3
      call jvalut (' Points in map :',1,bobo1)
      bobo2 = bobo1/20
      bobo3 = bobo2
c
      do 100 k=1,exta3
      do 100 j=1,exta2
      do 100 i=1,exta1
        mapb(i,j,k) = 0.
          ngjk = ngjk + 1
c
          if (ngjk .eq. bobo3) then
            xdum = 100.0 * float(ngjk-1)/float(bobo1)
            call fvalut (' Progress (% map) :',1,xdum)
            bobo3 = bobo3 + bobo2
          end if
c
          x(1) = (i-1+ orgna(1))*grid(1)/cell(1)
          x(2) = (j-1+ orgna(2))*grid(2)/cell(2)
          x(3) = (k-1+ orgna(3))*grid(3)/cell(3)
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
                  call fvalut (' Map point :',3,x2)
                  call fvalut (' NCS point :',3,x1)
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
              avy (loop)   = avy(loop) + value
              avysq (loop) = avysq(loop) + value*value
              avxy (loop)  = avxy(loop) + value*val1
              avxmy (loop) = avxmy(loop) + abs(value-val1)
c
            else
              if (nerr2 .lt. 10) then
                call prompt (' Interpolation error !')
                write (*,'(1x,a,1x,3f8.2,1x,3i4)')
     +            'x & ext',x,exta1,exta2,exta3
              else if (nerr2 .eq. 10) then
                call prompt (
     +        ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
c
            end if
120       continue
100   continue
c
      ierr = 0
      call jvalut (' Calls to BLDBIT      :',1,nbit)
      call jvalut (' Severe FRCSYM errors :',1,nerr1)
      call jvalut (' Interpolation errors :',1,nerr2)
c
c ... print correlation coefficients for the various operators
c
      call jvalut (' Nr of map points :',1,ngjk)
      f = float (ngjk)
      do loop=1,ctrt
        cc = (avxy(loop)/f- avx*avy(loop)/(f*f))/ 
     +    (sqrt(avxsq/f- (avx/f)**2) *
     +     sqrt(avysq(loop)/f- (avy(loop)/f)**2))
        write (*,'(a,i3,a,f8.5)') ' Corr. coeff. for operator ',
     +    loop,' = ',cc
        cc = avxmy(loop) / avax
        write (*,'(a,i3,a,f8.5)') ' R-factor for operator ',
     +    loop,' w.r.t. operator 1 = ',cc
      end do
c
c --- Normalize by number of additions
c
        if (ctrt .gt. 1) then
          xdum = 1.0 / float(ctrt)
          do 200 k=1,exta3
          do 200 j=1,exta2
          do 200 i=1,exta1
200         mapb(i,j,k) = mapb(i,j,k) * xdum
        end if
c
      return
      end
