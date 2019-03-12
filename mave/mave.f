      program mave
c
c ... do averaging in different spacegroups
c
c     930317 - Gerard Kleywegt - first attempt; average works;
c                                expand does NOT; based on
c                                Alwyn's standard averaging code
c     930318 - complete rewrite; seems to work okay now
c     930329 - minor changes
c     930330 - implementation of RT_ortho_improve
c     930331 - minor changes
c     951106 - implement EZ skewing
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'MAVE', vers = '080901/4.1.3')
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
      call domave (%val(iaptr),%val(ibptr),%val(icptr),
     +             mapsize, masksize)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
c
   99 continue
      call gkquit ()
c
      end
c
c
c
      subroutine domave (mapa, mapb, maskb, maxsiz, maxmsk)
c
c ... do averaging in different spacegroups
c
c     930317 - Gerard Kleywegt - first attempt; average works;
c                                expand does NOT; based on
c                                Alwyn's standard averaging code
c     930318 - complete rewrite; seems to work okay now
c     930329 - minor changes
c     930330 - implementation of RT_ortho_improve
c     930331 - minor changes
c     951106 - implement EZ skewing
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz,maxsym,maxncs,maxmsk
c      parameter (maxsiz = maxgk1)
      parameter (maxsym = maxgk2)
      parameter (maxncs = maxgk2)
c
      integer orgna(3),exta(3),grida(3),uvwa(3)
      real mapa(maxsiz),cella(6),spaca(3)
c
      integer orgnb(3),extb(3),gridb(3),maskb(maxmsk)
      real mapb(maxsiz),cellb(6),spacb(3)
c
      integer orgnr(3),extr(3),gridr(3),uvwr(3),spgrpr
      real cellr(6)
c
      real artatob(12,maxncs),artbtoa(12,maxncs)
      real rtsyma(12,maxsym),rtsymb(12,maxsym),or2or(12)
      real avy(maxncs),avysq(maxncs),avxy(maxncs),avxmy(maxncs)
      real f2ca(3,3),c2fa(3,3),f2cb(3,3),c2fb(3,3)
      real x(3),x1(3),x2(3),x3(3),x4(3)
c
      character file*80,fmt*80,newmap*80,par*25,partyp*1
      character task*1,fmapa*80,fmapb*80,fmask*80,ncspar*25
c
      integer ctrta,ctrtb,ctsyma,ctsymb,errcod,i,j,spgrpa,ierr
      integer length,k,l,gmask,ndum,margin
c
      real total,user,sys
c
      logical xinter,check,lposit,lback
c
code ...
c
c      call gainit (prognm,vers)
c
c ... 950118 - check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
c      call jvalut (' Max size of maps and masks    :',1,maxsiz)
      call jvalut (' Max nr of spacegroup symm ops :',1,maxsym)
      call jvalut (' Max nr of NCS symm ops        :',1,maxncs)
c
      task = 'A'
      ctrta = 0
      ctrtb = 0
      ctsyma = 0
      ctsymb = 0
      fmapa = ' '
      fmapb = ' '
      fmask = ' '
c
      lposit = .false.
      lback  = .true.
c
      write (*,*)
      write (*,*) 'Which task ?'
      write (*,*) ' A = Average map on mask in reference crystal'
      write (*,*) ' P = average & enforce positivity'
      write (*,*) ' E = Expand map from mask in reference crystal'
      write (*,*) ' F = expand & keep zero background'
      write (*,*) ' I = Improve RT-operator from reference to target'
      write (*,*) ' S = ez Skewing'
      write (*,*) ' Q = Quit right now'
c
      write (*,*)
      call textin (' Task ?',task)
      call upcase (task)
      call textut (' Task :',task)
c
      if (task .eq. 'P') then
        task = 'A'
        lposit = .true.
      else if (task .eq. 'F') then
        task = 'E'
        lback = .false.
      end if
c
c ... QUIT
c
      if (task .eq. 'Q') return
c
      if (task .ne. 'A' .and.
     +    task .ne. 'E' .and.
     +    task .ne. 'S' .and.
     +    task .ne. 'I') then
        call errcon ('Invalid task')
        return
      end if
c
      if (task .eq. 'A') then
        call prompt (' ==> Average density')
      else if (task .eq. 'E') then
        call prompt (' ==> Expand density')
      else if (task .eq. 'S') then
        call prompt (' ==> EZ Skewing')
      end if
c
      if (task .eq. 'A') then
        if (lposit) then
          call prompt (' ==> Enforce positivity')
        else
          call prompt (' ==> No positivity constraint')
        end if
      end if
c
      if (task .eq. 'E') then
        if (lback) then
          call prompt (' ==> Set background level')
        else
          call prompt (' ==> Keep zero background')
        end if
      end if
c
      if (task .ne. 'I' .and. task .ne. 'S') then
        write (*,*)
        write (*,*) 'IMPORTANT NOTE ...'
        write (*,*)
     +    'The REFERENCE and TARGET crystals may be identical !'
      end if
c
c ... get map A
c
      write (*,*)
      if (task .eq. 'A') then
        call textin (' Map to be averaged of TARGET crystal ?',fmapa)
        call textut (' Map to be averaged of TARGET crystal :',fmapa)
      else if (task .eq. 'E') then
        call textin (
     +    ' Averaged map in mask of REFERENCE crystal ?',fmapa)
        call textut (
     +    ' Averaged map in mask of REFERENCE crystal :',fmapa)
      else if (task .eq. 'S') then
        call textin (
     +    ' Input map of REFERENCE crystal ?',fmapa)
        call textut (
     +    ' Input map of REFERENCE crystal :',fmapa)
      else
        call textin (
     +    ' Map of TARGET crystal ?',fmapa)
        call textut (
     +    ' Map of TARGET crystal :',fmapa)
      end if
c
      if (length(fmapa) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
      call edin (fmapa, 1, mapa,  orgna, exta, grida, uvwa,
     +  cella, spgrpa, maxsiz)
      close (1)
      call telmap (grida,orgna,exta,cella)
c
c ... get mask
c
      write (*,*)
      if (task .eq. 'S') then
        write (*,*) 'Mask is *OPTIONAL* for skewing !'
        write (*,*) 'If you do not provide a mask name the entire'
        write (*,*) 'map will be skewed.'
      end if
      call textin (' Mask in REFERENCE crystal ?',fmask)
      call textut (' Input mask :',fmask)
c
      if (task .eq. 'S' .and. length(fmask) .lt. 1) then
        write (*,*) 'Skew entire map; set entire mask to "1"'
        margin = 3
        call ivalut (' Except for a border margin :',1,margin)
c
        i=exta(1)*exta(2)*exta(3)
        if (i .gt. maxmsk) then
          call errcon ('Requested mask too big!')
          call jvalut (' Requested mask size :',1,i)
          call jvalut (' Available mask size :',1,maxmsk)
          call prompt (' Rerun program with more memory for mask')
          call errstp ('Not enough memory for mask')
        end if
        call fillem (maskb,exta(1),exta(2),exta(3),margin)
c
c        do i=1,maxmsk
c          maskb (i) = 1
c        end do
c
        do i=1,3
          orgnb(i) = orgna(i)
          gridb(i) = grida(i)
          extb(i) = exta(i)
          cellb(i) = cella(i)
          cellb(i+3) = cella(i+3)
        end do
      else
        if (length(fmask) .lt. 1) then
          call errcon ('No file name provided')
          return
        end if
        call xopxoa (3,fmask,xinter(),errcod)
        if (errcod .ne. 0) then
          call errcon ('While opening mask')
          return
        end if
        call maskin 
     $    (3, maskb, orgnb, extb, gridb, cellb, maxmsk, errcod)
        if (errcod .ne. 0) then
          call errcon ('While reading mask')
          return
        end if
        close (3)
        call telmap (gridb,orgnb,extb,cellb)
      end if
c
c ... check
c
      if (task .eq. 'E' .or. task .eq. 'S') then
        check = .true.
        do i=1,3
          if (task .eq. 'E')
     +      check = (check .and. (orgna(i) .eq. orgnb(i)))
          check = (check .and. (grida(i) .eq. gridb(i)))
          if (task .eq. 'E')
     +      check = (check .and. (exta(i)  .eq. extb(i)))
          check = (check .and. (abs(cella(i)-cellb(i)) .le. 0.01))
          check = (check .and.
     +             (abs(cella(i+3)-cellb(i+3)) .le. 0.01))
        end do
        if (.not. check) then
          call errcon (
     +      'Map should have same grid etc. as mask')
          return
        end if
      end if
c
c ... bring mask to same grid as input map
c
      if (task .eq. 'S' .and. length(fmask) .ge. 1) then
        check = .true.
        do i=1,3
          if (orgna(i) .ne. orgnb(i)) check = .false.
          if (exta(i) .ne. extb(i))   check = .false.
        end do
        if (check) then
          call prompt (
     +      ' Mask has same origin and extent as input map')
        else
          call prompt (
     +      ' Bring mask to same origin and extent as input map')
          call fixmsk (mapb,maskb,exta(1),exta(2),exta(3),
     +      extb(1),extb(2),extb(3),orgna,orgnb)
        end if
        do i=1,3
          orgnb(i) = orgna(i)
          gridb(i) = grida(i)
          extb(i) = exta(i)
          cellb(i) = cella(i)
          cellb(i+3) = cella(i+3)
        end do
      end if
c
c ... get RT operator relating the two molecules from O datablock
c
      file = ' '
      write (*,*)
      call textin (
     +  ' File with RT FROM reference TO target crystal ?',file)
      call textut (
     +  ' File with RT FROM reference TO target crystal :',file)
      call opoodb (1,file,par,partyp,j,fmt,errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening file')
        return
      end if
      if (j .ne. 12) then
        call errcon ('Not an RT operator')
        return
      end if
      read (1,fmt,err=9999,end=9999) (or2or(i),i=1,j)
      close (1)
      if (task .eq. 'I') ncspar = par
c
c ... sometimes we need SYMMOPS of reference crystal
c
      if (task .eq. 'E' .or. task .eq. 'I' .or. task .eq. 'S') then
        file = ' '
        write (*,*)
        call textin (
     +    ' File with symmetry operators for REFERENCE crystal ?',file)
        call textut (
     +    ' File with symmetry operators for REFERENCE crystal :',file)
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
          call errcon ('While opening file')
          return
        end if
        ctsymb = j/12
        call jvalut (' Nr of symmetry operators :',1,ctsymb)
        if (ctsymb .lt. 1 .or. ctsymb .gt. maxsym) then
          call errcon ('Invalid number of symmetry operators')
          return
        end if
        if ((12*ctsymb) .ne. j) then
          call errcon ('Invalid nr of elements; must be N*12')
        return
      end if
        read (1,fmt,err=9999,end=9999)
     +    ((rtsymb(i, j),i=1,12),j=1,ctsymb)
        do j=1,ctsymb
          call fratra (rtsymb(10,j))
        end do
        close (1)
c
        call prompt ('0Reference crystal symm-ops:')
        call anasgs (ctsymb,rtsymb,.true.,ierr)
        if (ierr .ne. 0) then
          call errcon ('In spacegroup symmetry operators')
          return
        end if
c
      end if
c
      if (task .eq. 'S') then
        ctsyma = 1
        do i=1,12
          rtsyma(i,1) = 0.0
        end do
        rtsyma (1,1) = 1.0
        rtsyma (5,1) = 1.0
        rtsyma (9,1) = 1.0
      else
        file = ' '
        write (*,*)
        call textin (
     +    ' File with symmetry operators for TARGET crystal ?',file)
        call textut (
     +    ' File with symmetry operators for TARGET crystal :',file)
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
          call errcon ('While opening file')
          return
        end if
        ctsyma = j/12
        call jvalut (' Nr of symmetry operators :',1,ctsyma)
        if (ctsyma .lt. 1 .or. ctsyma .gt. maxsym) then
          call errcon ('Invalid number of symmetry operators')
          return
        end if
        if ((12*ctsyma) .ne. j) then
          call errcon ('Invalid nr of elements; must be N*12')
          return
        end if
        read (1,fmt,err=9999,end=9999)
     +    ((rtsyma(i, j),i=1,12),j=1,ctsyma)
        do j=1,ctsyma
          call fratra (rtsyma(10,j))
        end do
        close (1)
      end if
c
      call prompt ('0Target crystal symm-ops:')
      call anasgs (ctsyma,rtsyma,.true.,ierr)
      if (ierr .ne. 0) then
        call errcon ('In spacegroup symmetry operators')
        return
      end if
c
c ... handle IMPROVE option here
c
      if (task .eq. 'I') then
        file = ' '
        write (*,*)
        call textin (
     +    ' Map in REFERENCE crystal ?',file)
        call textut (
     +    ' Map in REFERENCE crystal :',file)
        call edin (file,1,mapb,orgnr,extr,gridr,uvwr,cellr,spgrpr,
     +             maxsiz)
        call telmap (gridr,orgnr,extr,cellr)
        check = .true.
        do i=1,3
          check = (check .and. (gridr(i) .eq. gridb(i)))
          check = (check .and. (abs(cellr(i)-cellb(i)) .le. 0.01))
          check = (check .and.
     +             (abs(cellr(i+3)-cellb(i+3)) .le. 0.01))
        end do
        if (.not. check) then
          call errcon ('Map should have same grid and cell as mask')
          return
        end if
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call flusho (6)
c
        call mimp (mapa,exta(1),exta(2),exta(3),grida,cella,orgna,
     +       maskb,extb(1),extb(2),extb(3),gridb,cellb,orgnb,
     +       mapb,extr(1),extr(2),extr(3),orgnr,
     +       or2or,rtsyma,ctsyma,rtsymb,ctsymb,ncspar)
c
        goto 9000
c
      end if
c
c ... handle SKEW option here
c
      if (task .eq. 'S') then
        newmap = ' '
        write (*,*)
        call textin (
     +    ' Output map in TARGET crystal ?',newmap)
        call textut (
     +    ' Output map in TARGET crystal :',newmap)
c
        do i=1,3
          uvwr  (i) = i
          cellr (i) = 100.0
          cellr (i+3) = 90.0
        end do
c
        call fvalin (' Cell parameters of TARGET crystal ?',6,cellr)
        call fvalut (' Cell parameters of TARGET crystal :',6,cellr)
c
        do i=1,3
          if (cellr(i) .le. 1.0 .or. cellr(i+3) .le. 20.0 .or.
     +        cellr(i+3) .ge. 160.0) then
            call errcon ('Invalid cell parameters')
            return
           end if
c
          spaca(i) = cella(i)/float(grida(i))
          gridr(i) = nint(cellr(i)/spaca(i))
        end do
c
        spgrpr = 1
        call jvalin (' Spacegroup NUMBER ?',1,spgrpr)
        call jvalut (' Spacegroup NUMBER :',1,spgrpr)
        if (spgrpr .lt. 1 .or. spgrpr .gt. 230) then
          call errcon ('Invalid spacegroup number')
          return
        end if
c
        call jvalin (' Grid of TARGET crystal ?',3,gridr)
        call jvalut (' Grid of TARGET crystal :',3,gridr)
        do i=1,3
          spacb(i) = cellr(i)/float(gridr(i))
          if (spacb(i) .le. 0.01 .or. spacb(i) .ge. 10.0) then
            call errcon ('Invalid grid spacing')
            return
          end if
          orgnr (i) =  999999
          extr  (i) = -999999
        end do
c
        call orthog (cellr, f2ca, 0)
        call orthog (cellr, c2fa, 1)
        call orthog (cellb, f2cb, 0)
        call orthog (cellb, c2fb, 1)
c
        ndum = 0
        do i=1,exta(1)
          x(1) = float(orgna(1) + i - 1) / 
     +           float(grida(1))
          do j=1,exta(2)
            x(2) = float(orgna(2) + j - 1) / 
     +             float(grida(2))
            do k=1,exta(3)
              l = gmask (i,j,k,
     +                   extb(1),extb(2),extb(3),maskb)
c
ccc              l = gmask (i+orgna(1)-orgnb(1),
ccc     +                   j+orgna(2)-orgnb(2),
ccc     +                   k+orgna(3)-orgnb(3),
ccc     +                   extb(1),extb(2),extb(3),maskb)
c
              if (l .eq. 1) then
                ndum = ndum + 1
                x(3) = float(orgna(3) + k - 1) / 
     +                 float(grida(3))
c
c                write (*,*)
c                call fvalut (' Frac R :',3,x)
                call mulmtx (f2cb,x,x1,3,3,1)
c                call fvalut (' Cart R :',3,x1)
                call vecrtv (x1,x2,1,or2or(1),or2or(10))
c                call fvalut (' Cart T :',3,x2)
                call mulmtx (c2fa,x2,x3,3,3,1)
c                call fvalut (' Frac T :',3,x3)
c
                do l=1,3
                  x4(l)=x3(l)*gridr(l)
                  if ( ndum .eq. 1) then
                    orgnr(l) = nint(x4(l))
                    extr (l) = nint(x4(l))
                  else
                    orgnr(l) = min(orgnr(l),nint(x4(l)))
                    extr(l)  = max(extr(l), nint(x4(l)))
                  end if
                end do
c                call fvalut (' Indx T :',3,x4)
              end if
c
            end do
          end do
        end do
c
        call jvalut (' Min mask indices in new map :',3,orgnr)
        call jvalut (' Max mask indices in new map :',3,extr)
c
        do i=1,3
          j = extr(i) - orgnr(i) + 1
          orgnr (i) = orgnr (i) - 5
          extr  (i) = j + 10
        end do
c
        call jvalut (' Origin of new map :',3,orgnr)
        call jvalut (' Extent of new map :',3,extr)
c
c ... check if map not too big
c
        check = .true.
 2413   continue
        k = extr(1)*extr(2)*extr(3)
        call jvalut (' Size of new map :',1,k)
        if (k .gt. maxsiz) then
          if (check) call errcon ('Map too big - reducing')
          do i=1,3
            orgnr (i) = orgnr (i) + 2
            extr  (i) = extr  (i) - 4
          end do
          check = .false.
          goto 2413
        end if
        if (.not. check) then
          call jvalut (' Origin of new map :',3,orgnr)
          call jvalut (' Extent of new map :',3,extr)
        end if
c
c ... no NCS
c
        ctrta = 1
        do i=1,12
          artatob (i,1) = 0
        end do
        artatob (1,1) = 1.0
        artatob (5,1) = 1.0
        artatob (9,1) = 1.0
        do i=1,12
          artbtoa (i,1) = artatob (i,1)
        end do
c
c ... keep zero background
c
        lback = .false.
c
c ... do the skewing
c
        call sub02x (
     +    mapb,extr(1),extr(2),extr(3),orgnr,spacb,cellr,gridr,
     +    mapa,exta(1),exta(2),exta(3),orgna,spaca,cella,
     +    maskb, artatob, artbtoa, ctrta, rtsyma, ctsyma,
     +    rtsymb, ctsymb, or2or, lback, ierr, 'S')
c
c ... write map
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call flusho (6)
c
        if (ierr .ne. 0) then
          call errcon ('While skewing')
          return
        end if
c
        if (length(newmap) .ge. 1) then
          call edout (newmap,4,mapb,orgnr,extr,gridr,uvwr,
     +                cellr,spgrpr)
          call prompt (' Map written')
        else
          call prompt (' No map written')
        end if
c
        goto 9000
c
      end if
c
  210 continue
      file = ' '
      write (*,*)
      call textin (
     +  ' File with NCS operator(s) for TARGET crystal ?',file)
      call textut (
     +  ' File with NCS operator(s) for TARGET crystal :',file)
      if (length(file) .lt. 1) goto 200
c
      call rdoncs (1,file,ctrta,maxncs,artatob,errcod)
      if (errcod .ne. 0) then
        call errcon ('While reading NCS operators')
        return
      end if
c
cc      call opoodb (1,file,par,partyp,j,fmt,errcod)
cc      if (errcod .ne. 0) call errstp ('While opening file')
cc      if (j .ne. 12) call errstp ('Not an RT operator')
cc      ctrta = ctrta + 1
cc      if (ctrta .gt. maxncs) call errstp (
cc     +  'Too many NCS operators')
cc      read (1,fmt,err=9999,end=9999) (artatob(i, ctrta),i=1,12)
cc      close (1)
c
      goto 210
c
  200 continue
      if (ctrta .lt. 1) then
        call errcon ('No NCS operators')
        return
      end if
c
cc      call anancs (ctrta,artatob,.true.,ierr)
cc      if (ierr .ne. 0) call errstp ('In NCS operators')
c
c ... task = AVERAGE
c
      if (task .eq. 'A') then
c
c ... need UVW and SPGRP of reference crystal
c
  220   continue
        write (*,*)
        write (*,*) 'I need parameters from the averaged map'
        write (*,*) 'around the mask in the REFERENCE crystal'
        file = ' '
        call textin (
     +    ' Averaged map in REFERENCE crystal ?',file)
        call textut (
     +    ' Averaged map in REFERENCE crystal :',file)
        call edhdr (file,1,orgnr,extr,gridr,uvwr,cellr,spgrpr,
     +              mapb,maxsiz)
        call telmap (gridr,orgnr,extr,cellr)
        check = .true.
        do i=1,3
cc          check = (check .and. (orgnr(i) .eq. orgnb(i)))
          check = (check .and. (gridr(i) .eq. gridb(i)))
cc          check = (check .and. (extr(i)  .eq. extb(i)))
          check = (check .and. (abs(cellr(i)-cellb(i)) .le. 0.01))
          check = (check .and.
     +             (abs(cellr(i+3)-cellb(i+3)) .le. 0.01))
        end do
        if (.not. check) then
          call errcon ('Map should have same grid etc. as mask')
          goto 220
        end if
        call prompt (' Thanks !')
c
        newmap = ' '
        write (*,*)
        call textin (
     +    ' Averaged map on mask in REFERENCE crystal ?',newmap)
        call textut (
     +    ' Averaged map on mask in REFERENCE crystal :',newmap)
c
        do i=1,3
          spaca(i) = cella(i)/float(grida(i))
          spacb(i) = cellb(i)/float(gridb(i))
        end do
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call flusho (6)
c
        call prompt (' Start averaging ...')
c
        call sub01x (
     +    mapa,exta(1),exta(2),exta(3),orgna,spaca,cella,
     +    mapb,extb(1),extb(2),extb(3),orgnb,spacb,cellb,
     +    maskb, artatob, ctrta, rtsyma, ctsyma, or2or, ierr,
     +    avy,avysq,avxy,avxmy,lposit)
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call flusho (6)
c
        if (ierr .ne. 0) then
          call errcon ('While averaging')
          return
        end if
c
        if (length(newmap) .ge. 1) then
          call edout (newmap,4,mapb,orgnb,extb,gridb,uvwr,
     +                cellb,spgrpr)
          call prompt (' Map written')
        else
          call prompt (' No map written')
        end if
c
c ... TASK = EXPAND
c
      else
c
c ... need UVW and SPGRP of reference crystal
c
        write (*,*)
        write (*,*) 'I need parameters from a map in the asymmetric'
        write (*,*) 'unit of the TARGET crystal'
        file = ' '
        call textin (
     +    ' Asymmetric unit map of TARGET crystal ?',file)
        call textut (
     +    ' Asymmetric unit map of TARGET crystal :',file)
        call edhdr (file,1,orgnr,extr,gridr,uvwr,cellr,spgrpr,
     +              mapb,maxsiz)
        call telmap (gridr,orgnr,extr,cellr)
        call prompt (' Thanks !')
c
        newmap = ' '
        write (*,*)
        call textin (
     +    ' Expanded map in TARGET crystal ?',newmap)
        call textut (
     +    ' Expanded map in TARGET crystal :',newmap)
        call textut (' Output map :',newmap)
c
        do i=1,3
          spaca(i) = cella(i)/float(grida(i))
          spacb(i) = cellr(i)/float(gridr(i))
        end do
c
c ... bring map from mask in REFERENCE crystal to the asymm unit
c     of the TARGET crystal
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call flusho (6)
c
        call prompt (' Start expanding ...')
c
        call sub02x (
     +    mapb,extr(1),extr(2),extr(3),orgnr,spacb,cellr,gridr,
     +    mapa,exta(1),exta(2),exta(3),orgna,spaca,cella,
     +    maskb, artatob, artbtoa, ctrta, rtsyma, ctsyma,
     +    rtsymb, ctsymb, or2or, lback, ierr, 'E')
c
        call gkdcpu (total,user,sys)
        write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
        call flusho (6)
c
        if (ierr .ne. 0) then
          call errcon ('While expanding')
          return
        end if
c
        if (length(newmap) .ge. 1) then
          call edout (newmap,4,mapb,orgnr,extr,gridr,uvwr,
     +                cellr,spgrpr)
          call prompt (' Map written')
        else
          call prompt (' No map written')
        end if
c
      end if
c
c ... done
c
 9000 continue
      return
c
 9999 continue
      call errcon ('While reading file')
      return
c
      end
c
c
c
      subroutine sub01x (
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
      integer bobo1,bobo2,bobo3,nerr1,nerr2
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
      call fvalut (' FORGNA :',3,forgna)
      call fvalut (' FEXTA  :',3,fexta)
      call fvalut (' GEXTA  :',3,gexta)
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
      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
      call jvalut (' Points in mask :',1,bobo1)
      bobo2 = bobo1/10
      bobo3 = bobo2
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
        if (ngjk .eq. bobo3) then
          xdum = 100.0 * float(ngjk-1)/float(bobo1)
          call fvalut (' Progress (% mask) :',1,xdum)
          bobo3 = bobo3 + bobo2
        end if
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
     +            ' NOTE: further FRCSYM errors but not listed!!!')
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
      call jvalut (' Calls to BLDBIT      :',1,nbit)
      call jvalut (' Severe FRCSYM errors :',1,nerr1)
      call jvalut (' Interpolation errors :',1,nerr2)
c
c ... print correlation coefficients for the various operators
c
      call jvalut (' Nr of mask points :',1,ngjk)
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
      subroutine sub02x (
     +    mapa,exta1,exta2,exta3,orgna,spaca,cella,grida,
     +    mapb,extb1,extb2,extb3,orgnb,spacb,cellb,
     +    maskb, rtbtoa, rtatob, ctrt, rtsyma, ctsyma,
     +    rtsymb, ctsymb, or2or, lback, ierr, mytask)
c
c ... expand into a (possibly) different spacegroup
c
c ... MAPB and MASKB are on the same grid etc.
c     MAPA has to be filled
c
c ... Gerard Kleywegt, 930317,18
c     based on code by Alwyn Jones, 1991
c
      implicit none
c
      integer orgna(3), exta1, exta2, exta3, grida(3)
      real mapa(exta1, exta2, exta3),spaca(3),cella(6)
c
      integer orgnb(3), extb1, extb2, extb3
      real mapb(extb1, extb2, extb3),spacb(3),cellb(6)
      integer maskb(extb1, extb2, extb3)
c
      real rtbtoa(12,*),rtatob(12,*),rtsyma(12,*),rtsymb(12,*)
      real or2or(12),ro2ro(12)
      integer ctrt,ctsyma,ctsymb,ierr,nbit
c
      integer ctin, errcod, i, ina(3), j, k, l, loop, nasu
      integer i1, j1, k1, ijk1(3), ijk(3)
      integer exta(3), i2, j2, k2, ijk2(3),extb(3)
      integer m,ngk1,ngk2,bobo1,bobo2,bobo3
c
      integer kk1,kk2,kk3,nmask,nout,ntwice,ijk3(3)
      integer nerr1,nerr2,nerr3,nerr4
c
      real avrho, back, value, x1(3), x2(3), mapbit(4,4,4)
      real f2ca(3,3),c2fa(3,3),f2cb(3,3),c2fb(3,3), x(3)
      real fmsk(3),cmsk(3),cmska(3),ncmska(3),fcmska(3),xdum
      real fcmapa(3),ncmapa(3),cmapa(3),cmapb(3),fmapb(3)
      real fbori(3),fbext(3),gbext(3),faori(3),faext(3),gaext(3)
c
      logical lback
c
      character task*1,mytask*1
c
      equivalence (ijk(1),  i),  (ijk(2),  j),  (ijk(3),  k)
      equivalence (ijk1(1), i1), (ijk1(2), j1), (ijk1(3), k1)
      equivalence (ijk2(1), i2), (ijk2(2), j2), (ijk2(3), k2)
c
code ...
c
      ierr = -1
      nbit = 0
c
      task = mytask
      call upcase (task)
      if (task .ne. 'S') task = 'E'
c
      exta(1) = exta1
      exta(2) = exta2
      exta(3) = exta3
c
      extb(1) = extb1
      extb(2) = extb2
      extb(3) = extb3
c
c      call fvalut (' Cell A    :',6,cella)
c      call fvalut (' Cell B    :',6,cellb)
c      call fvalut (' Spacing A :',3,spaca)
c      call fvalut (' Spacing B :',3,spacb)
c      call ivalut (' Grid    A :',3,grida)
c      call ivalut (' Origin  A :',3,orgna)
c      call ivalut (' Origin  B :',3,orgnb)
c      call ivalut (' Extent  A :',3,exta)
c      call ivalut (' Extent  B :',3,extb)
c
c      do i=1,ctrt
c        call fvalut (' NCS :',12,rtbtoa(1,i))
c      end do
c
      call orthog (cella, f2ca, 0)
      call orthog (cella, c2fa, 1)
      call orthog (cellb, f2cb, 0)
      call orthog (cellb, c2fb, 1)
c
c ... generate edges of both maps in fractionals
c
      do i=1,3
        fbori (i) = float(orgnb(i))          *spacb(i)/cellb(i)
        fbext (i) = float(extb(i)+orgnb(i)-1)*spacb(i)/cellb(i)
        gbext (i) = float(extb(i)+orgnb(i)-2)*spacb(i)/cellb(i)
        faori (i) = float(orgna(i))          *spaca(i)/cella(i)
        faext (i) = float(exta(i)+orgna(i)-1)*spaca(i)/cella(i)
        gaext (i) = float(exta(i)+orgna(i)-2)*spaca(i)/cella(i)
      end do
c
c      call fvalut (' FORGNA :',3,faori)
c      call fvalut (' FEXTA  :',3,faext)
c      call fvalut (' GEXTA  :',3,gaext)
c      call fvalut (' FORGNB :',3,fbori)
c      call fvalut (' FEXTB  :',3,fbext)
c      call fvalut (' GEXTB  :',3,gbext)
c
c ... get inverse RT operators
c
      do i=1,ctrt
        do j=1,9
          rtatob (j,i) = rtbtoa (j,i)
        end do
        call matinv (rtatob(1,i), 3, x, x1, x2)
        call mulmtx (rtatob(1,i), rtbtoa(10,i), rtatob(10,i), 3, 3, 1)
        do j=10,12
          rtatob (j,i) = - rtatob(j,i)
        end do
c        call fvalut (' Inverted NCS :',12,rtatob(1,i))
      end do
c
      do j=1,9
        ro2ro(j) = or2or(j)
      end do
      call matinv (ro2ro(1), 3, x, x1, x2)
      call mulmtx (ro2ro(1), or2or(10), ro2ro(10), 3,3,1)
      do j=10,12
        ro2ro(j) = -ro2ro(j)
      end do
      call fvalut (' REF -> TARGET :',12,or2or)
      call fvalut (' TARGET -> REF :',12,ro2ro)
c
c .. find out how many asymmetric units map A comprises
c
      i = grida(1)*grida(2)*grida(3)
      k = i/ctsyma
      j = exta1*exta2*exta3
      xdum = float(j)/float(k)
      nasu = nint (xdum)
      if (abs(xdum-float(nasu)) .gt. 0.2) then
        call jvalut (' Nr of points in unit cell  :',1,i)
        call jvalut (' Nr of symmetry operators   :',1,ctsyma)
        call jvalut (' Nr of points in asymm unit :',1,k)
        call jvalut (' Nr of points in map        :',1,j)
        call fvalut (' Nr of asymm units in map   :',1,xdum)
        call jvalut (' Nearest integer            :',1,nasu)
        call prompt (' WARNING - Not an integer number !!')
c        return
      end if
      call fvalut (' Nr of asymm units in map   :',1,xdum)
      call ivalut (' Nr of asymm units in map   :',1,nasu)
c
      if (xdum .lt. 1.0) then
        call errcon ('Less than one asymmetric unit')
        if (task .eq. 'E') return
      end if
c
      if (task .eq. 'S') then
        call prompt (' Skewing - create only one copy of the map')
        nasu = 1
      end if
c
      ngk1 = 0
      ngk2 = 0
      nerr1 = 0
      nerr2 = 0
      nerr3 = 0
      nerr4 = 0
      nout = 0
      ntwice = 0
c
c ... blank map A
c
      do 330 k=1,exta3
      do 330 j=1,exta2
      do 330 i=1,exta1
        mapa(i,j,k) = 0.0
  330 continue
c
      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
      call jvalut (' Points in mask :',1,bobo1)
      bobo2 = bobo1/10
      bobo3 = bobo2
c
c ... loop over the points in the mask
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
c
        if (maskb(i,j,k) .ne. 1) goto 100
        ngk1 = ngk1 + 1
c
        if (ngk1 .eq. bobo3) then
          xdum = 100.0 * float(ngk1-1)/float(bobo1)
          call fvalut (' Progress (% mask) :',1,xdum)
          bobo3 = bobo3 + bobo2
        end if
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
c                call ivalut (' Mask point in reference  :',3,ijk)
c                call fvalut (' Fractional               :',3,fmsk)
c                call fvalut (' Cartesian                :',3,cmsk)
c                call fvalut (' Cartesian in second xtal :',3,cmska)
c                call fvalut (' Cartesian of NCS mate    :',3,ncmska)
c                call fvalut (' Fractional               :',3,fcmska)
c
          do l=1,3
            x(l) = fcmska(l)*cella(l)/spaca(l)- float(orgna(l))+ 1
            ina(l) = x(l)
          end do
c
c ... set nearest eight points through interpolation
c
          do 150 k1=ina(3),ina(3)+1
          do 150 j1=ina(2),ina(2)+1
          do 150 i1=ina(1),ina(1)+1
c
c ... check if this point has been set previously
c
            call  frctrn (ijk1, ijk2, orgna, grida, faori, faext, 
     +                    rtsyma, ctsyma, errcod)
            if (errcod .ne. 0) then
              if (nerr1 .lt. 10) then
                call errcon ('FRCTRN error')
                call ivalut (' Mask indices :',3,ijk)
                call ivalut (' Map  indices :',3,ijk1)
              else if (nerr1 .eq. 10) then
                call prompt (
     +              ' NOTE: further FRCTRN errors but not listed!!!')
              end if
              nerr1 = nerr1 + 1
              goto 150
            end if
c
            if (mapa(i2,j2,k2) .ne. 0.0) goto 150
c
            do l=1,3
              fcmapa(l) = (ijk1(l)-1+ orgna(l))*spaca(l)/cella(l)
            end do
c
c ... FCMAPA = fract coords of NCS mate in A
c     NCMAPA = Cartesian coords of same
c     CMAPA  = Cartesian coords of reference molecule
c     CMAPB  = Cartesian coords in map B
c     FMAPB  = fract coords in B; should be close to FMSK !!!
c
            call mulmtx (f2ca,fcmapa,ncmapa,3,3,1)
            call vecrtv (ncmapa,cmapa,1,rtatob(1,loop),rtatob(10,loop))
            call vecrtv (cmapa,cmapb,1,ro2ro(1),ro2ro(10))
            call mulmtx (c2fb,cmapb,fmapb,3,3,1)
c
c            call fvalut (' Mask point  :',3,fmsk)
c            call fvalut (' FCMSKA      :',3,fcmska)
c            call fvalut (' FCMAPA      :',3,fcmapa)
c            call fvalut (' Transformed :',3,fmapb)
c
c ... force it into the mask by applying B's symmops
c
            call frcsym (fmapb,fbori,gbext,rtsymb,ctsymb,errcod)
c
c ... if error, then close to edge; build 4x4x4 map
c
            if (errcod .ne. 0) then
c
              nbit = nbit + 1
              call bldbit (mapb,extb1,extb2,extb3,orgnb,spacb,
     +          cellb,fbori,fbext,rtsymb,ctsymb,mapbit,fmapb,errcod)
c
              if (errcod .ne. 0) then
                if (nerr4 .lt. 10) then
                  call errcon ('Serious FRCSYM error')
                  call ivalut (' Mask point in reference  :',3,ijk)
                  call fvalut (' Fractional               :',3,fmsk)
                  call fvalut (' Cartesian                :',3,cmsk)
                  call fvalut (' Cartesian in second xtal :',3,cmska)
                  call fvalut (' Cartesian of NCS mate    :',3,ncmska)
                  call fvalut (' Fractional               :',3,fcmska)
                  call ivalut (' Map point in target      :',3,ijk1)
                  call fvalut (' Fractional               :',3,fcmapa)
                  call fvalut (' Cartesian of NCS mate    :',3,ncmapa)
                  call fvalut (' Cartesian in second xtal :',3,cmapa)
                  call fvalut (' Cartesian                :',3,cmapb)
                  call fvalut (' Fractional               :',3,fmapb)
                else if (nerr4 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr4 = nerr4 + 1
                goto 150
              end if
c
              do l=1,3
                x(l) = fmapb(l)*cellb(l)/spacb(l)- float(orgnb(l))+ 1
                m = x(l)
                x(l) = x(l) - float(m-1)
              end do
c
              call intrpl (mapbit, 4, 4, 4, x, value, errcod)
c
            else
c
c ... if FRCSYM worked, just interpolate
c
              do l=1,3
                x(l) = fmapb(l)*cellb(l)/spacb(l)- float(orgnb(l))+ 1
              end do
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
              call intrpl (mapb, extb1, extb2, extb3, x, value, errcod)
c
            end if
c
            if (errcod .eq. 0) then
              call frcval (ijk1,mapa,value,exta1,exta2,exta3,
     +          orgna,grida,faori,faext,rtsyma,ctsyma,
     +          nasu,ntwice,errcod)
              if (errcod .eq. 0) then
                ngk2 = ngk2 + 1
              else
                if (nerr3 .lt. 10) then
                  call errcon ('Serious FRCVAL error')
                  call ivalut (' Mask point in reference  :',3,ijk)
                  call fvalut (' Fractional               :',3,fmsk)
                  call fvalut (' Cartesian                :',3,cmsk)
                  call fvalut (' Cartesian in second xtal :',3,cmska)
                  call fvalut (' Cartesian of NCS mate    :',3,ncmska)
                  call fvalut (' Fractional               :',3,fcmska)
                  call ivalut (' Map point in target      :',3,ijk1)
                  call fvalut (' Fractional               :',3,fcmapa)
                  call fvalut (' Cartesian of NCS mate    :',3,ncmapa)
                  call fvalut (' Cartesian in second xtal :',3,cmapa)
                  call fvalut (' Cartesian                :',3,cmapb)
                  call fvalut (' Fractional               :',3,fmapb)
                else if (nerr3 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCVAL errors but not listed!!!')
                end if
                nerr3 = nerr3 + 1
              end if
            else
              if (nerr2 .lt. 10) then
                call errcon ('Interpolation error')
                call ivalut (' Mask point in reference  :',3,ijk)
                call fvalut (' Fractional               :',3,fmsk)
                call fvalut (' Cartesian                :',3,cmsk)
                call fvalut (' Cartesian in second xtal :',3,cmska)
                call fvalut (' Cartesian of NCS mate    :',3,ncmska)
                call fvalut (' Fractional               :',3,fcmska)
                call ivalut (' Map point in target      :',3,ijk1)
                call fvalut (' Fractional               :',3,fcmapa)
                call fvalut (' Cartesian of NCS mate    :',3,ncmapa)
                call fvalut (' Cartesian in second xtal :',3,cmapa)
                call fvalut (' Cartesian                :',3,cmapb)
                call fvalut (' Fractional               :',3,fmapb)
              else if (nerr2 .eq. 10) then
                call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
            end if
  150     continue
  120   continue
  100 continue
c
c ... now set background such that average density is zero
c
      call jvalut (' Calls to BLDBIT   :',1,nbit)
      call jvalut (' Points in mask    :',1,ngk1)
      call jvalut (' Set in asymm unit :',1,ngk2)
      call jvalut (' Total points set  :',1,(ngk2*nasu))
c
      call jvalut (' FRCTRN errors        :',1,nerr1)
      call jvalut (' FRCSYM errors        :',1,nerr4)
      call jvalut (' Interpolation errors :',1,nerr2)
      call jvalut (' FRCVAL errors        :',1,nerr3)
c
      call jvalut (' Points outside mask  :',1,nout)
      call jvalut (' Points set > 1 *     :',1,ntwice)
c
      if (.not. lback) then
        call prompt (' Keeping background at zero')
      else
        call prompt (' Calculating background ...')
        avrho = 0.0
        ctin = 0
        do 500 k=1,exta3
        do 500 j=1,exta2
        do 500 i=1,exta1
          if (mapa(i,j,k) .ne. 0.0) then
            avrho = avrho + mapa(i,j,k)
            ctin = ctin+1
          end if
  500   continue
        back  = -avrho/float(exta1*exta2*exta3-ctin)
        avrho =  avrho/float(ctin)
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
  510   continue
c
        avrho = avrho/float(exta1*exta2*exta3)
        call rvalut (' Average density overall      :',1,avrho)
c
      end if
c
      ierr = 0
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
250         i2(j) = nint(x1(j)*float(grid(j)))- orgn(j)+ 1
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
      subroutine mimp (mapa,exta1,exta2,exta3,grida,cella,orgna,
     +           maskb,extb1,extb2,extb3,gridb,cellb,orgnb,
     +           mapb,extr1,extr2,extr3,orgnr,
     +           rtbtoa,rtsym,ctsym,rtsymb,ctsymb,odbnam)
c
c ... MIMP - improve RT-operator
c
c ... Gerard Kleywegt @ 930330
c
      implicit none
c
      integer exta1,exta2,exta3,grida(3),orgna(3)
      real mapa(exta1,exta2,exta3),cella(6)
c
      integer extb1,extb2,extb3,gridb(3),orgnb(3)
      integer maskb(extb1,extb2,extb3)
      integer extr1,extr2,extr3,orgnr(3)
      real mapb(extr1,extr2,extr3),cellb(6)
c
      integer ctsym,ctsymb
      real rtbtoa(12),rtsym(12,ctsym),rtsymb(12,ctsymb)
c
      real agrid(3),bgrid(3),tshift,tend,rshift,rend,cswap,ccmin,shift
c
      integer i,ndens,length,errcod,leng1
c
      logical xinter
c
      character t*1,line*80,odbnam*(*)
c
code ...
c
      do i=1,3
        agrid (i) = cella(i)/float(grida(i))
        bgrid (i) = cellb(i)/float(gridb(i))
      end do
c
      t = 'A'
      tshift = 0.5
      tend = 0.02
      rshift = 2.0
      rend = 0.1
      cswap = 0.002
      ccmin = 0.00001
      ndens = 4
      shift = 0.5
c
c ... main event loop
c
c
200   continue
      write (*,*)
      call prompt (' Select one of the options Q(uit),')
      call prompt (' T(ranslation), R(otation) or A(uto)')
      call textin (' Search type ?',t)
      call textut (' Search type :',t)
      call upcase (t)
c
      if (t .eq. 'Q') then
        line = ' '
        call textin (' File for new operator ?',line)
        call textut (' File for new operator :',line)
        if (length(line) .lt. 1) return
        call xopxua (1,line,xinter(),errcod)
        if (errcod .ne. 0) then
          call errcon ('While opening O datablock file')
          return
        end if
        write (line,*) odbnam(1:leng1(odbnam)),' R 12 (3f15.6)'
        write (1,'(a)') line(1:leng1(line))
        write (1,'(3f15.6)') rtbtoa
        return
      end if
c
      if (t .ne. 'R' .and. t .ne. 'A' .and. t .ne. 'T') goto 200
c
      if (t .eq. 'A') then
c
        call fvalin (' Initial translation step size ?',1,tshift)
        call fvalut (' Initial translation step size :',1,tshift)
        if (tshift .le. 0.)  tshift = 0.5
        call fvalin (' Convergence when step size    ?',1,tend)
        call fvalut (' Convergence when step size    :',1,tend)
        if (tend .le. 0.)  tend = 0.02
        call fvalin (' Initial rotation step size    ?',1,rshift)
        call fvalut (' Initial rotation step size    :',1,rshift)
        if (rshift .le. 0.)  rshift = 2.0
        call fvalin (' Convergence when step size    ?',1,rend)
        call fvalut (' Convergence when step size    :',1,rend)
        if (rend .le. 0.)  rend = 0.1
        call rvalin (' Change for swapping R <-> T   ?',1,cswap)
        call rvalut (' Change for swapping R <-> T   :',1,cswap)
        if (ccmin .le. 0.)  cswap = 0.002
        call rvalin (' Change for cc convergence     ?',1,ccmin)
        call rvalut (' Change for cc convergence     :',1,ccmin)
        if (ccmin .le. 0.)  ccmin = 0.0001
c
      else
c
        if (t .eq. 'T') then
          call fvalin (' Step size ?',1,shift)
          call fvalut (' Step size :',1,shift)
          if (shift .le. 0.)  shift = 0.5
        else
          call fvalin (' Step size ?',1,shift)
          call fvalut (' Step size :',1,shift)
          if (shift .le. 0.)  shift = 0.5
        end if
c
      end if
c
c .. ask for sample density
c
      call prompt (' Enter the sample density; 1 means: use all')
      call prompt (' points, 3 means: use every third point etc.')
      call prompt (' Higher value -> faster (but less accurate)')
      call ivalin (' Sample density ?',1,ndens)
      call ivalut (' Sample density :',1,ndens)
      if (ndens .le. 0) then
        ndens = 4
      else if (ndens .gt. 100) then
        call ivalin (' Value too high; re-enter :',1,ndens)
        call ivalut (' Sample density :',1,ndens)
        if (ndens .le. 0 .or. ndens .gt. 100) ndens = 4
      end if
c
c ---	Now go do everything. This is to make subscripting easier
c
      call sub03x (
     +  mapa,exta1,exta2,exta3,orgna,cella,grida,agrid,
     +  maskb,extb1,extb2,extb3,orgnb,cellb,gridb,bgrid,
     +  mapb,extr1,extr2,extr3,orgnr,
     +  rtbtoa,rtsym,ctsym,rtsymb,ctsymb,
     +  t, shift,ndens,tshift,tend,rshift,rend,ccmin,cswap)
c
      goto 200
c
      end
c
c
c
      subroutine sub03x (
     +  mapa,exta1,exta2,exta3,orgna,cella,grida,agrid,
     +  maskb,extb1,extb2,extb3,orgnb,cellb,gridb,bgrid,
     +  mapb,extr1,extr2,extr3,orgnr,
     +  or2or,rtsyma,ctsyma,rtsymb,ctsymb,
     +  orit,shift,ndens,tshift,tend,rshift,rend,ccmin,cswap)
c
      implicit none
c
c ---	Map A data structures
      integer orgna(3), exta1, exta2, exta3, grida(3)
      real mapa(exta1, exta2, exta3),cella(6),agrid(3)
c
c ---	Map B data structure
      integer orgnb(3), extb1, extb2, extb3, gridb(3)
      integer maskb(extb1, extb2, extb3)
      integer extr1,extr2,extr3,orgnr(3)
      real mapb(extr1, extr2, extr3),bgrid(3),cellb(6)
c
      integer ctsyma,ctsymb
      real or2or(12),rtsyma(12,ctsyma),rtsymb(12,ctsymb)
c
      integer ndens
      character orit*1,t*1
      real shift
c
      integer ctmask,errcod,i,i1,i2,i3,j,k,l,ncnt,ncyc
      integer i1max,i2max,i3max,m,nerr1,nerr2
c
      real x(3), x1(3), x2(3),mapbit(4,4,4),cswap,old
      real avx,avxsq,avxy,avy,avysq,ccoef,f,peak,peakrt(12),rt(12)
      real del(3),q(9),qq(9),xmo(3),xp(3),xo(3),fgk(27),ogk(3),oextr
      real sgk(3),fextr,tshift,tend,rshift,rend,ccmin
      real f2cb(3,3),c2fb(3,3),f2ca(3,3),c2fa(3,3),ro2ro(12)
      integer exta(3),extb(3)
      real aforgn(3),bforgn(3),afext(3),bfext(3),agext(3),bgext(3)
      real avalue,bvalue,total,user,sys
c
      integer sumx, sumy, sumz,ngk
c
      logical init,tdone,rdone,first
c
      data init /.false./, first /.true./
c
      save init,first
c
code ...
c
c --- get conversion matrices fractional <-> cartesian
c
      call orthog (cella, f2ca, 0)
      call orthog (cella, c2fa, 1)
      call orthog (cellb, f2cb, 0)
      call orthog (cellb, c2fb, 1)
c
      exta(1) = exta1
      exta(2) = exta2
      exta(3) = exta3
      extb(1) = extr1
      extb(2) = extr2
      extb(3) = extr3
c
c ---	Generate envelope edges in fractional cell coords
      do i=1,3
        aforgn(i) = float(orgna(i))          *agrid(i)/cella(i)
        afext(i)  = float(exta(i)+orgna(i)-1)*agrid(i)/cella(i)
        agext(i)  = float(exta(i)+orgna(i)-2)*agrid(i)/cella(i)
        bforgn(i) = float(orgnr(i))          *bgrid(i)/cellb(i)
        bfext(i)  = float(extb(i)+orgnr(i)-1)*bgrid(i)/cellb(i)
        bgext(i)  = float(extb(i)+orgnr(i)-2)*bgrid(i)/cellb(i)
      end do
c
      if (first) then
        first = .false.
        write (*,*)
        call fvalut (' Cell A    :',6,cella)
        call ivalut (' Origin A  :',3,orgna)
        call ivalut (' Extent A  :',3,exta)
        call ivalut (' Grid A    :',3,grida)
        call fvalut (' Spacing A :',3,agrid)
        call fvalut (' AFORGN    :',3,aforgn)
        call fvalut (' AFEXT     :',3,afext)
        call fvalut (' AGEXT     :',3,agext)
        call fvalut (' Cell B    :',6,cellb)
        call ivalut (' Origin B  :',3,orgnr)
        call ivalut (' Extent B  :',3,extb)
        call ivalut (' Grid B    :',3,gridb)
        call fvalut (' Spacing B :',3,bgrid)
        call fvalut (' BFORGN    :',3,bforgn)
        call fvalut (' BFEXT     :',3,bfext)
        call fvalut (' BGEXT     :',3,bgext)
        write (*,*)
      end if
c
      if (orit .eq. 'A') then
        t = 'T'
        tdone = .false.
        rdone = .false.
        ncyc = 0
      else
        t = orit
      end if
c
      old = -2.0
c
      write (*,15) or2or
c
c ... main loop in case of automatic improvement
c
 6996 continue
c
      if (orit .eq. 'A') then
        ncyc = ncyc + 1
        write (*,6903) ncyc,tshift,tend,rshift,rend
c
 6903 format (/' Automatic Improvement ... Cycle ',i3/
     +  ' Translation step, target ',2f12.6/
     +  ' Rotation    step, target ',2f12.6)
c
        if (t .eq. 'T') then
          shift = tshift
          write (*,*) 'This cycle: Translation refinement'
        else
          shift = rshift
          write (*,*)  'This cycle: Rotation refinement'
        end if
      end if
c
      do 510 j=1,9
510     ro2ro(j) = or2or(j)
      call matinv (ro2ro, 3, x, x1, x2)
      call mulmtx (ro2ro, or2or(10), ro2ro(10), 3, 3, 1)
      do 520 j=10,12
520     xo(j-9) = -ro2ro(j)
      call fvalut (' Rotation origin :',3,xo)
c
      ctmask = 0
      sumx = 0
      sumy = 0
      sumz = 0
c
      do 600 k=1,extb3
      do 600 j=1,extb2
      do 600 i=1,extb1
        if (maskb(i,j,k) .eq. 1) then
          sumx = sumx+ i
          sumy = sumy+ j
          sumz = sumz+ k
          ctmask = ctmask+1
        end if
600   continue
      x(1) = (sumx/ctmask- 1+ orgnb(1))*bgrid(1)/cellb(1)
      x(2) = (sumy/ctmask- 1+ orgnb(2))*bgrid(2)/cellb(2)
      x(3) = (sumz/ctmask- 1+ orgnb(3))*bgrid(3)/cellb(3)
      call mulmtx (f2cb, x, xmo, 3, 3, 1)
      call fvalut (' Centre of mask in REFERENCE :',3,xmo)
      call mulmtx (or2or, xmo, xp, 3, 3, 1)
      do 610 i=1,3
610     xp(i) = xp(i) + or2or(i+9)
      call fvalut (' Centre of mask in TARGET    :',3,xp)
      call jvalut (' Nr of points in mask        :',1,ctmask)
c
      peak = -999.
      i1max = 0
      i2max = 0
      i3max = 0
c
      ngk = 0
      call flusho (6)
c
      do 300 i1=1,3
      do 300 i2=1,3
      do 300 i3=1,3
        do 310 i=1,12
310       rt(i) = or2or(i)
        del(1) = (i1-2)*shift
        del(2) = (i2-2)*shift
        del(3) = (i3-2)*shift
c
        if (t .eq. 'T') then
          do 370 j=10,12
370         rt(j) = rt(j)+ del(j-9)
        else
          do 350 i=1,3
            call matrot (i, del(i), q)
            call mulmtx (q, rt, qq, 3, 3, 3)
            do 360 j=1,9
360           rt(j) = qq(j)
350       continue
          call mulmtx (rt, xmo, rt(10), 3, 3, 1)
          do 380 j=10,12
380         rt(j) = -rt(j)+ xp(j-9)
        end if
c
c ---	Set sums to zero
c
        avx = 0.
        avy = 0.
        avxy = 0.
        avxsq = 0.
        avysq = 0.
        ctmask = 0
        ncnt = 0
        nerr1 = 0
        nerr2 = 0
c
c ---	Loop over mask looking for something
c
        do 100 k=1,extb3
        do 100 j=1,extb2
        do 100 i=1,extb1
c
          if (maskb(i,j,k) .eq. 1) then
c
c ... check if this point needs to be used
c
            if (ndens .gt. 1) then
              ncnt = ncnt + 1
              if (ncnt .lt. ndens) goto 100
              ncnt = 0
            end if
c
            x1(1) = (i-1+ orgnb(1))*bgrid(1)/cellb(1)
            x1(2) = (j-1+ orgnb(2))*bgrid(2)/cellb(2)
            x1(3) = (k-1+ orgnb(3))*bgrid(3)/cellb(3)
c
c ... X1(3) = fractional coordinates in reference map B;
c             get electron density value (BVALUE)
c
            call frcsym (x1, bforgn, bgext, rtsymb, ctsymb, errcod)
            if (errcod .ne. 0) then
	      call bldbit (mapb, extr1, extr2, extr3, orgnr, bgrid,
     $	        cellb, bforgn, bfext, rtsymb, ctsymb, mapbit, x1, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call prompt (' FRCSYM error REF B')
                  call fvalut (' Coordinates :',3,x1)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 100
              end if
	      do 225 l=1,3
	        x(l) = x1(l)*cellb(l)/bgrid(l)- float(orgnr(l))+ 1
	        m = x(l)
	        x(l) = x(l)- float(m-1)
225           continue
	      call intrpl (mapbit, 4, 4, 4, x, bvalue, errcod)
            else
              do 210 l=1,3
210             x(l) = x1(l)*cellb(l)/bgrid(l)- float(orgnr(l))+ 1
c ---	      x is now set of indices to the mapb array
              call intrpl 
     $          (mapb, extr1, extr2, extr3, x, bvalue, errcod)
            end if
            if (errcod .ne. 0) then
              if (nerr2 .lt. 10) then
                call errcon ('Interpolation error REF B')
                call fvalut (' Coordinates 1 :',3,x1)
                call fvalut (' Coordinates 2 :',3,x)
              else if (nerr2 .eq. 10) then
                call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
              goto 100
            end if
c
c ... now convert X1 to fractional coordinates in target map A
c
c ... X1 = fractional map B
c     X2 = cartesian map B
c     X  = cartesian map A
c     X1 = fractional map A
c
            call mulmtx (f2cb, x1, x2, 3, 3, 1)
            call vecrtv (x2,x,1,rt(1),rt(10))
            call mulmtx (c2fa,x,x1,3,3,1)
c
c ---	      Is it in the envelope, if not make it (if possible)
c
            call frcsym (x1, aforgn, agext, rtsyma, ctsyma, errcod)
            if (errcod .ne. 0) then
	      call bldbit (mapa, exta1, exta2, exta3, orgna, agrid,
     $	        cella, aforgn, afext, rtsyma, ctsyma, mapbit, x1, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call prompt (' FRCSYM error TAR A')
                  call fvalut (' Coordinates :',3,x1)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 100
              end if
	      do 125 l=1,3
	        x(l) = x1(l)*cella(l)/agrid(l)- float(orgna(l))+ 1
	        m = x(l)
	        x(l) = x(l)- float(m-1)
125           continue
	      call intrpl (mapbit, 4, 4, 4, x, avalue, errcod)
            else
              do 110 l=1,3
110             x(l) = x1(l)*cella(l)/agrid(l)- float(orgna(l))+ 1
c ---	      x is now set of indices to the mapa array
              call intrpl 
     $          (mapa, exta1, exta2, exta3, x, avalue, errcod)
            end if
            if (errcod .eq. 0) then
              avx = avx+ bvalue
              avy = avy+ avalue
              avxsq = avxsq+ bvalue**2
              avysq = avysq+ avalue**2
              avxy = avxy+ bvalue*avalue
              ctmask = ctmask+ 1
	      else
              if (nerr2 .lt. 10) then
                call prompt (' Interpolation error TAR A')
              else if (nerr2 .eq. 10) then
                call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
            end if
          else
          end if
100     continue
c
        f = float (ctmask)
        ccoef = (avxy/f- avx*avy/(f*f))/ 
     $    (sqrt(avxsq/f- (avx/f)**2)* sqrt(avysq/f- (avy/f)**2))
        write (*,10) del, ccoef
c
        ngk = ngk + 1
        fgk (ngk) = ccoef
c
        if (ccoef .gt. peak) then
          do 330 l=1,12
330         peakrt(l) = rt(l)
          peak = ccoef
          i1max = i1
          i2max = i2
          i3max = i3
	end if
        call flusho (6)
300   continue
c
      write (*,*) 'Nr of mask points checked :',ctmask
      if (nerr1 .gt. 0) then
        write (*,*) 'ERROR - Nr of FRCSYM errors :',nerr1
      end if
      if (nerr2 .gt. 0) then
        write (*,*) 'ERROR - Nr of interpolation errors :',nerr2
      end if
c
      call gkdcpu (total,user,sys)
      write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
      call flusho (6)
c
c ... if max not at 0,0,0 accept it as is
c
      if (i1max.ne.2 .or. i2max.ne.2 .or. i3max.ne.2) then
        do 1400 i=1,12
 1400     or2or(i) = peakrt(i)
        write (*,20) peakrt, peak
        if (orit .ne. 'A') return
        if (abs(old-peak) .le. cswap) then
          if (t .eq. 'T' .and. (.not. rdone)) then
            t = 'R'
          else if (t .eq. 'R' .and. (.not. tdone)) then
            t = 'T'
          end if
        end if
        old = peak
        goto 6996
      end if
c
c  ... else do interpolation
c
      oextr = fgk(14)
      do i=1,3
        ogk (i) = 0.0
      end do
c
      call svdmx3 (fgk,ogk,oextr,sgk,fextr,init)
      write (*,69) 'Central value ',oextr,' at ',ogk
      write (*,69) 'Interpolated  ',fextr,' at ',(shift*sgk(i),i=1,3)
   69 format (1x,a,f10.5,a,3f10.5)
c
      if (fextr .lt. oextr) then
        write (*,*) 'Oops - corr coeff goes up; abort search'
        write (*,20) or2or,oextr
        return
      end if
c
      do 1310 i=1,12
1310    rt(i) = or2or(i)
      del(1) = sgk(1)*shift
      del(2) = sgk(2)*shift
      del(3) = sgk(3)*shift
      if (t .eq. 'T') then
        do 1370 j=10,12
1370      rt(j) = rt(j)+ del(j-9)
      else
        do 1350 i=1,3
          call matrot (i, del(i), q)
          call mulmtx (q, rt, qq, 3, 3, 3)
          do 1360 j=1,9
1360        rt(j) = qq(j)
1350    continue
        call mulmtx (rt, xmo, rt(10), 3, 3, 1)
        do 1380 j=10,12
1380      rt(j) = -rt(j)+ xp(j-9)
      end if
c
      do 400 i=1,12
400     or2or(i) = rt(i)
c
      write (*,25) or2or,fextr
c
c ... decide what to do for the automatic option
c
      old = fextr
c
      if (orit .ne. 'A') return
c
c ... check if corr coeff has changed sufficiently
c
      if ( (fextr-oextr) .lt. ccmin) then
        write (*,*) 'Correlation coefficient converged'
        return
      end if
c
      if (t .eq. 'T') then
c
c ... set the new translation step size to twice the maximum
c      value of the shifts determined in this cycle
c
        tshift = 2.0 * max (abs(del(1)),abs(del(2)),abs(del(3)))
        if (tshift .lt. tend) then
          tdone = .true.
          write (*,*) 'Translation refinement converged'
        end if
        if (tdone .and. rdone) return
        if (.not. rdone) t = 'R'
      else
c
c ... similarly, for the rotation
c
        rshift = 2.0 * max (abs(del(1)),abs(del(2)),abs(del(3)))
        if (rshift .lt. rend) then
          rdone = .true.
          write (*,*) 'Rotation refinement converged'
        end if
        if (tdone .and. rdone) return
        if (.not. tdone) t = 'T'
      end if
c
      goto 6996
c
10    format (' Shift=',3f12.6, ' | Corr coeff=',f12.6)
15    format (' Start Rotation Matrix'/3(1x,3f12.6/),
     +  ' Start Translation'/1x,3f12.6)
20    format (' Best Rotation Matrix'/3(1x,3f12.6/),
     +  ' Best Translation'/1x,3f12.6/
     +  ' Correlation Coefficient = ',f12.6)
25    format (' Best Rotation Matrix'/3(1x,3f12.6/),
     +  ' Best Translation'/1x,3f12.6/
     +  ' Predicted Correlation Coefficient = ',f12.6)
c
      end
c
c
c
      integer function gmask (i,j,k,ext1,ext2,ext3,mask)
c
      implicit none
c
      integer i,j,k,ext1,ext2,ext3
      integer mask(ext1,ext2,ext3)
c
code ...
c
      gmask = - 1
      if (i .le. 0 .or. i .gt. ext1) return
      if (j .le. 0 .or. j .gt. ext2) return
      if (k .le. 0 .or. k .gt. ext3) return
      gmask = mask (i,j,k)
c
      return
      end
c
c
c
      subroutine fixmsk (map,mask,exta1,exta2,exta3,
     +      extb1,extb2,extb3,orgna,orgnb)
c
      implicit NONE
c
      integer exta1,exta2,exta3,extb1,extb2,extb3
c
      real map(*)
c
      integer mask(*)
c
      integer orgna(3),orgnb(3)
c
code ...
c
      call zeromp (map,exta1*exta2*exta3)
c
      call msk2mp (map,mask,exta1,exta2,exta3,
     +      extb1,extb2,extb3,orgna,orgnb)
c
      call mp2msk (map,mask,exta1,exta2,exta3)
c
      return
      end
c
c
c
      subroutine zeromp (map,nn)
c
      implicit NONE
c
      integer nn
c
      real map(nn)
c
      integer i
c
code ...
c
      do i=1,nn
        map (i) = 0.0
      end do
c
      return
      end
c
c
c
      subroutine msk2mp (map,mask,exta1,exta2,exta3,
     +      extb1,extb2,extb3,orgna,orgnb)
c
      implicit NONE
c
      integer exta1,exta2,exta3,extb1,extb2,extb3
c
      real map(exta1,exta2,exta3)
c
      integer mask(extb1,extb2,extb3)
c
      integer orgna(3),orgnb(3)
      integer i,j,k,ii,jj,kk
c
code ...
c
      do i=1,extb1
        ii = i + orgnb(1) - orgna(1)
        if (ii .lt. 1 .or. ii .gt. exta1) goto 1000
        do j=1,extb2
          jj = j + orgnb(2) - orgna(2)
          if (jj .lt. 1 .or. jj .gt. exta2) goto 2000
          do k=1,extb3
            kk = k + orgnb(3) - orgna(3)
            if (kk .lt. 1 .or. kk .gt. exta3) goto 3000
            if (mask(i,j,k) .eq. 1) map(ii,jj,kk) = 1.0
 3000       continue
          end do
 2000     continue
        end do
 1000   continue
      end do
c
      return
      end
c
c
c
      subroutine mp2msk (map,mask,exta1,exta2,exta3)
c
      implicit NONE
c
      integer exta1,exta2,exta3
c
      real map(exta1,exta2,exta3)
c
      integer mask(exta1,exta2,exta3)
c
      integer i,j,k
c
code ...
c
      do i=1,exta1
        do j=1,exta2
          do k=1,exta3
            if (map(i,j,k) .gt. 0.5) then
              mask (i,j,k) = 1
            else
              mask (i,j,k) = 0
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
      subroutine fillem (mask,nx,ny,nz,margin)
c
      implicit none
c
      integer nx,ny,nz,margin
      integer mask(nx,ny,nz)
      integer i,j,k,my,nn
c
code ...
c
      my = max(0,margin)
      do i=1,nx
        do j=1,ny
          do k=1,nz
            mask(i,j,k) = 0
          end do
        end do
      end do
c
      nn = 0
      do i=my+1,nx-my
        do j=my+1,ny-my
          do k=my+1,nz-my
            mask(i,j,k) = 1
            nn = nn + 1
          end do
        end do
      end do
c
      call jvalut (' Mask points total :',1,(nx*ny*nz))
      call jvalut (' Mask points set=1 :',1,nn)
      call jvalut (' Mask points set=0 :',1,(nx*ny*nz - nn))
c
      return
      end
