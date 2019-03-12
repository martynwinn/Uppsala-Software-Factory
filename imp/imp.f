      program imp
c
c ---	Improve the transformation operator
c ---	Assumes both maps are on same grid.
c ---	Alwyn Jones, 2-Oct-91
c
c ... Altered by Gerard Kleywegt @ 930219
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'IMP', vers = '080626/3.0.2')
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
c      integer maskb(1)
c      integer malloc
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
      call doimpr (%val(iaptr), %val(ibptr), %val(icptr),
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
      subroutine doimpr (mapa, mapb, maskb, maxsiz, maxmsk)
c
c ---	Improve the transformation operator
c ---	Assumes both maps are on same grid.
c ---	Alwyn Jones, 2-Oct-91
c
c ... Altered by Gerard Kleywegt @ 930219
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz,maxsym,maxmsk
c      parameter (maxsiz = 9 * maxgk1 / 10)
      parameter (maxsym=maxgk2)
c
c ---	Map A data structures
      integer orgna(3), exta(3), grida(3), uvwa(3)
      real mapa(maxsiz), cella(6)
c
c ---	Map B data structure
      integer orgnb(3), extb(3), gridb(3), maskb(maxmsk)
      real mapb(maxsiz), cellb(6)
c
      real cell(6), grid(3), rtbtoa(12), rtsym(12,maxsym), shift
      real tshift,tend,rshift,rend,ccmin,cswap,total,user,sys
      character file*80, fmt*80, line*80, par*25, partyp*1
      character t*1,ityp*1,mapnam*80,f927*1
      integer ctrt,ctsym,errcod,i,j,spgrp,ndens,length,niter,ierr
      integer leng1
      logical xinter,check
c
code ...
c
c      call gainit (prognm,vers)
c
c      call jvalut (' Max size of map & mask :',1,maxsiz)
      call jvalut (' Max nr of symmetry ops :',1,maxsym)
c
      write (*,*)
      ctrt = 0
      ndens = 1
      f927 = 'C'
c
      file = ' '
      call textin (' Map file ?',file)
      call textut (' Map file :',file)
      if (length(file) .lt. 1) then
        call errcon ('No filename provided')
        return
      end if
c
c ---	read in map a
c
      call edin   (file, 1, mapa,  orgna, exta, grida, uvwa,
     $             cella, spgrp, maxsiz)
      close (1)
      call telmap (grida,orgna,exta,cella)
      mapnam = file
c
c ---	the mask comes from map b
c
      write (*,*)
      file = ' '
      call textin (' Mask file ?',file)
      call textut (' Mask file :',file)
      if (length(file) .lt. 1) then
        call errcon ('No filename provided')
        return
      end if
c
      call xopxoa (3,file,xinter(),errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening mask')
        return
      end if
      call maskin 
     $  (3, maskb, orgnb, extb, gridb, cellb, maxmsk, errcod)
      if (errcod .ne. 0) then
        call errcon ('While reading mask')
        return
      end if
      close (3)
      call telmap (gridb,orgnb,extb,cellb)
c
c ... check cell and grid
c
      check = .true.
      do i=1,3
        check = (check .and. (grida(i) .eq. gridb(i)))
        check = (check .and. (abs(cella(i)-cellb(i)) .le. 0.01))
        check = (check .and.
     +           (abs(cella(i+3)-cellb(i+3)) .le. 0.01))
      end do
      if (.not. check) then
        call errcon (
     +    'Map should have same grid and cell as mask')
        return
      end if
c
      do 100 i=1,3
100     grid(i) = cella(i)/float(grida(i))
c
      do 110 i=1,6
110     cell(i) = cella(i)
c
      file = ' '
      write (*,*)
      call textin (
     +  ' File with symmetry operators (O-style) ?',file)
      call textut (
     +  ' File with symmetry operators (O-style) :',file)
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
        call errcon ('Too many space group operators')
        return
      end if
      if ((12*ctsym) .ne. j) then
        call errcon ('Invalid nr of elements; must be N*12')
        return
      end if
      do j=1,ctsym
        read (1, fmt,err=190,end=190) (rtsym(i,j),i=1,12)
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
      goto 210
c
 190  continue
      call errcon ('While reading O datablock file')
      return
c
c ---	Get rt a->b
c
210   continue
      file = ' '
      write (*,*)
      call textin (
     +  ' File with NCS operator (O style) ?',file)
      call textut (
     +  ' File with NCS operator (O style) :',file)
      if (length(file) .lt. 1) then
        call errcon ('No filename supplied')
        return
      end if
c
      call opoodb (1,file,par,partyp,j,fmt,errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening O datablock file')
        return
      end if
c
      j = min (12,j)
      ctrt = j/12
      if (ctrt .ne. 1) then
	call errcon ('Not an NCS operator')
        return
      end if
      read (1, fmt,err=190,end=190) (rtbtoa(i),i=1,j)
      close (1)
c
      call anancs (1,rtbtoa,.true.,ierr)
      if (ierr .ne. 0) then
        call errcon ('In NCS operators')
        return
      end if
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
      ityp = 'P'
      niter = 0
c
200   continue
      niter = niter + 1
      write (*,*)
      call gkdcpu (total,user,sys)
      write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      write (*,'(1x,a)') ' ',
     +  'Select one of the following options:',
     +  'Q = save current operator and quit',
     +  'L = list settings',
     +  'N = save current operator and read in a new one',
     +  'T = do a 3D translation search',
     +  'R = do a 3D rotation search',
     +  'A = do an automatic R/T-search',
     +  'P = do a pseudo-6D R/T-search',
     +  '6 = do a complete 6D R/T-search',
     +  ' '
c
      call textin (' Option ?',t)
      call textut (' Option :',t)
      call upcase (t)
c
      if (index('QNTRAP6L',t) .le. 0) t = 'L'
c
c ... LIST
c
      if (t .eq. 'L') then
c
        write (*,*)
        call prompt (' STATUS:')
        call textut (' Current map :',mapnam)
        call fvalut (' Shift for options T and R  :',1,shift)
        call fvalut (' Translation step (A/P/6)   :',1,tshift)
        call fvalut (' Rotation    step (A/P/6)   :',1,rshift)
        call ivalut (' Sample density             :',1,ndens)
        call textut (' Interpolation method (P/Q) :',ityp)
        call textut (' Sampling method (C/Q)      :',f927)
        call anancs (1,rtbtoa,.true.,ierr)
        write (*,*)
c
        goto 200
      end if
c
c ... QUIT - save improved operator
c
      if (t .eq. 'Q') then
        write (*,*)
        line = 'rt_improved.o'
        call textin (' File for new operator ?',line)
        call textut (' File for new operator :',line)
        write (*,*)
        if (line .eq. ' ') return
c
        call xopxua (1,line,xinter(),errcod)
        if (errcod .ne. 0) then
          call errcon ('While opening O datablock file')
          return
        end if
        write (1,'(a,1x,a)',err=299) par(1:leng1(par)),
     +    'R 12 (3f15.8)'
        write (1,'(3f15.8)',err=299) rtbtoa
        close (1)
 299    return
      end if
c
c ... NEW operator (save old one)
c
      if (t .eq. 'N') then
        write (*,*)
        line = 'rt_improved.o'
        call textin (' File to store CURRENT operator ?',line)
        call textut (' File to store CURRENT operator :',line)
        if (line .eq. ' ') goto 210
        call xopxua (1,line,xinter(),errcod)
        if (errcod .ne. 0) then
          call errcon ('While opening O datablock file')
          return
        end if
        write (1,'(a,1x,a)',err=299) par(1:leng1(par)),
     +    'R 12 (3f15.8)'
        write (1,'(3f15.8)',err=299) rtbtoa
        close (1)
        goto 210
      end if
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
      else if (t .eq. 'P' .or. t .eq. '6') then
c
        call fvalin (' Translation step size ?',1,tshift)
        call fvalut (' Translation step size :',1,tshift)
        if (tshift .le. 0.)  tshift = 0.5
        call fvalin (' Rotation step size    ?',1,rshift)
        call fvalut (' Rotation step size    :',1,rshift)
        if (rshift .le. 0.)  rshift = 2.0
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
      if (niter .eq. 1) then
        call prompt (' Enter the sample density; 1 means: use all')
        call prompt (' points, 3 means: use every third point etc.')
        call prompt (' Higher value -> faster (but less accurate)')
      end if
c
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
      if (niter .eq. 1) then
        call prompt (' Proper interpolation = 64-point spline')
        call prompt (' Q-n-D  interpolation = nearest grid point')
      end if
      call textin (
     +  ' Proper or Quick-n-Dirty interpolation (P/Q) ?',ityp)
      call textut (
     +  ' Proper or Quick-n-Dirty interpolation (P/Q) :',ityp)
      call upcase (ityp)
      if (ityp .ne. 'Q') ityp = 'P'
c
      if (t .eq. 'T' .or. t .eq. 'R' .or. t .eq. 'A') then
        call prompt (' Complete = try all 27 orientations')
        call prompt (' Quick    = try only 8 orientations')
        call textin (' Complete or Quick (C/Q) ?',f927)
        call textut (' Complete or Quick (C/Q) :',f927)
        call upcase (f927)
        if (f927 .ne. 'Q') f927 = 'C'
      else
        f927 = 'C'
      end if
c
c ---	Now go do everything. This is to make subscripting easier
c
      call sub03x (mapa, exta(1), exta(2), exta(3), orgna, 
     $  mapb, maskb, extb(1), extb(2), extb(3), orgnb,
     $  cell, grid, rtbtoa, ctrt, rtsym, ctsym,
     $  t, shift,ndens,tshift,tend,rshift,rend,ccmin,cswap,ityp,
     $  f927)
c
      goto 200
c
      end
c
c
c
      subroutine sub03x 
     $ (mapa, exta1, exta2, exta3, orgna, 
     $ mapb, maskb, extb1, extb2, extb3, orgnb, 
     $ cell, grid, rtbtoa, ctrt, rtsym, ctsym,
     $ orit, shift,ndens,tshift,tend,rshift,rend,ccmin,cswap,ityp,
     $ f927)
c
c ---	Improve the operator by improving the correlation coefficient between
c	the two maps
c ---	Assumes both maps are on same grid.
c ---	Alwyn Jones, 2-Oct-91
c
      implicit none
c
c ---	Map A data structures
      integer orgna(3), exta1, exta2, exta3
      real mapa(exta1, exta2, exta3)
c
c ---	Map B data structure
      integer orgnb(3), extb1, extb2, extb3
      integer maskb(extb1, extb2, extb3)
      real mapb(extb1, extb2, extb3)
c
      real cell(6), grid(3), rtbtoa(12), rtsym(12,*)
      integer ctrt, ctsym,ndens
      character orit*1,t*1,ityp*1,f927*1
      real shift,bvalue,delt(3),delr(3)
c
      integer ctmask,errcod,i,i1,i2,i3,j1,j2,j3,j,k,l,ix(3)
      integer ncnt,ncyc,nerr1,nerr2
      integer i1max,i2max,i3max,j1max,j2max,j3max,m,ext(3)
      real a(3,3), b(3,3),forgn(3),fext(3),rtatob(12)
      real avalue, x(3), x1(3), x2(3),mapbit(4,4,4),cswap,old
      real avx,avxsq,avxy,avy,avysq,ccoef,f,peak,peakrt(12),rt(12)
      real del(3),q(9),qq(9),xmo(3),xp(3),xo(3),fgk(27),ogk(3),oextr
      real sgk(3),fextr,tshift,tend,rshift,rend,ccmin,gext(3)
      integer sumx, sumy, sumz,ngk
      logical init,tdone,rdone,first
c
      data init,first /.false.,.true./
c
      save init,first
c
code ...
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
c ---	A = frac to ang, B=Ang to frac
c
      call orthog (cell, a, 0)
      call orthog (cell, b, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
c
      nerr1 = 0
      nerr2 = 0
c
c ---	Generate envelope edges in fractional cell coords
c
      do 130 i=1,3
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
130     gext(i)  = float(ext(i)+orgna(i)-2)*grid(i)/cell(i)
c
c ... precompute densities inside mask only once !
c
      if (first) then
        first = .false.
        call prompt (' Precomputing density inside mask')
c
        call fvalut (' FORGN :',3,forgn)
        call fvalut (' FEXT  :',3,fext)
        call fvalut (' GEXT  :',3,gext)
        ctmask = 0
c
        do 1600 k=1,extb3
        do 1600 j=1,extb2
        do 1600 i=1,extb1
          mapb(i,j,k) = -1.0E5
          if (maskb(i,j,k) .ne. 1) goto 1600
c
            x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
            x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
            x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
            call mulmtx (a, x, x2, 3, 3, 1)
            call mulmtx (b, x2, x1, 3, 3, 1)
c
c ---	      Is it in the envelope, if not make it (if possible)
c
            call frcsym (x1, forgn, fext, rtsym, ctsym, errcod)
            if (errcod .ne. 0) then
              if (nerr1 .lt. 10) then
                call errcon ('Severe FRCSYM error 1')
                call fvalut (' Mask point :',3,x2)
                call fvalut (' NCS  point :',3,x1)
              else if (nerr1 .eq. 10) then
                call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
              end if
              nerr1 = nerr1 + 1
              goto 1600
            end if
            do l=1,3
              ix (l) = nint(x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1)
            end do
            mapb (i,j,k) = mapa(ix(1),ix(2),ix(3))
            ctmask = ctmask+1
1600    continue
        call jvalut (' Nr of mask points with density :',1,ctmask)
      end if
c
      do 510 j=1,9
510     rtatob(j) = rtbtoa(j)
      call matinv (rtatob, 3, x, x1, x2)
      call mulmtx (rtatob, rtbtoa(10), rtatob(10), 3, 3, 1)
      do 520 j=10,12
520     xo(j-9) = -rtatob(j)
      write (*,30) xo
30    format (' Rotation origin at ',3f10.3)
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
c
      x(1) = (sumx/ctmask- 1+ orgnb(1))*grid(1)/cell(1)
      x(2) = (sumy/ctmask- 1+ orgnb(2))*grid(2)/cell(2)
      x(3) = (sumz/ctmask- 1+ orgnb(3))*grid(3)/cell(3)
      call mulmtx (a, x, xmo, 3, 3, 1)
      write (*,40) xmo
40    format (' Centre of mask ', 3f10.3)
      call mulmtx (rtbtoa, xmo, xp, 3, 3, 1)
      do 610 i=1,3
610     xp(i) = xp(i)+ rtbtoa(i+9)
c
      peak = -999.
      i1max = 0
      i2max = 0
      i3max = 0
      j1max = 0
      j2max = 0
      j3max = 0
c
      ngk = 0
      call flusho (6)
c
      if (orit .eq. '6') goto 3000
c
      do 300 i1=1,3
      do 300 i2=1,3
      do 300 i3=1,3
c
c ... check complete or quick
c
        if (f927 .eq. 'Q' .and. index ('R/T/A',orit) .gt. 0) then
          if (i1.eq.2 .or. i2.eq.2 .or. i3.eq.2) then
            if (.not.(i1.eq.2.and.i2.eq.2.and.i3.eq.2)) goto 300
          end if
        end if
c
        do i=1,12
          rt(i) = rtbtoa(i)
        end do
c
        if (orit .eq. 'P') then
c
c ... translation part
c
          delt(1) = (i1-2)*tshift
          delt(2) = (i2-2)*tshift
          delt(3) = (i3-2)*tshift
          do j=10,12
            rt(j) = rt(j)+ delt(j-9)
          end do
c
c ... rotation part
c
          delr(1) = (i1-2)*rshift
          delr(2) = (i2-2)*rshift
          delr(3) = (i3-2)*rshift
          do i=1,3
            call matrot (i, delr(i), q)
            call mulmtx (q, rt, qq, 3, 3, 3)
            do j=1,9
              rt(j) = qq(j)
            end do
          end do
          call mulmtx (rt, xmo, rt(10), 3, 3, 1)
          do j=10,12
            rt(j) = -rt(j)+ xp(j-9)
          end do
c
        else
          del(1) = (i1-2)*shift
          del(2) = (i2-2)*shift
          del(3) = (i3-2)*shift
          if (t .eq. 'T') then
            do j=10,12
              rt(j) = rt(j)+ del(j-9)
            end do
          else
            do i=1,3
              call matrot (i, del(i), q)
              call mulmtx (q, rt, qq, 3, 3, 3)
              do j=1,9
                rt(j) = qq(j)
              end do
            end do
            call mulmtx (rt, xmo, rt(10), 3, 3, 1)
            do j=10,12
              rt(j) = -rt(j)+ xp(j-9)
            end do
          end if
        end if
c ---	Set sums to zero
        avx = 0.
        avy = 0.
        avxy = 0.
        avxsq = 0.
        avysq = 0.
        ctmask = 0
        ncnt = 0
c ---	Loop over mask looking for something
        do 100 k=1,extb3
        do 100 j=1,extb2
        do 100 i=1,extb1
          if (maskb(i,j,k) .eq. 1) then
c
            if (mapb(i,j,k) .lt. -0.99E5) goto 100
c
c ... check if this point needs to be used (sample density)
c
            if (ndens .gt. 1) then
              ncnt = ncnt + 1
              if (ncnt .lt. ndens) goto 100
              ncnt = 0
            end if
c
            x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
            x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
            x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
c
            call mulmtx (a, x, x2, 3, 3, 1)
c
c ... generate NCS point
c
            call vecrtv (x2, x, 1, rt(1), rt(10))
            call mulmtx (b, x, x1, 3, 3, 1)
c
c ---	    Is it in the envelope, if not make it (if possible)
c
            call frcsym (x1, forgn, gext, rtsym, ctsym, errcod)
            if (errcod .ne. 0) then
	      call bldbit (mapa, exta1, exta2, exta3, orgna, grid,
     $	        cell, forgn, fext, rtsym, ctsym, mapbit, x1, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call errcon ('Severe FRCSYM error 2')
                  call fvalut (' Mask point :',3,x2)
                  call fvalut (' NCS  point :',3,x1)
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
              if (ityp .eq. 'P') then
  	        call intrpl (mapbit, 4, 4, 4, x, bvalue, errcod)
              else
  	        call intrp1 (mapbit, 4, 4, 4, x, bvalue, errcod)
              end if
            else
              do 110 l=1,3
110             x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c
c ---	      x is now set of indices to the mapa array
c
              if (ityp .eq. 'P') then
                call intrpl 
     $            (mapa, exta1, exta2, exta3, x, bvalue, errcod)
              else
                call intrp1 
     $            (mapa, exta1, exta2, exta3, x, bvalue, errcod)
              end if
            end if
c
            if (errcod .eq. 0) then
              avalue = mapb(i,j,k)
              avx = avx+ avalue
              avy = avy+ bvalue
              avxsq = avxsq+ avalue**2
              avysq = avysq+ bvalue**2
              avxy = avxy+ avalue*bvalue
              ctmask = ctmask+ 1
  	      else
              if (nerr2 .lt. 10) then
                call errcon ('Interpolation error')
              else if (nerr2 .eq. 10) then
                call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
            end if
          end if
100     continue
c
        f = float (ctmask)
        ccoef = (avxy/f- avx*avy/(f*f))/ 
     $    (sqrt(avxsq/f- (avx/f)**2)* sqrt(avysq/f- (avy/f)**2))
c
        if (orit .eq. 'P') then
          write (*,11) delt,delr,ccoef
        else
          write (*,10) del, ccoef
        end if
c
c ... 950208 - changed for the case of a quick scan
cc        ngk = ngk + 1
        ngk = 9*(i1-1) + 3*(i2-1) + i3
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
c
c ... if max not at 0,0,0 accept it as is, also if 
c     CC >= 0.999
c
      if (i1max.ne.2 .or. i2max.ne.2 .or. i3max.ne.2 .or.
     +    peak .ge. 0.999) then
        do 1400 i=1,12
 1400     rtbtoa(i) = peakrt(i)
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
c ... else do interpolation, UNLESS 6D or Quick search
c
      if (orit .eq. 'P') then
        write (*,20) rtbtoa,peak
        call prompt (' Interpolation not reliable for 6D search !!!')
        return
      end if
c
      if (f927 .eq. 'Q') then
        write (*,20) rtbtoa,peak
        call prompt (' Interpolation not possible for Quick search !!!')
        return
      end if
c
      oextr = fgk(14)
      do i=1,3
        ogk (i) = 0.0
      end do
c
      call svdmx3 (fgk,ogk,oextr,sgk,fextr,init)
      write (*,69) 'Central value ',oextr,' at ',ogk
      write (*,69) 'Interpolated  ',fextr,' at ',(shift*sgk(i),i=1,3)
c
   69 format (1x,a,f10.5,a,3f10.5)
c
      if (fextr .lt. oextr) then
        write (*,*) 'Oops - corr coeff goes up; abort search'
        return
      end if
c
      do 1310 i=1,12
1310    rt(i) = rtbtoa(i)
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
400     rtbtoa(i) = rt(i)
c
      write (*,25) rtbtoa,fextr
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
c ... 6D search
c
 3000 continue
c
      do 1300 i1=1,3
      do 1300 i2=1,3
      do 1300 i3=1,3
      do 1300 j1=1,3
      do 1300 j2=1,3
      do 1300 j3=1,3
c
        do i=1,12
          rt(i) = rtbtoa(i)
        end do
c
        ngk = ngk + 1
c
c ... translation part
c
          delt(1) = (i1-2)*tshift
          delt(2) = (i2-2)*tshift
          delt(3) = (i3-2)*tshift
          do j=10,12
            rt(j) = rt(j)+ delt(j-9)
          end do
c
c ... rotation part
c
          delr(1) = (j1-2)*rshift
          delr(2) = (j2-2)*rshift
          delr(3) = (j3-2)*rshift
          do i=1,3
            call matrot (i, delr(i), q)
            call mulmtx (q, rt, qq, 3, 3, 3)
            do j=1,9
              rt(j) = qq(j)
            end do
          end do
          call mulmtx (rt, xmo, rt(10), 3, 3, 1)
          do j=10,12
            rt(j) = -rt(j)+ xp(j-9)
          end do
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
c ---	Loop over mask looking for something
        do 1100 k=1,extb3
        do 1100 j=1,extb2
        do 1100 i=1,extb1
          if (maskb(i,j,k) .eq. 1) then
c
            if (mapb(i,j,k) .lt. -0.99E5) goto 1100
c
c ... check if this point needs to be used (sample density)
c
            if (ndens .gt. 1) then
              ncnt = ncnt + 1
              if (ncnt .lt. ndens) goto 1100
              ncnt = 0
            end if
c
            x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
            x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
            x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
c
            call mulmtx (a, x, x2, 3, 3, 1)
c
c ... generate NCS point
c
            call vecrtv (x2, x, 1, rt(1), rt(10))
            call mulmtx (b, x, x1, 3, 3, 1)
c
c ---	    Is it in the envelope, if not make it (if possible)
c
            call frcsym (x1, forgn, gext, rtsym, ctsym, errcod)
            if (errcod .ne. 0) then
	      call bldbit (mapa, exta1, exta2, exta3, orgna, grid,
     $	        cell, forgn, fext, rtsym, ctsym, mapbit, x1, errcod)
              if (errcod .ne. 0) then
                if (nerr1 .lt. 10) then
                  call errcon ('Severe FRCSYM error 2')
                  call fvalut (' Mask point :',3,x2)
                  call fvalut (' NCS  point :',3,x1)
                else if (nerr1 .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nerr1 = nerr1 + 1
                goto 1100
              end if
	      do 1125 l=1,3
	        x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
	        m = x(l)
	        x(l) = x(l)- float(m-1)
1125          continue
              if (ityp .eq. 'P') then
  	        call intrpl (mapbit, 4, 4, 4, x, bvalue, errcod)
              else
  	        call intrp1 (mapbit, 4, 4, 4, x, bvalue, errcod)
              end if
            else
              do 1110 l=1,3
1110            x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c
c ---	      x is now set of indices to the mapa array
c
              if (ityp .eq. 'P') then
                call intrpl 
     $            (mapa, exta1, exta2, exta3, x, bvalue, errcod)
              else
                call intrp1 
     $            (mapa, exta1, exta2, exta3, x, bvalue, errcod)
              end if
            end if
c
            if (errcod .eq. 0) then
              avalue = mapb(i,j,k)
              avx = avx+ avalue
              avy = avy+ bvalue
              avxsq = avxsq+ avalue**2
              avysq = avysq+ bvalue**2
              avxy = avxy+ avalue*bvalue
              ctmask = ctmask+ 1
	      else
              if (nerr2 .lt. 10) then
                call errcon ('Interpolation error')
              else if (nerr2 .eq. 10) then
                call prompt (
     +          ' NOTE: further interpolation errors but not listed!!!')
              end if
              nerr2 = nerr2 + 1
            end if
          end if
1100    continue
c
        f = float (ctmask)
        ccoef = (avxy/f- avx*avy/(f*f))/ 
     $    (sqrt(avxsq/f- (avx/f)**2)* sqrt(avysq/f- (avy/f)**2))
c
        if (ccoef .gt. peak) then
          do 1330 l=1,12
1330        peakrt(l) = rt(l)
          peak = ccoef
          write (*,11) delt,delr,peak
          i1max = i1
          i2max = i2
          i3max = i3
          j1max = j1
          j2max = j2
          j3max = j3
        else if (i1.eq.0 .and. i2.eq.0 .and. i3.eq.0) then
          if (j1.eq.0 .and. j2.eq.0 .and. j3.eq.0) then
            write (*,11) delt,delr,ccoef
          end if
	end if
        call flusho (6)
1300   continue
c
      write (*,*) 'Nr of mask points checked :',ctmask
      write (*,*) 'Nr of orientations tried  :',ngk
c
c ... if max not at 0,0,0,0,0,0 accept it as is, also if 
c     CC >= 0.999
c
      if (i1max.ne.2 .or. i2max.ne.2 .or. i3max.ne.2 .or.
     +    j1max.ne.2 .or. j2max.ne.2 .or. j3max.ne.2 .or.
     +    peak .ge. 0.999) then
        do 1500 i=1,12
 1500     rtbtoa(i) = peakrt(i)
        write (*,20) peakrt, peak
      else
        write (*,20) rtbtoa,peak
        call prompt (' Interpolation not done for 6D search !!!')
        return
      end if
c
      return
c
10    format (' Shift=',3f12.6, ' | Corr coeff=',f12.6)
11    format (' T=',3f8.4,' R=',3f8.4,' | CC=',f10.6)
20    format (' Best Rotation Matrix'/3(1x,3f12.6/),
     +  ' Best Translation'/1x,3f12.6/
     +  ' Correlation Coefficient = ',f12.6)
25    format (' Best Rotation Matrix'/3(1x,3f12.6/),
     +  ' Best Translation'/1x,3f12.6/
     +  ' Predicted Correlation Coefficient = ',f12.6)
c
      end
