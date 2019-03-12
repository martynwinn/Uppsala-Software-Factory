      program comdem
c
c ---	COMbine Multiple Domain Extended Maps
c ---	Must have same grid, origin, extent and cell
c ...   Gerard Kleywegt @ 931213
c
      implicit none
c
      include 'maxdim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'COMDEM', vers = '040701/3.0.1')
c
      integer maxsiz
      parameter (maxsiz = maxgk1)
c
c      pointer (iaptr,mapa)
c      pointer (ibptr,mapb)
c      pointer (icptr,mapc)
c
c      real mapa(1), mapb(1), mapc(1)
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
      integer nb, mapsize
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
c ... WRDBYT accounts for 4 or 8 bytes per word
c
      nb = wrdbyt*mapsize
      iaptr = fmalloc (nb)
      ibptr = fmalloc (nb)
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
      call docomd (%val(iaptr),%val(ibptr),%val(icptr), mapsize)
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
      subroutine docomd (mapa, mapb, mapc, maxsiz)
c
c ---	COMbine Multiple Domain Extended Maps
c ---	Must have same grid, origin, extent and cell
c ...   Gerard Kleywegt @ 931213
c
      implicit none
c
      include 'maxdim.incl'
c
      integer maxsiz
c      parameter (maxsiz = maxgk1)
c
c ---	Map A data structures
      integer orgna(3), exta(3), grida(3), uvwa(3)
      real mapa(maxsiz), cella(6)
c
c ---	Map B data structure
      integer orgnb(3), extb(3), gridb(3), uvwb(3)
      real mapb(maxsiz), cellb(6)
c
c ---   Map C
      real mapc(maxsiz)
c
      character file*80,newmap*80
c
      integer i,spgrp,spgrpb,nmap,length,nsiz
c
      real prod,wgta,wgtb,sumwgt,x1,x2,x3,x4,x5
      real avea,vara,sdva,aveb,varb,sdvb,xmin,xmax,xtot
      real x6,x7,x8,x9
c
code ...
c
c      call gainit (prognm,vers)
c
c ... 950118 - check if CCP4_OPEN is defined; crash if not
c
      call gkccp4 (.true.)
c
c      call jvalut (' Max size of maps :',1,maxsiz)
c
      sumwgt = 0.0
c
      write (*,*)
      file = ' '
      call textin (' First domain map ?',file)
      call textut (' First domain map :',file)
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
c ---	read in map a
      call edin (file, 1, mapa,  orgna, exta, grida, uvwa,
     $           cella, spgrp, maxsiz)
      close (1)
      call telmap (grida,orgna,exta,cella)
c
      nmap = 1
      wgta = 1.0
      nsiz = exta(1)*exta(2)*exta(3)
      call xstats (mapa,nsiz,avea,sdva,xmin,xmax,xtot)
      vara = sdva*sdva
      write (*,6010) nmap,nsiz,xmin,xmax,xtot,avea,vara,sdva
      call fvalin (' Weight ?',1,wgta)
      call fvalut (' Weight :',1,wgta)
      sumwgt = sumwgt + wgta
c
      do i=1,nsiz
        mapc (i) = wgta * mapa(i)
      end do
c
 6010 format (' Map nr ',i3,5x,' Size ',i10/
     +        ' ED min, max, total ',1p,3e12.4/
     +        ' ED ave, var, stdev ',3e12.4/)
c
c ... get file name of next map
c
 1000 continue
      write (*,*)
      call ivalut (' Map nr :',1,(nmap+1))
      file = ' '
      call textin (' Domain map file name ?',file)
      call textut (' Domain map file name :',file)
c
      if (length(file) .lt. 1) then
        write (*,*) 'No more maps to add'
        goto 2000
      end if
c
c ---	read in map b
      call edin (file, 2, mapb,  orgnb, extb, gridb, uvwb,
     $           cellb, spgrpb, maxsiz)
      close (2)
      call telmap (gridb,orgnb,extb,cellb)
c
      do 100 i=1,3
        if (exta (i) .ne. extb (i)) goto 110
        if (orgna(i) .ne. orgnb(i)) goto 110
        if (grida(i) .ne. gridb(i)) goto 110
        if (abs(cella(i)-cellb(i)) .gt. 0.01) goto 110
        if (abs(cella(i+3)-cellb(i+3)) .gt. 0.01) goto 110
100   continue
c
      nmap = nmap + 1
      call xstats (mapb,nsiz,aveb,sdvb,xmin,xmax,xtot)
      varb = sdvb*sdvb
      write (*,6010) nmap,nsiz,xmin,xmax,xtot,aveb,varb,sdvb
c
      wgtb = 1.0
      call fvalin (' Weight ?',1,wgtb)
      call fvalut (' Weight :',1,wgtb)
c
      sumwgt = sumwgt + wgtb
c
      call xystat (mapa,mapb,nsiz,x1,x2,x3,x4,x5,x6,x7,x8,x9)
      call rvalut (' RMSD map A / B   :',1,x1)
      call rvalut (' R-factor (A)     :',1,x4)
      call rvalut (' Ditto, B scaled  :',1,x6)
      call rvalut (' Scale for B      :',1,x8)
      call rvalut (' R-factor (B)     :',1,x5)
      call rvalut (' Ditto, A scaled  :',1,x7)
      call rvalut (' Scale for A      :',1,x9)
      call fvalut (' Corr coeff       :',1,x3)
      call fvalut (' Shape similarity :',1,x2)
c
c ---	Now go do everything. This is to make subscripting easier
c
      write (*,*) 'Busy ...'
      call sub04x (mapc, mapb, exta(1), exta(2), exta(3), wgtb)
      goto 1000
c
c ---	When finished, mapc needs outputting
c
 2000 continue
c
      write (*,*)
c
      if (nmap .lt. 2) then
        call errcon ('You must supply AT LEAST TWO maps !')
        return
      end if
c
      prod = float(nmap)/sumwgt
c
      call rvalut (' Sum of weights / Nr of maps :',1,
     +  (sumwgt/float(nmap)))
      call prompt (' Dividing by this number ...')
      do i=1,nsiz
        mapc (i) = mapc (i) * prod
      end do
c
c ... reset background
c
      call sub05x (mapc,exta(1),exta(2),exta(3))
c
      call xstats (mapc,nsiz,aveb,sdvb,xmin,xmax,xtot)
      varb = sdvb*sdvb
      write (*,6010) 0,nsiz,xmin,xmax,xtot,aveb,varb,sdvb
c
      newmap = ' '
      call textin (' New summed CCP4 map ?',newmap)
      call textut (' New summed CCP4 map :',newmap)
c
      if (length(newmap) .gt. 0) then
        call edout 
     $    (newmap, 3, mapc, orgna, exta, grida, uvwa, cella, spgrp)
      else
        call prompt (' No filename; no map written')
      end if
c
      return
c
110   continue
      call errcon ('Different origin, extent, grid, cell')
      call errcon ('This map will NOT be added')
      goto 1000
c
      end
c
c
c
      subroutine sub04x (mapa, mapb, exta1, exta2, exta3, prod)
c
c ---	Add 2 identical sized maps
c ---	Alwyn Jones, 5-Nov-91
c
      implicit none
c
c ---	Map A data structures
      integer exta1, exta2, exta3
      real mapa(exta1, exta2, exta3)
c
c ---	Map B data structure
      real mapb(exta1, exta2, exta3)
      real prod
c
      integer i, j, k
c
code ...
c
      do 100 k=1,exta3
      do 100 j=1,exta2
      do 100 i=1,exta1
100     mapa(i,j,k) = mapa(i,j,k) + prod*mapb(i,j,k)

      return
      end
c
c
c
      subroutine sub05x (mapa, exta1, exta2, exta3)
c
c ---	set background
c ---	Alwyn Jones, 5-Nov-91
c
      implicit none
c
c ---	Map A data structures
      integer exta1, exta2, exta3
      real mapa(exta1, exta2, exta3)
c
      real totden,avemsk,avebck
c
      integer i, j, k, n1, n2
c
code ...
c
      n1 = 0
      n2 = 0
      totden = 0.0
c
      call prompt (' Calculating background density ...')
      do 100 k=1,exta3
      do 100 j=1,exta2
      do 100 i=1,exta1
        if (mapa(i,j,k) .ne. 0.0) then
          n1 = n1 + 1
          totden = totden + mapa(i,j,k)
        else
          n2 = n2 + 1
        end if
  100 continue
c
      call jvalut (' Nr of points set         :',1,n1)
      call jvalut (' Nr of background points  :',1,n2)
c
      if (n1 .le. 0 .or. n2 .le. 0) then
        call errcon ('Can not calculate background')
        if (n1 .eq. 0) call prompt (
     +    ' No map points with density ???')
        if (n2 .eq. 0) call prompt (
     +    ' No background points in map ???')
        return
      end if
c
      avemsk = totden / float (n1)
      avebck = -1.0 * totden / float(n2)
      call rvalut (' Sum of density set       :',1,totden)
      call rvalut (' Average density in masks :',1,avemsk)
      call rvalut (' Average background level :',1,avebck)
c
      call prompt (' Setting background density ...')
      do 200 k=1,exta3
      do 200 j=1,exta2
      do 200 i=1,exta1
        if (mapa(i,j,k) .eq. 0.0) mapa(i,j,k) = avebck
  200 continue
c
      return
      end
