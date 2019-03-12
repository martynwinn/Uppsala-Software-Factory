c
c ... mask_subs.f - subroutines for MAMA
c
c *** 961126 - this file now contains all subroutines which DO NOT
c              use the MASKMAN.INCL file (this was easier while
c              implementing dynamic memory allocation)
c
      subroutine sub06m (maska, exta1, exta2, exta3, orgna, grida,
     +                   maskb, extb1, extb2, extb3, orgnb, gridb,
     +                   cell,ierr)
c
c ... copy maskA into maskB, where both have the same cell, but
c     (possibly) different origin, extent and grid
c ... written by Alwyn
c
      implicit none
c
      integer exta1, exta2, exta3, orgna(3), grida(3)
      integer extb1, extb2, extb3, orgnb(3), gridb(3)
      integer maska(exta1, exta2, exta3)
      integer maskb(extb1, extb2, extb3)
      integer i, ii, iii, i0, j, jj, jjj, k, kk, kkk, ierr
c
      real cell,c(3)
c
code ...
c
      ierr = 0
c
      do 100 k=1,extb3
        do 100 j=1,extb2
          do 100 i=1,extb1
            maskb(i,j,k) = 0
100   continue
c
      do i=1,3
        c(i) = float(gridb(i))/float(grida(i))
      end do
      if (max(c(1),c(2),c(3)) .gt. 1.7) then
        call errcon ('Grids are too different !')
        call prompt (' Use the NEw UNit_cell command instead')
        ierr = 1
        return
      end if
c
      if (grida(1) .eq. gridb(1) .and.
     +    grida(2) .eq. gridb(2) .and.
     +    grida(3) .eq. gridb(3)) then
	  i0 = 0
      else
        i0 = 1
      end if
c
      do 110 k=1,exta3
        if (i0 .eq. 1) then
          kk = int (float(k-1+ orgna(3))* c(3)- float(orgnb(3)-1))
        else
          kk = nint (float(k-1+ orgna(3))* c(3)- float(orgnb(3)-1))
        end if
        do 110 j=1,exta2
          if (i0 .eq. 1) then
            jj = int (float(j-1+ orgna(2))* c(2)- float(orgnb(2)-1))
          else
            jj = nint (float(j-1+ orgna(2))* c(2)- float(orgnb(2)-1))
          end if
          do 110 i=1,exta1
            if (maska(i,j,k) .eq. 0) goto 110
            if (i0 .eq. 1) then
              ii = int (float(i-1+ orgna(1))* c(1)- float(orgnb(1)-1))
            else
              ii = nint (float(i-1+ orgna(1))* c(1)- float(orgnb(1)-1))
            end if
            do 120 kkk=kk,kk+i0
              do 120 jjj=jj,jj+i0
                do 120 iii=ii,ii+i0
                  if (kkk .lt. 1) goto 120
                  if (jjj .lt. 1) goto 120
                  if (iii .lt. 1) goto 120
                  if (kkk .gt. extb3) goto 120
                  if (jjj .gt. extb2) goto 120
                  if (iii .gt. extb1) goto 120
                  maskb(iii, jjj, kkk) = 1
120         continue
110   continue
c
      return
      end
c
c
c
      subroutine sub05x 
     $ (imapa, exta1, exta2, exta3, orgna, 
     $  maskb, extb1, extb2, extb3, orgnb, 
     $  cell, grid, rtbtoa, ctrt, rtsym, ctsym,nbad,ierr)
c
c --- Overlap computation using LABELling of asymm unit points
c
c --- Alwyn Jones, 21-Feb-91
c
      implicit none
c
      integer maxcnt,mxm1
      parameter (maxcnt = 1001)
      parameter (mxm1 = maxcnt - 1)
c
c ---	Map A data structures
      integer orgna(3), exta1, exta2, exta3
      integer imapa(exta1, exta2, exta3)
c ---	Map B data structure
      integer orgnb(3), extb1, extb2, extb3
      integer maskb(extb1, extb2, extb3)
      real cell(6), grid(3), rtbtoa(12,*), rtsym(12,*)
      integer ctrt, ctsym, ierr, nerr
c
      integer ct, ctup(maxcnt), errcod, i, j, k, l, loop,nbad,ll
      integer i1, j1, k1, ijk1(3), ina(3), nnot
      integer ext(3), grida(3), v, value(31)
      real a(3,3), b(3,3), forgn(3), fext(3),x(3), x1(3), x2(3)
c
      equivalence (ijk1(1), i1), (ijk1(2), j1), (ijk1(3), k1)
c
code ...
c
      ierr = -1
      nerr = 0
c
      call ivalut (' Nr of symmetry operators :',1,ctsym)
      call ivalut (' Nr of NCS RT   operators :',1,ctrt)
c
c ---	value is the masking value associated with each rt
c
      value(1) = 1
      do i=2,31
        value(i) = value(i-1)*2
      end do
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
c ---	Generate envelope edges in fractional cell coords
c
      do i=1,3
        grida(i) = nint(cell(i)/grid(i))
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
      end do
c
      write (*,*) 'Checking the following part of the unit cell :'
      call fvalut (' Start in fractional coordinates :',3,forgn)
      call fvalut (' End   in fractional coordinates :',3,fext)
      write (*,*) '... Busy ...'
c
c ---	Zero the A array
c
      do 330 k=1,exta3
      do 330 j=1,exta2
      do 330 i=1,exta1
330     imapa(i,j,k) = 0
c
c ---	Loop over mask looking for something
c
      call prompt (' Expanding mask ...')
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        if (maskb(i,j,k) .eq. 1) then
          x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
          x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
          x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
          call mulmtx (a, x, x2, 3, 3, 1)
c
c ---	    Now loop over the number of operators
c
          do 120 loop=1, ctrt
c
            if (loop .gt. 31) then
              ll = mod (loop, 31)
              if (ll .le. 0) ll = 1
            else
              ll = loop
            end if
c
            call vecrtv (x2, x, 1, rtbtoa(1,loop), rtbtoa(10,loop))
            call mulmtx (b, x, x1, 3, 3, 1)
            do 110 l=1,3
110           x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c
c ---	      x is now set of indices to the mapa array
c
            do 130 l=1,3
130           ina(l) = x(l)
c
            do 150 k1=ina(3),ina(3)+1
            do 150 j1=ina(2),ina(2)+1
            do 150 i1=ina(1),ina(1)+1
c
              call  frcv2 
     $            (ijk1,  imapa, exta1, exta2, exta3,
     $            orgna, grida, forgn, fext, 
     $            rtsym, ctsym, value(ll), errcod)
              nerr = nerr + errcod
c
150         continue
120       continue
        end if
100   continue
c
      call jvalut (' Nr of errors in FRC2V :',1,nerr)
      if (nerr .ne. 0) then
        call errcon ('Probably NOT an asymmetric unit !!!')
        return
      end if
c
      do i=1,maxcnt
        ctup(i) = 0
      end do
c
      nnot = 0
      ll = min(31,ctrt)
      call prompt (' Counting ...')
c
      do 500 k=1,exta3
      do 500 j=1,exta2
      do 500 i=1,exta1
c
        ct = 0
        if (imapa(i,j,k) .gt. 0) then
          do l=1,ll
            v = iand(imapa(i,j,k), value(l))
            if (v .gt. 0) ct = ct+1
          end do
        end if
c
        imapa(i,j,k) = ct
c
        if (imapa(i,j,k) .gt. 0 .and. imapa(i,j,k) .le. mxm1) then 
          ctup(imapa(i,j,k)) = ctup(imapa(i,j,k))+ 1
        else if (imapa(i,j,k) .gt. mxm1) then
          ctup (maxcnt) = ctup (maxcnt) + 1
        else
          nnot = nnot + 1
        end if
c
500   continue
c
      write (*,*) 'Counts for asymmetric unit :'
      nbad = 0
      if (ctup(1) .gt. 0) write (*, 12) ctup(1)
      do i=2,mxm1
        if (ctup(i) .gt. 0) write (*, 10) i, ctup(i)
        nbad = nbad + ctup(i)
      end do
      if (ctup(maxcnt) .gt. 0) write (*, 11) maxcnt, ctup(maxcnt)
      nbad = nbad + ctup(maxcnt)
c
      write (*,13) nnot
      write (*,14) 100.0*float(nnot)/(float(exta1)*
     +             float(exta2)*float(exta3))
c
      ierr = 0
c
      return
c
10    format (' Nr of points set ',i4,' times is ........... ', i10)
11    format (' Nr of points set ',i4,' times or more is ... ', i10)
12    format (' Nr of points set only once is .......... ', i10)
13    format (' Nr of points in solvent areas .......... ', i10)
14    format (' Solvent areas as %-age of asymm. unit .. ', f10.2)
c
      end
c
c
c
      subroutine sub15 
     $ (imapa, exta1, exta2, exta3, orgna, 
     $  maskb, extb1, extb2, extb3, orgnb, 
     $  cell, grid, rtbtoa, ctrt, rtsym, ctsym,nbad,ierr)
c
c --- Overlap computation using COUNTing of asymm unit points
c
      implicit none
c
      integer maxcnt,mxm1
      parameter (maxcnt=1001)
      parameter (mxm1=maxcnt-1)
c
c ---	Map A data structures
      integer orgna(3), exta1, exta2, exta3
      real imapa(exta1, exta2, exta3)
c ---	Map B data structure
      integer orgnb(3), extb1, extb2, extb3
      integer maskb(extb1, extb2, extb3)
      real cell(6), grid(3), rtbtoa(12,*), rtsym(12,*)
      integer ctrt, ctsym, ierr, nerr
c
      integer ct,ctup(maxcnt), errcod, i, j, k, l, loop,nbad
      integer i1, j1, k1, ijk1(3), ina(3), nnot
      integer ext(3), grida(3)
      real a(3,3),b(3,3),forgn(3),fext(3),x(3),x1(3),x2(3)
c
      equivalence (ijk1(1), i1), (ijk1(2), j1), (ijk1(3), k1)
c
code ...
c
      ierr = -1
      nerr = 0
c
      call ivalut (' Nr of symmetry operators :',1,ctsym)
      call ivalut (' Nr of NCS RT   operators :',1,ctrt)
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
c ---	Generate envelope edges in fractional cell coords
c
      do i=1,3
        grida(i) = nint(cell(i)/grid(i))
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
      end do
c
      write (*,*) 'Checking the following part of the unit cell :'
      call fvalut (' Start in fractional coordinates :',3,forgn)
      call fvalut (' End   in fractional coordinates :',3,fext)
      write (*,*) '... Busy ...'
c
c ---	Zero the A array
c
      do 330 k=1,exta3
      do 330 j=1,exta2
      do 330 i=1,exta1
330     imapa(i,j,k) = 0.0
c
c ---	Loop over mask looking for something
c
      call prompt (' Expanding mask ...')
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        if (maskb(i,j,k) .eq. 1) then
          x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
          x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
          x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
          call mulmtx (a, x, x2, 3, 3, 1)
c
c ---	    Now loop over the number of operators
c
          do 120 loop=1, ctrt
c
            call vecrtv (x2, x, 1, rtbtoa(1,loop), rtbtoa(10,loop))
            call mulmtx (b, x, x1, 3, 3, 1)
            do 110 l=1,3
110           x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c
c ---	      x is now set of indices to the mapa array
c
            do 130 l=1,3
130           ina(l) = x(l)
c
            do 150 k1=ina(3),ina(3)+1
            do 150 j1=ina(2),ina(2)+1
            do 150 i1=ina(1),ina(1)+1
c
              call  frcv3 
     $            (ijk1,  imapa, exta1, exta2, exta3,
     $            orgna, grida, forgn, fext, 
     $            rtsym, ctsym, errcod)
              nerr = nerr + errcod
c
150         continue
120       continue
        end if
100   continue
c
      call jvalut (' Nr of errors in FRC2V :',1,nerr)
      if (nerr .ne. 0) then
        call errcon ('Probably NOT an asymmetric unit !!!')
        return
      end if
c
      do i=1,maxcnt
        ctup(i) = 0
      end do
c
      nnot = 0
      call prompt (' Counting ...')
c
c ... divide all accumulated counts by 8
c
      do 500 k=1,exta3
      do 500 j=1,exta2
      do 500 i=1,exta1
c
        imapa(i,j,k) = 0.1250 * imapa(i,j,k)
        ct = int(imapa(i,j,k) + 0.01)
c
        if (imapa(i,j,k) .gt. 0.0 .and. ct .lt. mxm1) then 
          ctup(ct+1) = ctup(ct+1)+ 1
        else if (ct .ge. mxm1) then
          ctup (maxcnt) = ctup (maxcnt) + 1
        else
          nnot = nnot + 1
        end if
c
500   continue
c
      write (*,*) 'Counts for asymmetric unit :'
      nbad = 0
      do i=1,mxm1
        if (ctup(i) .gt. 0) then
          write (*,20) i-1,i,ctup(i)
          if (i. gt. 2) nbad = nbad + ctup(i)
        end if
      end do
      if (ctup(maxcnt) .gt. 0) then
        write (*,21) ctup(maxcnt)
        nbad = nbad + ctup(maxcnt)
      end if
      write (*,22) nbad
c
      write (*,13) nnot
      write (*,14) 100.0*float(nnot)/(float(exta1)*
     +             float(exta2)*float(exta3))
c
      ierr = 0
c
      return
c
20    format (' Nr of points set ',i4,' to ',i4,' times : ',i10)
21    format (' Nr of points set more often : ',i10)
22    format (' Nr of suspect points        : ',i10)
13    format (' Nr of points in solvent areas .......... ', i10)
14    format (' Solvent areas as %-age of asymm. unit .. ', f10.2)
c
      end
c
c
c
      subroutine frcv2 (i1, map, e1, e2, e3, 
     $      orgn, grid, min, max, 
     $      rt, ct, v, errcod)
c
c ---	FoRCe_VALue into the mask or the asymm unit
c ---	Given the grid index, force via symops
c	into the envelope. Errcod .ne. 0 if impossible
c ---	Alwyn Jones
c
      implicit none
c
      integer ct, e1, e2, e3, errcod, i1(3), grid(3), orgn(3), v
      integer map(e1, e2, e3)
c
      real min(3), max(3), rt(12,*)
c
      integer i, i2(3), j, nok
c
      real x(3), x1(3)
c
code ...
c
      nok = 0
      errcod = 1
c
      do 100 i=1,3
100     x(i) = float(i1(i)+ orgn(i)- 1)/float(grid(i))
c
      do 200 i=1,ct
        call vecrtv (x, x1, 1, rt(1,i), rt(10,i))
c
        do 210 j=1,3
c
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
          if (x1(j) .lt. min(j)) goto 200
c
210     continue
c
          do 250 j=1,3
250         i2(j) = nint(x1(j)*float(grid(j)))- orgn(j) + 1
          if (i2(1) .le. 0 .or. i2(2) .le. 0 .or. i2(3) .le. 0) then
            call errcon ('In FRCV2 (lower bounds)')
            write (*,*) i1, i2, e1, e2, e3
            goto 200
          end if
          if (i2(1) .gt. e1 .or. i2(2) .gt. e2 .or. i2(3) .gt. e3) then
            call errcon ('In FRCV2 (upper bounds)')
            write (*,*) i1, i2, e1, e2, e3
            goto 200
          end if
          map(i2(1), i2(2), i2(3)) = ior (map(i2(1), i2(2), i2(3)), v)
          nok = nok + 1
c
200   continue
c
      if (nok .ge. 1) errcod = 0
c
      return
      end
c
c
c
      subroutine frcv3 (i1, map, e1, e2, e3, 
     $      orgn, grid, min, max, 
     $      rt, ct, errcod)
c
c ---	FoRCe_VALue into the mask or the asymm unit
c ---	Given the grid index, force via symops
c	into the envelope. Errcod .ne. 0 if impossible
c ---	Alwyn Jones
c
      implicit none
c
      integer ct, e1, e2, e3, errcod, i1(3), grid(3), orgn(3)
      real map(e1, e2, e3)
c
      real min(3), max(3), rt(12,*)
c
      integer i, i2(3), j, nok
c
      real x(3), x1(3)
c
code ...
c
      nok = 0
      errcod = 1
c
      do 100 i=1,3
100     x(i) = float(i1(i)+ orgn(i)- 1)/float(grid(i))
c
      do 200 i=1,ct
        call vecrtv (x, x1, 1, rt(1,i), rt(10,i))
c
        do 210 j=1,3
c
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
          if (x1(j) .lt. min(j)) goto 200
c
210     continue
c
          do 250 j=1,3
250         i2(j) = nint(x1(j)*float(grid(j)))- orgn(j) + 1
          if (i2(1) .le. 0 .or. i2(2) .le. 0 .or. i2(3) .le. 0) then
            call errcon ('In FRCV2 (lower bounds)')
            write (*,*) i1, i2, e1, e2, e3
            goto 200
          end if
          if (i2(1) .gt. e1 .or. i2(2) .gt. e2 .or. i2(3) .gt. e3) then
            call errcon ('In FRCV2 (upper bounds)')
            write (*,*) i1, i2, e1, e2, e3
            goto 200
          end if
          map(i2(1), i2(2), i2(3)) = map(i2(1), i2(2), i2(3)) + 1.0
          nok = nok + 1
c
200   continue
c
      if (nok .ge. 1) errcod = 0
c
      return
      end
c
c
c
      subroutine trimov (mask,over,nbrs,na,nb,nc,minov,minbr,ierr)
c
c ... remove mask points which are occupied by at least MINOV
c     NCS-related mates and which have at least MINBR neighbours
c     outside the mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,minov,minbr,ncnt,ctup(10),idum,ierr,i
      integer mask(na,nb,nc),over(na,nb,nc),nbrs(na,nb,nc)
c
code ...
c
      ierr = -1
c
      if (minov .le. 1) return
      if (minbr .lt. 0 .or. minbr .gt. 27) return
c
      do i=1,10
        ctup (i) = 0
      end do
c
      ncnt = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. 1) then
              idum = over(i1,i2,i3)
              if (idum .gt. 0 .and. idum .lt. 10) then
                ctup (idum) = ctup (idum) + 1
              else if (idum .ge. 10) then
                ctup (10) = ctup (10) + 1
              else
                ncnt = ncnt + 1
              end if
            end if
          end do
        end do
      end do
c
      if (ncnt .gt. 0) then
        call errcon ('There are points in the mask with ZERO count !!')
        call jvalut (' Number of them :',1,ncnt)
        return
      end if
c
      write (*,6000) ctup(1)
      ncnt = 0
      do i=2,9
        if (ctup(i) .gt. 0) write (*,6010) (i-1),ctup(i)
        ncnt = ncnt + ctup(i)
      end do
      if (ctup(10) .gt. 0) write (*,6020) ctup(10)
      ncnt = ncnt + ctup(10)
c
      call jvalut (' Total nr of suspicious mask points :',1,ncnt)
      if (ncnt .le. 0) return
c
 6000 format (' Nr of mask points without overlap    = ',i10)
 6010 format (' Nr of mask points overlapping ',i2,' *   = ',i10)
 6020 format (' Nr of mask points overlapping more   = ',i10)
c
      ncnt = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. 1) then
              if (over(i1,i2,i3) .ge. minov) then
                if (nbrs(i1,i2,i3) .ge. minbr) then
                  mask (i1,i2,i3) = 0
                  ncnt = ncnt + 1
                end if
              end if
            end if
          end do
        end do
      end do
c
      call jvalut (' Nr of mask points removed :',1,ncnt)
      ierr = 0
c
      return
      end
c
c
c
      subroutine sub01x (
     $  mapa, exta1, exta2, exta3, orgna, 
     $  mapb, maskb, extb1, extb2, extb3, orgnb, 
     $  cell, grid, rtbtoa, ctrt, rtsym, ctsym, ierr)
c
c ---	Build up a map. Given map a, map b and mask b, and rt_b_to_a
c	For every point in mask b, transform to map a, interpolate
c	onto map a . Add mapa value to map b.
c ---	Assumes both maps are on same grid.
c ---	Does not require a unit operator, i.e. uses map b value.
c ---	Alwyn Jones, 21-Feb-91
c
      implicit none
c
c ---	Map A data structures
c
      integer orgna(3), exta1, exta2, exta3
      real mapa(exta1, exta2, exta3)
c
c ---	Map B data structure
c
      integer orgnb(3), extb1, extb2, extb3
      integer maskb(extb1, extb2, extb3)
      real mapb(extb1, extb2, extb3)
      real cell(6), grid(3), rtbtoa(12,*), rtsym(12,*)
      integer ctrt, ctsym,ierr
c
      integer errcod, i, j, k, l, loop, m, nerr, nsev, nintp
      real a(3,3), b(3,3), forgn(3), fext(3), mapbit(4,4,4)
      real value, x(3), x1(3), x2(3), gext(3)
      integer ext(3),extb(3)
c
      logical xinter,isitok
c
code ...
c
      ierr = -1
c
c ---	A = frac to ang, B=Ang to frac
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
c ---	Generate envelope edges in fractional cell coords
c
      do i=1,3
        forgn(i) = float(orgna(i))         *grid(i)/cell(i)
        fext(i)  = float(ext(i)+orgna(i)-1)*grid(i)/cell(i)
        gext(i)  = float(ext(i)+orgna(i)-2)*grid(i)/cell(i)
      end do
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
      write (*,*) '... Busy ...'
c
c ---	Loop over mask looking for something
c
      nerr = 0
      nsev = 0
      nintp = 0
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        mapb(i,j,k) = 0.0
        if (maskb(i,j,k) .eq. 1) then
          x(1) = (i-1+ orgnb(1))*grid(1)/cell(1)
          x(2) = (j-1+ orgnb(2))*grid(2)/cell(2)
          x(3) = (k-1+ orgnb(3))*grid(3)/cell(3)
          call mulmtx (a, x, x2, 3, 3, 1)
c
c ---	    Now loop over the number of operators
c
          do 120 loop=1, ctrt
            call vecrtv (x2, x, 1, rtbtoa(1,loop), rtbtoa(10,loop))
            call mulmtx (b, x, x1, 3, 3, 1)
c
c ---	      Is it in the envelope, if not make it (if possible)
c
            call frcsym (x1, forgn, gext, rtsym, ctsym, errcod)
            if (errcod .ne. 0) then
c
c ---	Unable to do it. Probably around the edge so build a 4x4x4 map
c
              nerr = nerr + 1
              call bldbit (mapa, exta1, exta2, exta3, orgna, grid,
     $          cell, forgn, fext, rtsym, ctsym, mapbit, x1, errcod)
              if (errcod .ne. 0) then
c
c ---	Still unable to do it.
c
                if (nsev .lt. 10) then
                  call errcon (' Severe FRCSYM error')
                  call fvalut (' Mask point :',3,x2)
                  call fvalut (' NCS point  :',3,x1)
                else if (nsev .eq. 10) then
                  call prompt (
     +              ' NOTE: further FRCSYM errors but not listed!!!')
                end if
                nsev = nsev + 1
                goto 120
              end if
              do 125 l=1,3
                x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
                m = x(l)
                x(l) = x(l)- float(m-1)
125           continue
              call intrpl (mapbit, 4, 4, 4, x, value, errcod)
            else
              do 110 l=1,3
110             x(l) = x1(l)*cell(l)/grid(l)- float(orgna(l))+ 1
c
c ---	        x is now set of indices to the mapa array
c
              call intrpl 
     $          (mapa, exta1, exta2, exta3, x, value, errcod)
            end if
            if (errcod .eq. 0) then
              mapb(i,j,k) =  mapb(i,j,k) + value
            else
              if (nintp .lt. 10) then
                call prompt (' Interpolation error !')
                write (*,'(1x,a,1x,3f8.2,1x,3i4)')
     +            'x & ext',x,exta1,exta2,exta3
              else if (nintp .eq. 10) then
                call prompt (
     +        ' NOTE: further interpolation errors but not listed!!!')
              end if
              nintp = nintp + 1
            end if
120       continue
        else
        end if
100   continue
c
      call jvalut (' Nr of boundary points :',1,nerr)
      call jvalut (' Severe FRCSYM errors  :',1,nsev)
      if (nsev .gt. 0) then
        call errcon ('Probably NOT an asymmetric unit !')
        if (xinter()) then
          if (.not. isitok(' Continue anyway ?')) return
        else
          return
        end if
      end if
c
      if (nintp .gt. 0) then
        call errcon ('There were interpolation errors !')
        call jvalut (' Number of errors :',1,nintp)
        if (xinter()) then
          if (.not. isitok(' Continue anyway ?')) return
        else
          return
        end if
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine trrmov (mask,over,nbrs,na,nb,nc,xminov,minbr,ierr)
c
c ... remove mask points which are occupied by at least XMINOV
c     NCS-related mates and which have at least MINBR neighbours
c     outside the mask
c
      implicit none
c
      integer maxcnt,mxm1
      parameter (maxcnt = 1001)
      parameter (mxm1 = maxcnt - 1)
c
      integer na,nb,nc,i1,i2,i3,minbr,ncnt,ctup(maxcnt),nbad
      integer mask(na,nb,nc),nbrs(na,nb,nc),ntot,ierr,i,idum
c
      real xminov,over(na,nb,nc),rdum,perc
c
code ...
c
      ierr = -1
c
      if (xminov .le. 1.0) return
      if (minbr .lt. 0 .or. minbr .gt. 27) return
c
      do i=1,maxcnt
        ctup (i) = 0
      end do
c
      ncnt = 0
      nbad = 0
      ntot = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. 1) then
              ntot = ntot + 1
              rdum = over(i1,i2,i3)
              idum = int (rdum)
              if (rdum .ge. xminov) nbad = nbad + 1
              if (rdum .gt. 0.0001 .and. rdum .lt. mxm1) then
                ctup (idum+1) = ctup (idum+1) + 1
              else if (rdum .ge. mxm1) then
                ctup (maxcnt) = ctup (maxcnt) + 1
              else
                ncnt = ncnt + 1
              end if
            end if
          end do
        end do
      end do
c
      do i=1,mxm1
        if (ctup(i) .gt. 0) write (*,6010) (i-1),i,ctup(i)
      end do
      if (ctup(maxcnt) .gt. 0) write (*,6020) ctup(maxcnt)
c
      call jvalut (' Total nr of mask points :',1,ntot)
      call jvalut (' Suspicious mask points  :',1,nbad)
      perc = 100.0 * float (nbad) / float (ntot)
      call fvalut (' Percentage suspicious   :',1,perc)
c
      if (nbad .le. 0) return
c
      if (ncnt .gt. 0) then
        call errcon ('There are points in the mask with count < 1 !!')
        call jvalut (' Number of them :',1,ncnt)
      end if
c
 6010 format (' Points with count ',i4,' - ',i4,' = ',i10)
 6020 format (' Points with higher count    = ',i10)
c
      ncnt = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. 1) then
              if (over(i1,i2,i3) .ge. xminov) then
                if (nbrs(i1,i2,i3) .ge. minbr) then
                  mask (i1,i2,i3) = 0
                  ncnt = ncnt + 1
                end if
              end if
            end if
          end do
        end do
      end do
c
      call jvalut (' Nr of mask points removed :',1,ncnt)
      perc = 100.0 * float (ncnt) / float (ntot)
      call fvalut (' Percentage removed        :',1,perc)
      ierr = 0
c
      return
      end
c
c
c
      subroutine lmokay (mask,exta,extb,extc,origin,grid,cell,
     +                   rad,iunit,junit,kunit)
c
      implicit none
c
      real cell(6),rad,b(3,3),g(3),x(3),fake(6),a(3,3)
      real c(3,3),x1(3),xp(3),x2(3),distce
c
      integer exta,extb,extc,origin(3),grid(3),iunit,off,leng1
      integer i,j,k,l,i1,i2,i3,natoms,n1,n2,n3,nok,junit,kunit
      integer mask (exta,extb,extc)
c
      character line*120,atmnam*14
c
code ...
c
c ... foreplay
c
      do i=1,3
        fake(i) = 1.0
        fake(i+3) = cell(i+3)
      end do
      call orthog (fake, a, 1)
      call orthog (fake, c, 0)
      call orthog (cell, b, 1)
      do i=1,3
        g(i) = cell(i)/float(grid(i))
      end do
      off = 1 + int(0.1 + (rad/(min(g(1),g(2),g(3)))))
      natoms = 0
      n1 = 0
      n2 = 0
      n3 = 0
      nok = 0
c
      call fvalut (' Checking atoms with radius :',1,rad)
c
c ... read/check loop
c
   10 continue
      read (iunit,'(a)',end=100,err=9999) line
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') then
        if (junit.gt.0) write (junit,'(a)') line(1:leng1(line))
        if (kunit.gt.0) write (kunit,'(a)') line(1:leng1(line))
        goto 10
      end if
c
      i = natoms + 1
      call upcase (line)
      read (line,6000,err=9999) atmnam,(x(j),j=1,3)
c
 6000 format (12x,a14,4x,3f8.3)
c
      natoms = i
      call mulmtx (a, x, x1, 3, 3, 1)
c
      do 130 j=-off,off	
        do 130 k=-off,off	
          do 130 l=-off,off
            xp(1) = x1(1)+ l*g(1)
            xp(2) = x1(2)+ k*g(2)
            xp(3) = x1(3)+ j*g(3)
            call mulmtx (c, xp, x2, 3, 3, 1)
            if (distce(x,x2) .gt. rad) goto 130
            i1 = nint(xp(1)/g(1))- origin(1) +1
            i2 = nint(xp(2)/g(2))- origin(2) +1
            i3 = nint(xp(3)/g(3))- origin(3) +1
            if (i1 .le. 0) goto 199
            if (i2 .le. 0) goto 199
            if (i3 .le. 0) goto 199
            if (i1 .gt. exta) goto 198
            if (i2 .gt. extb) goto 198
            if (i3 .gt. extc) goto 198
            if (mask(i1,i2,i3) .ne. 1) goto 197
130   continue
c
c ... no error
c
      if (junit.gt.0) write (junit,'(a)') line(1:leng1(line))
      nok = nok + 1
      goto 10
c
c ... below lower grid bounds
c
  199 continue
      if (kunit.gt.0) write (kunit,'(a)') line(1:leng1(line))
      write (*,6010) atmnam
      n1 = n1 + 1
 6010 format (' ERROR - ',a14,' below lower grid bounds')
      goto 10
c
c ... above upper grid bounds
c
  198 continue
      if (kunit.gt.0) write (kunit,'(a)') line(1:leng1(line))
      write (*,6020) atmnam
      n2 = n2 + 1
 6020 format (' ERROR - ',a14,' above upper grid bounds')
      goto 10
c
c ... atom + radius outside current mask
c
  197 continue
      if (kunit.gt.0) write (kunit,'(a)') line(1:leng1(line))
      write (*,6030) atmnam
      n3 = n3 + 1
 6030 format (' ERROR - ',a14,' not covered by the mask')
      goto 10
c
 9999 continue
      call errcon ('While reading PDB file')
c
  100 continue
      call jvalut (' Nr of atoms read        :',1,natoms)
      call jvalut (' Nr covered by the mask  :',1,nok)
c
      call jvalut (' Below lower grid bounds :',1,n1)
      call jvalut (' Above upper grid bounds :',1,n2)
      call jvalut (' Nr not covered by mask  :',1,n3)
c
      if (n3 .gt. 0) then
        if (n1 .gt. 0 .or. n2 .gt. 0) then
          call errcon ('Grid AND mask too tight !')
        else
          call errcon ('Mask too tight !')
        end if
      end if
c
      return
      end
c
c
c
      subroutine island (mask,na,nb,nc)
c
c ... erase isolated "islands" in a mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,from,to,temp,nchn,j1,j2,j3,ncyc
      integer inow,ibest,nbest,nzap,ntot,nleft
      integer mask(na,nb,nc)
c
c      integer i
c
code ...
c
      inow = 0
      nbest = 0
      ibest = 1
c
      from = 1
      ncyc = 0
c
      ntot = 0
      do i1=1,na
        do i2=2,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. 1) ntot = ntot + 1
          end do
        end do
      end do
      nleft = ntot
c
   20 continue
c
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. from) then
c
c ... we found a new mask "drop"; zap it
c
              inow = inow + 1
              to = -inow
              temp = to - 1
              mask (i1,i2,i3) = to
              ncyc = 0
              nzap = 1
              call prompt (' ...')
              call jvalut (' Zapping island nr :',1,inow)
              goto 10
            end if
          end do
        end do
      end do
c
      write (*,*) 'No more islands ...'
      goto 30
c
c ... "zap" outside world
c
   10 continue
      ncyc = ncyc + 1
      nchn = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. to) then
              do j1=i1-1,i1+1
                do j2=i2-1,i2+1
                  do j3=i3-1,i3+1
                    if (mask(j1,j2,j3) .eq. from) then
                      nchn = nchn + 1
                      nzap = nzap + 1
                      mask(j1,j2,j3) = to
                    end if
                  end do
                end do
              end do
              mask (i1,i2,i3) = temp
            end if
          end do
        end do
      end do
c
      if (nchn .gt. 0) goto 10
c
      write (*,6206) ncyc,nzap
 6206 format (' ',i4,' zap cycles - ',i8,' mask points in island')
c
c      call jvalut (' Nr of zap cycles  :',1,ncyc)
c      call jvalut (' Nr of mask points :',1,nzap)
c
      if (nzap .gt. nbest) then
        ibest = inow
        nbest = nzap
c
        if (2*nzap .gt. ntot) then
          call prompt (
     +  ' This island contains more than 50% of all the mask')
          call prompt (
     +  ' points, so it must be the core mask ! Ha !')
          call alterm (mask,na,nb,nc,temp,to)
          call alterm (mask,na,nb,nc,1,0)
          goto 30
        end if
c
      end if
c
      nleft = nleft - nzap
      if (nleft .lt. nbest) then
        call prompt (
     +  ' The largest island found so far contains more mask')
        call prompt (
     +  ' points than there are unexamined mask points left,')
        call prompt (
     +  ' so it must be the core mask ! Ha !')
        call alterm (mask,na,nb,nc,temp,to)
        call alterm (mask,na,nb,nc,1,0)
        goto 30
      end if
c
      call alterm (mask,na,nb,nc,temp,to)
c
      goto 20
c
c ... let's see what happened
c
   30 continue
      call prompt (' ...')
c
      if (inow .lt. 1) then
        call errcon ('Unexpected error - no mask points found')
        return
      end if
c
c ... if one big connected mask, then reset mask points to +1
c
      if (inow .eq. 1) then
        write (*,*) 'No isolated islands found'
        call alterm (mask,na,nb,nc,-1,1)
        return
      end if
c
      call jvalut (' Nr of the core mask :',1,ibest)
      call jvalut (' Nr of points in it  :',1,nbest)
c
      call prompt (' Resetting core mask ...')
      do i1=1,na
        do i2=2,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. -ibest) then
              mask (i1,i2,i3) = 1
            else
              mask (i1,i2,i3) = 0
            end if
          end do
        end do
      end do
c
c      do i=1,inow
c        if (i.eq.ibest) then
c          call jvalut (' Resetting core mask :',1,i)
c          call alterm (mask,na,nb,nc,-i,1)
c        else
c          call jvalut (' Zeroing island      :',1,i)
c          call alterm (mask,na,nb,nc,-i,0)
c        end if
c      end do
c
      call prompt (' ...')
c
      return
      end
c
c
c
      subroutine border (mask,na,nb,nc)
c
c ... check if mask touches borders
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,nbad(6),k(6),ndum,i
      integer mask(na,nb,nc)
c
      character which(6)*7
c
      data which /'Lower X','Upper X','Lower Y','Upper Y',
     +            'Lower Z','Upper Z'/
c
code ...
c
c ... find current padding to all six borders
c
      ndum = max (na+nb+nc,9999)
      do i=1,6
        nbad (i) = ndum
      end do
c
      write (*,'(a,3i5)') ' Extent : ',na,nb,nc
c
      do i1=1,na
        k(1) = i1 - 1
        k(2) = na - i1
        do i2=1,nb
          k(3) = i2 - 1
          k(4) = nb - i2
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. 1) then
              k(5) = i3 - 1
              k(6) = nc - i3
              do i=1,6
                nbad (i) = min (nbad(i),k(i))
              end do
            end if
          end do
        end do
      end do
c
      ndum = 0
      do i=1,6
        write (*,6000) which(i),nbad(i)
        if (nbad(i) .le. 2) then
          write (*,6010) which (i)
          ndum = ndum + 1
        else if (nbad(i) .le. 5) then
          write (*,6020) which (i)
        end if
      end do
c
 6000 format (1x,a,' padding = ',i6,' grid points')
 6010 format (1x,a,' padding MUST be extended !!!')
 6020 format (1x,a,' padding could be extended a bit')
 6030 format (' Your mask touches ',i1,' borders of the grid !')
c
      if (ndum .gt. 0) then
        write (*,6030) ndum
      end if
c
      return
      end
c
c
c
      subroutine fillm (mask,na,nb,nc)
c
c ... fill voids in a mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,from,to,temp,nchn,j1,j2,j3,ncyc
      integer mask(na,nb,nc)
c
code ...
c
c ... test if mask not intersecting borders
c
      i1 = 1
      do i2=1,nb
        do i3=1,nc
          if (mask(i1,i2,i3) .eq. 1) goto 5
        end do
      end do
c
      i1 = na
      do i2=1,nb
        do i3=1,nc
          if (mask(i1,i2,i3) .eq. 1) goto 5
        end do
      end do
c
      i2 = 1
      do i1=1,na
        do i3=1,nc
          if (mask(i1,i2,i3) .eq. 1) goto 5
        end do
      end do
c
      i2 = nb
      do i1=1,na
        do i3=1,nc
          if (mask(i1,i2,i3) .eq. 1) goto 5
        end do
      end do
c
      i3 = 1
      do i1=1,na
        do i2=1,nb
          if (mask(i1,i2,i3) .eq. 1) goto 5
        end do
      end do
c
      i3 = nc
      do i1=1,na
        do i2=1,nb
          if (mask(i1,i2,i3) .eq. 1) goto 5
        end do
      end do
c
      goto 7
c
    5 continue
      call errcon ('Mask not contained within borders; enlarge !')
      write (*,'(a,3i6)') ' Mask set at border point : ',i1,i2,i3
      return
c
c ... set outside world from '0' to '2'
c
    7 continue
      from = 0
      to   = 2
      temp = 3
      ncyc = 0
c
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. from) then
              mask (i1,i2,i3) = to
              goto 10
            end if
          end do
        end do
      end do
c
      call errcon ('Unexpected error in FILLM - contact Gerard')
      return
c
c ... "zap" outside world
c
   10 continue
      ncyc = ncyc + 1
      nchn = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. to) then
              do j1=i1-1,i1+1
                do j2=i2-1,i2+1
                  do j3=i3-1,i3+1
                    if (mask(j1,j2,j3) .eq. from) then
                      nchn = nchn + 1
                      mask(j1,j2,j3) = to
                    end if
                  end do
                end do
              end do
              mask (i1,i2,i3) = temp
            end if
          end do
        end do
      end do
c
      if (nchn .gt. 0) goto 10
c
      call jvalut (' Nr of zap cycles :',1,ncyc)
c
c ... every mask point which is still 0 is inside the protein
c
c ... now set all remaining 0's to 1
c
      call alterm (mask,na,nb,nc,0,1)
c
c ... now reset the outside world
c
      call alterm (mask,na,nb,nc,to,0)
      call alterm (mask,na,nb,nc,temp,0)
c
      return
      end
c
c
c
      subroutine nbrcnt (mask,shadow,na,nb,nc)
c
c ... list counts of nr of nbrs
c
      implicit none
c
      real pi,one3rd
      parameter (pi = 0.50 * 6.2831853071796)
      parameter (one3rd = 1.0/3.0)
c
      integer na,nb,nc,i1,i2,i3,ns,cnt(0:27),ntot,i,j,nsur
      integer mask(na,nb,nc),shadow(na,nb,nc)
c
      real perc,frac,srad,surf,sfra,smoo
c
code ...
c
      call cntmsk (mask,na,nb,nc,ns)
c
c ... first count non-mask nbrs for mask points
c
      call cntnbr (mask,shadow,na,nb,nc,1,0)
c
      do i=0,27
        cnt (i) = 0
      end do
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. 1) then
              j = shadow (i1,i2,i3)
              cnt (j) = cnt (j) + 1
            end if
          end do
        end do
      end do
c
      write (*,'(/a/a)')
     +  ' Count of nr of non-mask neighbours for mask points',
     +  ' Count     Number      Total  %MASK_ACC'
c
      ntot = 0
      do i=26,0,-1
        ntot = ntot + cnt (i)
        if (cnt(i) .gt. 0) then
          perc = 100.0 * float(ntot) / float(ns)
          write (*,'(1x,i5,1x,i10,1x,i10,1x,f10.2)')
     +      i,cnt(i),ntot,perc
        end if
      end do
c
      nsur = ntot - cnt(0)
      frac = 100.0 * float(nsur) / float(ntot)
      srad = (0.75 * float(ntot) / pi) ** one3rd
      surf = 4.0 * pi * (srad**2)
      sfra = 100.0 * surf / float(ntot)
      smoo = frac / sfra
c
      write (*,6000) float(ntot),srad,surf,sfra,frac,smoo
c
c ... now count mask nbrs for non-mask points
c
      call cntnbr (mask,shadow,na,nb,nc,0,1)
c
      do i=0,27
        cnt (i) = 0
      end do
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. 0) then
              j = shadow (i1,i2,i3)
              cnt (j) = cnt (j) + 1
            end if
          end do
        end do
      end do
c
      write (*,'(/a/a)')
     +  ' Count of nr of mask neighbours for non-mask points',
     +  ' Count     Number      Total  %MASK_ACC'
c
      ntot = 0
      do i=26,0,-1
        ntot = ntot + cnt (i)
        if (cnt(i) .gt. 0) then
          perc = 100.0 * float(ntot) / float(ns)
          write (*,'(1x,i5,1x,i10,1x,i10,1x,f10.2)')
     +      i,cnt(i),ntot,perc
        end if
      end do
c
 6000 format (/
     +  ' Volume of mask (voxels)         = ',f15.3/
     +  ' Corresponds to sphere of radius = ',f15.3/
     +  ' This has a surface (voxels)     = ',f15.3/
     +  ' Fraction SURFACE/VOLUME sphere  = ',f15.3/
     +  ' Fraction SURFACE/VOLUME mask    = ',f15.3/
     +  ' Roughness of mask (smooth=1.0)  = ',f15.3)
c
      return
      end
c
c
c
      subroutine newunt (
     +   maska, exta1, exta2, exta3, orgna, grida, cella,
     +   maskb, extb1, extb2, extb3, orgnb, gridb, cellb,
     +   scratch,shadow,rt,pad,maxsiz,ierr)
c
c ... copy maskA into maskB, where both MAY have
c     different origin, extent, grid AND CELL !!!
c ... written by Gerard using Alwyn's SUB06M
c
      implicit none
c
c ... ACCU = required accuracy of actual mask volume when
c            going to another grid in PERCENT units
c ... CLOS = allowable deviation for considering two cell
c            constants to be the same or an RT-operator
c            element to be equal to 1.0 or 0.0 or a grid
c            point to be integer
c
      real accu,clos
      parameter (accu = 0.1, clos = 0.01)
c
      integer maxcnt
      parameter (maxcnt = 1001)
c
      integer exta1, exta2, exta3, orgna(3), grida(3)
      integer extb1, extb2, extb3, orgnb(3), gridb(3)
      integer pad(3),maxsiz,ierr,exta(3),extb(3)
      integer maska(exta1, exta2, exta3)
      integer maskb(maxsiz),scratch(maxsiz),shadow(maxsiz)
      integer i,ii,iii,i0,j,jjj,k,kkk
      integer ndum,i1,i2,i3,j0,k0,nxy,ctr(0:maxcnt),ntar,ncnt
      integer nmask,nbrs(0:26),iimax,ilo,ihi,jlo,jhi,klo,khi
c
      real cella(6),cellb(6),af2c(3,3),ac2f(3,3)
      real bf2c(3,3),bc2f(3,3),xmin(3),xmax(3),x1(3),x2(3)
      real spaca(3),xdum(3),rt(12),voxva,voxvb,volume
      real maskva,maskvb,diff
c
      logical cut
c
code ...
c
      ierr = -1
c
c ... foreplay
c
   69 continue
c
      exta(1) = exta1
      exta(2) = exta2
      exta(3) = exta3
c
      do i=1,3
        spaca(i)  = 1.0/float(grida(i))
      end do
c
c ... set up conversion matrices Fract-to-Cart etc.
c
      call orthog (cella,af2c,0)
      call orthog (cella,ac2f,1)
      call orthog (cellb,bf2c,0)
      call orthog (cellb,bc2f,1)
c
c ... find limits of mask A in Cartesian space after RT
c
      do i=1,3
        xmin (i) =  99999.9
        xmax (i) = -99999.9
      end do
c
c ... get min and max fractional coordinates of actual mask
c     points for the new cell
c
c ... note: previous bug: used min/max Cartesian and converted
c     those to fractional -> is wrong when not orthorhombic !!!
c
      ndum = 0
      do i1=1,exta1
        do i2=1,exta2
          do i3=1,exta3
            if (maska(i1,i2,i3) .eq. 1) then
              ndum = ndum + 1
              xdum (1) = float (i1 + orgna(1) - 1) * spaca(1)
              xdum (2) = float (i2 + orgna(2) - 1) * spaca(2)
              xdum (3) = float (i3 + orgna(3) - 1) * spaca(3)
              call mulmtx (af2c,xdum,x1,3,3,1)
              call vecrtv (x1,x2,1,rt(1),rt(10))
              call mulmtx (bc2f,x2,x1,3,3,1)
              do i=1,3
                xmin(i)=min(xmin(i),x1(i))
                xmax(i)=max(xmax(i),x1(i))
              end do
            end if
          end do
        end do
      end do
c
      call jvalut (' Nr of points in old mask :',1,ndum)
      if (ndum .le. 0) then
        call errcon ('This mask is empty !')
      end if
c
      call ivalut (' Old origin  :',3,orgna)
      call ivalut (' Old extent  :',3,exta)
      call ivalut (' Old grid    :',3,grida)
      call fvalut (' Old cell    :',6,cella)
      call fvalut (' RT operator :',3,rt(1))
      call fvalut (' RT operator :',3,rt(4))
      call fvalut (' RT operator :',3,rt(7))
      call fvalut (' RT operator :',3,rt(10))
c
      call voxvol (cella,grida,volume,voxva)
      maskva = float(ndum) * voxva
c
      call rvalut (' Old unit-cell volume  :',1,volume)
      call rvalut (' Old voxel volume      :',1,voxva)
      call rvalut (' Volume of actual mask :',1,maskva)
c
      do i=1,3
        x1(i)=xmin(i)
        x2(i)=xmax(i)
      end do
c
      call ivalut (' New grid    :',3,gridb)
      call fvalut (' New cell    :',6,cellb)
      call fvalut (' Fractional mask bottom :',3,x1)
      call fvalut (' Fractional mask top    :',3,x2)
      call ivalut (' Padding     :',3,pad)
c
      do i=1,3
        orgnb (i) = nint (x1(i)*float(gridb(i))) - 1
        ndum = nint (x2(i)*float(gridb(i))) + 1
        extb (i) = ndum - orgnb (i) + 1 + 2*max(5,pad (i)) + 1
        orgnb (i) = orgnb (i) - max(5,pad (i))
      end do
c
      call ivalut (' New origin  :',3,orgnb)
      call ivalut (' New extent  :',3,extb)
c
      extb1 = extb(1)
      extb2 = extb(2)
      extb3 = extb(3)
      ndum = extb1*extb2*extb3
      if (ndum .gt. maxsiz) then
        call errcon ('Mask too big !')
        call jvalut (' Required size :',1,ndum)
        call jvalut (' Maximum  size :',1,maxsiz)
        return
      end if
c
c ... now go do it
c
      write (*,*) '... Busy ...'
c
      ndum = 0
      nxy = extb1*extb2
      call inimsk (maskb,extb1,extb2,extb3,0)
      call inimsk (scratch,extb1,extb2,extb3,0)
c
c ... only if cells or grids differ or if RT-operator is
c     not the indentity operator do we have to mark 8
c     points instead of 1
c
      ii = 0
      do i=1,3
        if (abs(cella(i)-cellb(i))     .gt. clos) goto 1110
        if (abs(cella(i+3)-cellb(i+3)) .gt. clos) goto 1110
c
        if (grida(i) .ne. gridb(i)) goto 1110
c
        i1 = 1 + 4*(i-1)
        if (abs(rt(i1)-1.0) .gt. clos) goto 1110
        if (abs(rt(i1+1))   .gt. clos) goto 1110
        if (abs(rt(i1+2))   .gt. clos) goto 1110
        if (abs(rt(i1+3))   .gt. clos) goto 1110
      end do
      write (*,*) 'Same cell and grid; identity operator'
      goto 1120
c
 1110 continue
      ii = 1
      i0 = max(0,int( (cella(1)/float(grida(1))) /
     +                  (cellb(1)/float(gridb(1)))))
      j0 = max(0,int( (cella(2)/float(grida(2))) /
     +                  (cellb(2)/float(gridb(2)))))
      k0 = max(0,int( (cella(3)/float(grida(3))) /
     +                  (cellb(3)/float(gridb(3)))))
      write (*,'(1x,a,3i2)') 'Zap points within : ',i0,j0,k0
c
 1120 continue
cc      call ivalut (' Extra points :',1,ii)
c
      do 112 k=1,exta3
        x1 (3) = float (k + orgna(3) - 1) * spaca(3)
        do 111 j=1,exta2
          x1 (2) = float (j + orgna(2) - 1) * spaca(2)
          do 110 i=1,exta1
            if (maska(i,j,k) .eq. 0) goto 110
            x1 (1) = float (i + orgna(1) - 1) * spaca(1)
c
c ... X1 = fractional coordinates in cell A
c     X2 = fractional coordinates in cell B
c
            call mulmtx (af2c,x1,x2,3,3,1)
            call vecrtv (x2,xdum,1,rt(1),rt(10))
            call mulmtx (bc2f,xdum,x2,3,3,1)
c
            xdum (1) = x2(1)*gridb(1) - orgnb(1) + 1.0
            xdum (2) = x2(2)*gridb(2) - orgnb(2) + 1.0
            xdum (3) = x2(3)*gridb(3) - orgnb(3) + 1.0
c
            if (ii .gt. 0) then
c
              ilo = int(xdum(1)) - i0
              ihi = ilo + 2*i0 + 1
              jlo = int(xdum(2)) - j0
              jhi = jlo + 2*j0 + 1
              klo = int(xdum(3)) - k0
              khi = klo + 2*k0 + 1
c
            else
c
              ilo = int(xdum(1)+clos)
              ihi = ilo
              if (abs(float(ilo)-xdum(1)) .gt. clos) ihi=ihi+1
c
              jlo = int(xdum(2)+clos)
              jhi = jlo
              if (abs(float(jlo)-xdum(2)) .gt. clos) jhi=jhi+1
c
              klo = int(xdum(3)+clos)
              khi = klo
              if (abs(float(klo)-xdum(3)) .gt. clos) khi=khi+1
c
            end if
c
            do 140 kkk=klo,khi
              do 130 jjj=jlo,jhi
                do 120 iii=ilo,ihi
                  call mskset (scratch,extb1,extb2,extb3,iii,jjj,kkk,1)
120             continue
130           continue
140         continue
c
110       continue
111     continue
112   continue
c
      do i=0,maxcnt
        ctr (i) = 0
      end do
      ncnt = extb1*extb2*extb3
      do i=1,ncnt
        i1 = scratch (i)
        if (i1 .lt. maxcnt) then
          ctr (i1) = ctr (i1) + 1
        else
          ctr (maxcnt) = ctr (maxcnt) + 1
        end if
      end do
c
      call voxvol (cellb,gridb,volume,voxvb)
c
      call rvalut (' New unit-cell volume  :',1,volume)
      call rvalut (' New voxel volume      :',1,voxvb)
      ntar = nint (maskva/voxvb)
      call jvalut (' Required mask points  :',1,ntar)
      call fvalut (' Required accuracy (%) :',1,accu)
c
      ndum = 0
      i1 = 999999
      i2 = 0
      i0 = 999999
      j0 = 0
      iimax = 0
c
      do i=maxcnt,1,-1
c
        ndum = ndum + ctr(i)
c
        i3 = abs (ndum-ntar)
        if (i3 .lt. i1) then
          i1 = i3
          i2 = i
        end if
c
        k0 = ndum-ntar
        if (k0 .ge. 0 .and. k0 .lt. i0) then
          i0 = k0
          j0 = i
        end if
c
        if (ctr(i) .gt. 0) then
          if (iimax .eq. 0) iimax = i
          maskvb = float(ndum) * voxvb
          write (*,6010) i,ctr(i),ndum,maskvb
        end if
c
      end do
c
      cut = .false.
c
      if (i2 .gt. 0) then
c
        diff = 100.0 * float(i1) / float(ntar)
        call jvalut (' Best count     :',1,i2)
        call fvalut (' Difference (%) :',1,diff)
c
c ... check if closest is accurate enough
c
        if (abs(diff) .le. accu) then
          write (*,*) 'Accepted'
          ii = i2
          goto 9000
        end if
c
c ... is there a value which gives a slightly bigger mask ?
c
        if (j0 .le. 0) then
          ii = 1
          goto 9000
        end if
c
c ... if here, there is such a value; select it
c
        call jvalut (' First slightly bigger count :',1,j0)
        diff = 100.0 * float(i0) / float(ntar)
        call fvalut (' Difference (%)              :',1,diff)
        cut = .true.
        ii = j0
        goto 9000
c
c ... if here, don't know what to do
c
      else
c
        ii = 1
c
      end if
c
c ... set all mask points with count >= II
c
 9000 continue
      call ivalut (' Setting all points with count >=',1,ii)
      ndum = 0
      do i=1,ncnt
        if (scratch(i) .ge. ii) then
          maskb (i) = 1
          ndum = ndum + 1
        end if
      end do
c
      maskvb = float(ndum) * voxvb
      call rvalut (' Volume of actual mask   :',1,maskvb)
c
      ierr = 0
c
      if (.not. cut) then
        call rvalut (' Volume of original mask :',1,maskva)
        diff = 100.0 * (maskvb-maskva) / maskva
        call fvalut (' Difference (%)          :',1,diff)
        if (abs(diff) .gt. accu) then
          write (*,*) 'That is the best I can do ...'
          if (diff .gt. 0.0) then
            write (*,*) 'CUT/CONTRACT mask yourself'
          else
            write (*,*) 'SMOOTH/EXPAND mask yourself'
          end if
        end if
        return
      end if
c
c ... make mask a bit smaller
c
      write (*,*) 'Cutting a bit off your mask'
      nmask = ndum
c
 1000 continue
c
c ... count nr of non-mask nbrs for each mask point
c
      call cntnbr (maskb,shadow,extb1,extb2,extb3,1,0)
c
      do i=0,26
        nbrs (i) = 0
      end do
c
c ... count how many have >= N nbrs outside mask
c
      do i=1,ncnt
        if (maskb(i) .eq. 1) then
          if (scratch(i) .eq. ii) then
            i1 = shadow(i)
            nbrs (i1) = nbrs (i1) + 1
          end if
        end if
      end do
c
      ndum = 0
      i1 = 1
      do i=26,1,-1
        ndum = ndum + nbrs(i)
        if ((nmask-ndum) .le. ntar) then
          i1 = i
          goto 1100
        end if
      end do
c
 1100 continue
      write (*,6020) ndum,ii,i1
c
c ... remove them
c
      if (ndum .gt. 0) then
        do i=1,ncnt
          if (maskb(i) .eq. 1) then
            if (scratch(i) .eq. ii) then
              if (shadow(i) .ge. i1) then
                maskb (i) = 0
                nmask = nmask - 1
              end if
            end if
          end if
        end do
      end if
c
      call jvalut (' Mask points left :',1,nmask)
c
      if (nmask .gt. ntar .and. ii .lt. iimax) then
        ii = ii + 1
        goto 1000
      end if
c
      maskvb = float(nmask) * voxvb
      call rvalut (' Volume of actual   mask :',1,maskvb)
      call rvalut (' Volume of original mask :',1,maskva)
      diff = 100.0 * (maskvb-maskva) / maskva
      call fvalut (' Difference (%)          :',1,diff)
      if (abs(diff) .gt. accu) then
        write (*,*) 'That is the best I can do ...'
        if (diff .gt. 0.0) then
          write (*,*) 'CUT/CONTRACT mask yourself'
        else
          write (*,*) 'SMOOTH/EXPAND mask yourself'
        end if
      end if
c
      return
c
 6010 format (' Points with count ',i4,' = ',i8,' Cumul = ',
     +  i8,' Mask volume = ',1pe10.3)
 6020 format (' Removing ',i8,' surface points with count = ',
     +  i4,' and >= ',i2,' non-mask nbrs')
c
      end
c
c
c
      subroutine asucpy (m1,r1,mode,e1,e2,e3,m2,f1,f2,f3)
c
c ... asymmetric unit copy minus the borders
c     FROM mask M1 TO mask M2
c
      implicit NONE
c
      integer e1,e2,e3,f1,f2,f3,i,j,k,mode
      real    r1(e1,e2,e3)
      integer m1(e1,e2,e3),m2(f1,f2,f3)
c
code
c
      if (mode .eq. 1) then
c
c ... use the integer array
c
        do i=1,f1
          do j=1,f2
            do k=1,f3
              if (m1(i,j,k) .gt. 0) then
                m2(i,j,k) = 1
              else
                m2(i,j,k) = 0
              end if
            end do
          end do
        end do
c
      else if (mode .eq. 2) then
c
c ... use the real array
c
        do i=1,f1
          do j=1,f2
            do k=1,f3
              if (r1(i,j,k) .ge. 0.5) then
                m2(i,j,k) = 1
              else
                m2(i,j,k) = 0
              end if
            end do
          end do
        end do
c
      else
c
c ... this shouldn't happen
c
        call errcon ('Unknown mode in ASUCPY - call Gerard')
c
      end if
c
      return
      end
c
c
c
      subroutine blob (mask,na,nb,nc,icut,jcut)
c
c ... erase small isolated "blobs" in a mask
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,from,to,temp,nchn,j1,j2,j3,ncyc
      integer inow,nzap,i,icut,j,k,jcut
      integer mask(na,nb,nc)
c
code ...
c
      inow = 0
      from = 1
c
   20 continue
c
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. from) then
c
c ... we found a new mask blob; zap it
c
              inow = inow + 1
              to = -inow
              temp = to - 1
              mask (i1,i2,i3) = to
              ncyc = 0
              nzap = 1
              call jvalut (' Zapping blob nr :',1,inow)
              goto 10
            end if
          end do
        end do
      end do
c
      write (*,*) 'No more blobs ...'
      goto 30
c
c ... "zap" blob
c
   10 continue
      ncyc = ncyc + 1
      nchn = 0
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            if (mask(i1,i2,i3) .eq. to) then
              do j1=i1-1,i1+1
                do j2=i2-1,i2+1
                  do j3=i3-1,i3+1
                    if (mask(j1,j2,j3) .eq. from) then
                      nchn = nchn + 1
                      nzap = nzap + 1
                      mask(j1,j2,j3) = to
                    end if
                  end do
                end do
              end do
              mask (i1,i2,i3) = temp
            end if
          end do
        end do
      end do
c
      if (nchn .gt. 0) goto 10
c
      call jvalut (' Nr of zap cycles  :',1,ncyc)
      call jvalut (' Nr of mask points :',1,nzap)
c
      if (nzap .lt. icut) then
        call prompt (' Blob too small - erasing it ...')
        call alterm (mask,na,nb,nc,temp,0)
      else if (jcut. gt. icut .and. nzap .gt. jcut) then
        call prompt (' Blob too big - erasing it ...')
        call alterm (mask,na,nb,nc,temp,0)
      else
        call prompt (' Blob size okay - keeping it ...')
        call alterm (mask,na,nb,nc,temp,to)
      end if
c
      goto 20
c
c ... let's see what happened
c
   30 continue
c
      if (inow .lt. 1) then
        call errcon ('Unexpected error - no mask points found')
        return
      end if
c
c ... reset everything < 0 to +1
c
      call prompt ('Resetting erased blobs ...')
      do i=1,na
        do j=1,nb
          do k=1,nc
            if (mask(i,j,k) .gt. 0) then
              call errcon ('Oops - mask point left?!?!')
            end if
            if (mask(i,j,k) .lt. 0) mask(i,j,k) = 1
          end do
        end do
      end do
c
      return
      end
c
c
c
      subroutine dottum (mask,shadow,na,nb,nc,grid,origin,cell,
     +  a,sphrad,iunit,ierr)
c
      implicit none
c
      integer na,nb,nc,grid(3),origin(3),iunit,ierr
      integer mask(na,nb,nc),shadow(na,nb,nc)
c
      real cell(6),a(3,3),x(3),y(3),sphrad
c
      integer i,j,k,l,leng1,np,nd
c
      logical ldodot
c
      character line*80
c
code ...
c
 6000 format (' dot ',3f15.3)
 6100 format (' sphere_xyz ',4f15.3)
c
      ldodot = (sphrad .le. 0.0)
c
      if (ldodot) then
        write (iunit,'(a)',err=9999) 'begin mdot'
        write (iunit,'(a)',err=9999) 'colour magenta'
      else
        write (iunit,'(a)',err=9999) 'begin sdot'
        write (iunit,'(a)',err=9999) 'colour cyan'
      end if
c
      np = 0
      nd = 0
c
      do i=1,na
        x(1) = float(origin(1)+i-1) *
     +         (cell(1)/float(grid(1))) / cell(1)
        do j=1,nb
          x(2) = float(origin(2)+j-1) *
     +           (cell(2)/float(grid(2))) / cell(2)
          do k=1,nc
            if (mask(i,j,k) .gt. 0) then
              np = np + 1
              if (shadow(i,j,k) .gt. 0) then
                nd = nd + 1
                x(3) = float(origin(3)+k-1) *
     +                 (cell(3)/float(grid(3))) / cell(3)
                call mulmtx (a,x,y,3,3,1)
                if (ldodot) then
                  write (line,6000,err=9999) (y(l),l=1,3)
                else
                  write (line,6100,err=9999) (y(l),l=1,3),sphrad
                end if
                call pretty (line(2:))
                write (iunit,'(a)',err=9999) line(1:leng1(line))
              end if
            end if
          end do
        end do
      end do
c
      write (iunit,'(a)',err=9999) 'end_object'
      close (iunit)
      write (*,*) 'ODL file written'
c
      call jvalut (' Nr of mask points  :',1,np)
      call jvalut (' Nr of surface dots :',1,nd)
c
      ierr = 0
      return
c
 9999 continue
      call errcon ('While writing ODL file')
      close (iunit)
      ierr = -1
      return
c
      end
c
c
c
      subroutine combim (exoper,mask1,na1,nb1,nc1,
     +                   mask2,na2,nb2,nc2,lim1,lim2)
c
c ... combine masks
c
      implicit none
c
      integer na1,nb1,nc1,na2,nb2,nc2,lim1(2,3),lim2(2,3)
      integer i1,i2,i3,io1,io2,io3,j1,j2,n11,n22,n12,ntot,m1,m2
      integer mask1(na1,nb1,nc1),mask2(na2,nb2,nc2)
c
      real simil
c
      character exoper*(*),oper*(10)
c
code ...
c
      oper = exoper
c
      io1 = lim2(1,1) - lim1(1,1)
      io2 = lim2(1,2) - lim1(1,2)
      io3 = lim2(1,3) - lim1(1,3)
c
c ... SIMILARITY
c
      if (oper(1:3) .eq. 'SIM') then
        n11 = 0
        n12 = 0
        n22 = 0
        ntot = 0
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              ntot = ntot + 1
              m1 = mask1(i1,i2,i3)
              m2 = mask2(i1+io1,i2+io2,i3+io3)
              if (m1.eq.1 .and. m2 .eq. 1) then
                n12 = n12 + 1
                n22 = n22 + 1
                n11 = n11 + 1
              else if (m1.eq.1) then
                n11 = n11 + 1
              else if (m2.eq.1) then
                n22 = n22 + 1
              end if
            end do
          end do
        end do
c
        call jvalut (' Nr of points in common grid      :',1,ntot)
        call jvalut (' Nr of points set in mask 1 = N1  :',1,n11)
        call jvalut (' Nr of points set in mask 2 = N2  :',1,n22)
        call jvalut (' Nr of points set in both   = N12 :',1,n12)
c
        if (n11 .gt. 0 .and. n22 .gt. 0) then
          write (*,*)
          simil = float(n12)/(float(n11)+float(n22)-float(n12))
          call fvalut (
     +      ' Tanimoto index   = N12 / (N1+N2-N12) :',1,simil)
          simil = float(n12)/( sqrt(float(n11)) * sqrt(float(n22)) )
          call fvalut (
     +      ' Similarity index = N12 / SQRT(N1*N2) :',1,simil)
          simil = 2.0*float(n12)/(float(n11)+float(n22))
          call fvalut (
     +      ' Sim index 2      = 2*N12 / (N1+N2)   :',1,simil)
          simil = float(n12)/float(n11)
          call fvalut (
     +      ' Fraction common mask 1 = N12 / N1    :',1,simil)
          simil = float(n12)/float(n22)
          call fvalut (
     +      ' Fraction common mask 2 = N12 / N2    :',1,simil)
        end if
c
        write (*,*)
        call jvalut (
     +    ' AND    mask1 mask2 = N12         :',1,n12)
        call jvalut (
     +    ' OR     mask1 mask2 = N1+N2-N12   :',1,(n11+n22-n12))
        call jvalut (
     +    ' BUTNOT mask1 mask2 = N1-N12      :',1,(n11-n12))
        call jvalut (
     +    ' BUTNOT mask2 mask1 = N2-N12      :',1,(n22-n12))
        call jvalut (
     +    ' XOR    mask1 mask2 = N1+N2-2*N12 :',1,(n11+n22-2*n12))
c
c ... AND
c
      else if (oper(1:3) .eq. 'AND') then
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              if (mask1(i1,i2,i3).eq.1 .and.
     +            mask2(i1+io1,i2+io2,i3+io3) .eq. 1) then
                mask1(i1,i2,i3) = 1
              else
                mask1(i1,i2,i3) = 0
              end if
            end do
          end do
        end do
c
c ... OR
c
      else if (oper(1:2) .eq. 'OR') then
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              if (mask1(i1,i2,i3).eq.1 .or.
     +            mask2(i1+io1,i2+io2,i3+io3) .eq. 1) then
                mask1(i1,i2,i3) = 1
              else
                mask1(i1,i2,i3) = 0
              end if
            end do
          end do
        end do
c
c ... XOR
c
      else if (oper(1:3) .eq. 'XOR') then
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              j1 = mask1(i1,i2,i3)
              j2 = mask2(i1+io1,i2+io2,i3+io3)
              if ( (j1.eq.1 .and. j2.eq.0) .or.
     +             (j1.eq.0 .and. j2.eq.1) ) then
                mask1(i1,i2,i3) = 1
              else
                mask1(i1,i2,i3) = 0
              end if
            end do
          end do
        end do
c
c ... BUTNOT
c
      else if (oper(1:6) .eq. 'BUTNOT') then
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              if (mask1(i1,i2,i3).eq.1 .and.
     +            mask2(i1+io1,i2+io2,i3+io3) .eq. 0) then
                mask1(i1,i2,i3) = 1
              else
                mask1(i1,i2,i3) = 0
              end if
            end do
          end do
        end do
c
c ... INVALID OPERAND
c
      else
        call errcon ('Invalid logical operand :'//oper)
      end if
c
      return
      end
c
c
c
      subroutine pdbmsk (lun,cell,grid,pad,rad,rt,size,
     +                   mask,origin,extent,ierr)
c
c --- Given a PDB file, build up the AVGSYS mask
c --- Alwyn Jones 11-Nov-90, for bone_mask
c     3-Oct-91 for pdb version
c     930226 - gjk - turned into subroutine for use with MAMA
c                  - made alternative entry BONESM to prevent
c                    code duplication
c
      implicit none
c
      integer maxatm,maxopt
      parameter (maxatm=1500000,maxopt=10)
c
      character line*200,optpar(maxopt)*80
c
      integer i,j,k,l,ctmask,extent(3),grid(3),origin(3),pad(3)
      integer ct,i1,i2,i3,off,lun,ierr,size,ndum,inx,nxy,jerr
      integer nhet,nopt
      integer mask(size)
c
      real cell(6),g(3),distce,rad,rt(12),x(3),xp(3),xyzr(4,maxatm)
      real a(3,3),b(3,3),c(3,3),fake(6),ming(3),minr,maxg(3),maxr
      real x1(3),x2(3),gmin,xball,yball,zball,rball,minc(3),maxc(3)
c
      logical lcube,lbon
c
      data lcube /.false./, lbon /.false./
c
code ...
c
      ierr = -1
c
c --- Read in the atoms
c
      lcube = .false.
      ctmask = 0
      nhet = 0
200   continue
        read (lun,'(a)',err=6900,end=210) line
        if (line(1:6) .eq. 'HETATM') then
          nhet = nhet + 1
          goto 200
        end if
        if (line(1:6) .ne. 'ATOM  ') goto 200
        ctmask = ctmask+1
        if (ctmask .gt. maxatm) then
          call errcon ('Too many atoms - operation aborted')
          call jvalut (' Max nr of atoms :',1,maxatm)
          return
        end if
        read (line, '(30x,3f8.3)')  (xyzr(i,ctmask),i=1,3)
        xyzr(4,ctmask) = rad
      goto 200
c
210   continue
      call jvalut (' Number of atoms :',1,ctmask)
      if (nhet .gt. 0) then
        call jvalut (' WARNING - HETATM cards skipped :',1,nhet)
      end if
      close (lun)
      goto 6969
c
c --- Alternative entry CUBMSK : centre PLUS extent gives CUBE
c
      entry cubmsk (cell,grid,
     +      pad,rt,size,xball,yball,zball,rball,
     +      mask,origin,extent,ierr)
c
      lcube = .true.
      ctmask = 1
      xyzr(1,1) = xball
      xyzr(2,1) = yball
      xyzr(3,1) = zball
      xyzr(4,1) = rball
c
      call jvalut (' Number of atoms :',1,ctmask)
      goto 6969
c
c --- Alternative entry BALMSK : one centre PLUS a radius
c
      entry balmsk (cell,grid,
     +      pad,rt,size,xball,yball,zball,rball,
     +      mask,origin,extent,ierr)
c
      lcube = .false.
      ctmask = 1
      xyzr(1,1) = xball
      xyzr(2,1) = yball
      xyzr(3,1) = zball
      xyzr(4,1) = rball
c
      call jvalut (' Number of atoms :',1,ctmask)
      goto 6969
c
c --- Alternative entry BONESM : bones PLUS a radius
c
      entry bonesm (lun,cell,grid,pad,rad,rt,size,
     +              mask,origin,extent,ierr)
c
      lcube = .false.
      lbon = .true.
      ierr = -1
      jerr = 0
c
      read (lun,'(a)',err=6910) line
      call extrop (line,nopt,maxopt,optpar,jerr)
      if (jerr .ne. 0) goto 6910
      if (nopt .lt. 4) goto 6910
c
      call textut (' Data block name :',optpar(1))
      call textut (' Data block type :',optpar(2))
      call textut (' Nr of values    :',optpar(3))
      call textut (' Format          :',optpar(4))
c
      call str2i (optpar(3),ctmask,jerr)
      if (jerr .ne. 0) goto 6910
c
      ctmask = ctmask / 4
      if (ctmask .gt. maxatm) then
        call errcon ('Too many atoms - operation aborted')
        call jvalut (' Max nr of atoms :',1,maxatm)
        return
      end if
c
      read (lun,optpar(4),err=6910) ((xyzr(i,j),i=1,4), j=1,ctmask)
c
      call jvalut (' Number of atoms :',1,ctmask)
      close (lun)
      goto 6904
c
c --- Alternative entry BONESO : bones WITHOUT a radius
c
      entry boneso (lun,cell,grid,pad,rad,rt,size,
     +              mask,origin,extent,ierr)
c
      lcube = .false.
      ierr = -1
      jerr = 0
c
      read (lun,'(a)',err=6910) line
      call extrop (line,nopt,maxopt,optpar,jerr)
      if (jerr .ne. 0) goto 6910
      if (nopt .lt. 4) goto 6910
c
      call textut (' Data block name :',optpar(1))
      call textut (' Data block type :',optpar(2))
      call textut (' Nr of values    :',optpar(3))
      call textut (' Format          :',optpar(4))
c
      call str2i (optpar(3),ctmask,jerr)
      if (jerr .ne. 0) goto 6910
c
      ctmask = ctmask / 3
      if (ctmask .gt. maxatm) then
        call errcon ('Too many atoms - operation aborted')
        call jvalut (' Max nr of atoms :',1,maxatm)
        return
      end if
c
      read (lun,optpar(4),err=6910) ((xyzr(i,j),i=1,3), j=1,ctmask)
c
      do j=1,ctmask
        xyzr (4,j) = rad
      end do
c
      call jvalut (' Number of atoms :',1,ctmask)
      close (lun)
c
c ... get rid of any BONES at (1500,1500,1500)
c
 6904 continue
      write (*,*) 'Removing BONES at (1500,1500,1500) ...'
      do i=ctmask,1,-1
        if ( abs(xyzr(1,i)-1500.0) .gt. 0.1 .or.
     +       abs(xyzr(2,i)-1500.0) .gt. 0.1 .or.
     +       abs(xyzr(3,i)-1500.0) .gt. 0.1 ) then
          ctmask = i
          goto 6965
        end if
      end do
c
 6965 call jvalut (' Nr of atoms now :',1,ctmask)
c
      goto 6969
c
c ... START OF COMMON CODE FOR ALL ENTRIES
c
c --- Find the layout
c
 6969 continue
      if (ctmask .lt. 1) then
        call errcon ('No atoms')
        return
      end if
c
      do i=1,3
        fake(i) = 1.0
        fake(i+3) = cell(i+3)
      end do
      call orthog (fake, a, 1)
      call orthog (fake, c, 0)
c
c --- Transform the input coords
c
      do i=1,ctmask
        call vecrtv (xyzr(1,i), xyzr(1,i), 1, rt, rt(10))
      end do
c
c --- Check the coords to see the required origin and extent
c
      call orthog (cell, b, 1)
      do i=1,3
        ming(i) =  9999999.0
        maxg(i) = -9999999.0
        minc(i) =  9999999.0
        maxc(i) = -9999999.0
      end do
      minr =  9999999.0
      maxr = -9999999.0
c
      do 150 i=1,ctmask
        do 160 j=1,3
160       x1(j) = xyzr(j,i)
        call mulmtx (b, x1, x2, 3, 3, 1)
        do 170 j=1,3
c
          minc (j) = min (minc(j),xyzr(j,i))
          maxc (j) = max (maxc(j),xyzr(j,i))
c
          x1(j) = x2(j)*float(grid(j))
          ming (j) = min (ming(j),x1(j))
          maxg (j) = max (maxg(j),x1(j))
c
170     continue
        minr = min (minr,xyzr(4,i))
        maxr = max (maxr,xyzr(4,i))
150   continue
c
      call fvalut (' Lower bounds (coordinates) :',3,minc)
      call fvalut (' Upper bounds (coordinates) :',3,maxc)
      call fvalut (' Lower bounds (grid points) :',3,ming)
      call fvalut (' Upper bounds (grid points) :',3,maxg)
      call fvalut (' Smallest radius :',1,minr)
      call fvalut (' Largest  radius :',1,maxr)
c
c ... check if reasonable values for radii
c
      if (lbon .and. (minr .lt. 0.0 .or. maxr .gt. 9.99)) then
        call errcon ('Unreasonable radii !!! Use NEw OLd_bones ?')
        return
      end if
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
      do i=1,3
        extent(i) = extent(i) + 2*pad(i)
        origin(i) = origin(i) - pad(i)
      end do
c
      ndum = extent(1)*extent(2)*extent(3)
      if (ndum .gt. size) then
        call errcon ('Required mask is too big')
        call jvalut (' Padded origin :',3,origin)
        call jvalut (' Padded extent :',3,extent)
        call jvalut (' Maximum  nr of points :',1,size)
        call jvalut (' Required nr of points :',1,ndum)
        return
      end if
c
c --- Blank mask
c
      ct = 0
      do 110 k=1,extent(3)
        do 110 j=1,extent(2)
          do 110 i=1,extent(1)
            ct = ct+1
            mask(ct) = 0
110   continue
c
c --- Report
c
      call jvalut (' Mask origin :',3,origin)
      call jvalut (' Mask extent :',3,extent)
      call jvalut (' Grid points :',1,ct)
      call jvalut (' Mask grid   :',3,grid)
      call fvalut (' Mask cell   :',6,cell)
      call fvalut (' RT operator :',3,rt(1))
      call fvalut (' RT operator :',3,rt(4))
      call fvalut (' RT operator :',3,rt(7))
      call fvalut (' RT operator :',3,rt(10))
c
c --- Do the actual mask job
c
      nxy = extent(1)*extent(2)
      gmin = min (g(1),g(2),g(3))
c
      if (lcube) then
          do j=1,3
            x(j) = xyzr(j,1)
          end do
          call mulmtx (a, x, x1, 3, 3, 1)
          off = 1 + int (0.1 + (xyzr(4,1)/gmin))
          do 230 j=-off,off	
            do 230 k=-off,off	
              do 230 l=-off,off
                xp(1) = x1(1)+ l*g(1)
                xp(2) = x1(2)+ k*g(2)
                xp(3) = x1(3)+ j*g(3)
                call mulmtx (c, xp, x2, 3, 3, 1)
c
                if (abs(x(1)-x2(1)) .gt. xyzr(4,1)) goto 230
                if (abs(x(2)-x2(2)) .gt. xyzr(4,1)) goto 230
                if (abs(x(3)-x2(3)) .gt. xyzr(4,1)) goto 230
c
                i1 = nint(xp(1)/g(1))- origin(1) +1
                i2 = nint(xp(2)/g(2))- origin(2) +1
                i3 = nint(xp(3)/g(3))- origin(3) +1
c
                if (i1 .le. 0) goto 230
                if (i2 .le. 0) goto 230
                if (i3 .le. 0) goto 230
                if (i1 .gt. extent(1)) goto 230
                if (i2 .gt. extent(2)) goto 230
                if (i3 .gt. extent(3)) goto 230
c
                inx = (i3-1)*nxy + (i2-1)*extent(1) + i1
                mask(inx) = 1
230       continue
      else
        do 120 i=1,ctmask
          do j=1,3
            x(j) = xyzr(j,i)
          end do
          call mulmtx (a, x, x1, 3, 3, 1)
          off = 1 + int(0.1 + (xyzr(4,i)/gmin))
          do 130 j=-off,off	
            do 130 k=-off,off	
              do 130 l=-off,off
                xp(1) = x1(1)+ l*g(1)
                xp(2) = x1(2)+ k*g(2)
                xp(3) = x1(3)+ j*g(3)
                call mulmtx (c, xp, x2, 3, 3, 1)
                if (distce(x,x2) .gt. xyzr(4,i)) goto 130
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
130       continue
120     continue
      end if
c
c --- DONE
c
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
      subroutine dovrml (mask,shadow,na,nb,nc,grid,origin,cell,
     +  a,iunit,ierr)
c
      implicit none
c
      integer na,nb,nc,grid(3),origin(3),iunit,ierr
      integer mask(na,nb,nc),shadow(na,nb,nc)
c
      real cell(6),a(3,3),x(3),y(3)
c
      integer i,j,k,l,leng1,nat
c
      character line*80
c
code ...
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  PointSet { startIndex 0 '/
     +  '    numPoints ',i8,' } }')
c
      write (iunit,6000,err=9999)
c
      nat = 0
      do i=1,na
        x(1) = float(origin(1)+i-1) *
     +         (cell(1)/float(grid(1))) / cell(1)
        do j=1,nb
          x(2) = float(origin(2)+j-1) *
     +           (cell(2)/float(grid(2))) / cell(2)
          do k=1,nc
            if (mask(i,j,k) .gt. 0) then
              if (shadow(i,j,k) .gt. 0) then
                x(3) = float(origin(3)+k-1) *
     +                 (cell(3)/float(grid(3))) / cell(3)
                call mulmtx (a,x,y,3,3,1)
                write (line,6010,err=9999) (y(l),l=1,3)
                call pretty (line)
                write (iunit,'(a)',err=9999) line(1:leng1(line))
                nat = nat + 1
              end if
            end if
          end do
        end do
      end do
c
      write (iunit,6020,err=9999)
      write (iunit,6030,err=9999) nat
c
      call jvalut (' Nr of points written :',1,nat)
      ierr = 0
      return
c
 9999 continue
      call errcon ('While writing VRML file')
      ierr = -1
      return
c
      end
c
c
c
      subroutine mskset (maskb,extb1,extb2,extb3,i,j,k,ival)
c
      implicit none
c
      integer extb1,extb2,extb3,i,j,k,ival
      integer maskb(extb1,extb2,extb3)
c
code ...
c
      if (i .ge. 1 .and. i .le. extb1) then
        if (j .ge. 1 .and. j .le. extb2) then
          if (k .ge. 1 .and. k .le. extb3) then
            maskb (i,j,k) = ival
          end if
        end if
      end if
c
      return
      end
c
c
c
      subroutine cogmsk (mask,na,nb,nc,ns,xn,yn,zn)
c
c ... calculate centre-of-gravity of a mask
c
      implicit none
c
c ... MAXSIZ = maximum dimension of mask in any direction
c
      integer maxsiz
      parameter (maxsiz=10000)
c
      real*8 sx,sy,sz,ss
      real xn,yn,zn
c
      integer na,nb,nc,i1,i2,i3,i,ns
      integer nx(maxsiz),ny(maxsiz),nz(maxsiz)
      integer mask(na,nb,nc)
c
code ...
c
      do i=1,maxsiz
        nx (i) = 0
        ny (i) = 0
        nz (i) = 0
      end do
c
      ns = 0
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            if (mask(i1,i2,i3) .eq. 1) then
              ns = ns + 1
              nx(i1) = nx(i1) + 1
              ny(i2) = ny(i2) + 1
              nz(i3) = nz(i3) + 1
            end if
          end do
        end do
      end do
c
      ss = dble(ns)
      sx = 0.0D0
      sy = 0.0D0
      sz = 0.0D0
c
      do i=1,maxsiz
        if (nx(i).gt.0) sx = sx + (dble(nx(i))*dble(i)/ss)
        if (ny(i).gt.0) sy = sy + (dble(ny(i))*dble(i)/ss)
        if (nz(i).gt.0) sz = sz + (dble(nz(i))*dble(i)/ss)
      end do
c
      xn = sx
      yn = sy
      zn = sz
c
ccc      write (*,*) 'COG - NS = ',ns
ccc      write (*,*) 'COG - DP = ',sx,sy,sz
ccc      write (*,*) 'COG - SP = ',xn,yn,zn
c
      return
      end
