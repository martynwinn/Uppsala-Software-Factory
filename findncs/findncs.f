      program findncs
c
c ... FINDNCS - detect NCS from arbitrary set of points
c
c ... GJ Kleywegt @ 960821,22,23
c
c ... f77 -o FINDNCS findncs.f ../gklib/6d_kleylib; strip FINDNCS
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'FINDNCS', vers = '060503/0.4')
c
      integer maxatm
      parameter (maxatm=2000)
c
      real dismat(maxatm,maxatm),xyzrt(3,maxatm)
      real xyz(3,maxatm),test(3,maxatm),try(3,maxatm),newxyz(3),rtx(12)
      real dist,angle,tangle,cut,dummy,xyztol,rms3,rms,qq,dismin,distce
      real disint,bestrms,disira,rmstol,q1,q2
c
      integer na,nb,i,j,k,l,length,natmin,ii,jj,kk,ll,numncs,ierr,nnow
      integer itot,bestii,bestjj,bestkk,jtot,nfound,idum,mm,nnear,iunit
c
      logical ltest(maxatm),ltry(maxatm),ldone(maxatm),lnear(maxatm)
      logical linter,xinter
c
      character atmnam(maxatm)*15,line*80,nowres*10,file*128
c
code ...
c
cATOM    876  N   ARG   111      15.669  15.293  31.727  1.00  6.72      1CBS1098
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      call gkinit (prognm,vers)
c
      call jvalut (' Maximum number of atoms :',1,maxatm)
c
      iunit = 10
      linter = xinter()
c
      write (*,*)
      file = 'user.pdb'
      call textin (' Input PDB file ?',file)
      call textut (' Input PDB file :',file)
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file')
        goto 9999
      end if
c
      write (*,*)
      na = 0
      nowres = '#!$@%^&@^@%!'
c
c ... read ATOMs and HETATMs
c
  110 continue
      read (iunit,'(A)',end=99) line
c
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') goto 110
      na = na + 1
      if (na.gt.maxatm) call errstp (' Too many atoms !!!')
      atmnam(na) = line(13:26)
c
      read (line(31:54),*) (xyz(i,na),i=1,3)
      goto 110
c
c ... done reading
c
   99 continue
      close (iunit)
      write (*,*)
      call jvalut (' Nr of atoms read :',1,na)
      if (na .lt. 2) goto 999
c
      do i=1,na
        ldone (i) = .false.
      end do
c
c ... set parameters for now
c
      natmin = 3
      numncs = 2
      xyztol = 2.0
      rms3   = 2.0
      dismin = 1.0
      disint = 2.0
      disira = 100.0
      rmstol = 0.001
c
      write (*,*)
      call jvalin (' Min nr of atoms to match  ?',1,natmin)
      call jvalin (' Nr of RT-operators total  ?',1,numncs)
      call fvalin (' Max initial triple RMSD   ?',1,rms3)
      call fvalin (' Min initial triple RMSD   ?',1,rmstol)
      call fvalin (' Atom-matching tolerance   ?',1,xyztol)
      call fvalin (' Min dist RT-related atoms ?',1,dismin)
      call fvalin (' Min dist intra-mol atoms  ?',1,disint)
      call fvalin (' Max dist intra-mol atoms  ?',1,disira)
      write (*,*)
c
      write (*,*)
      call jvalut (' Min nr of atoms to match  :',1,natmin)
      call jvalut (' Nr of RT-operators total  :',1,numncs)
      call fvalut (' Max initial triple RMSD   :',1,rms3)
      call fvalut (' Min initial triple RMSD   :',1,rmstol)
      call fvalut (' Atom-matching tolerance   :',1,xyztol)
      call fvalut (' Min dist RT-related atoms :',1,dismin)
      call fvalut (' Min dist intra-mol atoms  :',1,disint)
      call fvalut (' Max dist intra-mol atoms  :',1,disira)
      write (*,*)
c
      itot = 0
      nfound = 1
c
c ... set up distance matrix
c
      do i=1,na-1
        do j=i+1,na
          qq = dist(i,j,xyz)
          dismat (i,j) = qq
          dismat (j,i) = qq
        end do
        dismat (i,i) = 0.0
      end do
      dismat (na,na) = 0.0
c
      if (na .le. 15) then
        write (*,6310) (i,i=1,na)
        do i=1,na
          write (*,6320) i,(dismat(j,i),j=1,na)
        end do
      end if
c
 6310 format (/4x,15i7)
 6320 format (1x,i3,15(1x,f6.1))
c
      do i=1,na-2
        if (ldone(i)) goto 70
        do j=i+1,na-1
          if (dismat(i,j) .lt. disint) goto 80
          if (dismat(i,j) .gt. disira) goto 80
          if (ldone(j)) goto 80
          do k=j+1,na
            if (ldone(k)) goto 90
            if (dismat(i,k) .lt. disint) goto 90
            if (dismat(j,k) .lt. disint) goto 90
            if (dismat(i,k) .gt. disira) goto 90
            if (dismat(j,k) .gt. disira) goto 90
c
            itot = itot + 1
            write (*,*)
            call jvalut (' Triple :',1,itot)
c
            write (*,*) '|',atmnam(i),'|',atmnam(j),'|',atmnam(k),'|'
            write (*,'(3f6.2)') dismat(i,j),dismat(i,k),dismat(j,k)
c
            do l=1,3
              test(l,1) = xyz(l,i)
              test(l,2) = xyz(l,j)
              test(l,3) = xyz(l,k)
            end do
c
            do ll=1,na
              ltest(ll) = .false.
              ltry (ll) = .false.
            end do
c
c ... set up a test array of atoms which are "near" I, J and K
c
            nnear = 0
            do ll=1,na
              lnear (ll) = .false.
              if (ldone(ll)) goto 600
              if (ll.eq.i.or.ll.eq.j.or.ll.eq.k) goto 600
c
ccc              print *,'TEST ',ll,dismat(i,ll),
ccc     +          dismat(j,ll),dismat(k,ll)
c
              if (dismat(i,ll) .gt. disira) goto 600
              if (dismat(j,ll) .gt. disira) goto 600
              if (dismat(k,ll) .gt. disira) goto 600
              if (dismat(i,ll) .lt. disint .and.
     +            dismat(j,ll) .lt. disint .and.
     +            dismat(k,ll) .lt. disint) goto 600
              lnear (ll) = .true.
              nnear = nnear + 1
  600         continue
            end do
c
c ... skip if no nearby atoms
c
            call jvalut (' Nr of nearby atoms :',1,nnear)
            if (nnear .lt. 4) goto 90
c
c ... test I,J,K
c
  111       continue
c
            bestrms = rms3
            bestii = -1
            bestjj = -1
            bestkk = -1
c
            jtot = 0
            do ii=1,na
              if (ii.eq.i .or. ii.eq.j .or. ii.eq.k) goto 12
              if (ldone(ii)) goto 12
              if (dismat(i,ii) .lt. dismin) goto 12
              do jj=1,na
                if (jj.eq.i .or. jj.eq.j .or. jj.eq.k) goto 11
                if (ldone(jj)) goto 11
                if (dismat(j,jj)  .lt. dismin) goto 11
                if (dismat(ii,jj) .lt. disint) goto 11
                if (dismat(ii,jj) .gt. disira) goto 11
                do kk=1,na
                  if (ii .eq. jj .or. ii .eq. kk) goto 10
                  if (jj .eq. kk) goto 10
                  if (kk.eq.i .or. kk.eq.j .or. kk.eq.k) goto 10
                  if (ldone(kk)) goto 10
                  if (dismat(k,kk)  .lt. dismin) goto 10
                  if (dismat(ii,kk) .lt. disint) goto 10
                  if (dismat(jj,kk) .lt. disint) goto 10
                  if (dismat(ii,kk) .gt. disira) goto 10
                  if (dismat(jj,kk) .gt. disira) goto 10
c
                  jtot = jtot + 1
c
                  do ll=1,3
                    try(ll,1) = xyz(ll,ii)
                    try(ll,2) = xyz(ll,jj)
                    try(ll,3) = xyz(ll,kk)
                  end do
c
                  call lsqgjk (try,test,3,rms,rtx,ierr)
c
                  if (rms .gt. rms3) goto 10
c
                  call vecrtv (xyz,xyzrt,na,rtx(1),rtx(10))
                  nnow = 3
                  q1 = xyztol*xyztol
                  do ll=1,na
                    if (.not. lnear(ll)) goto 400
                    do mm=1,na
                      if (mm.eq.i.or.mm.eq.j.or.mm.eq.k) goto 410
                      if (mm.eq.ii.or.mm.eq.jj.or.mm.eq.kk) goto 410
                      if (ldone(mm).or.mm.eq.ll) goto 410
                      if (dismat(ii,mm) .gt. disira) goto 410
                      if (dismat(jj,mm) .gt. disira) goto 410
                      if (dismat(kk,mm) .gt. disira) goto 410
c
                      q2 = (xyz(1,mm)-xyzrt(1,ll))**2
                      if (q2 .gt. xyztol) goto 410
                      q2 = q2 + (xyz(2,mm)-xyzrt(2,ll))**2
                      if (q2 .gt. xyztol) goto 410
                      q2 = q2 + (xyz(3,mm)-xyzrt(3,ll))**2
                      if (q2 .gt. xyztol) goto 410
                      nnow = nnow + 1
                      if (nnow .ge. (nnear/2)) goto 490
c
  410                 continue
                    end do
c
  400               continue
                  end do
  490             continue
c
ccc                  if (nnow .ge. 10) write (*,*) 
ccc     +              ' Trial, RMSD, #/# : ',jtot,rms,nnow,nnear
c
                  if (nnow .ge. (nnear/2)) then
                    bestrms = rms
                    bestii = ii
                    bestjj = jj
                    bestkk = kk
c
ccc                    write (*,*) 
ccc     +              ' Trial, RMSD, #/# : ',jtot,rms,nnow,nnear
c
                    goto 16
                  end if
c
   10             continue
                end do
   11           continue
              end do
   12         continue
            end do
c
   16       continue
            write (*,*)
            call jvalut (' Triple :',1,itot)
            call jvalut (' Trials :',1,jtot)
            if (bestii .gt. 0) then
              ii = bestii
              jj = bestjj
              kk = bestkk
              rms = bestrms
              goto 20
            end if
c
            call prompt (' Sorry - no triples found')
            goto 90
c
c ... try to extend this triple
c
   20       continue
            nnow = 3
            ltest (i) = .true.
            ltest (j) = .true.
            ltest (k) = .true.
            ltry (ii) = .true.
            ltry (jj) = .true.
            ltry (kk) = .true.
c
            do ll=1,3
              try(ll,1) = xyz(ll,ii)
              try(ll,2) = xyz(ll,jj)
              try(ll,3) = xyz(ll,kk)
            end do
c
            call lsqgjk (try,test,3,rms,rtx,ierr)
c
            call jvalut (' Best match for triple :',1,itot)
            call fvalut (' RMSD (A) :',1,rms)
            call fvalut (' Operator :',12,rtx)
            write (*,6000) atmnam(i),atmnam(ii)
            write (*,6000) atmnam(j),atmnam(jj)
            write (*,6000) atmnam(k),atmnam(kk)
c
c ... try adding an atom
c
            do ii=1,na
              if (ltest(ii) .or. ltry(ii) .or. ldone(ii)) goto 200
              call vecrtv (xyz(1,ii),newxyz,1,rtx(1),rtx(10))
              do jj=1,na
                if (ii .eq. jj) goto 210
                if (ltest(jj) .or. ltry(jj) .or. ldone(jj)) goto 210
c                if (distce(newxyz,xyz(1,jj)) .lt. dismin) goto 210
                qq = distce (xyz(1,jj),newxyz)
                if (qq .le. xyztol) then
                  write (*,6000) atmnam(ii),atmnam(jj),qq
                  ltry (jj)   = .true.
                  ltest  (ii) = .true.
                  nnow = nnow + 1
                  do ll=1,3
                    test (ll,nnow) = xyz(ll,ii)
                    try  (ll,nnow) = xyz(ll,jj)
                  end do
                  call lsqgjk (try,test,nnow,rms,rtx,ierr)
                  call fvalut (' RMSD now :',1,rms)
c
c                        if (rms .gt. rms3) then
c                          call prompt (' REJECT !!!')
c                          nnow = nnow - 1
c                          ltest (jj) = .false.
c                          ltry (ii) = .false.
c                        end if
c
                end if
  210           continue
              end do
  200         continue
            end do
c
            call jvalut (' Nr of matched atoms :',1,nnow)
            if (nnow .ge. natmin) then
              call prompt (' Found one !')
              call anancs (1,rtx,.true.,ierr)
              nfound = nfound + 1
              idum = 0
              do ll=1,na
                if (ltry(ll)) ldone(ll) = .true.
                if (ldone(ll)) idum = idum + 1
              end do
              idum = na - idum
              call jvalut (' Nr of atoms left :',1,idum)
              if (nfound .ge. numncs) goto 999
              if (idum .lt. 2*natmin) then
                call jvalut (' Min nr to match  :',1,2*natmin)
                call prompt (' Not enough atoms left for another RT')
                goto 999
              end if
              goto 111
            end if
c
   90       continue
          end do
   80     continue
        end do
   70   continue
      end do
c
 6000 format (' Match |',a,'| ... |',a,'| ',:,' Distance ',f6.2,' A')
c
  999 continue
      write (*,*)
      call jvalut (' Nr of operators found :',1,nfound-1)
c
 9999 continue
c
      call gkquit ()
c
      end
