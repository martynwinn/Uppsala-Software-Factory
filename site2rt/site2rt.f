      program site2rt
c
c ... SITE2RT - find RT operators between sites in different
c               crystal forms
c
c ... G J Kleywegt @ 971125
c
c ... f77 -o SITE2RT site2rt.f ../gklib/6d_kleylib; strip SITE2RT
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'SITE2RT', vers = '041001/0.4')
c
      integer maxatm
      parameter (maxatm = 1000)
c
      real dis1(maxatm,maxatm),dis2(maxatm,maxatm)
      real xyz1(3,maxatm),xyz2(3,maxatm),xyz2rt(3,maxatm)
      real dum1(3,maxatm),dum2(3,maxatm)
      real rtx(12)
      real dist,dtoler,p12,p13,p23,q12,q13,q23,rms,ptoler
      real xnear,q,distce,rbest
c
      integer nuse1(maxatm),nuse2(maxatm)
      integer nb1(maxatm),nb2(maxatm)
      integer iunit,junit,na1,na2,i,j,k,i2,j2,k2,m,ierr,msites
      integer nok,nold,lnear,l2,l,nhit,length,nbest,mbest
c
      logical luse1(maxatm),luse2(maxatm)
      logical ldone1(maxatm),ldone2(maxatm)
      logical xinter,linter,lex,lfull,lncs
c
      character atm1(maxatm)*15,atm2(maxatm)*15
      character line*256,file1*256,file2*256,answer*1
c
code ...
cATOM    876  N   ARG   111      15.669  15.293  31.727  1.00  6.72      1CBS1098
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      call gkinit (prognm,vers)
c
      iunit = 10
      junit = 11
c
      linter = xinter()
c
      dtoler = 3.0
      ptoler = 3.0
      msites = 5
      lex = .true.
      lfull = .false.
c
      call jvalut (' Max nr of atoms :',1,maxatm)
      write (*,*)
c
      file1 = 'm1.pdb'
      call textin (' First PDB file (smallest nr of atoms) ?',file1)
c
      file2 = 'm2.pdb'
      call textin (' Second PDB file (largest nr of atoms) ?',file2)
c
      call xopxoa (iunit,file1,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file 1')
        goto 9999
      end if
c
      call xopxoa (junit,file2,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file 2')
        goto 9999
      end if
c
      na1 = 0
c
c ... read ATOMs and HETATMs
c
  110 continue
      read (iunit,'(A)',end=99) line
c
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') goto 110
      na1 = na1 + 1
      if (na1 .gt. maxatm) then
        call errcon ('Too many atoms !!!  ')
        goto 9999
      end if
      atm1(na1) = line(13:26)
c
      read (line(31:54),*) (xyz1(i,na1),i=1,3)
      goto 110
c
c ... done reading
c
   99 continue
      write (*,*)
      call jvalut (' Nr of atoms read file 1 :',1,na1)
      if (na1 .lt. 3) then
        call errcon ('Fewer than 3 atoms in file 1')
        goto 9999
      end if
c
      na2 = 0
c
c ... read ATOMs and HETATMs
c
  210 continue
      read (junit,'(A)',end=199) line
c
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') goto 210
      na2 = na2 + 1
      if (na2 .gt. maxatm) then
        call errcon ('Too many atoms !!!  ')
        goto 9999
      end if
      atm2(na2) = line(13:26)
c
      read (line(31:54),*) (xyz2(i,na2),i=1,3)
      goto 210
c
c ... done reading
c
  199 continue
      write (*,*)
      call jvalut (' Nr of atoms read file 2 :',1,na2)
      if (na2 .lt. 3) then
        call errcon ('Fewer than 3 atoms in file 2')
        goto 9999
      end if
c
      close (iunit)
      close (junit)
c
      if (na2 .lt. na1) then
        call prompt (
     +    ' WARNING - file 2 has fewer atoms than file 1 !')
      end if
c
      msites = (2*min(na1,na2))/3
c
      if (min(na1,na2) .gt. 20) lex = .false.
c
      if (na1 .eq. na2) then
        lncs = .false.
        do i=1,na1
          do j=1,3
            if (abs(xyz1(j,i)-xyz2(j,i)) .gt. 0.01) goto 123
          end do
        end do
        lncs = .true.
        call prompt ('0Mol 1 == Mol 2 ==> Pure NCS search !!!')
  123   continue
      end if
c
c ... calc distance matrices
c
      do i=1,na1-1
        do j=i+1,na1
          dis1(i,j) = dist(i,j,xyz1)
          dis1(j,i) = dis1(i,j)
        end do
        dis1(i,i) = 0.0
      end do
      dis1(na1,na1) = 0.0
c
      if (na1 .le. 15) then
        write (*,6310) 1,(i,i=1,na1)
        do i=1,na1
          write (*,6320) i,(dis1(j,i),j=1,na1)
        end do
      end if
c
 6310 format (/' Distance matrix Mol ',i1/4x,15i7)
 6320 format (1x,i3,15(1x,f6.1))
c
      do i=1,na2-1
        do j=i+1,na2
          dis2(i,j) = dist(i,j,xyz2)
          dis2(j,i) = dis2(i,j)
        end do
        dis2(i,i) = 0.0
      end do
      dis2(na2,na2) = 0.0
c
      if (na2 .le. 15) then
        write (*,6310) 2,(i,i=1,na2)
        do i=1,na2
          write (*,6320) i,(dis2(j,i),j=1,na2)
        end do
      end if
c
 1234 continue
      write (*,*)
      nhit = 0
      nbest = -1
      rbest = 9999.99
      mbest = 0
c
      do i=1,max(na1,na2)
        ldone1(i) = .false.
        ldone2(i) = .false.
      end do
c
      call fvalin (' Triple distance tolerance  ?',1,dtoler)
c
      call fvalin (' Superpositioning tolerance ?',1,ptoler)
c
      call jvalin (' Min nr of matched sites    ?',1,msites)
c
      write (*,*)
      answer = 'N'
      if (lex) answer = 'Y'
      call textin (' Exhaustive search (Y/N)    ?',answer)
      call upcase (answer)
      lex = (answer .ne. 'N')
c
      write (*,*)
      answer = 'N'
      if (lfull) answer = 'Y'
      call textin (' Output ALL solutions (Y/N) ?',answer)
      call upcase (answer)
      lfull = (answer .ne. 'N')
c
      call prompt ('0Start search - may take some time !')
c
c ... find sets of three with similar distances
c
      do i=1,na1-2
ccc        if (ldone1(i)) goto 4000
        do i2=1,na2
          if (ldone2(i2)) goto 3900
          if (lncs .and. (i.eq.i2)) goto 3900
          do j=i+1,na1-1
ccc            if (ldone1(j)) goto 3000
ccc            if (i.eq.j) goto 3000
            if (lncs .and. (j.eq.i2)) goto 3000
            q12 = dis1(i,j)
            do j2=1,na2
              if (i2 .eq. j2) goto 2900
              if (ldone2(j2)) goto 2900
              if (lncs .and. (i.eq.j2)) goto 2900
              if (lncs .and. (j.eq.j2)) goto 2900
              p12 = dis2(i2,j2)
              if (abs(p12-q12) .le. dtoler) then
                do k=j+1,na1
ccc                  if (ldone1(k)) goto 2000
ccc                  if (i.eq.k) goto 2000
ccc                  if (j.eq.k) goto 2000
                  q13 = dis1(i,k)
                  q23 = dis1(j,k)
                  do k2=1,na2
                    if (i2 .eq. k2) goto 1900
                    if (j2 .eq. k2) goto 1900
                    if (ldone2(k2)) goto 1900
                    if (lncs .and. (i.eq.k2)) goto 1900
                    if (lncs .and. (j.eq.k2)) goto 1900
                    if (lncs .and. (k.eq.k2)) goto 1900
                    p13 = dis2(i2,k2)
                    p23 = dis2(j2,k2)
                    if (abs(p13-q13) .le. dtoler .and.
     +                  abs(p23-q23) .le. dtoler) then
c
c                      write (*,*)
c                      write (*,'(4a,3f8.2)') ' MOL 1 ',
c     +                  atm1(i),atm1(j),atm1(k),q12,q13,q23
c                      write (*,'(4a,3f8.2)') ' MOL 2 ',
c     +                  atm2(i2),atm2(j2),atm2(k2),p12,p13,p23
c
                      do m=1,3
                        dum1(m,1)=xyz1(m,i)
                        dum1(m,2)=xyz1(m,j)
                        dum1(m,3)=xyz1(m,k)
                        dum2(m,1)=xyz2(m,i2)
                        dum2(m,2)=xyz2(m,j2)
                        dum2(m,3)=xyz2(m,k2)
                      end do
c
                      call lsqgjk (dum1,dum2,3,rms,rtx,ierr)
c                      call fvalut (' RMSD :',1,rms)
c                      call fvalut (' RTX :',12,rtx)
c
                      do m=1,na1
                        luse1(m) = .false.
                      end do
                      luse1(i) = .true.
                      luse1(j) = .true.
                      luse1(k) = .true.
                      nuse1(1) = i
                      nuse1(2) = j
                      nuse1(3) = k
c
                      do m=1,na2
                        luse2(m) = .false.
                      end do
                      luse2(i2) = .true.
                      luse2(j2) = .true.
                      luse2(k2) = .true.
                      nuse2(1) = i2
                      nuse2(2) = j2
                      nuse2(3) = k2
c
                      call vecrtv (xyz2,xyz2rt,na2,rtx(1),rtx(10))
c
                      nok = 3
c
 5000                 continue
                      nold = nok
c
                      do l2=1,na2
                        if (.not. (luse2(l2).or.ldone2(l2))) then
                          xnear = 9999.99
                          lnear = -1
                          do l=1,na1
ccc                            if (.not. (luse1(l).or.ldone1(l))) then
                            if (.not. (luse1(l))) then
                              if (lncs .and. (l.eq.l2)) goto 900
                              if (lncs .and. luse1(l2)) goto 900
                              if (lncs .and. luse2(l)) goto 900
                              q = distce(xyz1(1,l),xyz2rt(1,l2))
                              if (q .lt. xnear) then
                                xnear = q
                                lnear = l
                              end if
  900                         continue
                            end if
                          end do
                          if (xnear .le. ptoler) then
                            nok = nok + 1
                            luse1(lnear) = .true.
                            luse2(l2) = .true.
                            nuse1(nok) = lnear
                            nuse2(nok) = l2
                            do m=1,3
                              dum1(m,nok)=xyz1(m,lnear)
                              dum2(m,nok)=xyz2(m,l2)
                            end do
                          end if
                        end if
                      end do
c
                      if (nok .gt. nold) then
                        call lsqgjk (dum1,dum2,nok,rms,rtx,ierr)
                        call vecrtv (xyz2,xyz2rt,na2,rtx(1),rtx(10))
                        goto 5000
                      end if
c
                      if (nok .ge. msites) then
                        nhit = nhit + 1
                        if (lfull) then
                          call prompt ('0-----------------------------')
                          write (*,*)
                          call jvalut (' Yeah ! Hit nr :',1,nhit)
                          call jvalut (' Matched sites :',1,nok)
                          call fvalut (' RMSD :',1,rms)
                        else
                          write (*,6100) nhit,nok,rms
                        end if
c
                        if (nok .ge. mbest) then
                          if ((nok.gt.mbest) .or. (rms.lt.rbest)) then
                            nbest = nhit
                            rbest = rms
                            mbest = nok
                            do m=1,nok
                              nb1(m)=nuse1(m)
                              nb2(m)=nuse2(m)
                            end do
                          end if
                        end if
c
                        if (lfull) then
                          write (*,6030) 'Mol 1','Mol 2'
                          do m=1,nok
                            q = distce(xyz1(1,nuse1(m)),
     +                                 xyz2rt(1,nuse2(m)))
                            write (*,6020) atm1(nuse1(m)),
     +                        atm2(nuse2(m)),q
                          end do
c                        write (*,'(6a)') ' MOL 1 ',
c     +                    (atm1(nuse1(m)),m=1,nok)
c                        write (*,'(6a)') ' MOL 2 ',
c     +                    (atm2(nuse2(m)),m=1,nok)
c
                        call prompt ('0 RT Mol 2 --> Mol 1 :')
                          write (line,*) '.LSQ_RT_SITE2RT_',nhit
                          call remspa (line)
                          write (*,6010) line(1:length(line)),rtx
ccc                        call fvalut (' RTX :',12,rtx)
c                        call anancs (1,rtx,.true.,ierr)
c
                          call lsqgjk (dum2,dum1,nok,rms,rtx,ierr)
                        call prompt ('0 RT Mol 1 --> Mol 2 :')
                          write (line,*) '.LSQ_RT_INVERSE_',nhit
                          call remspa (line)
                          write (*,6010) line(1:length(line)),rtx
c                        call anancs (1,rtx,.true.,ierr)
                        end if
c
                        if (.not. lex) then
                          do m=1,nok
ccc                          ldone1(nuse1(m)) = .true.
                            ldone2(nuse2(m)) = .true.
                          end do
                        end if
                      end if
c
                    end if
 1900               continue
                  end do
 2000             continue
                end do
              end if
 2900         continue
            end do
 3000       continue
          end do
 3900     continue
        end do
 4000   continue
      end do
c
 6010 format (/
     +  a,'   R   12   (3f16.8)'/4(3f16.8,:,/))
c 6010 format (/
c     +  ' ----- 8< ----- RT Operator in O format ----- 8< -----'/
c     +  a,'   R   12   (3f16.8)'/4(3f16.8/),
c     +  ' ----- 8< ----- RT Operator in O format ----- 8< -----')
 6020 format (1x,a15,' <--> ',a15,' = ',f8.3,' A')
 6030 format (/1x,a15,6x,a15,3x,'Distance'/
     +        1x,15('-'),6x,15('-'),3x,'--------')
 6100 format (' Hit nr ',i8,' *** Nr matched sites = ',i6,
     +  ' *** RMSD (A) = ',f8.3)
c
      call prompt ('0-----------------------------')
      write (*,*)
      call jvalut (' Nr of hits :',1,nhit)
c
      if (nhit .gt. 0 .and. nbest .gt. 0) then
        call prompt ('0-----------------------------')
        write (*,*)
        call jvalut (' Best hit was nr :',1,nbest)
        call jvalut (' Nr of matching sites :',1,mbest)
        do m=1,mbest
          do l=1,3
            dum1(l,m) = xyz1(l,nb1(m))
            dum2(l,m) = xyz2(l,nb2(m))
          end do
        end do
        call lsqgjk (dum1,dum2,mbest,rms,rtx,ierr)
        call vecrtv (xyz2,xyz2rt,na2,rtx(1),rtx(10))
        call fvalut (' With RMSD (A) :',1,rms)
c
        write (*,6030) 'Mol 1','Mol 2'
        do m=1,mbest
          q = distce(xyz1(1,nb1(m)),xyz2rt(1,nb2(m)))
          write (*,6020) atm1(nb1(m)),atm2(nb2(m)),q
        end do
c
        call prompt ('0 RT Mol 2 --> Mol 1 :')
        write (*,6010) '.LSQ_RT_SITE2RT_BEST',rtx
        call anancs (1,rtx,.true.,ierr)
c
        call lsqgjk (dum2,dum1,mbest,rms,rtx,ierr)
        call prompt ('0 RT Mol 1 --> Mol 2 :')
        write (*,6010) '.LSQ_RT_INVERSE_BEST',rtx
        call anancs (1,rtx,.true.,ierr)
c
      end if
c
      write (*,*)
      answer = 'N'
      if (nhit .le. 0) answer = 'Y'
      call textin (' Re-run with different parameters ?',answer)
      call upcase (answer)
      if (answer .ne. 'N') goto 1234
c
c ... all done
c
 9999 continue
      call gkquit
c
      end


