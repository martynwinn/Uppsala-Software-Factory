      program ssencs
c
c ... SSENCS - try to find NCS operators from a set of SSEs
c              (secondary-structure elements) in DEJAVU format
c
c ... G J Kleywegt @ 010222/020214
c
c ... f77 -o SSENCS ssencs.f ../../gklib/6d_kleylib; strip SSENCS
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'SSENCS', vers = '060503/0.3')
c
      integer maxsse, maxncs, maxpos, maxsol
      parameter (maxsse = 1000, maxncs = 60, maxpos=maxsse/5)
      parameter (maxsol = 10*maxsse)
c
      real sseco1(3,maxsse),sseco2(3,maxsse),sselen(maxsse)
      real ssecoc(3,maxsse),ssevec(3,maxsse)
      real ssert1(3,maxsse),ssert2(3,maxsse),ssertc(3,maxsse)
      real dismat(maxsse,maxsse)
      real solrt(12,maxsol),solpnt(3,maxsol),solrms(maxsol)
      real three1(3,maxsse),three2(3,maxsse)
      real rtlsq(12),xyzdum(3),xtest(3)
      real xdum,distce,mislen,disnbr,dismin,rmsd,maxini,max3rd
      real rmsd2,soldis,xbest,maxrtd,dismax,x1,x2,y1,y2
c
      integer nphit(maxsse),phits(maxpos,maxsse)
      integer ssetyp(maxsse),sseres(maxsse),ssenbr(maxsse)
      integer solcnt(maxsol),solind(maxsol),solmat(maxsol)
      integer sortem(maxsol),sse1(maxsse),sse2(maxsse)
      integer iunit,junit,i,j,k,l,m,nncs,nsse,ierr,nmat,idum
      integer misres,ibest,jbest,nokay,nsol,lbest,nmatch,maxlis
      integer nold,niter
c
      logical totry(maxsse)
      logical linter,xinter
c
      character label(maxsse)*10,res1(maxsse)*10,resn(maxsse)*10
      character file*256,line*256
c
code ...
c
      call gkinit (prognm,vers)
c
      iunit = 10
      junit = 11
c
      linter = xinter()
c
      write (*,*)
      call jvalut (' Max nr of SSEs in total :',1,maxsse)
      call jvalut (' Max nr of NCS molecules :',1,maxncs)
c
      write (*,*)
      file = 'user.sse'
      call textin (' Input SSE file ?',file)
      call textut (' Input SSE file :',file)
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening SSE file')
        goto 9999
      end if
c
      nsse = 0
c
   10 continue
      read (iunit,'(a)',err=90,end=100) line
      call textut (' >',line)
      if (line(1:6) .ne. 'ALPHA ' .and.
     +    line(1:6) .ne. 'BETA  ') goto 10
c
      if (nsse .eq. maxsse) then
        call errcon (' Too many SSEs - rest skipped !')
        goto 100
      end if
c
      nsse = nsse + 1
      ssetyp(nsse) = 1
      if (line(1:6) .eq. 'BETA  ') ssetyp(nsse) = 2
c
      read (line(7:),*,err=90,end=90)
     +  label(nsse),res1(nsse),resn(nsse),sseres(nsse),
     +  (sseco1(i,nsse),i=1,3),(sseco2(i,nsse),i=1,3)
c
      sselen (nsse) = distce (sseco1(1,nsse),sseco2(1,nsse))
      ssevec (1,nsse) = sseco2 (1,nsse) - sseco1 (1,nsse)
      ssecoc (1,nsse) = sseco1 (1,nsse) + 0.5 * ssevec (1,nsse)
      ssevec (2,nsse) = sseco2 (2,nsse) - sseco1 (2,nsse)
      ssecoc (2,nsse) = sseco1 (2,nsse) + 0.5 * ssevec (2,nsse)
      ssevec (3,nsse) = sseco2 (3,nsse) - sseco1 (3,nsse)
      ssecoc (3,nsse) = sseco1 (3,nsse) + 0.5 * ssevec (3,nsse)
c
      goto 10
c
   90 continue
      call errcon ('While reading SSE file')
      close (iunit)
      goto 9999
c
  100 continue
      close (iunit)
      write (*,*)
      call jvalut (' Nr of SSEs read :',1,nsse)
c
      write (*,*)
      nncs = 2
      call jvalut (' Min nr of molecules to look for :',1,nncs)
      call jvalut (' Max nr of molecules to look for :',1,maxncs)
      call jvalin (' Nr of molecules to look for ?',1,nncs)
      nncs = max (2, min (nncs, maxncs))
      call jvalut (' Nr of molecules to look for :',1,nncs)
      xdum = float(nsse)/float(nncs)
      call fvalut (' Average nr of SSEs per molecule :',1,xdum)
      if ((2*nncs) .gt. nsse) then
        call errcon ('Too few SSEs for so many molecules')
        goto 9999
      end if
c
      write (*,*)
      nmat = max (3, min ( int(xdum), 5))
      call prompt (' A core is a set of SSEs that matches another')
      call prompt (' set after applying an RT operator. This cut-off')
      call prompt (' only affects the output of RT operators.')
      call jvalin (' Min nr of SSEs in core ?',1,nmat)
      call jvalut (' Min nr of SSEs in core :',1,nmat)
      idum = nncs * nmat
      call jvalut (' Total nr of SSEs in core :',1,idum)
      if ( idum .gt. nsse) then
        call errcon ('Too few SSEs for such a big core')
        goto 9999
      end if
c
      write (*,*)
      maxlis = 10*nncs
      call jvalin (' Max nr of solutions to print ?',1,maxlis)
      maxlis = max (nncs, maxlis)
      call jvalut (' Max nr of solutions to print :',1,maxlis)
c
      write (*,*)
      misres = 3
      call prompt (' Potentially matchable SSEs may contain')
      call prompt (' different numbers of residues.')
      call jvalin (' Mismatch nr of residues ?',1,misres)
      misres = max (0, misres)
      call jvalut (' Mismatch nr of residues :',1,misres)
c
      write (*,*)
      mislen = misres*3.0
      call prompt (' Potentially matchable SSEs may have')
      call prompt (' different lengths.')
      call fvalin (' Mismatch length (A) ?',1,mislen)
      mislen = max (0.01, mislen)
      call fvalut (' Mismatch length (A) :',1,mislen)
c
      write (*,*)
      disnbr = 8.0
      call prompt (' SSE pairs to test should be close in space')
      call prompt (' (centre-of-gravity distance).')
      call fvalin (' Max nbr distance (A) ?',1,disnbr)
      disnbr = max (0.01, disnbr)
      call fvalut (' Max nbr distance (A) :',1,disnbr)
c
      write (*,*)
      maxini = 3.0
      call prompt (' For two pairs of SSEs to be considered')
      call prompt (' matchable, the RMSD of the 4 end-points')
      call prompt (' should not be too high (all 4 combinations')
      call prompt (' of matching them will be tried).')
      call fvalin (' Max initial RMSD (A) ?',1,maxini)
      maxini = max (0.01, maxini)
      call fvalut (' Max initial RMSD (A) :',1,maxini)
c
      write (*,*)
      max3rd = 3.0
      call prompt (' To extend a matching pair into a triple,')
      call prompt (' their superimposed C-of-Gs should be close.')
      call fvalin (' Max deviation 3rd SSE (A) ?',1,max3rd)
      max3rd = max (0.01, max3rd)
      call fvalut (' Max deviation 3rd SSE (A) :',1,max3rd)
c
      write (*,*)
      soldis = 3.0
      call prompt (' Distance between RT(test_vector) to decide')
      call prompt (' if two operators are essentially identical.')
      call fvalin (' Max projection distance (A) ?',1,soldis)
      soldis = max (0.01, soldis)
      call fvalut (' Max projection distance (A) :',1,soldis)
c
      write (*,*)
      maxrtd = 3.0
      call prompt (' Distance cut-off to decide if SSEs obey an')
      call prompt (' RT operator in the evaluation step.')
      call fvalin (' Max RT(SSE) distance (A) ?',1,maxrtd)
      maxrtd = max (0.01, maxrtd)
      call fvalut (' Max RT(SSE) distance (A) :',1,maxrtd)
c
      write (*,*)
      xtest(1) = 100.0
      xtest(2) = 100.0
      xtest(3) = 100.0
      call prompt (' Point for testing equivalence of operators')
      call fvalin (' Test vector ?',3,xtest)
      call fvalut (' Test vector :',3,xtest)
c
c ... do the real work
c
      write (*,*)
      call prompt (' Looking for seed SSEs ...')
      idum = 0
      do i=1,nsse
        totry (i) = .false.
        nphit (i) = 0
        ssenbr (i) = -1
        dismin = disnbr
        dismat (i,i) = 0.0
        do j=1,nsse
          if (i .eq. j) goto 200
          xdum = distce(ssecoc(1,i),ssecoc(1,j))
          dismat (i,j) = xdum
          if (xdum .le. dismin) then
            ssenbr (i) = j
            dismin = xdum
          end if
          if (ssetyp(i) .ne. ssetyp(j)) goto 200
          if (abs(sseres(i)-sseres(j)) .gt. misres) goto 200
          if (abs(sselen(i)-sselen(j)) .gt. mislen) goto 200
          if (nphit(i) .ge. maxpos) goto 210
          nphit(i) = nphit(i) + 1
          phits (nphit(i),i) = j
  200     continue
        end do
  210   continue
c
c ... maybe be more lenient here to generate more operators ?
c
        if (nphit(i) .ge. nncs .and. ssenbr(i) .gt. 0) then
          totry(i) = .true.
          idum = idum + 1
        end if
      end do
      call jvalut (' Nr of seed SSEs to try :',1,idum)
      if (idum .lt. 1) then
        call errcon (' No seed SSEs found - sorry !')
        goto 9999
      end if
c
      write (*,*)
      call prompt (' Generating RT operators and SSE triples ...')
      nokay = 0
      nsol = 0
c
      do i=1,nsse
c
        if (.not. totry(i)) goto 400
c
        do j=1,nphit(i)
          k = phits (j,i)
          if (.not. totry(k)) goto 300
          if (ssenbr(k) .le. 0) goto 300
          if (ssenbr(i).eq.k) goto 300
          if (ssenbr(k).eq.i) goto 300
          if (ssetyp(ssenbr(i)) .ne. ssetyp(ssenbr(k))) goto 300
c
          call getfit (
     +      sseco1(1,i),sseco2(1,i),
     +      sseco1(1,ssenbr(i)),sseco2(1,ssenbr(i)),
     +      sseco1(1,k),sseco2(1,k),
     +      sseco1(1,ssenbr(k)),sseco2(1,ssenbr(k)),
     +      rtlsq,rmsd,idum)
          if (rmsd .gt. maxini) goto 300
c
c ... try to find one more nearby SSE that matches
c
          dismin = max3rd
          ibest = -1
          jbest = -1
          do l=1,nsse
            if (l .ne. i .and. l. ne. k .and.
     +          l .ne. ssenbr(i) .and. l .ne. ssenbr(k)) then
              if (dismat (k,l) .le. disnbr) then
                call vecrtv (ssecoc(1,l),xyzdum,1,rtlsq(1),rtlsq(10))
                do m=1,nsse
                  if (m .ne. i .and. m. ne. k .and. m. ne. l .and.
     +                m .ne. ssenbr(i) .and. m .ne. ssenbr(k)) then
                    xdum = distce(xyzdum,ssecoc(1,m))
                    if (xdum .lt. dismin) then
                      dismin = xdum
                      ibest = l
                      jbest = m
                    end if
                  end if
                end do
              end if
            end if
          end do
c
          if (ibest .le. 0 .or. jbest .le. 0) goto 300
c
          do l=1,3
            three1 (l,1) = sseco1(l,i)
            three1 (l,2) = sseco2(l,i)
            three1 (l,3) = sseco1(l,ssenbr(i))
            three1 (l,4) = sseco2(l,ssenbr(i))
            three1 (l,5) = sseco1(l,jbest)
            three1 (l,6) = sseco2(l,jbest)
            three2 (l,1) = sseco1(l,k)
            three2 (l,2) = sseco2(l,k)
            three2 (l,3) = sseco1(l,ssenbr(k))
            three2 (l,4) = sseco2(l,ssenbr(k))
            three2 (l,5) = sseco1(l,ibest)
            three2 (l,6) = sseco2(l,ibest)
          end do
c
          call fit3rd (three1,three2,sseco1(1,ibest),
     +      sseco2(1,ibest),rtlsq,rmsd2,idum)
          if (rmsd2 .gt. maxini) goto 300
c
          nokay = nokay + 1
c
          call vecrtv (xtest,xyzdum,1,rtlsq(1),rtlsq(10))
c
          if (nsol .eq. 0) then
            nsol = nsol + 1
            solind (nsol) = nsol
            solcnt (nsol) = 1
            do l=1,12
              solrt (l,nsol) = rtlsq (l)
            end do
            do l=1,3
              solpnt (l,nsol) = xyzdum (l)
            end do
c
c          print *
c          print *,i,ssenbr(i),jbest,' <-> ',k,ssenbr(k),ibest
c          print *,rmsd,dismin,rmsd2
c          write (*,'(3f15.6)') rtlsq
c          write (*,'(a,3f15.6)') 'TEST ',xyzdum
c
          else
            lbest = -1
            xbest = soldis
            do l=1,nsol
              xdum = distce(solpnt(1,l),xyzdum)
              if (xdum .le. xbest) then
                lbest = l
                xbest = xdum
              end if
            end do
            if (lbest .gt. 0) then
              solcnt (lbest) = solcnt(lbest) + 1
            else
c
c          print *
c          print *,i,ssenbr(i),jbest,' <-> ',k,ssenbr(k),ibest
c          print *,rmsd,dismin,rmsd2
c          write (*,'(3f15.6)') rtlsq
c          write (*,'(a,3f15.6)') 'TEST ',xyzdum
c
              if (nsol .eq. maxsol) then
                call errcon ('Too many solutions - skipped')
              else
                nsol = nsol + 1
                solind (nsol) = nsol
                solcnt (nsol) = 1
                do l=1,12
                  solrt (l,nsol) = rtlsq (l)
                end do
                do l=1,3
                  solpnt (l,nsol) = xyzdum (l)
                end do
              end if
            end if
          end if
c
  300     continue
        end do
c
  400   continue
      end do
c
      call jvalut (' Nr of operators generated :',1,nokay)
      call jvalut (' Nr of unique solutions    :',1,nsol)
      if (nsol .lt. 1) goto 9999
c
ccc      write (*,*)
ccc      call prompt (' Sorting RT operators ...')
ccc      call shell (solcnt,solind,nsol)
c
ccc      do i=nsol,1,-1
ccc        write (*,'(a,i5,a,i5,a,i5)')
ccc     +    ' sol # ',nsol-i+1,' = ',solind(i),
ccc     +    ' multiplicity = ',solcnt(i)
ccc      end do
c
c ... now: - apply each operator to all SSEs
c          - figure out how many SSEs are matched with an RT(SSE)
c          - use these to compute a better operator
c          - sort by nr of matched SSEs
c          - output if > NMAT
c
      write (*,*)
      call prompt (' Evaluating RT operators ...')
c
      do i=1,nsol
        solind(i) = i
        j = i
        nold = 0
        niter = 0
c
 2024   continue
        niter = niter + 1
c
        call vecrtv (ssecoc(1,1),ssertc(1,1),nsse,
     +               solrt(1,j),solrt(10,j))
        call vecrtv (sseco1(1,1),ssert1(1,1),nsse,
     +               solrt(1,j),solrt(10,j))
        call vecrtv (sseco2(1,1),ssert2(1,1),nsse,
     +               solrt(1,j),solrt(10,j))
        nmatch = 0
c
        do k=1,nsse
          dismax = maxrtd
          ibest = -1
          do l=1,nsse
            if (ssetyp(k) .eq. ssetyp(l)) then
              xdum = distce(ssecoc(1,k),ssertc(1,l))
              if (xdum .lt. dismax) then
                dismax = xdum
                ibest = l
              end if
            end if
          end do
          if (ibest .gt. 0) then
            x1 = distce(sseco1(1,k),ssert1(1,ibest))
            x2 = distce(sseco1(1,k),ssert2(1,ibest))
            y1 = distce(sseco2(1,k),ssert1(1,ibest))
            y2 = distce(sseco2(1,k),ssert2(1,ibest))
            if (x1 .le. maxrtd .and. y2 .le. maxrtd) then
              nmatch = nmatch + 2
              do l=1,3
                three1(l,nmatch-1) = sseco1(l,k)
                three2(l,nmatch-1) = sseco1(l,ibest)
                three1(l,nmatch) = sseco2(l,k)
                three2(l,nmatch) = sseco2(l,ibest)
              end do
              sse1 (nmatch/2) = k
              sse2 (nmatch/2) = ibest
            else if (x2 .le. maxrtd .and. y1 .le. maxrtd) then
              nmatch = nmatch + 2
              do l=1,3
                three1(l,nmatch-1) = sseco1(l,k)
                three2(l,nmatch-1) = sseco2(l,ibest)
                three1(l,nmatch) = sseco2(l,k)
                three2(l,nmatch) = sseco1(l,ibest)
              end do
              sse1 (nmatch/2) = k
              sse2 (nmatch/2) = ibest
            end if
          end if
        end do
c
        solmat (j) = nmatch/2
        solrms (j) = 999.99
ccc        if (2*solmat (j) .gt. 0) then
        if (2*solmat (j) .gt. 2) then
          call lsqgjk (three1,three2,nmatch,solrms(j),solrt(1,j),ierr)
        end if
        call vecrtv (xtest,solpnt(1,j),1,solrt(1,j),solrt(10,j))
c
c ... iterate to see if we can catch more SSEs
c
        if (nmatch .gt. nold) then
          nold = nmatch
          goto 2024
        end if
c
        write (*,6060) i,j,niter,solcnt(j),solmat(j),solrms(j)
        sortem (j) = solmat (j)
c
        if (solmat(j) .ge. nmat) then
          m = min (15,solmat(j))
          write (*,6070) 'SSE',(label(sse1(l)),l=1,m)
          write (*,6070) 'RT(SSE)',(label(sse2(l)),l=1,m)
          write (*,*)
        end if
c
      end do
c
 6060 format (' # ',i5,' Operator ',i5,' Iterations ',i3,
     +  ' Multiplicity ',i3,' Matched SSEs ',i5,' RMSD (A) ',f8.2)
 6070 format (1x,a8,1x,15a6)
c
      write (*,*)
      call prompt (' Sorting RT operators ...')
c
ccc      call ivalut (' solmat before :',5,solmat)
ccc      call ivalut (' solind before :',5,solind)
c
      call shell (sortem,solind,nsol)
c
ccc      call ivalut (' solmat after  :',5,solmat)
ccc      call ivalut (' solind after  :',5,solind)
c
      l = 0
      m = 0
      do i=nsol,1,-1
        j = solind(i)
ccc        write (*,6060) i,j,l,solcnt(j),solmat(j),solrms(j)
        if ( solmat(j) .lt. nmat ) goto 500
        l = l + 1
        if (l .gt. maxlis) goto 500
        write (*,*)
        write (*,'(a,i5,a,i5,a,i5,a,f8.3,a)')
     +    ' RT Sol # ',nsol-i+1,' = ',j,' Matched SSEs = ',
     +    solmat(j),' RMSD = ',solrms(j),' A'
        write (*,'(a)') '.LSQ_RT_SSENCS R 12 (3f15.8)'
        write (*,'(3f15.8)') (solrt(k,j),k=1,12)
        m = m + 1
      end do
  500 continue
c
      write (*,*)
      call jvalut (' Nr of RT operators listed :',1,m)
c
 9999 continue
c
      call gkquit ()
      end
c
c
c
      subroutine getfit (x1,x2,x3,x4,y1,y2,y3,y4,lsqrt,rmsd,ibest)
c
      implicit none
c
      real x1(3),x2(3),x3(3),x4(3),y1(3),y2(3),y3(3),y4(3)
      real lsqrt(12),rt(12,4)
      real rmsd,rms(4)
c
      integer ibest,i
c
code ...
c
      call fittem (x1,x2,x3,x4,y1,y2,y3,y4,rt(1,1),rms(1))
      ibest = 1
c
      call fittem (x1,x2,x3,x4,y2,y1,y3,y4,rt(1,2),rms(2))
      if (rms(2) .lt. rms(ibest)) ibest = 2
c
      call fittem (x1,x2,x3,x4,y1,y2,y4,y3,rt(1,3),rms(3))
      if (rms(3) .lt. rms(ibest)) ibest = 3
c
      call fittem (x1,x2,x3,x4,y2,y1,y4,y3,rt(1,4),rms(4))
      if (rms(4) .lt. rms(ibest)) ibest = 4
c
      rmsd = rms(ibest)
      do i=1,12
        lsqrt(i) = rt(i,ibest)
      end do
c
      if (ibest .eq. 2 .or. ibest .eq. 4) call swap3d (y1,y2)
      if (ibest .eq. 3 .or. ibest .eq. 4) call swap3d (y3,y4)
c
      return
      end
c
c
c
      subroutine fittem (x1,x2,x3,x4,y1,y2,y3,y4,lsqrt,rmsd)
c
      implicit none
c
      real xx(12),yy(12)
      real x1(3),x2(3),x3(3),x4(3),y1(3),y2(3),y3(3),y4(3)
      real lsqrt(12)
      real rmsd
c
      integer i,ierr
c
code ...
c
      do i=1,3
        xx(i) = x1(i)
        xx(i+3) = x2(i)
        xx(i+6) = x3(i)
        xx(i+9) = x4(i)
        yy(i) = y1(i)
        yy(i+3) = y2(i)
        yy(i+6) = y3(i)
        yy(i+9) = y4(i)
      end do
c
      call lsqgjk (xx, yy, 4, rmsd, lsqrt, ierr)
c
      return
      end
c
c
c
      subroutine swap3d (x,y)
c
      implicit none
c
      real x(3),y(3)
      real xdum
c
      integer i
c
code ...
c
      do i=1,3
        xdum  = x(i)
        x (i) = y(i)
        y (i) = xdum
      end do
c
      return
      end
c
c
c
      subroutine fit3rd (xx,yy,z1,z2,lsqrt,rmsd,ibest)
c
      implicit none
c
      real xx(3,6),yy(3,6),z1(3),z2(3)
      real lsqrt(12),rt(12,2)
      real rmsd,rms(2)
c
      integer ibest,i,ierr
c
code ...
c
      do i=1,3
        yy(i,5) = z1(i)
        yy(i,6) = z2(i)
      end do
c
      call lsqgjk (xx, yy, 6, rms(1), rt(1,1), ierr)
      ibest = 1
c
      do i=1,3
        yy(i,5) = z2(i)
        yy(i,6) = z1(i)
      end do
c
      call lsqgjk (xx, yy, 6, rms(2), rt(1,2), ierr)
      if (rms(2) .lt. rms(1)) ibest = 2
c
      rmsd = rms(ibest)
      do i=1,12
        lsqrt(i) = rt(i,ibest)
      end do
c
      if (ibest .eq. 1) call swap3d (yy(1,5),yy(1,6))
      if (ibest .eq. 2) call swap3d (z1,z2)
c
      return
      end

