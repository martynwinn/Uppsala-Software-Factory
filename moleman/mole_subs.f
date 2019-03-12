c
c ... subroutines for MOLEMAN
c
c
c
c
      subroutine calcpp (natoms,allxyz,atmnam,iresid,
     +    nres,resptr,phipsi)
c
c ... CALCPP - calculate PHI and PSI angles
c
      implicit none
c
      include 'mole_dim.incl'
c
cc      integer maxres
cc      parameter (maxres=10000)
c
      real allxyz(3,*),phipsi(2,*),tangle,dist
c      real angle
c
      integer natoms,nres,iresid(*),resptr(*),nowres,ndum
      integer i,j,inx
c
      logical aminoa(maxres)
c
      character atmnam(*)*4
c
code ...
c
      nres = 0
      nowres = -1
      ndum = 0
c
      do i=1,natoms
c
c ... a new residue ?
c
        if (iresid(i) .ne. nowres) then
          if (nres .gt. 0) then
            aminoa (nres) = (ndum .eq. 3)
          end if
          nres = nres + 1
          resptr (nres) = i
          nowres = iresid (i)
          ndum = 0
        end if
c
c ... get N, CA, C coordinates
c
        if (atmnam(i) .eq. ' CA ') then
          ndum = ndum + 1
          inx = 3*(nres-1) + 2
          do j=1,3
            allxyz (j,inx) = allxyz (j,10+i)
          end do
        else if (atmnam(i) .eq. ' C  ') then
          ndum = ndum + 1
          inx = 3*(nres-1) + 3
          do j=1,3
            allxyz (j,inx) = allxyz (j,10+i)
          end do
        else if (atmnam(i) .eq. ' N  ') then
          ndum = ndum + 1
          inx = 3*(nres-1) + 1
          do j=1,3
            allxyz (j,inx) = allxyz (j,10+i)
          end do
        end if
      end do
c
      do i=1,nres
        phipsi (1,i) = 0.0
        phipsi (2,i) = 0.0
      end do
c
c ... calculate PHI an PSI
c
      do i=1,nres
c
c ... PHI
c
        if (i .gt. 1) then
          if (aminoa(i) .and. aminoa(i-1)) then
            inx = 3*(i-1)
            if (dist(inx,inx+1,allxyz) .le. 2.0) then
              phipsi (1,i) = tangle (inx+3,inx+2,inx+1,inx,allxyz)
            end if
c
c ... TAU
c
cc            phipsi(3,i) = angle(inx+3,inx+2,inx+1,allxyz)
          end if
        end if
c
c ... PSI
c
        if (i .lt. nres) then
          if (aminoa(i) .and. aminoa(i+1)) then
            inx = 3*(i-1) + 1
            if (dist(inx+3,inx+2,allxyz) .le. 2.0) then
              phipsi (2,i) = tangle (inx+3,inx+2,inx+1,inx,allxyz)
            end if
c
c ... TAU
c
cc            phipsi(3,i) = angle(inx+2,inx+1,inx,allxyz)
          end if
        end if
c
      end do
c
c ... discard waters etc which follow the amino acid sequence
c
      ndum = nres
      do i=nres,1,-1
        if (aminoa(i) .and. abs(phipsi(1,i)) .ge. 0.001 .and.
     +      abs(phipsi(2,i)) .ge. 0.001) then
          ndum = i
          goto 7403
        end if
      end do
 7403 continue
      nres = ndum
c
      return
      end
c
c
c
      subroutine calcca (natoms,allxyz,atmnam,iresid,
     +    nres,resptr,phipsi,iflag)
c
c ... CALCCA - calculate CA-CA-CA angles and CA-CA-CA-CA dihedrals
c
      implicit none
c
      include 'mole_dim.incl'
c
cc      integer maxres
cc      parameter (maxres=10000)
c
      real allxyz(3,*),phipsi(2,*),dist,angle,tangle
      real cadi(maxres),ave,sdv,xmin,xmax,xtot,xx,x,y
c
      integer natoms,nres,iresid(*),resptr(*),nowres,ndum
      integer i,j,inx,n,iflag,k,cnts(5)
c
      logical aminoa(maxres)
c
      character atmnam(*)*4
c
code ...
c
      nres = 0
      nowres = -1
      ndum = 0
      do i=1,natoms
c
c ... a new residue ?
c
        if (iresid(i) .ne. nowres) then
          if (nres .gt. 0) then
            aminoa (nres) = (ndum .eq. 1)
          end if
          nres = nres + 1
          resptr (nres) = i
          nowres = iresid (i)
          ndum = 0
        end if
c
c ... get CA coordinates
c
        if (atmnam(i) .eq. ' CA ') then
          ndum = ndum + 1
          inx = nres
          do j=1,3
            allxyz (j,inx) = allxyz (j,i)
          end do
        end if
c
      end do
c
      do i=1,nres
        phipsi (1,i) = -999.0
        phipsi (2,i) = -999.0
      end do
c
c ... calculate
c
      if (iflag .eq. 1) then
        n = 0
        do k=1,5
          cnts(k) = 0
        end do
        do k=2,nres-1
          if (.not. aminoa(k)) goto 7766
          if (.not. aminoa(k-1)) goto 7766
          xx = dist(k-1,k,allxyz)
          if (xx .gt. 5.0) goto 7766
          n = n + 1
          cadi (n) = xx
c
          if (xx .le. 2.80) then
            cnts(1)=cnts(1)+1
          else if (xx .le. 3.00) then
            cnts(2)=cnts(2)+1
          else if (xx .le. 3.70) then
            cnts(3)=cnts(3)+1
          else if (xx .le. 3.90) then
            cnts(4)=cnts(4)+1
          else
            cnts(5)=cnts(5)+1
          end if
c
 7766     continue
        end do
        if (n .lt. 2) then
          call errcon ('Fewer than 2 CA-CA distances')
        else
c
          call xstats (cadi,n,ave,sdv,xmin,xmax,xtot)
          write (*,6000) n,ave,sdv,xmin,xmax
c
          x = 100.0*float(cnts(1))/float(n)
          y = abs(x-0.008)/0.079
          write (*,6010) 'Short (<= 2.8 A) ',
     +      cnts(1),x,y
          if (y .gt. 3.0) call prompt (
     +      ' WARNING - 3 SIGMA deviant !')
c
          x = 100.0*float(cnts(2))/float(n)
          y = abs (x-0.240)/0.458
          write (*,6010) 'CIS   (<= 3.0 A) ',
     +      cnts(2),x,y
          if (y .gt. 3.0) call prompt (
     +      ' WARNING - 3 SIGMA deviant !')
c
          x = 100.0*float(cnts(3))/float(n)
          y = abs (x-1.517)/3.492
          write (*,6010) 'Poor  (<= 3.7 A) ',
     +      cnts(3),x,y
          if (y .gt. 3.0) call prompt (
     +      ' WARNING - 3 SIGMA deviant !')
c
          x = 100.0*float(cnts(4))/float(n)
          y = abs (x-96.818)/7.044
          write (*,6010) 'TRANS (<= 3.9 A) ',
     +      cnts(4),x,y
          if (y .gt. 3.0) call prompt (
     +      ' WARNING - 3 SIGMA deviant !')
c
          x = 100.0*float(cnts(5))/float(n)
          y = abs (x-1.416)/4.004
          write (*,6010) 'Long  (>  3.9 A) ',
     +      cnts(5),x,y
          if (y .gt. 3.0) call prompt (
     +      ' WARNING - 3 SIGMA deviant !')
c
          call r5caca (n,cnts(1))
c
        end if
        return
      end if
c
 6000 format (1x,i7,' CA-CA distances'/
     +  ' Average CA-CA distance = ',f8.3,' Sigma = ',f8.3/
     +  ' Minimum CA-CA distance = ',f8.3,' Maxim = ',f8.3)
 6010 format (1x,a,' CA-CA dists : ',i6,' res = ',f8.2,' % = ',
     +  f4.1,' SIGMA from mean')
c
      do i=2,nres-1
        if (.not. aminoa(i)) goto 7777
        if (.not. aminoa(i-1)) goto 7777
        if (dist(i-1,i,allxyz) .gt. 5.0) goto 7777
        if (.not. aminoa(i+1)) goto 7777
        if (dist(i,i+1,allxyz) .gt. 5.0) goto 7777
        inx = i
        phipsi (1,i) = angle (inx-1,inx,inx+1,allxyz)
c
        if (.not. aminoa(i+2)) goto 7777
        if (dist(i+1,i+2,allxyz) .gt. 5.0) goto 7777
        phipsi (2,i) = tangle (inx-1,inx,inx+1,inx+2,allxyz)
c
 7777   continue
      end do
c
c ... discard waters etc which follow the amino acid sequence
c
      ndum = nres
      do i=nres,1,-1
        if (aminoa(i) .and. phipsi(1,i) .ge. -900.0 .and.
     +      phipsi(2,i) .ge. -900.0) then
          ndum = i
          goto 7403
        end if
      end do
 7403 continue
      nres = ndum
c
      return
      end
c
c
c
      subroutine dorama (natoms,iresid,
     +    nres,resptr,phipsi,resnam)
c
c ... DORAMA - do Ramachandran plot
c
      implicit none
c
      real phipsi(2,*)
c
      integer natoms,nres,iresid(*),resptr(*),iunit,length
      integer ierr,lunhp,lenmol,i,leng1
c
      logical xinter
c
      character filer*80,molnam*10,resnam(*)*3
c
code ...
c
      iunit = 11
c
      write (*,*) 
      write (*,*) 'In the following, hit RETURN if you do NOT'
      write (*,*) 'want to produce the file the programs asks for'
      write (*,*) 
c
c ... ASCII file
c
      filer = ' '
      call textin (' Text file with PHI-PSIs ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          do i=1,nres
            if (abs(phipsi(1,i)) .ge. 0.001 .and.
     +          abs(phipsi(2,i)) .ge. 0.001) then
              write (filer,*) resnam(resptr(i)),iresid(resptr(i))
              call remspa(filer)
              write (iunit,'(a10,1x,3f10.2)')
     +          filer(1:leng1(filer)),phipsi(1,i),phipsi(2,i)
cccc     +          phipsi(3,i)
            end if
          end do
        end if
      end if
c
c ... O2D file
c
      filer = ' '
      call textin (' O2D Ramachandran plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          lunhp = iunit
          call stamp (filer)
          write (lunhp,8000) 'REMARK','Ramachandran plot'
          write (lunhp,8000) 'REMARK',(' '//filer(1:leng1(filer)))
          write (lunhp,8010) 'COLOUR',4
          write (lunhp,8010) 'MRKTYP',1
          write (lunhp,8020) 'MRKSIZ',2.0,2.0
          write (lunhp,8020) 'XYVIEW',-180.0,180.0,-180.0,180.0
          write (lunhp,8000) 'XLABEL','Phi'
          write (lunhp,8000) 'YLABEL','Psi'
          write (lunhp,8010) 'NPOINT',nres
          write (lunhp,8000) 'VALABS','(2f10.2,1x,a)'
          do i=1,nres
            write (filer,*) resnam(resptr(i)),iresid(resptr(i))
            call remspa(filer)
            write (iunit,'(2f10.2,1x,a)')
     +        phipsi(1,i),phipsi(2,i),filer(1:leng1(filer))
          end do
          write (lunhp,8000) 'END   '
        end if
      end if
c
 8000 format (a6,1x,a)
 8010 format (a6,1x,10i6)
 8020 format (a6,1x,6f10.3)
c
c ... O datablock file
c
      filer = ' '
      call textin (' O PHI-PSI datablock file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          lunhp = iunit
          molnam = 'M1'
          call textin (' Current molecule name in O ?',molnam)
          lenmol = length(molnam)
          write (lunhp,'(a,a,1x,a1,1x,i6,1x,a)')
     +      molnam(1:lenmol),'_RESIDUE_PHI','R',nres,
     +      '(7f10.3)'
          write (lunhp,'(7f10.3)') (phipsi(1,i),i=1,nres)
          write (lunhp,'(a,a,1x,a1,1x,i6,1x,a)')
     +      molnam(1:lenmol),'_RESIDUE_PSI','R',nres,
     +      '(7f10.3)'
          write (lunhp,'(7f10.3)') (phipsi(2,i),i=1,nres)
        end if
      end if
c
c ... Ramachandran HPGL file
c
      filer = ' '
      call textin (' HPGL Ramachandran plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          call hprama (iunit,nres,phipsi,resptr,resnam)
        end if
      end if
c
c ... Ramachandran PostScript file
c
      filer = ' '
      call textin (' PostScript Ramachandran plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call psrama (iunit,nres,phipsi,resptr,resnam,filer,'C')
      end if
c
c ... POLAR Ramachandran PostScript file
c
      filer = ' '
      call textin (' PostScript POLAR Ramachandran plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call psrama (iunit,nres,phipsi,resptr,resnam,filer,'P')
      end if
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine docara (natoms,iresid,
     +    nres,resptr,phipsi,resnam,caqual)
c
c ... DOCARA - do "CA"-Ramachandran plot
c
      implicit none
c
      real phipsi(2,*)
c
      integer caqual (0:60,-60:60)
      integer natoms,nres,iresid(*),resptr(*),iunit,length
      integer ierr,lunhp,i,i1,i2,leng1
c
      logical xinter
c
      character filer*80,resnam(*)*3,qual(0:3)*20
c
      data qual / 'DISALLOWED !','Generous','Additional','Core'/
c
code ...
c
      iunit = 11
c
      write (*,*) 
      write (*,*) 'In the following, hit RETURN if you do NOT'
      write (*,*) 'want to produce the file the programs asks for'
      write (*,*) 
c
c ... ASCII file
c
      filer = ' '
      call textin (' Text file with CA angles/dihedrals ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          do i=1,nres
            write (filer,*) resnam(resptr(i)),iresid(resptr(i))
            call remspa(filer)
            if (phipsi(1,i) .ge. -181.0 .and.
     +          phipsi(2,i) .ge. -181.0) then
              if (resnam(resptr(i)) .ne. 'GLY') then
                i1 = max (0, min (60, nint (phipsi(1,i)/3.0) ) )
                i2 = max (-60, min (60, nint (phipsi(2,i)/3.0) ) )
                write (iunit,'(a10,1x,2f10.2,2x,a)')
     +            filer(1:leng1(filer)),phipsi(1,i),phipsi(2,i),
     +            qual(caqual(i1,i2))
              else
                write (iunit,'(a10,1x,2f10.2,2x,a)')
     +            filer(1:leng1(filer)),phipsi(1,i),phipsi(2,i),
     +            'Glycine'
              end if
            else
              write (iunit,'(a10,2x,a)')
     +          filer(1:leng1(filer)),'Near terminus or break'
            end if
          end do
        end if
      end if
c
c ... O2D file
c
      filer = ' '
      call textin (' O2D CA angles/dihedrals plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          lunhp = iunit
          call stamp (filer)
          write (lunhp,8000) 'REMARK','CA pseudo-Ramachandran plot'
          write (lunhp,8000) 'REMARK',(' '//filer(1:leng1(filer)))
          write (lunhp,8010) 'COLOUR',4
          write (lunhp,8010) 'MRKTYP',1
          write (lunhp,8020) 'MRKSIZ',2.0,2.0
          write (lunhp,8020) 'XYVIEW',0.0,180.0,-180.0,180.0
          write (lunhp,8000) 'XLABEL','CA-CA-CA angle'
          write (lunhp,8000) 'YLABEL','CA-CA-CA-CA dihedral'
          write (lunhp,8010) 'NPOINT',nres
          write (lunhp,8000) 'VALABS','(2f10.2,1x,a)'
          do i=1,nres
            write (filer,*) resnam(resptr(i)),iresid(resptr(i))
            call remspa(filer)
            write (iunit,'(2f10.2,1x,a)')
     +        phipsi(1,i),phipsi(2,i),filer(1:leng1(filer))
          end do
          write (lunhp,8000) 'END   '
        end if
      end if
c
 8000 format (a6,1x,a)
 8010 format (a6,1x,10i6)
 8020 format (a6,1x,6f10.3)
c
c ... CA-Ramachandran PostScript file
c
      filer = ' '
      call textin (' PostScript CA angles/dihedrals plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call pscara (iunit,nres,phipsi,resptr,resnam,filer,
     +               caqual(0,-60))
      end if
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine hprama (lunhp,nres,phipsi,resptr,resnam)
c
      implicit none
c
      integer lunhp,nres,resptr(*),i,myleng,length
      integer plowl(2),pupr(2),psiz
c
      real phipsi(2,*),x1(2,42)
      real ulowl(2),usiz,phi,psi
c
      character resnam(*)*3,mytext*80
c
      data x1 /
     $	-176.0,-67.1 ,-176.0,-20.,
     $	-147.9,-35.3,-129.5,-31.5,
     $	-122.1,-19.5,-129.4,-13.2,
     $	-176.,40.4,-176.,180.,
     $	-53.,180.,-39.8,151.6,
     $	-38.7,96.9,-68.8,77.1,
     $	-86.,66.2,-84.8,38.1,
     $	-57.0,-1.1,-35.,-18.3,
     $	-35.5,-64.8,-176.,-67.,
     +	34.7,87.4,58.2,110.4,
     $	58.2,11.2,35.0,27.5,
     $	34.7,87.4,
     +	-175.4,-180.,-175.4,-169.5,
     $	-61.6,-169.5,-57.0,-180.,
     $	-175.4,-180.,
     +	-161.,146.,-161.,98.,
     $	-59.,98.,-59.,177.,
     $	-138.,177.,-138.,160.,
     $	-161.,146.,
     +	-161.6,-52.7,-61.9,-51.,
     $	-61.9,-34.1,-111.2,-34.1,
     $	-127.5,-43.3,-161.6,-43.3,
     $	-161.6,-52.7/
c
code ...
c
      plowl(1) = 0
      plowl(2) = 0
      pupr(1) = 15999
      pupr(2) = 11399
      psiz = 10000
      ulowl(1) = -300.
      ulowl(2) = -300.
      usiz = 700.
c
      call xhp_init(lunhp)
      call xhp_window(ulowl,plowl,pupr,usiz,psiz)
      call xhp_penclr(1)
c
c ---	first draw box
c
      call xhp_moveto(-180.,-180.)
      call xhp_lineto(+180.,-180.)
      call xhp_lineto( 180., 180.)
      call xhp_lineto(-180., 180.)
      call xhp_lineto(-180.,-180.)
      call xhp_moveto(   0.,-180.)
      call xhp_lineto(   0., 180.)
      call xhp_moveto(-180.,   0.)
      call xhp_lineto( 180.,   0.)
c
c ---	phi axis ticks
c
      do 300 i=1,11
        phi = -180.+i*30.
        call xhp_moveto(phi,-180.)
        call xhp_lineto(phi,-175.)
300   continue
c
c ---	psi axis ticks
c
      do 310 i=1,11
        psi = -180+i*30.
        call xhp_moveto(-180.,psi)
        call xhp_lineto(-175.,psi)
310   continue
c
c ---	text
c
      call xhp_txtsiz(4)
      call xhp_moveto(-190.,-190.)
      call xhp_text('-180',4)
      call xhp_moveto(-100.,-190.)
      call xhp_text('-90',3)
      call xhp_moveto(0.,-190.)
      call xhp_text('0',1)
      call xhp_moveto(80.,-190.)
      call xhp_text('90',2)
      call xhp_moveto(170.,-190.)
      call xhp_text('+180',4)
      call xhp_txtsiz(4)
      call xhp_moveto(-10.,-200.)
      call xhp_text('PHI',3)
c
      call xhp_txtsiz(4)
      call xhp_moveto(-210.,-180.)
      call xhp_text('-180',4)
      call xhp_moveto(-210.,-90.)
      call xhp_text('-90',3)
      call xhp_moveto(-190.,0.)
      call xhp_text('0',1)
      call xhp_moveto(-210.,90.)
      call xhp_text(' 90',3)
      call xhp_moveto(-210.,180.)
      call xhp_text('+180',4)
      call xhp_txtsiz(4)
      call xhp_moveto(-230.,0.)
      call xhp_text('PSI',3)
c
      call xhp_moveto (-180.0,-250.0)
      mytext = 'Ramachandran plot'
      call textin (' Label ?',mytext)
      myleng = length (mytext)
      if (myleng.gt.0) call xhp_text (mytext,myleng)
c
      call xhp_stamp (2,-180.0,-280.0,'MOLEMAN')
c
c ---	boundary
c
      call xhp_lintyp(1,2)
      call xhp_moveto(x1(1,1),x1(2,1))
      do 340 i=2,18
340      call xhp_lineto(x1(1,i),x1(2,i))
      call xhp_moveto(x1(1,19),x1(2,19))
      do 350 i=20,23
350      call xhp_lineto(x1(1,i),x1(2,i))
      call xhp_moveto(x1(1,24),x1(2,24))
      do 360 i=25,28
360      call xhp_lineto(x1(1,i),x1(2,i))
      call xhp_lintyp(1,1)
      call xhp_moveto(x1(1,29),x1(2,29))
      do 370 i=30,35
370      call xhp_lineto(x1(1,i),x1(2,i))
      call xhp_moveto(x1(1,36),x1(2,36))
      do 380 i=37,42
380      call xhp_lineto(x1(1,i),x1(2,i))
c
c ---	points now
c
      do 330 i=1,nres
c
        if (abs(phipsi(1,i)) .lt. 0.001 .or.
     +      abs(phipsi(2,i)) .lt. 0.001) goto 330
c
        if (resnam(resptr(i)) .eq. 'GLY' ) then
          call xhp_moveto(phipsi(1,i)-2.,phipsi(2,i)-2.)
          call xhp_lineto(phipsi(1,i)+2.,phipsi(2,i)-2.)
          call xhp_lineto(phipsi(1,i)+2.,phipsi(2,i)+2.)
          call xhp_lineto(phipsi(1,i)-2.,phipsi(2,i)+2.)
          call xhp_lineto(phipsi(1,i)-2.,phipsi(2,i)-2.)
        else
          call xhp_moveto(phipsi(1,i)-2.,phipsi(2,i)   )
          call xhp_lineto(phipsi(1,i)+2.,phipsi(2,i)   )
          call xhp_moveto(phipsi(1,i)   ,phipsi(2,i)-2.)
          call xhp_lineto(phipsi(1,i)   ,phipsi(2,i)+2.)
        end if
c
330   continue
c
      close (lunhp)
c
      return
      end
c
c
c
      subroutine psrama (lunhp,nres,phipsi,resptr,resnam,psfile,how)
c
      implicit none
c
      integer lunhp,nres,resptr(*),i,j
      integer coregn(37,37),nongly,noutlr,iphi,ipsi
c
      real phipsi(2,*),x1(2,42),phi,psi,dx,dy
c
      character resnam(*)*3,labx*40,laby*40,psfile*(*),how*1,line*128
c
      data x1 /
     $	-176.0,-67.1 ,-176.0,-20.,
     $	-147.9,-35.3,-129.5,-31.5,
     $	-122.1,-19.5,-129.4,-13.2,
     $	-176.,40.4,-176.,180.,
     $	-53.,180.,-39.8,151.6,
     $	-38.7,96.9,-68.8,77.1,
     $	-86.,66.2,-84.8,38.1,
     $	-57.0,-1.1,-35.,-18.3,
     $	-35.5,-64.8,-176.,-67.,
     +	34.7,87.4,58.2,110.4,
     $	58.2,11.2,35.0,27.5,
     $	34.7,87.4,
     +	-175.4,-180.,-175.4,-169.5,
     $	-61.6,-169.5,-57.0,-180.,
     $	-175.4,-180.,
     +	-161.,146.,-161.,98.,
     $	-59.,98.,-59.,177.,
     $	-138.,177.,-138.,160.,
     $	-161.,146.,
     +	-161.6,-52.7,-61.9,-51.,
     $	-61.9,-34.1,-111.2,-34.1,
     $	-127.5,-43.3,-161.6,-43.3,
     $	-161.6,-52.7/
c
code ...
c
      call xps_init ()
      call xps_open (lunhp,psfile,'MOLEMAN')
c
      call defcor (coregn)
c
c      dx = 0
c      do i=1,36
c        do j=1,36
c          if (coregn(i,j) .eq. 1) dx = dx + 1.0
c        end do
c      end do
c      dx = 100.0 * dx / (36.0 * 36.0)
c      print *,' CORE = ',dx,' %'
c
      if (how .eq. 'P') then
c
        call xps_polar (0.0,360.0)
c
        do i=1,42
          call fix360 (x1(1,i))
        end do
        do i=1,nres
          call fix360 (phipsi(1,i))
        end do
c
c        do i=1,36
c          do j=1,36
c            if (coregn(i,j) .eq. 1) then
c              phi = float(i-1)*10.0 - 180.0
c              psi = float(j-1)*10.0 - 180.0
c              if (phi .lt. 0.) phi = phi + 360.0
c              if (psi .lt. 0.) psi = psi + 360.0
c              call xps_grey_box (phi,phi+10.0,psi,psi+10.0)
c            end if
c          end do
c        end do
c
      else
c
        call xps_scale (-180.0,180.0,-180.0,180.0)
c
        do i=1,42
          call fixang (x1(1,i))
        end do
        do i=1,nres
          call fixang (phipsi(1,i))
        end do
c
        call xps_stroke ()
c
        do i=1,36
          do j=1,36
            if (coregn(i,j) .eq. 1) then
              phi = float(i-1)*10.0 - 180.0
              psi = float(j-1)*10.0 - 180.0
              call xps_grey_box (phi,phi+10.0,psi,psi+10.0)
            end if
          end do
        end do
c
        call xps_move (-180.,-180.)
        call xps_draw (+180.,-180.)
        call xps_draw ( 180., 180.)
        call xps_draw (-180., 180.)
        call xps_draw (-180.,-180.)
c
c      call xps_move (   0.,-180.)
c      call xps_draw (   0., 180.)
c      call xps_move (-180.,   0.)
c      call xps_draw ( 180.,   0.)
c
c ---	phi axis ticks
c
        do 300 i=1,11
          phi = -180.+i*30.
          call xps_move (phi,-180.)
          call xps_draw (phi,-175.)
          call xps_move (phi, 180.)
          call xps_draw (phi, 175.)
300     continue
c
c ---	psi axis ticks
c
        do 310 i=1,11
          psi = -180+i*30.
          call xps_move (-180.,psi)
          call xps_draw (-175.,psi)
          call xps_move ( 180.,psi)
          call xps_draw ( 175.,psi)
310     continue
c
      end if
c
c ---	text
c
      if (how .eq. 'P') then
        labx = 'PHI mapped to [0,360>'
        laby = 'PSI mapped to [0,360>'
      else
        labx = 'PHI mapped to [-180,180>'
        laby = 'PSI mapped to [-180,180>'
      end if
c      call textin (' Label for PHI axis ?',labx)
c      call textin (' Label for PSI axis ?',laby)
      call xps_label (labx,laby)
c
c ---	boundary
c
      if (how .eq. 'P') then
c
        call xps_dash ()
        call xps_dash ()
c
        call xps_move (x1(1,1),x1(2,1))
        do 340 i=2,18
340       call xps_draw (x1(1,i),x1(2,i))
        call xps_move (x1(1,19),x1(2,19))
        do 350 i=20,23
350       call xps_draw (x1(1,i),x1(2,i))
        call xps_move (x1(1,24),x1(2,24))
        do 360 i=25,28
360       call xps_draw (x1(1,i),x1(2,i))
        call xps_move (x1(1,29),x1(2,29))
        do 370 i=30,35
370       call xps_draw (x1(1,i),x1(2,i))
        call xps_move (x1(1,36),x1(2,36))
        do 380 i=37,42
380       call xps_draw (x1(1,i),x1(2,i))
c
        call xps_solid ()
c
      end if
c
      if (how .eq. 'C') then
        dx = 2.0
        dy = 2.0
      else
        dx = 2.0
        dy = 1.0
      end if
c
      nongly = 0
      noutlr = 0
c
c ---	points now
c
      do 330 i=1,nres
c
        if (abs(phipsi(1,i)) .lt. 0.001 .or.
     +      abs(phipsi(2,i)) .lt. 0.001) goto 330
c
        if (resnam(resptr(i)) .eq. 'GLY' ) then
          call xps_move (phipsi(1,i)-dx,phipsi(2,i)-dy)
          call xps_draw (phipsi(1,i)+dx,phipsi(2,i)-dy)
          call xps_draw (phipsi(1,i)+dx,phipsi(2,i)+dy)
          call xps_draw (phipsi(1,i)-dx,phipsi(2,i)+dy)
          call xps_draw (phipsi(1,i)-dx,phipsi(2,i)-dy)
        else
          call xps_move (phipsi(1,i)-dx,phipsi(2,i)   )
          call xps_draw (phipsi(1,i)+dx,phipsi(2,i)   )
          call xps_move (phipsi(1,i)   ,phipsi(2,i)-dy)
          call xps_draw (phipsi(1,i)   ,phipsi(2,i)+dy)
c
          if (how .ne. 'P') then
c
            iphi = int ( (180.0+phipsi(1,i)) / 10.0 ) + 1
            ipsi = int ( (180.0+phipsi(2,i)) / 10.0 ) + 1
c
            if (iphi .ge. 1 .and. iphi .le. 37 .and.
     +          ipsi .ge. 1 .and. ipsi .le. 37) then
              nongly = nongly + 1
              if (coregn(iphi,ipsi) .ne. 1) then
                noutlr = noutlr + 1
                call xps_move (phipsi(1,i)-dx,phipsi(2,i)-dy)
                call xps_draw (phipsi(1,i)+dx,phipsi(2,i)+dy)
                call xps_move (phipsi(1,i)+dx,phipsi(2,i)-dy)
                call xps_draw (phipsi(1,i)-dx,phipsi(2,i)+dy)
              end if
            end if
          end if
        end if
c
330   continue
c
      if (how .eq. 'P') goto 555
c
      call jvalut (' RAMA - Nr of non-Gly residues :',1,nongly)
      call jvalut (' RAMA - Nr of outliers (98%)   :',1,noutlr)
      if (nongly .gt. 0 .and. noutlr .gt. 0) then
        call fvalut (' RAMA - % Outliers             :',1,
     +      100.0*float(noutlr)/float(nongly))
      else
        call fvalut (' RAMA - % Outliers             :',1,0.0)
      end if
c
      write (line,*) ' Nr of non-Gly residues :',nongly
      call xps_legend (line)
      write (line,*) ' Nr of outliers (98%)   :',noutlr
      call xps_legend (line)
c
      if (nongly .gt. 0 .and. noutlr .gt. 0) then
        write (line,*) ' % Outliers             :',
     +    100.0*float(noutlr)/float(nongly)
      else
        write (line,*) ' % Outliers             :',0.0
      end if
      call xps_legend (line)
c
      line = ' Outliers are defined as residues lying outside the '//
     +  'core areas in which 98% of'
      call prompt (line)
      call xps_legend (line)
c
      line = ' all non-Gly residues were found to lie in a sample '//
     +  'of 403 PDB structures of'
      call prompt (line)
      call xps_legend (line)
c
      line = ' <= 95% homologous structures (Hobohm & Sander list) '//
     +  'at 2.0 A resolution or'
      call prompt (line)
      call xps_legend (line)
c
      line = ' better (GJK, October 1995, unpublished results; data '//
     +  'for 74,893 residues).'
      call prompt (line)
      call xps_legend (line)
c
      line = ' Outliers are shown as asterisks in the plot; core '//
     +  'areas are coloured grey.'
      call prompt (line)
      call xps_legend (line)
c
      line = ' The core areas occupy 19.7 % of the total plot area.'
      call prompt (line)
      call xps_legend (line)
c
  555 continue
      call xps_close ()
c
      close (lunhp)
c
      return
      end
c
c
c
      subroutine pscara (lunhp,nres,phipsi,resptr,resnam,psfile,
     +                   caqual)
c
      implicit none
c
      integer lunhp,nres,resptr(*),i
c
      real phipsi(2,*),x1(2,52),phi,psi,x,y
c
      integer caqual(0:60,-60:60),nn(0:3),ntot,i1,i2,ii
c
      character resnam(*)*3,labx*40,laby*40,psfile*(*),line*128
c
      data x1 /
     +  80,50,  89,76,  102,45,  90,28,  80,50,
     +  106,-180,  106,-135,  133,-135,  133,-180,  106,-180,
     +  110,175,  110,180,  121,180,  121,175,  110,175,
     +  79,51,  88,85,  90,124,  98,124,  98,72,  105,43,
     +          101,18,  90,18,  79,51,
     +  98,-180,  112,-93,  139,-93,  139,-180,  98,-180,
     +  103,165,  103,180,  131,180,  114,157,  103,165,
     +  84,180,  145,180,  138,140,  119,123,  119,57,
     +           141,45,  141,-6,  97,-65,  123,-50,  143,-50,
     +           149,-102,  149,-180,  83,-180,  84,-55,
     +           90,-25,  73,51,  82,104,  84,180/
c
code ...
c
      do i=0,3
        nn(i) = 0
      end do
      ntot = 0
c
      call xps_init ()
      call xps_open (lunhp,psfile,'MOLEMAN')
      call xps_scale (0.0,180.0,-180.0,180.0)
c
      call xps_move (   0.,-180.)
      call xps_draw ( 180.,-180.)
      call xps_draw ( 180., 180.)
      call xps_draw (   0., 180.)
      call xps_draw (   0.,-180.)
c
c      call xps_move (   0.,-180.)
c      call xps_draw (   0., 180.)
c      call xps_move (   0.,   0.)
c      call xps_draw ( 180.,   0.)
c
c ---	phi axis ticks
c
      do 300 i=1,10
        phi = (i-1)*20.
        call xps_move (phi,-180.)
        call xps_draw (phi,-175.)
        call xps_move (phi, 180.)
        call xps_draw (phi, 175.)
300   continue
c
c ---	psi axis ticks
c
      do 310 i=1,19
        psi = -180.+(i-1)*20.
        call xps_move (  0.,psi)
        call xps_draw (  5.,psi)
        call xps_move (180.,psi)
        call xps_draw (175.,psi)
310   continue
c
c ---	text
c
      labx = 'CA-CA-CA angle'
      laby = 'CA-CA-CA-CA dihedral'
c      call textin (' Label for X axis ?',labx)
c      call textin (' Label for Y axis ?',laby)
      call xps_label (labx,laby)
c
c ---	boundary
c
      call xps_colour (2)
      call xps_move (x1(1,1),x1(2,1))
      do i=2,5
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_move (x1(1,6),x1(2,6))
      do i=7,10
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_move (x1(1,11),x1(2,11))
      do i=12,15
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_colour (4)
      call xps_move (x1(1,16),x1(2,16))
      do i=17,24
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_move (x1(1,25),x1(2,25))
      do i=26,29
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_move (x1(1,30),x1(2,30))
      do i=31,34
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_colour (1)
      call xps_move (x1(1,35),x1(2,35))
      do i=36,52
        call xps_draw (x1(1,i),x1(2,i))
      end do
c
      call xps_colour (0)
      call xps_text (50.,  50.,14.,'A')
      call xps_text (50.,-170.,14.,'B')
      call xps_text (50., 170.,14.,'B')
c
c ---	points now
c
      do 330 i=1,nres
c
        if (phipsi(1,i) .le. -181.0 .or.
     +      phipsi(2,i) .le. -181.0) goto 330
c
        if (resnam(resptr(i)) .eq. 'GLY' ) then
          call xps_colour (0)
          call xps_move (phipsi(1,i)-1.,phipsi(2,i)-2.)
          call xps_draw (phipsi(1,i)+1.,phipsi(2,i)-2.)
          call xps_draw (phipsi(1,i)+1.,phipsi(2,i)+2.)
          call xps_draw (phipsi(1,i)-1.,phipsi(2,i)+2.)
          call xps_draw (phipsi(1,i)-1.,phipsi(2,i)-2.)
        else
          ntot = ntot + 1
          i1 = max (0, min (60, nint (phipsi(1,i)/3.0) ) )
          i2 = max (-60, min (60, nint (phipsi(2,i)/3.0) ) )
          ii = caqual(i1,i2)
          nn (ii) = nn (ii) + 1
          call xps_colour (0)
          if (ii .eq. 0) call xps_colour (1)
          call xps_move (phipsi(1,i)-1.,phipsi(2,i)   )
          call xps_draw (phipsi(1,i)+1.,phipsi(2,i)   )
          call xps_move (phipsi(1,i)   ,phipsi(2,i)-2.)
          call xps_draw (phipsi(1,i)   ,phipsi(2,i)+2.)
        end if
c
330   continue
c
      call jvalut (' Nr of non-GLY residues :',1,ntot)
      if (ntot .gt. 0) then
c
        x=100.0*float(nn(3))/float(ntot)
        y=abs(x-72.8)/8.9
        write (line,1010) 'CORE',nn(3),x,y
        call xps_legend (line)
        call prompt (line)
        if (y .gt. 3.0) then
          line = ' WARNING - 3 SIGMA deviant !'
ccc          call xps_legend (line)
          call prompt (line)
        end if
c
        write (line,1000) 'Additional',nn(2),
     +    100.0*float(nn(2))/float(ntot)
        call xps_legend (line)
        call prompt (line)
c
        write (line,1000) 'Generous',nn(1),
     +    100.0*float(nn(1))/float(ntot)
        call xps_legend (line)
        call prompt (line)
c
        x=100.0*float(nn(0))/float(ntot)
        y=abs(x-3.1)/2.2
        write (line,1010) 'DISALLOWED',nn(0),x,y
        call xps_legend (line)
        call prompt (line)
        if (y .gt. 3.0) then
          line = ' WARNING - 3 SIGMA deviant !'
ccc          call xps_legend (line)
          call prompt (line)
        end if
c
        line = ' For <= 2.0 A structures :'
        call xps_legend (line)
        call prompt (line)
        line = ' Core       -  7.1 % area - 72.8 % (8.9) residues'
        call xps_legend (line)
        call prompt (line)
        line = ' Additional -  5.2 % area - 12.8 % (4.0) residues'
        call xps_legend (line)
        call prompt (line)
        line = ' Generous   - 15.0 % area - 11.3 % (4.6) residues'
        call xps_legend (line)
        call prompt (line)
        line = ' Disallowed - 72.6 % area -  3.1 % (2.2) residues'
        call xps_legend (line)
        call prompt (line)
c
        call r5cara (ntot,nn(0)) 
c
      end if
c
 1000 format (1x,a10,' res : ',i6,' = ',f8.2,' %')
 1010 format (1x,a10,' res : ',i6,' = ',f8.2,' % = ',f4.1,
     +  ' SIGMA from mean')
c
      call xps_close ()
c
      close (lunhp)
c
      return
      end
c
c
c
      subroutine dobala (natoms,iresid,
     +    nres,resptr,phipsi,resnam)
c
c ... DOBALA - do Balasubramanian plot
c
      implicit none
c
      real phipsi(2,*)
c
      integer natoms,nres,iresid(*),resptr(*),iunit,length
      integer ierr,lunhp,lenmol,i,leng1
c
      logical xinter
c
      character filer*80,molnam*10,resnam(*)*3
c
code ...
c
      iunit = 11
c
      write (*,*) 
      write (*,*) 'In the following, hit RETURN if you do NOT'
      write (*,*) 'want to produce the file the programs asks for'
      write (*,*) 
c
c ... ASCII file
c
      filer = ' '
      call textin (' Text file with PHI-PSIs ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          do i=1,nres
            if (abs(phipsi(1,i)) .ge. 0.001 .or.
     +          abs(phipsi(2,i)) .ge. 0.001) then
              write (filer,*) resnam(resptr(i)),iresid(resptr(i))
              call remspa(filer)
              write (iunit,'(a10,1x,2f10.2)')
     +          filer(1:leng1(filer)),phipsi(1,i),phipsi(2,i)
            end if
          end do
        end if
      end if
c
c ... O datablock file
c
      filer = ' '
      call textin (' O PHI-PSI datablock file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          lunhp = iunit
          molnam = 'M1'
          call textin (' Current molecule name in O ?',molnam)
          lenmol = length(molnam)
          write (lunhp,'(a,a,1x,a1,1x,i6,1x,a)')
     +      molnam(1:lenmol),'_RESIDUE_PHI','R',nres,
     +      '(7f10.3)'
          write (lunhp,'(7f10.3)') (phipsi(1,i),i=1,nres)
          write (lunhp,'(a,a,1x,a1,1x,i6,1x,a)')
     +      molnam(1:lenmol),'_RESIDUE_PSI','R',nres,
     +      '(7f10.3)'
          write (lunhp,'(7f10.3)') (phipsi(2,i),i=1,nres)
        end if
      end if
c
c ... Balasubramanian HPGL file
c
      filer = ' '
      call textin (' HPGL Balasubramanian plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call xopxua (iunit,filer,xinter(),ierr)
        if (ierr .eq. 0) then
          call hpbala (iunit,nres,phipsi,resptr,resnam,iresid)
        end if
      end if
c
c ... Balasubramanian PostScript file
c
      filer = ' '
      call textin (' PostScript Balasubramanian plot file ?',filer)
      if (length(filer) .gt. 0) then
        close (iunit)
        call psbala (iunit,nres,phipsi,resptr,resnam,filer)
      end if
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine hpbala (lunhp,nres,phipsi,resptr,resnam,iresid)
c
      implicit none
c
      integer lunhp,nres,resptr(*),i,myleng,length,iresid(*)
      integer plowl(2),pupr(2),psiz,iax1,iax2
c
      real phipsi(2,*),ulowl(2),usiz,phi,psi,times,pos
c
      character resnam(*)*3,mytext*80
c
code ...
c
      plowl(1) = 0
      plowl(2) = 0
      pupr(1) = 15999
      pupr(2) = 11399
      psiz = 10000
      ulowl(1) = -120.
      ulowl(2) = -300.
      usiz = 700.
c
      call xhp_init(lunhp)
      call xhp_window(ulowl,plowl,pupr,usiz,psiz)
      call xhp_penclr(1)
c
      iax1 = 0
      iax2 = 20 * ((nres+19)/20)
      times = 2.0
c
c ---	first draw box
c
      call xhp_moveto(times*iax1,-180.)
      call xhp_lineto(times*iax2,-180.)
      call xhp_lineto(times*iax2, 180.)
      call xhp_lineto(times*iax1, 180.)
      call xhp_lineto(times*iax1,-180.)
c
c ---	phi/psi axis ticks
c
      do 420 i=1,11
        psi = -180+i*30.
        call xhp_moveto(times*iax1,psi)
        call xhp_lineto(times*iax1+5,psi)
420   continue
c
c ---	text
c
      call xhp_txtsiz(4)
      do 430 i=10,iax2,20
        write (mytext,*) iresid(resptr(i))
        call remspa (mytext)
        myleng = length (mytext)
        pos = i*times-10.
        call xhp_moveto(pos,-190.)
        call xhp_text (mytext,myleng)
        phi = i*times
        call xhp_moveto(phi,-180.)
        call xhp_lineto(phi,-175.)
430   continue
c
      pos = times*(iax1+iax2-20)*0.5
      call xhp_moveto(pos,-210.)
      call xhp_text('Residue',7)
c
      call xhp_txtsiz(4)
      pos = times*iax1
      call xhp_moveto(pos-30.,-180.)
      call xhp_text('-180',4)
      call xhp_moveto(pos-30.,-90.)
      call xhp_text('-90',3)
      call xhp_moveto(pos-10.,0.)
      call xhp_text('0',1)
      call xhp_moveto(pos-30.,90.)
      call xhp_text(' 90',3)
      call xhp_moveto(pos-30.,180.)
      call xhp_text('+180',4)
      call xhp_txtsiz(4)
      call xhp_moveto(pos-50.,0.)
      call xhp_text('PSI',3)
      call xhp_moveto(pos-50.,15.)
      call xhp_text('PHI',3)
c
      call xhp_moveto (pos,-250.0)
      mytext = 'Balasubramanian plot'
      call textin (' Label ?',mytext)
      myleng = length (mytext)
      if (myleng.gt.0) call xhp_text (mytext,myleng)
c
      call xhp_stamp (2,pos,-280.0,'MOLEMAN')
c
c ---	points now
c
      do 440 i=1,nres
        pos = i*times
c
        if (abs(phipsi(1,i)) .lt. 0.001 .or.
     +      abs(phipsi(2,i)) .lt. 0.001) goto 440
c
        if (resnam(resptr(i)) .eq. 'GLY' ) then
          call xhp_moveto(pos-2.,phipsi(1,i))
          call xhp_lineto(pos+2.,phipsi(1,i))
          call xhp_moveto(pos-2.,phipsi(1,i)-2.)
          call xhp_lineto(pos-2.,phipsi(1,i)+2.)
          call xhp_lineto(pos+2.,phipsi(1,i)+2.)
          call xhp_lineto(pos+2.,phipsi(1,i)-2.)
          call xhp_lineto(pos-2.,phipsi(1,i)-2.)
          call xhp_moveto(pos,phipsi(1,i))
          call xhp_lineto(pos,phipsi(2,i))
        else
          call xhp_moveto(pos-2.,phipsi(1,i))
          call xhp_lineto(pos+2.,phipsi(1,i))
          call xhp_moveto(pos,phipsi(1,i))
          call xhp_lineto(pos,phipsi(2,i))
        end if
440   continue
c
      close (lunhp)
c
      return
      end
c
c
c
      subroutine psbala (lunhp,nres,phipsi,resptr,resnam,psfile)
c
      implicit none
c
      integer lunhp,nres,resptr(*),i,iax1,iax2
c
      real phipsi(2,*),psi,times,pos
c
      character resnam(*)*3,labx*40,laby*40,psfile*(*)
c
code ...
c
      call xps_init ()
      call xps_open (lunhp,psfile,'MOLEMAN')
      call xps_scale (0.0,float(nres),-180.0,180.0)
c
      iax1 = 0
      iax2 = 20 * ((nres+19)/20)
      times = 1.0
c
c ---	first draw box
c
      call xps_move(times*iax1,-180.)
      call xps_draw(times*iax2,-180.)
      call xps_draw(times*iax2, 180.)
      call xps_draw(times*iax1, 180.)
      call xps_draw(times*iax1,-180.)
c
c ---	phi/psi axis ticks
c
      do 420 i=1,11
        psi = -180+i*30.
        call xps_move(times*iax1,psi)
        call xps_draw(times*iax1+5,psi)
420   continue
c
c ---	text
c
      labx = 'Residue'
      laby = 'PHI/PSI'
      call textin (' Label for residue axis ?',labx)
      call textin (' Label for PHI-PSI axis ?',laby)
      call xps_label (labx,laby)
c
c ---	points now
c
      do 440 i=1,nres
        pos = i*times
c
        if (abs(phipsi(1,i)) .lt. 0.001 .or.
     +      abs(phipsi(2,i)) .lt. 0.001) goto 440
c
        if (resnam(resptr(i)) .eq. 'GLY' ) then
          call xps_move(pos-times,phipsi(1,i))
          call xps_draw(pos+times,phipsi(1,i))
          call xps_move(pos-times,phipsi(1,i)-times)
          call xps_draw(pos-times,phipsi(1,i)+times)
          call xps_draw(pos+times,phipsi(1,i)+times)
          call xps_draw(pos+times,phipsi(1,i)-times)
          call xps_draw(pos-times,phipsi(1,i)-times)
          call xps_move(pos,phipsi(1,i))
          call xps_draw(pos,phipsi(2,i))
        else
          call xps_move(pos-times,phipsi(1,i))
          call xps_draw(pos+times,phipsi(1,i))
          call xps_move(pos,phipsi(1,i))
          call xps_draw(pos,phipsi(2,i))
        end if
440   continue
c
      call xps_close ()
c
      close (lunhp)
c
      return
      end
c
c
c
      subroutine mereul (ang,amat)
c
c.  the rotation matrices have been pinched from the MERLOT manual.
c.  This is used  once you solve the rotation function using
c.  merlot , you can rotate your molecule using merlot itself and 
c.  write out a merlot coordinate file - which is not a pdb file
c.  then write a code to rewrite it to the pdb format.  But this is 
c.  way of just extracting the angels is better as you can rotate
c.   other molecules too without bothering about converting them to merlot
c.. format too.
c...
c...........................................................................
c.. the matrix  is the transpose of the ccp4 euler angle matrix
c.. hence call that first and then do the transpose
c
      implicit NONE
c
      real ang(3),amat(3,3),bmat(3,3)
c
      integer i,j
c
code ...
c
      call ccpeul(ang,bmat)
      do i = 1,3
        do j=1,3
          amat(i,j)=bmat(j,i)
        enddo
      enddo
c
      return
      end
c
c
c
      subroutine merpol (ang,amat)
c..
c..  As in the earlier cases this matrix was also pinched from
c..  the manual.
c.
c.  rotation 1 (Phi ) about Z     -----   called F
c.. rotation 2 (Psi )  about the new Y  ------ called P
c.   rotation 3 (Kappa ) about the new Z ------ called K
cc.  rotation 4 (Psi -1) about the new Y
c...  Rotation 5 (Phi-1)  about the new Z
c.. the MATRIX
c..
c   +CK+(CF*CF*SP*SP(1-CK))  +CP*SK+(CF*SF*SP*SP(1-CK) -SF*SP*SK+(CF*CP*SP(1-CK))
C.. -CP*SK+(CF*SF*SP*SP(1-CK)  +CK+(SF*SF*SP*SP(1-CK)) +CF*SP*SK+(SF*CP*SP(1-CK))
C.. +SF*SP*SK+(CF*CP*SP(1-CK)) -CF*SP*SK+(SF*CP*SP(1-CK) +CK+(CP*CP(1-CK))
c..    now make the matrices and from the angles
cc.. remember ang(1) is  Phi ang(2) is Psi and ang(3) is Kappa
cc.
c
        real degtor
        parameter (degtor = 6.2831853071796/360.0)
c
        real ang(3),amat(3,3)
        real cf,cf2,sf,sf2,cp,cp2,sp,sp2,ck,sk,cck
c
code ...
c
	cf=cos(degtor*ang(1))
	cf2=cf*cf
	sf=sin(degtor*ang(1))
	sf2=sf*sf
	cp=cos(degtor*ang(2))
        cp2=cp*cp
	sp=sin(degtor*ang(2))
	sp2=sp*sp
	ck=cos(degtor*ang(3)) 
	sk=sin(degtor*ang(3))
	cck=1.0-ck
cc...
c..
c.. the bad code for making the matrix.  Note many things are calculated many times.
c.. a wonderful example of a bad code which works any way.
c...
	amat(1,1)=ck+(cf2*sp2*cck)
	amat(1,2)=cp*sk+(cf*sf*sp2*cck)
	amat(1,3)=-sf*sp*sk+(cf*sp*cp*cck)
	amat(2,1)=-cp*sk+(cf*sf*sp2*cck)
	amat(2,2)=ck+(sf2*sp2*cck)
	amat(2,3)=cf*sp*sk+(sf*cp*sp*cck)
	amat(3,1)=sf*sp*sk+(cf*cp*sp*cck)
	amat(3,2)=-cf*sp*sk+(sf*cp*sp*cck)
	amat(3,3)=ck+(cp2*cck)
c
	return
	end
c
c
c
      subroutine matana (r)
c
      implicit none
c
      real twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      real r(3,3),det,det3,theta,dtheta,v(3),vlen,p(3),e(3)
      real xsign,qqq
c
      integer i
c
code ...
c
      write (*,*)
      write (*,'(a,3f10.6)') ' Rot mat ',r(1,1),r(1,2),r(1,3)
      write (*,'(a,3f10.6)') ' Rot mat ',r(2,1),r(2,2),r(2,3)
      write (*,'(a,3f10.6)') ' Rot mat ',r(3,1),r(3,2),r(3,3)
c
      det = det3(r)
      write (*,'(a,f10.6)')
     +  ' Determinant of rotation matrix = ',det
c
c ... find POLAR angles from matrix
c
      qqq = max (-1.0, min (1.0, (r(1,1)+r(2,2)+r(3,3)-1.0)*0.5))
      theta = acos ( qqq )
      dtheta = theta*rtodeg
c
c ... changed 940316
c
      xsign = 1.0
      if (dtheta .gt. 180.0) then
        dtheta = 360.0 - dtheta
        xsign = -1.0
      end if
cc      if (dtheta .gt. 180.0) dtheta = 360.0 - dtheta
c
      write (*,'(a,f10.6)')
     +  ' Rotation angle                 = ',dtheta
c
      if (dtheta .eq. 0.0) then
        write (*,'(a)')
     +    ' ERROR - indeterminate direction cosines'
        goto 200
      else if (dtheta .eq. 180.0) then
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
cc        vlen = -0.5 / sin(theta)
cc        v(1) = r(2,3)-r(3,2)
cc        v(2) = r(1,3)-r(3,1)
cc        v(3) = r(2,1)-r(1,2)
        do i=1,3
          v(i) = v(i)*vlen
        end do
        write (*,'(a,3f10.4)')
     +    ' Direction cosines of rotation axis : ',
     +    (v(i),i=1,3)
      end if
c
c
c ... changed 940316
c
      det = sqrt( abs(1.0-v(3)*v(3)) )
      p(1) = atan2 (det,v(3)) * rtodeg
      p(2) = atan2 (v(2),v(1)) * rtodeg
      p(3) = dtheta
cc      p(3) = dtheta
cc      p(1) = acos (v(3)) * rtodeg
cc      p(2) = atan2 (v(2),v(1)) * rtodeg
c
      write (*,'(a,3f10.4)')
     +  ' ALMN Polar angles Omega/Phi/Chi    : ',(p(i),i=1,3)
c
c ... find EULER angles from matrix
c
  200 continue
c
      e(1) = atan2 (r(2,3),r(1,3)) * rtodeg
      qqq = max (-1.0, min (1.0, r(3,3)))
      e(2) = acos ( qqq ) * rtodeg
      e(3) = atan2 (r(3,2),-r(3,1)) * rtodeg
c
      write (*,'(a,3f10.4)')
     +  ' ALMN Euler angles Alpha/Beta/Gamma : ',(e(i),i=1,3)
c
      write (*,*)
c
      return
      end
c
c
c
      subroutine bondbs (natoms,bndcut,batom,b1,b2,xyz,atmnam)
c
      implicit none
c
      real batom(*),b1(*),b2(*),xyz(3,*)
      real bndcut,rmsd,shap,corr,q1,q2,q3,q4,q5,q6
      real dist
c
      integer natoms,nbo,i,j,k,nonh
c
      logical lhydro
c
      character atmnam(*)*4
c
code ...
c
      if (natoms .lt. 2) then
        call errcon ('No atoms')
        return
      end if
c
      nbo = 0
      nonh = 0
      do i=1,natoms-1
        if (lhydro(atmnam(i))) goto 6592
        nonh = nonh + 1
        k=min(natoms,i+50)
        do j=i+1,k
          if (.not. lhydro(atmnam(j))) then
            if (dist(i,j,xyz) .le. bndcut) then
              nbo = nbo + 1
              b1 (nbo) = batom(i)
              b2 (nbo) = batom(j)
            end if
          end if
        end do
 6592   continue
      end do
c
      if (.not. lhydro(atmnam(natoms))) nonh = nonh + 1
c
      if (nbo .le. 5) then
        rmsd = -1.0
        corr = -1.0
      else
        call xystat (b1,b2,nbo,rmsd,shap,corr,q1,q2,q3,q4,q5,q6)
      end if
c
      write (*,*)
      call jvalut (' Nr of atoms    :',1,natoms)
      call jvalut (' Nr of non-Hs   :',1,nonh)
      call jvalut (' Nr of bonds    :',1,nbo)
      call fvalut (' RMSD bonded Bs :',1,rmsd)
      call fvalut (' Corr. coeff.   :',1,corr)
      call prompt (' Hydrogen atoms have been ignored')
      call prompt (' Ignore the following line in case of grouped Bs')
c
      if (rmsd .le. 5.0) then
        call prompt (' Bonded Bs properly restrained')
      else if (rmsd .le. 10.0) then
        call prompt (' Use stronger restraints on bonded Bs')
      else
        call prompt (' Use restraints on bonded Bs !!!')
      end if
c
      call r5bond (nonh,bndcut,nbo,rmsd,corr)
c
      return
      end
c
c
c
      subroutine cpcoor (natoms,xf,yf,zf,xxyyzz)
c
      implicit none
c
      real xf(*),yf(*),zf(*),xxyyzz(3,*)
c
      integer natoms,i
c
code ...
c
      do i=1,natoms
        xxyyzz (1,i) = xf (i)
        xxyyzz (2,i) = yf (i)
        xxyyzz (3,i) = zf (i)
      end do
c
      return
      end
c
c
c
      subroutine nonbbs (natoms,bndcut,nobcut,batom,
     +  b1,b2,xyz,atmnam,inb,usem)
c
      implicit none
c
      real batom(*),b1(*),b2(*),xyz(3,*)
      real bndcut,nobcut,rmsd,shap,corr,q1,q2,q3,q4,q5,q6
      real dist,dx,dy,dz
c
      integer natoms,nbo,i,j,nonh,inb
c
      logical usem(*),lhydro
c
      character atmnam(*)*4
c
code ...
c
      if (natoms .lt. 2) then
        call errcon ('No atoms')
        return
      end if
c
      do i=1,natoms
        usem (i) = .false.
        if (lhydro(atmnam(i))) goto 69
        if (inb .eq. 2) then
          if (atmnam(i)(1:2) .eq. ' N') usem(i) = .true.
          if (atmnam(i)(1:2) .eq. ' O') usem(i) = .true.
        else if (inb .eq. 3) then
          if (atmnam(i)(1:2) .eq. ' N') goto 69
          if (atmnam(i)(1:2) .eq. ' O') goto 69
          usem (i) = .true.
        else if (inb .eq. 1) then
          usem (i) = .true.
        end if
   69   continue
      end do
c
      nbo = 0
      nonh = 0
c
      do i=1,natoms-1
        if (.not. usem(i)) goto 6592
        nonh = nonh + 1
        do j=i+1,natoms
          if (.not. usem(i)) goto 6594
c
          dx = abs(xyz(1,i)-xyz(1,j))
          if (dx .gt. nobcut) goto 6594
c
          dy = abs(xyz(2,i)-xyz(2,j))
          if (dy .gt. nobcut) goto 6594
c
          dz = abs(xyz(3,i)-xyz(3,j))
          if (dz .gt. nobcut) goto 6594
c
          dist = sqrt (dx*dx + dy*dy + dz*dz)
          if (dist .ge. bndcut .and. dist .le. nobcut) then
              nbo = nbo + 1
              b1 (nbo) = batom(i)
              b2 (nbo) = batom(j)
          end if
 6594     continue
        end do
 6592   continue
      end do
c
      if (usem(natoms)) nonh = nonh + 1
c
      if (nbo .le. 5) then
        rmsd = -1.0
        corr = -1.0
      else
        call xystat (b1,b2,nbo,rmsd,shap,corr,q1,q2,q3,q4,q5,q6)
      end if
c
      write (*,*)
      call jvalut (' Number of atoms           :',1,natoms)
      call jvalut (' Number of atoms used      :',1,nonh)
      call jvalut (' Nr of non-bonded contacts :',1,nbo)
      call fvalut (' RMSD non-bonded Bs        :',1,rmsd)
      call fvalut (' Correlation coefficient   :',1,corr)
      call prompt (' Hydrogen atoms have been ignored')
c
      if (nbo .le. 5) return
c
      if (inb .eq. 1) then
        call r5non1 (bndcut,nobcut,nonh,nbo,rmsd,corr)
      else if (inb .eq. 2) then
        call r5non2 (bndcut,nobcut,nonh,nbo,rmsd,corr)
      else if (inb .eq. 3) then
        call r5non3 (bndcut,nobcut,nonh,nbo,rmsd,corr)
      end if
c
      return
      end
c
c
c
      subroutine calchi (natoms,allxyz,atmnam,resnam,iresid)
c
c ... CALCHI - calculate CHI*-DIHEDRAL angles
c
      implicit none
c
      include 'mole_dim.incl'
c
cc      integer maxres
cc      parameter (maxres=10000)
c
      real allxyz(3,*),xyz(3,9)
c
      integer natoms,nres,iresid(*),ndum,nowres
      integer i,j,inx
c
      logical know(9)
c
      character atmnam(*)*4,resnam(*)*3
c
code ...
c
      nres = 0
      nowres = -1
      ndum = 0
c
      write (*,6000) 'Residue','CHI1','CHI2','CHI3','CHI4','CHI5'
 6000 format (1x,a10,5(1x,a10))
c
      do i=1,natoms
c
c ... a new residue ?
c
        if (iresid(i) .ne. nowres) then
          if (nres .gt. 0) then
            if (ndum .ge. 3)
     +        call dochis (xyz,know,resnam(i-1),iresid(i-1))
          end if
          ndum = 0
          nowres = iresid (i)
          nres = nres + 1
          do j=1,9
            know(j)=.false.
          end do
        end if
c
c ... get coordinates
c
        if (atmnam(i) .eq. ' CA ') then
          ndum = ndum + 1
          inx = 1
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' C  ') then
          ndum = ndum + 1
          inx = 2
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' N  ') then
          ndum = ndum + 1
          inx = 3
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' CB ') then
          ndum = ndum + 1
          inx = 4
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' CG ' .or.
     +           atmnam(i) .eq. ' CG1' .or.
     +           atmnam(i) .eq. ' SG ' .or.
     +           atmnam(i) .eq. ' OG ' .or.
     +           atmnam(i) .eq. ' OG1') then
          ndum = ndum + 1
          inx = 5
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' CD ' .or.
     +           atmnam(i) .eq. ' CD1' .or.
     +           atmnam(i) .eq. ' OD1' .or.
     +           atmnam(i) .eq. ' ND1' .or.
     +           atmnam(i) .eq. ' SD ') then
          ndum = ndum + 1
          inx = 6
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' CE ' .or.
     +           atmnam(i) .eq. ' OE1' .or.
     +           atmnam(i) .eq. ' NE ') then
          ndum = ndum + 1
          inx = 7
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' CZ ' .or.
     +           atmnam(i) .eq. ' NZ ') then
          ndum = ndum + 1
          inx = 8
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        else if (atmnam(i) .eq. ' NH1') then
          ndum = ndum + 1
          inx = 9
          do j=1,3
            xyz (j,inx) = allxyz (j,i)
          end do
          know (inx) = .true.
        end if
      end do
c
      if (ndum .ge. 3) then
        call dochis (xyz,know,resnam(natoms),iresid(natoms))
      end if
c
      return
      end
c
c
c
      subroutine dochis (xyz,know,rnam,ires)
c
      implicit none
c
      real xyz(3,9),xx,tangle
c
      integer ires,i,leng1
c
      logical know(9)
c
      character rnam*3,chis(5)*10,line*80
c
code ...
c
      do i=1,5
        chis (i) = ' '
      end do
c
      if (know(3).and.know(1).and.know(4).and.know(5)) then
        xx = tangle (3,1,4,5,xyz)
        write (chis(1),'(f10.2)') xx
      else
        goto 10
      end if
c
      if (know(1).and.know(4).and.know(5).and.know(6)) then
        xx = tangle (1,4,5,6,xyz)
        write (chis(2),'(f10.2)') xx
      else
        goto 10
      end if
c
      if (know(4).and.know(5).and.know(6).and.know(7)) then
        xx = tangle (4,5,6,7,xyz)
        write (chis(3),'(f10.2)') xx
      else
        goto 10
      end if
c
      if (know(5).and.know(6).and.know(7).and.know(8)) then
        xx = tangle (5,6,7,8,xyz)
        write (chis(4),'(f10.2)') xx
      else
        goto 10
      end if
c
      if (know(6).and.know(7).and.know(8).and.know(9)) then
        xx = tangle (6,7,8,9,xyz)
        write (chis(5),'(f10.2)') xx
      else
        goto 10
      end if
c
   10 continue
      write (line,6000) rnam,ires,(chis(i),i=1,5)
      write (*,'(a)') line(1:leng1(line))
c
 6000 format (1x,a3,1x,i6,5(1x,a10))
c
      return
      end
c
c
c
      subroutine smoobs (natoms,bndcut,batom,b1,b2,xyz,atmnam)
c
      implicit none
c
      real batom(*),b1(*),b2(*),xyz(3,*)
      real bndcut,dist
c
      integer natoms,nbo,i,j,k,l,nonh
c
      logical lhydro
c
      character atmnam(*)*4
c
code ...
c
      if (natoms .lt. 2) then
        call errcon ('No atoms')
        return
      end if
c
      call prompt (' Smoothing temperature factors ...')
      nbo = 0
      nonh = 0
      do i=1,natoms
        if (lhydro(atmnam(i))) goto 6592
        nonh = nonh + 1
        k=max(1,i-50)
        l=min(natoms,i+50)
        b1 (i) = 1.0
        b2 (i) = batom(i)
        do j=k,l
          if (.not. lhydro(atmnam(j))) then
            if (i.ne.j) then
              if (dist(i,j,xyz) .le. bndcut) then
                nbo = nbo + 1
                b1 (i) = b1 (i) + 1
                b2 (i) = b2 (i) + batom(j)
              end if
            end if
          end if
        end do
        batom (i) = b2(i)/b1(i)
 6592   continue
      end do
c
      write (*,*)
      call jvalut (' Nr of atoms    :',1,natoms)
      call jvalut (' Nr of non-Hs   :',1,nonh)
      call jvalut (' Nr of bonds    :',1,(nbo/2))
      call prompt (' Hydrogen atoms have been ignored')
c
      return
      end
c
c
c
      subroutine counte (natoms,atmnam)
c
      implicit none
c
      integer maxelm
      parameter (maxelm=250)
c
      integer counts(maxelm),natoms,nelm,i,j
c
      logical lhydro
c
      character names(maxelm)*2,atmnam(*)*4,line*80
c
code ...
c
      if (natoms .lt. 1) then
        call errcon ('No atoms')
        return
      end if
c
      do i=1,maxelm
        counts (i) = 0
        names (i) = '  '
      end do
      nelm = 6
      names (1) = ' C'
      names (2) = ' H'
      names (3) = ' N'
      names (4) = ' O'
      names (5) = ' S'
      names (6) = ' P'
c
      do i=1,natoms
        if (lhydro(atmnam(i))) then
          counts (2) = counts (2) + 1
          goto 9999
        end if
        do j=1,nelm
          if (atmnam(i)(1:2) .eq. names(j)) then
            counts (j) = counts (j) + 1
            goto 9999
          end if
        end do
        nelm = nelm + 1
        names (nelm) = atmnam(i)(1:2)
        counts (nelm) = 1
c
 9999   continue
      end do
c
      if (counts(1) .gt. 0) then
        call jvalut (' Nr of Carbon atoms      :',1,counts(1))
      else
        call prompt (' There are *NO* Carbon atoms')
      end if
c
      if (counts(2) .gt. 0) then
        call jvalut (' Nr of Hydrogen atoms    :',1,counts(2))
        call prompt (' (Includes any NMR "pseudo-atoms")')
      else
        call prompt (' There are *NO* Hydrogen atoms')
      end if
c
      if (counts(3) .gt. 0) then
        call jvalut (' Nr of Nitrogen atoms    :',1,counts(3))
      else
        call prompt (' There are *NO* Nitrogen atoms')
      end if
c
      if (counts(4) .gt. 0) then
        call jvalut (' Nr of Oxygen atoms      :',1,counts(4))
      else
        call prompt (' There are *NO* Oxygen atoms')
      end if
c
      if (counts(5) .gt. 0) then
        call jvalut (' Nr of Sulfur atoms      :',1,counts(5))
      else
        call prompt (' There are *NO* Sulfur atoms')
      end if
c
      if (counts(6) .gt. 0) then
        call jvalut (' Nr of Phosphorous atoms :',1,counts(6))
      else
        call prompt (' There are *NO* Phosphorous atoms')
      end if
c
      if (nelm .gt. 6) then
        call prompt (' Additional elements:')
        do i=7,nelm
          line = ' Nr of "'//names(i)//'" atoms        :'
          call jvalut (line,1,counts(i))
        end do
      end if
c
      call jvalut (' Total nr of atoms       :',1,natoms)
c
      return
      end
c
c
c
      subroutine buried (natoms,allxyz,atmnam,resnam,iresid)
c
c ... BURIED - find buried charges (protein only)
c
      implicit none
c
      include 'mole_dim.incl'
c
cc      integer maxres
cc      parameter (maxres=10000)
c
      real allxyz(3,*)
      real dist,cutdis,x,cc
c
      integer natoms,nres,iresid(*),nowres,ityp
      integer i,j,nnit,noxy,nwat,ncar,nsul,npos,nneg,ntot,noth
c
      logical lcheck,lwater,lhydro
c
      character atmnam(*)*4,resnam(*)*3
c
code ...
c
      cutdis = 3.5
      nres = 0
      nowres = -1
      ityp = 0
c
 6000 format (/' ==> ',a3,' - ',i5)
 6010 format ( ' Neighbours for ',a4,' ...')
 6020 format ( ' ... ',a3,' - ',i5,' - ',a4,' @ ',f6.2,' A')
 6100 format ( ' #C ',i2,' | #N ',i2,' | #O ',i2,' | #S ',i2,
     +  ' | #other ',i2)
 6110 format ( ' #water ',i2,' | #positive ',i2,' | #negative ',i2,
     +  ' | #total ',i2)
c
      do i=1,natoms
c
c ... a new residue ?
c
        if (iresid(i) .ne. nowres) then
c
          if (nres .gt. 0 .and. ityp .gt. 0) then
            write (*,6100) ncar,nnit,noxy,nsul,noth
            write (*,6110) nwat,npos,nneg,ntot
c
            if (ntot .le. 1 .or. (nwat+2) .ge. ntot .or.
     +          nwat .ge. 4) then
              call prompt (' > Appears to be solvent exposed')
            else
              call prompt (' > Appears not to be solvent exposed')
              if (ityp .eq. 1 .or. ityp .eq. 2) then
                if (npos .gt. 0)
     +            call prompt (' > Positive charge(s) nearby !!')
                if (nneg .le. 0)
     +            call prompt (' > No negative charge(s) nearby !!')
                if ((noxy+nsul) .le. 2)
     +            call prompt (' > Not many O/S nearby !!')
                if (ncar .ge. noxy)
     +            call prompt (' > Many carbons nearby !!')
              else
                if (nneg .gt. 0)
     +            call prompt (' > Negative charge(s) nearby !!')
                if (npos .le. 0)
     +            call prompt (' > No positive charge(s) nearby !!')
                if (nnit .le. 1)
     +            call prompt (' > Not many nitrogens nearby !!')
                if (npos .le. 0 .and. nnit .le. 1 .and. nwat .le. 1)
     +            call prompt (' > Not many waters nearby !!')
                if ((ncar+nsul) .ge. nnit)
     +            call prompt (' > Many C/S nearby !!')
              end if
            end if
          end if
c
          ityp = 0
          if (resnam(i) .eq. 'ARG') then
            ityp = 1
          else if (resnam(i) .eq. 'LYS') then
            ityp = 2
          else if (resnam(i) .eq. 'ASP') then
            ityp = 3
          else if (resnam(i) .eq. 'GLU') then
            ityp = 4
          end if
c
          npos = 0
          nneg = 0
          nnit = 0
          noxy = 0
          nsul = 0
          ncar = 0
          noth = 0
          ntot = 0
          nwat = 0
c
          if (ityp .ne. 0) write (*,6000) resnam(i),iresid(i)
c
          nowres = iresid (i)
          nres = nres + 1
        end if
c
c ... check atom
c
        lcheck = .false.
        if (ityp.eq.1) then
          if (atmnam(i) .eq. ' NE ' .or.
     +       atmnam(i) .eq. ' NH1' .or.
     +       atmnam(i) .eq. ' NH2') lcheck = .true.
        else if (ityp .eq. 2) then
          if (atmnam(i) .eq. ' NZ ') lcheck = .true.
        else if (ityp .eq. 3) then
          if (atmnam(i) .eq. ' OD1' .or.
     +       atmnam(i) .eq. ' OD2') lcheck = .true.
        else if (ityp .eq. 4) then
          if (atmnam(i) .eq. ' OE1' .or.
     +       atmnam(i) .eq. ' OE2') lcheck = .true.
        end if
c
        if (lcheck) then
          write (*,6010) atmnam(i)
          do j=1,natoms
            if (lhydro(atmnam(j))) goto 1619
            if (iresid(j) .eq. nowres) goto 1619
            x = dist (i,j,allxyz)
            if (x .gt. cutdis) goto 1619
c
            ntot = ntot + 1
            write (*,6020) resnam(j),iresid(j),atmnam(j),x
c
            if (atmnam(j)(1:2) .eq. ' N') then
              nnit = nnit + 1
            else if (atmnam(j)(1:2) .eq. ' O') then
              noxy = noxy + 1
              if (lwater(resnam(j))) nwat = nwat + 1
            else if (atmnam(j)(1:2) .eq. ' S') then
              nsul = nsul + 1
            else if (atmnam(j)(1:2) .eq. ' C') then
              ncar = ncar + 1
            else
              noth = noth + 1
            end if
c
            call charge(resnam(j),atmnam(j),cc)
            if (cc .gt. 0.001) then
              npos = npos + 1
            else if (cc .le. -0.001) then
              nneg = nneg + 1
            end if
c
 1619       continue
c
          end do
        end if
c
      end do
c
      write (*,*)
c
      return
      end
c
c
c
      subroutine disdis (nat,xyz,atmnam)
c
      implicit none
c
      integer nat,maxdis
      parameter (maxdis = 1000)
c
      real xyz(3,nat),dist,x,y,z
c
      integer discnt(0:maxdis),ndis,i,j,k,noh
c
      logical lhydro
c
      character atmnam(nat)*4
c
code ...
c
      do i=0,maxdis
        discnt(i) = 0
      end do
c
      call prompt (' Counting distances (Hs ignored) ...')
      ndis = 0
      noh = 0
      do i=1,nat-1
        if (.not. lhydro(atmnam(i))) then
          noh = noh + 1
          do j=i+1,nat
            if (.not. lhydro(atmnam(j))) then
              x = dist (i,j,xyz)
              ndis = ndis + 1
              k = max (0, min (maxdis, int (x/0.5)))
              discnt (k) = discnt (k) + 1
            end if
          end do
        end if
      end do
      if (.not. lhydro(atmnam(nat))) noh = noh + 1
c
      call jvalut (' Nr of non-hydrogen atoms     :',1,noh)
      call jvalut (' Nr of inter-atomic distances :',1,ndis)
      if (ndis .lt. 3) return
c
      x = 100.0 / float(ndis)
      k = 0
      write (*,6000) 'Lower(A)','Upper(A)','     Count','       %',
     +  '     Cumul','  Cumul%'
      do i=0,maxdis
        if (discnt(i) .gt. 0) then
          k = k + discnt(i)
          y = float(discnt(i)) * x
          z = float (k) * x
          write (*,6100) float(i)*0.5,float(i+1)*0.5,discnt(i),y,k,z
        end if
      end do
c
 6000 format (1x,a8,1x,a8,1x,a10,1x,a8,1x,a10,1x,a8)
 6100 format (1x,f8.2,1x,f8.2,1x,i10,1x,f8.3,1x,i10,1x,f8.3)
c
      return
      end
c
c
c
      subroutine charge (resnam,atmnam,cc)
c
      implicit none
c
      real cc
c
      character resnam*3,atmnam*4
c
code ...
c
      cc = 0.0
c
      if (resnam .eq. 'ARG') then
        if (atmnam .eq. ' NE ' .or. atmnam .eq. ' NH1' .or.
     +      atmnam .eq. ' NH2') cc = 0.33
        return
      else if (resnam .eq. 'LYS') then
        if (atmnam .eq. ' NZ ') cc = 1.0
        return
      else if (resnam .eq. 'ASP') then
        if (atmnam .eq. ' OD1' .or. atmnam .eq. ' OD2') cc = -1.0
        return
      else if (resnam .eq. 'GLU') then
        if (atmnam .eq. ' OE1' .or. atmnam .eq. ' OE2') cc = -1.0
        return
      end if
c
      if (atmnam(3:3) .eq. '+' .or. atmnam(4:4) .eq. '+') cc = 1.0
      if (atmnam(3:3) .eq. '-' .or. atmnam(4:4) .eq. '-') cc = -1.0
c
      return
      end
c
c
c
      subroutine chknom (natoms,allxyz,atmnam,resnam,iresid,lchang)
c
      implicit none
c
      include 'mole_dim.incl'
c
cc      integer natoms,maxres
cc      parameter (maxres=10000)
c
      integer natoms
c
      real allxyz(3,natoms)
c
      integer iresid(natoms),resptr(maxres),nchk(5),nerr(5)
      integer i,ires,j,oldres,nres
c
      logical lchang,lerror
c
      character resnam(natoms)*3,atmnam(natoms)*4
c
code ...
c
      oldres = -99999
      nres = 0
c
      call jvalut (' Nr of atoms :',1,natoms)
      do i=1,natoms
        ires = iresid (i)
        if (oldres .ne. ires) then
          nres = nres + 1
          oldres = ires
          resptr (nres) = i
        end if
      end do
      resptr (nres+1) = natoms + 1
c
      call jvalut (' Nr of residues :',1,nres)
      if (nres .lt. 1) return
c
      do i=1,5
        nchk (i) = 0
        nerr (i) = 0
      end do
c
 6000 format (' Checking ',a3,1x,i6,' ...')
 6100 format (' # of ',a3,' checked : ',i5,' # errors : ',i5)
c
      write (*,*)
      do i=1,nres
        j = resptr(i)
        if (resnam(j) .eq. 'PHE') then
          if (.not.lchang) write (*,6000) resnam(j),iresid(j)
          call chkres (resnam(j),j,resptr(i+1)-1,allxyz,atmnam,
     +      iresid(j),natoms,lchang,lerror)
          nchk(1) = nchk(1) + 1
          if (lerror) nerr(1) = nerr(1) + 1
        else if (resnam(j) .eq. 'TYR') then
          if (.not.lchang) write (*,6000) resnam(j),iresid(j)
          call chkres (resnam(j),j,resptr(i+1)-1,allxyz,atmnam,
     +      iresid(j),natoms,lchang,lerror)
          nchk(2) = nchk(2) + 1
          if (lerror) nerr(2) = nerr(2) + 1
        else if (resnam(j) .eq. 'ASP') then
          if (.not.lchang) write (*,6000) resnam(j),iresid(j)
          call chkres (resnam(j),j,resptr(i+1)-1,allxyz,atmnam,
     +      iresid(j),natoms,lchang,lerror)
          nchk(3) = nchk(3) + 1
          if (lerror) nerr(3) = nerr(3) + 1
        else if (resnam(j) .eq. 'GLU') then
          if (.not.lchang) write (*,6000) resnam(j),iresid(j)
          call chkres (resnam(j),j,resptr(i+1)-1,allxyz,atmnam,
     +      iresid(j),natoms,lchang,lerror)
          nchk(4) = nchk(4) + 1
          if (lerror) nerr(4) = nerr(4) + 1
        else if (resnam(j) .eq. 'ARG') then
          if (.not.lchang) write (*,6000) resnam(j),iresid(j)
          call chkres (resnam(j),j,resptr(i+1)-1,allxyz,atmnam,
     +      iresid(j),natoms,lchang,lerror)
          nchk(5) = nchk(5) + 1
          if (lerror) nerr(5) = nerr(5) + 1
        end if
      end do
c
      write (*,*)
      write (*,6100) 'PHE',nchk(1),nerr(1)
      write (*,6100) 'TYR',nchk(2),nerr(2)
      write (*,6100) 'ASP',nchk(3),nerr(3)
      write (*,6100) 'GLU',nchk(4),nerr(4)
      write (*,6100) 'ARG',nchk(5),nerr(5)
c
      i=nerr(1)+nerr(2)+nerr(3)+nerr(4)+nerr(5)
      if (i .gt. 0) then
        if (lchang) then
          call prompt (
     +      ' WARNING - any attached hydrogens NOT renamed')
        end if
      else
        call prompt (' No problem, mon !')
      end if
c
      return
      end
c
c
c
      subroutine chkres (resnam,i1,i2,allxyz,atmnam,id,natoms,
     +                   lchang,lerror)
c
      implicit none
c
      integer natoms
c
      real allxyz(3,natoms),tangle,x1,x2
c
      integer i1,i2,ica,icb,icg,icd,icd1,icd2,iod1,iod2,ine
      integer ice1,ice2,ioe1,ioe2,icz,inh1,inh2,i,id
c
      logical lchang,lerror
c
      character atmnam(natoms)*4,resnam*3
c
      data ica,icb,icg,icd,icd1,icd2,iod1,iod2,ine /9*-1/
      data ice1,ice2,ioe1,ioe2,icz,inh1,inh2 /7*-1/
c
code ...
c
      lerror = .true.
c
      do i=i1,i2
        if (atmnam(i) .eq. ' CA ') then
          ica = i
        else if (atmnam(i) .eq. ' CB ') then
          icb = i
        else if (atmnam(i) .eq. ' CG ') then
          icg = i
        else if (atmnam(i) .eq. ' CD ') then
          icd = i
        else if (atmnam(i) .eq. ' CD1') then
          icd1 = i
        else if (atmnam(i) .eq. ' CD2') then
          icd2 = i
        else if (atmnam(i) .eq. ' OD1') then
          iod1 = i
        else if (atmnam(i) .eq. ' OD2') then
          iod2 = i
        else if (atmnam(i) .eq. ' NE ') then
          ine = i
        else if (atmnam(i) .eq. ' CE1') then
          ice1 = i
        else if (atmnam(i) .eq. ' CE2') then
          ice2 = i
        else if (atmnam(i) .eq. ' OE1') then
          ioe1 = i
        else if (atmnam(i) .eq. ' OE2') then
          ioe2 = i
        else if (atmnam(i) .eq. ' CZ ') then
          icz = i
        else if (atmnam(i) .eq. ' NH1') then
          inh1 = i
        else if (atmnam(i) .eq. ' NH2') then
          inh2 = i
        end if
      end do
c
 6000 format (' Error in ',a3,1x,i6,' ...')
c
      if (resnam .eq. 'PHE' .or. resnam .eq. 'TYR') then
        if (min(ica,icb,icg,icd1,icd2) .le. 0) then
          call errcon ('Missing atom(s) !')
          return
        end if
        x1 = tangle(ica,icb,icg,icd1,allxyz)
        x2 = tangle(ica,icb,icg,icd2,allxyz)
        if (abs(x1) .gt. abs(x2)) then
          if (lchang) then
            write (*,6000) resnam,id
            atmnam (icd1) = ' CD2'
            atmnam (icd2) = ' CD1'
            atmnam (ice1) = ' CE2'
            atmnam (ice2) = ' CE1'
            call prompt (' Swapped CD1/2 and CE1/2')
          else
            call errcon ('Wrong names for CD1/2 and CE1/2')
            call fvalut (' Torsion CD1 :',1,x1)
            call fvalut (' Torsion CD2 :',1,x2)
          end if
          write (*,*)
          return
        end if
        lerror = .false.
        return
c
      else if (resnam .eq. 'ASP') then
        if (min(ica,icb,icg,iod1,iod2) .le. 0) then
          call errcon ('Missing atom(s) !')
          return
        end if
        x1 = tangle(ica,icb,icg,iod1,allxyz)
        x2 = tangle(ica,icb,icg,iod2,allxyz)
        if (abs(x1) .gt. abs(x2)) then
          if (lchang) then
            write (*,6000) resnam,id
            atmnam (iod1) = ' OD2'
            atmnam (iod2) = ' OD1'
            call prompt (' Swapped OD1/2')
          else
            call errcon ('Wrong names for OD1/2')
            call fvalut (' Torsion OD1 :',1,x1)
            call fvalut (' Torsion OD2 :',1,x2)
          end if
          write (*,*)
          return
        end if
        lerror = .false.
        return
c
      else if (resnam .eq. 'GLU') then
        if (min(icb,icg,icd,ioe1,ioe2) .le. 0) then
          call errcon ('Missing atom(s) !')
          return
        end if
        x1 = tangle(icb,icg,icd,ioe1,allxyz)
        x2 = tangle(icb,icg,icd,ioe2,allxyz)
        if (abs(x1) .gt. abs(x2)) then
          if (lchang) then
            write (*,6000) resnam,id
            atmnam (ioe1) = ' OE2'
            atmnam (ioe2) = ' OE1'
            call prompt (' Swapped OE1/2')
          else
            call errcon ('Wrong names for OE1/2')
            call fvalut (' Torsion OE1 :',1,x1)
            call fvalut (' Torsion OE2 :',1,x2)
          end if
          write (*,*)
          return
        end if
        lerror = .false.
        return
c
      else if (resnam .eq. 'ARG') then
        if (min(icd,ine,icz,inh1,inh2) .le. 0) then
          call errcon ('Missing atom(s) !')
          return
        end if
        x1 = tangle(icd,ine,icz,inh1,allxyz)
        x2 = tangle(icd,ine,icz,inh2,allxyz)
        if (abs(x1) .gt. abs(x2)) then
          if (lchang) then
            write (*,6000) resnam,id
            atmnam (inh1) = ' NH2'
            atmnam (inh2) = ' NH1'
            call prompt (' Swapped NH1/2')
          else
            call errcon ('Wrong names for NH1/2')
            call fvalut (' Torsion NH1 :',1,x1)
            call fvalut (' Torsion NH2 :',1,x2)
          end if
          write (*,*)
          return
        end if
        lerror = .false.
        return
      end if
c
      return
      end
c
c
c
      subroutine badman (
     +  option,f7,file7,natoms,maxapr,maxbuf,iresid,atmnam,
     +  iot,names,sorted,allxyz,xf,yf,zf,altloc,resnam,achain,
     +  insert,qatom,batom,ibuff)
c
c ... Bonds/Angles/Dihedrals files
c
      implicit none
c
      integer maxapr,natoms,maxbuf
c
      real allxyz(3,natoms),xf(natoms),yf(natoms),zf(natoms)
      real qatom(natoms),batom(natoms),dmin,dummy
      real dist,angle,tangle
c
      integer f7,ierr,nwrit,ires,napr,i,j,kk,l,k,ii,jj
      integer iot(maxapr),sorted(natoms),insert(natoms),iref
      integer ibuff(maxbuf),iresid(natoms),jref,ll,leng1
c
      logical xinter,lhydro
c
      character option*(*),file7*(*),atmnam(natoms)*4
      character names(0:maxapr)*4,altloc(natoms)*1,line*256
      character resnam(natoms)*3,achain(natoms)*2,xnote*40
c
code ...
c
 4711 format (a6,a)
 4713 format (i5,1x,a4,a1,a3,a2,i4,a1,3x,3f8.3,2f6.2,a40)
c
      if (option(1:4) .eq. 'EXPO') then
c
        write (*,*)
        call textin (
     +    ' Output Bonds/Angles/Dihedrals file ?',file7)
        close (f7)
        call xopxna (f7,file7,xinter(),ierr)
        if (ierr .ne. 0) goto 898
c
        nwrit =0
        ires = 0
        napr = 0
c
        i = 0
 3020   continue
        i = i + 1
        if (i .gt. natoms) goto 3021
c
        if (iresid(i) .eq. ires) then
          napr = napr + 1
          iot (napr) = i
          names (napr) = atmnam (i)
          if (i .eq. natoms) goto 3021
          goto 3020
        end if
c
        if (napr .gt. 0) goto 3021
c
 3022   continue
        ires = iresid(i)
        napr = 1
        iot (1) = i
        names (1) = atmnam(i)
c
        goto 3020
c
 3021   continue
        if (napr .le. 0) goto 3029
c
cc        print *,' NAPR = ',napr
c
c ... sort if more than 2 atoms (put N, CA, C, O, OT1, OT2, CB at the top)
c
        if (napr .eq. 2) then
c
          if (names(2) .lt. names(1)) then
            names(0) = names(2)
            names(2) = names(1)
            names(1) = names(0)
            l = iot(2)
            iot(2) = iot(1)
            iot(1) = l
          end if
c
        else if (napr .gt. 2) then
c
          j = 1
c
          do k=1,napr
            if (names(k) .eq. ' N  ') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8300
            end if
          end do
 8300     continue
c
          do k=1,napr
            if (names(k) .eq. ' CA ') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8301
            end if
          end do
 8301     continue
c
          do k=1,napr
            if (names(k) .eq. ' C  ') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8302
            end if
          end do
 8302     continue
c
          do k=1,napr
            if (names(k) .eq. ' O  ') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8303
            end if
          end do
 8303     continue
c
          do k=1,napr
            if (names(k) .eq. ' OT1') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8304
            end if
          end do
 8304     continue
c
          do k=1,napr
            if (names(k) .eq. ' OT2') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8305
            end if
          end do
 8305     continue
c
          do k=1,napr
            if (names(k) .eq. ' OTX') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8306
            end if
          end do
 8306     continue
c
          do k=1,napr
            if (names(k) .eq. ' OXT') then
              names(0) = names(k)
              names(k) = names(j)
              names(j) = names(0)
              l = iot(k)
              iot(k) = iot(j)
              iot(j) = l
              j = j + 1
              goto 8307
            end if
          end do
 8307     continue
c
          if (j .lt. napr) call nasort((napr-j+1),names(j),iot(j))
c
        end if
c
        do k=1,napr
          nwrit = nwrit + 1
          sorted (nwrit) = iot(k)
          allxyz (1,nwrit) = xf(iot(k))
          allxyz (2,nwrit) = yf(iot(k))
          allxyz (3,nwrit) = zf(iot(k))
        end do
c
        napr = 0
c
        if (nwrit .eq. natoms) goto 3029
        if (i .le. natoms) goto 3022
c
 3029   continue
c
        call jvalut (' Nr of atoms sorted :',1,nwrit)
        write (*,*) '... I am toiling ...'
c
        kk = sorted (1)
        write (xnote,'(3i6)') 0,0,0
        ibuff (1) = 0
        ibuff (2) = 0
        ibuff (3) = 0
        write (line,4713) 1,atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    0.0,0.0,0.0,
     +    qatom(kk),batom(kk),xnote
        write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
c
        kk = sorted(2)
        ibuff (4) = 1
        ibuff (5) = 0
        ibuff (6) = 0
        write (xnote,'(3i6)') 1,0,0
        write (line,4713) 2,atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    dist(2,1,allxyz),0.0,0.0,
     +    qatom(kk),batom(kk),xnote
        write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
c
        kk = sorted(3)
        ibuff (7) = 2
        ibuff (8) = 1
        ibuff (9) = 0
        write (xnote,'(3i6)') 2,1,0
        write (line,4713) 3,atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    dist(3,2,allxyz),
     +    angle(3,2,1,allxyz),0.0,
     +    qatom(kk),batom(kk),xnote
        write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
c
        do i=4,natoms
          ii = sorted(i)
          dmin = 9999.99
          iref = i-1
          if (dist(i,i-1,allxyz) .lt. 1.8 .and.
     +        .not. lhydro(atmnam(sorted(i-1)))) goto 3026
c
          jref = max(3,i-10)
          l = 1
 3028     continue
          do j=i-1,jref,-1
            if (.not. lhydro(atmnam(sorted(j)))) then
              dummy = dist(i,j,allxyz)
              if (dummy .lt. dmin) then
                dmin = dummy
                iref = j
              end if
            end if
          end do
          if (dmin .gt. 1.8 .and. jref .gt. 3) then
            jref = max (3, jref - l*10)
            l = l * 5
            goto 3028
          end if
c
 3026     continue
          jref = 1 + 3*(iref-1)
          jj = iref
          kk = ibuff (jref)
          ll = ibuff (jref+1)
          jref = (i-1)*3 + 1
          ibuff (jref) = jj
          ibuff (jref+1) = kk
          ibuff (jref+2) = ll
          write (xnote,'(3i6)') jj,kk,ll
          write (line,4713) i,atmnam(ii),altloc(ii),
     +      resnam(ii),achain(ii),iresid(ii),insert(ii),
     +      dist   (i,jj,      allxyz),
     +      angle  (i,jj,kk,   allxyz),
     +      tangle (i,jj,kk,ll,allxyz),
     +      qatom(ii),batom(ii),xnote
          write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
        end do
c
        close (f7)
        write (*,*)
c
      else if (option(1:4) .eq. 'SAME') then
c
        write (*,*)
        call textin (
     +    ' Output Bonds/Angles/Dihedrals file ?',file7)
        close (f7)
        call xopxna (f7,file7,xinter(),ierr)
        if (ierr .ne. 0) goto 898
c
        do i=1,natoms
          k = sorted (i)
          allxyz (1,i) = xf(k)
          allxyz (2,i) = yf(k)
          allxyz (3,i) = zf(k)
        end do
c
        kk = sorted (1)
        write (xnote,'(3i6)') 0,0,0
        write (line,4713) 1,atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    0.0,0.0,0.0,
     +    qatom(kk),batom(kk),xnote
        write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
c
        kk = sorted(2)
        write (xnote,'(3i6)') 1,0,0
        write (line,4713) 2,atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    dist(2,1,allxyz),0.0,0.0,
     +    qatom(kk),batom(kk),xnote
        write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
c
        kk = sorted(3)
        write (xnote,'(3i6)') 2,1,0
        write (line,4713) 3,atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    dist(3,2,allxyz),
     +    angle(3,2,1,allxyz),0.0,
     +    qatom(kk),batom(kk),xnote
        write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
c
        do i=4,natoms
          ii = sorted(i)
          jref = (i-1)*3 + 1
          jj = ibuff (jref)
          kk = ibuff (jref+1)
          ll = ibuff (jref+2)
          write (xnote,'(3i6)') jj,kk,ll
          write (line,4713) i,atmnam(ii),altloc(ii),
     +      resnam(ii),achain(ii),iresid(ii),insert(ii),
     +      dist   (i,jj,      allxyz),
     +      angle  (i,jj,kk,   allxyz),
     +      tangle (i,jj,kk,ll,allxyz),
     +      qatom(ii),batom(ii),xnote
          write (f7,4711,err=996) 'BAD   ',(line(1:leng1(line)))
        end do
c
        close (f7)
        write (*,*)
c
      end if
c
      return
c
  898 write (*,*) ' *** ERROR - while opening output BAD file'
      return
c
  996 write (*,*) ' *** ERROR - while writing output PDB file'
      return
c
      end
c
c
c
      subroutine odicts (
     +  option,natoms,ires,bndcut,iot,iresid,atmnam,xf,yf,zf,
     +  altloc,resnam,achain,insert,qatom,batom,inote,otxyz,
     +  f1,isbond,okbond,dismat,deftor,valtor,afftor,tornam)
c
      implicit none
c
      include 'mole_dim.incl'
c
      real xf(maxatm),yf(maxatm),zf(maxatm),qatom(maxatm),bndcut
      real batom(maxatm),otxyz(3,maxapr),dismat(maxapr,maxapr)
      real dist,tangle,dummy,xx,valtor(maxtor),cuthea,mass,radius
c
      integer natoms,iopt,ires,napr,iot(maxapr),i,iresid(maxatm)
      integer insert(maxatm),f1,i1,nl,length,j,j2,j1,k,ncnt,l
      integer deftor(4,maxtor),ierr,leng1
      integer chemic(maxapr)
c
      logical xinter,lhydro,isbond(maxapr,maxapr)
      logical okbond(maxapr,maxapr),afftor(maxapr,maxtor)
c
      character line*(6*maxapr),fulnam*80
      character option*(*),line2*256,junky*256,myres*3
      character atmnam(maxatm)*4,altloc(maxatm)*1,achain(maxatm)*1
      character inote(maxatm)*40,resnam(maxatm)*3,tornam(maxtor)*6
c
      data cuthea /2.3/
c
      save cuthea
c
code ...
c
 4713 format (i5,1x,a4,a1,a3,a2,i4,a1,3x,3f8.3,2f6.2,a40)
c
        iopt = 1
        if (option(1:4) .eq. 'RSR_') iopt = 2
        if (option(1:4) .eq. 'CONN') iopt = 3
        if (option(1:4) .eq. 'TORS') iopt = 4
c
        call ivalin (' Which residue number ?',1,ires)
        if (ires .le. 0) then
          call errcon ('Residue number must be positive !')
          return
        end if
c
        if (iopt .eq. 3 .or. iopt .eq. 4) then
          call prompt (
     +      ' Two cut-offs determine if two atoms are bonded:')
          call prompt (
     +      ' - "Light", if both in {H,He,Li,Be,B,C,N,O,F,Ne}')
          call prompt (
     +      ' - "Heavy", if either in {Na,Mg,Al,Si,P,S,Cl,...}')
          call fvalin (' "Light" cut-off ?',1,bndcut)
          bndcut = max (0.1, bndcut)
          cuthea = max (cuthea,bndcut+0.1)
          call fvalin (' "Heavy" cut-off ?',1,cuthea)
        end if
c
        write (*,*)
        napr = 0
        line = ' '
        do i=1,natoms
          if (iresid(i) .eq. ires) then
c
            if (napr .ge. maxapr) then
              call errcon ('Too many atoms')
              call ivalut (' Maximum :',1,maxapr)
              goto 2206
            end if
c
            napr = napr + 1
            iot (napr) = i
            if (napr .eq. 1) then
              myres = resnam(i)
            end if
            line (1+5*(napr-1):4+5*(napr-1)) = atmnam(i)
            otxyz (1,napr) = xf (i)
            otxyz (2,napr) = yf (i)
            otxyz (3,napr) = zf (i)
            write (junky,4713) i,atmnam(i),altloc(i),
     +        resnam(i),achain(i),iresid(i),insert(i),xf(i),
     +        yf(i),zf(i),qatom(i),batom(i),inote(i)
            write (*,'(a)') junky(1:leng1(junky))
c
c ... get chemical element nr
c
            call elinfo (atmnam(i)(1:2),fulnam,chemic(napr),
     +        mass,radius,.false.)
c
          else
            if (napr .gt. 0) goto 2206
          end if
        end do              
c
 2206   continue
        write (*,*)
        call ivalut (' Nr of atoms found :',1,napr)
        if (napr .le. 0) then
          call errcon ('No atoms means nothing to do')
          return
        end if
c
        if (iopt .eq. 4 .and. natoms .lt. 4) then
          call errcon (' Fewer than 4 atoms means no torsions')
          return
        end if
c
        call textut (' Residue type :',myres)
        call textut (' Atom types   :',line)
c
        if (iopt .eq.1) then
          junky = 'rsfit_'//myres//'.odb'
        else if (iopt .eq. 2) then
          junky = 'rsr_dict_'//myres//'.odb'
        else if (iopt .eq. 3) then
          junky = 'connect_'//myres//'.dat'
        else if (iopt .eq. 4) then
          junky = 'torsion_'//myres//'.dat'
        end if
        call locase (junky)
        call remspa (junky)
        if (iopt .eq. 1 .or. iopt .eq. 2 .or. iopt .eq. 4) then
          call textin (' Datablock file ?',junky)
        else if (iopt .eq. 3) then
          call textin (' Connect file ?',junky)
        end if
        close (f1)
        call xopxua (f1,junky,xinter(),ierr)
        if (ierr .ne. 0) return
c
        if (iopt .eq. 3 .or. iopt .eq. 4) goto 2313
c
c ... rs-fit or rsr datablock file
c
        call stamp (line2)
        write (f1,'(a1,1x,a)') '!',line2(1:leng1(line2))
        write (f1,'(a1)') '!'
        if (iopt .eq.1) then        
          write (f1,'(a,1x,a)') '! RS-fit datablock for',myres
        else if (iopt .eq. 2) then
          write (f1,'(a,1x,a)') '! RSR datablock for',myres
        end if
        write (f1,'(a1)') '!'
c
        if (iopt .eq.1) then        
          junky = 'rsfit_'//myres
        else if (iopt .eq. 2) then
          junky = 'RSR_dict_'//myres
        end if
        call remspa (junky)
        call textut (' Datablock name :',junky)
        nl = 1 + ((napr*5)/71)
        write (line2,'(a,1x,a,1x,i3,1x,i3)')
     +    junky(1:leng1(junky)),'T',nl,70
        call pretty (line2)
        write (f1,'(a)') line2(1:leng1(line2))
c
        i1 = 0
        do i=1,napr,14
          j = min (napr,i+13)
          write (f1,'(a)') line(1+5*(i-1):4+5*(j-1))
          i1 = i1 + 1
        end do
c
        if (i1 .lt. nl) then
          do i=i1+1,nl
            write (f1,*)
          end do
        end if
c
        close (f1)
        call prompt (' Datablock written')
        return
c
c ... connect file or torsion datablock
c
 2313   continue
c
        nl = 0
        do i=1,napr
          do j=i+1,napr
            isbond (i,j) = .false.
            isbond (j,i) = .false.
            okbond (i,j) = .false.
            okbond (j,i) = .false.
c
            dismat (i,j) = dist (i,j,otxyz)
            dismat (j,i) = dismat (i,j)
c
            if (chemic(i).le.10 .and. chemic(j).le.10) then
              if (dismat(i,j) .le. bndcut) then
                if (.not. (lhydro(atmnam(iot(i))) .and.
     +              lhydro(atmnam(iot(j))) )) then
                  isbond (i,j) = .true.
                  isbond (j,i) = .true.
                  nl = nl + 1
                end if
              end if
            else
              if (dismat(i,j) .le. cuthea) then
                if (.not. (lhydro(atmnam(iot(i))) .and.
     +              lhydro(atmnam(iot(j))) )) then
                  isbond (i,j) = .true.
                  isbond (j,i) = .true.
                  nl = nl + 1
                end if
              end if
            end if
          end do
          isbond (i,i) = .false.
          okbond (i,i) = .false.
          dismat (i,i) = 0.0
        end do
        call jvalut (' Nr of bonds :',1,nl)
c
        if (iopt .eq. 4) then
          if (nl .le. 0) then
            call errcon (' No bonds means no torsions')
            return
          end if
          goto 2452
        end if
c
c ... connect file
c
        j2 = 0
        write (f1,'(a3)') myres
        do i=1,napr,10
          j = min (napr,i+9)
          write (f1,'(a4,1x,a)') 'ATOM',line(1+5*(i-1):4+5*(j-1))
          j2 = j2 + 1
        end do
c
        if (nl .gt. 0) then
c
          j1 = 0
c
 2336     continue
c
          i1 = 0
c
          do i=1,napr
            do j=1,napr
              if ( isbond(i,j) .and. (.not. okbond(i,j)) ) then
                j1 = j1 + 1
                if (j1 .eq. 1) then
                  line2 = 'CONNECT - '//line(1+5*(i-1):4+5*(i-1))//
     +                    ' '//line(1+5*(j-1):4+5*(j-1))
                else
                  line2 = 'CONNECT '//line(1+5*(i-1):4+5*(i-1))//
     +                    ' '//line(1+5*(j-1):4+5*(j-1))
                end if
                i1 = 2
                k = j
                okbond (i,j) = .true.
                okbond (j,i) = .true.
                goto 2339
              end if
            end do
          end do
          goto 2341
c
 2339     continue
          if (i1 .ge. 12) goto 2343
          do i=1,napr
            if (isbond(k,i) .and. (.not. okbond(k,i)) ) then
              line2 = line2(1:leng1(line2))//' '//
     +                line(1+5*(i-1):4+5*(i-1))
              i1 = i1 + 1
              okbond (i,k) = .true.
              okbond (k,i) = .true.
              k = i
              goto 2339
            end if
          end do
c
 2343     continue
          if (j1 .eq. 1) then
            line2 = line2(1:leng1(line2))//' +'
          end if
          write (f1,'(a)') line2(1:leng1(line2))
          goto 2336
c
 2341     continue
c
c ... check orphan atoms without any bonds
c
          do i=1,napr
            do j=1,napr
              if (isbond(i,j)) goto 2424
            end do
            write (f1,'(a7,1x,a,1x,a)') 'CONNECT',
     +        line(1+5*(i-1):4+5*(i-1)),line(1+5*(i-1):4+5*(i-1))
            j1 = j1 + 1
 2424       continue
          end do
c
        else
c
c ... no bonds at all; just "connect" individual atoms with themselves
c
          do i=1,napr
            write (f1,'(a7,1x,a,1x,a)') 'CONNECT',
     +        line(1+5*(i-1):4+5*(i-1)),line(1+5*(i-1):4+5*(i-1))
          end do
          j1 = napr
        end if
c
        close (f1)
        call jvalut (' Nr of ATOM records :',1,j2)
        call jvalut (' Nr of CONNECT ,,   :',1,j1)
        call prompt (' Connect file written (append to all.dat)')
        return
c
c ... torsion datablock
c
 2452   continue
c
        dummy = 5.0
c
        ncnt = 0
        do i=1,napr
          do j=1,napr
            if (i.eq.j) goto 2465
            if (chemic(i).le.10.and.chemic(j).le.10) then
              if (dismat(i,j) .gt. bndcut) goto 2465
            else
              if (dismat(i,j) .gt. cuthea) goto 2465
            end if
            do k=1,napr
              if (k.eq.j .or. k.eq.i) goto 2466
              if (chemic(i).le.10.and.chemic(k).le.10) then
                if (dismat(i,k) .gt. bndcut) goto 2466
              else
                if (dismat(i,k) .gt. cuthea) goto 2466
              end if
              do l=1,napr
                if (l.eq.j .or. l.eq.i .or. l.eq.k) goto 2467
                if (chemic(l).le.10.and.chemic(j).le.10) then
                  if (dismat(l,j) .gt. bndcut) goto 2467
                else
                  if (dismat(l,j) .gt. cuthea) goto 2467
                end if
c
c ... we have found a dihedral
c
                xx = tangle(k,i,j,l,otxyz)
                call fixang (xx)
                write (*,'(/1x,a8,4(1x,a4,i4),1x,f8.2)')
     +            'DIHEDRAL',atmnam(iot(k)),iresid(iot(k)),
     +            atmnam(iot(i)),iresid(iot(i)),
     +            atmnam(iot(j)),iresid(iot(j)),
     +            atmnam(iot(l)),iresid(iot(l)),xx
c
c ... skip if torsion near 0 or +/- 180 degrees
c
                if (abs (xx) .le. dummy) then
                  call prompt (' Skip -> torsion ~ 0')
                  goto 2467
                else if (abs(xx-180) .le. dummy .or.
     +                   abs(xx+180) .le. dummy) then
                  call prompt (' Skip -> torsion ~ 180')
                  goto 2467
                end if
c
c ... get all atoms affected by the torsion
c
                do i1=1,napr
                  do j1=1,napr
                    okbond (i1,j1) = .false.
                  end do
                end do
c
c ... flag all atoms connected to atom "J", except "I", as affected
c
                do i1=1,napr
                  if (i1.ne.i) then
                    if (isbond(j,i1)) okbond (i1,i1) = .true.
                  end if
                end do
c
 2514           continue
                j2 = 0
                do i1=1,napr
                  if (okbond(i1,i1)) then
                    do j1=1,napr
                      if (j1.ne.j) then
                        if (isbond(i1,j1) .and.
     +                      (.not. okbond(j1,j1))) then
                           okbond (j1,j1) = .true.
                           j2 = j2 + 1
                        end if
                      end if
                    end do
                  end if
                end do
c
c ... if I or K are flagged, we probably have a ring torsion
c
                if (okbond(i,i) .or. okbond(k,k)) then

                  call prompt (' Skip -> ring torsion')
                  goto 2467
                end if
c
                if (j2 .gt. 0) goto 2514
c
                j2 = 1
                line2 = atmnam (iot(l))
                do i1=1,napr
                  if (i1 .ne. l .and. okbond(i1,i1)) then
                    line2 = line2(1:leng1(line2))//' '//
     +                      atmnam(iot(i1))
                    j2 = j2 + 1
                  end if
                end do
c
c ... if almost all atoms affected, better to use it the other
c     way around
c
                call textut (' Affected atoms :',line2)
                if (j2 .ge. (napr-4) ) then
                  call prompt (' Skip -> too many affected atoms')
                  goto 2467
                end if
c
c ... check if it is a permutation of a previous torsion
c
                if (ncnt .ge. 1) then
                  do i1=1,ncnt
                    if (
     +           (i.eq.deftor(2,i1) .and. j.eq.deftor(3,i1))) then
ccccc     +           (j.eq.deftor(2,i1) .and. i.eq.deftor(3,i1))) then
                      do j1=1,napr
                        if (okbond(j1,j1) .neqv. afftor(j1,i1))
     +                    goto 2532
                      end do
                      call textut (' Permutation of :',tornam(i1))
                      call prompt (' Skip -> permutation')
                      goto 2467
 2532                 continue
                    end if
                  end do
                end if
c
                if (ncnt .ge. maxtor) then
                  call jvalut (' Max nr of torsions :',1,maxtor)
                  call errcon (' Too many torsion angles !!!')
                  goto 2467
                end if
c
c ... we have a new torsion !!!
c
                ncnt = ncnt + 1
                write (tornam(ncnt),'(a3,i3)') 'TOR',ncnt
                call remspa (tornam(ncnt))
                call textut (' OKAY ==> new torsion :',tornam(ncnt))
c
                deftor (1,ncnt) = k
                deftor (2,ncnt) = i
                deftor (3,ncnt) = j
                deftor (4,ncnt) = l
c
                valtor (ncnt) = float(nint(xx))
c
                do i1=1,napr
                  afftor(i1,ncnt) = okbond(i1,i1)
                end do
c
 2467           continue
              end do
 2466         continue
            end do
 2465       continue
          end do
        end do
c
        write (*,*)
        call ivalut (' Nr of unique rotatable torsions :',1,ncnt)
c
c ... now write the file
c
        write (f1,'(a,1x,a)') 'RESIDUE',myres
        nl = 1
c
        do i=1,ncnt
          write (line2,'(a,1x,a,1x,f8.0,4(1x,a4))') 'TORSION',
     +      tornam(i),valtor(i),atmnam(iot(deftor(1,i))),
     +      atmnam(iot(deftor(2,i))),atmnam(iot(deftor(3,i))),
     +      atmnam(iot(deftor(4,i)))
          call pretty (line2)
c
          do j=1,napr
            if (afftor(j,i)) then
              if (length(line2) .ge. 65) then
                write (f1,'(a,1x,a1)') line2(1:leng1(line2)),'\\'
                nl = nl + 1
                line2 = '   '//atmnam(iot(j))
              else
                line2 = line2(1:leng1(line2))//' '//atmnam(iot(j))
              end if
              call pretty (line2)
            end if
          end do
          write (f1,'(a)') line2(1:leng1(line2))
          nl = nl + 1
        end do
c
        close (f1)
        call jvalut (' Nr of lines written :',1,nl)
        call prompt (' Torsion file written (append to torsion.o)')
        return
c
      end
c
c
c
      subroutine shorty (natoms,badcut,xyz,
     +    atmnam,resnam,iresid,achain)
c
      implicit none
c
      include 'mole_dim.incl'
c
      real xyz(3,maxatm),badcut,cog(3,maxres),rad(maxres)
      real dist,dummy,dumdum
c
      integer natoms,i,j,iresid(maxatm),resptr(maxres)
      integer k,l,nres,nowres
c
      logical lhydro
c
      character atmnam(maxatm)*4,achain(maxatm)*2
      character resnam(maxatm)*3
c
code ...
c
c
      nres = 0
      nowres = -1
c
      do i=1,natoms
        if (iresid(i) .ne. nowres) then
          nres = nres + 1
          resptr (nres) = i
          nowres = iresid (i)
          do j=1,3
            cog(j,nres) = xyz(j,i)
          end do
        else
          do j=1,3
            cog(j,nres) = cog(j,nres) + xyz(j,i)
          end do
        end if
      end do
      resptr (nres+1) = natoms + 1
c
      do i=1,nres
        nowres = resptr(i+1) - resptr(i)
        do j=1,3
          cog(j,i) = cog(j,i)/float(nowres)
        end do
        rad (i) = 0.0
        do j=resptr(i),resptr(i+1)-1
          dummy =  (xyz(1,j)-cog(1,i))**2 + (xyz(2,j)-cog(2,i))**2 +
     +             (xyz(3,j)-cog(3,i))**2
          rad (i) = max (rad(i), sqrt(dummy))
        end do
c
cc      write (*,'(1x,a3,1x,i5,1x,3f8.3,1x,f8.2,1x,i5)')
cc     +  resnam(resptr(i)),iresid(resptr(i)),
cc     +  (cog(j,i),j=1,3),rad(i),nowres
c
      end do
c
      nowres = 0
      do i=1,nres
        do j=i,nres
          dumdum = dist (i,j,cog)
          if (dumdum .le. (rad(i)+rad(j)+badcut+1.0) ) then
            do k=resptr(i),resptr(i+1)-1
c
              if (lhydro(atmnam(k))) goto 2354
              do l=resptr(j),resptr(j+1)-1
c
                if (lhydro(atmnam(l))) goto 2356
                if (k .eq. l) goto 2354
                dummy = dist(k,l,xyz)
c
c ... skip disulfide links
c
                if (atmnam(k).eq.' SG ' .and.
     +              atmnam(l).eq.' SG ') goto 2356
c
                if (abs(iresid(k)-iresid(l)) .le. 1) then
                  if (dummy .le. 1.1) then
                    write (*,6000) resnam(k),achain(k),iresid(k),
     +                atmnam(k),resnam(l),achain(l),iresid(l),
     +                atmnam(l),dummy
                    nowres = nowres + 1
                  end if
                else
                  if (dummy .le. badcut) then
                    write (*,6000) resnam(k),achain(k),iresid(k),
     +                atmnam(k),resnam(l),achain(l),iresid(l),
     +                atmnam(l),dummy
                    nowres = nowres + 1
                  end if
                end if
c
 2356           continue
              end do
 2354         continue
            end do
          end if
        end do
      end do
c
 6000 format (' Dist ',a3,'-',a2,i5,'-',a4,' TO ',a3,'-',a2,i5,
     +  '-',a4,' = ',f8.2,' A')
c
      call jvalut (' Short contacts :',1,nowres)
      call prompt (' Any hydrogen atoms have been ignored')
      call prompt (
     +  ' Short contacts between I-1,I,I+1 residues ignored')
      call prompt (
     +  ' (unless they are shorter than 1.1 A)')
c
      call r5shor (nowres,badcut)
c
      return
      end
c
c
c
      subroutine extinc (nat,atmnam,resnam)
c
c ... estimate extinction coefficient at 280 nm
c
      implicit none
c
      real e
c
      integer nat,i,nr,ny,nw,nc
c
      character atmnam(nat)*4,resnam(nat)*3
c
code ...
c
      call prompt (' Reference: Gill & Von Hippel,')
      call prompt (' Anal. Biochem. 182: 319-326 (1989).')
      call prompt (' e280 = (1280*TYR + 5690*TRP + 120*CYS)')
      call prompt (' accuracy roughly 5%; units M-1 cm-1')
c
      nr = 0
      ny = 0
      nw = 0
      nc = 0
c
      do i=1,nat
        if (atmnam(i) .eq. ' CA ') then
          nr = nr + 1
          if (resnam(i) .eq. 'TYR') then
            ny = ny + 1
          else if (resnam(i) .eq. 'TRP') then
            nw = nw + 1
          else if (resnam(i) .eq. 'CYS') then
            nc = nc + 1
          end if
        end if
      end do
c
      call ivalut (' Nr of residues :',1,nr)
      call ivalut (' Nr of TYR      :',1,ny)
      call ivalut (' Nr of TRP      :',1,nw)
      call ivalut (' Nr of CYS      :',1,nc)
c
      e = 1280.0*ny + 5690.0*nw + 120.0*nc
      call rvalut (' e280 (M-1 cm-1) ~',1,e)
c
      return
      end
c
c
c
      subroutine minclu (natoms,xf,yf,zf,buffer)
c
c ... MINCLU - subroutine to get least-cluttered view of
c              a set of NATOMS atoms with coordinates
c              (xf,yf,zf)
c
c ... GJ Kleywegt @ 950705
c
c assumption: we look along the Z-axis; no perspective
c
c least-cluttered view means: the sum of the Z-components
c of all inter-atomic difference vectors I,J (J>I,I=1,NATOMS-1)
c is at a minimum (e.g.: three atoms, least cluttered if they
c are all in the XY plane)
c
c i.e.: SUM (d.z) = minimal, where Z = unit vector (0 0 1)
c and the sum extends over all difference vectors
c d = Rx, where R is the unknown rotation matrix and x are
c the difference vectors before the rotation
c
c use Euler angles, alpha, beta, gamma
c since we look along Z, the alpha rotation does not matter
c -> choose alpha = 0.0
c
c now the rotation matrix elements are (cb=cos(beta) etc):
c
c     ( cb cg		-cb sg		sb  )
c R = ( sg		cg		0   )
c     ( -sb cg		sb sg		cb  )
c
c if x = diff. vector (dx dy dz) then the sum to be minimised is:
c SUM (d.z) = SUM (Rx).(0 0 1) = ... =
c SUM { -sb cd dx + sb sg dy + cb dz } = minimal
c
c minimum means d()/dbeta = d()/dgamma = 0
c
c -> gamma: sb sg SUM(dx) + sb cg SUM(dy) = 0
c    beta may not be zero/180 (then we would only rotate around z)
c    => sb may not be zero
c    => sg SUM(dx) + cg SUM(dy) = 0 (1)
c
c -> beta : -cb cg SUM(dx) + cb sg SUM(dy) - sb SUM(dz) = 0 (2)
c
c from (1) -> tan(gamma) = -SUM(dy)/SUM(dx)
c this gives gamma; substitute into (2) ->
c tan(beta) = (-cg SUM(dx) + sg SUM(dy)) / SUM(dz)
c this gives beta
c
c the difference vector sums can be evaluated very efficiently, e.g.:
c
c SUM(dx) = SUM     SUM     (Xj-Xi) =
c           i=1,N-1 j=i+1,N
c
c = ... = SUM     { (SUM   Xi ) - (SUM   Xj) - (N-i+2)*Xi }
c         i=1,N-1    i=1,N         j=1,i
c
c and similarly for Y and Z
c
c afterwards: apply the Euler rotations (or the matrix) to the
c original coordinates to get the least-cluttered view along the
c Z axis
c
c ... this subroutine calculates the Euler angles
c
      implicit none
c
      real twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      integer natoms
c
      real xf(natoms),yf(natoms),zf(natoms),buffer(natoms)
      real sum,sumdx,sumdy,sumdz,al,be,ga,dum
c
      integer i
c
code ...
c
c ... get sum of X part of all difference vectors I,J, J>I
c
      sum = 0.0
      do i=1,natoms
        sum = sum + xf(i)
        buffer(i) = sum
      end do
c
c ... now: SUM = SUM (Xi), i=1,N
c          BUFFER (I) = SUM (Xj), j=1,I
c
      sumdx = (natoms-1)*sum
      do i=1,natoms-1
        sumdx = sumdx - buffer(i) - (natoms-i+2)*xf(i)
      end do
c
c ... get sum of Y part of all difference vectors I,J, J>I
c
      sum = 0.0
      do i=1,natoms
        sum = sum + yf(i)
        buffer(i) = sum
      end do
      sumdy = (natoms-1)*sum
      do i=1,natoms-1
        sumdy = sumdy - buffer(i) - (natoms-i+2)*yf(i)
      end do
c
c ... get sum of Z part of all difference vectors I,J, J>I
c
      sum = 0.0
      do i=1,natoms
        sum = sum + zf(i)
        buffer(i) = sum
      end do
      sumdz = (natoms-1)*sum
      do i=1,natoms-1
        sumdz = sumdz - buffer(i) - (natoms-i+2)*zf(i)
      end do
c
ccc      write (*,'(1x,a15,3f15.1)') 'SUMDX/Y/Z',sumdx,sumdy,sumdz
c
c ... alpha = 0.0
c
      al = 0.0
c
c ... tan gamma = -SUMdy / SUMdx
c
      ga = atan2 (-sumdy,sumdx)
c
c ... tan beta = (-cos gamma SUMdx + sin gamma SUMdy) / SUMdz
c
      dum = -cos(ga)*sumdx + sin(ga)*sumdy
      be = atan2 (dum,sumdz)
c
      call fvalut (' Euler ALPHA :',1,rtodeg*al)
      call fvalut (' Euler BETA  :',1,rtodeg*be)
      call fvalut (' Euler GAMMA :',1,rtodeg*ga)
c
c ... check that the derivatives are really zero (within round-off)
c
ccc      print *,' SG*X + CG*Y = 0 ? = ',sin(ga)*sumdx+cos(ga)*sumdy
ccc      print *,' -CB*CG*X + CB*SG*Y - SB*Z = 0 ? = ',
ccc     +  -cos(be)*cos(ga)*sumdx + cos(be)*sin(ga)*sumdy
ccc     +  - sin(be)*sumdz
c
      return
      end
c
c
c
      subroutine bradpl (iunit,filnam,nat,xyz,b)
c
      implicit none
c
      integer maxbin
      parameter (maxbin=100)
c
      integer iunit,nat,i,j,k,kmax,ierr,leng1
c
      real xyz(3,nat),b(nat)
c
      integer nbin(0:maxbin)
c
      real bbin(0:maxbin),cog(3),d,bmax
c
      logical xinter
c
      character filnam*(*),line*128
c
code ...
c
      cog(1)=0.0
      cog(2)=0.0
      cog(3)=0.0
      do i=1,nat
        do j=1,3
          cog(j)=cog(j)+xyz(j,i)
        end do
      end do
      cog(1)=cog(1)/float(nat)
      cog(2)=cog(2)/float(nat)
      cog(3)=cog(3)/float(nat)
c
      do i=0,maxbin
        nbin(i)=0
        bbin(i)=0.0
      end do
c
      kmax = 0
      do i=1,nat
        d = 0.0
        do j=1,3
          d = d + (xyz(j,i)-cog(j))**2
        end do
        d = sqrt(d)
        k = max (0, min (maxbin,int(d/2.0)))
        if (k.gt.kmax) kmax = k
        nbin (k) = nbin(k) + 1
        bbin (k) = bbin(k) + b(i)
      end do
c
 6000 format (' Shell ',f6.1,' - ',f6.1,' A - ',i6,' atoms; <B> =',
     +  f6.2,' A**2')
      bmax = 0.0
      do i=0,maxbin
        if (nbin(i).gt.0) then
          bbin(i)=bbin(i)/float(nbin(i))
          write (*,6000) 2.0*i,2.0*(i+1),nbin(i),bbin(i)
        end if
        if (bbin(i).gt.bmax) bmax = bbin(i)
      end do
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
c
 5000 format (a6,1x,a)
 5010 format (a6,1x,12i6)
 5020 format (a6,1x,6f12.4)
c
      write (iunit,5000) 'REMARK',' Radial temperature factor plot'
      write (iunit,5000) 'REMARK',(' '//line(1:leng1(line)))
      write (iunit,5000) 'REMARK'
      write (iunit,5010) 'NPOINT',kmax+1
      write (iunit,5010) 'COLOUR',4
      write (iunit,5000) 'XLABEL','Distance from Centre-of-Gravity (A)'
      write (iunit,5000) 'YLABEL','Average B in 2 A bins (A**2)'
      write (iunit,5020) 'XYVIEW',0.0,float(3*kmax+1),0.0,bmax+1.0
      write (iunit,5020) 'XLIMIT',0.0,3.0
      write (iunit,5000) 'YVALUE','(8f10.4)'
c
      write (iunit,'(8f10.4)') (bbin(i),i=0,kmax)
c
      write (iunit,5000) 'END   '
      close (iunit)
c
      call prompt (' Plot file written')
c
      return
      end
c
c
c
      subroutine helstr (type,start,end,natoms,
     +  atomnr,atmnam,altloc,resnam,achain,iresid,insert,
     +  xf,yf,zf,allxyz,qatom,batom,inote,xxyyzz,ibuff,maxatm)
c
      implicit none
c
      integer maxatm,natoms
c
      real start(3),end(3),xf(maxatm),yf(maxatm),zf(maxatm)
      real allxyz(3,maxatm),xxyyzz(3,maxatm),qatom(maxatm)
      real batom(maxatm)
      real try(3,3),act(3,3),rt(12),rms,x1(3),x2(3)
      real d,d1,d2,q1,q2,qq
      real dist
c
      integer atomnr(maxatm),iresid(maxatm),insert(maxatm)
      integer ibuff(3,maxatm),i,j,ierr,nres,ires
c
      character type*1,atmnam(maxatm)*4,altloc(maxatm)*1
      character resnam(maxatm)*3,achain(maxatm)*2
      character inote(maxatm)*40,dumstr*1
c
code ...
c
c ... get desired length of SSE
c
      d = sqrt ( (start(1)-end(1))**2 + (start(2)-end(2))**2 + 
     +           (start(3)-end(3))**2 )
      call fvalut (' Length (A) :',1,d)
c
c ... helix rises 1.46 A per residue
c
      if (type .eq. 'A') then
        nres = 1 + int(d/1.46)
      else
c
c ... strand rises 3.32 A per residue
c
        nres = 1 + int(d/3.32)
      end if
c
c ... can't get too many atoms
c
      if (5*nres .gt. maxatm) nres = maxatm/5
      call jvalut (' Nr of residues :',1,nres)
c
c ... poly-ALA
c
      natoms = 5 * nres
c
ccc      call prompt (' Generating internal coordinates')
c
      if (type .eq. 'A') then
c
c ... ideal helix geometry (internal coordinates)
c
cBAD      11  N   ALA     3       1.317 117.032 -50.524  1.00 20.00     8     7     6
cBAD      12  CA  ALA     3       1.483 121.842 177.987  1.00 20.00    11     8     7
cBAD      13  C   ALA     3       1.534 111.564 -53.333  1.00 20.00    12    11     8
cBAD      14  O   ALA     3       1.270 119.208 133.221  1.00 20.00    13    12    11
cBAD      15  CB  ALA     3       1.549 109.843 178.252  1.00 20.00    12    11     8
c
        ires = 0
        do i=1,natoms-1,5
c
          ires = ires + 1
c
c ... N
c
          j = i
          xxyyzz(1,j) = 1.333
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 117.1
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = -50.0
          ibuff (3,j) = max (0,j-5)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' N  '
c
c ... CA
c
          j = i+1
          xxyyzz(1,j) = 1.480
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 121.5
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 180.0
          ibuff (3,j) = max (0,j-5)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' CA '
c
c ... C
c
          j = i+2
          xxyyzz(1,j) = 1.520
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 110.0
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = -60.0
          ibuff (3,j) = max (0,j-5)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' C  '
c
c ... O
c
          j = i+3
          xxyyzz(1,j) = 1.245
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 120.0
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = 133.0
          ibuff (3,j) = max (0,j-3)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' O  '
c
c ... CB
c
          j = i+4
          xxyyzz(1,j) = 1.548
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 110.8
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 180.0
          ibuff (3,j) = max (0,j-7)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' CB '
c
c ... CB first residue different
c
          if (j .eq. 5) then
            xxyyzz(1,j) = 2.568
            ibuff (1,j) = 3
            xxyyzz(2,j) = 32.68
            ibuff (2,j) = 2
            xxyyzz(3,j) = -123.0
            ibuff (3,j) = 1
          end if
c
        end do
c
      else
c
c
c ... ideal strand geometry (internal coordinates)
c
cBAD      11  N   ALA     3       1.305 115.737 147.052  1.00 20.00     8     7     6
cBAD      12  CA  ALA     3       1.466 123.431-179.357  1.00 20.00    11     8     7
cBAD      13  C   ALA     3       1.532 108.628-131.194  1.00 20.00    12    11     8
cBAD      14  O   ALA     3       1.282 115.511 -52.183  1.00 20.00    13    12    11
cBAD      15  CB  ALA     3       1.551 108.713 105.733  1.00 20.00    12    11     8
c
        ires = 0
        do i=1,natoms-1,5
c
          ires = ires + 1
c
c ... N
c
          j = i
          xxyyzz(1,j) = 1.32
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 115.7
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 120.0
          ibuff (3,j) = max (0,j-5)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' N  '
c
c ... CA
c
          j = i+1
          xxyyzz(1,j) = 1.466
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 121.9
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 180.0
          ibuff (3,j) = max (0,j-5)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' CA '
c
c ... C
c
          j = i+2
          xxyyzz(1,j) = 1.530
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 108.6
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = -125.0
          ibuff (3,j) = max (0,j-5)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' C  '
c
c ... O
c
          j = i+3
          xxyyzz(1,j) = 1.255
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 120.0
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = -40.0
          ibuff (3,j) = max (0,j-3)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' O  '
c
c ... CB
c
          j = i+4
          xxyyzz(1,j) = 1.548
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 111.3
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 108.9
          ibuff (3,j) = max (0,j-7)
          atomnr (j) = j
          iresid (j) = ires
          atmnam (j) = ' CB '
c
c ... CB first residue different
c
          if (j .eq. 5) then
            xxyyzz(1,j) = 2.487
            ibuff (1,j) = 3
            xxyyzz(2,j) = 35.99
            ibuff (2,j) = 2
            xxyyzz(3,j) = -121.3
            ibuff (3,j) = 1
          end if
c
        end do
c
      end if
c
ccc 6000 format (i6,3f8.2,3i6)
ccc      do i=1,natoms
ccc        write (*,6000) i,(xxyyzz(j,i),j=1,3),(ibuff(j,i),j=1,3)
ccc      end do
c
c ... convert internal to Cartesian coordinates
c
ccc      call prompt (' Converting to Cartesian coordinates')
c
      call c2cart (natoms,allxyz,xxyyzz,ibuff)
c
c ... now set rest of attributes
c
ccc      call prompt (' Setting other attributes to defaults')
c
      dumstr = ' '
      do i=1,natoms
        qatom (i) = 1.0
        batom (i) = 20.0
        read (dumstr,'(a1)') insert(i)
        resnam(i) = 'ALA'
        achain(i) = ' Z'
        inote (i) = ' '
        altloc(i) = ' '
      end do
c
c ... get actual distance between first and last CA
c
      d = dist (2,(nres-1)*5+2,allxyz)
      call fvalut (' Distance terminal CAs (A) :',1,d)
c
c ... replace "END" by a point on the line START -> END which has
c     distance D to START and is closest to the old END
c
c ... line START -> END has equation:
c
c         ( Xstart )       (Xend - Xstart)
c ... x = ( Ystart ) + L * (Yend - Ystart)
c         ( Zstart )       (Zend - Zstart)
c
c ... there are two points on this line which have distance D to
c     the point N; these occur for values of L = +/- SQRT (P),
c     with P = D**2 / { (Xend - Xstart)**2 + ... }
c
c ... try both roots and keep the one which is closest to the old
c     point END
c
      qq = d*d / ( (start(1)-end(1))**2 + (start(2)-end(2))**2 +
     +             (start(3)-end(3))**2 )
c
c ... + root
c
      q1 = sqrt(qq)
      x1(1) = start(1) + q1*(end(1)-start(1))
      x1(2) = start(2) + q1*(end(2)-start(2))
      x1(3) = start(3) + q1*(end(3)-start(3))
      d1 = sqrt ( (x1(1)-end(1))**2 + (x1(2)-end(2))**2 +
     +            (x1(3)-end(3))**2 )
c
ccc      print *,' + root = ',q1,' dist = ',d1
ccc      call fvalut (' + point :',3,x1)
c
c ... - root
c
      q2 = -q1
      x2(1) = start(1) + q2*(end(1)-start(1))
      x2(2) = start(2) + q2*(end(2)-start(2))
      x2(3) = start(3) + q2*(end(3)-start(3))
      d2 = sqrt ( (x2(1)-end(1))**2 + (x2(2)-end(2))**2 +
     +            (x2(3)-end(3))**2 )
c
ccc      print *,' - root = ',q2,' dist = ',d2
ccc      call fvalut (' - point :',3,x2)
c
      if (d1 .lt. d2) then
        try(1,3) = x1(1)
        try(2,3) = x1(2)
        try(3,3) = x1(3)
      else
        try(1,3) = x2(1)
        try(2,3) = x2(2)
        try(3,3) = x2(3)
      end if
c
ccc      call fvalut (' Actual end point (CA) :',3,try(1,3))
c
      try(1,1) = start(1)
      try(2,1) = start(2)
      try(3,1) = start(3)
c
      try(1,2) = 0.5*(try(1,1)+try(1,3))
      try(2,2) = 0.5*(try(2,1)+try(2,3))
      try(3,2) = 0.5*(try(3,1)+try(3,3))
c
      act(1,1) = allxyz(1,2)
      act(2,1) = allxyz(2,2)
      act(3,1) = allxyz(3,2)
c
      act(1,3) = allxyz(1,(nres-1)*5+2)
      act(2,3) = allxyz(2,(nres-1)*5+2)
      act(3,3) = allxyz(3,(nres-1)*5+2)
c
      act(1,2) = 0.5*(act(1,1)+act(1,3))
      act(2,2) = 0.5*(act(2,1)+act(2,3))
      act(3,2) = 0.5*(act(3,1)+act(3,3))
c
c ... add some "random noise" to the user's desired positions
c     to make sure that the LSQ routine works (otherwise we
c     have three points exactly co-linear -> LSQ fails)
c
      do i=1,3
        do j=1,3
          call gkrand (qq,-0.05,0.05,0)
          try(i,j)=try(i,j)+qq
        end do
      end do
c
c ... get operator to position SSE according to the user's
c     target coordinates
c
      call lsqgjk (try,act,3,rms,rt,ierr)
c
ccc      call rvalut (' Rotation    :',9,rt)
ccc      call fvalut (' Translation :',3,rt(10))
ccc      print *,' rms = ',rms,' error = ',ierr
c
ccc      call anancs (1,rt,.true.,ierr)
c
c ... apply operator to all atoms
c
      do i=1,natoms
        call vecrtv (allxyz(1,i),x1,1,rt(1),rt(10))
        xf(i) = x1(1)
        yf(i) = x1(2)
        zf(i) = x1(3)
      end do
c
      call prompt (' All done ...')
c
      return
      end
c
c
c
      subroutine chkrem (lpres,string)
c
c ... check if a certain REMARK (etc.) card is already present
c
      implicit none
c
      include 'mole_dim.incl'
c
      integer i,length,ll
c
      logical lpres
c
      character string*(*)
c
code ...
c
      lpres = .false.
      if (nrem .le. 0) return
c
      ll = length(string)
c
      do i=1,nrem
        if (index(remark(i),string(1:ll)) .gt. 0) then
          lpres = .true.
          return
        end if
      end do
c
      return
      end
c
c
c
      subroutine remrk5 ()
c
c ... add REMARK 5 cards
c
      implicit none
c
      include 'mole_dim.incl'
c
      integer cnts(*),n,nonh,nbo,m
c
      real bstats(*),bmax(*),bndcut,nobcut,x,cut,rmsd,corr
c
      logical ldorem,lkeeph
c
      character whichn*2
c
code ...
c
 6600 format ('REMARK   5',a,i8)
 6605 format ('REMARK   5',10(1x,a))
 6610 format ('REMARK   5',a,f8.2)
c
      return
c
c ... CA-CA DISTANCE DICTRIBUTION
c
      entry r5caca (n,cnts)
c
      call chkrem (ldorem,' CA-CA DISTANCE DISTRIBUTION')
      if (ldorem) return
      if (nrem.gt.(maxcom-14)) return
      if (n .lt. 2) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' CA-CA DISTANCE DISTRIBUTION'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF DISTANCES      :',n
      write (remark(nrem+5),6600)
     +  '   SHORT (<= 2.8 A)    (NR) :',cnts(1)
      x = 100.0*float(cnts(1))/float(n)
      write (remark(nrem+6),6610)
     +  '   SHORT (<= 2.8 A)     (%) :',x
      write (remark(nrem+7),6600)
     +  '   CIS (2.8 - 3.0 A)   (NR) :',cnts(2)
      x = 100.0*float(cnts(2))/float(n)
      write (remark(nrem+8),6610)
     +  '   CIS (2.8 - 3.0 A)    (%) :',x
      write (remark(nrem+9),6600)
     +  '   POOR (3.0 - 3.7 A)  (NR) :',cnts(3)
      x = 100.0*float(cnts(3))/float(n)
      write (remark(nrem+10),6610)
     +  '   POOR (3.0 - 3.7 A)   (%) :',x
      write (remark(nrem+11),6600)
     +  '   TRANS (3.7 - 3.9 A) (NR) :',cnts(4)
      x = 100.0*float(cnts(4))/float(n)
      write (remark(nrem+12),6610)
     +  '   TRANS (3.7 - 3.9 A)  (%) :',x
      write (remark(nrem+13),6600)
     +  '   LONG (> 3.9 A)      (NR) :',cnts(5)
      x = 100.0*float(cnts(5))/float(n)
      write (remark(nrem+14),6610)
     +  '   LONG (> 3.9 A)       (%) :',x
c
      nrem = nrem + 14
c
      return
c
c ... CA RAMACHANDRAN PLOT
c
      entry r5cara (n,cnts)
c
      call chkrem (ldorem,' CA RAMACHANDRAN PLOT')
      if (ldorem) return
      if (nrem.gt.(maxcom-12)) return
      if (n .lt. 2) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' CA RAMACHANDRAN PLOT'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF RESIDUES       :',n
      write (remark(nrem+5),6600)
     +  '   CORE REGIONS        (NR) :',cnts(4)
      x = 100.0*float(cnts(4))/float(n)
      write (remark(nrem+6),6610)
     +  '   CORE REGIONS         (%) :',x
      write (remark(nrem+7),6600)
     +  '   ADDITIONAL REGIONS  (NR) :',cnts(3)
      x = 100.0*float(cnts(3))/float(n)
      write (remark(nrem+8),6610)
     +  '   ADDITIONAL REGIONS   (%) :',x
      write (remark(nrem+9),6600)
     +  '   GENEROUS REGIONS    (NR) :',cnts(2)
      x = 100.0*float(cnts(2))/float(n)
      write (remark(nrem+10),6610)
     +  '   GENEROUS REGIONS     (%) :',x
      write (remark(nrem+11),6600)
     +  '   DISALLOWED REGIONS  (NR) :',cnts(1)
      x = 100.0*float(cnts(1))/float(n)
      write (remark(nrem+12),6610)
     +  '   DISALLOWED REGIONS   (%) :',x
c
      nrem = nrem + 12
c
      return
c
c ... SHORT CONTACTS
c
      entry r5shor (n,cut)
c
      x = cut
c
      call chkrem (ldorem,' SHORT CONTACTS')
      if (ldorem) return
      if (nrem.gt.(maxcom-5)) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' SHORT CONTACTS (EXCL. HYDROGENS AND SYMMETRY)'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF SHORT CONTACTS :',n
      write (remark(nrem+5),6610)
     +  '   DISTANCE CUTOFF      (A) :',x
c
      nrem = nrem + 5
c
      return
c
c ... BONDED B-FACTORS
c
      entry r5bond (nonh,bndcut,nbo,rmsd,corr)
c
      call chkrem (ldorem,' BONDED B-FACTORS')
      if (ldorem) return
      if (nrem.gt.(maxcom-8)) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' BONDED B-FACTORS (EXCL. HYDROGENS)'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF NON-H ATOMS    :',nonh
      write (remark(nrem+5),6610)
     +  '   BOND DISTANCE CUTOFF (A) :',bndcut
      write (remark(nrem+6),6600)
     +  '   NUMBER OF BONDS          :',nbo
      write (remark(nrem+7),6610)
     +  '   RMSD B-BONDED     (A**2) :',rmsd
      write (remark(nrem+8),6610)
     +  '   CORR. COEFF. B-BONDED    :',corr
c
      nrem = nrem + 8
c
      return
c
c ... NON-BONDED B-FACTORS (EXCL. HYDROGENS)
c
      entry  r5non1 (bndcut,nobcut,nonh,nbo,rmsd,corr)
c
      call chkrem (ldorem,' NON-BONDED B-FACTORS (EXCL. HYDR')
      if (ldorem) return
      if (nrem.gt.(maxcom-9)) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' NON-BONDED B-FACTORS (EXCL. HYDROGENS)'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF NON-H ATOMS    :',nonh
      write (remark(nrem+5),6610)
     +  '   BOND DISTANCE CUTOFF (A) :',bndcut
      write (remark(nrem+6),6610)
     +  '   NON-BONDED CUTOFF    (A) :',nobcut
      write (remark(nrem+7),6600)
     +  '   NUMBER OF INTERACTIONS   :',nbo
      write (remark(nrem+8),6610)
     +  '   RMSD B-NON-BONDED (A**2) :',rmsd
      write (remark(nrem+9),6610)
     +  '   CORR. COEFF. B-NON-BONDED:',corr
c
      nrem = nrem + 9
c
      return
c
c ... NON-BONDED B-FACTORS (ONLY N AND O)
c
      entry  r5non2 (bndcut,nobcut,nonh,nbo,rmsd,corr)
c
      call chkrem (ldorem,' NON-BONDED B-FACTORS (ONLY N')
      if (ldorem) return
      if (nrem.gt.(maxcom-9)) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' NON-BONDED B-FACTORS (ONLY N AND O)'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF N AND O ATOMS  :',nonh
      write (remark(nrem+5),6610)
     +  '   BOND DISTANCE CUTOFF (A) :',bndcut
      write (remark(nrem+6),6610)
     +  '   NON-BONDED CUTOFF    (A) :',nobcut
      write (remark(nrem+7),6600)
     +  '   NUMBER OF INTERACTIONS   :',nbo
      write (remark(nrem+8),6610)
     +  '   RMSD B-NON-BONDED (A**2) :',rmsd
      write (remark(nrem+9),6610)
     +  '   CORR. COEFF. B-NON-BONDED:',corr
c
      nrem = nrem + 9
c
      return
c
c ... NON-BONDED B-FACTORS (EXCL. H, N, O)
c
      entry  r5non3 (bndcut,nobcut,nonh,nbo,rmsd,corr)
c
      call chkrem (ldorem,' NON-BONDED B-FACTORS (EXCL. H,')
      if (ldorem) return
      if (nrem.gt.(maxcom-9)) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      write (remark(nrem+2),6600)
     +  ' NON-BONDED B-FACTORS (EXCL. H, N, O)'
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      write (remark(nrem+4),6600)
     +  '   NUMBER OF NON-H/N/O ATOMS:',nonh
      write (remark(nrem+5),6610)
     +  '   BOND DISTANCE CUTOFF (A) :',bndcut
      write (remark(nrem+6),6610)
     +  '   NON-BONDED CUTOFF    (A) :',nobcut
      write (remark(nrem+7),6600)
     +  '   NUMBER OF INTERACTIONS   :',nbo
      write (remark(nrem+8),6610)
     +  '   RMSD B-NON-BONDED (A**2) :',rmsd
      write (remark(nrem+9),6610)
     +  '   CORR. COEFF. B-NON-BONDED:',corr
c
      nrem = nrem + 9
c
      return
c
c ... B-FACTOR STATISTICS
c
      entry  r5bq (whichn,lkeeph,cnts,bstats,bmax)
c
      if (nrem.gt.(maxcom-31)) return
c
      call prompt (' Generating REMARK records ...')
c
      write (remark(nrem+1),6600)
      if (lkeeph) then
        write (remark(nrem+2),6605)
     +    'B-FACTOR STATISTICS FOR CHAIN',whichn,
     +    '(INCL. HYDROGENS)'
      else
        write (remark(nrem+2),6605)
     +    'B-FACTOR STATISTICS FOR CHAIN',whichn,
     +    '(EXCL. HYDROGENS)'
      end if
      write (remark(nrem+3),6600)
     +  '   PROGRAM USED: MOLEMAN'
      nrem = nrem + 3
c
      m = 0
c
      if (cnts(1) .gt. 0) then
        write (remark(nrem+1),6600)
     +    '   PROTEIN MAIN CHAIN'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(1)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(1)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(1)
        nrem = nrem + 4
      end if
c
      if (cnts(2) .gt. 0) then
        write (remark(nrem+1),6600)
     +    '   PROTEIN SIDE CHAIN'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(2)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(2)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(2)
        nrem = nrem + 4
      end if
c
      if (cnts(6) .gt. 0) then
        write (remark(nrem+1),6600)
     +    '   PROTEIN ALL ATOMS'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(6)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(6)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(6)
        nrem = nrem + 4
        m = m + 1
      end if
c
      if (cnts(7) .gt. 0) then
        write (remark(nrem+1),6600)
     +    '   LIGAND/SUBSTRATE'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(7)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(7)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(7)
        nrem = nrem + 4
        m = m + 1
      end if
c
      if (cnts(4) .gt. 0) then
        write (remark(nrem+1),6600)
     +    '   WATER MOLECULES'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(4)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(4)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(4)
        nrem = nrem + 4
        m = m + 1
      end if
c
      if (cnts(3) .gt. 0) then
        write (remark(nrem+1),6600)
     +    '   OTHER ENTITIES'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(3)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(3)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(3)
        nrem = nrem + 4
        m = m + 1
      end if
c
      if (cnts(5) .gt. 0 .and. m .gt. 1) then
        write (remark(nrem+1),6600)
     +    '   ALL ATOMS'
        write (remark(nrem+2),6600)
     +    '     NUMBER OF ATOMS        :',cnts(5)
        write (remark(nrem+3),6610)
     +    '     AVERAGE B       (A**2) :',bstats(5)
        write (remark(nrem+4),6610)
     +    '     MAXIMUM B       (A**2) :',bmax(5)
        nrem = nrem + 4
      end if
c
      return
c
c ... end of subroutine
c
      end
