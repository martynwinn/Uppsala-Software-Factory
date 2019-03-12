      program xpand
c
c ... program XPAND
c
c ... NCS and SGS expansion of molecules
c
c ... Gerard Kleywegt @ 940512
c
c ... to do:
c     - S = expand into shell (use mask)
c
      include 'xpand.incl'
c
      real rtinv(12,maxncs)
c
      integer length,ierr,i,j,lunit,munit
c
      logical linter,xinter
c
      character*80 pdbin,pdbut,ncsfil,line,valid,par,fmt,symfil,ncsut
      character*80 rtout,cnsut,spgrp
      character*1  task,partyp,askme
c
      data cpoint /0.0,0.0,0.0/
c
code ...
c
      call gkinit (prognm,vers)
c
      call jvalut (' Max nr of input atoms   = ',1,maxatm)
      call jvalut (' Max nr of SGS operators = ',1,maxsym)
      call jvalut (' Max nr of NCS operators = ',1,maxncs)
c
      write (*,*)
      call prompt (' ==> SGS = Space Group Symmetry')
      call prompt (' ==> NCS = Non-Crystallographic Symmetry')
c
      linter = xinter()
      task = 'N'
      pdbin = 'in.pdb'
      pdbut = 'out.pdb'
      ncsut = 'xplor_ncs.include'
      cnsut = 'cns_ncs.def_blurb'
      rtout  = 'inverted.rt'
      spgrp = ' '
      symfil = 'p1.sym'
      ncsfil = ' '
      shell = 5.0
      do i=1,3
        cell(i) = 100.0
        cell(i+3) = 90.0
        fralo(i) = 0.0
        frahi(i) = 1.0
      end do
c
      j = ichar('A') - 1
      do i=1,maxncs
        j = j + 1
        if (j .gt. ichar('Z')) j = ichar('A')
        chains (i) = char (j)
      end do
c
   10 continue
      valid = 'NEPFAXCBWMORID?Q'
      write (*,'(20(1x,a,:,/))') ' ',
     +  'N = expand PDB file under NCS',
     +  'E = expand into a sphere under SGS',
     +  'P = expand PDB file under SGS',
     +  'F = expand fractional under SGS',
     +  'A = expand around a point under SGS',
     +  'X = generate complete X-PLOR (X)NCS file',
     +  'C = generate CNS (X)NCS blurb for "ncs.def" file',
     +  'B = generate non-bonded contacts (NCS+SGS)',
     +  'W = water scrutinizer',
     +  'M = generate PDB MTRIX records',
     +  'O = convert PDB MTRIX records to O-style RT operator(s)',
     +  'R = convert PDB REMARK 350 records to O-style RT operator(s)',
     +  'I = invert O-style RT operator(s)',
     +  'D = get default unit cell and spacegroup from PDB file',
     +  '? = print this list of options',
     +  'Q = quit'
c
   20 continue
      write (*,*)
      call textin (' Task ?',task)
      call upcase (task)
      if (index(valid,task) .le. 0) goto 10
      if (task .eq. '?') goto 10
      if (task .eq. 'Q') goto 999
c
c ... NCS EXPANSION
c     MTRIX RECORDS
c     INVERT RT
c
      if (task .eq. 'N' .or. task .eq. 'M' .or. task .eq. 'I') then
c
        if (task .eq. 'N') then
c
          call textin (' Input  PDB file ?',pdbin)
          call xopxoa (iunit,pdbin,linter,ierr)
          if (ierr .ne. 0) goto 1990
c
          call textin (' Output PDB file ?',pdbut)
          call xopxua (junit,pdbut,linter,ierr)
          if (ierr .ne. 0) goto 1990
c
        else if (task .eq. 'M') then
c
          call textin (' Output MTRIX file ?',pdbut)
          call xopxua (junit,pdbut,linter,ierr)
          if (ierr .ne. 0) goto 1990
c
        else if (task .eq. 'I') then
c
          call textin (' Output RT file ?',rtout)
          call xopxua (junit,rtout,linter,ierr)
          if (ierr .ne. 0) goto 1990
          close (junit)
c
        end if
c
        write (*,'(20(1x,a,:,/))') ' ',
     +    'Enter the names of the O-style NCS operator files',
     +    'including the unit operator; finish with <CR>'
c
        nncs = 0
  946   continue
        ncsfil = ' '
        call textin (' File with NCS operator(s) ?',ncsfil)
        if (length(ncsfil) .lt. 1) goto 1000
c
        call rdoncs (kunit,ncsfil,nncs,maxncs,rtncs,ierr)
        if (ierr .ne. 0) then
          call errcon (' While reading NCS operator file')
          goto 1990
        end if
c
        goto 946
c
 1000   continue
        write (*,*)
        if (nncs .lt. 1) then
          call errcon ('No NCS operators provided')
          goto 1990
        end if
c
        if (task .eq. 'N') then
          askme = 'A'
          call prompt (' Select:')
          call prompt (' A = automatic generation of chain IDs')
          call prompt (' M = enter chain IDs manually')
          call textin (' Option (A/M) ?',askme)
          call upcase (askme)
          if (askme .ne. 'M') askme = 'A'
          if (askme .eq. 'M') then
            do i=1,nncs
              write (line,'(1x,a,i3,a)') 'Chain ID for operator ',
     +          i,' ?'
              call textin (line,chains(i))
              call upcase (chains(i))
            end do
          else
            j = ichar('A') - 1
            do i=1,nncs
              j = j + 1
              if (j .gt. ichar('Z')) j = ichar('A')
              chains (i) = char (j)
            end do
          end if
        end if
c
cc        call anancs (nncs,rtncs,.true.,ierr)
cc        if (ierr .ne. 0) then
cc          call errcon ('In NCS operators')
cc          goto 1990
cc        end if
c
        if (task .eq. 'N') then
          call expncs (ierr)
        else if (task .eq. 'M') then
          call gmtrix (ierr)
        else if (task .eq. 'I') then
          call invncs (nncs,rtncs,rtinv)
          call wroncs (junit,rtout,nncs,maxncs,rtinv,ierr)
        end if
c
        if (ierr .ne. 0) goto 1990
c
        goto 20
c
c ... error
c
 1990   continue
        call errcon ('Sorry !')
        close (iunit)
        close (junit)
        close (kunit)
        goto 20
c
c ... CONVERT MTRIX
c     CONVERT REMARK 350
c
      else if (task .eq. 'O' .or. task .eq. 'R') then
c
        call textin (' Input PDB file ?',pdbin)
        call xopxoa (iunit,pdbin,linter,ierr)
        if (ierr .ne. 0) goto 1990
c
        ncsfil = ' '
        call textin (' Output file with NCS operator(s) ?',ncsfil)
        if (length(ncsfil) .lt. 1) goto 1000
        call xopxua (junit,ncsfil,linter,ierr)
c
        call pdb2o (task,ierr)
c
        if (ierr .ne. 0) then
          call prompt (' Error occurred during conversion')
        end if
c
        close (iunit)
        close (junit)
c
        goto 20
c
c ... GET DEFAULTS
c
      else if (task .eq. 'D') then
c
        call textin (' Input  PDB file ?',pdbin)
        call xopxoa (iunit,pdbin,linter,ierr)
        if (ierr .ne. 0) goto 1990
c
        call symcel (spgrp,ierr)
        close (iunit)
        if (ierr .ne. 0) goto 1990
c
        call fvalut (' Cell axes (A)     :',3,cell(1))
        call fvalut (' Cell angles (deg) :',3,cell(4))
c
        if (length(spgrp) .gt. 0) then
          call textut (' Spacegroup name   :',spgrp)
          call remspa (spgrp)
          call locase (spgrp)
          symfil = spgrp(1:length(spgrp)) // '.sym'
          call textut (' SGS operator file :',symfil)
        else
          call prompt (' No spacegroup name found')
        end if
c
        goto 20
c
c ... WATER SCRUTINIZER
c
      else if (task .eq. 'W') then
c
        call textin (' Input PDB file ?',pdbin)
        call xopxoa (iunit,pdbin,linter,ierr)
        if (ierr .ne. 0) goto 2990
c
        line = ' '
        call textin (' Log file (RETURN for screen output) ?',line)
        if (length(line) .lt. 1) then
          lunit = 6
        else
          lunit = 21
          call xopxua (lunit,line,linter,ierr)
          if (ierr .ne. 0) goto 2990
        end if
c
        line = ' '
        call textin (' Output PDB file (RETURN to skip) ?',line)
        if (length(line) .lt. 1) then
          munit = -1
        else
          munit = 22
          call xopxua (munit,line,linter,ierr)
          if (ierr .ne. 0) goto 2990
        end if
c
        call waters (lunit,munit)
c
        if (lunit .ne. 6) close (lunit)
        if (munit .gt. 0) close (munit)
c
        goto 20
c
c ... SPHERE EXPANSION
c     P1 EXPANSION
c     FRACTIONAL EXPANSION
c     AROUND-A-POINT EXPANSION
c
      else if (task .eq. 'E' .or. task .eq. 'P' .or.
     +         task .eq. 'F' .or. task .eq. 'A') then
c
        call textin (' Input  PDB file ?',pdbin)
        call xopxoa (iunit,pdbin,linter,ierr)
        if (ierr .ne. 0) goto 2990
c
        if (task .eq. 'E' .or. task .eq. 'P' .or.
     +      task .eq. 'F' .or. task .eq. 'A') then
          call textin (' Output PDB file ?',pdbut)
          call xopxua (junit,pdbut,linter,ierr)
          if (ierr .ne. 0) goto 2990
        end if
c
        call textin (' SGS operator file ?',symfil)
        call osymop (kunit,symfil,ierr)
        if (ierr .ne. 0) then
          call errcon ( 'While opening O datablock file')
          goto 2990
        end if
        close (kunit)
c
        call opoodb (kunit,symfil,par,partyp,j,fmt,ierr)
        call upcase (partyp)
        nsym = j/12
        i = nint (12.0 * (float(j)/12.0))
        if (ierr.ne.0 .or. i.ne.j .or. partyp.ne.'R') then
          call errcon ('In O datablock file')
          goto 2990
        end if
        call jvalut (' Nr of SGS operators :',1,nsym)
        if (nsym .gt. maxsym) then
          call errcon ('Too many SGS operators')
          goto 2990
        end if
        read (kunit,fmt,err=2990,end=2990)
     +    ((rtsym(j,i),j=1,12),i=1,nsym)
        close (kunit)
        do j=1,nsym
          call fratra (rtsym(10,j))
        end do
c
        call anasgs (nsym,rtsym,.true.,ierr)
        if (ierr .ne. 0) then
          call errcon ('In SGS operators')
          goto 2990
        end if
c
        call fvalin (' Cell constants ?',6,cell)
        call orthog (cell,fra2ca,0)
        call orthog (cell,ca2fra,1)
c
        if (task .eq. 'A') then
          call fvalin (' Fractional point to expand around ?',3,fpoint)
          call mulmtx (fra2ca,fpoint,cpoint,3,3,1)
          call fvalin (' Cartesian point to expand around ?',3,cpoint)
          call mulmtx (ca2fra,cpoint,fpoint,3,3,1)
          call fvalut (' Fractional :',3,fpoint)
          call fvalut (' Cartesian  :',3,cpoint)
          if (fpoint(1) .lt. -1.0 .or. fpoint(1) .gt. 1.0 .or.
     +        fpoint(2) .lt. -1.0 .or. fpoint(2) .gt. 1.0 .or.
     +        fpoint(3) .lt. -1.0 .or. fpoint(3) .gt. 1.0) then
            call errcon ('Fractional point not in range [-1,+1]')
            goto 2990
          end if
        end if
c
        if (task .eq. 'E') then
          call fvalin (' Shell radius (A) ?',1,shell)
          if (shell .le. 0.1 .or. shell .gt. 99.9) then
            call errcon ('Inappropriate shell radius')
            goto 2990
          end if
        end if
c
        if (task .eq. 'F') then
c
          askme = 'M'
          call prompt (' Select option:')
          call prompt (' M = apply to the Molecule as a whole')
          call prompt (' A = apply to each Atom individually')
          call textin (' Option (M/A) ?',askme)
          call upcase (askme)
          if (askme .ne. 'A') askme = 'M'
c
          if (askme .eq. 'M') then
            call prompt (' Max allowed fractional range [0,+1>')
          else
            call prompt (' Max allowed fractional range [-1,+2>')
          end if
          call fvalin (' Fractional lower limits ?',3,fralo)
          call fvalin (' Fractional upper limits ?',3,frahi)
          do i=1,3
            call rlohi(fralo(i),frahi(i))
            if (fralo(i) .eq. frahi(i)) then
              call errcon ('Empty fractional range')
              goto 2990
            end if
            if (askme .eq. 'M') then
              if (fralo(i) .lt. 0.0 .or. frahi(i) .gt. 1.0) then
                call errcon ('Outside allowed fractional range')
                goto 2990
              end if
            else
              if (fralo(i) .lt. -1.0 .or. frahi(i) .gt. 2.0) then
                call errcon ('Outside allowed fractional range')
                goto 2990
              end if
            end if
          end do
        end if
c
        if (task .eq. 'E') then
          call expsph (ierr)
          if (ierr .ne. 0) goto 2990
        else if (task .eq. 'P') then
          call expp1 (ierr)
          if (ierr .ne. 0) goto 2990
        else if (task .eq. 'F') then
          if (askme .eq. 'M') then
            call expfra (ierr)
          else
            call expatm (ierr)
          end if
          if (ierr .ne. 0) goto 2990
        else if (task .eq. 'A') then
          call expapt (ierr)
          if (ierr .ne. 0) goto 2990
        end if
c
        goto 20
c
c ... error
c
 2990   continue
        call errcon ('Sorry !')
        close (iunit)
        close (junit)
        close (kunit)
        goto 20
c
c ... X-PLOR NCS INCLUDE FILE GENERATION
c ... BONDAGE
c
      else if (task .eq. 'X' .or.
     +         task .eq. 'C' .or.
     +         task .eq. 'B') then
c
        call textin (' Input  PDB file ?',pdbin)
        call xopxoa (iunit,pdbin,linter,ierr)
        if (ierr .ne. 0) goto 2990
c
        if (task .eq. 'X') then
          call textin (' Output NCS file ?',ncsut)
          call xopxua (junit,ncsut,linter,ierr)
          if (ierr .ne. 0) goto 2990
        else if (task .eq. 'C') then
          call textin (' Output NCS file ?',cnsut)
          call xopxua (junit,cnsut,linter,ierr)
          if (ierr .ne. 0) goto 2990
        end if
c
        write (*,'(20(1x,a,:,/))') ' ',
     +    'Enter the names of the O-style NCS operator files',
     +    'including the unit operator; finish with <CR>'
c
        nncs = 0
 1946   continue
        ncsfil = ' '
        call textin (' File with NCS operator(s) ?',ncsfil)
        if (length(ncsfil) .lt. 1) goto 3000
c
        call rdoncs (kunit,ncsfil,nncs,maxncs,rtncs,ierr)
        if (ierr .ne. 0) then
          call errcon (' While reading NCS operator file')
          goto 1990
        end if
c
        goto 1946
c
 3000   continue
        write (*,*)
        call anancs (nncs,rtncs,.true.,ierr)
        if (ierr .ne. 0) then
          call errcon ('In NCS operators')
          goto 3990
        end if
c
        write (*,*)
        call textin (' SGS operator file ?',symfil)
        call osymop (kunit,symfil,ierr)
        if (ierr .ne. 0) then
          call errcon ( 'While opening O datablock file')
          goto 3990
        end if
        close (kunit)
c
        call opoodb (kunit,symfil,par,partyp,j,fmt,ierr)
        call upcase (partyp)
        nsym = j/12
        i = nint (12.0 * (float(j)/12.0))
        if (ierr.ne.0 .or. i.ne.j .or. partyp.ne.'R') then
          call errcon ('In O datablock file')
          goto 3990
        end if
        call jvalut (' Nr of SGS operators :',1,nsym)
        if (nsym .gt. maxsym) then
          call errcon ('Too many SGS operators')
          goto 3990
        end if
        read (kunit,fmt,err=3990,end=3990)
     +    ((rtsym(j,i),j=1,12),i=1,nsym)
        close (kunit)
        do j=1,nsym
          call fratra (rtsym(10,j))
        end do
c
        call anasgs (nsym,rtsym,.true.,ierr)
        if (ierr .ne. 0) then
          call errcon ('In SGS operators')
          goto 3990
        end if
c
        write (*,*)
        call fvalin (' Cell constants ?',6,cell)
        call orthog (cell,fra2ca,0)
        call orthog (cell,ca2fra,1)
c
        write (*,*)
        if (task .eq. 'X') then
          call fvalin (' Shell radius (A) ?',1,shell)
          if (shell .lt. 0.0 .or. shell .gt. 99.9) then
            call errcon ('Inappropriate contact distance')
            goto 3990
          end if
          call explor (ierr)
          if (ierr .ne. 0) goto 3990
        else if (task .eq. 'C') then
          call fvalin (' Shell radius (A) ?',1,shell)
          if (shell .lt. 0.0 .or. shell .gt. 99.9) then
            call errcon ('Inappropriate contact distance')
            goto 3990
          end if
          call ncscns (ierr)
          if (ierr .ne. 0) goto 3990
        else if (task .eq. 'B') then
          call prompt (
     +    ' Enter maximum distance for NCS+SGS interactions')
          call prompt (
     +    ' E.g., 2.4 A for bad contacts, 3.6 A for hydrogen')
          call prompt (
     +    ' bonds and salt links, or 4.5 A to include')
          call prompt (
     +    ' hydrophobic interactions')
          call fvalin (' Max contact distance (A) ?',1,shell)
          if (shell .lt. 0.0 .or. shell .gt. 99.9) then
            call errcon ('Inappropriate contact distance')
            goto 3990
          end if
          call bondag (ierr)
          if (ierr .ne. 0) goto 3990
        end if
c
        goto 20
c
c ... error
c
 3990   continue
        call errcon ('Sorry !')
        close (iunit)
        close (junit)
        close (kunit)
        goto 20
c
c ... invalid option
c
      else
        call errcon ('Invalid task')
        goto 10
      end if
c
      goto 20
c
c ... end of program
c
  999 continue
      call gkquit ()
c
      end
c
c
c
      subroutine expncs (ierr)
c
      include 'xpand.incl'
c
      real xut(3)
c
      integer nut,length,i,j,k,iat,ierr
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .false.
      lrem = .true.
      lxyz = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
 6000 format ('REMARK ',a,i2)
 6010 format ('REMARK X" = ',f10.6,' * X + ',f10.6,' * Y + ',
     +  f10.6,' * Z + ',f10.3)
 6020 format ('REMARK Y" = ',f10.6,' * X + ',f10.6,' * Y + ',
     +  f10.6,' * Z + ',f10.3)
 6030 format ('REMARK Z" = ',f10.6,' * X + ',f10.6,' * Y + ',
     +  f10.6,' * Z + ',f10.3)
c
c ... now do the NCS expansion
c
      iat = 0
      write (*,*)
c
 6110 format (' NCS expansion for operator # ',i6,
     +        ' with chain ID |',a1,'|')
c
      do i=1,nncs
        write (*,6110) i,chains(i)
c
ccc        call ivalut (' NCS expansion for operator :',1,i)
c
        write (junit,6000,err=9910)
        write (junit,6000,err=9910) 'NCS Operator nr ',i
        write (junit,6000,err=9910) 'Chain ID '//chains(i)
        write (junit,6010,err=9910) (rtncs(k,i),k=1,10,3)
        write (junit,6020,err=9910) (rtncs(k,i),k=2,11,3)
        write (junit,6030,err=9910) (rtncs(k,i),k=3,12,3)
        write (junit,6000,err=9910)
        nut = nut + 7
c
        do j=1,natoms
          iat = iat + 1
          line = atline(j)
          call vecrtv (atmxyz(1,j),xut,1,rtncs(1,i),rtncs(10,i))
          write (line(31:54),'(3f8.3)',err=9920) (xut(k),k=1,3)
          write (line(7:11),'(i5)',err=9920) iat
          line (22:22) = chains(i)
          write (junit,'(a)',err=9910) line(1:length(line))
          nut = nut + 1
        end do
c
      end do
c
      write (junit,'(a)',err=9910) 'END   '
      nut = nut + 1
c
      write (*,*)
      call jvalut (' Nr of atoms generated :',1,iat)
      call jvalut (' Nr of lines written   :',1,nut)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      goto 9999
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine expsph (ierr)
c
      include 'xpand.incl'
c
      real cogxyz(3),minxyz(3),maxxyz(3),radmax
      real cogfra(3),x1(3),x2(3),x3(3),d2max
      real dis,xeff(3),y1(3),y2(3),d3max
c
      integer nut,length,i,j,k,iat,ierr,k1,k2,k3,inew
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .false.
      lrem = .true.
      lxyz = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
      call getcog (cogxyz)
c
      call getlim (cogxyz,cogfra,minxyz,maxxyz,radmax,d2max,d3max)
c
 6000 format ('REMARK ',a,i2)
c
      iat = 0
c
c ... loop over the symmetry operators
c
      do i=1,nsym
c
        write (*,*)
        call ivalut (' Checking operator :',1,i)
c
c ... apply operator to FRACTIONAL centre-of-gravity
c
        call vecrtv (cogfra,x2,1,rtsym(1,i),rtsym(10,i))
c
c ... loop over 7*7*7 unit cells
c
        do k1=-3,3
          do k2=-3,3
            do k3=-3,3
c
c ... skip the Identity operator
c
              if (i .eq. 1) then
                if (k1.eq.0.and.k2.eq.0.and.k3.eq.0) goto 6969
              end if
c
c ... translate molecule
c
              x1(1) = x2(1) + float(k1)
              x1(2) = x2(2) + float(k2)
              x1(3) = x2(3) + float(k3)
c
c ... get CoG in Cartesian coordinates
c
              call mulmtx (fra2ca,x1,x3,3,3,1)
c
              dis = (cogxyz(1)-x3(1))**2 +
     +              (cogxyz(2)-x3(2))**2 +
     +              (cogxyz(3)-x3(3))**2
c
c ... check if close enough
c
              if (dis .le. d2max) then
c
                write (*,*)
                write (*,*) 'Translation ',k1,k2,k3
                call fvalut (' Distance :',1,sqrt(dis))
                inew = 0
c
c ... generate the atoms explicitly and check distance
c
                xeff(1) = float(k1)
                xeff(2) = float(k2)
                xeff(3) = float(k3)
c
                do j=1,natoms
c
c ... fractionalise, apply operator & translation, orthogonalise
c
                  call mulmtx (ca2fra,atmxyz(1,j),y2,3,3,1)
                  call vecrtv (y2,y1,1,rtsym(1,i),rtsym(10,i))
                  do k=1,3
                    y1(k)=y1(k)+xeff(k)
                  end do
                  call mulmtx (fra2ca,y1,y2,3,3,1)
c
                  dis = (cogxyz(1)-y2(1))**2 +
     +                  (cogxyz(2)-y2(2))**2 +
     +                  (cogxyz(3)-y2(3))**2
c
                  if (dis .le. d3max) then
                    iat = iat + 1
                    inew = inew + 1
c
                    if (inew .eq. 1) then
                      write (junit,6000,err=9910)
                      write (junit,6000,err=9910)
     +                  'SGS Operator nr ',i
                      write (junit,6000,err=9910)
     +                  'Translation A ',nint(xeff(1))
                      write (junit,6000,err=9910)
     +                  'Translation B ',nint(xeff(2))
                      write (junit,6000,err=9910)
     +                  'Translation C ',nint(xeff(3))
                      write (junit,6000,err=9910)
                      nut = nut + 6
                    end if
c
                    line = atline(j)
                    write (line(31:54),'(3f8.3)',err=9920)
     +                (y2(k),k=1,3)
c
                    write (junit,'(a)',err=9910)
     +                line(1:length(line))
                    nut = nut + 1
                  end if
                end do
c
                call jvalut (' Nr of atoms found :',1,inew)
c
              end if
c
 6969         continue
c
            end do
          end do
        end do
c
      end do
c
      write (junit,'(a)',err=9910) 'END   '
      nut = nut + 1
c
      write (*,*)
      call jvalut (' Nr of atoms written :',1,iat)
      call jvalut (' Nr of lines written :',1,nut)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      goto 9999
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine expp1 (ierr)
c
      include 'xpand.incl'
c
      real y1(3),y2(3)
c
      integer nut,length,i,j,k,iat,ierr
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .true.
      lxyz = .true.
      lrem = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
 6000 format ('REMARK ',a,i2)
c
      iat = natoms
c
c ... loop over the symmetry operators
c     (BUT SKIP THE IDENTITY OPERATOR (NR 1) !!!)
c
      write (*,*)
      call prompt (' NOTE - Assuming that the first symm-op')
      call prompt ('        is the identity operator !!!')
c
      do i=2,nsym
c
        write (*,*)
        call ivalut (' Generating operator :',1,i)
c
        write (junit,6000,err=9910)
        write (junit,6000,err=9910) 'SGS Operator nr ',i
        write (junit,6000,err=9910)
        nut = nut + 3
c
        do j=1,natoms
c
c ... fractionalise, apply operator, orthogonalise
c
          call mulmtx (ca2fra,atmxyz(1,j),y2,3,3,1)
          call vecrtv (y2,y1,1,rtsym(1,i),rtsym(10,i))
          call mulmtx (fra2ca,y1,y2,3,3,1)
c
c ... write ATOM/HETATM card
c
          iat = iat + 1
          line = atline(j)
          write (line(31:54),'(3f8.3)',err=9920) (y2(k),k=1,3)
          write (junit,'(a)',err=9910) line(1:length(line))
          nut = nut + 1
c
        end do
c
      end do
c
      write (junit,'(a)',err=9910) 'END   '
      nut = nut + 1
c
      write (*,*)
      call jvalut (' Nr of atoms written :',1,iat)
      call jvalut (' Nr of lines written :',1,nut)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      goto 9999
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine getpdb (lrem,lecho,lxyz,nut,ierr)
c
c ... GETPDB
c
c     if LREM  then echo REMARK etc. cards to output file
c     if LECHO then echo ATOM/HETATM cards to output file
c     if LXYZ  then store XYZ coordinates
c
      include 'xpand.incl'
c
      integer nin,nut,length,k,ierr
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nin = 0
      nut = 0
      ierr = 0
c
      natoms = 0
c
      write (*,*)
      call prompt (' Reading PDB file ...')
c
      if (lrem) then
        call stamp (line)
        write (junit,6000,err=9910) line(1:length(line))
      end if
c
 6000 format ('REMARK ',a,i2)
c
   10 continue
      read (iunit,'(a)',err=9900,end=1000) line
      nin = nin + 1
      call upcase (line)
      if (line(1:6) .eq. 'ATOM  ') goto 20
      if (line(1:6) .eq. 'HETATM') goto 20
      if (line(1:6) .eq. 'END   ') goto 1000
c
c ... HEADER, REMARK etc. cards
c
      if (lrem) then
        write (junit,'(a)',err=9910) line(1:length(line))
        nut = nut + 1
      end if
c
      goto 10
c
c ... ATOM or HETATM card
c
   20 continue
      natoms = natoms + 1
      if (natoms .gt. maxatm) then
        call errcon ('Too many atoms; rest skipped')
        natoms = maxatm
        goto 1000
      end if
      atline (natoms) = line
c
      if (lxyz) then
        read (line(31:54),'(3f8.3)',err=9920)
     +    (atmxyz(k,natoms),k=1,3)
      end if
c
      if (lecho) then
        write (junit,'(a)',err=9910) line(1:length(line))
        nut = nut + 1
      end if
c
      goto 10
c
c ... end of file (or too many atoms)
c
 1000 continue
      write (*,*)
      call jvalut (' Nr of lines read :',1,nin)
      call jvalut (' Nr of atoms read :',1,natoms)
      if (natoms .lt. 1) goto 9930
c
      return
c
c ... errors
c
 9900 continue
      write (*,*)
      call errcon ('While reading PDB file')
      write (*,*)
      ierr = -1
      return
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      return
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      return
c
 9930 continue
      write (*,*)
      call errcon ('No atoms found')
      write (*,*)
      ierr = -1
      return
c
      end
c
c
c
      subroutine getcog (cogxyz)
c
      include 'xpand.incl'
c
      real cogxyz(3)
c
      integer i,j
c
code ...
c
      do i=1,3
        cogxyz (i) = 0.0
      end do
c
      do j=1,natoms
        do i=1,3
          cogxyz (i) = cogxyz (i) + atmxyz (i,j)
        end do
      end do
c
      do i=1,3
        cogxyz(i) = cogxyz(i) / float(natoms)
      end do
c
      return
      end
c
c
c
      subroutine getlim (cogxyz,cogfra,minxyz,maxxyz,
     +  radmax,cogmx2,atmmx2)
c
      include 'xpand.incl'
c
      real cogxyz(3),minxyz(3),maxxyz(3),radmax,cogmx2,atmmx2
      real rad,dismax,cogfra(3)
c
      integer i,j,ifar
c
code ...
c
      do i=1,3
        minxyz(i) = atmxyz(i,1)
        maxxyz(i) = atmxyz(i,1)
      end do
c
      radmax = 0.0
      ifar = 1
      do i=1,natoms
        rad = 0.0
        do j=1,3
          rad = rad + (cogxyz (j) - atmxyz (j,i))**2
          minxyz(j) = min (minxyz(j),atmxyz(j,i))
          maxxyz(j) = max (maxxyz(j),atmxyz(j,i))
        end do
        if (i .eq. 1 .or. rad .gt. radmax) then
          radmax = rad
          ifar = i
        end if
      end do
c
      write (*,*)
      call fvalut (' Min XYZ :',3,minxyz)
      call fvalut (' Max XYZ :',3,maxxyz)
      call fvalut (' Centre  :',3,cogxyz)
      call mulmtx (ca2fra,cogxyz,cogfra,3,3,1)
      call fvalut (' Fractl  :',3,cogfra)
      call textut (' Farthest atom :',atline(ifar)(7:54))
      call fvalut (' Distance :',1,sqrt(radmax))
c
c ... get max allowable distance between centres-of-gravity
c     and square (to avoid square roots)
c
      dismax = 2.0 * sqrt(radmax) + shell
      call fvalut (' Max centre distance :',1,dismax)
      cogmx2 = dismax*dismax
      atmmx2 = (sqrt(radmax)+shell)**2
c
      return
      end
c
c
c
      subroutine expfra (ierr)
c
      include 'xpand.incl'
c
      real cogxyz(3),minxyz(3),maxxyz(3),radmax
      real cogfra(3),d2max,d3max
      real xtra(3),y1(3),y2(3),x2(3)
c
      integer nut,length,i,j,k,iat,ierr
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .false.
      lxyz = .true.
      lrem = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
 6000 format ('REMARK ',a,i2)
c
      call getcog (cogxyz)
c
      call getlim (cogxyz,cogfra,minxyz,maxxyz,radmax,d2max,d3max)
c
      iat = 0
c
c ... loop over the symmetry operators
c     (INCLUDING THE IDENTITY OPERATOR (NR 1) !!!)
c
      do i=1,nsym
c
        write (*,*)
        call ivalut (' Checking operator :',1,i)
c
c ... apply operator to FRACTIONAL centre-of-gravity
c
        call vecrtv (cogfra,x2,1,rtsym(1,i),rtsym(10,i))
c
c ... force inside the unit cell at (0,0,0)/(1,1,1)
c
        do k=1,3
          xtra(k) = 0.0
  220     if (x2(k) .lt. fralo(k)) then
            x2(k) = x2(k) + 1.0
            xtra(k)=xtra(k)+1.0
            goto 220
          end if
  230     if (x2(k) .ge. frahi(k)) then
            x2(k) = x2(k) - 1.0
            xtra(k)=xtra(k)-1.0
            goto 230
          end if
          if (x2(k) .lt. fralo(k)) then
            call prompt (' Does not fit into fractional range')
            goto 1999 
          end if
        end do
c
        call fvalut (' Translation applied :',3,xtra)
c
        write (junit,6000,err=9910)
        write (junit,6000,err=9910) 'SGS Operator nr ',i
        write (junit,6000,err=9910) 'Translation A ',nint(xtra(1))
        write (junit,6000,err=9910) 'Translation B ',nint(xtra(2))
        write (junit,6000,err=9910) 'Translation C ',nint(xtra(3))
        write (junit,6000,err=9910)
        nut = nut + 6
c
        do j=1,natoms
c
c ... fractionalise, apply operator & translation, orthogonalise
c
          call mulmtx (ca2fra,atmxyz(1,j),y2,3,3,1)
          call vecrtv (y2,y1,1,rtsym(1,i),rtsym(10,i))
          do k=1,3
            y1(k)=y1(k)+xtra(k)
          end do
          call mulmtx (fra2ca,y1,y2,3,3,1)
c
c ... write ATOM/HETATM card
c
          iat = iat + 1
          line = atline(j)
          write (line(31:54),'(3f8.3)',err=9920) (y2(k),k=1,3)
          write (junit,'(a)',err=9910) line(1:length(line))
          nut = nut + 1
c
        end do
c
 1999   continue
      end do
c
      write (junit,'(a)',err=9910) 'END   '
      nut = nut + 1
c
      write (*,*)
      call jvalut (' Nr of atoms written :',1,iat)
      call jvalut (' Nr of lines written :',1,nut)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      goto 9999
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine explor (ierr)
c
      include 'xpand.incl'
c
      integer maxinv
      parameter (maxinv=1000)
c
      real cogxyz(3),minxyz(3),maxxyz(3),radmax
      real cogfra(3),d2max,d3max,dis,rms
      real x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),rtx(12)
      real y2(3),y3(3),y4(3),y5(3),z1(3),z2(3),z3(3)
      real dum1(3,5),dum2(3,5)
      real rtinv(12,maxinv),rtnorm(12,maxinv)
c
      integer nut,length,i,j,k,ierr,inew,iii,jjj,k1,k2,k3
      integer itry,ninv
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      ierr = 0
      lecho = .false.
      lrem = .false.
      lxyz = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
      call stamp (line)
      write (junit,'(a6,1x,a)',err=9910) 'REMARK',line(1:length(line))
c
      write (junit,*,err=9910)
      write (junit,*,err=9910)
     +  '{ Invoke strict non-crystallographic symmetry }'
      write (junit,*,err=9910)
      write (junit,*,err=9910)
     +  'ncs strict'
      write (junit,*,err=9910)
      write (junit,*,err=9910)
     +  '  { ==> Assuming identity skew matrix }'
      write (junit,*,err=9910)
      write (junit,*,err=9910)
     +  '  skew'
      write (junit,*,err=9910)
     +  '    matrix = ( 1.000000 0.000000 0.000000 )'
      write (junit,*,err=9910)
     +  '             ( 0.000000 1.000000 0.000000 )'
      write (junit,*,err=9910)
     +  '             ( 0.000000 0.000000 1.000000 )'
      write (junit,*,err=9910)
     +  '    translation = ( 0.0000 0.0000 0.0000 )'
      write (junit,*,err=9910)
     +  '  end'
      write (junit,*,err=9910)
c
      do i=1,nncs
        write (junit,'(1x,a,i3,a)',err=9910)
     +    '  xncsrel { #',i,' }'
        write (junit,'(1x,a,3f13.6,a)',err=9910)
     +    '    matrix      = ( ',
     +    rtncs(1,i),rtncs(4,i),rtncs(7,i),' )'
        write (junit,'(1x,a,3f13.6,a)',err=9910)
     +    '                = ( ',
     +    rtncs(2,i),rtncs(5,i),rtncs(8,i),' )'
        write (junit,'(1x,a,3f13.6,a)',err=9910)
     +    '                = ( ',
     +    rtncs(3,i),rtncs(6,i),rtncs(9,i),' )'
        write (junit,'(1x,a,3f13.4,a)',err=9910)
     +    '    translation = ( ',
     +    rtncs(10,i),rtncs(11,i),rtncs(12,i),' )'
        write (junit,'(1x,a,i3,a)',err=9910)
     +    '  end { xncsrel  #',i,' }'
        write (junit,*,err=9910)
      end do
      write (junit,*,err=9910)
     +  '  { xncsrel end }'
      write (junit,*,err=9910)
c
c ... set up inverse operators
c
      ninv = 0
      do i=1,nncs
        ninv = ninv + 1
        do j=1,12
          rtnorm (j,ninv) = rtncs (j,i)
        end do
        do j=1,9
          rtinv (j,ninv) = rtncs (j,i)
        end do
        call matinv (rtinv(1,ninv), 3, z1, z2, z3)
        call mulmtx (rtinv(1,ninv),rtncs(10,i),rtinv(10,ninv),3,3,1)
        do j=10,12
          rtinv (j,ninv) = - rtinv (j,i)
        end do
      end do
c
      call getcog (cogxyz)
c
      call getlim (cogxyz,cogfra,minxyz,maxxyz,radmax,d2max,d3max)
c
      inew = 0
      itry = 0
c
c ... LOOP OVER ALL MOLECULES (NCS)
c
      do i=1,nncs
c
c ... get CoG of this (NCS) molecule
c
        call vecrtv (cogxyz,x1,1,rtncs(1,i),rtncs(10,i))
c
ccc        print *
ccc        print *,' TRY NCS # ',i
ccc        call fvalut (' CoG :',3,x1)
c
c ... get Cartesian coords of first five atoms
c
        do iii=1,5
          call vecrtv (atmxyz(1,iii),dum1(1,iii),1,
     +      rtncs(1,i),rtncs(10,i))
c
ccc          line = atline(iii)
ccc          write (line(31:54),'(3f8.3)',err=9920) (dum1(k,iii),k=1,3)
ccc          write (kunit,'(a)',err=9910) line(1:length(line))
c
        end do
c
c ... now loop over NCS-related molecules
c
        do j=1,nncs
c
c ... get fractional CoG of this NCS molecule
c
          call vecrtv (cogxyz,x2,1,rtncs(1,j),rtncs(10,j))
          call mulmtx (ca2fra,x2,x3,3,3,1)
c
ccc        print *,'   NOW NCS # ',j
ccc        call fvalut (' CoG :',3,x2)
ccc        call fvalut (' Fra :',3,x3)
c
c ... loop over SGS operators
c
          do k=1,nsym
c
c ... apply SGS operator to fractional CoG
c
            call vecrtv (x3,x4,1,rtsym(1,k),rtsym(10,k))
c
ccc        print *,'   NOW SGS # ',k
ccc        call fvalut (' Fra :',3,x4)
c
c ... loop over 7*7*7 unit cells
c
            do k1=-3,3
              do k2=-3,3
                do k3=-3,3
c
                  itry = itry + 1
c
c ... translate fractional CoG
c
                  x5(1) = x4(1) + float(k1)
                  x5(2) = x4(2) + float(k2)
                  x5(3) = x4(3) + float(k3)
c
c ... get CoG in Cartesian coordinates
c
                  call mulmtx (fra2ca,x5,x6,3,3,1)
c
                  dis = (x1(1)-x6(1))**2 +
     +                  (x1(2)-x6(2))**2 +
     +                  (x1(3)-x6(3))**2
c
ccc        print *,'   NOW TRA ',k1,k2,k3
ccc        call fvalut (' Fra :',3,x5)
ccc        call fvalut (' CoG :',3,x6)
ccc        print *,'   DIS ',sqrt(dis)
c
c ... check if close enough
c
                  if (dis .le. d2max .and. dis .gt. 1.0) then
c
                    write (*,*)
                    call fvalut (' Distance :',1,sqrt(dis))
c
ccc                    call fvalut (' NCS :',3,x1)
ccc                    call fvalut (' NEW :',3,x6)
c
c ... generate Cartesian coordinates of the first five atoms
c
                    do iii=1,5
                      call vecrtv (atmxyz(1,iii),y2,1,
     +                  rtncs(1,j),rtncs(10,j))
                      call mulmtx (ca2fra,y2,y3,3,3,1)
                      call vecrtv (y3,y4,1,rtsym(1,k),rtsym(10,k))
                      y5(1) = y4(1) + float(k1)
                      y5(2) = y4(2) + float(k2)
                      y5(3) = y4(3) + float(k3)
                      call mulmtx (fra2ca,y5,dum2(1,iii),3,3,1)
                    end do
c
c ... get the operator NCS(i) ---> T(k1,k2,k3) SGS(k) NCS(j)
c     use Wolfgang Kabsch's least-squares subroutine to do it
c
                    call lsqgjk (dum2,dum1,5,rms,rtx,ierr)
c
c ... check if it's not an existing (X)NCS operator
c
                    do iii=1,ninv
                      do jjj=1,9
                        if (abs(rtnorm(jjj,iii)-rtx(jjj)) .ge.
     +                      0.001) goto 6968
                      end do
                      do jjj=10,12
                        if (abs(rtnorm(jjj,iii)-rtx(jjj)) .ge.
     +                      0.01) goto 6968
                      end do
                      goto 6967
 6968                 continue
                    end do
                    goto 6969
c
c ... if here, NCS === XNCS
c
 6967               continue
                    call fvalut (' NCS = (X)NCS :',12,rtx)
                    goto 6970
c
c ... check if it's not the INVERSE of an existing (X)NCS operator
c
 6969               continue
                    do iii=1,ninv
                      do jjj=1,9
                        if (abs(rtinv(jjj,iii)-rtx(jjj)) .ge.
     +                      0.001) goto 6966
                      end do
                      do jjj=10,12
                        if (abs(rtinv(jjj,iii)-rtx(jjj)) .ge.
     +                      0.01) goto 6966
                      end do
                      goto 6967
 6966                 continue
                    end do
c
c ... new operator is okay and new; store its inverse
c
                    ninv = ninv + 1
                    do jjj=1,12
                      rtnorm (jjj,ninv) = rtx (jjj)
                    end do
                    do jjj=1,9
                      rtinv (jjj,ninv) = rtx (jjj)
                    end do
                    call matinv (rtinv(1,ninv), 3, z1, z2, z3)
                    call mulmtx (rtinv(1,ninv),rtx(10),
     +                rtinv(10,ninv),3,3,1)
                    do jjj=10,12
                      rtinv (jjj,ninv) = - rtinv (jjj,ninv)
                    end do
c
ccc        do iii=1,5
ccc          line = atline(iii)
ccc          write (line(31:54),'(3f8.3)',err=9920) (dum2(k,iii),k=1,3)
ccc          write (kunit,'(a)',err=9910) line(1:length(line))
ccc        end do
c
                    call fvalut (' New NCS op :',12,rtx)
                    write (*,'(1x,6(a,i3),a,f8.1,a)')
     +    'From NCS #',i,' -> NCS #',j,' SGS #',k,
     +    ' T=(',k1,',',k2,',',k3,') DIST = ',sqrt(dis),' A'
c
                    call anancs (1,rtx,.true.,ierr)
c
                    inew = inew + 1
c
                    write (junit,'(1x,a,i3,a)',err=9910)
     +    '  ncsrel { #',inew,' }'
                    write (junit,'(1x,6(a,i3),a,f8.1,a)',err=9910)
     +    '    { from NCS #',i,' -> NCS #',j,' SGS #',k,
     +    ' T=(',k1,',',k2,',',k3,') DIST = ',sqrt(dis),' A }'
                    write (junit,'(1x,a,3f13.6,a)',err=9910)
     +    '    matrix      = ( ',
     +    rtx(1),rtx(4),rtx(7),' )'
                    write (junit,'(1x,a,3f13.6,a)',err=9910)
     +    '                = ( ',
     +    rtx(2),rtx(5),rtx(8),' )'
                    write (junit,'(1x,a,3f13.6,a)',err=9910)
     +    '                = ( ',
     +    rtx(3),rtx(6),rtx(9),' )'
                    write (junit,'(1x,a,3f13.4,a)',err=9910)
     +    '    translation = ( ',
     +    rtx(10),rtx(11),rtx(12),' )'
                    write (junit,'(1x,a,i3,a)',err=9910)
     +    '  end { ncsrel  #',inew,' }'
                    write (junit,*,err=9910)
c
 6970               continue
c
                  end if
c
                end do
              end do
            end do
c
          end do
c
        end do
c
      end do
c
      write (junit,*,err=9910) '  { ncsrel end }'
      write (junit,*,err=9910)
      write (junit,*,err=9910) '  ?'
      write (junit,*,err=9910)
      write (junit,*,err=9910) 'end {ncs strict}'
c
      write (*,*)
      call jvalut (' Nr of NCSRel generated :',1,inew)
      call jvalut (' Nr of operators tested :',1,itry)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing NCS file')
      write (*,*)
      ierr = -1
      goto 9999
c
ccc 9920 continue
ccc      write (*,*)
ccc      call errcon ('During internal I/O operation')
ccc      write (*,*)
ccc      ierr = -1
ccc      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine ncscns (ierr)
c
      include 'xpand.incl'
c
      integer maxinv
      parameter (maxinv=1000)
c
      real cogxyz(3),minxyz(3),maxxyz(3),radmax
      real cogfra(3),d2max,d3max,dis,rms
      real x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),rtx(12)
      real y2(3),y3(3),y4(3),y5(3),z1(3),z2(3),z3(3)
      real dum1(3,5),dum2(3,5)
      real rtinv(12,maxinv),rtnorm(12,maxinv)
c
      integer nut,length,i,j,k,ierr,inew,iii,jjj,k1,k2,k3
      integer itry,ninv
c
      logical lecho,lxyz,lrem
c
      character line*256,dum*20
c
code ...
c
      ierr = 0
      lecho = .false.
      lrem = .false.
      lxyz = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
      call stamp (line)
      write (junit,'(/3a/)',err=9910)
     +  '{* ',line(1:length(line)),' *}'
c
      write (junit,'(a,a)',err=9910)
     +  '{=============== strict NCS (X-ray and',
     +  ' nonbonded terms) ==============}'
      write (junit,'(/a,a)',err=9910)
     +  '{* these operators are used to generate',
     +  ' the complete asymmetric unit for'
      write (junit,'(a)',err=9910)
     +  '   structure factor and nonbonded energy calculations *}'
c
 6000 format (/'{* NCS operator ',i4,' *}'//
     +  '{* use this operator *}'/
     +  '{+ choice: true false +}')
 6010 format ('{===>} ',a,'=true;'//
     +  '{* real-space rotation matrix *}')
 6020 format ('{===>} ',a20,' = (',3(1x,f13.6),' )'/
     + '                         (',3(1x,f13.6),' )'/
     + '                         (',3(1x,f13.6),' );')
 6030 format (/'{* real-space translation vector *}'/
     +  '{===>} ',a,' = (',3(1x,f13.4),' );')
c
      do i=1,nncs
        write (junit,6000,err=9910) i
        write (dum,'(a16,i4)') 'ncs_op_',i
        call remspa (dum)
        write (junit,6010,err=9910) dum(1:length(dum))
        write (dum,'(a16,i4)') 'ncs_matrix_',i
        call remspa (dum)
        write (junit,6020,err=9910) dum(1:length(dum)),
     +    rtncs(1,i),rtncs(4,i),rtncs(7,i),
     +    rtncs(2,i),rtncs(5,i),rtncs(8,i),
     +    rtncs(3,i),rtncs(6,i),rtncs(9,i)
        write (dum,'(a16,i4)') 'ncs_vector_',i
        call remspa (dum)
        write (junit,6030,err=9910) dum(1:length(dum)),
     +    rtncs(10,i),rtncs(11,i),rtncs(12,i)
      end do
      write (junit,*,err=9910)
c
c ... set up inverse operators
c
      ninv = 0
      do i=1,nncs
        ninv = ninv + 1
        do j=1,12
          rtnorm (j,ninv) = rtncs (j,i)
        end do
        do j=1,9
          rtinv (j,ninv) = rtncs (j,i)
        end do
        call matinv (rtinv(1,ninv), 3, z1, z2, z3)
        call mulmtx (rtinv(1,ninv),rtncs(10,i),rtinv(10,ninv),3,3,1)
        do j=10,12
          rtinv (j,ninv) = - rtinv (j,i)
        end do
      end do
c
      write (junit,'(/a,a)',err=9910)
     +  '{=============== strict NCS (nonbonded',
     +  ' terms only) ==============}'
      write (junit,'(/a/a)',err=9910)
     +  '{* these operators can be used to generate the full',
     +  '   nonbonded environment around the protomer *}'
c
 6100 format (/'{* NCS nonbonded operator ',i4,' *}'/)
 6105 format (/'{* use this operator *}'/
     +  '{+ choice: true false +}')
 6110 format ('{===>} ',a,'=true;'//
     +  '{* real-space rotation matrix *}')
 6120 format ('{===>} ',a20,' = (',3(1x,f13.6),' )'/
     + '                         (',3(1x,f13.6),' )'/
     + '                         (',3(1x,f13.6),' );')
 6130 format (/'{* real-space translation vector *}'/
     +  '{===>} ',a,' = (',3(1x,f13.4),' );')
c
      call getcog (cogxyz)
c
      call getlim (cogxyz,cogfra,minxyz,maxxyz,radmax,d2max,d3max)
c
      inew = 0
      itry = 0
c
c ... LOOP OVER ALL MOLECULES (NCS)
c
      do i=1,nncs
c
c ... get CoG of this (NCS) molecule
c
        call vecrtv (cogxyz,x1,1,rtncs(1,i),rtncs(10,i))
c
ccc        print *
ccc        print *,' TRY NCS # ',i
ccc        call fvalut (' CoG :',3,x1)
c
c ... get Cartesian coords of first five atoms
c
        do iii=1,5
          call vecrtv (atmxyz(1,iii),dum1(1,iii),1,
     +      rtncs(1,i),rtncs(10,i))
c
ccc          line = atline(iii)
ccc          write (line(31:54),'(3f8.3)',err=9920) (dum1(k,iii),k=1,3)
ccc          write (kunit,'(a)',err=9910) line(1:length(line))
c
        end do
c
c ... now loop over NCS-related molecules
c
        do j=1,nncs
c
c ... get fractional CoG of this NCS molecule
c
          call vecrtv (cogxyz,x2,1,rtncs(1,j),rtncs(10,j))
          call mulmtx (ca2fra,x2,x3,3,3,1)
c
ccc        print *,'   NOW NCS # ',j
ccc        call fvalut (' CoG :',3,x2)
ccc        call fvalut (' Fra :',3,x3)
c
c ... loop over SGS operators
c
          do k=1,nsym
c
c ... apply SGS operator to fractional CoG
c
            call vecrtv (x3,x4,1,rtsym(1,k),rtsym(10,k))
c
ccc        print *,'   NOW SGS # ',k
ccc        call fvalut (' Fra :',3,x4)
c
c ... loop over 7*7*7 unit cells
c
            do k1=-3,3
              do k2=-3,3
                do k3=-3,3
c
                  itry = itry + 1
c
c ... translate fractional CoG
c
                  x5(1) = x4(1) + float(k1)
                  x5(2) = x4(2) + float(k2)
                  x5(3) = x4(3) + float(k3)
c
c ... get CoG in Cartesian coordinates
c
                  call mulmtx (fra2ca,x5,x6,3,3,1)
c
                  dis = (x1(1)-x6(1))**2 +
     +                  (x1(2)-x6(2))**2 +
     +                  (x1(3)-x6(3))**2
c
ccc        print *,'   NOW TRA ',k1,k2,k3
ccc        call fvalut (' Fra :',3,x5)
ccc        call fvalut (' CoG :',3,x6)
ccc        print *,'   DIS ',sqrt(dis)
c
c ... check if close enough
c
                  if (dis .le. d2max .and. dis .gt. 1.0) then
c
                    write (*,*)
                    call fvalut (' Distance :',1,sqrt(dis))
c
ccc                    call fvalut (' NCS :',3,x1)
ccc                    call fvalut (' NEW :',3,x6)
c
c ... generate Cartesian coordinates of the first five atoms
c
                    do iii=1,5
                      call vecrtv (atmxyz(1,iii),y2,1,
     +                  rtncs(1,j),rtncs(10,j))
                      call mulmtx (ca2fra,y2,y3,3,3,1)
                      call vecrtv (y3,y4,1,rtsym(1,k),rtsym(10,k))
                      y5(1) = y4(1) + float(k1)
                      y5(2) = y4(2) + float(k2)
                      y5(3) = y4(3) + float(k3)
                      call mulmtx (fra2ca,y5,dum2(1,iii),3,3,1)
                    end do
c
c ... get the operator NCS(i) ---> T(k1,k2,k3) SGS(k) NCS(j)
c     use Wolfgang Kabsch's least-squares subroutine to do it
c
                    call lsqgjk (dum2,dum1,5,rms,rtx,ierr)
c
c ... check if it's not an existing (X)NCS operator
c
                    do iii=1,ninv
                      do jjj=1,9
                        if (abs(rtnorm(jjj,iii)-rtx(jjj)) .ge.
     +                      0.001) goto 6968
                      end do
                      do jjj=10,12
                        if (abs(rtnorm(jjj,iii)-rtx(jjj)) .ge.
     +                      0.01) goto 6968
                      end do
                      goto 6967
 6968                 continue
                    end do
                    goto 6969
c
c ... if here, NCS === XNCS
c
 6967               continue
                    call fvalut (' NCS = (X)NCS :',12,rtx)
                    goto 6970
c
c ... check if it's not the INVERSE of an existing (X)NCS operator
c
 6969               continue
                    do iii=1,ninv
                      do jjj=1,9
                        if (abs(rtinv(jjj,iii)-rtx(jjj)) .ge.
     +                      0.001) goto 6966
                      end do
                      do jjj=10,12
                        if (abs(rtinv(jjj,iii)-rtx(jjj)) .ge.
     +                      0.01) goto 6966
                      end do
                      goto 6967
 6966                 continue
                    end do
c
c ... new operator is okay and new; store its inverse
c
                    ninv = ninv + 1
                    do jjj=1,12
                      rtnorm (jjj,ninv) = rtx (jjj)
                    end do
                    do jjj=1,9
                      rtinv (jjj,ninv) = rtx (jjj)
                    end do
                    call matinv (rtinv(1,ninv), 3, z1, z2, z3)
                    call mulmtx (rtinv(1,ninv),rtx(10),
     +                rtinv(10,ninv),3,3,1)
                    do jjj=10,12
                      rtinv (jjj,ninv) = - rtinv (jjj,ninv)
                    end do
c
ccc        do iii=1,5
ccc          line = atline(iii)
ccc          write (line(31:54),'(3f8.3)',err=9920) (dum2(k,iii),k=1,3)
ccc          write (kunit,'(a)',err=9910) line(1:length(line))
ccc        end do
c
                    call fvalut (' New NCS op :',12,rtx)
                    write (*,'(1x,6(a,i3),a,f8.1,a)')
     +    'From NCS #',i,' -> NCS #',j,' SGS #',k,
     +    ' T=(',k1,',',k2,',',k3,') DIST = ',sqrt(dis),' A'
c
                    call anancs (1,rtx,.true.,ierr)
c
                    inew = inew + 1
c
                    write (junit,6100,err=9910) inew
                    write (junit,'(6(a,i3),a,f8.1,a)',err=9910)
     +    '{* from NCS #',i,' -> NCS #',j,' SGS #',k,
     +    ' T=(',k1,',',k2,',',k3,') DIST = ',sqrt(dis),' A *}'
                    write (junit,6105,err=9910)
                    write (dum,'(a16,i4)')
     +                'nb_ncs_op_',inew
                    call remspa (dum)
                    write (junit,6110,err=9910)
     +                dum(1:length(dum))
                    write (dum,'(a16,i4)')
     +                'nb_ncs_matrix_',inew
                    call remspa (dum)
                    write (junit,6120,err=9910)
     +                dum(1:length(dum)),
     +                rtx(1),rtx(4),rtx(7),
     +                rtx(2),rtx(5),rtx(8),
     +                rtx(3),rtx(6),rtx(9)
                    write (dum,'(a16,i4)')
     +                'nb_ncs_vector_',inew
                    call remspa (dum)
                    write (junit,6130,err=9910)
     +                dum(1:length(dum)),
     +                rtx(10),rtx(11),rtx(12)
c
 6970               continue
c
                  end if
c
                end do
              end do
            end do
c
          end do
c
        end do
c
      end do
      write (junit,*,err=9910)
c
      write (*,*)
      call jvalut (' Nr of NCSRel generated :',1,inew)
      call jvalut (' Nr of operators tested :',1,itry)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing NCS file')
      write (*,*)
      ierr = -1
      goto 9999
c
ccc 9920 continue
ccc      write (*,*)
ccc      call errcon ('During internal I/O operation')
ccc      write (*,*)
ccc      ierr = -1
ccc      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine bondag (ierr)
c
      include 'xpand.incl'
c
      integer maxinv
      parameter (maxinv=1000)
c
      real cogxyz(3),minxyz(3),maxxyz(3),radmax
      real cogfra(3),d2max,d3max,dis,sh2,adis,cutmax
      real x1(3),x2(3),x3(3),x4(3),x5(3),x6(3)
      real y2(3),y3(3),y4(3),y5(3)
c
      integer nut,i,j,k,ierr,iii,jjj,k1,k2,k3,icnt,ipart
      integer i1,i2,i3,iseek,idi,j1,j2,j3,itot,length
c
      logical lecho,lxyz,lrem
c
      character line*80
c
code ...
c
      ierr = 0
      lecho = .false.
      lrem = .false.
      lxyz = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) return
c
      call getcog (cogxyz)
c
      call getlim (cogxyz,cogfra,minxyz,maxxyz,radmax,d2max,d3max)
c
      call boxem (cogxyz,0)
c
      sh2 = shell*shell
      idi = int (shell/5.0) + 1
      itot = 0
      cutmax = (sqrt(radmax)+shell+1.0)**2
c
c ... LOOP OVER ALL MOLECULES (NCS)
c
      do i=1,nncs
c
c ... get CoG of this (NCS) molecule
c
        call vecrtv (cogxyz,x1,1,rtncs(1,i),rtncs(10,i))
c
c ... get Cartesian coords of NCS atoms
c
        do iii=1,natoms
          call vecrtv (atmxyz(1,iii),ncsxyz(1,iii),1,
     +      rtncs(1,i),rtncs(10,i))
        end do
c
        call boxem (cogxyz,1)
c
ccc        print *
ccc        print *,' TRY NCS # ',i
ccc        call fvalut (' CoG :',3,x1)
c
c ... now loop over NCS-related molecules
c
        do j=1,nncs
c
c ... get fractional CoG of this NCS molecule
c
          call vecrtv (cogxyz,x2,1,rtncs(1,j),rtncs(10,j))
          call mulmtx (ca2fra,x2,x3,3,3,1)
c
ccc        print *,'   NOW NCS # ',j
ccc        call fvalut (' CoG :',3,x2)
ccc        call fvalut (' Fra :',3,x3)
c
c ... loop over SGS operators
c
          do k=1,nsym
c
c ... apply SGS operator to fractional CoG
c
            call vecrtv (x3,x4,1,rtsym(1,k),rtsym(10,k))
c
ccc        print *,'   NOW SGS # ',k
ccc        call fvalut (' Fra :',3,x4)
c
c ... loop over 7*7*7 unit cells
c
            do k1=-3,3
              do k2=-3,3
                do k3=-3,3
c
c ... translate fractional CoG
c
                  x5(1) = x4(1) + float(k1)
                  x5(2) = x4(2) + float(k2)
                  x5(3) = x4(3) + float(k3)
c
c ... get CoG in Cartesian coordinates
c
                  call mulmtx (fra2ca,x5,x6,3,3,1)
c
                  dis = (x1(1)-x6(1))**2 +
     +                  (x1(2)-x6(2))**2 +
     +                  (x1(3)-x6(3))**2
c
ccc        print *,'   NOW TRA ',k1,k2,k3
ccc        call fvalut (' Fra :',3,x5)
ccc        call fvalut (' CoG :',3,x6)
ccc        print *,'   DIS ',sqrt(dis)
c
c ... check if close enough
c
                  if (dis .le. d2max .and. dis .gt. 1.0) then
c
                    write (*,*)
                    write (*,'(1x,6(a,i3),a,f8.1,a)')
     +    'From NCS #',i,' -> NCS #',j,' SGS #',k,
     +    ' T=(',k1,',',k2,',',k3,') DIST = ',sqrt(dis),' A'
c
                    if (i .eq. j) then
                      print *,'Pure SGS interaction'
                    else
                      print *,'Mixed NCS and SGS interaction'
                    end if
c
ccc      print *,'Generating symm coords'
                    do iii=1,natoms
                      call vecrtv (atmxyz(1,iii),y2,1,
     +                  rtncs(1,j),rtncs(10,j))
                      call mulmtx (ca2fra,y2,y3,3,3,1)
                      call vecrtv (y3,y4,1,rtsym(1,k),rtsym(10,k))
                      y5(1) = y4(1) + float(k1)
                      y5(2) = y4(2) + float(k2)
                      y5(3) = y4(3) + float(k3)
                      call mulmtx (fra2ca,y5,symxyz(1,iii),3,3,1)
                    end do
c
ccc      print *,'Boxing'
                    call boxem (cogxyz,2)
c
                    icnt = 0
c
                    do iii=1,natoms
c
c ... check if atom near symm-related molecule
c
                      adis = (ncsxyz(1,iii)-x6(1))**2 +
     +                       (ncsxyz(2,iii)-x6(2))**2 +
     +                       (ncsxyz(3,iii)-x6(3))**2
                      if (adis .gt. cutmax) goto 6969
c
                      call packut (j1,j2,j3,ipart,ncsbox(iii))
ccc      if (iii.eq.1) print *,'Atom 1 ',j1,j2,j3,ncsbox(iii),idi
                      do i1=j1-idi,j1+idi
                        do i2=j2-idi,j2+idi
                          do i3=j3-idi,j3+idi
                            call packin (i1,i2,i3,0,iseek)
                            do jjj=1,natoms
                              if (symbox(jjj).ne.iseek) goto 6970
ccc                              if (abs(ncsxyz(1,iii)-symxyz(1,jjj)) .gt.
ccc     +                            shell) goto 6970
ccc                              if (abs(ncsxyz(2,iii)-symxyz(2,jjj)) .gt.
ccc     +                            shell) goto 6970
ccc                              if (abs(ncsxyz(3,iii)-symxyz(3,jjj)) .gt.
ccc     +                            shell) goto 6970
                              adis = (ncsxyz(1,iii)-symxyz(1,jjj))**2 +
     +                               (ncsxyz(2,iii)-symxyz(2,jjj))**2 +
     +                               (ncsxyz(3,iii)-symxyz(3,jjj))**2
ccc      print *,'TRY ',jjj,adis
                              if (adis .gt. sh2) goto 6970
c
                              adis = sqrt (adis)
                              call intera(iii,jjj,adis,line)
                              write (*,'(4(1x,a),1x,f8.2,a,1x,a)')
     +                          atline(iii)(13:27),' ---> ',
     +                          atline(jjj)(13:27),' = ',
     +                          adis,' A',line(1:length(line))
                              icnt = icnt + 1
c
 6970                         continue
                            end do
                          end do
                        end do
                      end do
 6969                 continue
                    end do
c
                    itot = itot + icnt
                    call jvalut (' Nr of contacts :',1,icnt)
c
                  end if
c
                end do
              end do
            end do
c
          end do
c
        end do
c
      end do
c
      write (*,*)
      call jvalut (' Total nr of contacts :',1,itot)
      call prompt (' Every contact should occur twice !')
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine boxem (cog,jflag)
c
      include 'xpand.incl'
c
      integer boxsiz
      parameter (boxsiz = 5.0)
c
      real cog(3)
c
      integer jflag,iflag,i,i1,i2,i3
c
code ...
c
      iflag = max (0,min(2,jflag))
c
      if (iflag .eq. 0) then
        do i=1,natoms
          i1 = int ( (atmxyz(1,i)-cog(1)) / boxsiz )
          i2 = int ( (atmxyz(2,i)-cog(2)) / boxsiz )
          i3 = int ( (atmxyz(3,i)-cog(3)) / boxsiz )
          call packin (i1,i2,i3,0,atmbox(i))
        end do
      else if (iflag .eq. 1) then
        do i=1,natoms
          i1 = int ( (ncsxyz(1,i)-cog(1)) / boxsiz )
          i2 = int ( (ncsxyz(2,i)-cog(2)) / boxsiz )
          i3 = int ( (ncsxyz(3,i)-cog(3)) / boxsiz )
          call packin (i1,i2,i3,0,ncsbox(i))
        end do
      else if (iflag .eq. 2) then
        do i=1,natoms
          i1 = int ( (symxyz(1,i)-cog(1)) / boxsiz )
          i2 = int ( (symxyz(2,i)-cog(2)) / boxsiz )
          i3 = int ( (symxyz(3,i)-cog(3)) / boxsiz )
          call packin (i1,i2,i3,0,symbox(i))
        end do
      end if
c
      return
      end
c
c
c
      subroutine intera (i,j,dis,line)
c
      include 'xpand.incl'
c
      real dis,badcon,hydro1,hydro2,hbond
c
      integer i,j
c
      logical plus1,plus2,neg1,neg2
c
      character line*(*),at1*8,at2*8
c
code ...
c
      badcon = 2.4
      hbond = 3.6
      hydro1 = 3.5
      hydro2 = 4.5
c
      line = '?'
c
      at1 = atline(i)(13:20)
      at2 = atline(j)(13:20)
c
      if (dis .lt. badcon) then
        if (at1 .eq. ' SG  CYS' .and.
     +      at2 .eq. ' SG  CYS') then
          line = '==> inter-molecular SS-bond ?'
          goto 10
        else
          line = '==> BAD contact !'
          goto 10
        end if
      end if
c
      if (dis .gt. hydro2) then
        line = '-'
        goto 10
      end if
c
      if (at1(1:2) .ne. ' C' .and. at1(1:2) .ne. ' N' .and.
     +    at1(1:2) .ne. ' O' .and. at1(1:2) .ne. ' S') goto 10
c
      if (at2(1:2) .ne. ' C' .and. at2(1:2) .ne. ' N' .and.
     +    at2(1:2) .ne. ' O' .and. at2(1:2) .ne. ' S') goto 10
c
      if ( (at1(1:2) .eq. ' C' .or. at1(1:2) .eq. ' S') .and.
     +     (at2(1:2) .eq. ' C' .or. at2(1:2) .eq. ' S') ) then
        if (dis .lt. hydro1) then
          line = '==> (too ?) close contact !'
          goto 10
        else if (dis .gt. hydro2) then
          line = '-'
          goto 10
        else
          line = 'hydrophobic interaction'
          goto 10
        end if
      end if
c
      if ( (at1(1:2) .eq. ' N' .or. at1(1:2) .eq. ' O'
     +         .or. at1(1:2) .eq. ' S') .and.
     +     (at2(1:2) .eq. ' N' .or. at2(1:2) .eq. ' O'
     +         .or. at2(1:2) .eq. ' S') ) then
        if (dis .gt. hbond) then
          line = '-'
          goto 10
        end if
        plus1 = (at1.eq.' NH1 ARG'.or.at1.eq.' NH2 ARG'.or.
     +           at1.eq.' NZ  LYS'.or.at1.eq.' ND1 HIS'.or.
     +           at1.eq.' NE2 HIS')
        plus2 = (at2.eq.' NH1 ARG'.or.at2.eq.' NH2 ARG'.or.
     +           at2.eq.' NZ  LYS'.or.at2.eq.' ND1 HIS'.or.
     +           at2.eq.' NE2 HIS')
        neg1  = (at1.eq.' OD1 ASP'.or.at1.eq.' OD2 ASP'.or.
     +           at1.eq.' OE1 GLU'.or.at1.eq.' OE2 GLU')
        neg2  = (at2.eq.' OD1 ASP'.or.at2.eq.' OD2 ASP'.or.
     +           at2.eq.' OE1 GLU'.or.at2.eq.' OE2 GLU')
        if ( (plus1.and.neg2) .or. (plus2.and.neg1) ) then
          line = 'salt link'
          goto 10
        else if (plus1 .and. plus2) then
          line = '==> (too ?) close positive charges ?'
          goto 10
        else if (neg1 .and. neg2) then
          line = '==> (too ?) close negative charges ?'
          goto 10
        end if
c
        line = 'potential hydrogen bond ?'
        goto 10
      end if
c
      if (dis .le. hbond) then
        line = '==> (too ?) close contact !'
        goto 10
      end if
c
   10 continue
c
      if (atline(i)(18:27) .eq. atline(j)(18:27)) then
        line = 'SELF '//line
      end if
c
      return
      end
c
c
c
      subroutine waters (lunit,munit)
c
      include 'xpand.incl'
c
      real dist,x,chbo,ccha,cbad,ccok
c
      integer nut,i,j,ierr,nwat,npok,nmok,lunit,munit
      integer nplus,nminu,nhbdo,nhbac,ndoac,ncyss,nmets
      integer nnh4,nh2o,nf,nna,nmg,nnot,length
c
      logical lecho,lxyz,lrem,lnh4,lh2o,lf,lna,lmg,lhydro,lwater
c
      character rtyp*3,anam*4,rnam*6,rtyp2*3,anam2*4,rnam2*6,line*128
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .false.
      lrem = .false.
      lxyz = .true.
c
      nnh4 = 0
      nh2o = 0
      nf   = 0
      nna  = 0
      nmg  = 0
      nnot = 0
c
      chbo = 3.5
      ccha = 3.5
      ccok = 2.8
      cbad = 1.8
c
      write (*,*)
      call fvalut (' Max interaction distance :',1,chbo)
      call fvalut (' Min bad-charge  distance :',1,ccha)
      call fvalut (' Max good-charge distance :',1,ccok)
      call fvalut (' Min contact     distance :',1,cbad)
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      close (iunit)
      if (ierr .ne. 0) goto 9999
c
      if (munit .gt. 0) then
        call stamp (line)
        write (munit,'(a6,1x,a)',err=33) 'REMARK',line(1:length(line))
   33   continue
      end if
c
cATOM   1110  O1  HOH A 300       1.742  13.554   0.117  1.00 21.92
c1234567890123456789012345678901234567890123456789012345678901234567890
c
      nwat = 0
      do i=1,natoms
        anam = atline(i)(13:16)
        if (lhydro(anam)) goto 100
        rtyp = atline(i)(18:20)
        if (.not.lwater(rtyp)) goto 100
        rnam = atline(i)(21:26)
        nwat = nwat + 1
        write (lunit,6000) nwat,rtyp,rnam,anam
c
        nplus = 0
        npok  = 0
        nhbdo = 0
        nminu = 0
        nmok  = 0
        nhbac = 0
        ndoac = 0
        ncyss = 0
        nmets = 0
c
        lnh4 = .true.
        lh2o = .true.
        lf   = .true.
        lna  = .true.
        lmg  = .true.
c
        do j=1,natoms
          if (i .eq. j) goto 200
          anam2 =  atline(j)(13:16)
          if (lhydro(anam2)) goto 200
          rtyp2 = atline(j)(18:20)
          rnam2 = atline(j)(21:26)
c
          x = dist (i,j,atmxyz)
c
          if (x .le. chbo) then
c
            write (lunit,6010) rtyp2,rnam2,anam2,x
c
            if (anam2(1:2) .eq. ' C') then
              write (lunit,'(a)') '    BAD contact !'
              goto 200
            end if
c
            if (x .le. cbad) then
              write (lunit,'(a)') '    Too close (bad ?) contact !'
              goto 200
            end if
c
            if (index ('| N| O| S|',anam2(1:2)) .lt. 1) then
              write (lunit,'(a)') '    Can only handle C/N/O/S'
              goto 200
            end if
c
            if (rtyp2 .eq. 'ARG') then
              if (anam2 .eq. ' NH1' .or.
     +            anam2 .eq. ' NH2') then
                nplus = nplus + 1
                nhbdo = nhbdo + 1
                if (x .le. ccha) then
                  lnh4 = .false.
                  lna  = .false.
                  lmg  = .false.
                end if
                if (x .le. ccok) then
                  npok = npok + 1
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'LYS') then
              if (anam2 .eq. ' NE ') then
                nplus = nplus + 1
                nhbdo = nhbdo + 1
                if (x .le. ccha) then
                  lnh4 = .false.
                  lna  = .false.
                  lmg  = .false.
                end if
                if (x .le. ccok) then
                  npok = npok + 1
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'HIS') then
              if (anam2 .eq. ' ND1' .or.
     +            anam2 .eq. ' NE2') then
                nplus = nplus + 1
                nhbdo = nhbdo + 1
                if (x .le. ccha) then
                  lnh4 = .false.
                  lna  = .false.
                  lmg  = .false.
                end if
                if (x .le. ccok) then
                  npok = npok + 1
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'ASP') then
              if (anam2 .eq. ' OD1' .or.
     +            anam2 .eq. ' OD2') then
                nminu = nminu + 1
                nhbac = nhbac + 1
                if (x .le. ccha) then
                  lf = .false.
                end if
                if (x .le. ccok) then
                  nmok = nmok + 1
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'GLU') then
              if (anam2 .eq. ' OE1' .or.
     +            anam2 .eq. ' OE2') then
                nminu = nminu + 1
                nhbac = nhbac + 1
                if (x .le. ccha) then
                  lf = .false.
                end if
                if (x .le. ccok) then
                  nmok = nmok + 1
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'CYS') then
              if (anam2 .eq. ' SG ') then
                ncyss = ncyss + 1
                ndoac = ndoac + 1
                if (x .le. ccha) then
                  lf = .false.
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'MET') then
              if (anam2 .eq. ' SD ') then
                nmets = nmets + 1
                nhbac = nhbac + 1
                if (x .le. ccha) then
                  lf = .false.
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'ASN') then
              if (anam2 .eq. ' OD1') then
                nhbac = nhbac + 1
                goto 200
              else if (anam2 .eq. ' ND2') then
                nhbdo = nhbdo + 1
                if (x .le. ccha) then
                  lnh4 = .false.
                  lna  = .false.
                  lmg  = .false.
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'GLN') then
              if (anam2 .eq. ' OE1') then
                nhbac = nhbac + 1
                goto 200
              else if (anam2 .eq. ' NE2') then
                nhbdo = nhbdo + 1
                if (x .le. ccha) then
                  lnh4 = .false.
                  lna  = .false.
                  lmg  = .false.
                end if
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'THR') then
              if (anam2 .eq. ' OG1') then
                ndoac = ndoac + 1
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'SER') then
              if (anam2 .eq. ' OG ') then
                ndoac = ndoac + 1
                goto 200
              end if
            end if
c
            if (rtyp2 .eq. 'TYR') then
              if (anam2 .eq. ' OH ') then
                ndoac = ndoac + 1
                goto 200
              end if
            end if
c
            if (lwater(rtyp2)) then
              ndoac = ndoac + 1
              goto 200
            end if
c
            if (anam2 .eq. ' O  ') then
              nhbac = nhbac + 1
              goto 200
            end if
c
            if (anam2 .eq. ' N  ') then
              nhbdo = nhbdo + 1
              if (x .le. ccha) then
                lnh4 = .false.
                lna  = .false.
                lmg  = .false.
              end if
              goto 200
            end if
c
            if (anam2(1:2) .eq. ' O') then
              ndoac = ndoac + 1
              goto 200
            end if
c
            if (anam2(1:2) .eq. ' N') then
              ndoac = ndoac + 1
              goto 200
            end if
c
            if (anam2(1:2) .eq. ' S') then
              ndoac = ndoac + 1
              goto 200
            end if
c
          end if
c
  200     continue
        end do
c
        write (lunit,6020) nplus,nminu,nhbdo,nhbac,ndoac,ncyss,nmets
        if (lf) then
          if (nplus .lt. 1) lf = .false.
          if (npok .lt. 1) lf = .false.
        end if
        if (lnh4) then
          if (nminu .lt. 1) lnh4 = .false.
          if (nhbdo .gt. 1) lnh4 = .false.
          if ((nhbac+ndoac) .lt. 3) lnh4 = .false.
          if (nmok .lt. 1) lnh4 = .false.
        end if
        if (lna) then
          if (nminu .lt. 1) lna = .false.
          if (nhbdo .gt. 2) lna = .false.
          if (nmok .lt. 1) lna = .false.
        end if
        if (lmg) then
          if (nminu .lt. 2) lmg = .false.
          if (nhbdo .gt. 2) lmg = .false.
          if (nmok .lt. 1) lmg = .false.
        end if
        if (lh2o) then
          if ((nhbac+ndoac) .lt. 1) lh2o = .false.
          if ((nhbdo+ndoac) .lt. 1) lh2o = .false.
        end if
c
        if (lnh4) then
          write (lunit,6200) 'NH4+ (or other small cation with H)'
          nnh4 = nnh4 + 1
        end if
c
        if (lh2o) then
          write (lunit,6200) 'HOH'
          nh2o = nh2o + 1
        end if
c
        if (lf) then
          write (lunit,6200) 'F- (or other small anion)'
          nf = nf + 1
        end if
c
        if (lna) then
          write (lunit,6200) 'Na+ (or other small monovalent cation)'
          nna = nna + 1
        end if
c
        if (lmg) then
          write (lunit,6200) 'Mg++ (or other small divalent cation)'
          nmg = nmg + 1
        end if
c
        if (.not. (lnh4 .or. lh2o .or. lf .or. lna .or. lmg)) then
          write (lunit,'(a)')
     +      '    Does *not* look like a tightly bound water !'
          nnot = nnot + 1
c
        else if (munit .gt. 0) then
c
          write (munit,'(a)',err=100)
     +      atline(i)(1:length(atline(i)))
c
        end if
c
cATOM   1110  O1  HOH A 300       1.742  13.554   0.117  1.00 21.92
c1234567890123456789012345678901234567890123456789012345678901234567890
c
  100   continue
      end do
c
      if (munit .gt. 0) then
        write (munit,'(a6)',err=101) 'END   '
  101   continue
      end if
c
      write (*,*)
      call jvalut (' Nr of waters scrutinized :',1,nwat)
      call jvalut (' # ? NH4+ (or other small cation with H)     :',
     +  1,nnh4)
      call jvalut (' # ? F-   (or other small anion)             :',
     +  1,nf)
      call jvalut (' # ? Na+  (or other small monovalent cation) :',
     +  1,nna)
      call jvalut (' # ? Mg++ (or other small divalent cation)   :',
     +  1,nmg)
      call jvalut (' # ? HOH  (with reasonable neighbours)       :',
     +  1,nh2o)
      call jvalut (' # non-tightly bound water molecules         :',
     +  1,nnot)
c
 6000 format (/' => Water # ',i5,' = ',a3,'-',a6,' atom ',a4)
 6010 format ( '    Nbr ',a3,'-',a6,'-',a4,' @ ',f6.2,' A')
 6020 format ( '    Nr + charges ',i3,' | Nr - charges ',i3,/
     +         '    Nr Hb donors ',i3,' | Nr Hb accept ',i3,
     +         ' | Nr don/acc   ',i3/
     +         '    Nr CYS S     ',i3,' | Nr MET S     ',i3)
 6200 format ( '    This could be ',a)
c
c ... end
c
 9999 continue
c
      return
      end
c
c
c
      subroutine pdb2o (task,ierr)
c
      include 'xpand.incl'
c
      integer ierr,i,k,nline
      integer length
c
      logical lmtrix
c
      character task*1,line*100,ofmt*20,rtname*20
c
      data ofmt / '(3f15.7)' /
c
code ...
c
      ierr = -1
c
      lmtrix = (task .eq. 'O')
c
      nncs = 0
      nline = 0
c
  100 continue
      read (iunit,'(a)',err=996,end=900) line
      nline = nline + 1
c
cREMARK 350   BIOMT1  54  0.309013  0.755757  0.577344        0.00000            
cREMARK 350   BIOMT2  54 -0.755763  0.563657 -0.333326        0.00000            
cREMARK 350   BIOMT3  54 -0.577363 -0.333334  0.745361        0.00000            
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      if (line(1:19) .eq. 'REMARK 350   BIOMT1') then
        if (.not. lmtrix) then
          read (line(20:23),'(i4)') i
          nncs = nncs + 1
          call jvalut (' BioMT nr :',1,i)
ccc          call jvalut (' O RT nr  :',1,nncs)
          read (line(24:),'(3f10.7,f15.5)')
     +      (rtncs(k,nncs),k=1,10,3)
  200     continue
          read (iunit,'(a)',err=998,end=998) line
          nline = nline + 1
          if (line(1:19) .ne. 'REMARK 350   BIOMT2') goto 200
          read (line(24:),'(3f10.7,f15.5)')
     +      (rtncs(k,nncs),k=2,11,3)
  300     continue
          read (iunit,'(a)',err=998,end=998) line
          nline = nline + 1
          if (line(1:19) .ne. 'REMARK 350   BIOMT3') goto 300
          read (line(24:),'(3f10.7,f15.5)')
     +      (rtncs(k,nncs),k=3,12,3)
        end if
      end if
c
cMTRIX1   1  1.000000  0.000000  0.000000        0.00000    1                    
cMTRIX2   1  0.000000  1.000000  0.000000        0.00000    1                    
cMTRIX3   1  0.000000  0.000000  1.000000        0.00000    1                    
cMTRIX1   2 -0.000007 -0.934176  0.356813        0.00000                         
cMTRIX2   2 -0.356827  0.333330  0.872672        0.00000                         
cMTRIX3   2 -0.934180 -0.127331 -0.333324        0.00000                         
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      if (line(1:6) .eq. 'MTRIX1') then
        if (lmtrix) then
          read (line(7:10),'(i4)') i
          nncs = nncs + 1
          call jvalut (' MTRIX nr :',1,i)
ccc          call jvalut (' O RT nr  :',1,nncs)
          read (line(11:),'(3f10.6,f15.5)')
     +      (rtncs(k,nncs),k=1,10,3)
  400     continue
          read (iunit,'(a)',err=997,end=997) line
          nline = nline + 1
          if (line(1:6) .ne. 'MTRIX2') goto 400
          read (line(11:),'(3f10.6,f15.5)')
     +      (rtncs(k,nncs),k=2,11,3)
  500     continue
          read (iunit,'(a)',err=997,end=997) line
          nline = nline + 1
          if (line(1:6) .ne. 'MTRIX3') goto 500
          read (line(11:),'(3f10.6,f15.5)')
     +      (rtncs(k,nncs),k=3,12,3)
        end if
      end if
c
      goto 100
c
c ... now write O-style operator file
c
  900 continue
      call jvalut (' Nr of lines read :',1,nline)
      call jvalut (' Operators read   :',1,nncs)
      if (nncs .lt. 1) goto 999
c
      call stamp (line)
      write (junit,'(a,a)') '! ',line(1:length(line))
c
      do i=1,nncs
        write (rtname,'(a,i3)') '.LSQ_RT_',i
        call remspa (rtname)
        call textut (' Writing :',rtname)
        write (junit,'(a)') '!'
        write (line,'(a,1x,a,1x,i4,1x,a)') rtname,'R',12,ofmt
        call pretty (line)
        write (junit,'(a)') line(1:length(line))
        write (junit,fmt=ofmt) (rtncs(k,i),k=1,12)
      end do
c
      ierr = 0
c
      goto 999
c
  996 continue
      call errcon ('While reading PDB file')
      goto 999
c
  997 continue
      call errcon ('While reading MTRIX line')
      nncs = nncs - 1
      goto 999
c
  998 continue
      call errcon ('While reading REMARK 350 line')
      nncs = nncs - 1
      goto 999
c
  999 continue
c
      return
      end
c
c
c
      subroutine gmtrix (ierr)
c
      include 'xpand.incl'
c
      integer ierr,i,k
c
      character line*24
c
code ...
c
      ierr = -1
c
 6000 format ('REMARK ',10(a,1x))
 6010 format ('MTRIX',i1,1x,i3,3f10.6,5x,f10.5,5x)
c
      line = ' '
      call gkdate (line(1:24))
c
      write (junit,6000,err=9910)
      write (junit,6000,err=9910) 'NCS MTRX CARDS'
      write (junit,6000,err=9910) 'GENERATED BY :',prognm
      write (junit,6000,err=9910) 'VERSION      :',vers
      write (junit,6000,err=9910) 'DATE         :',line(1:24)
      write (junit,6000,err=9910)
c
      do i=1,nncs
        write (junit,6010,err=9910) 1,i,(rtncs(k,i),k=1,10,3)
        write (junit,6010,err=9910) 2,i,(rtncs(k,i),k=2,11,3)
        write (junit,6010,err=9910) 3,i,(rtncs(k,i),k=3,12,3)
      end do
c
      ierr = 0
c
      goto 9920
c
 9910 continue
      call errcon ('While writing MTRIX file')
      ierr = -1
c
 9920 continue
      close (junit)
c
      return
      end
c
c
c
      subroutine expapt (ierr)
c
      include 'xpand.incl'
c
      real cogxyz(3),cogfra(3),x1(3),x2(3),x3(3)
      real xeff(3),y1(3),y2(3),dis,dmin
c
      integer nut,length,i,j,k,iat,ierr,k1,k2,k3,m1,m2,m3
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .false.
      lrem = .true.
      lxyz = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
      call getcog (cogxyz)
      call fvalut (' Centre  :',3,cogxyz)
      call mulmtx (ca2fra,cogxyz,cogfra,3,3,1)
      call fvalut (' Fractl  :',3,cogfra)
c
      if (cogfra(1) .lt. -1.0 .or. cogfra(1) .gt. 1.0 .or.
     +    cogfra(2) .lt. -1.0 .or. cogfra(2) .gt. 1.0 .or.
     +    cogfra(3) .lt. -1.0 .or. cogfra(3) .gt. 1.0) then
        call errcon ('C-of-G fractional not in [-1,+1] range')
        call prompt (
     +    ' First translate your input model an integer number')
        call prompt (
     +    ' of unit cells until its centre of gravity in fractional')
        call prompt (
     +    ' coordinates lies in the [-1,+1] range (X/Y/Z); then')
        call prompt (
     +    ' re-run this program.')
        ierr = -1
        return
      end if
c
 6000 format ('REMARK ',a,i2)
 6010 format ('REMARK ',a,f10.1)
c
      iat = 0
c
c ... loop over the symmetry operators
c
      do i=1,nsym
c
        write (*,*)
        call ivalut (' Checking operator :',1,i)
c
c ... apply operator to FRACTIONAL centre-of-gravity
c
        call vecrtv (cogfra,x2,1,rtsym(1,i),rtsym(10,i))
c
c ... loop over 7*7*7 unit cells
c
        dmin = 10.0**30
        m1 = 999
        m2 = 999
        m3 = 999
c
        do k1=-3,3
          do k2=-3,3
            do k3=-3,3
c
c ... translate molecule
c
              x1(1) = x2(1) + float(k1)
              x1(2) = x2(2) + float(k2)
              x1(3) = x2(3) + float(k3)
c
c ... get CoG in Cartesian coordinates
c
              call mulmtx (fra2ca,x1,x3,3,3,1)
c
              dis = (cpoint(1)-x3(1))**2 +
     +              (cpoint(2)-x3(2))**2 +
     +              (cpoint(3)-x3(3))**2
c
              if (dis .lt. dmin) then
ccc                print *,' K123, DIS ',k1,k2,k3,sqrt(dis)
                dmin = dis
                m1 = k1
                m2 = k2
                m3 = k3
              end if
            end do
          end do
        end do
c
        if (m1 .gt. 900) then
          call errcon ('No suitable translation found !?')
          ierr = -1
          return
        end if
c
        dmin = sqrt(dmin)
        call fvalut (' Minimal distance (A) :',1,dmin)
        write (*,'(3(a,i2))') ' t1=',m1,' t2=',m2,' t3=',m3
c
        k1 = m1
        k2 = m2
        k3 = m3
c
c ... generate the atoms explicitly and check distance
c
        xeff(1) = float(k1)
        xeff(2) = float(k2)
        xeff(3) = float(k3)
c
        do j=1,natoms
c
c ... fractionalise, apply operator & translation, orthogonalise
c
          call mulmtx (ca2fra,atmxyz(1,j),y2,3,3,1)
          call vecrtv (y2,y1,1,rtsym(1,i),rtsym(10,i))
          do k=1,3
            y1(k)=y1(k)+xeff(k)
          end do
          call mulmtx (fra2ca,y1,y2,3,3,1)
c
          iat = iat + 1
c
          if (j .eq. 1) then
            write (junit,6000,err=9910)
            write (junit,6000,err=9910)
     +         'SGS Operator nr ',i
            write (junit,6000,err=9910)
     +        'Translation A ',nint(xeff(1))
            write (junit,6000,err=9910)
     +        'Translation B ',nint(xeff(2))
            write (junit,6000,err=9910)
     +        'Translation C ',nint(xeff(3))
            write (junit,6010,err=9910)
     +        'Distance (A)  ',dmin
            write (junit,6000,err=9910)
            nut = nut + 7
          end if
c
          line = atline(j)
          write (line(31:54),'(3f8.3)',err=9920)
     +      (y2(k),k=1,3)
c
          write (junit,'(a)',err=9910)
     +      line(1:length(line))
          nut = nut + 1
c
        end do
c
      end do
c
      write (junit,'(a)',err=9910) 'END   '
      nut = nut + 1
c
      write (*,*)
      call jvalut (' Nr of atoms written :',1,iat)
      call jvalut (' Nr of lines written :',1,nut)
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      goto 9999
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine expatm (ierr)
c
      include 'xpand.incl'
c
      real atmfra(3),y2(3),x2(3),x3(3)
c
      integer nut,length,i,j,k,iat,ierr,noff,ii1,ii2,ii3
c
      logical lecho,lxyz,lrem
c
      character*128 line
c
code ...
c
      nut = 0
      ierr = 0
      lecho = .false.
      lxyz = .true.
      lrem = .true.
c
      call getpdb (lrem,lecho,lxyz,nut,ierr)
      if (ierr .ne. 0) goto 9999
c
 6000 format ('REMARK ',a,i2)
c
      iat = 0
c
      do i=1,natoms
        symbox (i) = 0
        ncsbox (i) = 0
      end do
      noff = 0
c
c ... loop over atoms
c
      do j=1,natoms
c
c ... get fractional coordinates
c
        call mulmtx (ca2fra,atmxyz(1,j),atmfra,3,3,1)
c
c ... loop over the symmetry operators
c     (INCLUDING THE IDENTITY OPERATOR (NR 1) !!!)
c
        do i=1,nsym
c
c ... apply operator to FRACTIONAL coordinates
c
          call vecrtv (atmfra,x2,1,rtsym(1,i),rtsym(10,i))
c
c ... force inside the unit cell at (0,0,0)/(1,1,1)
c
          do k=1,3
  220       if (x2(k) .lt. 0.0) then
              x2(k) = x2(k) + 1.0
              goto 220
            end if
  230       if (x2(k) .ge. 1.0) then
              x2(k) = x2(k) - 1.0
              goto 230
            end if
            if (x2(k) .lt. 0.0) then
              call errcon ('Unexpected force-symm problem')
              goto 1999 
            end if
          end do
c
c ... now try all 9 unit cells around the one from 0->1
c
          do ii1=-1,1
            x3(1) = x2(1) + float(ii1)
            if (x3(1).ge.fralo(1) .and.
     +          x3(1).lt.frahi(1)) then
              do ii2=-1,1
                x3(2) = x2(2) + float(ii2)
                if (x3(2).ge.fralo(2) .and.
     +              x3(2).lt.frahi(2)) then
                  do ii3=-1,1
                    x3(3) = x2(3) + float(ii3)
                    if (x3(3).ge.fralo(3) .and.
     +                  x3(3).lt.frahi(3)) then
c
                      call mulmtx (fra2ca,x3,y2,3,3,1)
c
c ... write ATOM/HETATM card
c
                      iat = iat + 1
                      line = atline(j)
                      write (line(31:54),'(3f8.3)',err=9920)
     +                  (y2(k),k=1,3)
                      write (junit,'(a)',err=9910)
     +                  line(1:length(line))
                      nut = nut + 1
                      ncsbox (j) = ncsbox (j) + 1
c
                    end if
                  end do
                end if
              end do
            end if
          end do
c
 1999     continue
c
        end do
c
        if (ncsbox(j) .le. 0) then
          noff = noff + 1
        else
          symbox (ncsbox(j)) = symbox (ncsbox(j)) + 1
        end if
c
      end do
c
      write (junit,'(a)',err=9910) 'END   '
      nut = nut + 1
c
      write (*,*)
      call jvalut (' Nr of atoms written :',1,iat)
      call jvalut (' Nr of lines written :',1,nut)
c
c ... histogram
c
      write (*,*)
      call prompt (' Off-spring counts :')
c
 6392 format (' Nr of atoms with ',i6,' off-spring = ',i6)
c
      do i=natoms,1,-1
        if (symbox(i) .gt. 0) then
          write (*,6392) i,symbox(i)
        end if
      end do
      write (*,6392) 0,noff
c
      goto 9999
c
c ... errors
c
 9910 continue
      write (*,*)
      call errcon ('While writing PDB file')
      write (*,*)
      ierr = -1
      goto 9999
c
 9920 continue
      write (*,*)
      call errcon ('During internal I/O operation')
      write (*,*)
      ierr = -1
      goto 9999
c
c ... end
c
 9999 continue
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine symcel (spgrp,ierr)
c
c ... SYMCEL
c
      include 'xpand.incl'
c
      integer ierr,i
c
      character spgrp*(*)
      character*128 line
c
code ...
c
      ierr = -1
      spgrp = ' '
c
      write (*,*)
      call prompt (' Reading PDB file ...')
c
c ... CRYST1   45.650   47.560   77.610  90.00  90.00  90.00 P 21 21 21    4  1CBS 216
c     1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c              1         2         3         4         5         6         7         8         9        10
c
   10 continue
      read (iunit,'(a)',err=9900,end=1000) line
      call upcase (line)
      if (line(1:6) .eq. 'CRYST1') then
        call textut (' CRYST1 card found :',line)
        read (line,'(6x,3f9.3,3f7.2)',err=9800) (cell(i),i=1,6)
        spgrp = line(56:66)
        ierr = 0
        return
      end if
c
      goto 10
c
 9800 continue
      call errcon ('While reading CRYST1 card')
      return
c
 9900 continue
      call errcon ('While reading PDB file')
      return
c
 1000 continue
      call errcon ('No CRYST1 card found')
      return
c
      end
