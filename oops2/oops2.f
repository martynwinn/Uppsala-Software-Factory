      program oops2
c
c ... find "bad" spots in a protein model
c
c ... Gerard Kleywegt @ 930304 (OOPS) / 990311 (OOPS2)
c
c ... changes OOPS2:
c     - interface
c     - default file names
c     - subjectivity stuff
c     - "learn" from previous subjective score
c
      include 'oops2.incl'
c
      integer maxopt,maxhis,maxcom
      parameter (maxopt = 10, maxhis = 25, maxcom = 50)
c
      real xdum(maxres),histo(maxhis),stats(6,8)
      real ave,sdv,xmin,xmax,xtot,perc,rdum
c
      integer length,iunit,i,j,ierr,nopt,ires,nbad,iuser
      integer nhisto,nhcnt(maxhis+1),nwat,junit,nmacro,npos
      integer npep,nrsf,nrsc,nram,nmas,nblo,nbhi,nrsm,nrss
      integer npla,nzet,nrus(maxuse),nbco,nokay,nqlo,nqhi
      integer nrsr,ndbb,nins,nmut,ndis,ntem,nocc,nphi,nchi
      integer nstats,iflip,irsfal,irsfmc,irsfsc,irsral,npro
      integer irsc,iplanp,ichirl,leng1,kunit,nome,ncis,ndaa
      integer nssd,nsst,idisu,k,nl,l,nbwi,iindex,idum,numcom
      integer nxgood,nipoor,nlef,moltyp,itype,cnttyp(4)
      integer cntbad(4)
c
      logical xinter,linter,baddy,badpep,badrsf,lfirst,lwater
      logical badram,listem,badmas,badblo,badbhi,badrsc
      logical badrsm,badrss,badpla,badzet,lstats,badcon
      logical baduse(maxuse),lplot,lalres,lchain,lmanam
      logical badqlo,badqhi,badwat,badrsr,lnotes,baddbb
      logical baddis,badtem,badocc,badphi,badchi,lprwat
      logical lhtml,badome,badcis,baddaa,badpro,badpos,badisu
      logical badssd,badsst,lpos,lpro,lworst,lrelax,badwif
      logical lsubj,lknow,lxgood,lipoor,badpor,lignore,llefth
      logical badlef,lgrab,laltyp(4)
c
      character line*200,optpar(maxopt)*80,reply*1,filnam*80
      character gktext*256,title*80,notfil*80,html*80,myline*80
      character filnext*80,xlines(80)*80,dbsub*25,comand*3
      character comlin(maxcom)*39,rnamdb*25,rtypdb*25
      character chjunk*6,chdoll*6,chtype*6,typnam(4)*20
c
      data npep,nrsf,nrsc,nram,nmas,nblo,nbhi,nrsm,nrss /9*0/
      data npla,nzet,nbco,nqlo,nqhi,nwat,nrsr,ndbb,nins /9*0/
      data nmut,ndis,ntem,nocc,nphi,nchi,nome,ncis,ndaa /9*0/
      data npro,npos,nssd,nsst,nbwi,nlef /6*0/
c
      data nrus /maxuse*0/, cntbad /4*0/
      data mcname /' C  ',' CA ',' N  ',' O  ',' CB '/
      data iflip,irsfal,irsfmc,irsfsc,irsral /5*-1/
      data irsc,iplanp,ichirl /3*-1/, cnttyp /4*0/
      data typnam /'Protein','Nucleic acid','Water','Heterogen'/
c
      data lpep,lrsfit,lrama,lmask,llob,lhib,lrsc /7*.false./
      data lrsm,lrss,lplan,lchir,lbadco,lloq,lhiq /7*.false./
      data lwat,lrsr,ldbb,lprev,ldisu,lwif /6*.false./
      data lpos,lpro,lworst,lrelax,llefth /5*.false./
      data lsubj /.true./, laltyp /4*.false./
c
code ...
c
c --- INITIALISATION
c
      call gkinit (prognm,vers)
c
      write (*,6000) maxres,maxbuf,maxatm,maxss,maxwif,maxuse
 6000 format (/' Max nr of residues              : ',i10/
     +         ' Max size of mask                : ',i10/
     +         ' Max nr of atoms                 : ',i10/
     +         ' Max nr of disulfide bridges     : ',i10/
     +         ' Max nr of WHAT IF messages      : ',i10/
     +         ' Max nr of user-defined criteria : ',i10/)
c
      linter = xinter()
c
      iunit  = 10
      junit  = 11
      kunit  = 12
c
c ... first check if directory "oops" exists
c
      line = 'oops/test'
      call xopxua (iunit,line,.false.,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot open file oops/test')
        call prompt (' Probably directory oops does not exist')
        call prompt (' Create it with: mkdir oops')
        call prompt (' And re-run OOPS')
        call errstp ('No (access to) oops directory')
      end if
      close (unit=iunit,status='DELETE')
c
ccc      call system ('rm oops/test')
c
c ... default cut-off values
c
      pepcut = 2.5
      rsfcut = 0.7
      rsfwat = 0.7
      rsmcut = 0.8
      rsscut = 0.6
      mskrad = 2.0
      lobfac = 5.0
      hibfac = 40.0
      hibwat = 40.0
      loqocc = 0.99
      hiqocc = 1.0
      rsccut = 1.0
      placut = 5.8
      omecut = 6.0
      chicut = 3.5
      rsrcut = 0.3
      rsrwat = 0.3
      dbbcut = 5.0
      mxrmsd = 0.5
      mxrmsb = 5.0
      mxrmsq = 0.1
      mxphps = 20.0
      mxchid = 20.0
      bmpcut = 0.2
c
c ... watcut default = 100*exp(-2), i.e. B = 8 * resoln^2
c     i.e., if resoln = 2.5 A -> Bmax = 50 (OCC = 1.0)
c                       2.0 A -> Bmax = 32
c                       1.5 A -> Bmax = 18
c                       1.2 A -> Bmax = 11.5
c                       1.0 A -> Bmax =  8
c
      watcut = 13.53
      nbadcu = 1
c
      resoln = 2.0
      mxcaca = 4.5
      mxatat = 2.0
c
c ... YASSPA cut-offs
c
      yalcut = 0.5
      ybecut = 0.8
      ylhcut = 0.3
c
      natoms = 0
      nres   = 0
      oline = 'bell   save   print DONE and SAVED'
      nwif = 0
c
      do i=1,maxres
        pepflip (i) = 0.0
        rsfit (i)   = 1.0
        rsfitm(i)   = 1.0
        rsfits(i)   = 1.0
        rsc (i)     = 0.0
        phi (i)     = 999.99
        psi (i)     = 999.99
        plangl (i)  = 999.99
        zeta (i)    = 999.99
        watqua (i)  = -1.0
        write (resnam(i),'(i6)') i
        call remspa (resnam(i))
        restyp (i) = 'UNK'
        count (i) = 0
        do j=1,maxuse
          usercr (i,j) = 0.0
        end do
        swater (i) = .false.
        rsrfac (i) = 0.0
        rmsdbb (i) = 0.0
        xxrmsd (i) = 0.0
        xxrmsb (i) = 0.0
        xxrmsq (i) = 0.0
        xxphps (i) = 0.0
        xxchid (i) = 0.0
        insert (i) = .false.
        mutate (i) = .false.
      end do
c
      nstats = 0
c
      do i=1,maxuse
        write (usernm (i),*) 'User-criterion nr ',i
        call pretty (usernm(i))
c
        write (userln (i),*) 'Bad user-criterion nr ',i
        call pretty (userln(i))
c
        usersn (i) = '>'
        usercu (i) = 999.999
c
      end do
c
      do i=1,maxcom
        comlin (i) = ' '
      end do
c
c ... PRINT STATISTICS AND HISTOGRAMS ?
c
      reply = 'Y'
      call textin (' Print statistics and histograms ?',reply)
      call upcase (reply)
      lstats = (reply .eq. 'Y')
c
      if (lstats) then
        reply = 'N'
        call textin (' Auto-generate (some) O2D plot files ?',reply)
        call upcase (reply)
        lplot = (reply .eq. 'Y')
      end if
c
c ... MOLECULE NAME
c
      molnam ='M1'
      write (*,*)
      call textin (' Molecule name in O ?',molnam)
      if (length(molnam) .lt. 1) then
        call errstp ('No molecule name given')
      end if
      call upcase (molnam)
c
c ... PDB file
c
      write (*,*)
      filnam = molnam(1:leng1(molnam))//'.pdb'
      call locase (filnam)
      call textin (' PDB file  ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening PDB file')
      end if
c
      call getpdb (iunit,ierr)
      if (ierr .ne. 0) then
        call errstp ('While reading PDB file')
      end if
c
c ... figure out which residues are waters
c
      nwater = 0
      do i=1,nres
        if (lwater(restyp(i))) then
          nwater = nwater + 1
          swater (i) = .true.
        end if
      end do
      write (*,*)
      call jvalut (' Nr of WATERs :',1,nwater)
c
c ... read old subjectivity table if it exists
c
      write (*,*)
      do i=1,nres
        oldsub (i) = -1
      end do
      lknow = .false.
      filnam = 'oops_subjectivity.table'
      call xexist (iunit,filnam,.false.,ierr)
      if (ierr .ne. 0) then
        call prompt (' Old subjective score table not found')
        call prompt (' All residues set to status: unexamined')
      else
        call prompt (' Reading old subjective score table')
        call ressub (iunit)
        close (iunit)
        lknow = .true.
      end if
c
c ... list available commands
c
 1900 continue
      write (*,*)
      call prompt (' POSSIBLE COMMANDS :')
c
      if (lpep) then
        comlin (1) = ' PEP = switch off pep-flip'
      else
        comlin (1) = ' PEP = include pep-flip'
      end if
c
      if (lrsc) then
        comlin (2) = ' RSC = switch off rotamer'
      else
        comlin (2) = ' RSC = include rotamer'
      end if
c
      if (lrama) then
        comlin (3) = ' RAM = switch off Ramachandran'
      else
        comlin (3) = ' RAM = include Ramachandran'
      end if
c
      if (lwif) then
        comlin (4) = ' WIF = switch off WHAT IF diagnostics'
      else
        comlin (4) = ' WIF = include WHAT IF diagnostics'
      end if
c
      if (lrsr) then
        comlin (5) = ' RSR = switch off real-space R (all)'
      else
        comlin (5) = ' RSR = include real-space R (all atoms)'
      end if
c
      if (lrsfit) then
        comlin (6) = ' CCA = switch off real-space CC (all)'
      else
        comlin (6) = ' CCA = include real-space CC (all)'
      end if
c
      if (lrsm) then
        comlin (7) = ' CCM = switch off real-space CC (main)'
      else
        comlin (7) = ' CCM = include real-space CC (main)'
      end if
c
      if (lrss) then
        comlin (8) = ' CCS = switch off real-space CC (side)'
      else
        comlin (8) = ' CCS = include real-space CC (side)'
      end if
c
      if (llob) then
        comlin (9) = ' BFA = switch off B-factors'
      else
        comlin (9) = ' BFA = include B-factors'
      end if
c
      if (lloq) then
        comlin (10) = ' OCC = switch off occupancies'
      else
        comlin (10) = ' OCC = include occupancies'
      end if
c
      if (lmask) then
        comlin (11) = ' MSK = switch off mask fit'
      else
        comlin (11) = ' MSK = include mask fit'
      end if
c
      if (ldisu) then
        comlin (12) = ' DIS = switch off disulfides'
      else
        comlin (12) = ' DIS = include disulfides'
      end if
c
      if (lchir) then
        comlin (13) = ' CAC = switch off CA chirality'
      else
        comlin (13) = ' CAC = include CA chirality'
      end if
c
      if (lplan) then
        comlin (14) = ' PLA = switch off peptide planarity'
      else
        comlin (14) = ' PLA = include peptide planarity'
      end if
c
      if (lprev) then
        comlin (15) = ' PRE = switch off previous model'
      else
        comlin (15) = ' PRE = include previous model'
      end if
c
      numcom = 15
c
c ... various program settings
c
      if ( 2*(numcom/2) .ne. numcom) then
        numcom = numcom + 1
        comlin (numcom) = ' '
      end if
      numcom = numcom + 1
      comlin (numcom) = ' ---------------------------------------'
      numcom = numcom + 1
      comlin (numcom) = ' ---------------------------------------'
      numcom = numcom + 1
      comlin (numcom) = ' GO  = get going !'
      numcom = numcom + 1
      comlin (numcom) = ' QUI = quit without doing anything'
c
      do i=1,numcom,2
        write (*,'(a,1x,a)') comlin(i),comlin(i+1)
      end do
c
c ... enter command loop here
c
 2000 continue
      comand = '?'
      write (*,*)
      call textin (' Next command ?',comand)
      call remspa (comand)
      call upcase (comand)
c
      if (comand .eq. '?  ') goto 1900
c
      if (comand .eq. 'PEP') then
        lpep = (.not. lpep)
        if (lpep) goto 6500
      else if (comand .eq. 'RSC') then
        lrsc = (.not. lrsc)
        if (lrsc) goto 6550
      else if (comand .eq. 'RAM') then
        lrama = (.not. lrama)
        if (lrama) goto 6620
      else if (comand .eq. 'WIF') then
        lwif = (.not. lwif)
        if (lwif) goto 6680
      else if (comand .eq. 'RSR') then
        lrsr = (.not. lrsr)
        if (lrsr) goto 6540
      else if (comand .eq. 'CCA') then
        lrsfit = (.not. lrsfit)
        if (lrsfit) goto 6510
      else if (comand .eq. 'CCM') then
        lrsm = (.not. lrsm)
        if (lrsm) goto 6520
      else if (comand .eq. 'CCS') then
        lrss = (.not. lrss)
        if (lrss) goto 6530
      else if (comand .eq. 'BFA') then
        llob = (.not. llob)
        if (llob) goto 6570
      else if (comand .eq. 'OCC') then
        lloq = (.not. lloq)
        if (lloq) goto 6600
      else if (comand .eq. 'MSK') then
        lmask = (.not. lmask)
        if (lmask) goto 6560
      else if (comand .eq. 'DIS') then
        ldisu = (.not. ldisu)
        if (ldisu) goto 6645
      else if (comand .eq. 'PRE') then
        lprev = (.not. lprev)
        if (lprev) goto 6650
      else if (comand .eq. 'CAC') then
        lchir = (.not. lchir)
        if (lchir) goto 6640
      else if (comand .eq. 'PLA') then
        lplan = (.not. lplan)
        if (lplan) goto 6630
      else if (comand .eq. 'GO ') then
        goto 3500
      else if (comand .eq. 'QUI') then
        reply = 'N'
        call textin (' QUIT - Are you sure (Y|N) ?',reply)
        call upcase (reply)
        if (reply .ne. 'Y') goto 2000
        goto 9900
      else
        call textut (' ERROR - Invalid command :',comand)
        goto 1900
      end if
c
      goto 2000






c
c ... PEP-FLIP VALUES
c
 6500 continue
c
      filnam = 'pepflip.o'
      call textin (
     +  ' O data block with pep-flip values ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening datablock file')
        goto 6501
      end if
      call skipem (iunit)
c
      read (iunit,'(a)',err=6501,end=6501) line
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (nopt .lt. 4 .or. ierr .ne. 0) goto 6501
c
      call prdb (optpar(1),optpar(2),optpar(3),optpar(4))
      call str2i (optpar(3),ires,ierr)
      if (ierr .ne. 0) goto 6501
c
      call odbchk (molnam,'_residue_pepflip',optpar(1),
     +             'r',optpar(2),nres,ires,lignore)
      if (lignore) goto 6501
c
      if (ires .gt. maxres) then
        call errcon ('Too many residues')
        goto 6501
      end if
c
      read (iunit,optpar(4),err=6501,end=6501)
     +    (pepflip(i),i=1,ires)
      close (iunit)
c
      if (lstats) then
        j=0
        do i=1,nres
          if (pepflip(i) .gt. 0.001) then
            j = j + 1
            xdum (j) = pepflip(i)
          end if
        end do
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 0.25 * float(i-1)
        end do
        call oopsts ('Pep-flip values (>0)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_pepflip.plt'
        title = 'Pep-flip values'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,j,xdum,0.0,1.05*xmax,j,ave,sdv,xmin,xmax)
c
      end if
c
      write (*,*)
      call fvalin (' Pep-flip cut-off ?',1,pepcut)
      if (lstats) call outlie (j,xdum,pepcut,0)
c
      if (lstats) then
c
        nstats = nstats + 1
        iflip = nstats
        stats (1,nstats) = pepcut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
 6501 continue
      call errcon ('Sorry - pep-flip NOT included')
      lpep = .false.
      goto 2000
c
c ... RSC VALUES
c
 6550 continue
c
      filnam = 'rsc.o'
      call textin (
     +  ' O data block with RSC values   ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening datablock file')
        goto 6551
      end if
      call skipem (iunit)
c
      read (iunit,'(a)',err=6551,end=6551) line
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (nopt .lt. 4 .or. ierr .ne. 0) goto 6551
c
      call prdb (optpar(1),optpar(2),optpar(3),optpar(4))
      call str2i (optpar(3),ires,ierr)
      if (ierr .ne. 0) goto 6551
c
      call odbchk (molnam,'_residue_rsc',optpar(1),
     +             'r',optpar(2),nres,ires,lignore)
      if (lignore) goto 6551
c
      if (ires .gt. maxres) then
        call errcon ('Too many residues')
        goto 6551
      end if
c
      read (iunit,optpar(4),err=6551,end=6551)
     +  (rsc(i),i=1,ires)
      close (iunit)
c
      if (lstats) then
        j=0
        do i=1,nres
          if (rsc(i) .gt. 0.001 .and.
     +        restyp(i) .ne. 'GLY' .and.
     +        restyp(i) .ne. 'ALA') then
            j = j + 1
            xdum (j) = rsc(i)
          end if
        end do
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 0.25 * float(i-1)
        end do
        call oopsts ('RSC-fit values (non-GLY/ALA) (>0)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_rsc.plt'
        title = 'RSC-fit values (non-Gly/Ala)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,nres,rsc,0.0,1.05*xmax,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' RSC cut-off ?',1,rsccut)
      if (lstats) call outlie (j,xdum,rsccut,0)
c
      if (lstats) then
c
        nstats = nstats + 1
        irsc = nstats
        stats (1,nstats) = rsccut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
 6551 continue
      call errcon ('Sorry - rotamer NOT included')
      lrsc = .false.
      goto 2000
c
c ... RAMACHANDRAN
c
 6620 continue
c
c ... initialise regions in PHI-PSI space which are okay
c
      write (*,*) 'Checking allowed PHI-PSI areas'
      call phipsi ()
      call newram ()
c
      reply = 'N'
      write (*,*)
      call textin (' Check unusual Pro PHI ?',reply)
      call upcase (reply)
      lpro = (reply .ne. 'N')
c
      reply = 'N'
      write (*,*)
      call textin (' Check positive non-Gly PHI ?',reply)
      call upcase (reply)
      lpos = (reply .ne. 'N')
c
      reply = 'N'
      write (*,*)
      call textin (' Check left-handed helical residues ?',reply)
      call upcase (reply)
      llefth = (reply .ne. 'N')
c
      goto 2000
c
c ... WHAT IF DIAGNOSTICS
c
 6680 continue
c
      filnam = 'pdbout.txt'
      call textin (
     +  ' WHAT IF report file ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening WHAT IF report file')
        goto 6681
      end if
c
c      reply = 'Y'
c      call textin (' Nit-picking mode (e.g., final model) ?',reply)
c      call upcase (reply)
c      lrelax = (reply .eq. 'N')
c
      lrelax = .true.
c
ccc      if (lrelax) then
        call fvalin (' Cut-off for short contacts (bumps) ?',1,bmpcut)
ccc      end if
c
      nwif = 0
      call whatif (iunit,lrelax)
c
      write (*,*)
      call jvalut (' Nr of WHAT IF diagnostics :',1,nwif)
c
      goto 2000
c
 6681 continue
      call errcon ('Sorry - WHAT IF diagnostics NOT included')
      lwif = .false.
      goto 2000
c
c ... RS R-factor (all atoms)
c
 6540 continue
c
      reply = 'O'
      call textin (' Read values from O or MAPMAN (O|M) ?',reply)
      call upcase (reply)
c
      if (reply .eq. 'M') then
c
        filnam = 'rs_fit.list'
        call textin (
     +    ' MAPMAN list file with RS-fit values   ?',filnam)
        call xopxoa (iunit,filnam,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening file')
          goto 6511
        end if
        call skipem (iunit)
c
        do i=1,nres
          rsrfac (i) = -9.99
        end do
c
        j = 0
 6543   continue
        read (iunit,'(a)',err=6541,end=6545) line
        if (line(1:1) .eq. '!') goto 6543
        read (line(27:32),*,err=6541) rdum
        chjunk = line(13:18)
        call remspa (chjunk)
        chdoll = '$' // chjunk
        chtype = line(9:11)
        call remspa (chtype)
        j = j + 1
        if ( (resnam(j) .eq. chjunk .or.
     +        resnam(j) .eq. chdoll) .and.
     +      restyp(j) .eq. chtype) then
          rsrfac (j) = rdum
          goto 6543
        end if
        do k=1,nres
          if ( (resnam(k) .eq. chjunk .or.
     +          resnam(k) .eq. chdoll) .and.
     +        restyp(k) .eq. chtype) then
            rsrfac (k) = rdum
            goto 6543
          end if
        end do
        call errcon ('Residue not found !')
        call textut (' >',line)
        goto 6543
c
 6545   continue
        close (iunit)
c
      else
c
        filnam = 'rsrfac_all.o'
        call textin (
     +    ' O data block with RS R-factors  ?',filnam)
        call xopxoa (iunit,filnam,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening datablock file')
          goto 6541
        end if
        call skipem (iunit)
c
        read (iunit,'(a)',err=6541,end=6541) line
        call extrop (line,nopt,maxopt,optpar,ierr)
        if (nopt .lt. 4 .or. ierr .ne. 0) goto 6541
c
        call prdb (optpar(1),optpar(2),optpar(3),optpar(4))
        call str2i (optpar(3),ires,ierr)
        if (ierr .ne. 0) goto 6541
c
        call odbchk (molnam,'_residue_rsfit',optpar(1),
     +               'r',optpar(2),nres,ires,lignore)
        if (lignore) goto 6541
c
        if (ires .gt. maxres) then
          call errcon ('Too many residues')
          goto 6541
        end if
c
        read (iunit,optpar(4),err=6541,end=6541)
     +    (rsrfac(i),i=1,ires)
        close (iunit)
c
      end if
c
      if (lstats) then
        j=0
        do i=1,nres
          if (rsrfac(i) .ge. 0.001) then
            j = j + 1
            xdum (j) = rsrfac(i)
          end if
        end do
        if (j .le. 0) then
          call errcon ('No residues with RS R-factor found')
          goto 6541
        end if
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 0.025 * float(i-1)
        end do
        call oopsts ('RS R-factors (all atoms) (>0)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        rsrcut = max(xmin+0.1,min(xmax-0.1,ave+2.0*sdv))
        rsrwat = max(xmin+0.1,min(xmax-0.1,ave+3.0*sdv))
c
        filnam = molnam(1:leng1(molnam))//'_rsrfac_all.plt'
        title = 'RS R-factors (all atoms)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,j,xdum,0.0,xmax,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' RS R-factor cut-off ?',1,rsrcut)
      rsrwat = rsrcut
      if (lstats) call outlie (j,xdum,rsrcut,0)
c
ccc      call fvalin (' RS R-factor cut-off WATERs ?',1,rsrwat)
c
      if (lstats) then
c
        nstats = nstats + 1
        irsral = nstats
        stats (1,nstats) = rsrcut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
 6541 continue
      call errcon ('Sorry - real-space R NOT included')
      lrsr = .false.
      goto 2000
c
c ... RS-FIT VALUES (all atoms)
c
 6510 continue
c
      reply = 'O'
      call textin (' Read values from O or MAPMAN (O|M) ?',reply)
      call upcase (reply)
c
      if (reply .eq. 'M') then
c
        filnam = 'rs_fit.list'
        call textin (
     +    ' MAPMAN list file with RS-fit values   ?',filnam)
        call xopxoa (iunit,filnam,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening file')
          goto 6511
        end if
        call skipem (iunit)
c
        do i=1,nres
          rsfit (i) = -9.99
        end do
c
        j = 0
 6513   continue
        read (iunit,'(a)',err=6511,end=6515) line
        if (line(1:1) .eq. '!') goto 6513
        read (line(21:26),*,err=6511) rdum
        chjunk = line(13:18)
        call remspa (chjunk)
        chdoll = '$' // chjunk
        chtype = line(9:11)
        call remspa (chtype)
        j = j + 1
        if ( (resnam(j) .eq. chjunk .or.
     +        resnam(j) .eq. chdoll) .and.
     +      restyp(j) .eq. chtype) then
          rsfit (j) = rdum
          goto 6513
        end if
        do k=1,nres
          if ( (resnam(k) .eq. chjunk .or.
     +          resnam(k) .eq. chdoll) .and.
     +        restyp(k) .eq. chtype) then
            rsfit (k) = rdum
            goto 6513
          end if
        end do
        call errcon ('Residue not found !')
        call textut (' >',line)
        goto 6513
c
 6515   continue
        close (iunit)
c
      else
c
        filnam = 'rsfit_all.o'
        call textin (
     +    ' O data block with RS-fit values   ?',filnam)
        call xopxoa (iunit,filnam,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening datablock file')
          goto 6511
        end if
        call skipem (iunit)
c
        read (iunit,'(a)',err=6511,end=6511) line
        call extrop (line,nopt,maxopt,optpar,ierr)
        if (nopt .lt. 4 .or. ierr .ne. 0) goto 6511
c
        call prdb (optpar(1),optpar(2),optpar(3),optpar(4))
        call str2i (optpar(3),ires,ierr)
        if (ierr .ne. 0) goto 6511
c
        call odbchk (molnam,'_residue_rsfit',optpar(1),
     +               'r',optpar(2),nres,ires,lignore)
        if (lignore) goto 6511
c
        if (ires .gt. maxres) then
          call errcon ('Too many residues')
          goto 6511
        end if
c
        read (iunit,optpar(4),err=6511,end=6511)
     +    (rsfit(i),i=1,ires)
        close (iunit)
c
      end if
c
      if (lstats) then
        j=0
        do i=1,nres
          if (rsfit(i) .ge. -1.0) then
            j = j + 1
            xdum (j) = rsfit(i)
          end if
        end do
        if (j .le. 0) then
          call errcon ('No residues with RS-fit values found')
          goto 6511
        end if
        nhisto = 20
        do i=1,nhisto
          histo(i) = 0.05 * float(i-1)
        end do
        call oopsts ('RS-fit values (all atoms)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        rsfcut = max(xmin+0.1,min(xmax-0.1,ave-2.0*sdv))
        rsfwat = max(xmin+0.1,min(xmax-0.1,ave-3.0*sdv))
c
        filnam = molnam(1:leng1(molnam))//'_rsfit_all.plt'
        title = 'RS-fit values (all atoms)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,j,xdum,min(0.0,xmin),1.0,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' RS-fit cut-off ?',1,rsfcut)
      rsfwat = rsfcut
      if (lstats) call outlie (j,xdum,rsfcut,1)
c
ccc      call fvalin (' RS-fit cut-off WATERs ?',1,rsfwat)
c
      if (lstats) then
c
        nstats = nstats + 1
        irsfal = nstats
        stats (1,nstats) = rsfcut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
 6511 continue
      call errcon ('Sorry - real-space CC (all) NOT included')
      lrsfit = .false.
      goto 2000
c
c ... RS-FIT VALUES (main chain atoms)
c
 6520 continue
c
      filnam = 'rsfit_mc.o'
      call textin (
     +  ' O data block with RS-fit values   ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening datablock file')
        goto 6521
      end if
      call skipem (iunit)
c
      read (iunit,'(a)',err=6521,end=6521) line
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (nopt .lt. 4 .or. ierr .ne. 0) goto 6521
c
      call prdb (optpar(1),optpar(2),optpar(3),optpar(4))
      call str2i (optpar(3),ires,ierr)
      if (ierr .ne. 0) goto 6521
c
      call odbchk (molnam,'_residue_rsfit',optpar(1),
     +             'r',optpar(2),nres,ires,lignore)
      if (lignore) goto 6521
c
      if (ires .gt. maxres) then
        call errcon ('Too many residues')
        goto 6521
      end if
c
      read (iunit,optpar(4),err=6521,end=6521)
     +  (rsfitm(i),i=1,ires)
      close (iunit)
c
      if (lstats) then
        j=0
        do i=1,nres
          if (rsfitm(i) .ge. -1.0 .and.
     +        (.not. swater(i)) .and.
     +        abs(rsfitm(i)) .ge. 0.0001) then
            j = j + 1
            xdum (j) = rsfitm(i)
          end if
        end do
        nhisto = 20
        do i=1,nhisto
          histo(i) = 0.05 * float(i-1)
        end do
        call oopsts ('RS-fit values (main-chain atoms)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        rsmcut = max(xmin+0.1,min(xmax-0.1,ave-2.0*sdv))
c
        filnam = molnam(1:leng1(molnam))//'_rsfit_mc.plt'
        title = 'RS-fit values (main-chain atoms)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,j,xdum,min(0.0,xmin),1.0,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' RS-fit cut-off ?',1,rsmcut)
      if (lstats) call outlie (j,xdum,rsmcut,1)
c
      if (lstats) then
c
        nstats = nstats + 1
        irsfmc = nstats
        stats (1,nstats) = rsmcut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
 6521 continue
      call errcon ('Sorry - real-space CC (main) NOT included')
      lrsm = .false.
      goto 2000
c
c ... RS-FIT VALUES (side chain atoms)
c
 6530 continue
c
      filnam = 'rsfit_sc.o'
      call textin (
     +  ' O data block with RS-fit values   ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening datablock file')
        goto 6531
      end if
      call skipem (iunit)
c
      read (iunit,'(a)',err=6531,end=6531) line
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (nopt .lt. 4 .or. ierr .ne. 0) goto 6531
c
      call prdb (optpar(1),optpar(2),optpar(3),optpar(4))
      call str2i (optpar(3),ires,ierr)
      if (ierr .ne. 0) goto 6531
c
      call odbchk (molnam,'_residue_rsfit',optpar(1),
     +             'r',optpar(2),nres,ires,lignore)
      if (lignore) goto 6531
c
      if (ires .gt. maxres) then
        call errcon ('Too many residues')
        goto 6531
      end if
c
      read (iunit,optpar(4),err=6531,end=6531)
     +    (rsfits(i),i=1,ires)
      close (iunit)
c
      if (lstats) then
        j=0
        do i=1,nres
          if (rsfits(i) .ge. -1.0 .and.
     +        restyp(i) .ne. 'GLY' .and.
     +        (.not. swater(i)) .and.
     +        abs(rsfits(i)) .ge. 0.0001) then
            j = j + 1
            xdum (j) = rsfits(i)
          end if
        end do
        nhisto = 20
        do i=1,nhisto
          histo(i) = 0.05 * float(i-1)
        end do
        call oopsts ('RS-fit values (side-chain; non-GLY)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        rsscut = max(xmin+0.1,min(xmax-0.1,ave-2.0*sdv))
c
        filnam = molnam(1:leng1(molnam))//'_rsfit_sc.plt'
        title = 'RS-fit values (side-chain atoms)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,j,xdum,min(0.0,xmin),1.0,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' RS-fit cut-off ?',1,rsscut)
      if (lstats) call outlie (j,xdum,rsscut,1)
c
      if (lstats) then
c
        nstats = nstats + 1
        irsfsc = nstats
        stats (1,nstats) = rsscut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
 6531 continue
      call errcon ('Sorry - real-space CC (side) NOT included')
      lrss = .false.
      goto 2000
c
c ... LOW/HIGH TEMPERATURE FACTORS
c
 6570 continue
c
      if (lstats) then
        nhisto = maxhis
        histo (1) = 1.0
        histo (2) = 2.0
        histo (3) = 5.0
        do i=4,nhisto
          histo(i) = 10.0 * float(i-3)
        end do
        call oopsts ('Temperature factors (all ATOMS !)',
     +               natoms,bfac,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
      end if
c
      write (*,*)
      call fvalin (' Threshold for low Bs  ?',1,lobfac)
      call fvalin (' Threshold for high Bs ?',1,hibfac)
      call fvalin (' Threshold for high Bs WATERs ?',1,hibwat)
c
      write (*,*) 'Checking low Bs ...'
      call lowbs ()
c
      write (*,*) 'Checking high Bs ...'
      call highbs ()
c
      goto 2000
c
c ... LOW/HIGH OCCUPANCIES
c
 6600 continue
c
      if (lstats) then
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 0.1 * float(i-1)
        end do
        call oopsts ('Occupancies (all ATOMS !)',
     +               natoms,qocc,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
      end if
c
      write (*,*)
      call fvalin (' Threshold for low Qs ?',1,loqocc)
      call fvalin (' Threshold for high Qs ?',1,hiqocc)
c
      write (*,*) 'Checking low Qs ...'
      call lowqs ()
c
      write (*,*) 'Checking high Qs ...'
      call highqs ()
c
      goto 2000
c
c ... MASK TOO TIGHT
c
 6560 continue
c
      filnam = molnam(1:leng1(molnam))//'.mask'
      call locase (filnam)
      call textin (' Mask file  ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening mask file')
        goto 6561
      end if
c
      call maskin (iunit,mask,origin,extent,grid,cell,maxbuf,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading mask file')
        goto 6561
      end if
c
      write (*,*)
      call fvalin (' Radius to check around atoms ?',1,mskrad)
c
      write (*,*) 'Checking if mask is too tight'
      call tightm ()
c
      goto 2000
c
 6561 continue
      call errcon ('Sorry - mask fit NOT included')
      lmask = .false.
      goto 2000
c
c ... DISULFIDE BRIDGES
c
 6645 continue
c
c ... check disulfides
c
      write (*,*) 'Checking disulfides'
      call disulf ()
c
      goto 2000
c
c ... COMPARE WITH PREVIOUS MODEL
c
 6650 continue
c
c ... calc Phi/Psi angles if not already done
c
      if (.not. lrama) call phipsi ()
c
c ... get PDB file of previous model
c
      write (*,*)
      filnam = molnam(1:leng1(molnam))//'_prev.pdb'
      call locase (filnam)
      call textin (
     +  ' PDB file of previous model ?',filnam)
      call xopxoa (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file')
        goto 6651
      end if
c
      call getold (iunit,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading PDB file')
        goto 6651
      end if
c
      write (*,*)
      prevnm = 'Prev Model'
      call textin (' Name of this previous model ?',prevnm)
      reply = 'N'
      call textin (' Include WATERs in comparison ?',reply)
      call upcase (reply)
      lprwat = (reply .ne. 'N')
c
      write (*,*)
      write (*,*) 'Comparing with previous model'
      call doprev (lprwat)
c
      if (lstats) then
c
c ... RMSD
c
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 0.1 * float(i-1)
        end do
        call oopsts ('RMSD current/previous model',
     +               nres,xxrmsd,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_rmsd.plt'
        title = 'RMSD current/previous model'
        optpar (1) = 'RMSD for each residue (A)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,nres,xxrmsd,0.0,1.05*xmax+0.01,nres,
     +    ave,sdv,xmin,xmax)
c
c ... RMS delta-B
c
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = float(i-1)
        end do
        call oopsts ('RMS delta-B current/previous model',
     +               nres,xxrmsb,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_rmsb.plt'
        title = 'RMS delta-B current/previous model (A**2)'
        optpar (1) = 'RMS delta-B for each residue (A**2)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,nres,xxrmsb,0.0,1.05*xmax+0.01,nres,
     +    ave,sdv,xmin,xmax)
c
c ... RMS delta-Q
c
        nhisto = 15
        do i=1,nhisto
          histo(i) = 0.1 * float(i-1)
        end do
        call oopsts ('RMS delta-Q current/previous model',
     +               nres,xxrmsq,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_rmsq.plt'
        title = 'RMS delta-Q current/previous model'
        optpar (1) = 'RMS delta-Q for each residue'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,nres,xxrmsq,0.0,1.05*xmax+0.01,nres,
     +    ave,sdv,xmin,xmax)
c
c ... PHI/PSI distance
c
        nhisto = 19
        do i=1,nhisto
          histo(i) = 10.0 * float(i-1)
        end do
        call oopsts ('Phi/Psi distance',
     +               nres,xxphps,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_phipsi.plt'
        title = 'Phi/Psi distance'
        optpar (1) = 'Phi/Psi distance for each residue'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,nres,xxphps,0.0,180.0,nres,
     +    ave,sdv,xmin,xmax)
c
c ... CHI1/2 distance
c
        nhisto = 19
        do i=1,nhisto
          histo(i) = 10.0 * float(i-1)
        end do
        call oopsts ('Chi1/2 distance',
     +               nres,xxchid,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_chi12.plt'
        title = 'Chi1/2 distance'
        optpar (1) = 'Chi1/2 distance for each residue'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,nres,xxchid,0.0,180.0,nres,
     +    ave,sdv,xmin,xmax)
c
      end if
c
      write (*,*)
      call fvalin (' Maximum RMSD (A)         ?',1,mxrmsd)
      call fvalin (' Maximum RMS delta-B (A2) ?',1,mxrmsb)
      call fvalin (' Maximum RMS delta-Q (A2) ?',1,mxrmsq)
      call fvalin (' Maximum Phi/Psi distance ?',1,mxphps)
      call fvalin (' Maximum Chi1/2 distance  ?',1,mxchid)
c
      reply = 'N'
      call textin (' Do you want an O datablock file ?',reply)
      call upcase (reply)
      if (reply .eq. 'Y') call prevdb (iunit)
c
      goto 2000
c
 6651 continue
      call errcon ('Sorry - previous model NOT included')
      lprev = .false.
      goto 2000
c
c ... C-ALPHA CHIRALITY
c
 6640 continue
c
      write (*,*) 'Calculating improper twist angle'
      write (*,*) 'CA(i) - N(i) - C(i) - CB(i)'
      write (*,*) 'This angle should be +33.9, sigma = 3.5,'
      write (*,*) 'for standard L-amino acids'
      call dozeta ()
c
      if (lstats) then
        j=0
        do i=1,nres
          if (zeta(i) .le. 999.) then
            j=j+1
            xdum (j) = zeta(i)
          end if
        end do
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 33.9 + 1.0 * float(i - maxhis/2 - 1)
        end do
        call oopsts ('C-alpha chirality (non-GLY/PRO)',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_ca_chir.plt'
        title = 'C-alpha chirality angles (ideal 33.9; sigma 3.5)'
        optpar (1) = 'Improper CA(i) - N(i) - C(i) - CB(i)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +      title,j,xdum,xmin-1.0,xmax+1.0,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' Maximum absolute deviation ?',1,chicut)
c
      if (lstats) then
c
        nstats = nstats + 1
        ichirl = nstats
        stats (1,nstats) = chicut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000
c
c ... PEPTIDE PLANARITY
c
 6630 continue
c
      write (*,*) 'Calculating improper twist angle'
      write (*,*) 'C(i) - CA(i) - N(i+1) - O(i)'
      write (*,*) 'This angle should be ZERO, sigma = 5.8'
      write (*,*)
      write (*,*) 'Also calculating OMEGA torsion angle'
      write (*,*) 'CA(i) - C(i) - N(i+1) - CA(i+1)'
      write (*,*) 'This angle should be +178, sigma = 6'
      write (*,*) 'for TRANS peptide links.'
      call plapep ()
c
      if (lstats) then
        j=0
        do i=1,nres
          if (plangl(i) .le. 999.) then
            j = j + 1
            xdum (j) = plangl(i)
          end if
        end do
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = 1.0 * float(i - maxhis/2 - 1)
        end do
        call oopsts ('Peptide planarity',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_pep_plan.plt'
        title = 'Peptide-planarity angles (ideal 0.0; sigma 5.8)'
        optpar (1) = 'Improper C(i) - CA(i) - N(i+1) - O(i)'
        if (lplot) call ooplot (iunit,filnam,optpar(1),
     +    title,j,xdum,xmin-1.0,xmax+1.0,j,ave,sdv,xmin,xmax)
      end if
c
      write (*,*)
      call fvalin (' Max deviation improper ?',1,placut)
      write (*,*)
      call prompt (' Suggested OMEGA cut-off for CNS 3-6 degrees')
      call prompt (' For REFMAC (and probably others) 10-12 degrees')
      call fvalin (' Max deviation OMEGA    ?',1,omecut)
c
      if (lstats) then
c
        nstats = nstats + 1
        iplanp = nstats
        stats (1,nstats) = placut
        stats (2,nstats) = ave
        stats (3,nstats) = sdv
        stats (4,nstats) = xmin
        stats (5,nstats) = xmax
        stats (6,nstats) = j
c
      end if
c
      goto 2000




c
c ... RMS TEMP. FACTORS BONDED ATOMS *** NOT USED ***
c
 6590 continue
      reply = 'N'
      write (*,*)
      call textin (' Analyse RMS delta-B bonded atoms ?',reply)
      call upcase (reply)
      ldbb = (reply .ne. 'N')
c
      if (ldbb) then
c
        call fvalin (' Maximum distance for bonded atoms ?',
     +    1,mxatat)
c
        call rdbb ()
c
        if (lstats) then
          j = nres
          nhisto = maxhis
          do i=1,nhisto
            histo(i) = 1.0 * float(i-1)
          end do
          call oopsts ('RMS delta-B bonded atoms',
     +                 nres,rmsdbb,nhisto,histo,nhcnt,
     +                 ave,sdv,xmin,xmax,xtot)
c
          filnam = molnam(1:leng1(molnam))//'_rmsdbb.plt'
          title = 'RMS delta-B bonded atoms'
          if (lplot) call ooplot (iunit,filnam,optpar(1),
     +      title,nres,rmsdbb,0.0,1.05*xmax,j,ave,sdv,xmin,xmax)
        end if
c
        write (*,*)
        call fvalin (' Threshold for RMS delta-B bonded atoms ?',
     +    1,dbbcut)
c
      end if
c
c ... QUALWAT VALUES *** NOT USED ***
c
 6660 continue
      reply = 'N'
      write (*,*)
      call textin (' Analyse QualWat values ?',reply)
      call upcase (reply)
      lwat = (reply .ne. 'N')
c
      if (lwat) then
c
        if (nwater .le. 0) then
          call errcon ('No water molecules found')
          lwat = .false.
          goto 6912
        end if
c
        call fvalin (' Resolution limit of your data (A) ?',
     +    1,resoln)
        resoln = max (resoln, 0.1)
c
        write (*,*) 'Calculating QualWat'
        write (*,*) ' = 100 * Q * EXP(-B/(4D^2))'
        write (*,*) 'This is 0 for absent and 100 for perfect'
        write (*,*) 'water molecules'
        call dowatq ()
c
        if (lstats) then
          j=0
          do i=1,nres
            if (swater(i)) then
              j = j + 1
              xdum (j) = watqua (i)
            end if
          end do
          nhisto = maxhis
          do i=1,nhisto
            histo(i) = 10.0 * float(i-1)
          end do
          call oopsts ('QualWat values (waters)',
     +                 j,xdum,nhisto,histo,nhcnt,
     +                 ave,sdv,xmin,xmax,xtot)
c
          filnam = molnam(1:leng1(molnam))//'_qualwat.plt'
          optpar (1) = 'QualWat = 100 * Q * EXP(-B/(4D^2))'
          title = 'QualWat (0 = absent, 100 = perfect water)'
          if (lplot) call ooplot (iunit,filnam,optpar(1),
     +      title,j,xdum,0.0,100.0,j,ave,sdv,xmin,xmax)
        end if
c
        write (*,*)
        call fvalin (' QualWat cut-off ?',1,watcut)
c
 6912   continue
c
      end if
c
c ... NR OF BAD CONTACTS *** NOT USED ***
c
 6670 continue
      reply = 'N'
      write (*,*)
      call textin (' Analyse nr of bad contacts ?',reply)
      call upcase (reply)
      lbadco = (reply .ne. 'N')
c
      if (lbadco) then
c
        filnam = 'badcon.o'
        call textin (
     +    ' O data block with contact counts  ?',filnam)
        call xopxoa (iunit,filnam,linter,ierr)
        if (ierr .ne. 0) then
          if (linter) goto 6670
          call errstp ('While opening datablock file')
        end if
        call skipem (iunit)
c
        read (iunit,'(a)',err=999,end=999) line
        call extrop (line,nopt,maxopt,optpar,ierr)
        if (nopt .lt. 4 .or. ierr .ne. 0) goto 999
c
        call textut (' Datablock :',optpar(1))
        call textut (' Data type :',optpar(2))
        call textut (' Number    :',optpar(3))
        call textut (' Format    :',optpar(4))
        call str2i (optpar(3),ires,ierr)
        if (ierr .ne. 0) goto 999
c
        call odbchk (molnam,'_residue_badcon',optpar(1),
     +               'i',optpar(2),nres,ires,lignore)
ccc      if (lignore) goto 6551
c
        if (ires .gt. maxres) then
          call errstp ('Too many residues')
        end if
c
        read (iunit,optpar(4),err=999,end=999)
     +    (nbadco(i),i=1,ires)
        close (iunit)
c
        if (lstats) then
          do i=1,nres
            xdum (i) = float (nbadco(i))
          end do
          nhisto = 11
          do i=1,nhisto
            histo(i) = float(i-1)
          end do
          call oopsts ('Nr of bad contacts',
     +                 nres,xdum,nhisto,histo,nhcnt,
     +                 ave,sdv,xmin,xmax,xtot)
c
          filnam = molnam(1:leng1(molnam))//'_badcon.plt'
          title = 'Nr of bad contacts'
          if (lplot) call ooplot (iunit,filnam,optpar(1),
     +      title,nres,xdum,0.0,1.05*xmax,j,ave,sdv,xmin,xmax)
        end if
c
        write (*,*)
        call ivalin (' Bad contacts cut-off ?',1,nbadcu)
c
      end if



c
c ... GO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
 3500 continue
c
c ... USER-DEFINABLE CRITERIA
c
      nuser = 0
      write (*,*)
      call prompt (' User-definable criteria')
      call ivalut (' Max number of them :',1,maxuse)
c
      do iuser=1,maxuse
c
 6990   continue
        write (*,*)
        filnam = ' '
        call textin (
     +    ' Enter file with user datablock (<CR> to stop):',
     +    filnam)
        if (length(filnam) .lt. 1) goto 3004
c
        call xopxoa (iunit,filnam,linter,ierr)
        if (ierr .ne. 0) then
          if (linter) goto 6990
          call errstp ('While opening datablock file')
        end if
c
        nuser = nuser + 1
c
        call skipem (iunit)
        read (iunit,'(a)',err=999,end=999) line
        call extrop (line,nopt,maxopt,optpar,ierr)
        if (nopt .lt. 4 .or. ierr .ne. 0) goto 999
c
        call textut (' Datablock :',optpar(1))
        call textut (' Data type :',optpar(2))
        call textut (' Number    :',optpar(3))
        call textut (' Format    :',optpar(4))
        call str2i (optpar(3),ires,ierr)
        if (ierr .ne. 0) goto 999
c
        if (ires .gt. maxres) then
          call errstp ('Too many residues')
        end if
c
        call upcase (optpar(2))
        if (optpar(2)(1:1) .eq. 'R') then
          read (iunit,optpar(4),err=999,end=999)
     +      (usercr(i,nuser),i=1,ires)
        else if (optpar(2)(1:1) .eq. 'I') then
          read (iunit,optpar(4),err=999,end=999)
     +      (mask(i),i=1,ires)
          do i=1,ires
            usercr (i,nuser) = float(mask(i))
          end do
        else
          call errstp ('Datablock MUST be of type I or R !')
        end if
c
        close (iunit)
c
        call textin (' Name of the property ?',usernm(nuser))
c
        userln (nuser) = ' Bad ' //
     +                   usernm(nuser)(1:leng1(usernm(nuser))) //
     +                   '; value = '
        call pretty (userln(nuser))
        call textin (
     +    ' Warning for bad residues in O ?',userln (nuser))
c
        if (lstats) then
c
          xdum (1) = usercr(1,nuser)
          xdum (2) = usercr(1,nuser)
          do i=2,ires
            xdum(1) = min ( xdum(1), usercr(i,nuser) )
            xdum(2) = max ( xdum(2), usercr(i,nuser) )
          end do
          xdum (3) = 0.1 * (xdum(2) - xdum(1))
          xdum (1) = xdum(1) - xdum(3)
          nhisto = 12
          do i=1,nhisto
            histo (i) = xdum(1) + float(i-1)*xdum(3)
          end do
          call oopsts (usernm(nuser),
     +                 ires,usercr(1,nuser),nhisto,histo,nhcnt,
     +                 ave,sdv,xmin,xmax,xtot)
c
        end if
c
        write (*,'(/a,a,a/a)')
     +    ' A residue has a BAD "',
     +    usernm(nuser)(1:leng1(usernm(nuser))),'" IF:',
     +    ' its value either > CUTOFF or < CUTOFF'
c
        call textin (' Should I use ">" or "<" ?',usersn(nuser))
        if (usersn(nuser) .ne. '<' .and.
     +      usersn(nuser) .ne. '>') then
          call errstp ('Invalid comparison operator !')
        end if
c
        if (lstats) then
          if (usersn(nuser) .eq. '>') then
            usercu(nuser) = ave + 2.0 * sdv
          else
            usercu(nuser) = ave - 2.0 * sdv
          end if
        end if
c
        call rvalin (' Cut-off value to use ?',1,usercu(nuser))
c
        call prompt (' Checking ...')
        if (usersn(nuser) .eq. '>') then
          do i=1,ires
            userba (i,nuser) = (usercr(i,nuser) .gt. usercu(nuser))
          end do
        else
          do i=1,ires
            userba (i,nuser) = (usercr(i,nuser) .lt. usercu(nuser))
          end do
        endif
c
      end do
c
c ... END OF USER-DEFINABLE CRITERIA
c
 3004 continue
      write (*,*)
      call ivalut (' Nr of user criteria :',1,nuser)
      write (*,*)
c
c ... ANALYSE THEM
c
      lfirst = .true.
      nbad = 0
      nokay = 0
      write (*,*)
c
c      reply = 'Y'
c      call prompt (' You may opt to get the details listed')
c      call prompt (' on the screen.')
c      call textin (' Do you want to see the details ?',reply)
c      call upcase (reply)
c      listem = (reply .eq. 'Y')
c      write (*,*)
c
      listem = .true.
c
c      reply = 'Y'
c      call prompt (' Also, you can get this list written to a file;')
c      call prompt (' this is very handy for electronic note-keeping:')
c      call prompt (' just edit the file as you rebuild your model &')
c      call prompt (' print it & stick it in your notebook !')
c      call textin (' Do you want to have a list file ?',reply)
c      call upcase (reply)
c      lnotes = (reply .eq. 'Y')
c
      lnotes = .true.
      if (lnotes) then
        notfil = molnam // '_rebuild.notes'
        call remspa (notfil)
        call locase (notfil)
ccc        call textin (' Name of the list file ?',notfil)
        call textut (' Electronic notebook file :',notfil)
        call xopxua (junit,notfil,.false.,ierr)
        if (ierr .ne. 0) lnotes = .false.
      end if
      write (*,*)
c
c      reply = 'Y'
c      call prompt (' In addition, you can get an HTML file with a')
c      call prompt (' residue-by-residue critique of your current')
c      call prompt (' model.')
c      call textin (' Do you want to have an HTML file ?',reply)
c      call upcase (reply)
c      lhtml = (reply .eq. 'Y')
c
      lhtml = .true.
      if (lhtml) then
        html = molnam // '_oops.html'
        call remspa (html)
        call locase (html)
ccc        call textin (' Name of the HTML file ?',html)
        call textut (' HTML model critique file :',html)
        call xopxua (kunit,html,.false.,ierr)
        if (ierr .ne. 0) then
          lhtml = .false.
        else
          myline = '<html><head><title>OOPS Critique of model ' //
     +      molnam // '</title>'
          call putlin (kunit,myline,.false.,ierr)
         myline = '</head><body bgcolor="#fefefe">'
          call putlin (kunit,myline,.false.,ierr)
           myline = '<h1>OOPS Critique of model ' //
     +      molnam // '</h1>'
          call putlin (kunit,myline,.false.,ierr)
          myline = '<table border=3 cellspacing=3> <tr> ' //
     +      '<th>Residue</th><th>Remarks</th></tr>'
          call putlin (kunit,myline,.false.,ierr)
        end if
      end if
      write (*,*)
c
c      if (lrama .and. lpep) then
c        reply = 'N'
c        call prompt (' You can also get a "pseudo-PDB file" with')
c        call prompt (' one CA atom for each residue which has:')
c        call prompt ('   X = Phi(residue)')
c        call prompt ('   Y = Psi')
c        call prompt ('   Z = Pep-flip * 100')
c        call prompt ('   B = Maximum temp. factor')
c        call prompt ('   Q = RSC')
c        call prompt (' Displaying and colouring this "molecule"')
c        call prompt (' in O may reveal poor areas and correlations')
c        call prompt (' between poor Phi/Psi and Pep-flip values')
c        call textin (' Create pseudo-PDB file ?',reply)
c        call upcase (reply)
c        if (reply .eq. 'Y') then
c          notfil = molnam // '_pseudo.pdb'
c          call remspa (notfil)
c          call locase (notfil)
c          call textin (' Name of the PDB file ?',notfil)
c          call xopxua (iunit,notfil,.false.,ierr)
c          if (ierr .eq. 0) then
c            call pseudo (iunit)
c            notfil = molnam // '_pseudo.odb'
c            call remspa (notfil)
c            call locase (notfil)
c            call textin (' Name of the O datablock file ?',notfil)
c            call xopxua (iunit,notfil,.false.,ierr)
c            if (ierr .eq. 0) then
c              call psodb (iunit)
c            end if
c          end if
c        end if
c        write (*,*)
c      end if
c
      call prompt (
     + ' You may enter one line of O commands which will be')
      call prompt (
     + ' executed by EVERY macro generated by OOPS.')
      call textin (' O command(s) to execute in every macro ?',oline)
      write (*,*)
c
      reply = 'N'
      call prompt (
     + ' You may opt to get macros for ALL or merely the')
      call prompt (
     + ' BAD residues.')
      call textin (' Do you want macros for ALL residues ?',reply)
      call upcase (reply)
      lalres = (reply .eq. 'Y')
      write (*,*)
c
      if (.not. lalres) then
        reply = 'N'
        call prompt (
     +    ' Alternatively, you may opt to get macros for all')
        call prompt (
     +    ' residues (good or bad) of a certain type (protein,')
        call prompt (
     +    ' nucleic acid, water, heterogen).')
        call textin (' Do you want macros for certain types ?',reply)
        call upcase (reply)
        laltyp (1) = .false.
        laltyp (2) = .false.
        laltyp (3) = .false.
        laltyp (4) = .false.
        if (reply .eq. 'Y') then
          reply = 'N'
          call textin (' Macros for all amino acids ?',reply)
          call upcase (reply)
          laltyp (1) = (reply .eq. 'Y')
c
          reply = 'N'
          call textin (' Macros for all nucleotides ?',reply)
          call upcase (reply)
          laltyp (2) = (reply .eq. 'Y')
c
          reply = 'N'
          call textin (' Macros for all waters      ?',reply)
          call upcase (reply)
          laltyp (3) = (reply .eq. 'Y')
c
          reply = 'N'
          call textin (' Macros for all heterogens  ?',reply)
          call upcase (reply)
          laltyp (4) = (reply .eq. 'Y')
        end if
        write (*,*)
      end if
c
      lxgood = .false.
      lipoor = .false.
c
      if ( (.not. lalres) .and. lknow) then
        reply = 'N'
        call prompt (
     + ' If you wish, residues that you have previously flagged')
        call prompt (
     + ' as "good density, good fit" can be skipped by OOPS,')
        call prompt (
     + ' even if they would normally be flagged as "bad" (e.g.,')
        call prompt (
     + ' due to an unusual side-chain conformation).')
        call textin (' Skip good residues ?',reply)
        call upcase (reply)
        lxgood = (reply .eq. 'Y')
        write (*,*)
        reply = 'N'
        call prompt (
     + ' Conversely, residues that you have previously flagged')
        call prompt (
     + ' as having poor density or a poor fit can be flagged by')
        call prompt (
     + ' OOPS, even if they are not flagged as "bad".')
        call textin (' Include poor residues ?',reply)
        call upcase (reply)
        lipoor = (reply .eq. 'Y')
        write (*,*)
      end if
c
      reply = 'Y'
      call prompt (
     + ' You may opt to get chained macros or individual macros')
      call prompt (
     + ' without instructions to put the next macro on the menu.')
      call textin (' Do you want CHAINED macros ?',reply)
      call upcase (reply)
      lchain = (reply .eq. 'Y')
      write (*,*)
c
      lworst = .false.
      if (lchain) then
        reply = 'N'
        call prompt (
     + ' The macros can be executed either from first residue to')
        call prompt (
     + ' last, or starting with the residue with the highest')
        call prompt (
     + ' number of possible problems ("worst-first").')
        call textin (' Do you want the macros "WORST-FIRST" ?',reply)
        call upcase (reply)
        lworst = (reply .eq. 'Y')
        write (*,*)
      end if
c
      reply = 'Y'
      call prompt (
     + ' From O version 9.0.7 on, you can opt to use the Grab_build')
      call prompt (
     + ' panel for every residue flagged by OOPS2.')
      call textin (' Do you want to use GRAB_BUILD ?',reply)
      call upcase (reply)
      lgrab = (reply .eq. 'Y')
      write (*,*)
c
c      reply = 'Y'
c      call prompt (' You may opt to give the macros the')
c      call prompt (' same name as the residue which they')
c      call prompt (' deal with (in lowercase), or to just')
c      call prompt (' give them sequential numbers (1, 2, ...).')
c      call textin (' Do you want macros named as RESIDUES ?',reply)
c      call upcase (reply)
c      lmanam = (reply .eq. 'Y')
c      write (*,*)
c
      lmanam = .true.
c
c ... subjective judgement datablock file
c
      if (lsubj) then
        dbsub = molnam
        call appstr (dbsub,'_residue_quality')
        rnamdb = molnam
        call appstr (rnamdb,'_residue_name')
        rtypdb = molnam
        call appstr (rtypdb,'_residue_type')
        line = 'oops/oops.odb'
        call xopxua (iunit,line,.false.,ierr)
        if (ierr .ne. 0) then
          lsubj = .false.
        else
          write (iunit,'(4(a/a,a,a/a/a/))')
     +    '@Good_fit t 3 60',
     +    'db_set_dat ',dbsub,' $oops_irc $oops_irc 0',
     +    'print QUALITY : Good fit, good density',
     +    'message QUALITY : Good fit, good density',
     +    '@Poor_fit t 3 60',
     +    'db_set_dat ',dbsub,' $oops_irc $oops_irc 1',
     +    'print QUALITY : Poor fit, good density',
     +    'message QUALITY : Poor fit, good density',
     +    '@Poor_dens t 3 60',
     +    'db_set_dat ',dbsub,' $oops_irc $oops_irc 2',
     +    'print QUALITY : Poor density, cannot fit better',
     +    'message QUALITY : Poor density, cannot fit better',
     +    '@No_dens t 3 60',
     +    'db_set_dat ',dbsub,' $oops_irc $oops_irc 3',
     +    'print QUALITY : No density, cannot fit',
     +    'message QUALITY : No density, cannot fit'
        end if
        close (iunit)
      end if
c
      nmacro = 0
      nxgood = 0
      nipoor = 0
c
c ... merged low & high Bs and Qs, so flags should be the same
c
      lhib = llob
      lhiq = lloq
c
      do i=1,nres
c
c ... check if protein=1, nucleic acid=2, water=3, heterogen=4
c
        moltyp = 1
        call nuctyp (restyp(i),itype)
        if (itype .gt. 0) then
          moltyp = 2
        else if (swater(i)) then
          moltyp = 3
        else if (sstype(i) .lt. 0) then
          moltyp = 4
        end if
        cnttyp (moltyp) = cnttyp (moltyp) + 1
c
        badpep = .false.
        badrsf = .false.
        badrsm = .false.
        badrss = .false.
        badram = .false.
        badpro = .false.
        badpos = .false.
        badlef = .false.
        badmas = .false.
        badblo = .false.
        badbhi = .false.
        badqlo = .false.
        badqhi = .false.
        badrsc = .false.
        badpla = .false.
        badome = .false.
        badcis = .false.
        badzet = .false.
        baddaa = .false.
        badisu = .false.
        badssd = .false.
        badsst = .false.
        badwat = .false.
        badcon = .false.
        badrsr = .false.
        baddbb = .false.
        baddis = .false.
        badtem = .false.
        badocc = .false.
        badphi = .false.
        badchi = .false.
        badwif = .false.
        badpor = .false.
c
        do j=1,nuser
          baduse(j) = .false.
        end do
c
c        write (*,'(1x,a6,2(2x,f10.3))') resnam(i),pepflip(i),
c     +    rsfit(i)
c
        if (lpep) then
          if (abs(pepflip(i)) .gt. 0.001) then
            if (pepflip(i) .ge. pepcut) then
              badpep = .true.
              count (i) = count (i) + 1
            end if
          end if
        end if
c
        if (lrsfit) then
          if (abs(rsfit(i)) .gt. 0.001) then
            if (swater(i)) then
              if (rsfit(i) .le. rsfwat) then
                badrsf = .true.
                count (i) = count (i) + 1
              end if
            else
              if (rsfit(i) .le. rsfcut) then
                badrsf = .true.
                count (i) = count (i) + 1
              end if
            end if
          end if
        end if
c
        if (lrsm) then
          if (abs(rsfitm(i)) .gt. 0.001 .and.
     +        (.not. swater(i)) ) then
            if (rsfitm(i) .le. rsmcut) then
              badrsm = .true.
              count (i) = count (i) + 1
            end if
          end if
        end if
c
        if (lrss) then
          if (abs(rsfits(i)) .gt. 0.001 .and.
     +        (.not. swater(i)) ) then
            if (rsfits(i) .le. rsscut .and.
     +          restyp(i) .ne. 'GLY') then
              badrss = .true.
              count (i) = count (i) + 1
            end if
          end if
        end if
c
        if (lrsr) then
          if (abs(rsrfac(i)) .gt. 0.001) then
            if (swater(i)) then
              if (rsrfac(i) .ge. rsrwat) then
                badrsr = .true.
                count (i) = count (i) + 1
              end if
            else
              if (rsrfac(i) .ge. rsrcut) then
                badrsr = .true.
                count (i) = count (i) + 1
              end if
            end if
          end if
        end if
c
        if (lrsc) then
          if (abs(rsc(i)) .gt. 0.001) then
            if (rsc(i) .ge. rsccut .and.
     +          restyp(i) .ne. 'GLY' .and.
     +          restyp(i) .ne. 'ALA') then
              badrsc = .true.
              count (i) = count (i) + 1
            end if
          end if
        end if
c
        if (lrama) then
          badram = badpp(i)
          if (badram) count (i) = count (i) + 1
c
          if (lpro .and. restyp(i) .eq. 'PRO' .and.
     +        phi(i) .lt. 999.0) then
            badpro = ( abs(phi(i)+65.4) .gt. 11.2 )
c
c ... don't count this as an error !
c
ccc            if (badpro .and. (.not. badram)) count(i)=count(i)+1
c
          end if
c
          if (lpos .and. restyp(i) .ne. 'GLY' .and.
     +        phi(i) .lt. 999.0) then
            badpos = ( phi(i) .gt. 0.0 )
c
c ... don't count this as an error !
c
ccc            if (badpos .and. (.not. badram)) count(i)=count(i)+1
c
          end if
c
          if (llefth) then
            badlef = ( sstype(i) .eq. 3 )
c
c ... don't count this as an error !
c
          end if
c
        end if
c
        if (lmask) then
          badmas = badmsk(i)
          if (badmas) count (i) = count (i) + 1
        end if
c
        if (llob) then
          badblo = badlob (i)
          if (badblo) count (i) = count (i) + 1
        end if
c
        if (lhib) then
          badbhi = badhib (i)
          if (badbhi) count (i) = count (i) + 1
        end if
c
        if (ldbb) then
          baddbb = (rmsdbb(i) .ge. dbbcut)
          if (baddbb) count (i) = count (i) + 1
        end if
c
        if (lloq) then
          badqlo = badloq (i)
          if (badqlo) count (i) = count (i) + 1
        end if
c
        if (lhiq) then
          badqhi = badhiq (i)
          if (badqhi) count (i) = count (i) + 1
        end if
c
        if (lplan) then
          if (plangl(i) .le. 999.) then
            badpla = (abs(plangl(i)) .ge. placut)
            if (badpla) count (i) = count (i) + 1
          end if
          if (omega(i) .le. 999.) then
            badome = ( min(abs(omega(i)),abs(omega(i)-180.0),
     +                     abs(omega(i)+180.0)) .ge. omecut)
            if (badome) count (i) = count (i) + 1
            badcis = (abs(omega(i)) .le. 90.0 .or.
     +                abs(omega(i)) .ge. 270.0)
          end if
        end if
c
        if (lchir) then
          if (zeta(i) .le. 999.) then
            badzet = (abs(zeta(i)-33.9) .ge. chicut)
            if (badzet) count (i) = count (i) + 1
            baddaa = (zeta(i) .lt. 0.0)
          end if
        end if
c
        if (ldisu .and. ndisu.gt. 0 .and.
     +      restyp(i) .eq. 'CYS') then
          idisu = -1
          do j=1,ndisu
            if (ssresi(1,j) .eq. i .or. ssresi(2,j) .eq. i) then
              badisu = .true.
              idisu = j
              badssd = (abs(ssdist(j)-2.06) .gt. 0.1)
              if (badssd) count(i) = count(i) + 1
              badsst = (abs(sstors(j)-96.8) .gt. 10.1 .and.
     +          abs(sstors(j)+85.8) .gt. 8.6)
              if (badsst) count(i) = count(i) + 1
              goto 7382
            end if
          end do
 7382     continue
        end if
c
        if (lprev) then
          if ( (.not.lprwat) .or.
     +         (.not. swater(i)) ) then
            if (.not. (insert(i).or.mutate(i))) then
              if (xxrmsd(i).ge.mxrmsd) then
                baddis = .true.
                count(i) = count(i) + 1
              end if
              if (xxrmsb(i).ge.mxrmsb) then
                badtem = .true.
                count(i) = count(i) + 1
              end if
              if (xxrmsq(i).ge.mxrmsq) then
                badocc = .true.
                count(i) = count(i) + 1
              end if
              if (xxphps(i).ge.mxphps) then
                badphi = .true.
                count(i) = count(i) + 1
              end if
              if (xxchid(i).ge.mxchid) then
                badchi = .true.
                count(i) = count(i) + 1
              end if
            end if
          end if
        end if
c
        if (lwat .and. swater(i)) then
          badwat = (watqua(i) .le. watcut)
          if (badwat) count (i) = count (i) + 1
        end if
c
        if (lbadco) then
          badcon = (nbadco(i) .ge. nbadcu)
          if (badcon) count (i) = count (i) + 1
        end if
c
        if (lwif .and. nwif .gt. 0) then
          j = iindex (i,0,nwif,wifptr)
          badwif = (j .gt. 0)
          if (badwif) then
            do j=1,nwif
              if (wifptr(j) .eq. i) count (i) = count (i) + 1
            end do
          end if
        end if
c
        if (nuser .gt. 0) then
          do j=1,nuser
            baduse(j) = userba(i,j)
            if (baduse(j)) count (i) = count (i) + 1
          end do
        end if
c
        if (lipoor .and. oldsub(i) .gt. 0) then
          count (i) = count (i) + 1
          badpor = .true.
          nipoor = nipoor + 1
        end if
c
        baddy = (badpep .or. badrsf .or. badram .or. badmas .or.
     +           badblo .or. badbhi .or. badrsc .or. badrsm .or.
     +           badrss .or. badpla .or. badzet .or. badcon .or.
     +           badqlo .or. badqhi .or. badwat .or. badrsr .or.
     +           baddbb .or. baddis .or. badocc .or. badtem .or.
     +           mutate(i) .or. insert(i) .or. badphi .or.
     +           badchi .or. badome .or. badcis .or. baddaa .or.
     +           badpro .or. badpos .or. badisu .or. badssd .or.
     +           badsst .or. badwif .or. badpor .or. badlef)
c
        if (.not. baddy) then
          if (nuser .gt. 0) then
            do j=1,nuser
              baddy = ( baddy .or. baduse(j) )
            end do
          end if
        end if
c
        if (baddy .and. lxgood .and. oldsub(i) .eq. 0) then
          baddy = .false.
          count (i) = 0
          nxgood = nxgood + 1
          gktext = (restyp(i)//' '//resnam(i))
          call pretty (gktext)
          call textut (' SKIP - (good dens, good fit) -',gktext)
        end if
c
        if (baddy .or. lalres .or. laltyp(moltyp) ) then
c
          cntbad (moltyp) = cntbad (moltyp) + 1
c
          gktext = (restyp(i)//' '//resnam(i)//' ['//
     +      ssenam(sstype(i))(1:leng1(ssenam(sstype(i))))//'] ['//
     +      typnam(moltyp)(1:leng1(typnam(moltyp)))//']')
c
          call pretty (gktext)
          if (baddy) then
ccc            write (*,*)
            call textut (' OOPS -',gktext)
          else
ccc            write (*,*)
            call textut (' OKAY -',gktext)
          end if
c
          nbad = nbad + 1
          if (lmanam) then
            filnam = 'oops/'//resnam(i)
ccc            filnam = 'oops/'//restyp(i)//'_'//resnam(i)
            call remspa (filnam)
            call locase (filnam)
          else
            write (filnam,*) 'oops/',nbad
            call remspa (filnam)
          end if
c
c ... replace $ (for HETATM residues) by underscore in macro file name
c
          call subchr (filnam,'$','_',idum)
c
c ... keep track of things for "worst-first" option
c
          if (lworst) then
            bindex (nbad) = nbad
            bcount (nbad) = count(i)
            bfile (nbad) = filnam
          end if
c
          if (lfirst) then
c
            lfirst = .false.
            line = 'oops.omac'
            call xopxua (iunit,line,.false.,ierr)
            if (ierr .ne. 0) goto 998
c
c ... subjectivity stuff
c
            if (lsubj) then
              write (iunit,'(a/a,a,1x,i6,1x,a/a,a,a,1x,i6,1x,a)')
     +          ' read oops/oops.odb',
     +          ' db_create ',dbsub,nres,' i',
     +          ' db_set_dat ',dbsub,' 1 ',nres,' -1'
              write (iunit,'(a/a/a/a)')
     +          ' menu @Good_fit on',
     +          ' menu @Poor_fit on',
     +          ' menu @Poor_dens on',
     +          ' menu @No_dens on'
            end if
c
            call stamp (line)
            call pretty (line)
            write (iunit,'(a,a)') ' print ',line(1:leng1(line))
            if (lnotes) write (junit,*)
            if (lnotes) write (junit,'(a,a)') ' ',line(1:leng1(line))
c
            write (iunit,*) 'print Molecule ',
     +        molnam(1:leng1(molnam))
            if (lnotes) write (junit,*) 'Molecule ',
     +        molnam(1:leng1(molnam))
c
            write (iunit,'(a,a)') ' print OOPS has checked:'
            if (lnotes) write (junit,'(a,a)') ' OOPS has checked:'
c
            if (lpep) write (iunit,*)
     +        'print Pep-flip values; cutoff = ',
     +        pepcut
            if (lpep .and. lnotes) write (junit,*)
     +        'Pep-flip values; cutoff = ',
     +        pepcut
c
            if (lrsfit) write (iunit,'(a,f8.3,a,f8.3)')
     +        ' print RS-fit (all atoms); cutoff = ',
     +        rsfcut,' ; WATERs = ',rsfwat
            if (lrsfit .and. lnotes) write (junit,'(a,f8.3,a,f8.3)')
     +        ' RS-fit (all atoms); cutoff = ',
     +        rsfcut,' ; WATERs = ',rsfwat
c
            if (lrsm) write (iunit,*)
     +        'print RS-fit (main-chain atoms); cutoff = ',
     +        rsmcut
            if (lrsm .and. lnotes) write (junit,*)
     +        'RS-fit (main-chain atoms); cutoff = ',
     +        rsmcut
c
            if (lrss) write (iunit,*)
     +        'print RS-fit (side-chain atoms); cutoff = ',
     +        rsscut
            if (lrss .and. lnotes) write (junit,*)
     +        'RS-fit (side-chain atoms); cutoff = ',
     +        rsscut
c
            if (lrsr) write (iunit,'(a,f8.3,a,f8.3)')
     +        ' print RS R-factor (all atoms); cutoff = ',
     +        rsrcut,' ; WATERs = ',rsrwat
            if (lrsr .and. lnotes) write (junit,'(a,f8.3,a,f8.3)')
     +        ' RS R-factor (all atoms); cutoff = ',
     +        rsrcut,' ; WATERs = ',rsrwat
c
            if (lrsc) write (iunit,*)
     +        'print RSC values; cutoff = ',
     +        rsccut
            if (lrsc .and. lnotes) write (junit,*)
     +        'RSC values; cutoff = ',
     +        rsccut
c
            if (lmask) write (iunit,*)
     +        'print Mask violations'
            if (lmask .and. lnotes) write (junit,*)
     +        'Mask violations'
c
            if (llob) write (iunit,*)
     +        'print Too low temperature factors; cutoff = ',
     +        lobfac
            if (llob .and. lnotes) write (junit,*)
     +        'Too low temperature factors; cutoff = ',
     +        lobfac
c
            if (lhib) write (iunit,'(a,f8.3,a,f8.3)')
     +        ' print Too high temperature factors; cutoff = ',
     +        hibfac,' ; WATERs = ',hibwat
            if (lhib .and. lnotes) write (junit,'(a,f8.3,a,f8.3)')
     +        ' Too high temperature factors; cutoff = ',
     +        hibfac,' ; WATERs = ',hibwat
c
            if (ldbb) write (iunit,*)
     +        'print Too high RMS delta-B bonded atoms; cutoff = ',
     +        dbbcut
            if (ldbb .and. lnotes) write (junit,*)
     +        'Too high RMS delta-B bonded atoms; cutoff = ',
     +        dbbcut
c
            if (lloq) write (iunit,*)
     +        'print Too low occupancies; cutoff = ',
     +        loqocc
            if (lloq .and. lnotes) write (junit,*)
     +        'Too low occupancies; cutoff = ',
     +        loqocc
c
            if (lhiq) write (iunit,*)
     +        'print Too high occupancies; cutoff = ',
     +        hiqocc
            if (lhiq .and. lnotes) write (junit,*)
     +        'Too high occupancies; cutoff = ',
     +        hiqocc
c
            if (lrama) write (iunit,*)
     +        'print Phi-Psi angle combinations (Ramachandran)'
            if (lrama .and. lnotes) write (junit,*)
     +        'Phi-Psi angle combinations (Ramachandran)'
c
            if (lplan) write (iunit,*)
     +        'print Peptide planarity; cutoff = ',
     +        placut
            if (lplan .and. lnotes) write (junit,*)
     +        'Peptide planarity; cutoff = ',
     +        placut
c
            if (lplan) write (iunit,*)
     +        'print Omega torsion; cutoff = ',
     +        omecut
            if (lplan .and. lnotes) write (junit,*)
     +        'Omega torsion; cutoff = ',
     +        omecut
c
            if (lchir) write (iunit,*)
     +        'print CA chirality; cutoff = ',
     +        chicut
            if (lchir .and. lnotes) write (junit,*)
     +        'CA chirality; cutoff = ',
     +        chicut
c
            if (ldisu) write (iunit,*)
     +        'print Disulfide bridges'
            if (ldisu .and. lnotes) write (junit,*)
     +        'Disulfide bridges'
c
            if (lprev) then
              write (iunit,*)
     +          'print Comparison with previous model ',
     +          prevnm
              write (iunit,*)
     +          'print RMSD; cutoff = ',mxrmsd
              write (iunit,*)
     +          'print RMS delta-B; cutoff = ',mxrmsb
              write (iunit,*)
     +          'print RMS delta-Q; cutoff = ',mxrmsq
              write (iunit,*)
     +          'print Phi/Psi distance; cutoff = ',mxphps
              write (iunit,*)
     +          'print Chi1/2 distance; cutoff = ',mxchid
              if (lprwat) then
                write (iunit,*)
     +            'print WATERs included in comparison'
              else
                write (iunit,*)
     +            'print WATERs NOT included in comparison'
              end if
c
              if (lnotes) then
                write (junit,*)
     +            'Comparison with previous model ',
     +            prevnm
                write (junit,*)
     +            'RMSD; cutoff = ',mxrmsd
                write (junit,*)
     +            'RMS delta-B; cutoff = ',mxrmsb
                write (junit,*)
     +            'RMS delta-Q; cutoff = ',mxrmsq
                write (junit,*)
     +            'Phi/Psi distance; cutoff = ',mxphps
                write (junit,*)
     +            'Chi1/2 distance; cutoff = ',mxchid
                if (lprwat) then
                  write (junit,*)
     +              'WATERs included in comparison'
                else
                  write (junit,*)
     +              'WATERs NOT included in comparison'
                end if
              end if
            end if
c
            if (lwat) write (iunit,*)
     +        'print QualWat; cutoff = ',
     +        watcut
            if (lwat .and. lnotes) write (junit,*)
     +        'QualWat; cutoff = ',
     +        watcut
c
            if (lbadco) write (iunit,*)
     +        'print Number of bad contacts; cutoff = ',
     +        nbadcu
            if (lbadco .and. lnotes) write (junit,*)
     +        'Number of bad contacts; cutoff = ',
     +        nbadcu
c
            if (lwif) write (iunit,*)
     +        'print WHAT IF diagnostics'
            if (lwif .and. lnotes) write (junit,*)
     +        'WHAT IF diagnostics'
c
c ... user-defined criteria
c
            if (nuser .gt. 0) then
              do j=1,nuser
                write (gktext,*) 'print ',
     +            usernm(j)(1:leng1(usernm(j))),
     +            ' cut-off = ',usercu(j)
                call pretty (gktext)
                write (iunit,'(1x,a)') gktext(1:leng1(gktext))
                if (lnotes) 
     +            write (junit,'(1x,a)') gktext(7:leng1(gktext))
              end do
            end if
c
            if (lalres) then
              write (iunit,*)
     +          'print Macros were generated for ALL residues'
            else
              write (iunit,*)
     +          'print Macros were generated for BAD residues'
              if (laltyp(1)) write (iunit,*)
     +          'print and all amino acids'
              if (laltyp(2)) write (iunit,*)
     +          'print and all nucleotides'
              if (laltyp(3)) write (iunit,*)
     +          'print and all water molecules'
              if (laltyp(4)) write (iunit,*)
     +          'print and all heterogens'
            end if
c
            if (lchain) then
              write (iunit,*)
     +          'print Macros were chained together'
            else
              write (iunit,*)
     +          'print Macros were unchained'
            end if
c
            if (lworst) then
              write (iunit,*)
     +          'print Macros were chained "worst-first"'
            else
              write (iunit,*)
     +          'print Macros were chained sequentially'
            end if
c
            if (lmanam) then
              write (iunit,*)
     +          'print Macros were given residue names'
            else
              write (iunit,*)
     +          'print Macros were named 1, 2, 3, ...'
            end if
c
            if (lgrab) then
              write (iunit,*)
     +          'print Using the Grab_build panel for rebuilding'
              write (iunit,*)
     +          'symbol oops_c_atom centre_build'
              write (iunit,*)
     +          'symbol oops_c_zone centre_build'
            else
              write (iunit,*)
     +          'print Not using the Grab_build panel for rebuilding'
              write (iunit,*)
     +          'symbol oops_c_atom centre_atom'
              write (iunit,*)
     +          'symbol oops_c_zone centre_zone'
            end if
c
            write (iunit,*) 'mol ',molnam(1:leng1(molnam))
c
            if (.not. lworst) write (iunit,*)
     +        '@',filnam(1:leng1(filnam))
            close (iunit)
c
            if (lnotes) write (junit,*)
c
          else
c
ccc            write (iunit,*) oline(1:leng1(oline))
c
            if (lchain .and. (.not. lworst)) then
              if (baddy) then
                write (iunit,*) 'print Hit or type "@',
     +            filnam(1:leng1(filnam)),'" for next baddy'
              else
                write (iunit,*) 'print Hit or type "@',
     +            filnam(1:leng1(filnam)),'" for next residue'
              end if
c
              write (iunit,*)
     +          'menu @',filnam(1:leng1(filnam)),' on'
c     +          'menu @',filnam(1:leng1(filnam)),' on on_off'
              write (iunit,*)
     +          'menu @',line(1:leng1(line)),' off'
c     +          'menu @',line(1:leng1(line)),' off on_off'
            end if
c
c ... tell them to take a break every now and then
c
            nmacro = nmacro + 1
            if (nmacro .eq. 50 .and. (.not. lworst)) then
              write (iunit,*)
     +   'print You checked another 50 residues ... take a break !'
              write (iunit,*)
     +   'message You checked another 50 residues ... take a break !'
              write (iunit,*) 'bell'
              nmacro = 0
            end if
c
            close (iunit)
c
          end if
c
          if (lnotes) then
            gktext = (restyp(i)//' '//resnam(i)//' ['//
     +                ssenam(sstype(i))//']')
            call pretty (gktext)
            if (baddy) then
              write (junit,*) 
              write (junit,*) 'OOPS - ',gktext(1:leng1(gktext))
            else
              write (junit,*) 
              write (junit,*) 'OKAY - ',gktext(1:leng1(gktext))
            end if
          end if
c
          if (lhtml) then
            gktext = (restyp(i)//' '//resnam(i)//' ['//
     +                ssenam(sstype(i))//']')
            call pretty (gktext)
            myline = '<tr><td>' // gktext(1:leng1(gktext)) //
     +               '</td><td>'
            call putlin (kunit,myline,.false.,ierr)
          end if
c
          line = filnam
c
          call xopxua (iunit,filnam,.false.,ierr)
          if (ierr .ne. 0) goto 998
c
c ... centre on this residue's CA if it exists, else on its centre-of-mass
c
          if (mcptr(ca,i) .gt. 0) then
            write (iunit,*) '$oops_c_atom ',molnam(1:leng1(molnam)),
     +        ' ',resnam(i)(1:leng1(resnam(i))),' CA '
          else
            write (iunit,*) '$oops_c_zone ',molnam(1:leng1(molnam)),
     +        ' ',resnam(i)(1:leng1(resnam(i))),' ;'
          end if
c
c ... 021108 - do user commands straight after centering
c
          write (iunit,*) oline(1:leng1(oline))
c
          write (iunit,*) 'print ..... '
          write (iunit,*) 'print Residue ',
     +      restyp(i)(1:leng1(restyp(i))),' ',
     +      resnam(i)(1:leng1(resnam(i))),' [',
     +      ssenam(sstype(i))(1:leng1(ssenam(sstype(i)))),'] [',
     +      typnam(moltyp)(1:leng1(typnam(moltyp))),']'
          if (baddy) then
            write (iunit,*) 'message OOPS - Residue ',
     +        restyp(i)(1:leng1(restyp(i))),' ',
     +        resnam(i)(1:leng1(resnam(i))),' [',
     +        ssenam(sstype(i))(1:length(ssenam(sstype(i)))),'] [',
     +        typnam(moltyp)(1:leng1(typnam(moltyp))),']'
          else
            write (iunit,*) 'message OKAY - Residue ',
     +        restyp(i)(1:leng1(restyp(i))),' ',
     +        resnam(i)(1:leng1(resnam(i))),' [',
     +        ssenam(sstype(i))(1:length(ssenam(sstype(i)))),'] [',
     +        typnam(moltyp)(1:leng1(typnam(moltyp))),']'
          end if
c
          if (lsubj) then
            write (iunit,'(a,i6)')
     +        ' symbol oops_irc ',i
          end if
c
          if (badpep) then
            npep = npep + 1
            write (gktext,'(a,f10.2)')
     +        'print Bad pep-flip = ',pepflip(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badrsf) then
            nrsf = nrsf + 1
            write (gktext,'(a,f10.3)')
     +        'print Bad RS-fit (all atoms) = ',rsfit(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem)
     +        write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badrsm) then
            nrsm = nrsm + 1
            write (gktext,'(a,f10.3)')
     +        'print Bad RS-fit (main chain) = ',rsfitm(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem)
     +        write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badrss) then
            nrss = nrss + 1
            write (gktext,'(a,f10.3)')
     +        'print Bad RS-fit (side chain) = ',rsfits(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem)
     +        write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badrsr) then
            nrsr = nrsr + 1
            write (gktext,'(a,f10.3)')
     +        'print Bad RS R-factor (all atoms) = ',rsrfac(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem)
     +        write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badrsc) then
            nrsc = nrsc + 1
            write (gktext,'(a,f10.2)')
     +        'print Bad RSC = ',rsc(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badram) then
            nram = nram + 1
            write (gktext,'(a,f10.2,1x,f10.2)')
     +        'print Bad Phi-Psi = ',phi(i),psi(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badpro) then
            npro = npro + 1
            write (gktext,'(a,f10.2,1x,f10.2)')
     +        'print NOTE - Pro PHI (not near -65.4 +/- 11.2) = ',phi(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badpos) then
            npos = npos + 1
            write (gktext,'(a,f10.2)')
     +        'print NOTE - non-Gly positive PHI = ',phi(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badlef) then
            nlef = nlef + 1
            write (gktext,'(a,a,f10.2,1x,f10.2)')
     +        'print NOTE - left-handed helical fragment; ',
     +        'Phi-Psi = ',phi(i),psi(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badmas) then
            nmas = nmas + 1
            write (iunit,*)
     +        'print Mask too tight'
            if (lnotes) write (junit,*)
     +        'Mask too tight'
            if (lhtml) write (kunit,'(1x,a,a)')
     +        'Mask too tight','<br>'
            if (listem) write (*,*) 'Mask too tight'
          end if
c
          if (badblo) then
            nblo = nblo + 1
            write (gktext,'(a,f10.2)')
     +        'print Too low temperature factor = ',lob(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badbhi) then
            nbhi = nbhi + 1
            write (gktext,'(a,f10.2)')
     +        'print Too high temperature factor = ',hib(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (baddbb) then
            ndbb = ndbb + 1
            write (gktext,'(a,f10.2)')
     +        'print Too high RMS delta-B bonded = ',rmsdbb(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badqlo) then
            nqlo = nqlo + 1
            write (gktext,'(a,f10.2)')
     +        'print Too low occupancy = ',loq(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badqhi) then
            nqhi = nqhi + 1
            write (gktext,'(a,f10.2)')
     +        'print Too high occupancy = ',hiq(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badpla) then
            npla = npla + 1
            write (gktext,'(a,f10.2)')
     +        'print Non-planar peptide; improper = ',plangl(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badome) then
            nome = nome + 1
            write (gktext,'(a,f10.2)')
     +        'print Unusual omega value; omega = ',omega(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badcis) then
            ncis = ncis + 1
            write (gktext,'(a,f10.2)')
     +        'print NOTE - cis-peptide; omega = ',omega(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badzet) then
            nzet = nzet + 1
            write (gktext,'(a,f10.2)')
     +        'print Bad CA chirality; zeta = ',zeta(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (baddaa) then
            ndaa = ndaa + 1
            write (gktext,'(a,f10.2)')
     +        'print NOTE - D-amino acid; zeta = ',zeta(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badisu) then
            if (sstors(idisu) .gt. 0.0) then
              write (gktext,'(5a,f6.2,a,f6.1)')
     +        'print NOTE - Right-handed disulfide bridge Cys ',
     +        resnam(ssresi(1,idisu)),' - ',resnam(ssresi(2,idisu)),
     +        '; S=S = ',ssdist(idisu),'; CB-S-S-CB = ',sstors(idisu)
            else
              write (gktext,'(5a,f6.2,a,f6.1)')
     +        'print NOTE - Left-handed disulfide bridge Cys ',
     +        resnam(ssresi(1,idisu)),' - ',resnam(ssresi(2,idisu)),
     +        '; S=S = ',ssdist(idisu),'; CB-S-S-CB = ',sstors(idisu)
            end if
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badssd) then
            nssd = nssd + 1
            write (gktext,'(a,f10.2)')
     +        'print Bad S=S disulfide distance = ',ssdist(idisu)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badsst) then
            nsst = nsst + 1
            write (gktext,'(a,f10.2)')
     +        'print Bad CB-S-S-CB disulfide torsion = ',sstors(idisu)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (insert(i)) then
            nins = nins + 1
            write (gktext,'(a)')
     +        'print Newly inserted residue in this model'
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (mutate(i)) then
            nmut = nmut + 1
            write (gktext,'(a)')
     +        'print Newly mutated residue in this model'
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (baddis) then
            ndis = ndis + 1
            write (gktext,'(a,f10.2)')
     +        'print Large positional shift; RMSD = ',xxrmsd(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badtem) then
            ndis = ndis + 1
            write (gktext,'(a,f10.2)')
     +        'print Large temp.-factor shift; RMSD B = ',xxrmsb(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badocc) then
            nocc = nocc + 1
            write (gktext,'(a,f10.2)')
     +        'print Large occupancy shift; RMSD Q = ',xxrmsq(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badphi) then
            nphi = nphi + 1
            write (gktext,'(a,f10.2)')
     +        'print Large Phi/Psi shift; RMSD = ',xxphps(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badchi) then
            nchi = nchi + 1
            write (gktext,'(a,f10.2)')
     +        'print Large Chi1/2 shift; RMSD = ',xxchid(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badwat) then
            nwat = nwat + 1
            write (gktext,'(a,i10)')
     +        'print Bad QualWat; value = ',watqua(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badcon) then
            nbco = nbco + 1
            write (gktext,'(a,i10)')
     +        'print Bad contact(s); count = ',nbadco(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
          if (badwif) then
            nbwi = nbwi + 1
            do j=1,nwif
              if (wifptr(j) .eq. i) then
                gktext = 'print ' // wifmes(j)
                call pretty (gktext)
                write (iunit,'(1x,a)') gktext(1:leng1(gktext))
                if (lnotes) write (junit,'(1x,a)')
     +            gktext(7:leng1(gktext))
                if (lhtml) write (kunit,'(1x,a,a)')
     +            gktext(7:leng1(gktext)),'<br>'
                if (listem) write (*,*) gktext(7:leng1(gktext))
              end if
            end do
          end if
c
          if (badpor) then
            write (gktext,'(a,i2)')
     +        'print Residue flagged by user; status = ',
     +        oldsub(i)
            call pretty (gktext)
            write (iunit,'(1x,a)') gktext(1:leng1(gktext))
            if (lnotes) write (junit,'(1x,a)') gktext(7:leng1(gktext))
            if (lhtml) write (kunit,'(1x,a,a)')
     +        gktext(7:leng1(gktext)),'<br>'
            if (listem) write (*,*) gktext(7:leng1(gktext))
          end if
c
c ... user-defined criteria
c
          if (nuser .gt. 0) then
            do j=1,nuser
              if (baduse(j)) then
                nrus(j) = nrus(j) + 1
                write (gktext,*) 'print ',
     +            userln(j)(1:leng1(userln(j))),
     +            ' ',usercr(i,j)
                call pretty (gktext)
                write (iunit,'(1x,a)') gktext(1:leng1(gktext))
                if (lnotes) 
     +            write (junit,'(1x,a)') gktext(7:leng1(gktext))
                if (lhtml) write (kunit,'(1x,a,a)')
     +            gktext(7:leng1(gktext)),'<br>'
                if (listem) write (*,*) gktext(7:leng1(gktext))
              end if
            end do
          end if
c
          if (.not. baddy) then
            nokay = nokay + 1
            write (iunit,'(1x,a)')
     +        'print ... Nothing wrong with this residue'
ccc            if (lnotes) write (junit,'(1x,a)')
ccc     +        '... Nothing wrong with this residue'
          end if
c
          write (iunit,*) 'print ..... '
c
          if (lnotes .and. baddy) then
            write (junit,*) 'COMMENTS/ACTION --> '
          end if
c
          if (lhtml) then
            write (kunit,'(a)') '&nbsp;</td></tr>'
          end if
c
        end if
c
      end do
c
      if (lhtml) then
        myline = ' '
        call gkdate (myline)
        myline = '</table><hr noshade>Generated at ' //
     +    myline(1:leng1(myline))
        call putlin (kunit,myline,.false.,ierr)
        myline = 'by ' // prognm // ' version ' // vers
        call putlin (kunit,myline,.false.,ierr)
        myline = '</body></html>'
        call putlin (kunit,myline,.false.,ierr)
        close (kunit)
      end if
c
      write (*,*)
ccc      if (lalres) then
      call jvalut (
     +  ' SUMMARY - Nr of macros generated :',1,nbad)
      call jvalut (
     +  ' SUMMARY - Nr of baddies          :',1,(nbad-nokay))
ccc      else
ccc        call jvalut (' Nr of baddies :',1,nbad)
ccc      end if
c
 6062 format (
     +  ' SUMMARY - Generated macros for ',i6,' out of ',i6,1x,a)
      write (*,6062) cntbad(1),cnttyp(1),'amino acids'
      write (*,6062) cntbad(2),cnttyp(2),'nucleotides'
      write (*,6062) cntbad(3),cnttyp(3),'water molecules'
      write (*,6062) cntbad(4),cnttyp(4),'heterogens'
c
      if (nbad .le. 0) then
        call errcon ('No macros generated - copping out')
        goto 9900
      end if
c
      if (.not. lworst) then
ccc      if (lalres .and. (.not. lworst)) then
        write (iunit,*) 'print No more residues'
ccc      else
ccc        write (iunit,*) 'print No more baddies'
      end if
c
      if (nbad .gt. 0 .and. lchain .and. (.not. lworst)) then
        write (iunit,*)
     +    'menu @',line(1:leng1(line)),' off'
c     +    'menu @',line(1:leng1(line)),' off on_off'
ccc        write (iunit,*) oline(1:leng1(oline))
      end if
c
      if (lchain .and. (.not. lworst)) then
ccc        if (lalres) then
          write (iunit,*) 'bell message No more residues'
ccc        else
ccc          write (iunit,*) 'bell message No more baddies'
ccc        end if
        write (iunit,*) 'bell print ...'
        write (iunit,*)
     +    'print Type the following to clean the OOPS directory:'
        write (iunit,*)
     +    'print $ rm oops/*'
        write (iunit,*) 'bell print ...'
        if (lsubj) then
          write (iunit,*)
     +      'print Subjective residue quality is stored in'
          write (iunit,*)
     +      'print datablock ',dbsub
          write (iunit,*)
     +      'print and file oops_subjectivity.table'
          write (iunit,*) 'bell print ...'
          write (iunit,*) 'pa_prop atom_z > 0 green'
          write (iunit,*) 'pa_prop residue_quality = -1 blue'
          write (iunit,*) 'pa_prop residue_quality = 1 red'
          write (iunit,*) 'pa_prop residue_quality = 2 magenta'
          write (iunit,*) 'pa_prop residue_quality = 3 white'
          write (iunit,*) 'obj qual ca ; end'
          write (iunit,*) 
     +  'pa_case atom_z 5 6 7 8 15 16 yellow blue red orange green'
          write (iunit,*) 'menu @Good_fit off'
          write (iunit,*) 'menu @Poor_fit off'
          write (iunit,*) 'menu @Poor_dens off'
          write (iunit,*) 'menu @No_dens off'
          write (iunit,*) 'write ',rtypdb,' qqq1 (a6)'
          write (iunit,*) 'write ',rnamdb,' qqq2 (a6)'
          write (iunit,*) 'write ',dbsub,' qqq3 (i6)'
          write (iunit,*)
     + '$ paste qqq1 qqq2 qqq3 | grep -v RESIDUE > qqq4'
          write (iunit,*)
     + '$ cp qqq4 oops_subjectivity.table'
ccc     + '$ cat qqq4 | grep -v "\-1" > oops_subjectivity.table'
          write (iunit,*) '$ rm qqq1 qqq2 qqq3 qqq4'
          write (iunit,*) 'bell print ... All done !'
          write (iunit,*) 'bell message ... All done !'
        end if
      end if
c
      close (iunit)
      if (lnotes) then
        write (junit,*)
        close (junit)
      end if
c
      if (nbad .gt. 0) then
        write (*,*)
        write (*,*) 'Start by typing @oops.omac in O !!!'
      end if
c
      call prompt ('0SUMMARY - Crystallographic indicators:')
      if (lrsfit) call jvalut (
     +  ' SUMMARY - Bad RS-fit (all atoms)       :',1,nrsf)
      if (lrsm) call jvalut (
     +  ' SUMMARY - Bad RS-fit (main chain)      :',1,nrsm)
      if (lrss) call jvalut (
     +  ' SUMMARY - Bad RS-fit (side chain)      :',1,nrss)
      if (lrsr) call jvalut (
     +  ' SUMMARY - Bad RS R-factor (all atoms)  :',1,nrsr)
      if (lmask) call jvalut (
     +  ' SUMMARY - Bad mask                     :',1,nmas)
      if (llob) call jvalut (
     +  ' SUMMARY - Bad B (low)                  :',1,nblo)
      if (lhib) call jvalut (
     +  ' SUMMARY - Bad B (high)                 :',1,nbhi)
      if (ldbb) call jvalut (
     +  ' SUMMARY - Bad RMS delta-B bonded atoms :',1,ndbb)
      if (lloq) call jvalut (
     +  ' SUMMARY - Bad Q (low)                  :',1,nqlo)
      if (lhiq) call jvalut (
     +  ' SUMMARY - Bad Q (high)                 :',1,nqhi)
      if (lwat) call jvalut (
     +  ' SUMMARY - Bad QualWat                  :',1,nwat)
c
      call prompt ('0SUMMARY - Geometric indicators:')
      if (lpep) call jvalut (
     +  ' SUMMARY - Bad pep-flip                 :',1,npep)
      if (lrsc) call jvalut (
     +  ' SUMMARY - Bad RSC                      :',1,nrsc)
      if (lrama) call jvalut (
     +  ' SUMMARY - Bad Phi/Psi                  :',1,nram)
      if (lpro) call jvalut (
     +  ' SUMMARY - Proline Phi not near -65.4   :',1,npro)
      if (lpos) call jvalut (
     +  ' SUMMARY - Positive Phi for non-Gly     :',1,npos)
      if (llefth) call jvalut (
     +  ' SUMMARY - Left-handed helical residues :',1,nlef)
      if (lplan) call jvalut (
     +  ' SUMMARY - Bad peptide planarity        :',1,npla)
      if (lplan) call jvalut (
     +  ' SUMMARY - Bad omega torsion            :',1,nome)
      if (lplan) call jvalut (
     +  ' SUMMARY - Cis-peptides                 :',1,ncis)
      if (lchir) call jvalut (
     +  ' SUMMARY - Bad C-alpha chirality        :',1,nzet)
      if (lchir) call jvalut (
     +  ' SUMMARY - D-amino acids                :',1,ndaa)
      if (ldisu) call jvalut (
     +  ' SUMMARY - Disulfide bridges            :',1,ndisu)
      if (ldisu) call jvalut (
     +  ' SUMMARY - Bad disulfide S=S distance   :',1,nssd)
      if (ldisu) call jvalut (
     +  ' SUMMARY - Bad disulfide CB-S-S-CB      :',1,nsst)
      if (lbadco) call jvalut (
     +  ' SUMMARY - Bad contact(s)               :',1,nbco)
      if (lwif) call jvalut (
     +  ' SUMMARY - Flagged by WHAT IF           :',1,nbwi)
c
      if (lprev) then
        call prompt ('0SUMMARY - Comparison to previous model:')
        call jvalut (
     +  ' SUMMARY - Newly inserted residues      :',1,nins)
        call jvalut (
     +  ' SUMMARY - Newly mutated residues       :',1,nmut)
        call jvalut (
     +  ' SUMMARY - Large positional shift       :',1,ndis)
        call jvalut (
     +  ' SUMMARY - Large temp.-factor shift     :',1,ntem)
        call jvalut (
     +  ' SUMMARY - Large occupancy shift        :',1,nocc)
        call jvalut (
     +  ' SUMMARY - Large Phi/Psi shift          :',1,nphi)
        call jvalut (
     +  ' SUMMARY - Large Chi1/2 shift           :',1,nchi)
      end if
c
      if (nuser .gt. 0) then
        call prompt ('0SUMMARY - User-defined criteria:')
        do i=1,nuser
          gktext = ' SUMMARY - Bad ' //
     +             usernm(i)(1:leng1(usernm(i))) // ' :'
          call jvalut (gktext,1,nrus(i))
        end do
      end if
c
      if (lxgood .or. lipoor) then
        write (*,*)
        call jvalut (' Skipped (good dens, good fit):',1,nxgood)
        call jvalut (' Flagged (poor dens or fit)   :',1,nipoor)
      end if
c
ccc      if (lalres) then
        write (*,*)
        call jvalut (
     +  ' SUMMARY - Residues without problems    :',1,nokay)
ccc      end if
      write (*,*)
c
c ... write O datablock with badcounts
c
      filnam = 'oops_badcounts.o'
      call xopxua (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening badcounts file')
        goto 9002
      end if
c
      call prompt (' Writing oops_badcounts.o ...')
      line = molnam // '_residue_badcounts'
      call remspa (line)
      write (iunit,'(a,1x,a,1x,i6,1x,a)')
     +  line(1:leng1(line)),'I',nres,'(25i3)'
      write (iunit,'(25i3)') (count(i),i=1,nres)
      close (iunit)
c
 9002 continue
c
      if (nbad .gt. 0) then
c
        write (*,*) 'Read oops_badcounts.o into O and use the'
        write (*,*) 'residue property badcounts to colour your'
        write (*,*) 'molecule, or plot this file, to reveal'
        write (*,*) 'areas where the structure may be poor.'
c
        j=nres
        do i=1,nres
          xdum (i) = float(count(i))
        end do
        nhisto = maxhis
        do i=1,nhisto
          histo(i) = float(i-1)
        end do
        call oopsts ('Bad counts',
     +               j,xdum,nhisto,histo,nhcnt,
     +               ave,sdv,xmin,xmax,xtot)
c
        filnam = molnam(1:leng1(molnam))//'_badcounts.plt'
        title = 'Bad counts'
        line = molnam // '_residue_badcounts'
        call remspa (line)
        if (lplot) call ooplot (iunit,filnam,line,
     +    title,nres,xdum,0.0,1.05*xmax,j,ave,sdv,xmin,xmax)
c
      end if
c
c ... write PDB REMARK include file
c
      filnam = 'oops_remarks.pdb'
      call xopxua (iunit,filnam,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB REMARK file')
        goto 9000
      end if
c
 6100 format ('REMARK   5 ',10(a,1x))
 6110 format ('REMARK   5 ',a,i8,a)
 6120 format ('REMARK   5 ',a,f8.2,a)
c
      write (*,*)
      call prompt (' Writing oops_remarks.pdb ...')
      line = ' '
      call gkdate (line(1:24))
      call upcase (line)
      write (iunit,6100)
      write (iunit,6100) 'CREATED BY :',prognm
      write (iunit,6100) 'VERSION    :',vers
      write (iunit,6100) 'MODEL NAME :',molnam
      write (iunit,6100) 'DATE       :',line(1:24)
      write (iunit,6100)
      write (iunit,6100) 'MODEL QUALITY'
c
      if (iflip .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'PEPTIDE-OXYGEN ORIENTATION'
        write (iunit,6100)
     +    '  PROGRAM USED: O'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,iflip))
        write (iunit,6120)
     +    '  MAXIMUM VALUE         (A) :',stats(5,iflip)
        write (iunit,6120)
     +    '  AVERAGE VALUE         (A) :',stats(2,iflip)
        write (iunit,6120)
     +    '  CUT-OFF VALUE         (A) :',stats(1,iflip)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',npep
        perc = 100.0*float(npep)/stats(6,iflip)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (irsc .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'ROTAMER SIDECHAIN ANALYSIS'
        write (iunit,6100)
     +    '  PROGRAM USED: O'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,irsc))
        write (iunit,6120)
     +    '  MAXIMUM VALUE         (A) :',stats(5,irsc)
        write (iunit,6120)
     +    '  AVERAGE VALUE         (A) :',stats(2,irsc)
        write (iunit,6120)
     +    '  CUT-OFF VALUE         (A) :',stats(1,irsc)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',nrsc
        perc = 100.0*float(nrsc)/stats(6,irsc)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (irsfal .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'REAL-SPACE FIT (CORR. COEFF./ALL ATOMS)'
        write (iunit,6100)
     +    '  PROGRAM USED: O'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,irsfal))
        write (iunit,6120)
     +    '  MINIMUM VALUE             :',stats(4,irsfal)
        write (iunit,6120)
     +    '  AVERAGE VALUE             :',stats(2,irsfal)
        write (iunit,6120)
     +    '  STANDARD DEVIATION        :',stats(3,irsfal)
        write (iunit,6120)
     +    '  CUT-OFF VALUE             :',stats(1,irsfal)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',nrsf
        perc = 100.0*float(nrsf)/stats(6,irsfal)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (irsral .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'REAL-SPACE FIT (R-FACTOR/ALL ATOMS)'
        write (iunit,6100)
     +    '  PROGRAM USED: O'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,irsral))
        write (iunit,6120)
     +    '  MAXIMUM VALUE             :',stats(5,irsral)
        write (iunit,6120)
     +    '  AVERAGE VALUE             :',stats(2,irsral)
        write (iunit,6120)
     +    '  STANDARD DEVIATION        :',stats(3,irsral)
        write (iunit,6120)
     +    '  CUT-OFF VALUE             :',stats(1,irsral)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',nrsr
        perc = 100.0*float(nrsr)/stats(6,irsral)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (irsfmc .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'REAL-SPACE FIT (CORR. COEFF./MAINCHAIN)'
        write (iunit,6100)
     +    '  PROGRAM USED: O'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,irsfmc))
        write (iunit,6120)
     +    '  MINIMUM VALUE             :',stats(4,irsfmc)
        write (iunit,6120)
     +    '  AVERAGE VALUE             :',stats(2,irsfmc)
        write (iunit,6120)
     +    '  STANDARD DEVIATION        :',stats(3,irsfmc)
        write (iunit,6120)
     +    '  CUT-OFF VALUE             :',stats(1,irsfmc)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',nrsm
        perc = 100.0*float(nrsm)/stats(6,irsfmc)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (irsfsc .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'REAL-SPACE FIT (CORR. COEFF./SIDECHAIN)'
        write (iunit,6100)
     +    '  PROGRAM USED: O'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,irsfsc))
        write (iunit,6120)
     +    '  MINIMUM VALUE             :',stats(4,irsfsc)
        write (iunit,6120)
     +    '  AVERAGE VALUE             :',stats(2,irsfsc)
        write (iunit,6120)
     +    '  STANDARD DEVIATION        :',stats(3,irsfsc)
        write (iunit,6120)
     +    '  CUT-OFF VALUE             :',stats(1,irsfsc)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',nrss
        perc = 100.0*float(nrss)/stats(6,irsfsc)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (iplanp .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'PEPTIDE PLANARITY IMPROPER TORSION'
        write (iunit,6100)
     +    '  PROGRAM USED: OOPS'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,iplanp))
        write (iunit,6120)
     +    '  MINIMUM VALUE       (DEG) :',stats(4,iplanp)
        write (iunit,6120)
     +    '  MAXIMUM VALUE       (DEG) :',stats(5,iplanp)
        write (iunit,6120)
     +    '  AVERAGE VALUE       (DEG) :',stats(2,iplanp)
        write (iunit,6120)
     +    '  STANDARD DEVIATION  (DEG) :',stats(3,iplanp)
        write (iunit,6120)
     +    '  CUT-OFF VALUE       (DEG) :',stats(1,iplanp)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',npla
        perc = 100.0*float(npla)/stats(6,iplanp)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      if (ichirl .gt. 0) then
        write (iunit,6100)
        write (iunit,6100) 'CA CHIRALITY IMPROPER TORSION'
        write (iunit,6100)
     +    '  PROGRAM USED: OOPS'
        write (iunit,6110)
     +    '  NR OF RESIDUES            :',nint(stats(6,ichirl))
        write (iunit,6120)
     +    '  MINIMUM VALUE       (DEG) :',stats(4,ichirl)
        write (iunit,6120)
     +    '  MAXIMUM VALUE       (DEG) :',stats(5,ichirl)
        write (iunit,6120)
     +    '  AVERAGE VALUE       (DEG) :',stats(2,ichirl)
        write (iunit,6120)
     +    '  STANDARD DEVIATION  (DEG) :',stats(3,ichirl)
        write (iunit,6120)
     +    '  CUT-OFF VALUE       (DEG) :',stats(1,ichirl)
        write (iunit,6110)
     +    '  OUTLIERS             (NR) :',nzet
        perc = 100.0*float(nzet)/stats(6,ichirl)
        write (iunit,6120)
     +    '  OUTLIERS              (%) :',perc
      end if
c
      write (iunit,6100)
c
      close (iunit)
c
c --- for "worst-first", sort out everything
c
 9000 continue
      if (.not. lworst) goto 9900
c
      write (*,*)
      call prompt (' Sorting "worst-first" ...')
      call flusho (6)
c
      j = bcount(1)
      do i=1,nbad
        if (bcount(i) .gt. j) j = bcount(i)
      end do
c
      k = 0
      do l=j,0,-1
        do i=1,nbad
          if (bcount(i) .eq. l) then
            k = k + 1
            bindex (k) = i
            if (k .eq. nbad) goto 9062
          end if
        end do
      end do
 9062 continue
c
c      call ivalut (' COUNTS BEFORE :',10,bcount)
c      call ivalut (' INDEX  BEFORE :',10,bindex)
c      call shell (bcount,bindex,nbad)
c      do i=1,nbad/2
c        j = nbad - i + 1
c        call iswap (bindex(i),bindex(j))
c        call iswap (bcount(i),bcount(j))
c      end do
c      call ivalut (' COUNTS AFTER  :',10,bcount)
c      call ivalut (' INDEX  AFTER  :',10,bindex)
c
 6206 format (' # ',i6,' = baddy # ',i6,' with ',i3,
     +  ' quality warnings in file : ',a)
c
      nmacro = 0
      do i=1,nbad
        nmacro = nmacro + 1
        j = bindex (i)
        k = bcount (j)
        filnam = bfile (j)
c
        write (*,6206) i,j,k,filnam(1:leng1(filnam))
c
ccc        print *,i,j,k,filnam
c
        if (i .eq. 1) then
          xlines(1) = ' @'//filnam(1:leng1(filnam))
          nl = 1
          call appfil (iunit,'oops.omac',nl,xlines)
        end if
c
        if (i .lt. nbad) then
          filnext = bfile ( bindex (i+1) )
          xlines (1) = ' print Hit or type "@' //
     +      filnext(1:leng1(filnext)) // '" for next residue'
          xlines (2) = ' menu @' // filnext(1:leng1(filnext))
     +      // ' on'
          nl = 2
          if (i .gt. 1) then
            nl = nl + 1
ccc            filnext = bfile ( bindex (i-1) )
            xlines (nl) = ' menu @' // filnam(1:leng1(filnam))
     +        // ' off'
          end if
          if (nmacro .eq. 50) then
            nl = nl + 1
            xlines (nl) = 
     +   ' print You checked another 50 residues ... take a break !'
            nl = nl + 1
            xlines (nl) = 
     +   ' message You checked another 50 residues ... take a break !'
            nmacro = 0
          end if
          nl = nl + 1
          write (xlines(nl),*) ' print ',i,' down, ',(nbad-i),' to go'
          call pretty (xlines(nl)(2:))
          call appfil (iunit,filnam,nl,xlines)
        end if
c
        if (i .eq. nbad) then
          xlines (1) = ' print No more residues'
          xlines (2) = ' bell message No more residues'
          xlines (3) = ' bell print ...'
          xlines (4) = 
     + ' print Type the following to clean the OOPS directory:'
          xlines (5) = ' print $ rm oops/*'
          xlines (6) = ' bell print ...'
          nl = 6
          if (lsubj) then
            xlines (7) =
     + ' print Subjective residue quality is stored in'
            xlines (8) = ' print datablock '//dbsub
            xlines (9) =
     + 'print and file oops_subjectivity.table'
            xlines (10) = ' bell print ...'
            xlines (11) = ' pa_prop atom_z > 0 green'
            xlines (12) = ' pa_prop residue_quality = -1 blue'
            xlines (13) = ' pa_prop residue_quality = 1 red'
            xlines (14) = ' pa_prop residue_quality = 2 magenta'
            xlines (15) = ' pa_prop residue_quality = 3 white'
            xlines (16) = ' obj qual ca ; end'
            xlines (17) =
     +  'pa_case atom_z 5 6 7 8 15 16 yellow blue red orange green'
            xlines (18) = ' menu @Good_fit off'
            xlines (19) = ' menu @Poor_fit off'
            xlines (20) = ' menu @Poor_dens off'
            xlines (21) = ' menu @No_dens off'
            xlines (22) = ' write '//rtypdb//' qqq1 (a6)'
            xlines (23) = ' write '//rnamdb//' qqq2 (a6)'
            xlines (24) = ' write '//dbsub//' qqq3 (i6)'
            xlines (25) = 
     + '$ paste qqq1 qqq2 qqq3 | grep -v RESIDUE > qqq4'
            xlines (26) = 
     + '$ cp qqq4 oops_subjectivity.table'
ccc     + '$ cat qqq4 | grep -v "\-1" > oops_subjectivity.table'
            xlines (27) = ' $ rm qqq1 qqq2 qqq3 qqq4'
            xlines (28) = ' bell print ... All done !'
            xlines (29) = ' bell message ... All done !'
            nl = 29
          end if
c
          if (i .gt. 1) then
            nl = nl + 1
            filnext = bfile ( bindex (i-1) )
            xlines (nl) = ' menu @' // filnext(1:leng1(filnext))
     +        // ' off'
          end if
c
          call appfil (iunit,filnam,nl,xlines)
c
        end if
c
      end do
c
 9900 continue
c
c --- END OF OOPS
c
      call gkquit ()
c
  998 continue
      call errstp ('While opening O macro file')
c
  999 continue
      call errstp ('While reading O datablock file')
c
      end
c
c
c
      subroutine newram
c
c ... set up areas which are allowed in PHI-PSI space
c     use our new definition of core regions
c
c ... Gerard Kleywegt @ 960801
c
      include 'oops2.incl'
c
      real pout
c
      integer coregn (37,37)
c
      integer i,ix,jx,nbad,ndef
c
code ...
c
      call defcor (coregn)
c
      ndef = 0
      nbad = 0
      badpp (1) = .false.
      badpp (nres) = .false.
c
      do i=2,nres-1
        if (restyp(i) .eq. 'GLY') then
          badpp (i) = .false.
        else if (phi(i) .ge. 999. .or.
     +           psi(i) .ge. 999.) then
          badpp (i) = .false.
        else
          ix = 1 + int( (180.0 + phi(i)) / 10.0)
          jx = 1 + int( (180.0 + psi(i)) / 10.0)
c
ccc          print *,i,phi(i),psi(i),ix,jx
c
          ndef = ndef + 1
          badpp (i) = (coregn(ix,jx) .ne. 1)
        end if
        if (badpp(i)) nbad = nbad + 1
      end do
c
      if (ndef .gt. 0) then
        pout = 100.0*float(nbad)/float(ndef)
      else
        pout = 0.0
      end if
c
      call jvalut (' Nr of residues considered :',1,ndef)
      call jvalut (' Nr of outliers            :',1,nbad)
      call fvalut (' Percentage outliers       :',1,pout)
c
      return
      end
c
c
c
      subroutine iniram
c
c ... set up areas which are allowed in PHI-PSI space
c
c ... Gerard Kleywegt @ 930305
c
      include 'oops2.incl'
c
      integer nrok
      parameter (nrok = 12)
c
      real area(6,nrok)
      real x1,x2,y1l,y1h,y2l,y2h,ymin,ymax
c
      integer i,ix,j,ixlo,ixhi
c
c ... define the 12 allowed areas
c
      data area/-176.0,-67.1,-20.0,  -147.9,-67.1,-35.3,
     + -147.9,-67.1,-35.3,  -129.5,-67.1,-31.5,
     + -129.5,-67.1,-31.5,  -122.1,-67.1,-19.5,
     + -122.1,-67.1,180.0,  -86.0,-67.1,180.0,
     + -176.0,40.4,180.0,   -122.1,-19.5,180.0,
     + -86.0,66.2,180.0,    -53.0,84.6,180.0,
     + -53.0,84.6,180.0,    -39.8,96.9,151.6,
     + -86.0,-67.1,38.1,    -57.0,-67.1,-1.1,
     + -57.0,-67.1,-1.1,    -35.5,-64.8,-18.3,
     + 34.7,27.5,87.4,      58.2,11.2,110.4,
     + -175.4,-180.0,-169.5,-61.6,-180.0,-169.5,
     + -61.6,-180.0,-169.5, -57.0,-180.0,-180.0/
c
code ...
c
      do i=-180,180
        nramok (i) = 0
      end do
c
      do i=1,nrok
        x1  = area(1,i)
        y1l = area(2,i)
        y1h = area(3,i)
        x2  = area(4,i)
        y2l = area(5,i)
        y2h = area(6,i)
        ixlo = max (-180,int(x1))
        ixhi = min (180,1+int(x2))
        do ix=ixlo,ixhi
          ymin = y1l + (float(ix)-x1)*(y2l-y1l)/(x2-x1)
          ymax = y1h + (float(ix)-x1)*(y2h-y1h)/(x2-x1)
          if (nramok(ix) .gt. 0) then
            do j=1,nramok(ix)
              if ((ymin .ge. ramaok(1,j,ix) .and.
     +             ymin .le. ramaok(2,j,ix)) .or.
     +            (ymax .ge. ramaok(1,j,ix) .and.
     +             ymax .le. ramaok(2,j,ix))) then
                ramaok(1,j,ix) = min (ymin,ramaok(1,j,ix))
                ramaok(2,j,ix) = max (ymax,ramaok(2,j,ix))
                goto 6900
              end if
            end do
          end if
c
          nramok(ix) = nramok(ix) + 1
          j = nramok(ix)
          ramaok (1,j,ix) = ymin
          ramaok (2,j,ix) = ymax
c
 6900     continue
c
        end do
      end do
c
      do i=1,nres
        if (restyp(i) .eq. 'GLY') then
          badpp (i) = .false.
        else if (phi(i) .ge. 999. .or.
     +           psi(i) .ge. 999.) then
          badpp (i) = .false.
        else
          ix = nint(phi(i))
          badpp (i) = .true.
          if (nramok(ix) .gt. 0) then
            do j=1,nramok(ix)
              if (psi(i) .ge. ramaok(1,j,ix) .and.
     +            psi(i) .le. ramaok(2,j,ix)) then
                badpp (i) = .false.
                goto 6800
              end if
            end do
          end if
 6800     continue
        end if
      end do
      badpp (1) = .false.
      badpp (nres) = .false.
c
      return
      end
c
c
c
      subroutine getpdb (iunit,ierr)
c
      include 'oops2.incl'
c
c      integer nwatnm
c      parameter (nwatnm = 7)
c
      integer iunit,ierr,i,j,ires,nhet,nhydro
c
      logical hetatm,lhydro
c
      character line*80,nowres*6
c
code ...
c
      ierr = -1
c
      nres = 0
      nowres = '??????'
c
      natoms = 0
      nhet = 0
      nhydro = 0
c
   10 continue
      read (iunit,'(a)',end=100,err=9999) line
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') then
        if (line(1:6) .ne. 'ANISOU') then
          line = ' '//line
          call textut (line(1:7),line(8:))
        end if
        goto 10
      end if
c
      hetatm = (line(1:6) .eq. 'HETATM')
c
      i = natoms + 1
      if (i .gt. maxatm) then
        call errcon ('Too many atoms; rest skipped')
        goto 100
      end if
c
      call upcase (line)
      read (line,6000,err=9999) atmnam(i),atmtyp(i),atnmnm(i),
     +  (pdbxyz(j,i),j=1,3),qocc(i),bfac(i)
c
      if (lhydro(atmnam(i))) then
        nhydro = nhydro + 1
        goto 10
      end if
c
 6000 format (12x,a4,1x,a3,1x,a6,3x,3f8.3,f6.2,f6.2)
c
      call remspa (atnmnm(i))
c
c ... O puts a "$" in front of HETATM residue names
c
      if (hetatm) then
        atnmnm(i) = '$'//atnmnm(i)
        nhet = nhet + 1
      end if
c
c ... new residue ???
c
      if (atnmnm(i) .ne. nowres) then
        if (nres .ge. maxres) then
          call errstp ('GETPDB - Too many residues; rest skipped')
          goto 100
        end if
        nres = nres + 1
        resnam (nres) = atnmnm(i)
        call upcase (resnam(nres))
        nowres = atnmnm(i)
        restyp (nres) = atmtyp(i)
        call upcase (restyp(nres))
        call remspa (restyp(nres))
      end if
c
      natoms = i
c
      goto 10
c
 9999 continue
      return
c
  100 continue
      close (iunit)
      call jvalut (' Nr of residues read :',1,nres)
      call jvalut (' Nr of atoms read    :',1,natoms)
      call jvalut (' Nr of HETATM cards  :',1,nhet)
      call jvalut (' Nr of Hs stripped   :',1,nhydro)
c
      if (natoms .lt. 10) then
        call errcon ('Too few atoms -- be serious !')
        return
      end if
c
      if (nres .lt. 3) then
        call errcon ('Too few residues -- be serious !')
        return
      end if
c
c ... generate pointer arrays (for each residue, pointer to
c     first and last atom; separately, pointers to N, CA, C,
c     O and CB)
c
      do i=1,maxres
        resptr (1,i) = 0
        resptr (2,i) = 0
        do j=1,nmc
          mcptr (j,i) = 0
        end do
      end do
c
      do i=1,natoms
        nowres = atnmnm(i)
c
c ... recognise it ?
c
        do j=1,nres
          if (resnam(j) .eq. nowres) then
            ires = j
            goto 6900
          end if
        end do
c
        call textut (' Unrecognised residue name :',nowres)
        ires = 0
c
 6900   continue
        if (ires .gt. 0) then
          if (resptr(1,ires) .eq. 0) then
            resptr(1,ires) = i
            resptr(2,ires) = i
          else
            if (i .eq. (resptr(2,ires)+1) ) then
              resptr (2,ires) = i
            else
              call errcon (
     +  'Atoms of residue scattered throughout the PDB file')
              call textut (' Residue :',nowres)
              return
            end if
          end if
c
c ... store pointers to main-chain atoms
c
          do j=1,nmc
            if (atmnam(i) .eq. mcname(j)) then
              mcptr(j,ires) = i
            end if
          end do
c
c ... special case for XPLOR OT1, OT2
c
          if (mcptr(o,ires) .lt. 1) then
            if (atmnam(i) .eq. ' OT1' .or.
     +          atmnam(i) .eq. ' OT2') mcptr(o,ires)=i
          end if
c
c ... special case for PROLSQ OT
c
          if (mcptr(o,ires) .lt. 1) then
            if (atmnam(i) .eq. ' OT ') mcptr(o,ires)=i
          end if
c
c ... special case for OTX
c
          if (mcptr(o,ires) .lt. 1) then
            if (atmnam(i) .eq. ' OTX') mcptr(o,ires)=i
          end if
c
c ... special case for PDB OXT
c
          if (mcptr(o,ires) .lt. 1) then
            if (atmnam(i) .eq. ' OXT') mcptr(o,ires)=i
          end if
c
        end if
c
      end do
c
      do i=1,nres
        if (resptr(1,i) .lt. 1) then
          call errcon ('Residue has no atoms')
          call textut (' Residue :',resnam(i))
          return
        end if
      end do
c
      do i=1,nmc
        ires = 0
        do j=1,nres
          if (mcptr(i,j) .gt. 0) ires=ires+1
        end do
        call jvalut (' Nr of residues with '//mcname(i)//' :',1,ires)
      end do
      call jvalut (' Total nr of residues     :',1,nres)
c
ccc      call fvalin (' Max CA-CA distance for neighbours ?',
ccc     +  1,mxcaca)
ccc      mxcaca = max (mxcaca,4.0)
c
      call fvalut (' Max CA-CA distance for neighbours :',
     +  1,mxcaca)
c
      call yasspa ()
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine tightm
c
      include 'oops2.incl'
c
      integer i,ires,istat,nbad
c
code ...
c
      nbad = 0
      do ires=1,nres
        badmsk (ires) = .false.
        do i=resptr(1,ires),resptr(2,ires)
          if (.not. badmsk(ires)) then
            call lmokay (mask,extent(1),extent(2),extent(3),
     +                   origin,grid,cell,pdbxyz(1,i),mskrad,istat)
            if (istat .ne. 0) badmsk(ires) = .true.
          end if
        end do
        if (badmsk(ires)) nbad = nbad + 1
      end do
      call jvalut (' Residues with mask problems :',1,nbad)
c
      return
      end
c
c
c
      subroutine lmokay (mask,exta,extb,extc,origin,grid,cell,
     +                   xyz,rad,istat)
c
      real cell(6),xyz(3),rad,b(3,3),g(3),x(3),fake(6),a(3,3)
      real c(3,3),x1(3),xp(3),x2(3),distce
c
      integer exta,extb,extc,origin(3),grid(3),istat,off
      integer i,j,k,l,i1,i2,i3
      integer mask (exta,extb,extc)
c
code ...
c
      istat = 0
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
c
      do j=1,3
        x(j) = xyz(j)
      end do
      call mulmtx (a, x, x1, 3, 3, 1)
      off = rad/g(1)+1
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
      return
c
c ... below lower grid bounds
c
  199 continue
      istat = -1
      return
c
c ... above upper grid bounds
c
  198 continue
      istat = -2
      return
c
c ... atom + radius outside current mask
c
  197 continue
      istat = -3
      return
c
      end
c
c
c
      subroutine lowbs
c
      include 'oops2.incl'
c
      integer i,ires,nbad
c
code ...
c
      nbad = 0
      do ires=1,nres
        badlob (ires) = .false.
        lob (ires) = 9999.99
        do i=resptr(1,ires),resptr(2,ires)
          lob (ires) = min (bfac(i),lob(ires))
          if (.not. badlob(ires)) then
            badlob (ires) = (bfac(i) .lt. lobfac)
          end if
        end do
        if (badlob(ires)) nbad = nbad + 1
      end do
      call jvalut (' Residues with low Bs :',1,nbad)
c
      return
      end
c
c
c
      subroutine lowqs
c
      include 'oops2.incl'
c
      integer i,ires,nbad
c
code ...
c
      nbad = 0
      do ires=1,nres
        badloq (ires) = .false.
        loq (ires) = 9999.99
        do i=resptr(1,ires),resptr(2,ires)
          loq (ires) = min (qocc(i),loq(ires))
          if (.not. badloq(ires)) then
            badloq (ires) = (qocc(i) .lt. loqocc)
          end if
        end do
        if (badloq(ires)) nbad = nbad + 1
      end do
      call jvalut (' Residues with low Qs :',1,nbad)
c
      return
      end
c
c
c
      subroutine highqs
c
      include 'oops2.incl'
c
      integer i,ires,nbad
c
code ...
c
      nbad = 0
      do ires=1,nres
        badhiq (ires) = .false.
        hiq (ires) = -9999.99
        do i=resptr(1,ires),resptr(2,ires)
          hiq (ires) = max (qocc(i),hiq(ires))
          if (.not. badhiq(ires)) then
            badhiq (ires) = (qocc(i) .gt. hiqocc)
          end if
        end do
        if (badhiq(ires)) nbad = nbad + 1
      end do
      call jvalut (' Residues with high Qs :',1,nbad)
c
      return
      end
c
c
c
      subroutine highbs
c
      include 'oops2.incl'
c
      integer i,ires,nbad
c
code ...
c
      nbad = 0
      do ires=1,nres
        badhib (ires) = .false.
        hib (ires) = -9999.99
        do i=resptr(1,ires),resptr(2,ires)
          hib (ires) = max (bfac(i),hib(ires))
          if (.not. badhib(ires)) then
            if (swater(ires)) then
              badhib (ires) = (bfac(i) .gt. hibwat)
            else
              badhib (ires) = (bfac(i) .gt. hibfac)
            end if
          end if
        end do
        if (badhib(ires)) nbad = nbad + 1
      end do
      call jvalut (' Residues with high Bs :',1,nbad)
c
      return
      end
c
c
c
      subroutine phipsi
c
      include 'oops2.incl'
c
      real tangle,dist
c
      integer ires,icnt,jcnt
c
code ...
c
      icnt = 0
      jcnt = 0
c
      do ires=1,nres
c
        phi (ires) = 999.99
        psi (ires) = 999.99
c
        if (mcptr(ca,ires)   .le. 0) goto 300
c
        if (ires .eq. 1) goto 200
c
        if (mcptr(ca,ires-1) .le. 0) goto 200
        if (dist(mcptr(ca,ires),mcptr(ca,ires-1),
     +      pdbxyz) .gt. mxcaca) goto 200
c
c        print *,ires,phi(ires),psi(ires)
c
        if (mcptr(c,ires-1) .gt. 0 .and.
     +      mcptr(n,ires)   .gt. 0 .and.
     +      mcptr(ca,ires)  .gt. 0 .and.
     +      mcptr(c,ires)   .gt. 0) then
          phi (ires) = tangle (mcptr(c,ires-1),
     +       mcptr(n,ires),mcptr(ca,ires),mcptr(c,ires),
     +       pdbxyz)
c
          call fixang (phi(ires))
c
          icnt = icnt + 1
        end if
c
  200   continue
c
        if (ires .eq. nres) goto 300
c
        if (mcptr(ca,ires+1) .le. 0) goto 300
        if (dist(mcptr(ca,ires),mcptr(ca,ires+1),
     +      pdbxyz) .gt. mxcaca) goto 300
c
        if (mcptr(n,ires)   .gt. 0 .and.
     +      mcptr(ca,ires)  .gt. 0 .and.
     +      mcptr(c,ires)   .gt. 0 .and.
     +      mcptr(n,ires+1) .gt. 0) then
          psi (ires) = tangle (mcptr(n,ires),
     +       mcptr(ca,ires),mcptr(c,ires),mcptr(n,ires+1),
     +       pdbxyz)
c
          call fixang (psi(ires))
c
          jcnt = jcnt + 1
        end if
c
  300   continue
c
c        print *,ires,phi(ires),psi(ires)
c
      end do
c
      phi (1) = 999.99
ccc      psi (1) = 999.99
ccc      phi (nres) = 999.99
      psi (nres) = 999.99
c
      call jvalut (' Nr of residues with defined PHI :',1,icnt)
      call jvalut (' Nr of residues with defined PSI :',1,jcnt)
c
      return
      end
c
c
c
      subroutine plapep
c
      include 'oops2.incl'
c
      real tangle,dist
c
      integer ires,icnt,jcnt
c
code ...
c
      icnt = 0
      jcnt = 0
c
      do ires=1,nres-1
c
        plangl (ires) = 999.99
        omega  (ires) = 999.99
c
        if (mcptr(ca,ires)   .le. 0) goto 200
        if (mcptr(ca,ires+1) .le. 0) goto 200
        if (dist(mcptr(ca,ires),mcptr(ca,ires+1),
     +      pdbxyz) .gt. mxcaca) goto 200
c
        if (mcptr(c,ires)   .gt. 0 .and.
     +      mcptr(ca,ires)  .gt. 0 .and.
     +      mcptr(n,ires+1) .gt. 0 .and.
     +      mcptr(o,ires)   .gt. 0) then
          plangl (ires) = tangle (mcptr(c,ires),
     +       mcptr(ca,ires),mcptr(n,ires+1),mcptr(o,ires),
     +       pdbxyz)
          icnt = icnt + 1
        end if
c
        if (mcptr(ca,ires)   .gt. 0 .and.
     +      mcptr(c,ires)  .gt. 0 .and.
     +      mcptr(n,ires+1) .gt. 0 .and.
     +      mcptr(ca,ires+1)   .gt. 0) then
          omega (ires) = tangle (mcptr(ca,ires),
     +       mcptr(c,ires),mcptr(n,ires+1),mcptr(ca,ires+1),
     +       pdbxyz)
          jcnt = jcnt + 1
        end if
c
  200   continue
c
      end do
c
      plangl (nres) = 999.99
      omega  (nres) = 999.99
c
      call jvalut (' Nr of residues with defined improper :',1,icnt)
      call jvalut (' Nr of residues with defined omega    :',1,jcnt)
c
      return
      end
c
c
c
      subroutine dozeta
c
      include 'oops2.incl'
c
      real tangle
c
      integer ires,icnt
c
code ...
c
      icnt = 0
c
      do ires=1,nres
c
        if (mcptr(ca,ires) .gt. 0 .and.
     +      mcptr(n,ires)  .gt. 0 .and.
     +      mcptr(c,ires)  .gt. 0 .and.
     +      mcptr(cb,ires) .gt. 0 .and.
     +      restyp(ires) .ne. 'PRO' ) then
          zeta (ires) = tangle (mcptr(ca,ires),
     +       mcptr(n,ires),mcptr(c,ires),mcptr(cb,ires),
     +       pdbxyz)
          icnt = icnt + 1
        else
          zeta (ires) = 999.99
        end if
c
      end do
c
      call jvalut (' Nr of residues with defined improper :',1,icnt)
c
      return
      end
c
c
c
      subroutine oopsts (text,n,x,nh,h,cnt,ave,sdv,xmin,xmax,xtot)
c
      implicit none
c
      real x(*),h(*)
      real ave,sdv,xmin,xmax,xtot
c
      integer cnt(*)
      integer n,nh
c
      character text*(*)
c
code ...
c
      write (*,6000) text
c
      if (n .le. 0) then
        call errcon ('No data points !')
        return
      end if
c
      call xstats (x,n,ave,sdv,xmin,xmax,xtot)
      write (*,6010) n,ave,sdv,xmin,xmax
c
      call histo (n,x,nh,h,cnt)
c
 6000 format (/1x,75('*')/1x,' Analysis of ',a/1x,75('*')/)
 6010 format (' Number of values .................... ',i20/
     +        ' Average value ....................... ',f20.3/
     +        ' Standard deviation .................. ',f20.3/
     +        ' Minimum value observed .............. ',f20.3/
     +        ' Maximum value observed .............. ',f20.3)
c
      return
      end
c
c
c
      subroutine ooplot (iunit,filnam,label,title,
     +                   nval,values,ymin,ymax,j,ave,sdv,xmin,xmax)
c
      implicit none
c
      real values(*),ymin,ymax,ave,sdv,xmin,xmax
c
      integer iunit,nval,ierr,i,j,leng1
c
      logical xinter
c
      character filnam*(*),label*(*),myline*128,title*(*)
c
code ...
c
      if (nval .lt. 2) return
c
      call remspa (filnam)
      call locase (filnam)
ccc      call textin (' O2D plot file ?',filnam)
      call textut (' O2D plot file :',filnam)
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening O2D plot file')
        return
      end if
c
 6900 format (a6,9i10)
 6910 format (a6,9f10.3)
c
      call stamp (myline)
      write (iunit,'(a,a)',err=9900)
     +  'REMARK ',myline(1:leng1(myline))
      write (iunit,'(a,a)',err=9900)
     +  'REMARK File ',filnam(1:leng1(filnam))
      write (iunit,'(a,a)',err=9900)
     +  'REMARK ',title(1:leng1(title))
      write (iunit,'(a,i8)',err=9900)
     +  'REMARK Number of observations = ',j
      write (iunit,'(a,f8.3)',err=9900)
     +  'REMARK Average value = ',ave
      write (iunit,'(a,f8.3)',err=9900)
     +  'REMARK Standard deviation = ',sdv
      write (iunit,'(a,f8.3)',err=9900)
     +  'REMARK Minimum observed value = ',xmin
      write (iunit,'(a,f8.3)',err=9900)
     +  'REMARK Maximum observed value = ',xmax
      write (iunit,6900,err=9900)  'NPOINT',nval
      write (iunit,'(a)',err=9900) 'XLABEL Residue number'
      write (iunit,'(a,a)',err=9900) 'YLABEL ',label(1:leng1(label))
      write (iunit,6900,err=9900)  'COLOUR',4
      write (iunit,6910,err=9900)  'XYVIEW',0.0,float(nval+1),ymin,ymax
      write (iunit,6910,err=9900)  'XLIMIT',1.0,1.0
      write (iunit,'(a)',err=9900) 'YVALUE *'
      write (iunit,'(1p,5e15.5)',err=9900) (values(i),i=1,nval)
      write (iunit,'(a)',err=9900) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
      return
c
 9900 continue
      call errcon ('While writing plot file')
      close (iunit)
c
      return
      end
c
c
c
c
      subroutine odbchk (molnam,extend,dbname,
     +             type,dbtype,nval,ival,lignore)
c
      implicit none
c
      integer nval,ival,length,l1,l2,leng1
c
      logical lerror,lignore
c
      character molnam*(*),extend*(*),dbname*(*),type*(*)
      character dbtype*(*),myline*128,oline*128,reply*1
c
code ...
c
      lerror = .false.
      lignore = .false.
c
      myline = molnam(1:leng1(molnam)) // 
     +         extend(1:leng1(extend))
      call locase (myline)
      call remspa (myline)
      oline = dbname
      call locase (oline)
      call remspa (oline)
c
      l1 = length(myline)
      l2 = length(oline)
      if (l1 .eq. l2) then
        if (myline(1:l1) .eq. oline(1:l1)) goto 100
      end if
c
      call prompt (' WARNING - inconsistent datablock name')
      call textut (' Expected    :',myline)
      call textut (' Encountered :',oline)
      lerror = .true.
c
  100 continue
c
      myline = type
      call locase (myline)
      call remspa (myline)
      oline = dbtype
      call locase (oline)
      call remspa (oline)
c
      l1 = length(myline)
      l2 = length(oline)
      if (l1 .eq. l2) then
        if (myline(1:l1) .eq. oline(1:l1)) goto 200
      end if
c
      call prompt (' WARNING - inconsistent datablock type')
      call textut (' Expected    :',myline)
      call textut (' Encountered :',oline)
      lerror = .true.
c
  200 continue
c
      if (nval .eq. ival) goto 300
c
      call prompt (
     +  ' WARNING - inconsistent number of datablock values')
      call jvalut (' Expected    :',1,nval)
      call jvalut (' Encountered :',1,ival)
      lerror = .true.
c
  300 continue
      if (.not. lerror) return
c
      reply = 'N'
      call textin (' Continue anyway ?',reply)
      call upcase (reply)
c
      if (reply .eq. 'Y') return
c
ccc      call errstp ('Inconsistent datablock')
c
      call errcon ('Ignoring inconsistent datablock')
      lignore = .true.
c
      return
      end
c
c
c
      subroutine dowatq
c
      include 'oops2.incl'
c
      real qqq
c
      integer ires,j
c
code ...
c
      qqq = -0.25/(resoln*resoln)
c
      do ires=1,nres
c
        if (swater(ires)) then
          j = resptr(1,ires)
          watqua (ires) = 100.0*qocc(j)*exp(qqq*bfac(j))
        end if
c
      end do
c
      return
      end
c
c
c
      subroutine rdbb
c
      include 'oops2.incl'
c
      integer maxbon
      parameter (maxbon = maxatm)
c
      real buf1(maxbon),buf2(maxbon)
      real dist,dd,q1,q2,q3,q4,q5,q6,q7,q8,q9
c
      integer i,ires,nbon,j
c
code ...
c
      do ires=1,nres
c
        nbon = 0
c
        do i=resptr(1,ires),resptr(2,ires)-1
          do j=i+1,resptr(2,ires)
            if (i .ne. j) then
              dd = dist (i,j,pdbxyz)
              if (dd .le. mxatat) then
                nbon = nbon + 1
                buf1 (nbon) = bfac(i)
                buf2 (nbon) = bfac(j)
              end if
            end if
          end do
        end do
c
        if (ires .gt. 1) then
          do i=resptr(1,ires),resptr(2,ires)
            do j=resptr(1,ires-1),resptr(2,ires-1)
              dd = dist (i,j,pdbxyz)
              if (dd .le. mxatat) then
                nbon = nbon + 1
                buf1 (nbon) = bfac(i)
                buf2 (nbon) = bfac(j)
              end if
            end do
          end do
        end if
c
        if (ires .lt. nres) then
          do i=resptr(1,ires),resptr(2,ires)
            do j=resptr(1,ires+1),resptr(2,ires+1)
              dd = dist (i,j,pdbxyz)
              if (dd .le. mxatat) then
                nbon = nbon + 1
                buf1 (nbon) = bfac(i)
                buf2 (nbon) = bfac(j)
              end if
            end do
          end do
        end if
c
        if (nbon .le. 0) then
          rmsdbb (ires) = 0.0
        else
          call xystat (buf1,buf2,nbon,q1,q2,q3,q4,q5,q6,q7,q8,q9)
          rmsdbb (ires) = q1
        end if
c
      end do
c
      return
      end
c
c
c
      subroutine getold (iunit,ierr)
c
      include 'oops2.incl'
c
      real xyz(3),bbb,qqq
c
      integer iunit,ierr,i,j,natold,nfound,n1,n2
c
      logical hetatm,lhydro
c
      character line*80,atmnax*4,atmtyx*6,atnmnx*6
c
code ...
c
      ierr = -1
c
      do i=1,natoms
        found (i) = .false.
      end do
c
      natold = 0
      nfound = 0
c
 6000 format (12x,a4,1x,a3,1x,a6,3x,3f8.3,f6.2,f6.2)
c
   10 continue
      read (iunit,'(a)',end=100,err=9999) line
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') then
        line = ' '//line
        call textut (line(1:7),line(8:))
        goto 10
      end if
c
      call upcase (line)
      read (line,6000,err=9999) atmnax,atmtyx,atnmnx,
     +  (xyz(j),j=1,3),qqq,bbb
c
      if (lhydro(atmnax)) goto 10
c
c ... O puts a "$" in front of HETATM residue names
c
      hetatm = (line(1:6) .eq. 'HETATM')
      call remspa (atnmnx)
      if (hetatm) atnmnx = '$'//atnmnx
c
      natold = natold + 1
c
ccc      write (*,'(1x,a4,1x,a6,1x,a6)') atmnax,atmtyx,atnmnx
c
c ... try to find the atom in the existing list of atoms
c
      do i=1,natoms
        if (.not. found(i)) then
          if (atnmnx .eq. atnmnm(i)) then
            if (atmnax .eq. atmnam(i)) then
ccc              if (atmtyx .eq. atmtyp(i)) then
c
                nfound = nfound + 1
                found (i) = .true.
                pdbold (1,i) = xyz(1)
                pdbold (2,i) = xyz(2)
                pdbold (3,i) = xyz(3)
                bold (i) = bbb
                qold (i) = qqq
c
                goto 10
c
ccc              end if
            end if
          end if
        end if
      end do
c
      write (*,6110) atmnax,atmtyx,atnmnx
 6110 format (' Not in current model : ',a4,1x,a6,1x,a6)
c
      goto 10
c
 9999 continue
      close (iunit)
      return
c
  100 continue
      close (iunit)
      call jvalut (' Nr of atoms read  :',1,natold)
      call jvalut (' Nr of atoms found :',1,nfound)
c
      if (nfound .lt. 10) then
        call errcon ('Too few atoms found -- wrong PDB file ?')
        return
      end if
c
c ... check for inserted and mutated residues
c
      write (*,*)
      do i=1,nres
        n1 = resptr(2,i) - resptr(1,i) + 1
        n2 = 0
        do j=resptr(1,i),resptr(2,i)
          if (found(j)) n2 = n2 + 1
        end do
        if (n2 .eq. 0) then
          insert (i) = .true.
          write (*,6100) 'Inserted',restyp(i),resnam(i)
        else if (n2 .ne. n1) then
          mutate (i) = .true.
          write (*,6100) 'Mutated',restyp(i),resnam(i)
        end if
      end do
c
 6100 format (1x,a10,' residue > ',2a6)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine doprev (usewat)
c
      include 'oops2.incl'
c
      integer i
c
      logical usewat
c
code ...
c
      do i=1,nres
c
c ... skip if water and user not interested in waters
c
        if ( (.not. usewat) .and. swater(i) ) goto 9000
c
c ... skip if inserted
c
        if (insert(i)) goto 9000
c
c ... if mutated, try to calc RMSD/B/Q and Phi/Psi distance
c
        if (mutate(i)) then
          call gkrmss (i)
          call gkphps (i)
c
c ... if "conserved", calc everything
c
        else
          call gkrmss (i)
c
          if (.not. swater(i)) then
            call gkphps (i)
            call gkchid (i)
          end if
        end if
c
 9000   continue
      end do
c
      return
      end
c
c
c
      subroutine gkrmss (ires)
c
      include 'oops2.incl'
c
      real dum,xd,xb,xq
c
      integer ires,i,nf
c
code ...
c
      nf = 0
      xd = 0.0
      xb = 0.0
      xq = 0.0
c
c ... if mutated, only use main chain atoms
c
      if (mutate(ires)) then
c
        if (mcptr(ca,ires) .le. 0) return
        if (.not.found(mcptr(ca,ires))) return
        if (mcptr(n,ires) .le. 0) return
        if (.not.found(mcptr(n,ires))) return
        if (mcptr(c,ires) .le. 0) return
        if (.not.found(mcptr(c,ires))) return
c
        i = mcptr(ca,ires)
        dum = (pdbxyz(1,i)-pdbold(1,i))**2 +
     +        (pdbxyz(2,i)-pdbold(2,i))**2 +
     +        (pdbxyz(3,i)-pdbold(3,i))**2
        xd = xd + dum
        xb = xb + ( ( bfac(i) - bold(i) ) **2 )
        xq = xq + ( ( qocc(i) - qold(i) ) **2 )
c
        i = mcptr(n,ires)
        dum = (pdbxyz(1,i)-pdbold(1,i))**2 +
     +        (pdbxyz(2,i)-pdbold(2,i))**2 +
     +        (pdbxyz(3,i)-pdbold(3,i))**2
        xd = xd + dum
        xb = xb + ( ( bfac(i) - bold(i) ) **2 )
        xq = xq + ( ( qocc(i) - qold(i) ) **2 )
c
        i = mcptr(c,ires)
        dum = (pdbxyz(1,i)-pdbold(1,i))**2 +
     +        (pdbxyz(2,i)-pdbold(2,i))**2 +
     +        (pdbxyz(3,i)-pdbold(3,i))**2
        xd = xd + dum
        xb = xb + ( ( bfac(i) - bold(i) ) **2 )
        xq = xq + ( ( qocc(i) - qold(i) ) **2 )
c
        nf = 3
c
        if (mcptr(o,ires) .gt. 0) then
          if (found(mcptr(o,ires))) then
            i = mcptr(o,ires)
            nf = nf + 1
            dum = (pdbxyz(1,i)-pdbold(1,i))**2 +
     +            (pdbxyz(2,i)-pdbold(2,i))**2 +
     +            (pdbxyz(3,i)-pdbold(3,i))**2
            xd = xd + dum
            xb = xb + ( ( bfac(i) - bold(i) ) **2 )
            xq = xq + ( ( qocc(i) - qold(i) ) **2 )
          end if
        end if
c
        if (mcptr(cb,ires) .gt. 0) then
          if (found(mcptr(cb,ires))) then
            i = mcptr(cb,ires)
            nf = nf + 1
            dum = (pdbxyz(1,i)-pdbold(1,i))**2 +
     +            (pdbxyz(2,i)-pdbold(2,i))**2 +
     +            (pdbxyz(3,i)-pdbold(3,i))**2
            xd = xd + dum
            xb = xb + ( ( bfac(i) - bold(i) ) **2 )
            xq = xq + ( ( qocc(i) - qold(i) ) **2 )
          end if
        end if
c
      else
c
c ... unchanged residue
c
        do i=resptr(1,ires),resptr(2,ires)
          if (found(i)) then
            nf = nf + 1
            dum = (pdbxyz(1,i)-pdbold(1,i))**2 +
     +            (pdbxyz(2,i)-pdbold(2,i))**2 +
     +            (pdbxyz(3,i)-pdbold(3,i))**2
            xd = xd + dum
            xb = xb + ( ( bfac(i) - bold(i) ) **2 )
            xq = xq + ( ( qocc(i) - qold(i) ) **2 )
          end if
        end do
c
      end if
c
      if (nf .le. 0) return
c
      dum = 1.0 / float(nf)
      xxrmsd (ires) = sqrt (xd * dum)
      xxrmsb (ires) = sqrt (xb * dum)
      xxrmsq (ires) = sqrt (xq * dum)
c
      return
      end
c
c
c
      subroutine gkphps (ires)
c
      include 'oops2.incl'
c
      real xd,xphi,xpsi,dum
      real dist,tangle
c
      integer ires,nc
c
code ...
c
c ... skip if not an amino acid or CA not found
c
      if (mcptr(ca,ires)   .le. 0) return
      if (.not. found(mcptr(ca,ires)))    return
c
      xd = 0.0
      nc = 0
c
c ... try to calc PHI
c
      if (phi(ires) .lt. 999.0) then
c
c ... check if the necessary atoms were found
c
        if (mcptr(ca,ires-1).le. 0) goto 200
        if (.not. found(mcptr(ca,ires-1))) goto 200
        if (dist(mcptr(ca,ires),mcptr(ca,ires-1),
     +      pdbold) .gt. mxcaca) goto 200
c
        if (mcptr(c,ires-1).le. 0) goto 200
        if (mcptr(c,ires).le. 0) goto 200
        if (mcptr(n,ires).le. 0) goto 200
        if (.not. found (mcptr(c,ires-1))) goto 200
        if (.not. found (mcptr(n,ires))) goto 200
        if (.not. found (mcptr(c,ires))) goto 200
c
        xphi = tangle (mcptr(c,ires-1),
     +     mcptr(n,ires),mcptr(ca,ires),mcptr(c,ires),
     +     pdbold)
c
        call fixang (xphi)
c
        call fixdif (phi(ires),xphi,dum)
        xd = xd + ( dum **2 )
        nc = nc + 1
ccc
ccc        print *,restyp(ires),resnam(ires)
ccc        print *,'PHI ',phi(ires),xphi,dum
c
      end if
c
c ... try to calc PSI
c
  200 continue
c
      if (psi(ires) .lt. 999.0) then
c
c ... check if the necessary atoms were found
c
        if (mcptr(ca,ires+1).le. 0) goto 200
        if (.not. found(mcptr(ca,ires+1))) goto 300
        if (dist(mcptr(ca,ires),mcptr(ca,ires+1),
     +      pdbold) .gt. mxcaca) goto 300
c
        if (mcptr(n,ires+1).le. 0) goto 300
        if (mcptr(c,ires).le. 0) goto 300
        if (mcptr(n,ires).le. 0) goto 300
        if (.not. found (mcptr(n,ires))) goto 300
        if (.not. found (mcptr(c,ires))) goto 300
        if (.not. found (mcptr(n,ires+1))) goto 300
c
        xpsi = tangle (mcptr(n,ires),
     +       mcptr(ca,ires),mcptr(c,ires),mcptr(n,ires+1),
     +       pdbold)
c
        call fixang (xpsi)
c
        call fixdif (psi(ires),xpsi,dum)
        xd = xd + ( dum **2 )
        nc = nc + 1
ccc
ccc        print *,'PSI ',psi(ires),xpsi,dum
c
      end if
c
c ... calc Phi/Psi distance
c
  300 continue
c
      if (nc .le. 0) return
c
      xxphps (ires) = sqrt (xd / float(nc))
ccc
ccc      print *,'PHI/PSI distance ',xxphps(ires)
c
      return
      end
c
c
c
      subroutine gkchid (ires)
c
      include 'oops2.incl'
c
      real xd,xchi1,xchi2,dum
      real tangle
c
      integer ires,iptr,jptr,nc
c
code ...
c
c ... skip if not an amino acid or CA not found
c
      if (mcptr(ca,ires) .le. 0) return
      if (.not. found(mcptr(ca,ires))) return
c
      xd = 0.0
      nc = 0
c
c ... try to calc CHI1
c
      if (mcptr(n,ires) .le. 0) goto 200
      if (.not. found(mcptr(n,ires))) goto 200
      if (mcptr(cb,ires) .le. 0) return
      if (.not. found(mcptr(cb,ires))) return
c
c ... try to find G(1) atom
c
      call getg1 (ires,iptr)
      if (iptr .le. 0) return
      if (.not. found (iptr)) return
c
      xchi1 = tangle (mcptr(n,ires),
     +       mcptr(ca,ires),mcptr(cb,ires),iptr,pdbxyz)
      call fixang (xchi1)
c
      xchi2 = tangle (mcptr(n,ires),
     +       mcptr(ca,ires),mcptr(cb,ires),iptr,pdbold)
      call fixang (xchi2)
c
      call fixdif (xchi1,xchi2,dum)
      xd = xd + ( dum **2 )
      nc = nc + 1
ccc
ccc      print *,restyp(ires),resnam(ires)
ccc      print *,'CHI1 ',xchi1,xchi2,dum
c
  200 continue
c
c ... try to find G(1) atom
c
      call getg1 (ires,iptr)
      if (iptr .le. 0) return
      if (.not. found (iptr)) return
c
c ... try to find G(2) atom
c
      call getg2 (ires,jptr)
      if (jptr .le. 0) goto 300
      if (.not. found (jptr)) goto 300
c
      xchi1 = tangle (mcptr(ca,ires),
     +       mcptr(cb,ires),iptr,jptr,pdbxyz)
      call fixang (xchi1)
c
      xchi2 = tangle (mcptr(ca,ires),
     +       mcptr(cb,ires),iptr,jptr,pdbold)
      call fixang (xchi2)
c
      call fixdif (xchi1,xchi2,dum)
      xd = xd + ( dum **2 )
      nc = nc + 1
ccc
ccc      print *,'CHI2 ',xchi1,xchi2,dum
c
  300 continue
c
      if (nc .le. 0) return
c
      xxchid (ires) = sqrt (xd / float(nc))
ccc
ccc      print *,'CHI distance ',xxchid(ires)
c
      return
      end
c
c
c
      subroutine getg1 (ires,iptr)
c
      include 'oops2.incl'
c
      integer ires,iptr,i
c
code ...
c
      iptr = -1
c
      do i=resptr(1,ires),resptr(2,ires)
        if (atmnam(i) .eq. ' CG ' .or.
     +      atmnam(i) .eq. ' CG1' .or.
     +      atmnam(i) .eq. ' SG ' .or.
     +      atmnam(i) .eq. ' OG ' .or.
     +      atmnam(i) .eq. ' OG1') then
          iptr = i
          return
        end if
      end do
c
      return
      end
c
c
c
      subroutine getg2 (ires,iptr)
c
      include 'oops2.incl'
c
      integer ires,iptr,i
c
code ...
c
      iptr = -1
c
      do i=resptr(1,ires),resptr(2,ires)
        if (atmnam(i) .eq. ' CD ' .or.
     +      atmnam(i) .eq. ' CD1' .or.
     +      atmnam(i) .eq. ' SD ' .or.
     +      atmnam(i) .eq. ' OD1' .or.
     +      atmnam(i) .eq. ' ND1') then
          iptr = i
          return
        end if
      end do
c
      return
      end
c
c
c
      subroutine pseudo (iunit)
c
      include 'oops2.incl'
c
      real x1,x2
c
      integer iunit,ires,ica,inr
c
      character achain*2
c
code ...
c
      do ires=1,nres
        ica = mcptr(ca,ires)
        if (ica .gt. 0) then
c
          x1 = phi(ires)
          if (x1.gt.999.0) x1 = 0.0
          x2 = psi(ires)
          if (x2.gt.999.0) x2 = 0.0
c
          call detaj (atnmnm(ica),achain,inr)
          write (iunit,6000,err=9990)
     +     ires,atmnam(ica),atmtyp(ica),achain,inr,x1,
     +     x2,100.0*pepflip(ires),rsc(ires),hib(ires)
c
        end if
      end do
c
 6000 format ('ATOM  ',i5,1x,a4,1x,a3,a2,i4,4x,3f8.3,2f6.2)
c
c ATOM   2202  CA  GLU B 136     -15.078 -12.989 -14.700  1.00 81.57   6
c 1234567890123456789012345678901234567890123456789012345678901234567890
c ATOM    137  CA  GLU   138     -61.507   0.000   0.000  0.00 37.31
c
      write (iunit,'(a3)',err=9990) 'END'
c
 9999 continue
      close (iunit)
      return
c
 9990 continue
      call errcon ('While writing pseudo-PDb file')
      goto 9999
c
      end
c
c
c
      subroutine psodb (iunit)
c
      implicit none
c
      real x1,x2,x3
c
      integer iunit
c
code ...
c
 6000 format (a,1x,3f8.2,1x,a)
c
      x1 = -180.0
      x2 =  180.0
      x3 =    0.0
c
      write (iunit,6000) '!'
      write (iunit,6000) '! Type: draw ramaflip'
      write (iunit,6000) '!'
      write (iunit,6000) '.GS_REAL R 27 (6(x,f8.5))'
      write (iunit,6000)
     +  '  0.03668  0.00000  0.00000  0.00000  0.04432 -0.00067'
      write (iunit,6000)
     +  ' -0.00330  0.00000  0.00123  0.04382  0.00683  0.00000'
      write (iunit,6000)
     +  '  0.00316 -0.00690  0.04371  0.00000  0.00000  0.00000'
      write (iunit,6000)
     +  '  0.00000  1.00000  0.09200  0.98200  0.00000  0.90000'
      write (iunit,6000)
     +  '  0.00000  0.00000  0.00000'
      write (iunit,6000) 'ramaflip t 18 72'
      write (iunit,6000) 'begin ramaflip'
      write (iunit,6000) ' colour white'
      write (iunit,6000) ' text_colour white'
      write (iunit,6000) ' move',x1,x1,x3
      write (iunit,6000) ' line',x2,x1,x3
      write (iunit,6000) ' line',x2,x2,x3
      write (iunit,6000) ' line',x1,x2,x3
      write (iunit,6000) ' line',x1,x1,x3
      write (iunit,6000) ' move',x1,x3,x3
      write (iunit,6000) ' line',x2,x3,x3
      write (iunit,6000) ' move',x3,x1,x3
      write (iunit,6000) ' line',x3,x2,x3
      write (iunit,6000) ' move',x3,x3,x3
      write (iunit,6000) ' line',x3,x3,300.0
      write (iunit,6000) ' text',190.0,0.0,0.0,'Phi'
      write (iunit,6000) ' text',0.0,190.0,0.0,'Psi'
      write (iunit,6000) ' text',0.0,0.0,300.0,'Pep_flip'
      write (iunit,6000) 'end_object'
      close (iunit)
c
      return
      end
c
c
c
      subroutine prevdb (iunit)
c
      include 'oops2.incl'
c
      integer i,iunit,ierr,leng1
c
      logical xinter
c
      character file*128,line*128,dbname*80
c
code ...
c
      file = molnam//'_prev_model.odb'
      call remspa (file)
      call locase (file)
      call textin (' O datablock file ?',file)
c
      call xopxua (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening O datablock file')
        return
      end if
c
 6000 format ('!',10(1x,a))
 6010 format (a,' R ',i6,' (9F8.3)')
 6020 format (9f8.3)
c
      call stamp (line)
      write (iunit,6000) line(1:leng1(line))
c
      dbname = molnam // '_residue_rmsd'
      call remspa (dbname)
      call upcase (dbname)
      write (iunit,6000) dbname(1:leng1(dbname))
      write (iunit,6010) dbname(1:leng1(dbname)),nres
      write (iunit,6020) (xxrmsd(i),i=1,nres)
c
      dbname = molnam // '_residue_rmsb'
      call remspa (dbname)
      call upcase (dbname)
      write (iunit,6000) dbname(1:leng1(dbname))
      write (iunit,6010) dbname(1:leng1(dbname)),nres
      write (iunit,6020) (xxrmsb(i),i=1,nres)
c
      dbname = molnam // '_residue_rmsq'
      call remspa (dbname)
      call upcase (dbname)
      write (iunit,6000) dbname(1:leng1(dbname))
      write (iunit,6010) dbname(1:leng1(dbname)),nres
      write (iunit,6020) (xxrmsq(i),i=1,nres)
c
      dbname = molnam // '_residue_phipsi'
      call remspa (dbname)
      call upcase (dbname)
      write (iunit,6000) dbname(1:leng1(dbname))
      write (iunit,6010) dbname(1:leng1(dbname)),nres
      write (iunit,6020) (xxphps(i),i=1,nres)
c
      dbname = molnam // '_residue_chi1chi2'
      call remspa (dbname)
      call upcase (dbname)
      write (iunit,6000) dbname(1:leng1(dbname))
      write (iunit,6010) dbname(1:leng1(dbname)),nres
      write (iunit,6020) (xxchid(i),i=1,nres)
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine outlie (n,val,cut,mode)
c
c ... OUTLIE - mode = 0 -> count residues >= cutoff
c              mode = 1 -> count residues <= cutoff
c
      implicit none
c
      integer n,mode,nout,i
c
      real val(n),cut,pout
c
code ...
c
      if (n .le. 0) return
c
      nout = 0
      do i=1,n
        if (mode .eq. 0) then
          if (val(i) .ge. cut) nout = nout + 1
        else
          if (val(i) .le. cut) nout = nout + 1
        end if
      end do
      pout = 100.0 * float(nout) / float(n)
c
      call ivalut (' Nr of outliers :',1,nout)
      call fvalut (' Percentage     :',1,pout)
c
      return
      end
c
c
c
      subroutine disulf ()
c
      include 'oops2.incl'
c
      real bridge(3,4),dd,dist,tangle
c
      integer i,j,ires,jres,k
c
      logical all4(4)
c
code ...
c
      ndisu = 0
c
      do ires=1,nres-1
        if (restyp(ires) .eq. 'CYS') then
          all4(1) = .false.
          all4(2) = .false.
          do i=resptr(1,ires),resptr(2,ires)
            if (atmnam(i) .eq. ' CB ') then
              do k=1,3
                bridge(k,1) = pdbxyz(k,i)
              end do
              all4(1) = .true.
            else if (atmnam(i) .eq. ' SG ') then
              do k=1,3
                bridge(k,2) = pdbxyz(k,i)
              end do
              all4(2) = .true.
            end if
            if (all4(1) .and. all4(2)) goto 10
          end do
          call errcon ('Missing atom(s) in Cys '//resnam(ires))
          goto 90
c
   10     continue
          do jres=ires+1,nres
            if (restyp(jres) .eq. 'CYS') then
              all4(3) = .false.
              all4(4) = .false.
              do j=resptr(1,jres),resptr(2,jres)
                if (atmnam(j) .eq. ' CB ') then
                  do k=1,3
                    bridge(k,4) = pdbxyz(k,j)
                  end do
                  all4(4) = .true.
                else if (atmnam(j) .eq. ' SG ') then
                  do k=1,3
                    bridge(k,3) = pdbxyz(k,j)
                  end do
                  all4(3) = .true.
                end if
                if (all4(3) .and. all4(4)) goto 20
              end do
              call errcon ('Missing atom(s) in Cys '//resnam(ires))
              goto 80
c
   20         continue
              dd = dist (2,3,bridge)
              if (dd .le. 3.2) then
                ndisu = ndisu + 1
                if (ndisu .gt. maxss) then
                  call errcon ('Too many disulfide bridges')
                  goto 99
                end if
                ssresi(1,ndisu) = ires
                ssresi(2,ndisu) = jres
                ssdist(ndisu) = dd
                sstors(ndisu) = tangle(1,2,3,4,bridge)
                write (*,6000) ndisu,resnam(ires),resnam(jres),dd
              end if
            end if
   80     continue
          end do
        end if
   90   continue
      end do
   99 continue
c
 6000 format (' Nr ',i3,' Cys ',a,' - ',a,' S=S (A) = ',f6.2)
c
      call jvalut (' Nr of disulfide bridges :',1,ndisu)
c
      return
      end
c
c
c
      subroutine appfil (iunit,myfile,nl,xl)
c
c ... append NL lines XL*(*) to an ASCII file of no more than
c     MAXLIN lines of 120 characters
c
      implicit none
c
      integer maxlin
      parameter (maxlin=1000)
c
      integer iunit,nl,ierr,nold,i,leng1
c
      character lines(maxlin)*120,filnam*120
      character*(*) myfile,xl(*)
c
code ...
c
      if (nl .le. 0) return
c
      filnam = myfile
      ierr = 0
      close (iunit)
      call xopxoa (iunit,filnam,.false.,ierr)
      nold = 0
c
   10 continue
      read (iunit,'(a)',end=20) lines(nold+1)
      nold = nold + 1
      goto 10
c
   20 continue
      close (iunit)
      call xopxua (iunit,filnam,.false.,ierr)
c
      if (nold .gt. 0) then
        do i=1,nold
          write (iunit,'(a)') lines(i)(1:leng1(lines(i)))
        end do
      end if
c
      do i=1,nl
        write (iunit,'(a)') xl(i)(1:leng1(xl(i)))
      end do
      close (iunit)
c
      return
      end
c
c
c
      subroutine whatif (iunit,lrelax)
c
c ... WHATIF - analyse WHAT IF diagnostics report
c
      include 'oops2.incl'
c
      integer iunit,ierr,nbah
c
      logical lrelax,missed,anymis
c
      character line*128
c
code ...
c
      anymis = .false.
c
c# 1 # Warning: Rounded coordinates detected
c 131 HOH  (HOH  ) V    254     1.000     1.000     1.000
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      write (*,*)
      call getlin (iunit,'Warning: Rounded coordinates',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No rounded coordinates detected')
        goto 10
      end if
      call textut (' Ouch :',line)
      line = ' Rounded coordinates : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   10 continue
c
c# 3 # Warning: Valine nomenclature problem
c   2 VAL  (   2 ) A
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      write (*,*)
      call getlin (iunit,'Warning: Valine nomenclature',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No incorrect Valine nomenclature')
        goto 20
      end if
      call textut (' Ouch :',line)
      line = ' Incorrect Valine nomenclature : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   20 continue
c
      write (*,*)
      call getlin (iunit,'Error: Threonine nomenclature',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No incorrect Threonine nomenclature')
        goto 30
      end if
      call textut (' Ouch :',line)
      line = ' Incorrect Threonine nomenclature : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   30 continue
c
      write (*,*)
      call getlin (iunit,'Error: Isoleucine nomenclature',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No incorrect Isoleucine nomenclature')
        goto 40
      end if
      call textut (' Ouch :',line)
      line = ' Incorrect Isoleucine nomenclature : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   40 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Leucine nomenclature',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No incorrect Leucine nomenclature')
        goto 50
      end if
      call textut (' Ouch :',line)
      line = ' Incorrect Leucine nomenclature : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   50 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Arginine nomenclature',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No incorrect Arginine nomenclature')
        goto 60
      end if
      call textut (' Ouch :',line)
      line = ' Incorrect Arginine nomenclature : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   60 continue
c
c# 13 # Warning: Chirality deviations detected
c   9 ALA  (   9 ) A    CA     -69.5    -94.3     34.4 Wrong hand
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      write (*,*)
      call getlin (iunit,'Warning: Chirality deviations',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No chirality deviations')
        goto 70
      end if
      call textut (' Ouch :',line)
      line = ' Chirality deviation : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   70 continue
c
      write (*,*)
      call getlin (iunit,'Error: Weights outside',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No impossible occupancies')
        goto 80
      end if
      call textut (' Ouch :',line)
      line = ' Impossible occupancy : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   80 continue
c
ccc      if (lrelax) goto 90
c
      write (*,*)
      call getlin (iunit,'Warning: Missing atoms',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No missing atoms')
        goto 90
      end if
      call textut (' Ouch :',line)
      line = ' Missing atoms : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
   90 continue
c
ccc      if (lrelax) goto 100
c
      write (*,*)
      call getlin (iunit,'Warning: Unusual bond lengths',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual bond lengths')
        goto 100
      end if
      call textut (' Ouch :',line)
      line = ' Unusual bond lengths : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  100 continue
c
ccc      if (lrelax) goto 110
c
      write (*,*)
      call getlin (iunit,'Warning: Unusual bond angles',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual bond angles')
        goto 110
      end if
      call textut (' Ouch :',line)
      line = ' Unusual bond angles : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  110 continue
c
      write (*,*)
      call getlin (iunit,'Error: Side chain planarity',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No non-planar side chains')
        goto 120
      end if
      call textut (' Ouch :',line)
      line = ' Non-planar side chains : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  120 continue
c
      write (*,*)
      call getlin (iunit,'Error: Connections to arom',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No atoms not co-planar with rings')
        goto 130
      end if
      call textut (' Ouch :',line)
      line = ' Atoms not co-planar with rings : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  130 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Torsion angle eval',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual torsion angles')
        goto 140
      end if
      call textut (' Ouch :',line)
      line = ' Unusual torsion angles : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  140 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Backbone torsion',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual backbone torsions')
        goto 150
      end if
      call textut (' Ouch :',line)
      line = ' Unusual backbone torsions : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  150 continue
c
      write (*,*)
      call getlin (iunit,'Error: Atoms too close to',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No atoms near symmetry axes')
        goto 160
      end if
      call textut (' Ouch :',line)
      line = ' Atoms near symmetry axes : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  160 continue
c
      write (*,*)
      call getlin (iunit,'Error: Abnormally short',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No short contacts')
        goto 170
      end if
      call textut (' Ouch :',line)
      line = ' Bumps : '
      call wifdo2 (iunit,line,nbah,lrelax,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  170 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Abnormal packing',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No abnormal packing environments')
        goto 180
      end if
      call textut (' Ouch :',line)
      line = ' Abnormal packing environments : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  180 continue
c
c ... PM - to do: stretch of abnormal packing (old and new) !!!
c
      write (*,*)
      call getlin (iunit,'Warning: Backbone oxygen',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual peptide oxygens')
        goto 190
      end if
      call textut (' Ouch :',line)
      line = ' Unusual peptide oxygens : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  190 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Unusual rotamers',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual rotamers')
        goto 200
      end if
      call textut (' Ouch :',line)
      line = ' Unusual rotamers : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  200 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Unusual backbone',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unusual backbone conformations')
        goto 210
      end if
      call textut (' Ouch :',line)
      line = ' Unusual backbone conformations : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  210 continue
c
      write (*,*)
      call getlin (iunit,'Error: Water clusters without',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No isolated water clusters')
        goto 220
      end if
      call textut (' Ouch :',line)
      line = ' Isolated water clusters : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  220 continue
c
      write (*,*)
      call getlin (iunit,'Warning: Water molecules need',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No waters in wrong asymmetric unit')
        goto 230
      end if
      call textut (' Ouch :',line)
      line = ' Waters in wrong asymmetric unit : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  230 continue
c
      write (*,*)
      call getlin (iunit,'Error: Water molecules without',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No water molecules without H-bonds')
        goto 240
      end if
      call textut (' Ouch :',line)
      line = ' Water molecules without H-bonds : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  240 continue
c
      write (*,*)
      call getlin (iunit,'Error: HIS, ASN, GLN side',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No H/N/Q side chain flips')
        goto 250
      end if
      call textut (' Ouch :',line)
      line = ' H/N/Q side chain flips : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  250 continue
c
      write (*,*)
      call getlin (iunit,'sfied hydrogen bond donors',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unsatisfied H-bond donors')
        goto 260
      end if
      call textut (' Ouch :',line)
      line = ' Unsatisfied H-bond donors : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  260 continue
c
      write (*,*)
      call getlin (iunit,'sfied hydrogen bond acceptors',line,ierr)
      if (ierr .ne. 0) then
        call prompt
     +    (' SUMMWIF - No unsatisfied H-bond acceptors')
        goto 270
      end if
      call textut (' Ouch :',line)
      line = ' Unsatisfied H-bond acceptors : '
      call wifdo (iunit,line,nbah,missed)
      if (missed) anymis = .true.
      call jvalut ((' SUMMWIF - '//line),1,nbah)
      if (nwif .ge. maxwif) return
  270 continue
c
      if (anymis) then
        write (*,'(99(a/))')
     +    ' ',
     +    ' --------------------------------------------------',
     +    ' WARNING !!! Not all diagnostics were printed by',
     +    ' What If/What_Check and therefore you may miss',
     +    ' important information during rebuilding.',
     +    ' ',
     +    ' To fix this problem, either issue the command',
     +    ' "setwif 593 100000" before running the "check"',
     +    ' command, or edit the What If/What_Check parameter',
     +    ' file PARAMS.FIG and change the value of parameter',
     +    ' number 593 to a large number (e.g., 100000).',
     +    ' ',
     +    ' Then re-run the check in What If/What_Check and',
     +    ' use the resulting "pdbout.txt" file as input for',
     +    ' this program',
     +    ' --------------------------------------------------',
     +    ' '
      end if
c
      return
      end
c
c
c
      subroutine wifdo (iunit,line,nn,missed)
c
      include 'oops2.incl'
c
      integer iunit,length,leng1,nn,i,idum
c
      logical missed
c
      character line*(*)
      character myline*128,resid*6,dollid*6,cins*1
c
code ...
c
c1405 PHE  ( 710C) B -3.3198
c2882 HOH  (HOH  ) W    797
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      nn = 0
      missed = .false.
c
   10 continue
      read (iunit,'(a)',end=100,err=100) myline
c
      if (myline(1:1) .eq. '#') then
        if (nn .eq. 0) call errcon (
     +         ' Could not find list of offenders ???')
        return
      end if
c
      if (myline(1:9) .eq. 'And so on') then
        call pretty (myline)
        myline = ' WARNING - List incomplete ! "' //
     +           myline(1:length(myline)) // '"'
        call prompt (myline)
        missed = .true.
      end if
c
      if (length(myline) .lt. 1) goto 10
      if (myline(11:11) .ne. '(') goto 10
      if (myline(17:17) .ne. ')') goto 10
c
      cins = ' '
      read (myline,'(18x,a1)',err=200) resid(1:1)
      read (myline,'(11x,i4)',err=20) idum
      cins = myline(16:16)
      goto 30
   20 continue
      read (myline,'(20x,i6)',err=200) idum
      cins = myline(27:27)
   30 continue
      write (resid(2:6),'(i5)') idum
      call remspa (resid)
      resid = resid(1:length(resid)) // cins
      call upcase (resid)
      dollid = '$' // resid
c
ccc      print *,'RESIDUE |',resid,'|'
c
      do i=1,nres
        if (resnam(i) .eq. resid .or.
     +      resnam(i) .eq. dollid) then
          call pretty (myline)
          nwif = nwif + 1
          if (nwif .gt. maxwif) then
            nwif = maxwif
            call errcon ('Too many WHAT IF diagnostics')
            call jvalut (' Maximum allowed :',1,maxwif)
            return
          end if
          nn = nn + 1
          wifptr (nwif) = i
          wifmes (nwif) = line(1:leng1(line)) // ' ' // myline
          call pretty (wifmes(nwif))
          goto 10
        end if
      end do
      call errcon (' Residue not recognised ???')
      call textut (' >',myline)
      goto 10
c
  100 continue
      call errcon ('While reading WHAT IF file')
      return
c
  200 continue
      call errcon ('While trying to recognise residue')
      call textut (' >',myline)
      goto 10
c
      end
c
c
c
      subroutine wifdo2 (iunit,line,nn,lrelax,missed)
c
      include 'oops2.incl'
c
      real xdum
c
      integer iunit,length,leng1,nn,i,idum,ierr
c
      logical lrelax,missed
c
      character line*(*)
      character myline*128,resid*6,resid2*6,dolid1*6,dolid2*6
      character cins*1
c
code ...
c
      nn = 0
      missed = .false.
c
c   8 LEU  (   8 ) A    O    --    9 ALA  (   9 ) A    CB     1.547   1.253 INTRA
c 133 HOH  (HOH  ) V    254  --  133 HOH  (HOH  ) V    262    0.986   1.414 INTRA
c 161 CYS  ( 245 ) A    SG   --  752 HOH  (HOH  ) C    600    0.290   2.560 INTRA HB
c 641 GLN  ( 362 ) K    NE2  --  758 HOH  (HOH  ) M    631    0.112   2.438 INTER BF
c 641 GLN  ( 362 ) K    NE2  --  758 HOH  (HOH  ) M    631    0.112   2.438 INTER B2
c 641 GLN  ( 362 ) K    NE2  --  758 HOH  (HOH  ) M    631    0.112   2.438 INTER B3
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c what_check:
c  12 GLU  (  13 ) A    CG    --   57 VAL  (  58 ) A    O       0.279   2.521 INTER
c what if:
c  12 GLU  (  13 ) A    CG   --   57 VAL  (  58 ) A    O      0.279   2.521 INTER
c
   10 continue
      read (iunit,'(a)',end=100,err=100) myline
c
      if (myline(1:1) .eq. '#') then
        if (nn .eq. 0) call errcon (
     +         ' Could not find list of offenders ???')
        return
      end if
c
      if (myline(1:9) .eq. 'And so on') then
        call pretty (myline)
        myline = ' WARNING - List incomplete ! "' //
     +           myline(1:length(myline)) // '"'
        call prompt (myline)
        missed = .true.
      end if
c
      if (length(myline) .lt. 1) goto 10
      if (myline(11:11) .ne. '(') goto 10
      if (myline(17:17) .ne. ')') goto 10
c
c ... @ CSHL @ 2003-10-26 - handle whatcheck (format partly shifted by one character)
c
      if (myline(43:43) .eq. '(') then
        myline = myline(1:28) // myline(30:)
ccc        print *,'SHIFTING 28-33 |',myline(28:33),'| 40-45 |',
ccc     +    myline(40:45),'|'
ccc      else
ccc        print *,'NOT SHIFTING 28-33 |',myline(28:33),'| 40-45 |',
ccc     +    myline(40:45),'|'
      end if
c
c ... skip small violations if relaxed mode
c
      if (lrelax .and. bmpcut .gt. 0.0) then
        call str2r (myline(60:66),xdum,ierr)
        if (ierr .eq. 0) then
          if (xdum .lt. bmpcut) goto 10
        end if
      end if
c
      idum = index (myline,'INTRA')
      if (idum .gt. 0) then
        myline = myline (1:idum-1) // myline (idum+5:)
      end if
c
      idum = index (myline,'INTER')
      if (idum .gt. 0) then
        myline = myline (1:idum-1) // 'SYMM' // myline (idum+5:)
      end if
c
      cins = ' '
      read (myline,'(18x,a1)',err=200) resid(1:1)
      read (myline,'(11x,i4)',err=20) idum
      cins = myline(16:16)
      goto 30
   20 continue
      read (myline,'(20x,i6)',err=200) idum
      cins = myline(27:27)
   30 continue
      write (resid(2:6),'(i5)') idum
      call remspa (resid)
      resid = resid(1:length(resid)) // cins
      call remspa (resid)
      call upcase (resid)
      dolid1 = '$' // resid
c
      cins = ' '
      read (myline,'(49x,a1)',err=200) resid2(1:1)
      read (myline,'(42x,i4)',err=120) idum
      cins = myline(47:47)
      goto 130
  120 continue
      read (myline,'(51x,i6)',err=200) idum
      cins = myline(58:58)
  130 continue
      write (resid2(2:6),'(i5)') idum
      call remspa (resid2)
      resid2 = resid2(1:length(resid2)) // cins
      call remspa (resid2)
      call upcase (resid2)
      dolid2 = '$' // resid2
c
ccc      print *,'RESIDUES |',resid,'|',resid2,'|'
c
      do i=1,nres
        if (resnam(i) .eq. resid .or.
     +      resnam(i) .eq. dolid1) then
          call pretty (myline)
          nwif = nwif + 1
          if (nwif .gt. maxwif) then
            nwif = maxwif
            call errcon ('Too many WHAT IF diagnostics')
            call jvalut (' Maximum allowed :',1,maxwif)
            return
          end if
          nn = nn + 1
          wifptr (nwif) = i
          wifmes (nwif) = line(1:leng1(line)) // ' ' // myline
          call pretty (wifmes(nwif))
          goto 40
        end if
      end do
      call errcon (' Residue 1 not recognised ???')
      call textut (' >',myline)
c
   40 continue
      do i=1,nres
        if (resnam(i) .eq. resid2 .or.
     +      resnam(i) .eq. dolid2) then
          call pretty (myline)
          nwif = nwif + 1
          if (nwif .gt. maxwif) then
            nwif = maxwif
            call errcon ('Too many WHAT IF diagnostics')
            call jvalut (' Maximum allowed :',1,maxwif)
            return
          end if
          nn = nn + 1
          wifptr (nwif) = i
          wifmes (nwif) = line(1:leng1(line)) // ' ' // myline
          call pretty (wifmes(nwif))
          goto 10
        end if
      end do
      call errcon (' Residue 2 not recognised ???')
      call textut (' >',myline)
      goto 10
c
  100 continue
      call errcon ('While reading WHAT IF file')
      return
c
  200 continue
      call errcon ('While trying to recognise residue')
      call textut (' >',myline)
      goto 10
c
      end
c
c
c
      subroutine getlin (iunit,target,result,ierr)
c
c ... GETLIN - find line in file (first rewind)
c
      implicit none
c
      integer iunit,ierr,lt,leng1,ll
c
      character target*(*),result*(*),line*128
c
code ...
c
      rewind (iunit)
      ierr = 0
c
      result = ' '
      lt = leng1 (target)
   10 continue
      read (iunit,'(a)',end=100,err=100) line
      ll = leng1 (line)
      if (ll .lt. lt) goto 10
      if (index(line(1:ll),target(1:lt)) .gt. 0) then
        result = line
        return
      end if
c
      goto 10
c
  100 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine gempty (iunit,ierr)
c
c ... GEMPTY - find next empty line in file
c
      implicit none
c
      integer iunit,ierr,length
c
      character line*10
c
code ...
c
      ierr = 0
c
   10 continue
      read (iunit,'(a)',end=100,err=100) line
      if (length(line) .lt. 1) return
      goto 10
c
  100 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine prdb (t1,t2,t3,t4)
c
      implicit none
c
      character*(*) t1,t2,t3,t4
c
code ...
c
      call textut (' Datablock :',t1)
      call textut (' Data type :',t2)
      call textut (' Number    :',t3)
      call textut (' Format    :',t4)
c
      return
      end
c
c
c
      subroutine ressub (iunit)
c
      include 'oops2.incl'
c
      integer iunit,istat,nopt,ierr,nl,nok,i,nn(-1:3)
c
      character line*80,optpar(3)*10,mytyp*6,mynam*6
c
code ...
c
      nl = 0
      nok = 0
   10 continue
      read (iunit,'(a)',end=99,err=99) line
      nl = nl + 1
      call extrop (line,nopt,3,optpar,ierr)      
      if (ierr .ne. 0) goto 10
      if (nopt .ne. 3) goto 10
      mytyp = optpar (1)
      mynam = optpar (2)
      call str2i (optpar(3),istat,ierr)
      if (ierr .ne. 0) goto 10
      if (istat .lt. -1 .or. istat .gt. 3) goto 10
      call upcase (mytyp)
      call remspa (mytyp)
      call upcase (mynam)
      call remspa (mynam)
      do i=1,nres
        if (resnam(i) .eq. mynam) then
          if (restyp(i) .eq. mytyp) then
            oldsub (i) = istat
            nok = nok + 1
            goto 10
          end if
        end if
      end do
      call pretty (line)
      call textut (' Residue not recognised :',line)
      goto 10
c
   99 continue
      call jvalut (' Lines read from table :',1,nl)
      call jvalut (' Residues recognised   :',1,nok)
      do i=-1,3
        nn (i) = 0
      end do
      do i=1,nres
        if (oldsub(i) .ge. -1 .and. oldsub(i) .le. 3) then
          nn(oldsub(i)) = nn(oldsub(i)) + 1
        end if
      end do
      call jvalut (' Unexamined            :',1,nn(-1))
      call jvalut (' Good dens, good fit   :',1,nn(0))
      call jvalut (' Good dens, poor fit   :',1,nn(1))
      call jvalut (' Poor dens, poor fit   :',1,nn(2))
      call jvalut (' No dens, no fit       :',1,nn(3))
c
      return
      end
c
c
c
      subroutine skipem (iunit)
c
c ... SKIPEM - skip header lines in ODB files that begin with
c              an Exclamation Mark
c
      implicit none
c
      integer iunit
c
      character ch*1
c
code ...
c
   10 continue
      read (iunit,'(a1)',end=99,err=99) ch
      if (ch .eq. '!') goto 10
c
   99 continue
      backspace (iunit)
c
      return
      end
c
c
c
      subroutine salpha (x)
c
c ... SALPHA - set up penta-residue ALPHA helix CA coordinates
c     (from O's alpha template)
c
      real x(3,5)
c
code ...
c
      x(1,1) = 3.633
      x(2,1) = 15.082
      x(3,1) = 31.410
c
      x(1,2) = 3.847
      x(2,2) = 16.332
      x(3,2) = 27.802
c
      x(1,3) = 4.017
      x(2,3) = 12.736
      x(3,3) = 26.460
c
      x(1,4) = 0.906
      x(2,4) = 11.724
      x(3,4) = 28.417
c
      x(1,5) = -1.006
      x(2,5) = 14.798
      x(3,5) = 27.217
c
      return
      end
c
c
c
      subroutine sbeta (x)
c
c ... SBETA - set up penta-residue BETA strand CA coordinates
c     (from O's beta template)
c
      real x(3,5)
c
code ...
c
      x(1,1) = 43.840
      x(2,1) = 28.467
      x(3,1) = -6.991
c
      x(1,2) = 45.504
      x(2,2) = 26.414
      x(3,2) = -4.284
c
      x(1,3) = 43.734
      x(2,3) = 23.423
      x(3,3) = -2.720
c
      x(1,4) = 44.607
      x(2,4) = 21.760
      x(3,4) = 0.604
c
      x(1,5) = 43.505
      x(2,5) = 18.374
      x(3,5) = 1.905
c
      return
      end
c
c
c
      subroutine slefth (x)
c
c ... SLEFTH - set up tetra-residue LEFT-Handed helix CA coordinates
c     (from GJK's lefth template)
c
      real x(3,4)
c
code ...
c
      x(1,1) = 0.000
      x(2,1) = 0.000
      x(3,1) = 0.000
c
      x(1,2) = -0.972
      x(2,2) = 0.114
      x(3,2) = -3.703
c
      x(1,3) = -4.115
      x(2,3) = 2.184
      x(3,3) = -2.993
c
      x(1,4) = -5.086
      x(2,4) = 0.000
      x(3,4) = 0.000
c
      return
      end
c
c
c
      subroutine yasspa ()
c
      include 'oops2.incl'
c
      real resxyz(3,maxres),calpha(3,5),cbeta(3,5),clefth(3,4)
      real rt(12),rmsd1,rmsd2,rmsd3,x
c
      integer cnt(-1:3),i,j,ierr
c
      character line*128
c
c ...
c
ccc      ssenam (-1) = 'Non-protein'
      ssenam (-1) = '-'
      ssenam (0)  = 'Loop or turn'
      ssenam (1)  = 'Alpha helix'
      ssenam (2)  = 'Beta strand'
      ssenam (3)  = 'Left-handed helix'
c
      call salpha (calpha)
      call sbeta  (cbeta)
      call slefth (clefth)
c
      do i=1,nres
        if (mcptr(ca,i) .gt. 0) then
          sstype (i) = 0
          j = mcptr(ca,i)
          resxyz (1,i) = pdbxyz (1,j)
          resxyz (2,i) = pdbxyz (2,j)
          resxyz (3,i) = pdbxyz (3,j)
        else
          sstype (i) = -1
          resxyz (1,i) = 999.99
          resxyz (2,i) = 999.99
          resxyz (3,i) = 999.99
        end if
      end do
c
      do i=1,nres
c
        if (mcptr (ca,i) .le. 0) goto 10
c
        if (i .lt. 3 .or. i .gt. (nres-2)) goto 10
c
        if (mcptr(ca,i-2) .le. 0) goto 10
        if (mcptr(ca,i-1) .le. 0) goto 10
        if (mcptr(ca,i+1) .le. 0) goto 10
        if (mcptr(ca,i+2) .le. 0) goto 10
c
        call lsqgjk (calpha,resxyz(1,i-2),5,rmsd1,rt,ierr)
        call lsqgjk (cbeta, resxyz(1,i-2),5,rmsd2,rt,ierr)
        call lsqgjk (clefth,resxyz(1,i-2),4,rmsd3,rt,ierr)
        if (rmsd1 .le. yalcut) then
          sstype (i-1) = 1
          sstype (i)   = 1
          sstype (i+1) = 1
        else if (rmsd2 .le. ybecut) then
          sstype (i-1) = 2
          sstype (i)   = 2
          sstype (i+1) = 2
        else if (rmsd3 .le. ylhcut) then
          sstype (i-1) = 3
          sstype (i)   = 3
        end if
c
   10   continue
c
      end do
c
      do i=-1,3
        cnt (i) = 0
      end do
c
      do i=1,nres
        j = sstype (i)
        cnt (j) = cnt (j) + 1
      end do
c
      call prompt ('0YASSPA results:')
      do i=-1,3
        line = ' Nr residues of type '//ssenam(i)//' : '
        call jvalut (line,1,cnt(i))
      end do
c
      j = cnt(0)+cnt(1)+cnt(2)+cnt(3)
      if (j .gt. 0) then
        x = 100.0 * float(cnt(1)) / float(j)
        call fvalut (' %-age alpha :',1,x)
        x = 100.0 * float(cnt(2)) / float(j)
        call fvalut (' %-age beta  :',1,x)
        x = 100.0 * float(cnt(3)) / float(j)
        call fvalut (' %-age lefth :',1,x)
      end if
c
      return
      end
