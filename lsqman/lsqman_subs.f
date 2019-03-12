c
c ... subroutines for LSQMAN - Gerard Kleywegt
c
      subroutine allocm (string,imsk,ierr)
c
c ... allocate a datamol
c
      include 'lsqman.incl'
c
      integer length,whichm
      integer i,imsk,ierr
c
      character string*(*)
c
code ...
c
      ierr = -1
c
      if (length(string) .lt. 1) return
c
      do i=1,maxmol
        if (.not. incore(i)) then
          imsk = i
          goto 910
        end if
      end do
c
      call errcon ('No more mols available')
      return
c
  910 continue
      name (imsk) = string
c
      call upcase (name(imsk))
      if (whichm(name(imsk)) .ne. imsk) then
        call errcon ('Invalid mol name (empty or not unique)')
        name (imsk) = '*&^$'
        return
      end if
c
      if (incore(imsk)) then
        call errcon ('Mol in use; DELETE it first')
        return
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      integer function whichm (nam)
c
c ... which datamol does the name "nam" correspond to ?
c
c ... if "*", then return 0, meaning ALL mols
c     if okay, return index of mol
c     otherwise:
c     -1 if length = 0
c     -2 if duplicate name
c     -3 if not found
c
      include 'lsqman.incl'
c
      integer i,ll,imsk,length
c
      character nam*(*)
c
code ...
c
      whichm = -1
      ll = length(nam)
      if (ll .le. 0) return
c
c ... is it the wildcard ?
c
      if (ll .eq. 1 .and. nam(1:1) .eq. '*') then
        whichm = 0
        return
      end if
c
      call upcase (nam)
      whichm = -2
      imsk = 0
      do i=1,maxmol
        if (nam(1:ll) .eq. name(i)(1:ll)) then
          if (imsk .ne. 0) return
          imsk = i
        end if
      end do
c
      whichm = -3
      if (imsk .eq. 0) return
c
      whichm = imsk
c
      return
      end
c
c
c
      subroutine selecm (nam,ierr)
c
c ... figure out which mol(s) are to be selected
c
      include 'lsqman.incl'
c
      integer ierr,i,imsk,whichm,nmask
c
      character nam*(*)
c
code ...
c
      ierr = -1
c
      nmask = 0
      do i=1,maxmol
        select (i) = .false.
        if (incore(i)) nmask = nmask + 1
      end do
c
      if (nmask .le. 0) then
        call errcon ('No mols in memory')
        return
      end if
c
      imsk = whichm(nam)
c
      if (imsk .lt. 0 .or. imsk .gt. maxmol) then
        call errcon ('Invalid mol name')
        return
      end if
c
      if (imsk .eq. 0) then
        do i=1,maxmol
          select (i) = incore (i)
        end do
        ierr = 0
      else
        if (incore(imsk)) then
          select (imsk) = .true.
          ierr = 0
        else
          call errcon ('Selected mol not in memory')
        end if
      end if
c
      return
      end
c
c
c
      subroutine molin (imol,wchn,watm,ierr)
c
      include 'lsqman.incl'
c
      integer imol,ierr,nlines,nat,kk,ichain,i,jchain,iprev
      integer nalts,nstrip,ns2,nwcstr,nwastr
c
      logical xinter,nmr,lhydro,lunder,lalchn,lalatm
c
      character key*6,line*128,chain*1,prev*1,segid*4,wchn*1,watm*4
c
code ...
c
      ierr = 0
c
c ... open files as READONLY
c
ccc      call xopenr (iunit,file(imol),'O','F',.true.,xinter(),ierr)
c
      call xopxoa (iunit,file(imol),xinter(),ierr)
      if (ierr .ne. 0) return
c
      ichain = ichar('A') - 1
      jchain = 0
      prev = '?'
      segid = '????'
      nmr = .false.
      iprev = 999999
      lunder = .false.
c
      lalchn = (wchn .eq. '*')
      lalatm = (watm .eq. '*')
c
c ... read the PDB file
c
 4711 format (a6,a)
 4713 format (i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2,6x,a4)
c
      nlines = 0
      nstrip = 0
      ns2 = 0
      nat = 0
      nalts = 0
      nwcstr = 0
      nwastr = 0
c
   20 continue
      read (iunit,4711,end=21,err=997) key,line
      call upcase (key)
      nlines = nlines + 1
c
      if (key .eq. 'ENDMDL' .and. jchain .gt. 0 .and. nmr) then
        if (.not. lnmral) then
          call prompt (' Skipping all but first NMR model')
          goto 21
        end if
      end if
c
      if (key .eq. 'MODEL ') then
        if (.not. nmr) call prompt (' Multiple NMR models')
        nmr = .true.
      end if
c
      if (key .eq. 'MODEL ' .and. chamod .eq. 'RE') then
        ichain = ichain + 1
        if (ichain .gt. ichar('Z')) then
          ichain = ichain - 1
          call errcon ('Too many MODELs; max = 26; rest skipped')
          goto 21
        end if
        chain = char(ichain)
        prev = chain
        write (*,'(1x,a,i3,a,a1)')
     +    'NMR model ',ichain - ichar('A') + 1,
     +    ' becomes chain ',chain
c
        jchain = jchain + 1
        chnptr (1,jchain,imol) = nat + 1
        chname (jchain,imol) = chain
      end if
c
c ... read cell from CRYST1 card if present
c
      if (key .eq. 'CRYST1') then
        read (line,*) (cell(i,imol),i=1,6)
ccc        call fvalut (' Cell :',6,cell(1,imol))
        call pdbinfo (key,line)
        goto 20
      end if
c
      if (key .ne. 'ATOM  ' .and.
     +    key .ne. 'HETATM') then
        call pdbinfo (key,line)
        goto 20
      end if
c
      if (key .eq. 'HETATM') then
        ns2 = ns2 + 1
        if (.not. ltkeep) goto 20
      end if
c
cATOM   1744  O   CYS B 217       4.341  58.615  38.325  1.00 14.18           O  
cATOM   1745  CB ACYS B 217       4.376  61.240  36.361  0.54 14.71           C  
cATOM   1746  CB BCYS B 217       4.545  61.115  36.319  0.45 15.22           C  
c      1234567890123456789012345678901234567890123456789012345678901234567890
c               1         2         3         4         5         6         7
c
      if (line(11:11) .ne. ' ' .and. line(11:11) .ne. 'A') then
        nalts = nalts + 1
        goto 20
      end if
c
      nat = nat + 1
      call upcase (line)
c
      if (nat .gt. maxatm) then
        call errcon ('Too many atoms -- rest skipped')
        call jvalut (' Max nr of atoms =',1,maxatm)
        nat = nat - 1
        goto 21
      end if
c
      kk = nat
      read (line,4713,err=5372) atomnr(kk,imol),
     +  atmnam(kk,imol),resnam(kk,imol),achain(kk,imol),
     +  iresid(kk,imol),atmxyz(1,kk,imol),atmxyz(2,kk,imol),
     +  atmxyz(3,kk,imol),qatom(kk,imol),batom(kk,imol),
     +  axplor(kk,imol)
c
c ... read a particular chain ?
c
      if (.not. lalchn) then
        if (achain(kk,imol) .ne. wchn) then
          nat = nat - 1
          nwcstr = nwcstr + 1
          goto 20
        end if
      end if
c
c ... read a particular atom type ?
c
      if (.not. lalatm) then
        if (atmnam(kk,imol) .ne. watm) then
          nat = nat - 1
          nwastr = nwastr + 1
          goto 20
        end if
      end if
c
      if (chamod .eq. 'NO') then
        if (achain(kk,imol) .eq. ' ') then
          achain(kk,imol) = '_'
          lunder = .true.
        end if
      end if
c
c ... strip hydrogens ?
c
      if (lhydro(atmnam(kk,imol))) then
        nstrip = nstrip + 1
        if (.not. lhkeep) then
          nat = nat - 1
          goto 20
        end if
      end if
c
      if (chamod .eq. 'RE') then
        if ((.not. nmr) .and. achain(kk,imol) .ne. prev) then
          ichain = ichain + 1
          if (ichain .gt. ichar('Z')) then
            ichain = ichain - 1
            call errcon ('Too many CHAINs; max = 26; rest skipped')
            goto 21
          end if
          chain = char(ichain)
          prev = achain (kk,imol)
          write (*,'(1x,a,a1,a,a1)')
     +      'Old chain |',prev,'| becomes chain ',chain
          jchain = jchain + 1
          chnptr (1,jchain,imol) = nat
          chname (jchain,imol) = chain
        end if
      else if (chamod .eq. 'OR' .or.
     +         chamod .eq. 'NO') then
        if (achain(kk,imol) .ne. prev) then
          ichain = ichain + 1
          if (ichain .gt. ichar('Z')) then
            ichain = ichain - 1
            call errcon ('Too many CHAINs; max = 26; rest skipped')
            goto 21
          end if
          chain = achain(kk,imol)
          prev = achain (kk,imol)
          write (*,'(1x,a,a1,a)')
     +      'Old chain name |',prev,'| kept'
          jchain = jchain + 1
          chnptr (1,jchain,imol) = nat
          chname (jchain,imol) = chain
        end if
      else if (chamod .eq. 'BR') then
        if (iresid(kk,imol) .ne. iprev) then
          if (iresid(kk,imol) .eq. (iprev+1) ) then
            iprev = iprev + 1
          else
            ichain = ichain + 1
            if (ichain .gt. ichar('Z')) then
              ichain = ichain - 1
              call errcon ('Too many CHAINs; max = 26; rest skipped')
              goto 21
            end if
            chain = char(ichain)
            iprev = iresid (kk,imol)
            write (*,'(1x,a,a1,a,a3,i6)')
     +         'New chain name |',chain,'| at residue ',
     +         resnam(kk,imol),iprev
            jchain = jchain + 1
            chnptr (1,jchain,imol) = nat
            chname (jchain,imol) = chain
          end if
        end if
      else if (chamod .eq. 'LO') then
        if (iresid(kk,imol) .ne. iprev) then
          if (iresid(kk,imol) .gt. iprev ) then
            iprev = iresid(kk,imol)
          else
            ichain = ichain + 1
            if (ichain .gt. ichar('Z')) then
              ichain = ichain - 1
              call errcon ('Too many CHAINs; max = 26; rest skipped')
              goto 21
            end if
            chain = char(ichain)
            iprev = iresid (kk,imol)
            write (*,'(1x,a,a1,a,a3,i6)')
     +         'New chain name |',chain,'| at residue ',
     +         resnam(kk,imol),iprev
            jchain = jchain + 1
            chnptr (1,jchain,imol) = nat
            chname (jchain,imol) = chain
          end if
        end if
      else if (chamod .eq. 'XP') then
        if (axplor(kk,imol) .ne. segid) then
          ichain = ichain + 1
          if (ichain .gt. ichar('Z')) then
            ichain = ichain - 1
            call errcon ('Too many CHAINs; max = 26; rest skipped')
            goto 21
          end if
          segid = axplor(kk,imol)
          chain = char(ichain)
          write (*,'(1x,a,a4,a,a1)')
     +      'XPLOR SEGId |',segid,'| becomes chain ',chain
          jchain = jchain + 1
          chnptr (1,jchain,imol) = nat
          chname (jchain,imol) = chain
        end if
      end if
      achain (kk,imol) = chain
c
      goto 20
c
c ... read from file error
c
  997 continue
      close (iunit)
      call errcon ('While reading PDB file')
      call jvalut (' Line nr :',1,nlines)
      call jvalut (' Atom nr :',1,nat)
      call textut (' Line :',line)
      ierr = -1
      return
c
c ... read from string error
c
 5372 continue
      close (iunit)
      call errcon ('While reading PDB file')
      call jvalut (' Line nr :',1,nlines)
      call jvalut (' Atom nr :',1,nat)
      call textut (' Line :',line)
      ierr = -1
      return
c
   21 continue
      call jvalut (' Nr of lines read from file :',1,nlines)
      call jvalut (' Nr of atoms in molecule    :',1,nat)
      if (nat .lt. 1) then
        call errcon ('No atoms found in PDB file !')
        close (iunit)
        ierr = -2
        return
      end if
      ichain = ichain - ichar('A') + 1
      call jvalut (' Nr of chains or models     :',1,ichain)
      if (lunder) then
        call prompt (' Blank chain IDs replaced by _underscores_')
      end if
c
      if (lhkeep) then
        call jvalut (' Nr of hydrogen atoms       :',1,nstrip)
      else
        call jvalut (' Stripped hydrogen atoms    :',1,nstrip)
      end if
c
      if (ltkeep) then
        call jvalut (' Nr of HETATMs              :',1,ns2)
      else
        call jvalut (' Stripped HETATMs           :',1,ns2)
      end if
c
      if (.not. lalchn) then
        call jvalut (' Stripped unwanted chains   :',1,nwcstr)
      end if
c
      if (.not. lalatm) then
        call jvalut (' Stripped unwanted atoms    :',1,nwastr)
      end if
c
      call jvalut (' Stripped alt. conf. atoms  :',1,nalts)
c
      natoms (imol) = nat
      mmodel (imol) = nmr
      nchain (imol) = ichain
c
      chnptr (2,ichain,imol) = nat
      if (ichain .gt. 1) then
        do i=1,ichain-1
          chnptr (2,i,imol) = chnptr (1,i+1,imol) - 1
        end do
      end if
c
ccc      do i=1,nchain(imol)
ccc        call jvalut (' Chain ptrs :',2,chnptr(1,i,imol))
ccc      end do
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine molut (imol,ierr,chain,i1,i2)
c
      include 'lsqman.incl'
c
      integer imol,ierr,kk,imod,nwr,leng1,i1,i2
c
      logical xinter,lhydro
c
      character key*6, line*128, prev*1, chain*1
c
code ...
c
      ierr = 0
c
      call xopxua (iunit,file(imol),xinter(),ierr)
      if (ierr .ne. 0) return
c
c ... write the PDB file
c
 4711 format (a6,a)
 4713 format (i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2,6x,a4)
c
      key = 'ATOM  '
      prev = '?'
      imod = 0
      nwr = 0
c
      do kk=1,natoms(imol)
c
        if (mmodel(imol)) then
          if (achain(kk,imol) .ne. prev) then
            imod = imod + 1
            prev = achain(kk,imol)
            write (iunit,'(a5,i10)',err=997) 'MODEL',imod
          end if
        end if
c
        if (.not. lhkeep) then
          if (lhydro(atmnam(kk,imol))) goto 1938
        end if
c
        if (chain .ne. '*') then
          if (chain .ne. achain(kk,imol)) goto 1938
        end if
c
c ... check residue range (if defined)
c
        if (i1 .gt. -9999) then
          if (iresid(kk,imol) .lt. i1) goto 1938
        end if
c
        if (i2 .gt. -9999) then
          if (iresid(kk,imol) .gt. i2) goto 1938
        end if
c
        write (line,4713,err=5372) atomnr(kk,imol),
     +    atmnam(kk,imol),resnam(kk,imol),achain(kk,imol),
     +    iresid(kk,imol),atmxyz(1,kk,imol),atmxyz(2,kk,imol),
     +    atmxyz(3,kk,imol),qatom(kk,imol),batom(kk,imol),
     +    axplor(kk,imol)
        write (iunit,4711,err=997) key,line(1:leng1(line))
        nwr = nwr + 1
c
 1938   continue
      end do
c
      write (iunit,'(a)',err=997) 'END'
      call jvalut (' Number of atoms written :',1,nwr)
c
      if (mmodel(imol)) then
        call jvalut (' Number of NMR models    :',1,imod)
      end if
c
      close (iunit)
      return
c
c ... write to file error
c
  997 continue
      close (iunit)
      call errcon ('While writing PDB file')
      call jvalut (' Atom nr :',1,kk)
      call textut (' Line :',line)
      ierr = -1
      return
c
c ... write to string error
c
 5372 continue
      close (iunit)
      call errcon ('While writing PDB file')
      call jvalut (' Atom nr :',1,kk)
      call textut (' Line :',line)
      ierr = -1
      return
c
      end
c
c
c
      subroutine lsqexp (imol,irange,jmol,jrange,mode,lprint,ierr)
c
      include 'lsqman.incl'
c
      integer maxopt
      parameter (maxopt = 100)
c
      real rtx(12),vec1(3),vec2(3)
      real alpha,rms,x1,x2,x3,x4,x5,x6,x7,x8,x9,usm,dsq,cr,rno
      real qd,qdmax,qdrms,distce,qqq,d0,dd,tmscor
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j,k,kk,kptr
      integer idum,jdum,kdum,ldum,nuse,iptr,jptr,ires,jres
      integer mode,what,nbskip,leng1,iii,nrms,mymode
c
      logical lhydro,lprint,mainch
c
      character irange*(*),jrange*(*),ichn*1,jchn*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20
c
code ...
c
      ierr = 0
c
      what = -1
      if (mode .gt. 0) what = 1
      mymode = mode
c
      if (what .eq. 1) then
        if (lprint) write (*,6000) name(imol)(1:leng1(name(imol))),
     +    irange(1:leng1(irange)),
     +    name(jmol)(1:leng1(name(jmol))),
     +    jrange(1:leng1(jrange)),
     +    (atypes(i),i=1,natype)
      else
        if (lprint) write (*,6001) name(imol)(1:leng1(name(imol))),
     +    irange(1:leng1(irange)),
     +    name(jmol)(1:leng1(name(jmol))),
     +    jrange(1:leng1(jrange)),
     +    (atypes(i),i=1,natype)
      end if
c
      if (lprint) write (*,6002) bcutlo,bcuthi
c
 6000 format (' Explicit fit of ',a,1x,a/
     +        ' And             ',a,1x,a/
     +        ' Atom types     |',15(a4,'|'))
 6001 format (' Calculate RMSD of ',a,1x,a/
     +        ' And               ',a,1x,a/
     +        ' Atom types       |',15(a4,'|'))
 6002 format (' B-factor range used: ',f8.2,' - ',f8.2,' A2')
c
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        if (lprint) call errcon ('In range1')
        goto 9900
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        if (lprint) call errcon ('In range2')
        goto 9900
      end if
c
      if (nopt1 .ne. nopt2) then
        if (lprint) call errcon ('Different nr of zones')
        goto 9900
      end if
c
      nuse = 0
      nbskip = 0
      iptr = 0
      jptr = 0
c
      do i=1,nopt1
c
        ichn = optpa1(i)(1:1)
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          if (lprint) call errcon (
ccc     +      'Invalid chain id in range1 : '//optpa1(i))
ccc          goto 9900
        end if
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
c
c ... 971111 - allow e.g. "A78" to mean "A78-78"
c
           k = length(optpa1(i))
           if (ichn .eq. ' ') then
             optpa1(i) = optpa1(i)(1:k) // '-' // optpa1(i)(1:k)
           else
             optpa1(i) = optpa1(i)(1:k) // '-' // optpa1(i)(2:k)
           end if
           j = k + 1
ccc          if (lprint) call errcon ('Not a zone in range1 : '//optpa1(i))
ccc          goto 9900
        end if
ccc      print *,optpa1(i)(2:j-1)
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
ccc      print *,optpa1(i)(j+1:)
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          if (lprint) call errcon (
     +      'Invalid zone in range1 : '//optpa1(i))
          goto 9900
        end if
c
        jchn = optpa2(i)(1:1)
        if (index('1234567890',jchn) .gt. 0) then
          jchn = ' '
ccc        else if (jchn .lt. 'A' .or. jchn .gt. 'Z') then
ccc          if (lprint) call errcon (
ccc     +      'Invalid chain id in range2 : '//optpa2(i))
ccc          goto 9900
        end if
        j = index (optpa2(i),'-')
        if (j .le. 0) j = index (optpa2(i),':')
        if (j .gt. 0) then
          if (lprint) call errcon (
     +      'Zone ignored in range2 : '//optpa2(i))
          optpa2(i) = optpa2(i)(1:j-1)
        end if
ccc      print *,optpa2(i)(2:)
        if (jchn .eq. ' ') then
          call str2i (optpa2(i)(1:),kdum,ierr)
        else
          call str2i (optpa2(i)(2:),kdum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        ldum = kdum + (jdum - idum)
c
        do ires=idum,jdum
          jres = kdum + (ires - idum)
          do j=1,natype
ccc      print *,'Fetch ',imol,ichn,ires,atypes(j)
ccc      print *,'and   ',jmol,jchn,jres
            call getptr (imol,ichn,ires,atypes(j),
     +                   jmol,jchn,jres,iptr,jptr,ierr)
ccc      print *,'Err/i/j ',ierr,iptr,jptr
            if (ierr .eq. 0) then
c
              nuse = nuse + 1
              do k=1,3
                buffi(k,nuse) = atmxyz(k,iptr,imol)
                buffj(k,nuse) = atmxyz(k,jptr,jmol)
ccc      print *,nuse,buffi(k,nuse),buffj(k,nuse)
              end do
c
c ... calc rms delta B
c
              buffb1(nuse) = batom(iptr,imol)
              buffb2(nuse) = batom(jptr,jmol)
c
c ... reject if one of atoms outside B-factor limits
c
              if (batom(iptr,imol) .lt. bcutlo .or.
     +            batom(iptr,imol) .gt. bcuthi .or.
     +            batom(jptr,jmol) .lt. bcutlo .or.
     +            batom(jptr,jmol) .gt. bcuthi) then
                nbskip = nbskip + 1
                nuse = nuse - 1
                goto 8001
              end if
c
c ... get others if ALL or NONH or SIDE or TRAC
c
              if ( natype .eq. 1 .and.
     +            (atypes(1).eq.'ALL ' .or.
     +             atypes(1).eq.'NONH' .or.
     +             atypes(1).eq.'SIDE' .or.
     +             atypes(1).eq.'TRAC') ) then
                do k=iptr+1,natoms(imol)
                  if (iresid(k,imol).ne.iresid(iptr,imol) .or.
     +                achain(k,imol).ne.achain(iptr,imol)) then
                    kptr = k - 1
                    goto 6511
                  else if (k .eq. natoms(imol)) then
                    kptr = k
                    goto 6511
                  end if
                end do
                kptr = iptr
 6511           continue
c
ccc      print *,'IRES,IPTR,KPTR,NUSE ',IRES,IPTR,KPTR,NUSE
c
                if (kptr .gt. iptr) then
                  do k=iptr+1,kptr
c
c ... 960729 - find same atom in other molecule
c
                    call getptr (imol,ichn,ires,atmnam(k,imol),
     +                           jmol,jchn,jres,iii,jptr,ierr)
c
                    if (ierr .ne. 0) goto 6436
c
ccc                    jptr = jptr + 1
ccc                    if ( atmnam(jptr,jmol) .eq.
ccc     +                   atmnam(k,imol) ) then
                      if (
     +   (atypes(1) .eq. 'ALL ')
     +      .or.
     +   (atypes(1) .eq. 'NONH' .and. (.not. lhydro(atmnam(k,imol))))
     +      .or.
     +   (atypes(1) .eq. 'SIDE' .and. (.not. mainch(atmnam(k,imol)))
     +                          .and. (.not. lhydro(atmnam(k,imol))))
     +      .or.
     +   (atypes(1) .eq. 'TRAC' .and. (.not. lhydro(atmnam(k,imol)))
     +       .and. (atmnam(k,imol) .eq. ' CA ' .or.
     +              (.not. mainch(atmnam(k,imol))) ))
     +      ) then
ccc     +                    atmnam(k,imol)(2:2) .ne. 'H') then
                        nuse = nuse + 1
                        do kk=1,3
                          buffi(kk,nuse) = atmxyz(kk,k,imol)
                          buffj(kk,nuse) = atmxyz(kk,jptr,jmol)
                        end do
c
c ... calc rms delta B
c
                        buffb1(nuse) = batom(k,imol)
                        buffb2(nuse) = batom(jptr,jmol)
c
                      end if
ccc                    end if
 6436               continue
                  end do
                end if
              end if
c
 8001         continue
            end if
          end do
        end do
c
      end do
c
      if (lprint) call jvalut (' Nr of atoms to match  :',1,nuse)
      if (lprint) call jvalut (' Nr skipped (B limits) :',1,nbskip)
      if (nuse .lt. 3) then
        if (lprint) call errcon ('Fewer than 3 atoms; cannot match')
        goto 9900
      end if
c
      if (what .eq. 1) then
        call lsqgjk (buffi,buffj,nuse,rms,rtx,ierr)
        if (ierr .ne. 0) goto 9900
      else
        do i=1,12
          rtx (i) = rtlsq (i,imol,jmol)
        end do
        call vecrtv (buffj,buffk,nuse,rtx(1),rtx(10))
c
        if (mymode .ge. 0) mymode = -1
        if (mymode .le. -25) then
          d0 = ( 1.24 * (float(-mymode)-15.0)**0.333 ) - 1.8
        else
          d0 = 1.0
        end if
ccc        print *,' *** d0 = ',d0
ccc        print *,' *** lt = ',-mymode
        d0 = d0*d0
        tmscor = 0.0
        rms = 0.0
c
        do i=1,nuse
          dd = 0.0
          do j=1,3
            dd = dd + (buffi(j,i)-buffk(j,i))**2
          end do
          rms = rms + dd
          tmscor = tmscor + (1.0 / (1.0 + (dd/d0)))
ccc          print *,' i,dd,contrib ',i,dd,(1.0 / (1.0 + (dd/d0)))
        end do
        tmscor = tmscor / float(-mymode)
        rms = sqrt ( rms / float(nuse) )
c
c ... analyse difference-distance matrix
c
        if (lprint) then
          qdmax = 0.0
          qdrms = 0.0
          nrms = 0
          do i=1,nuse-1
            do j=i+1,nuse
              qd = distce (buffi(1,i),buffi(1,j)) -
     +             distce (buffk(1,i),buffk(1,j))
              if (abs(qd) .gt. qdmax) qdmax = abs(qd)
              nrms = nrms + 1
              qdrms = qdrms + qd*qd
            end do
          end do
          qdrms = sqrt (qdrms / float(nrms))
        end if
c
      end if
c
c ... calc rms delta B
c
      call xystat (buffb1,buffb2,nuse,x1,x2,x3,x4,x5,x6,x7,x8,x9)
c
      if (lprint) write (*,6010) nuse,rms,x1,x3,
     +  rtx(1),rtx(4),rtx(7),rtx(2),rtx(5),rtx(8),
     +  rtx(3),rtx(6),rtx(9),rtx(10),rtx(11),rtx(12)
c
 6010 format (/
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' RMS delta B  = ',f8.3,' A2'/
     +  ' Corr. coeff. = ',f11.4/
     +  ' Rotation    : ',3f10.6/15x,3f10.6/15x,3f10.6/
     +  ' Translation : ',3f10.3)
c
 6013 format (
     +  ' Maiorov-Crippen RHO (0-2)            = ',f12.5/
     +  ' Estimated RMSD for 2 random proteins = ',f10.3,' A'/
     +  ' Relative RMSD                        = ',f12.5/
     +  ' Normalised RMSD (100)                = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A')
c
      if (what .eq. 1) then
c
        call mcrho (nuse,buffi,buffj,rms,
     +    cripp(imol,jmol),rrmsd(imol,jmol),normsd(imol,jmol),dsq)
c
        rmsdna (imol,jmol) = rms / float(nuse)
c
        if (lprint)
     +    write (*,6013) cripp(imol,jmol),dsq,rrmsd(imol,jmol),
     +                   normsd(imol,jmol),rmsdna(imol,jmol)
c
c ... store results
c
        do i=1,12
          rtlsq (i,imol,jmol) = rtx (i)
        end do
c
        nmatch (imol,jmol) = nuse
        rmsd (imol,jmol) = rms
        rmsb (imol,jmol) = x1
        corb (imol,jmol) = x3
        simind (imol,jmol) = rms
        matchi (imol,jmol) = 1.0 / (1.0+rms)
        qqq = 100.0/float(nuse)
        sas1 (imol,jmol) = rms * qqq
        sas2 (imol,jmol) = sas1 (imol,jmol) * qqq
        sas3 (imol,jmol) = sas2 (imol,jmol) * qqq
        sas4 (imol,jmol) = sas3 (imol,jmol) * qqq
c
      else
c
        if (lprint) then
c
          call mcrho (nuse,buffi,buffj,rms,cr,usm,rno,dsq)
          write (*,6013) cr,dsq,usm,rno,rms/float(nuse)
c
          if (mymode .le. -25) then
            write (*,*)
            call jvalut (
     +        ' Value of Ltarget               :',1,-mymode)
            call fvalut (
     +        ' Value of TM-score              :',1,tmscor)
            if (nuse .gt. -mymode) then
              call errcon (' Nr of atoms compared > Ltarget !!')
            end if
          end if
c
          write (*,*)
          call jvalut (
     +      ' Nr of unique elements in DDM   :',1,nrms)
          call fvalut (
     +      ' Max absolute DDM element (A)   :',1,qdmax)
          call fvalut (
     +      ' RMS of unique DDM elements (A) :',1,qdrms)
c
        end if
c
c ... calculate angle between vectors from first to last atom
c
        do i=1,3
          vec1 (i) = buffi(i,nuse)-buffi(i,1)
          vec2 (i) = buffk(i,nuse)-buffk(i,1)
        end do
c
        call vecang (vec1,vec2,alpha,ierr)
        if (ierr .eq. 0) then
          write (*,*)
          call prompt (
     +      ' Vectors between first and last selected atoms:')
          call fvalut (' Mol 1  :',3,vec1)
          call fvalut (' Mol 2" :',3,vec2)
          call fvalut (' Angle between them (deg) :',1,alpha)
        end if
c
      end if
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine getptr (imol,ichn,ires,atype,jmol,jchn,jres,
     +                   iptr,jptr,ierr)
c
      include 'lsqman.incl'
c
      integer ires,jres,iptr,jptr,imol,jmol,iold,jold,iok,jok
      integer i,ierr
c
      logical lhydro,mainch
c
      character ichn*1,jchn*1,atype*4
c
code ...
c
ccc      print *,'GETPTR ',imol,jmol
c
      if (imol .gt. 0) then
        iold = max (1, min (iptr-50,natoms(imol)-1))
        if (atype .eq. 'ALL ' .or. atype .eq. 'NONH' .or.
     +      atype .eq. 'SIDE' .or. atype .eq. 'TRAC') iold = 1
      else
        iold = 0
      end if
c
      if (jmol .gt. 0) then
        jold = max (1, min (jptr-50,natoms(jmol)-1))
        if (atype .eq. 'ALL ' .or. atype .eq. 'NONH' .or.
     +      atype .eq. 'SIDE' .or. atype .eq. 'TRAC') jold = 1
      else
        jold = 0
      end if
c
      iok = 0
      jok = 0
      ierr = 0
c
      if (imol .le. 0 .and. jmol .le. 0) then
        ierr = -1
        return
      end if
c
      if (imol .le. 0) goto 1000
c
ccc      print *,iold,jold,iptr,jptr
c
      do i=iold,natoms(imol)
        if (achain(i,imol) .eq. ichn) then
          if (iresid(i,imol) .eq. ires) then
            if (atype .eq. 'ALL ') then
              iok = i
              goto 1000
            else if (atype .eq. 'NONH') then
ccc              if (atmnam(i,imol)(2:2) .ne. 'H') then
              if (.not. lhydro (atmnam(i,imol))) then
                iok = i
                goto 1000
              end if
            else if (atype .eq. 'SIDE') then
              if (.not. mainch (atmnam(i,imol))) then
                iok = i
                goto 1000
              end if
            else if (atype .eq. 'TRAC') then
              if (.not. lhydro (atmnam(i,imol))) then
                if (atmnam(i,imol) .eq. ' CA ') then
                  iok = i
                  goto 1000
                else if (.not. mainch (atmnam(i,imol))) then
                  iok = i
                  goto 1000
                end if
              end if
            else if (atmnam(i,imol) .eq. atype) then
              iok = i
              goto 1000
            end if
          end if
        end if
      end do
c
      do i=1,iold
        if (achain(i,imol) .eq. ichn) then
          if (iresid(i,imol) .eq. ires) then
            if (atype .eq. 'ALL ') then
              iok = i
              goto 1000
            else if (atype .eq. 'NONH') then
ccc              if (atmnam(i,imol)(2:2) .ne. 'H') then
              if (.not. lhydro (atmnam(i,imol))) then
                iok = i
                goto 1000
              end if
            else if (atype .eq. 'SIDE') then
              if (.not. mainch (atmnam(i,imol))) then
                iok = i
                goto 1000
              end if
            else if (atype .eq. 'TRAC') then
              if (.not. lhydro (atmnam(i,imol))) then
                if (atmnam(i,imol) .eq. ' CA ') then
                  iok = i
                  goto 1000
                else if (.not. mainch (atmnam(i,imol))) then
                  iok = i
                  goto 1000
                end if
              end if
            else if (atmnam(i,imol) .eq. atype) then
              iok = i
              goto 1000
            end if
          end if
        end if
      end do
c
ccc      print *,'I not found'
      ierr = -1
      return
c
 1000 continue
c
      if (jmol .le. 0) goto 2000
c
      do i=jold,natoms(jmol)
        if (achain(i,jmol) .eq. jchn) then
          if (iresid(i,jmol) .eq. jres) then
            if (atype .eq. 'ALL ') then
              jok = i
              goto 2000
            else if (atype .eq. 'NONH') then
ccc              if (atmnam(i,jmol)(2:2) .ne. 'H') then
              if (.not. lhydro (atmnam(i,jmol))) then
                jok = i
                goto 2000
              end if
            else if (atype .eq. 'SIDE') then
              if (.not. mainch (atmnam(i,jmol))) then
                jok = i
                goto 2000
              end if
            else if (atype .eq. 'TRAC') then
              if (.not. lhydro (atmnam(i,jmol))) then
                if (atmnam(i,jmol) .eq. ' CA ') then
                  jok = i
                  goto 2000
                else if (.not. mainch (atmnam(i,jmol))) then
                  jok = i
                  goto 2000
                end if
              end if
            else if (atmnam(i,jmol) .eq. atype) then
              jok = i
              goto 2000
            end if
          end if
        end if
      end do
c
      do i=1,jold
        if (achain(i,jmol) .eq. jchn) then
          if (iresid(i,jmol) .eq. jres) then
            if (atype .eq. 'ALL ') then
              jok = i
              goto 2000
            else if (atype .eq. 'NONH') then
ccc              if (atmnam(i,jmol)(2:2) .ne. 'H') then
              if (.not. lhydro (atmnam(i,jmol))) then
                jok = i
                goto 2000
              end if
            else if (atype .eq. 'SIDE') then
              if (.not. mainch (atmnam(i,jmol))) then
                jok = i
                goto 2000
              end if
            else if (atype .eq. 'TRAC') then
              if (.not. lhydro (atmnam(i,jmol))) then
                if (atmnam(i,jmol) .eq. ' CA ') then
                  jok = i
                  goto 2000
                else if (.not. mainch (atmnam(i,jmol))) then
                  jok = i
                  goto 2000
                end if
              end if
            else if (atmnam(i,jmol) .eq. atype) then
              jok = i
              goto 2000
            end if
          end if
        end if
      end do
c
ccc      print *,'J not found'
      ierr = -1
      return
c
 2000 continue
      iptr = iok
      jptr = jok
      ierr = 0
c
ccc      print *,imol,ichn,ires,atype,jmol,jchn,jres,iptr,jptr,ierr
c
      return
      end
c
c
c
      subroutine selimp (imol,irange,jmol,jrange,lprint,ierr)
c
      include 'lsqman.incl'
c
      integer maxopt
      parameter (maxopt = 25)
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j,k
      integer idum,jdum,iptr,jptr,ires,jres
c
      logical lprint
c
      character irange*(*),jrange*(*),ichn*1,jchn*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20
c
code ...
c
      ierr = 0
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range1')
        goto 9900
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range2')
        goto 9900
      end if
c
      nuse1 = 0
      iptr = 0
c
      do i=1,natoms(imol)
        usei(i) = .false.
      end do
c
      do i=1,nopt1
c
        ichn = optpa1(i)(1:1)
c
c ... ALL MOLECULE
c
        if (ichn .eq. '*') then
          do j=1,natoms(imol)
            if (atmnam(j,imol) .eq. atypes(1)) then
              if (.not. usei(j)) then
                nuse1 = nuse1 + 1
                usei(j) = .true.
              end if
            end if
          end do
          goto 1000
        end if
c
c ... A WHOLE CHAIN
c
        if (optpa1(i)(2:2) .eq. '*' .or.
     +      length(optpa1(i)) .eq. 1) then
          do j=1,natoms(imol)
            if (achain(j,imol) .eq. ichn) then
              if (atmnam(j,imol) .eq. atypes(1)) then
                if (.not. usei(j)) then
                  nuse1 = nuse1 + 1
                  usei(j) = .true.
                end if
              end if
            end if
          end do
          goto 990
        end if
c
c ... A ZONE
c
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc          goto 1234
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range1 : '//optpa1(i))
ccc          goto 9900
        end if
c
ccc 1234   continue
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
          call errcon ('Not a zone in range1 : '//optpa1(i))
          goto 9900
        end if
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon ('Invalid zone in range1 : '//optpa1(i))
          goto 9900
        end if
c
        do ires=idum,jdum
          call getptr (imol,ichn,ires,atypes(1),
     +                 0,jchn,jres,iptr,jptr,ierr)
          if (ierr .eq. 0) then
            if (.not. usei(iptr)) then
              nuse1 = nuse1 + 1
              usei (iptr) = .true.
            end if
          end if
        end do
c
 990    continue
c
      end do
c
 1000 continue
c
      nuse1 = 0
      do i=1,natoms(imol)
        if (usei(i)) then
          nuse1 = nuse1 + 1
          ptri(nuse1) = i
          do k=1,3
            buffi(k,nuse1) = atmxyz(k,i,imol)
          end do
c
          buffb1 (nuse1) = batom (i,imol)
c
        end if
      end do
      if (lprint) call jvalut (' Nr of atoms in mol1 :',1,nuse1)
c
      if (nuse1 .lt. 3) then
        call errcon ('Fewer than 3 atoms - aborting')
        ierr = -1
        return
      end if
c
      do i=1,natoms(jmol)
        usej(i) = .false.
      end do
c
      nuse2 = 0
      iptr = 0
c
      do i=1,nopt2
c
        ichn = optpa2(i)(1:1)
c
c ... ALL MOLECULE
c
        if (ichn .eq. '*') then
          do j=1,natoms(jmol)
            if (atmnam(j,jmol) .eq. atypes(1)) then
              if (.not. usej(j)) then
                nuse2 = nuse2 + 1
                usej(j) = .true.
              end if
            end if
          end do
          goto 2000
        end if
c
c ... A WHOLE CHAIN
c
        if (optpa2(i)(2:2) .eq. '*' .or.
     +      length(optpa2(i)) .eq. 1) then
          do j=1,natoms(jmol)
            if (achain(j,jmol) .eq. ichn) then
              if (atmnam(j,jmol) .eq. atypes(1)) then
                if (.not. usej(j)) then
                  nuse2 = nuse2 + 1
                  usej(j) = .true.
                end if
              end if
            end if
          end do
          goto 1990
        end if
c
c ... A ZONE
c
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc          goto 2345
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range2 : '//optpa2(i))
ccc          goto 9900
        end if
c
ccc 2345   continue
        j = index (optpa2(i),'-')
        if (j .le. 1) j = index (optpa2(i),':')
        if (j .le. 1) then
          call errcon ('Not a zone in range2 : '//optpa2(i))
          goto 9900
        end if
        if (ichn .eq. ' ') then
          call str2i (optpa2(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa2(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        if (optpa2(i)(j+1:j+1) .eq. ichn)
     +    optpa2(i)(j+1:j+1) = ' '
        call str2i (optpa2(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon ('Invalid zone in range2 : '//optpa2(i))
          goto 9900
        end if
c
        do ires=idum,jdum
          call getptr (jmol,ichn,ires,atypes(1),
     +                 0,'*',0,iptr,jptr,ierr)
          if (ierr .eq. 0) then
            if (.not. usej(iptr)) then
              nuse2 = nuse2 + 1
              usej (iptr) = .true.
            end if
          end if
        end do
c
 1990   continue
c
      end do
c
 2000 continue
c
      nuse2 = 0
      do i=1,natoms(jmol)
        if (usej(i)) then
          nuse2 = nuse2 + 1
          ptrj(nuse2) = i
          do k=1,3
            buffj(k,nuse2) = atmxyz(k,i,jmol)
          end do
c
          buffb2 (nuse2) = batom (i,jmol)
c
        end if
      end do
      if (lprint) call jvalut (' Nr of atoms in mol2 :',1,nuse2)
c
      if (nuse2 .lt. 3) then
        call errcon ('Fewer than 3 atoms - aborting')
        ierr = -1
        return
      end if
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine lsqimp (imol,irange,jmol,jrange,lprint,ierr)
c
      include 'lsqman.incl'
c
      real rtx(12),rtnew(12)
      real rms,dmax,dmin,dist,si,ddec,mi,best,xdist,cr,usm,dsq
      real x1,x2,x3,x4,x5,x6,x7,x8,x9,rms0,rmsm1,rmsp1,dvalue
      real rno,rna,s1,s2,s3,s4,qqq
c
      integer imol,jmol,ierr,i,j,k,l,jmin,uselen,nsame,nali
      integer ncycle,nmat,imin,jold,kold,ilast,jlast,nres,leng1
c
      logical lprint
c
      character irange*(*),jrange*(*)
      character jchn*1,kchn*1,symbol*1
c
code ...
c
      ierr = 0
c
      if (lprint) write (*,6000) name(imol)(1:leng1(name(imol))),
     +  irange(1:leng1(irange)),
     +  name(jmol)(1:leng1(name(jmol))),
     +  jrange(1:leng1(jrange)),
     +  atypes(1)
c
 6000 format (' Improve fit of  ',a,1x,a/
     +        ' And             ',a,1x,a/
     +        ' Atom type      |',a4,'|')
c
      if (atypes(1) .eq. 'ALL ' .or. atypes(1) .eq. 'NONH' .or.
     +    atypes(1) .eq. 'SIDE' .or. atypes(1) .eq. 'TRAC') then
        call errcon ('Invalid atom type !')
        ierr = -1
        return
      end if
c
      call selimp (imol,irange,jmol,jrange,lprint,ierr)
      if (ierr .ne. 0) return
c
ccc      print *,'LSQIMP after SELIMP ',imol,jmol
c
      do i=1,12
        rtx (i) = rtlsq (i,imol,jmol)
      end do
c
      ncycle = 0
c
      dmax = dismax ** 2
      ddec = decay ** 2
c
      uselen = minlen
c
      if (optcri .eq. 'SI') then
        best = 999999.99
      else if (optcri .eq. 'MI') then
        best = -1.0
      else if (optcri .eq. 'S1') then
        best = 999999.99
      else if (optcri .eq. 'S2') then
        best = 999999.99
      else if (optcri .eq. 'S3') then
        best = 999999.99
      else if (optcri .eq. 'S4') then
        best = 999999.99
      else if (optcri .eq. 'RM') then
        best = 999999.99
      else if (optcri .eq. 'CR') then
        best = 999999.99
      else if (optcri .eq. 'RR') then
        best = 999999.99
      else if (optcri .eq. 'NR') then
        best = 999999.99
      else if (optcri .eq. 'NM') then
        best = -1.0
      else
        if (lprint) call errcon (
     +    'Invalid optimisation criterion : '//optcri)
        ierr = -1
        return
      end if
c
c ... major loop
c
   10 continue
      ncycle = ncycle + 1
c
c ... apply current operator
c
      call vecrtv (buffj,buffk,nuse2,rtx(1),rtx(10))
c
ccc      print *,'BUFFJ ',((buffj(i,j),i=1,3),j=1,3)
ccc      print *,'BUFFK ',((buffk(i,j),i=1,3),j=1,3)
ccc      print *,'RT ',rtx
ccc      print *,'NUSE2 ',nuse2
c
      do i=1,nuse2
        usej(i) = .false.
      end do
c
      do i=1,nuse1
        usei (i) = .false.
        ptrij (i) = 0
      end do
c
      nmat = 0
c
      if (uselen .le. 1) then
c
        jmin = 1
c
c ... CASE 1 - min fragment length <= 1
c
        do i=1,nuse1
c
          dmin = dmax
          imin = -1
c
          do j=jmin,nuse2
            if (.not. usej(j)) then
              dist = (buffi(1,i)-buffk(1,j))**2
              if (dist .lt. dmin) then
                dist = dist + (buffi(2,i)-buffk(2,j))**2
                if (dist .lt. dmin) then
                  dist = dist + (buffi(3,i)-buffk(3,j))**2
                  if (dist .lt. dmin) then
                    dmin = dist
                    imin = j
                  end if
                end if
              end if
            end if
          end do
c
          if (imin .gt. 0) then
            nmat = nmat + 1
            usei (i) = .true.
            ptrij (i) = imin
            if (ncycle .eq. 1) sptrij (i) = imin
            usej (imin) = .true.
            do k=1,3
              buffl (k,nmat) = buffi (k,i)
              buffm (k,nmat) = buffj (k,imin)
            end do
c
            buffb3 (nmat) = buffb1 (i)
            buffb4 (nmat) = buffb2 (imin)
c
            if (seqcri .eq. 'ON') jmin = imin + 1
c
          end if
c
        end do
c
      else
c
c ... CASE 2 - min fragment length > 1
c
        if (lprint) write (*,*)
c
        i = 1
        jmin = 1
c
   20   continue
c
          dmin = dmax
          imin = -1
c
ccc      print *,'JMIN,NUSE2 ',jmin,nuse2
c
          do j=jmin,nuse2
            if (.not. usej(j)) then
              dist = (buffi(1,i)-buffk(1,j))**2
              if (dist .lt. dmin) then
                dist = dist + (buffi(2,i)-buffk(2,j))**2
                if (dist .lt. dmin) then
                  dist = dist + (buffi(3,i)-buffk(3,j))**2
                  if (dist .lt. dmin) then
                    dmin = dist
                    imin = j
                  end if
                end if
              end if
            end if
          end do
c
ccc      print *,'I/IMIN ',i,imin,dmin
c
          if (imin .gt. 0) then
c
 1941       continue
c
ccc      print *,'PTRI/J ',ptri(i),ptrj(imin)
c
            jold = iresid(ptri(i),imol)
            kold = iresid(ptrj(imin),jmol)
            jchn = achain(ptri(i),imol)
            kchn = achain(ptrj(imin),jmol)
            ilast = i
            jlast = imin
c
ccc      print *,'IOLD  ETC. ',jold,kold,jchn,kchn
ccc      print *,'ILAST ETC. ',ilast,ptri(ilast+1),
ccc     +                      jlast,ptrj(jlast+1)
c
            do l=i+1,nuse1
c
c ... residues left ?
c
              if (jlast .ge. nuse2) goto 1999
c
              if (usej(jlast+1)) goto 1999
c
c ... same chains ?
c
              if ( jchn .ne. achain(ptri(ilast+1),imol)     ) goto 1999
              if ( kchn .ne. achain(ptrj(jlast+1),jmol)     ) goto 1999
c
ccc      print *,'INOW ',l,' CHAIN OKAY'
c
c ... consecutive residues ?
c
c
ccc      print *,jold+1,iresid(ptri(ilast+1),imol),
ccc     +        kold+1,iresid(ptrj(jlast+1),imol)
c
              if ( (jold+1) .ne. iresid(ptri(ilast+1),imol) ) goto 1999
              if ( (kold+1) .ne. iresid(ptrj(jlast+1),jmol) ) goto 1999
c
ccc      print *,'INOW ',l,' RESIDUE NUMBER OKAY'
c
c ... distance okay ?
c
              xdist = (buffi(1,ilast+1)-buffk(1,jlast+1)) ** 2 +
     +                (buffi(2,ilast+1)-buffk(2,jlast+1)) ** 2 +
     +                (buffi(3,ilast+1)-buffk(3,jlast+1)) ** 2
c
ccc      print *,'INOW ',l,' DISTANCE ',xdist
c
              if (xdist .gt. dmax) goto 1999
c
              jold = iresid(ptri(ilast+1),imol)
              kold = iresid(ptrj(jlast+1),jmol)
              jchn = achain(ptri(ilast+1),imol)
              kchn = achain(ptrj(jlast+1),jmol)
c
              ilast = ilast + 1
              jlast = jlast + 1
c
            end do
c
 1999       continue
            nres = (ilast-i+1)
c
ccc      print *,'NRES,MINLEN ',nres,uselen
c
c ... try frameshifts ?
c
c ... goto 1943 or goto 1941 ???
c     empirical test with 1CEL on to 1AYH shows that 1943 is better
c     IF you use 1941, use "3" instead of "USELEN" since fragments
c     may grow later
c
ccc            if (lshift .and. nres .ge. 3) then
            if (lshift .and. nres .ge. uselen) then
ccc              print *,'I, IMIN before :',i,imin
 1943         continue
c
              call lsqrms (nres,buffi(1,i),
     +                     buffk(1,imin),rms0)
c
              rmsm1 = rms0 + 1.0
              if (imin .gt. 1) then
                if (.not. usej(imin-1)) then
                if (kchn .eq. achain(ptrj(imin-1),jmol)) then
                if (iresid(ptrj(imin-1),jmol)+1  .eq.
     +              iresid(ptrj(imin),jmol) ) then
                  call lsqrms (nres,buffi(1,i),
     +                         buffk(1,imin-1),rmsm1)
                end if
                end if
                end if
              end if
c
              rmsp1 = rms0 + 1.0
              if ((imin+nres) .lt. nuse2) then
                if (.not. usej(imin+nres)) then
                if (kchn .eq. achain(ptrj(imin+nres),jmol)) then
                if (iresid(ptrj(imin+nres-1),jmol)+1  .eq.
     +              iresid(ptrj(imin+nres),jmol) ) then
                  call lsqrms (nres,buffi(1,i),
     +                         buffk(1,imin+1),rmsp1)
                end if
                end if
                end if
              end if
c
ccc              print *,'RMS 0/-1/+1 = ',rms0,rmsm1,rmsp1
              if (rmsm1 .lt. rms0 .and. rmsm1 .lt. rmsp1) then
                imin = imin - 1
                if (lprint) write (*,6935) '-',rmsm1,rms0,rmsp1
                goto 1943
              end if
c
              if (rmsp1 .lt. rms0 .and. rmsp1 .lt. rmsm1) then
                imin = imin + 1
                if (lprint) write (*,6935) '+',rmsm1,rms0,rmsp1
                goto 1943
              end if
ccc              print *,'I, IMIN after  :',i,imin
            end if
 6935 format (' Frameshift ',a1,'1 ; RMSD -1/0/+1 (A) :',3f8.3)
c
            if (nres .ge. uselen) then
c
c ... okay, fragment long enough
c
              if (lprint) call ivalut (
     +          ' Found fragment of length :',1,nres)
c
              do l=i,ilast
                nmat = nmat + 1
                usei (i) = .true.
                ptrij (i) = imin
                if (ncycle .eq. 1) sptrij (i) = imin
                usej (imin) = .true.
                do k=1,3
                  buffl (k,nmat) = buffi (k,i)
                  buffm (k,nmat) = buffj (k,imin)
                end do
c
                buffb3 (nmat) = buffb1 (i)
                buffb4 (nmat) = buffb2 (imin)
c
                imin = imin + 1
                i = i + 1
              end do
c
              if (seqcri .eq. 'ON') jmin = imin
c
            else
              i = i + 1
            end if
          else
            i = i + 1
          end if
c
        if (i .le. nuse1) goto 20
c
      end if
c
      if (lprint) then
        write (*,*)
        call jvalut (' Cycle :',1,ncycle)
        call fvalut (' Distance cut-off (A)      :',1,sqrt(dmax))
        call ivalut (' Min fragment length (res) :',1,uselen)
      end if
c
      if (nmat .lt. 3) then
        if (lprint) call jvalut (' Nr of matching residues :',1,nmat)
        if (lprint) call errcon ('Fewer than 3 atoms; cannot match')
        ierr = -1
        return
      end if
c
      call lsqgjk (buffl,buffm,nmat,rms,rtnew,ierr)
      if (ierr .ne. 0) return
c
      call xystat (buffb3,buffb4,nmat,x1,x2,x3,x4,x5,x6,x7,x8,x9)
c
      si = rms * float ( min (nuse1,nuse2) ) / float ( nmat )
      mi = (1.0+float(nmat)) /
     +     ( (1.0+rmswgt*rms) * (1.0+float(min (nuse1,nuse2))) )
      call mcrho (nmat,buffl,buffm,rms,cr,usm,rno,dsq)
      rna = rms / float (nmat)
      qqq = 100.0/float(nmat)
      s1 = rms * qqq
      s2 = s1 * qqq
      s3 = s2 * qqq
      s4 = s3 * qqq
c
      if (lprint) write (*,6010) nmat,rms,si,mi,
     +  cr,dsq,usm,rno,rna,s1,s2,s3,s4,x1,x3,
     +  rtnew(1),rtnew(4),rtnew(7),rtnew(2),rtnew(5),rtnew(8),
     +  rtnew(3),rtnew(6),rtnew(9),rtnew(10),rtnew(11),rtnew(12)
c
 6010 format (
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' SI = RMS * Nmin / Nmatch             = ',f12.5/
     +  ' MI = (1+Nmatch)/{(1+W*RMS)*(1+Nmin)} = ',f12.5/
     +  ' CR = Maiorov-Crippen RHO (0-2)       = ',f12.5/
     +  ' Estimated RMSD for 2 random proteins = ',f10.3,' A'/
     +  ' RR = Relative RMSD                   = ',f12.5/
     +  ' NR = Normalised RMSD (100)           = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A'/
     +  ' SAS(1) = cRMS * (100/Nmatch)         = ',f10.3,' A'/
     +  ' SAS(2) = cRMS * (100/Nmatch)^2       = ',f10.3,' A'/
     +  ' SAS(3) = cRMS * (100/Nmatch)^3       = ',f10.3,' A'/
     +  ' SAS(4) = cRMS * (100/Nmatch)^4       = ',f10.3,' A'/
     +  ' RMS delta B for matched atoms        = ',f9.3,' A2'/
     +  ' Corr. coefficient matched atom Bs    = ',f12.3/
     +  ' Rotation     : ',3f12.8/16x,3f12.8/16x,3f12.8/
     +  ' Translation  : ',3f12.4)
c
cxyz - check progress
c
      if (optcri .eq. 'SI') then
        if (si .lt. best) then
          best = si
          goto 6996
        end if
      else if (optcri .eq. 'MI') then
        if (mi .gt. best) then
          best = mi
          goto 6996
        end if
      else if (optcri .eq. 'S1') then
        if (s1 .lt. best) then
          best = s1
          goto 6996
        end if
      else if (optcri .eq. 'S2') then
        if (s2 .lt. best) then
          best = s2
          goto 6996
        end if
      else if (optcri .eq. 'S3') then
        if (s3 .lt. best) then
          best = s3
          goto 6996
        end if
      else if (optcri .eq. 'S4') then
        if (s4 .lt. best) then
          best = s4
          goto 6996
        end if
      else if (optcri .eq. 'RM') then
        if (rms .lt. best) then
          best = rms
          goto 6996
        end if
      else if (optcri .eq. 'CR') then
        if (cr .lt. best) then
          best = cr
          goto 6996
        end if
      else if (optcri .eq. 'RR') then
        if (usm .lt. best) then
          best = usm
          goto 6996
        end if
      else if (optcri .eq. 'NR') then
        if (rno .lt. best) then
          best = rno
          goto 6996
        end if
      else if (optcri .eq. 'NM') then
        if (float(nmat) .gt. best) then
          best = float(nmat)
          goto 6996
        end if
      end if
c
      if (lprint) then
        write (*,*)
ccc        call prompt (' Fit deteriorated in this cycle !')
        call prompt (' Fit did not improve in this cycle !')
        call prompt (' Alignment based on previous operator !')
      end if
      goto 7000
c
 6996 continue
c
      do i=1,nuse1
        sptrij (i) = ptrij (i)
      end do
c
      do i=1,12
        rtx (i) = rtnew (i)
      end do
c
      nmatch (imol,jmol) = nmat
      rmsd (imol,jmol) = rms
      rmsb (imol,jmol) = x1
      corb (imol,jmol) = x3
      simind (imol,jmol) = si
      matchi (imol,jmol) = mi
      cripp  (imol,jmol) = cr
      rrmsd  (imol,jmol) = usm
      normsd (imol,jmol) = rno
      rmsdna (imol,jmol) = rna
      sas1 (imol,jmol) = s1
      sas2 (imol,jmol) = s2
      sas3 (imol,jmol) = s3
      sas4 (imol,jmol) = s4
c
      dmax = dmax * ddec
c
      uselen = max (1, (uselen + lendec))
c
      if (ncycle .lt. maxcyc) goto 10
c
 7000 continue
c
      if (lprint) write (*,*)
c
c ... store improved operator
c
      do i=1,12
        rtlsq (i,imol,jmol) = rtx (i)
      end do
c
c ... apply current operator
c
      call vecrtv (buffj,buffk,nuse2,rtx(1),rtx(10))
c
c ... show alignment
c
      jold = -1
      kold = -1
      jchn = '?'
      kchn = '?'
c
      nsame = 0
      nali = 0
c
      do i=1,nuse1
        if (sptrij(i) .gt. 0) then
          j = ptri (i)
          k = ptrj (sptrij (i))
          l = sptrij (i)
          xdist = (buffi(1,i)-buffk(1,l)) ** 2 +
     +            (buffi(2,i)-buffk(2,l)) ** 2 +
     +            (buffi(3,i)-buffk(3,l)) ** 2
          xdist = sqrt (xdist)
c
          nali = nali + 1
          buffb3 (nali) = xdist
c
ccc          if (i .gt. 1) then 
            symbol = ' '
            if (resnam(j,imol) .eq. resnam(k,jmol)) symbol = '*'
            if ( (jold+1) .eq. iresid(j,imol) .and.
     +           (kold+1) .eq. iresid(k,jmol) .and.
     +           jchn     .eq. achain(j,imol) .and.
     +           kchn     .eq. achain(k,jmol) ) then
              if (lprint) write (*,6020)
     +          resnam(j,imol),achain(j,imol),iresid(j,imol),
     +          resnam(k,jmol),achain(k,jmol),iresid(k,jmol),
     +          xdist,symbol
            else
              if (lprint) write (*,6030)
     +          resnam(j,imol),achain(j,imol),iresid(j,imol),
     +          resnam(k,jmol),achain(k,jmol),iresid(k,jmol),
     +          xdist,symbol
            end if
ccc          end if
          jold = iresid(j,imol)
          kold = iresid(k,jmol)
          jchn = achain(j,imol)
          kchn = achain(k,jmol)
c
          if (resnam(j,imol) .eq. resnam(k,jmol)) nsame = nsame + 1
c
        end if
      end do
c
 6020 format ('          ',a3,'-',a1,i4,' <===> ',a3,'-',a1,i4,
     + ' @ ',f8.2,' A',1x,a1)
 6030 format (' Fragment ',a3,'-',a1,i4,' <===> ',a3,'-',a1,i4,
     + ' @ ',f8.2,' A',1x,a1)
c
      if (lprint) write (*,*)
      if (lprint) call ivalut (
     +  ' Nr of residues in mol1   :',1,nuse1)
      if (lprint) call ivalut (
     +  ' Nr of residues in mol2   :',1,nuse2)
      if (lprint) call ivalut (
     +  ' Nr of matched residues   :',1,nmatch(imol,jmol))
      if (lprint) call ivalut (
     +  ' Nr of identical residues :',1,nsame)
      x1 = 100.0 * float(nsame) / float(nmatch(imol,jmol))
      if (lprint) call fvalut (
     +  ' % identical of matched   :',1,x1)
c
      x1 = 100.0 * float(nmatch(imol,jmol)) / float(nuse1)
      x2 = 100.0 * float(nsame) / float(nuse1)
      if (lprint) call fvalut (
     +  ' % matched   of mol1      :',1,x1)
      if (lprint) call fvalut (
     +  ' % identical of mol1      :',1,x2)
      dvalue = 0.0001 * x1 * x2
      if (lprint) call fvalut (
     +  ' D-value    for mol1      :',1,dvalue)
c
      x1 = 100.0 * float(nmatch(imol,jmol)) / float(nuse2)
      x2 = 100.0 * float(nsame) / float(nuse2)
      if (lprint) call fvalut (
     +  ' % matched   of mol2      :',1,x1)
      if (lprint) call fvalut (
     +  ' % identical of mol2      :',1,x2)
      dvalue = 0.0001 * x1 * x2
      if (lprint) call fvalut (
     +  ' D-value    for mol2      :',1,dvalue)
c
      if (lprint) write (*,*)
c
      if (lprint) call distat (nali,buffb3)
c
      if (lprint) write (*,*)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine setdef (level)
c
      include 'lsqman.incl'
c
      integer level
c
code ...
c
      seqcri = 'ON'
      ashift = 'ON'
      lshift = .true.
c
c ... 0 = defaults
c
      if (level .lt. 1 .or. level .gt. 5) then
        dismax = 3.5
        decay  = 1.0
        rmswgt = 1.0
        minlen = 3
        lendec = 0
        maxcyc = 10
        optcri = 'CR'
c
c ... 1 = coarse 6 A fit
c
      else if (level .eq. 1) then
        dismax = 6.0
        decay  = 1.0
        rmswgt = 0.5
        minlen = 5
        lendec = 0
        maxcyc = 25
        optcri = 'NM'
c
c ... 2 = intermediate 4 A fit
c
      else if (level .eq. 2) then
        dismax = 4.0
        decay  = 0.975
        rmswgt = 0.5
        minlen = 4
        lendec = 0
        maxcyc = 10
        optcri = 'CR'
c
c ... 3 = fine-tune 3 A fit
c
      else if (level .eq. 3) then
        dismax = 3.0
        decay  = 0.975
        rmswgt = 0.5
        minlen = 5
        lendec = 0
        maxcyc = 10
        optcri = 'CR'
c
c ... 4 = 2 A fit of similar molecules
c
      else if (level .eq. 4) then
        dismax = 2.0
        decay  = 0.975
        rmswgt = 0.5
        minlen = 7
        lendec = 0
        maxcyc = 10
        optcri = 'CR'
c
c ... 5 = nucleic acids
c
      else if (level .eq. 5) then
        dismax = 4.0
        decay  = 1.0
        rmswgt = 0.5
        minlen = 3
        lendec = 0
        maxcyc = 10
        optcri = 'CR'
c
      end if
c
      return
      end
c
c
c
      subroutine seqnce (imol)
c
      include 'lsqman.incl'
c
      integer imol,previd,nout,kk
c
      character prevnm*3,prevch*1
c
code ...
c
      nout = 0
      previd = -1
      prevnm = '???'
      prevch = '?'
c
 4711 format (a6,a)
 4713 format (1x,2i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,
     +        3f8.3,2f6.2,6x,a4)
c
      write (*,*)
      do kk=1,natoms(imol)
        if (previd .ne. iresid(kk,imol) .or.
     +      prevnm .ne. resnam(kk,imol) .or.
     +      prevch .ne. achain(kk,imol)) then
          nout = nout + 1
          write (*,4713) nout,atomnr(kk,imol),
     +      atmnam(kk,imol),resnam(kk,imol),achain(kk,imol),
     +      iresid(kk,imol),atmxyz(1,kk,imol),atmxyz(2,kk,imol),
     +      atmxyz(3,kk,imol),qatom(kk,imol),batom(kk,imol),
     +      axplor(kk,imol)
          previd = iresid(kk,imol)
          prevnm = resnam(kk,imol)
          prevch = achain(kk,imol)
        end if
      end do
      write (*,*)
c
      return
      end
c
c
c
      subroutine plotem (what,imol,irange,jmol,jrange,filnam,
     +                   dum1,dum2,dum3,ierr)
c
      include 'lsqman.incl'
c
      integer maxopt
      parameter (maxopt = 25)
c
      real dist,angle,tangle,rms,pgt,ave,dmax,xmin,xmax,rtx(12)
      real d1,d2,d3,e1,e2,e3,q1,q2,q3,q4,q5,q6,q7,q8,corphi,corpsi
      real hiscut,bin,dum1,dum2,dum3,dmin,outcut
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j,ngt10,leng1
      integer idum,jdum,kdum,ldum,nuse,iptr,jptr,ires,jres,nmax
      integer inptr,jnptr,icptr,jcptr,ic1ptr,jc1ptr,in1ptr,jn1ptr
      integer i1ptr,i2ptr,i3ptr,i4ptr,j1ptr,j2ptr,j3ptr,j4ptr,nnot
      integer ngt25,ipptr,jpptr,ip1ptr,jp1ptr,ic11ptr,jc11ptr
c
      logical xinter
c
      character irange*(*),jrange*(*),ichn*1,jchn*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20
      character what*(*),filnam*(*)
      character line*256
c
code ...
c
      ierr = 0
c
      if (what(1:2) .eq. 'PH') then
        call prompt (' Delta-Phi/Delta-Psi plot')
        hiscut = 180.0
        bin = 10.0
        outcut = 10.0
      else if (what(1:2) .eq. 'ET') then
        call prompt (' Delta-Eta/Delta-Theta plot')
        hiscut = 180.0
        bin = 10.0
        outcut = 25.0
      else if (what(1:2) .eq. 'DI') then
        call prompt (' Central-atom distance plot')
        call textut (' Central atom type :',atypes(1))
        hiscut = 100.0
        bin = 0.5
        outcut = 1.0
      else if (what(1:2) .eq. 'DD') then
        call prompt (' Central-atom delta-dihedral plot')
        call textut (' Central atom type :',atypes(1))
        hiscut = 180.0
        bin = 10.0
        outcut = 10.0
      else if (what(1:2) .eq. 'D1') then
        call prompt (' Delta-1/Delta-2 plot')
        hiscut = 180.0
        bin = 10.0
        outcut = 10.0
      else
        call errcon ('PLOTEM - Bug ! Invalid option !')
        ierr = -1
        return
      end if
      if (dum1 .gt. 0.0 .and. dum1 .le. 180.0) hiscut = dum1
      if (dum2 .gt. 0.0 .and. dum2 .le. hiscut) bin = dum2
      if (dum3 .gt. 0.0 .and. dum3 .le. 180.0) outcut = dum3
c
      if (what(1:2) .eq. 'DI' .or. what(1:2) .eq. 'DD') then
        if (atypes(1) .eq. 'ALL ' .or. atypes(1) .eq. 'NONH' .or.
     +      atypes(1) .eq. 'SIDE' .or. atypes(1) .eq. 'TRAC') then
          call errcon ('Invalid atom type !')
          ierr = -1
          return
        end if
      end if
c
      write (*,6000) name(imol)(1:leng1(name(imol))),
     +    irange(1:leng1(irange)),
     +    name(jmol)(1:leng1(name(jmol))),
     +    jrange(1:leng1(jrange)),outcut,hiscut,bin
c
 6000 format (' Plot of ',a,1x,a/
     +        ' And     ',a,1x,a/
     +        ' Outlier list cut-off  : ',f8.2/
     +        ' Histogram upper limit : ',f8.2/
     +        ' Histogram bin size    : ',f8.2)
c
      do i=1,12
        rtx (i) = rtlsq (i,imol,jmol)
      end do
c
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range1')
        goto 9900
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range2')
        goto 9900
      end if
c
      if (nopt1 .ne. nopt2) then
        call errcon ('Different nr of zones')
        goto 9900
      end if
c
      nuse = 0
      iptr = 0
      jptr = 0
c
      dmax = 0.0
      dmin = 0.0
      nnot = 0
c
      do i=1,nopt1
c
        ichn = optpa1(i)(1:1)
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range1 : '//optpa1(i))
ccc          goto 9900
        end if
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
          call errcon ('Not a zone in range1 : '//optpa1(i))
          goto 9900
        end if
ccc      print *,optpa1(i)(2:j-1)
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
ccc      print *,optpa1(i)(j+1:)
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon ('Invalid zone in range1 : '//optpa1(i))
          goto 9900
        end if
c
        jchn = optpa2(i)(1:1)
        if (index('1234567890',jchn) .gt. 0) then
          jchn = ' '
ccc        else if (jchn .lt. 'A' .or. jchn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range2 : '//optpa2(i))
ccc          goto 9900
        end if
        j = index (optpa2(i),'-')
        if (j .le. 0) j = index (optpa2(i),':')
        if (j .gt. 0) then
          call errcon (
     +      'Zone ignored in range2 : '//optpa2(i))
          optpa2(i) = optpa2(i)(1:j-1)
        end if
ccc      print *,optpa2(i)(2:)
        if (jchn .eq. ' ') then
          call str2i (optpa2(i)(1:),kdum,ierr)
        else
          call str2i (optpa2(i)(2:),kdum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        ldum = kdum + (jdum - idum)
c
        do ires=idum,jdum
            jres = kdum + (ires - idum)
ccc      print *,imol,ichn,ires,atypes(j)
ccc      print *,jmol,jchn,jres
          if (what(1:2) .eq. 'DI') then
            call getptr (imol,ichn,ires,atypes(1),
     +                   -1,jchn,jres,iptr,jptr,ierr)
          else if (what(1:2) .eq. 'DD') then
            call getptr (imol,ichn,ires,atypes(1),
     +                   jmol,jchn,jres,iptr,jptr,ierr)
          else if (what(1:2) .eq. 'ET') then
            call getptr (imol,ichn,ires,' C4*',
     +                   jmol,jchn,jres,iptr,jptr,ierr)
          else
            call getptr (imol,ichn,ires,' CA ',
     +                   jmol,jchn,jres,iptr,jptr,ierr)
          end if
ccc      print *,ierr,iptr,jptr
            if (ierr .ne. 0) goto 6196
c
            buffi (1,nuse+1) = 0.0
            buffj (1,nuse+1) = 0.0
            buffk (1,nuse+1) = 0.0
            buffl (1,nuse+1) = 0.0
c
c ... PHI/PSI PLOT or D1/D2 PLOT
c
            if (what(1:2) .eq. 'PH' .or.
     +          what(1:2) .eq. 'D1') then
              call getptr (imol,ichn,ires,' N  ',
     +                   jmol,jchn,jres,inptr,jnptr,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' C  ',
     +                   jmol,jchn,jres,icptr,jcptr,ierr)
              if (ierr .ne. 0) goto 6996
c
              nuse = nuse + 1
c
              buffb1 (nuse) = float (ires)
              buffm (1,nuse) = float (inptr)
              buffm (2,nuse) = float (jnptr)
c
              call getptr (imol,ichn,ires-1,' C  ',
     +                   jmol,jchn,jres-1,ic1ptr,jc1ptr,ierr)
              if (ierr .ne. 0) goto 6986
c
              if (dist(ic1ptr,inptr,atmxyz(1,1,imol)) .gt. 3.0)
     +          goto 6986
              if (dist(jc1ptr,jnptr,atmxyz(1,1,jmol)) .gt. 3.0)
     +          goto 6986
c
              buffi (1,nuse) = tangle (ic1ptr,inptr,iptr,
     +          icptr,atmxyz(1,1,imol))
              buffj (1,nuse) = tangle (jc1ptr,jnptr,jptr,
     +          jcptr,atmxyz(1,1,jmol))
c
 6986         continue
              call getptr (imol,ichn,ires+1,' N  ',
     +                   jmol,jchn,jres+1,in1ptr,jn1ptr,ierr)
              if (ierr .ne. 0) goto 6976
c
              if (dist(icptr,in1ptr,atmxyz(1,1,imol)) .gt. 3.0)
     +          goto 6976
              if (dist(jcptr,jn1ptr,atmxyz(1,1,jmol)) .gt. 3.0)
     +          goto 6976
c
              buffk (1,nuse) = tangle (inptr,iptr,
     +          icptr,in1ptr,atmxyz(1,1,imol))
              buffl (1,nuse) = tangle (jnptr,jptr,
     +          jcptr,jn1ptr,atmxyz(1,1,jmol))
c
 6976         continue
c
 6996         continue
c
c ... ETA/THETA PLOT
c
            else if (what(1:2) .eq. 'ET') then
c
c ... IPTR/JPTR = C4* (i)
c     IPPTR/..  = P   (i)
c     IC1PTR/.. = C4* (i-1)
c     IP1PTR/.. = P   (i+1)
c     IC11PTR/..= C4* (i+1)
c
              call getptr (imol,ichn,ires,' P  ',
     +                   jmol,jchn,jres,ipptr,jpptr,ierr)
              if (ierr .ne. 0) goto 6916
c
              nuse = nuse + 1
c
              buffb1 (nuse)  = float (ires)
              buffm (1,nuse) = float (iptr)
              buffm (2,nuse) = float (jptr)
c
              call getptr (imol,ichn,ires-1,' C4*',
     +                   jmol,jchn,jres-1,ic1ptr,jc1ptr,ierr)
              if (ierr .ne. 0) goto 6916
              call getptr (imol,ichn,ires+1,' P  ',
     +                   jmol,jchn,jres+1,ip1ptr,jp1ptr,ierr)
              if (ierr .ne. 0) goto 6916
c
              if (dist(ic1ptr,ipptr,atmxyz(1,1,imol)) .gt. 8.0)
     +          goto 6916
              if (dist(jc1ptr,jpptr,atmxyz(1,1,jmol)) .gt. 8.0)
     +          goto 6916
              if (dist(iptr,ip1ptr,atmxyz(1,1,imol)) .gt. 8.0)
     +          goto 6916
              if (dist(jptr,jp1ptr,atmxyz(1,1,jmol)) .gt. 8.0)
     +          goto 6916
c
              buffi (1,nuse) = tangle (ic1ptr,ipptr,iptr,
     +          ip1ptr,atmxyz(1,1,imol))
              buffj (1,nuse) = tangle (jc1ptr,jpptr,jptr,
     +          jp1ptr,atmxyz(1,1,jmol))
c
              call getptr (imol,ichn,ires+1,' C4*',
     +                   jmol,jchn,jres+1,ic11ptr,jc11ptr,ierr)
              if (ierr .ne. 0) goto 6916
c
              buffk (1,nuse) = tangle (ipptr,iptr,
     +          ip1ptr,ic11ptr,atmxyz(1,1,imol))
              buffl (1,nuse) = tangle (jpptr,jptr,
     +          jp1ptr,jc11ptr,atmxyz(1,1,jmol))
c
 6916         continue
c
c ... DISTANCE PLOT
c
            else if (what(1:2) .eq. 'DI') then
c
              call getptr (imol,ichn,ires,atypes(1),
     +                     jmol,jchn,jres,inptr,jnptr,ierr)
c
ccc              if (ierr .ne. 0) goto 6896
c 
              if (ierr .eq. 0) then
c
                nuse = nuse + 1
c
                buffb1 (nuse) = float (ires)
                buffm (1,nuse) = float (inptr)
                buffm (2,nuse) = float (jnptr)
c
                do j=1,3
                  buffj(j,1) = atmxyz(j,inptr,imol)
                  buffk(j,1) = atmxyz(j,jnptr,jmol)
                end do
                call vecrtv (buffk,buffl,1,rtx(1),rtx(10))
c
                buffi (1,nuse) = sqrt (
     +            (buffj(1,1)-buffl(1,1))**2 +
     +            (buffj(2,1)-buffl(2,1))**2 +
     +            (buffj(3,1)-buffl(3,1))**2 )
c
                dmax = max (dmax,buffi(1,nuse))
              else
                nuse = nuse + 1
                buffb1 (nuse) = float(ires)
                buffm (1,nuse) = float (inptr)
                buffm (2,nuse) = float (jnptr)
                buffi (1,nuse) = -1.0
                dmin = -1.0
                nnot = nnot + 1
ccc                call prompt (' Missing residue !')
              end if
c
 6896         continue
c
c ... DELTA-DIHEDRAL PLOT
c
            else if (what(1:2) .eq. 'DD') then
c
              call getptr (imol,ichn,ires,atypes(1),
     +                     jmol,jchn,jres,i2ptr,j2ptr,ierr)
              if (ierr .ne. 0) goto 6796
c
              call getptr (imol,ichn,ires-1,atypes(1),
     +                     jmol,jchn,jres-1,i1ptr,j1ptr,ierr)
              if (ierr .ne. 0) goto 6796
c
              call getptr (imol,ichn,ires+1,atypes(1),
     +                     jmol,jchn,jres+1,i3ptr,j3ptr,ierr)
              if (ierr .ne. 0) goto 6796
c
              call getptr (imol,ichn,ires+2,atypes(1),
     +                     jmol,jchn,jres+2,i4ptr,j4ptr,ierr)
              if (ierr .ne. 0) goto 6796
c
              d1 = dist (i1ptr,i2ptr,atmxyz(1,1,imol))
              e1 = dist (j1ptr,j2ptr,atmxyz(1,1,jmol))
              d2 = dist (i2ptr,i3ptr,atmxyz(1,1,imol))
              e2 = dist (j2ptr,j3ptr,atmxyz(1,1,jmol))
              d3 = dist (i3ptr,i4ptr,atmxyz(1,1,imol))
              e3 = dist (j3ptr,j4ptr,atmxyz(1,1,jmol))
c
              xmin = min (d1,d2,d3,e1,e2,e3)
              xmax = max (d1,d2,d3,e1,e2,e3)
              if ( xmax .gt. (xmin+5.0) ) goto 6796
c
              nuse = nuse + 1
c
              buffb1 (nuse)  = float (ires)
              buffm (1,nuse) = float (i2ptr)
              buffm (2,nuse) = float (j2ptr)
c
              buffi (1,nuse) = tangle (i1ptr,i2ptr,
     +          i3ptr,i4ptr,atmxyz(1,1,imol))
              buffj (1,nuse) = tangle (j1ptr,j2ptr,
     +          j3ptr,j4ptr,atmxyz(1,1,jmol))
c
              buffk (1,nuse) = angle (i1ptr,i2ptr,
     +          i3ptr,atmxyz(1,1,imol))
              buffl (1,nuse) = angle (j1ptr,j2ptr,
     +          j3ptr,atmxyz(1,1,jmol))
c
 6796         continue
c
            end if
c
 6196       continue
c
        end do
c
      end do
c
      call jvalut (' Nr of residues found :',1,nuse)
      if (nuse .lt. 3) then
        call errcon ('Fewer than 3 residues; cannot plot')
        goto 9900
      end if
c
      nmax = 1 + int (hiscut/bin)
      if (nmax .gt. maxatm) nmax = maxatm
      do i=1,nmax
        buffb3 (i) = 0.0
        buffb4 (i) = 0.0
      end do
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Residue number'
c
      if (what(1:2) .eq. 'PH') then
        write (iunit,6200) 'YLABEL',
     +    'Delta PHI and Delta PSI'
        write (iunit,6200) 'REMARK',
     +    'Plot of delta PHI (solid blue) and delta PSI ',
     +    '(dotted red) as a function of residue nr'
cc        write (iunit,6200) 'REMARK',
cc     +    'Plot of delta PHI is the solid blue curve'
cc        write (iunit,6200) 'REMARK',
cc     +    'Plot of delta PSI is the dotted red curve'
        write (iunit,6200) 'REMARK',
     +    'Values are those of mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        write (iunit,6200) 'REMARK',
     +    '     minus those of mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
        write (iunit,6210) 'XYVIEW',
     +    nint(buffb1(1))-1, nint(buffb1(nuse))+1,-180,180
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DPHI R ',
     +    nuse,' (9f8.2)'
c
      else if (what(1:2) .eq. 'ET') then
        write (iunit,6200) 'YLABEL',
     +    'Delta ETA and Delta THETA'
        write (iunit,6200) 'REMARK',
     +    'Plot of delta ETA (solid blue) and delta THETA ',
     +    '(dotted red) as a function of residue nr'
cc        write (iunit,6200) 'REMARK',
cc     +    'Plot of delta ETA is the solid blue curve'
cc        write (iunit,6200) 'REMARK',
cc     +    'Plot of delta THETA is the dotted red curve'
        write (iunit,6200) 'REMARK',
     +    'Values are those of mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        write (iunit,6200) 'REMARK',
     +    '     minus those of mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
        write (iunit,6210) 'XYVIEW',
     +    nint(buffb1(1))-1, nint(buffb1(nuse))+1,-180,180
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DETA R ',
     +    nuse,' (9f8.2)'
c
      else if (what(1:2) .eq. 'D1') then
        write (iunit,6200) 'YLABEL',
     +    'Delta-1 and Delta-2'
        write (iunit,6200) 'REMARK',
     +    'Plot of Delta-1 (solid blue) and Delta-2 ',
     +    '(dotted red) as a function of residue nr'
        write (iunit,6200) 'REMARK',
     +    'Delta-1 = delta (PHI(i+1) - PSI(i)) mol 1 and 2'
        write (iunit,6200) 'REMARK',
     +    'Delta-2 = delta (PHI(i+1) + PSI(i)) mol 1 and 2'
        write (iunit,6200) 'REMARK',
     +    'Values are those of mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        write (iunit,6200) 'REMARK',
     +    '     minus those of mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
        write (iunit,6210) 'XYVIEW',
     +    nint(buffb1(1))-1, nint(buffb1(nuse))+1,-180,180
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DELTA1 R ',
     +    nuse,' (9f8.2)'
c
      else if (what(1:2) .eq. 'DI') then
c
        call ivalut (' Residues not matched in mol2 :',1,nnot)
        if (dmax .lt. 5.0*abs(dmin)) then
          dmin = min(-dmax/5.0,-0.001)
          do i=1,nuse
            if (buffi(1,i) .lt. 0.0) buffi(1,i) = dmin
          end do
          call fvalut (' Dummy distance set to :',1,dmin)
        end if
c
        write (iunit,6200) 'YLABEL',
     +    'Distance between corresponding ',atypes(1),
     +    ' atoms after superimposing'
        write (iunit,6200) 'REMARK',
     +    'Values are those of mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        write (iunit,6200) 'REMARK',
     +    '        compared to mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
        write (iunit,6213) 'XYVIEW',
     +    buffb1(1)-1.0,buffb1(nuse)+1.0,dmin-0.05,dmax+0.05
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DIST R ',
     +    nuse,' (9f8.2)'
c
      else if (what(1:2) .eq. 'DD') then
        write (iunit,6200) 'YLABEL',
     +    'Delta-angles/dihedrals of sequential ',atypes(1),
     +    ' atoms'
        write (iunit,6200) 'REMARK',
     +    'Delta X-X*-X-X dihedrals <-180,+180] solid blue curve'
        write (iunit,6200) 'REMARK',
     +    'Absolute Delta X-X*-X angles [0,+180] dotted red curve'
        write (iunit,6200) 'REMARK',
     +    'Values are those of mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        write (iunit,6200) 'REMARK',
     +    '        compared to mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
        write (iunit,6210) 'XYVIEW',
     +    nint(buffb1(1))-1, nint(buffb1(nuse))+1,-180,180
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DDIHE R ',
     +    nuse,' (9f8.2)'
c
      end if
c
      write (iunit,6210) 'NPOINT',nuse
      write (iunit,6210) 'COLOUR',4
c
      write (iunit,6210) 'LINE  ',
     +    nint(buffb1(1))-1,0,nint(buffb1(nuse))+1,0
c
      write (iunit,6200) 'XVALUE','*'
      write (iunit,6220) (nint(buffb1(i)),i=1,nuse)
c
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
      write (iunit,6200) 'YVALUE','*'
c
 6268 format (a,1x,a,f8.2,a,f8.2)
 6269 format (a,1x,a,i8,a,f8.2)
c
 6206 format (' >>> ',a,1x,a3,1x,a1,i4,' - ',
     +        a3,1x,a1,i4,' = ',f8.2)
 6208 format (' >>> ',a,1x,a3,1x,a1,i4,1x,a4,' - ',
     +        a3,1x,a1,i4,1x,a4,' = ',f8.2)
c
      write (*,*)
c
      if (what(1:2) .eq. 'PH') then
        rms = 0.0
        ngt10 = 0
        ave = 0
        do i=1,nuse
          call fixadf (buffi(1,i),buffj(1,i),buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 10.0) ngt10 = ngt10 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-PHI',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb3(j) = buffb3(j) + 1.0
        end do
c
        write (iunit,6230) (buffb2(i),i=1,nuse)
        call xystat (buffi,buffj,nuse,q1,q2,corphi,q3,q4,q5,q6,q7,q8)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt10) / float(nuse)
        call fvalut (' RMS delta PHI       :',1,rms)
        call fvalut (' Average |delta PHI| :',1,ave)
        call ivalut (' Nr |delta PHI| > 10 :',1,ngt10)
        call fvalut (' Percentage          :',1,pgt)
        call fvalut (' Corr. coeff. PHI    :',1,corphi)
        write (*,*)
        write (iunit,6268) 'REMARK','RMS delta PHI = ',rms,
     +    ' +++ Average |delta PHI| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta PHI| > 10 = ',
     +    ngt10,' +++ Percentage = ',pgt
        write (iunit,6268) 'REMARK','Corr. coeff. PHI(1)-PHI(2) = ',
     +    corphi
c
        write (iunit,6200) 'MORE  '
        write (iunit,6210) 'COLOUR',1
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DPSI R ',
     +    nuse,' (9f8.2)'
        call pretty (line)
        write (iunit,'(a)') line(1:leng1(line))
c
        write (iunit,6200) 'YVALUE','*'
c
        rms = 0.0
        ngt10 = 0
        ave = 0
        do i=1,nuse
          call fixadf (buffk(1,i),buffl(1,i),buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 10.0) ngt10 = ngt10 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-PSI',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb4(j) = buffb4(j) + 1.0
        end do
c
        write (iunit,6230) (buffb2(i),i=1,nuse)
        call xystat (buffk,buffl,nuse,q1,q2,corpsi,q3,q4,q5,q6,q7,q8)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt10) / float(nuse)
        call fvalut (' RMS delta PSI       :',1,rms)
        call fvalut (' Average |delta PSI| :',1,ave)
        call ivalut (' Nr |delta PSI| > 10 :',1,ngt10)
        call fvalut (' Percentage          :',1,pgt)
        call fvalut (' Corr. coeff. PSI    :',1,corpsi)
        write (iunit,6268) 'REMARK','RMS delta PSI = ',rms,
     +    ' +++ Average |delta PSI| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta PSI| > 10 = ',
     +    ngt10,' +++ Percentage = ',pgt
        write (iunit,6268) 'REMARK','Corr. coeff. PSI(1)-PSI(2) = ',
     +    corpsi
c
      else if (what(1:2) .eq. 'ET') then
        rms = 0.0
        ngt25 = 0
        ave = 0
        do i=1,nuse
          call fixadf (buffi(1,i),buffj(1,i),buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 25.0) ngt25 = ngt25 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-ETA',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb3(j) = buffb3(j) + 1.0
        end do
c
        write (iunit,6230) (buffb2(i),i=1,nuse)
        call xystat (buffi,buffj,nuse,q1,q2,corphi,q3,q4,q5,q6,q7,q8)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt25) / float(nuse)
        call fvalut (' RMS delta ETA       :',1,rms)
        call fvalut (' Average |delta ETA| :',1,ave)
        call ivalut (' Nr |delta ETA| > 25 :',1,ngt25)
        call fvalut (' Percentage          :',1,pgt)
        call fvalut (' Corr. coeff. ETA    :',1,corphi)
        write (*,*)
        write (iunit,6268) 'REMARK','RMS delta ETA = ',rms,
     +    ' +++ Average |delta ETA| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta ETA| > 25 = ',
     +    ngt25,' +++ Percentage = ',pgt
        write (iunit,6268) 'REMARK','Corr. coeff. ETA(1)-ETA(2) = ',
     +    corphi
c
        write (iunit,6200) 'MORE  '
        write (iunit,6210) 'COLOUR',1
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DTHETA R ',
     +    nuse,' (9f8.2)'
        call pretty (line)
        write (iunit,'(a)') line(1:leng1(line))
c
        write (iunit,6200) 'YVALUE','*'
c
        rms = 0.0
        ngt25 = 0
        ave = 0
        do i=1,nuse
          call fixadf (buffk(1,i),buffl(1,i),buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 25.0) ngt25 = ngt25 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-THETA',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb4(j) = buffb4(j) + 1.0
        end do
c
        write (iunit,6230) (buffb2(i),i=1,nuse)
        call xystat (buffk,buffl,nuse,q1,q2,corpsi,q3,q4,q5,q6,q7,q8)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt25) / float(nuse)
        call fvalut (' RMS delta THETA     :',1,rms)
        call fvalut (' Average |delta THE| :',1,ave)
        call ivalut (' Nr |delta THE| > 25 :',1,ngt25)
        call fvalut (' Percentage          :',1,pgt)
        call fvalut (' Corr. coeff. THETA  :',1,corpsi)
        write (iunit,6268) 'REMARK','RMS delta THETA = ',rms,
     +    ' +++ Average |delta THETA| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta THETA| > 25 = ',
     +    ngt25,' +++ Percentage = ',pgt
        write (iunit,6268) 'REMARK',
     +    'Corr. coeff. THETA(1)-THETA(2) = ',corpsi
c
      else if (what(1:2) .eq. 'D1') then
c
        rms = 0.0
        ngt10 = 0
        ave = 0
        do i=1,nuse-1
          call fixdif (buffi(1,i+1),-buffk(1,i),d1)
          call fixdif (buffj(1,i+1),-buffl(1,i),d2)
          call fixdif (d1,d2,buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 10.0) ngt10 = ngt10 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-1',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb3(j) = buffb3(j) + 1.0
        end do
c
        write (iunit,6230) (buffb2(i),i=1,nuse)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt10) / float(nuse)
        call fvalut (' RMS delta Delta-1   :',1,rms)
        call fvalut (' Average |delta D-1| :',1,ave)
        call ivalut (' Nr |delta D-1| > 10 :',1,ngt10)
        call fvalut (' Percentage          :',1,pgt)
        write (iunit,6268) 'REMARK','RMS delta Delta-1 = ',rms,
     +    ' +++ Average |delta Delta-1| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta Delta-1| > 10 = ',
     +    ngt10,' +++ Percentage = ',pgt
c
        write (iunit,6200) 'MORE  '
        write (iunit,6210) 'COLOUR',1
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DELTA2 R ',
     +    nuse,' (9f8.2)'
        call pretty (line)
        write (iunit,'(a)') line(1:leng1(line))
c
        write (iunit,6200) 'YVALUE','*'
c
        rms = 0.0
        ngt10 = 0
        ave = 0
        do i=1,nuse
          call fixdif (buffi(1,i+1),buffk(1,i),d1)
          call fixdif (buffj(1,i+1),buffl(1,i),d2)
          call fixdif (d1,d2,buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 10.0) ngt10 = ngt10 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-2',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb4(j) = buffb4(j) + 1.0
        end do
        write (iunit,6230) (buffb2(i),i=1,nuse)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt10) / float(nuse)
        call fvalut (' RMS delta Delta-2   :',1,rms)
        call fvalut (' Average |delta D-2| :',1,ave)
        call ivalut (' Nr |delta D-2| > 10 :',1,ngt10)
        call fvalut (' Percentage          :',1,pgt)
        write (iunit,6268) 'REMARK','RMS delta Delta-2 = ',rms,
     +    ' +++ Average |delta Delta-2| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta Delta-2| > 10 = ',
     +    ngt10,' +++ Percentage = ',pgt
c
      else if (what(1:2) .eq. 'DI') then
        write (iunit,6230) (buffi(1,i),i=1,nuse)
        ave = 0.0
        xmin = 999.99
        xmax = -999.9
        do i=1,nuse
          if (buffi(1,i) .gt. 0.0) then
            ave = ave + buffi(1,i)
            xmin = min (xmin,buffi(1,i))
            xmax = max (xmax,buffi(1,i))
c
          if (buffi(1,i) .gt. outcut) then
            write (*,6208) 'Distance',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        atmnam(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        atmnam(nint(buffm(2,i)),imol),
     +        buffi(1,i)
          end if
c
            j = 1 + int(buffi(1,i)/bin)
            if (j .gt. nmax) j = nmax
            buffb3(j) = buffb3(j) + 1.0
          end if
        end do
        ave = ave / float(nuse)
        call fvalut (' Average distance :',1,ave)
        call fvalut (' Minimum distance :',1,xmin)
        call fvalut (' Maximum distance :',1,xmax)
        write (iunit,6268) 'REMARK','Average distance = ',ave
        write (iunit,6268) 'REMARK','Minimum = ',xmin,
     +    ' +++ Maximum = ',xmax
        if (nnot .gt. 0) write (iunit,6269) 'REMARK',
     +    'Nr of residues not matched in mol2 : ',nnot
c
      else if (what(1:2) .eq. 'DD') then
        rms = 0.0
        ngt10 = 0
        ave = 0
        do i=1,nuse
          call fixdif (buffi(1,i),buffj(1,i),buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + abs(buffb2(i))
          if (abs(buffb2(i)) .gt. 10.0) ngt10 = ngt10 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-DIH',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb3(j) = buffb3(j) + 1.0
        end do
        write (iunit,6230) (buffb2(i),i=1,nuse)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt10) / float(nuse)
        call fvalut (' RMS delta DIH       :',1,rms)
        call fvalut (' Average |delta DIH| :',1,ave)
        call ivalut (' Nr |delta DIH| > 10 :',1,ngt10)
        call fvalut (' Percentage          :',1,pgt)
        write (iunit,6268) 'REMARK','RMS delta DIH = ',rms,
     +    ' +++ Average |delta DIH| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta DIH| > 10 = ',
     +    ngt10,' +++ Percentage = ',pgt
c
        write (iunit,6200) 'MORE  '
        write (iunit,6210) 'COLOUR',1
c
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_DANG R ',
     +    nuse,' (9f8.2)'
        call pretty (line)
        write (iunit,'(a)') line(1:leng1(line))
c
        write (iunit,6200) 'YVALUE','*'
c
        rms = 0.0
        ngt10 = 0
        ave = 0
        do i=1,nuse
          call fixdif (buffk(1,i),buffl(1,i),buffb2(i))
          buffb2(i) = abs(buffb2(i))
          rms = rms + buffb2(i)**2
          ave = ave + buffb2(i)
          if (buffb2(i) .gt. 5.0) ngt10 = ngt10 + 1
c
          if (abs(buffb2(i)) .gt. outcut) then
            write (*,6206) 'Delta-ANG',
     +        resnam(nint(buffm(1,i)),imol),
     +        achain(nint(buffm(1,i)),imol),
     +        iresid(nint(buffm(1,i)),imol),
     +        resnam(nint(buffm(2,i)),jmol),
     +        achain(nint(buffm(2,i)),jmol),
     +        iresid(nint(buffm(2,i)),jmol),
     +        buffb2(i)
          end if
c
          j = 1 + int(abs(buffb2(i))/bin)
          if (j .gt. nmax) j = nmax
          buffb4(j) = buffb4(j) + 1.0
        end do
        write (iunit,6230) (buffb2(i),i=1,nuse)
        rms = sqrt ( rms / float (nuse) )
        ave = ave / float(nuse)
        pgt = 100.0 * float(ngt10) / float(nuse)
        call fvalut (' RMS |delta ANG|     :',1,rms)
        call fvalut (' Average |delta ANG| :',1,ave)
        call ivalut (' Nr |delta ANG| > 5  :',1,ngt10)
        call fvalut (' Percentage          :',1,pgt)
        write (iunit,6268) 'REMARK','RMS |delta ANG| = ',rms,
     +    ' +++ Average |delta ANG| = ',ave
        write (iunit,6269) 'REMARK','Nr |delta ANG| > 5 = ',
     +    ngt10,' +++ Percentage = ',pgt
c
      end if
c
c re m1 /nfs/pdb/full/3sdp.pdb
c ph m1 a5-190 m1 b5 3sdp_phipsi.plt
c ex m1 a5-190 m1 b5
c di m1 a5-190 m1 b5 3sdp_cadist.plt
c dd m1 a5-190 m1 b5 3sdp_cadihe.plt
c
c re m2 /nfs/pdb/full/5rub.pdb
c ph m2 a2-457 m2 b2 5rub_phipsi.plt
c ex m2 a2-457 m2 b2
c di m2 a2-457 m2 b2 5rub_cadist.plt
c dd m2 a2-457 m2 b2 5rub_cadihe.plt
c
c re m3 /home/gerard/projects/lowres/comp/alex/alex.pdb
c ph m3 a1-221 m3 b1 alex_phipsi.plt
c ex m3 a1-221 m3 b1
c di m3 a1-221 m3 b1 alex_cadist.plt
c dd m3 a1-221 m3 b1 alex_cadihe.plt
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
c ... print histograms
c
 6020 format (1x,i8,' in [',f6.2,'-',f6.2,'> = ',f6.2,
     +  '% (Cumul ',f6.2,' %)')
c
      if (what(1:2) .eq. 'PH') then
        call prompt ('0Histogram of |delta PHI| values :')
      else if (what(1:2) .eq. 'ET') then
        call prompt ('0Histogram of |delta ETA| values :')
      else if (what(1:2) .eq. 'D1') then
        call prompt ('0Histogram of |delta Delta-1| values :')
      else if (what(1:2) .eq. 'DD') then
        call prompt ('0Histogram of |delta DIH| values :')
      else if (what(1:2) .eq. 'DI') then
        call prompt ('0Histogram of distances :')
      end if
c
      q4 = 0.0
      do i=1,nmax
        if (buffb3(i) .gt. 0.5) then
          q1 = float(i-1)*bin
          q2 = q1 + bin
          q3 = 100.0 * buffb3(i) / float(nuse)
          q4 = q4 + q3
          write (*,6020) nint(buffb3(i)),q1,q2,q3,q4
        end if
      end do
c
      if (what(1:2) .eq. 'PH') then
        call prompt ('0Histogram of |delta PSI| values :')
      else if (what(1:2) .eq. 'ET') then
        call prompt ('0Histogram of |delta THETA| values :')
      else if (what(1:2) .eq. 'D1') then
        call prompt ('0Histogram of |delta Delta-2| values :')
      else if (what(1:2) .eq. 'DD') then
        call prompt ('0Histogram of |delta ANG| values :')
      end if
c
      if (what(1:2) .ne. 'DI') then
        q4 = 0.0
        do i=1,nmax
          if (buffb4(i) .gt. 0.5) then
            q1 = float(i-1)*bin
            q2 = q1 + bin
            q3 = 100.0 * buffb4(i) / float(nuse)
            q4 = q4 + q3
            write (*,6020) nint(buffb4(i)),q1,q2,q3,q4
          end if
        end do
      end if
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine waters (imol,jmol,filnam)
c
      include 'lsqman.incl'
c
      integer maxbin
      parameter (maxbin = 25)
c
      real rtx(12),x,cw2,dmin,rms,rmb,abn,ab1,ab2
      real xave,xsdv,xmin,xmax,xtot,yave,ysdv,ymin,ymax,ytot
      real rr,ss,cc,q1,q2,q3,q4,q5,q6
c
      integer imol,jmol,i,j,ii,jj,nw1,nw2,nmat,nmin,ierr,leng1
c
      logical lwater,xinter
c
      character filnam*(*),line*256
c
code ...
c
      call textut (' Water comparison of mol 1 :',name(imol))
      call textut ('                 and mol 2 :',name(jmol))
c
c ... get waters in mol 1 (by residue name)
c
      nw1 = 0
      do i=1,natoms(imol)
        if (atmnam(i,imol)(1:2) .ne. ' O') goto 1001
        if (lwater(resnam(i,imol))) then
          nw1 = nw1 + 1
          buffi (1,nw1) = atmxyz (1,i,imol)
          buffi (2,nw1) = atmxyz (2,i,imol)
          buffi (3,nw1) = atmxyz (3,i,imol)
          buffb1 (nw1)  = batom (i,imol)
          ptri (nw1) = i
        end if
 1001   continue
      end do
      call jvalut (' Waters in mol 1 :',1,nw1)
      if (nw1 .lt. 3) then
        call errcon ('Too few waters')
        return
      end if
c
c ... get waters in mol 2 (by residue name)
c
      nw2 = 0
      do i=1,natoms(jmol)
        if (atmnam(i,jmol)(1:2) .ne. ' O') goto 1002
        if (lwater(resnam(i,jmol))) then
          nw2 = nw2 + 1
          buffk (1,nw2) = atmxyz (1,i,jmol)
          buffk (2,nw2) = atmxyz (2,i,jmol)
          buffk (3,nw2) = atmxyz (3,i,jmol)
          buffb2 (nw2)  = batom (i,jmol)
          ptrj (nw2) = i
        end if
 1002   continue
      end do
      call jvalut (' Waters in mol 2 :',1,nw2)
      if (nw2 .lt. 3) then
        call errcon ('Too few waters')
        return
      end if
c
c ... apply current operator to waters of mol 2
c
      do i=1,12
        rtx (i) = rtlsq (i,imol,jmol)
      end do
      call fvalut (' Applying current operator to mol 2 :',12,rtx)
      call vecrtv (buffk,buffj,nw2,rtx(1),rtx(10))
      write (*,*)
c
c ... match waters of mol 1
c
      cw2 = cutwat * cutwat
      nmat = 0
      rms = 0.0
      rmb = 0.0
      abn = 0.0
      ab1 = 0.0
      ab2 = 0.0
      do i=1,nw1
        dmin = cw2 + 0.001
        nmin = -1
        do j=1,nw2
          x = (buffi(1,i)-buffj(1,j))**2 +
     +        (buffi(2,i)-buffj(2,j))**2 +
     +        (buffi(3,i)-buffj(3,j))**2
          if (x .lt. dmin) then
            dmin = x
            nmin = j
          end if
        end do
        if (nmin .gt. 0) then
          nmat = nmat + 1
          rms = rms + dmin
          rmb = rmb + (buffb1(i)-buffb2(nmin))**2
          ab1 = ab1 + buffb1(i)
          ab2 = ab2 + buffb2(nmin)
          x = sqrt (dmin)
          buffb3 (nmat) = x
          buffb4 (nmat) = abs(buffb1(i)-buffb2(nmin))
          ii = ptri (i)
          jj = ptrj (nmin)
          write (*,6000) nmat,resnam(ii,imol),achain(ii,imol),
     +      iresid(ii,imol),resnam(jj,jmol),achain(jj,jmol),
     +      iresid(jj,jmol),x,batom(ii,imol),batom(jj,jmol)
        else
          abn = abn + buffb1(i)
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of matched waters   :',1,nmat)
      x = 100.0 * float(nmat) / float(nw1)
      call fvalut (' % matched waters mol 1 :',1,x)
      x = 100.0 * float(nmat) / float(nw2)
      call fvalut (' % matched waters mol 2 :',1,x)
      if (nmat .lt. nw1) then
        abn = abn / float (nw1 - nmat)
        call fvalut (
     +    ' Average B of non-matched waters in mol 1 :',1,abn)
      end if
      if (nmat .lt. 3) then
        call errcon ('Too few matched waters')
        return
      end if
      x = ab1 / float(nmat)
      call fvalut (
     +    ' Average B of matched waters in mol 1     :',1,x)
      x = ab2 / float(nmat)
      call fvalut (
     +    ' Average B of matched waters in mol 2     :',1,x)
c
      rms = sqrt (rms / float(nmat))
      rmb = sqrt (rmb / float(nmat))
      call fvalut (' RMS distance (A) :',1,rms)
      call fvalut (' RMS delta-B (A2) :',1,rmb)
c
      call xstats (buffb3,nmat,xave,xsdv,xmin,xmax,xtot)
      call xstats (buffb4,nmat,yave,ysdv,ymin,ymax,ytot)
      write (*,6010) 'Matching distances',xave,xsdv,xmin,xmax
      write (*,6010) 'Delta-B values',yave,ysdv,ymin,ymax
c
      call xystat (buffb3,buffb4,nmat,rr,ss,cc,
     +             q1,q2,q3,q4,q5,q6)
      call fvalut (' Correlation coefficient :',1,cc)
c
 6000 format (' ',i5,1x,a3,'-',a1,i4,' <-> ',a3,'-',a1,i4,
     +  ' | D ',f6.2,' | B ',2f7.2)
c
 6010 format (' ',A20,' :'/
     +  ' Average :',f8.2,'   St. dev. :',f8.2/
     +  ' Minimum :',f8.2,'   Maximum  :',f8.2)
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'LINFIT'
      write (iunit,6200) 'XLABEL',
     +  'Distance (A) between matched water molecules'
      write (iunit,6200) 'YLABEL',
     +  '|Delta-B| (A2) of matched water molecules'
      write (iunit,6200) 'REMARK',
     +  'Scatter plot of distance and abs(delta-B) of matched ',
     +  'water molecules in two models or chains'
c
      write (iunit,6200) 'REMARK',
     +  'Waters are those of mol 1 = ',name(imol)
     +  (1:leng1(name(imol))),' = file ',
     +  file(imol)(1:leng1(file(imol)))
      write (iunit,6200) 'REMARK',
     +  'Matched waters are from mol 2 = ',name(jmol)
     +  (1:leng1(name(jmol))),' = file ',
     +  file(jmol)(1:leng1(file(jmol)))
c
      write (iunit,6215) 'REMARK',
     +  'Distances ave/sdv/min/max = ',xave,xsdv,xmin,xmax
      write (iunit,6215) 'REMARK',
     +  '|Delta-B| ave/sdv/min/max = ',yave,ysdv,ymin,ymax
      write (iunit,6215) 'REMARK',
     +  'Correlation coefficient = ',cc
      write (iunit,6215) 'REMARK',
     +  'Cut-off distance used (A) = ',cutwat
      write (iunit,6217) 'REMARK',
     +  'Waters in mol 1, mol 2, matched ',nw1,nw2,nmat
c
      write (iunit,6213) 'XYVIEW',
     +    0.0,0.1*nint(10.0*(xmax+0.1)),0.0,float(nint(ymax+1.0))
c
      write (iunit,6210) 'NPOINT',nmat
      write (iunit,6210) 'COLOUR',4
c
      write (iunit,6200) 'XVALUE','*'
      write (iunit,6230) (buffb3(i),i=1,nmat)
c
      write (iunit,6200) 'YVALUE','*'
      write (iunit,6230) (buffb4(i),i=1,nmat)
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      return
      end
c
c
c
      subroutine closan (p1,p0)
c
c ... map angle P1 so it is closest to P0
c
      real p1,p0
c
code ...
c
      if (abs(p1-p0) .le. 180.0) return
c
      if (p1 .gt. p0) then
   10   p1 = p1 - 360.0
        if (abs(p1-p0) .le. 180.0) return
        goto 10
      else
   20   p1 = p1 + 360.0
        if (abs(p1-p0) .le. 180.0) return
        goto 20
      end if
c
      end
c
c
c
      subroutine muldih (imol,refchn,filnam,imode,how,cutoff,
     +                   ierr,prognm)
c
      include 'lsqman.incl'
c
      character*12 prognm
c
      real phi(0:maxchn),psi(0:maxchn)
      real dist,tangle
      real xave,xsdv,xmin,xmax,xtot,yave,ysdv,ymin,ymax,ytot,pmax
      real psize,cosphi,sinphi,cospsi,sinpsi,cvphi,cvpsi,qmax
      real sphiav,sphisd,sphimi,sphima,spsiav,spsisd,spsimi,spsima
      real rphiav,rphisd,rphimi,rphima,rpsiav,rpsisd,rpsimi,rpsima
      real cphiav,cphisd,cphimi,cphima,cpsiav,cpsisd,cpsimi,cpsima
      real cutoff
c
      integer imol,ierr,i,j,i1,i2,ir1,ir2,leng1
      integer nphi,npsi,iptr,jptr,ichn,ndone,mode,imode,nok,nlist
      integer inptr,jnptr,icptr,jcptr,ic1ptr,jc1ptr,in1ptr,jn1ptr
c
      logical xinter
c
      character refchn*(1),filnam*(*),jchn*1,kchn*1,how*1
      character line*256
c
code ...
c
      ierr = 0
c
      psize = 1.25
c
c ... mode = 1 -> SIGMA(phi), SIGMA(psi) plot
c            2 -> multiple Ramachandran plot
c            3 -> circular variance phi, psi plot
c
      mode = imode
c
      do i=1,nchain(imol)
        if (chname(i,imol) .eq. refchn) then
          ichn = i
          goto 10
        end if
      end do
      call errcon ('Reference chain not found')
      ierr = -1
      return
c
   10 continue
      call textut (' Reference chain :',chname(ichn,imol))
      kchn = chname(ichn,imol)
      i1 = chnptr(1,ichn,imol)
      i2 = chnptr(2,ichn,imol)
      ir1 = iresid(i1,imol)
      ir2 = ir1
      do i=i1,i2
        ir1 = min (ir1,iresid(i,imol))
        ir2 = max (ir2,iresid(i,imol))
      end do
      write (*,'(1x,a,i5,a,i5)') 'Residue range : ',ir1,' - ',ir2
c
      call fvalut (' Cut-off for printing :',1,cutoff)
c
      iptr = 1
      jptr = 1
c
      ndone = 0
      nlist = 0
c
c ... open PostScript file for multi Ramachandran plot
c
      if (mode .eq. 2) then
        if (how .ne. 'P') then
          call psrini (iunit,filnam,prognm)
        else
          call psrinp (iunit,filnam,prognm,how)
        end if
      end if
c
c ... loop over the residues
c
      write (*,6010) 'Residue','<phi>','S(phi)','Min','Max',
     +  '<psi>','S(psi)','Min','Max','#phi','#psi',
     +  'V(phi)','V(psi)'
c
      do i=ir1,ir2
c
        nphi = 0
        npsi = 0
        phi(1) = 0.0
        psi(1) = 0.0
        cosphi = 0.0
        sinphi = 0.0
        cospsi = 0.0
        sinpsi = 0.0
c
        call getptr (imol,kchn,i,' CA ',
     +               imol,kchn,i,iptr,jptr,ierr)
        if (ierr .ne. 0) goto 6896
c
        do j=1,nchain(imol)
          jchn = chname(j,imol)
          call getptr (imol,kchn,i,' CA ',
     +                 imol,jchn,i,iptr,jptr,ierr)
          if (ierr .ne. 0) goto 6996
c
          call getptr (imol,kchn,i,' N  ',
     +                 imol,jchn,i,inptr,jnptr,ierr)
          if (ierr .ne. 0) goto 6996
          call getptr (imol,kchn,i,' C  ',
     +                 imol,jchn,i,icptr,jcptr,ierr)
          if (ierr .ne. 0) goto 6996
c
c ... PHI ?
c
          call getptr (imol,kchn,i-1,' C  ',
     +                 imol,jchn,i-1,ic1ptr,jc1ptr,ierr)
          if (ierr .eq. 0) then
            if (dist(jc1ptr,jnptr,atmxyz(1,1,imol)) .le. 3.0) then
              nphi = nphi + 1
              phi (nphi) = tangle (jc1ptr,jnptr,jptr,jcptr,
     +                             atmxyz(1,1,imol))
c
              cosphi = cosphi + cos (degtor*phi(nphi))
              sinphi = sinphi + sin (degtor*phi(nphi))
c
              if (nphi .eq. 1) then
                phi (0) = tangle (ic1ptr,inptr,iptr,icptr,
     +                             atmxyz(1,1,imol))
              end if
ccc              print *,i,(phi(k),k=0,nphi)
              call closan (phi(nphi),phi(0))
            end if
          end if
c
c ... PSI ?
c
          call getptr (imol,kchn,i+1,' N  ',
     +                 imol,jchn,i+1,in1ptr,jn1ptr,ierr)
          if (ierr .eq. 0) then
            if (dist(jcptr,jn1ptr,atmxyz(1,1,imol)) .le. 3.0) then
              npsi = npsi + 1
              psi (npsi) = tangle (jnptr,jptr,jcptr,jn1ptr,
     +                             atmxyz(1,1,imol))
c
              cospsi = cospsi + cos (degtor*psi(npsi))
              sinpsi = sinpsi + sin (degtor*psi(npsi))
c
              if (npsi .eq. 1) then
                psi (0) = tangle (inptr,iptr,icptr,in1ptr,
     +                             atmxyz(1,1,imol))
              end if
              call closan (psi(npsi),psi(0))
            end if
          end if
c
 6996     continue
c
        end do
c
        xave = phi(1)
        yave = psi(1)
        xsdv = 0.0
        ysdv = 0.0
        xmin = phi(1)
        ymin = psi(1)
        xmax = phi(1)
        ymax = psi(1)
        cvphi = -1.0
        cvpsi = -1.0
c
        if (nphi .gt. 1) then
c
          cosphi = cosphi / float(nphi)
          sinphi = sinphi / float(nphi)
          xave = rtodeg * atan2 (sinphi,cosphi)
c
ccc          call xstats (phi(1),nphi,xave,xsdv,xmin,xmax,xtot)
 6910     continue
          cosphi = 0.0
          sinphi = 0.0
          do j=1,nphi
            call closan (phi(j),xave)
            cosphi = cosphi + cos (degtor*phi(j))
            sinphi = sinphi + sin (degtor*phi(j))
          end do
          call xstats (phi(1),nphi,xave,xsdv,xmin,xmax,xtot)
          cosphi = cosphi / float(nphi)
          sinphi = sinphi / float(nphi)
          xave = rtodeg * atan2 (sinphi,cosphi)
          call cirvar (phi(1),nphi,cvphi)
c
          if ( abs(xmax-xave) .gt. 180.0) goto 6910
          if ( abs(xmin-xave) .gt. 180.0) goto 6910
        end if
c
        if (npsi .gt. 1) then
c
          cospsi = cospsi / float(npsi)
          sinpsi = sinpsi / float(npsi)
          yave = rtodeg * atan2 (sinpsi,cospsi)
c
ccc          call xstats (psi(1),npsi,yave,ysdv,ymin,ymax,ytot)
 6920     continue
          cospsi = 0.0
          sinpsi = 0.0
          do j=1,npsi
            call closan (psi(j),yave)
            cospsi = cospsi + cos (degtor*psi(j))
            sinpsi = sinpsi + sin (degtor*psi(j))
          end do
          call xstats (psi(1),npsi,yave,ysdv,ymin,ymax,ytot)
          cospsi = cospsi / float(npsi)
          sinpsi = sinpsi / float(npsi)
          yave = rtodeg * atan2 (sinpsi,cospsi)
          call cirvar (psi(1),npsi,cvpsi)
c
          if ( abs(ymax-yave) .gt. 180.0) goto 6920
          if ( abs(ymin-yave) .gt. 180.0) goto 6920
        end if
c
        if ( (nphi+npsi) .gt. 0) then
c
          ndone = ndone + 1
          buffb1 (ndone)          = float(i)
          buffb1 (maxres+ndone)   = xsdv
          buffb1 (2*maxres+ndone) = abs (xmax-xmin)
          if (buffb1 (2*maxres+ndone) .gt. 360.0) then
            buffb1 (2*maxres+ndone) = buffb1 (2*maxres+ndone) - 360.0
          end if
          buffb1 (3*maxres+ndone) = ysdv
          buffb1 (4*maxres+ndone) = abs (ymax-ymin)
          if (buffb1 (4*maxres+ndone) .gt. 360.0) then
            buffb1 (4*maxres+ndone) = buffb1 (4*maxres+ndone) - 360.0
          end if
          buffb2 (ndone)          = cvphi
          buffb2 (maxres+ndone)   = cvpsi
c
          if (mode .eq. 1 .or. mode .eq. 2) then
            if (buffb1(2*maxres+ndone) .ge. cutoff .or.
     +          buffb1(4*maxres+ndone) .ge. cutoff) then
              write (*,6000) resnam(iptr,imol),i,xave,xsdv,
     +          xmin,xmax,yave,ysdv,ymin,ymax,nphi,npsi,cvphi,cvpsi
              nlist = nlist + 1
            end if
          else
            if (cvphi .ge. cutoff .or.
     +          cvpsi .ge. cutoff) then
              write (*,6000) resnam(iptr,imol),i,xave,xsdv,
     +          xmin,xmax,yave,ysdv,ymin,ymax,nphi,npsi,cvphi,cvpsi
              nlist = nlist + 1
            end if
          end if
c
c ... add point to multi Ramachandran plot
c
          if (mode .eq. 2) then
            if (nphi.gt.0 .and. npsi.gt.0 .and. nphi.eq.npsi) then
              if (how .eq. 'C') then
                call fixang (xave)
                call fixang (yave)
                call mmfold (xave,yave,nphi,phi(1),psi(1),1)
                do j=1,nphi
                  call fixang (phi(j))
                  call fixang (psi(j))
                end do
              else
                call fix360 (xave)
                call fix360 (yave)
                do j=1,nphi
                  call fix360 (phi(j))
                  call fix360 (psi(j))
                end do
              end if
              do j=1,nphi
                if (resnam(iptr,imol) .eq. 'GLY') then
                  call xps_move (phi(j)-psize,psi(j)-psize)
                  call xps_draw (phi(j)+psize,psi(j)-psize)
                  call xps_draw (phi(j)+psize,psi(j)+psize)
                  call xps_draw (phi(j)-psize,psi(j)+psize)
                  call xps_draw (phi(j)-psize,psi(j)-psize)
                else
                  call xps_move (phi(j)-psize,psi(j))
                  call xps_draw (phi(j)+psize,psi(j))
                  call xps_move (phi(j),psi(j)-psize)
                  call xps_draw (phi(j),psi(j)+psize)
                end if
                if (how .ne. 'C') then
                  call xps_move (xave,yave)
                  call xps_draw (phi(j),psi(j))
                end if
              end do
            end if
          end if
c
        end if
c
 6896   continue
c
      end do
c
      call jvalut (' Nr of residues found :',1,ndone)
      call jvalut (' Nr of residues shown :',1,nlist)
      if (ndone .lt. 3) return
c
      call xstats (buffb1(1),ndone,xave,xsdv,xmin,xmax,xtot)
      ir1 = nint(xmin)
      ir2 = nint(xmax)
c
 6110 format (1x,a,' Ave, Sdv, Min, Max, # : ',4f8.2,i8)
 6120 format (1x,a,' Ave, Sdv, Min, Max, # : ',4f8.3,i8)
c
      pmax = 90.0
      qmax = 0.0
c
      nok = 0
      do i=maxres+1,maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,sphiav,sphisd,sphimi,sphima,xtot)
      write (*,6110) 'SIGMA(phi)',sphiav,sphisd,sphimi,sphima,nok
      if (mode .eq. 2) then
        write (line,6110) 'SIGMA(phi)',sphiav,sphisd,sphimi,sphima,nok
        call xps_legend (line)
      end if
      pmax = max (pmax,sphima)
c
      nok = 0
      do i=2*maxres+1,2*maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,rphiav,rphisd,rphimi,rphima,xtot)
      write (*,6110) 'RANGE(phi)',rphiav,rphisd,rphimi,rphima,nok
      if (mode .eq. 2) then
        write (line,6110) 'RANGE(phi)',rphiav,rphisd,rphimi,rphima,nok
        call xps_legend (line)
      end if
c
      nok = 0
      do i=3*maxres+1,3*maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,spsiav,spsisd,spsimi,spsima,xtot)
      write (*,6110) 'SIGMA(psi)',spsiav,spsisd,spsimi,spsima,nok
      if (mode .eq. 2) then
        write (line,6110) 'SIGMA(psi)',spsiav,spsisd,spsimi,spsima,nok
        call xps_legend (line)
      end if
      pmax = max (pmax,spsima)
c
      nok = 0
      do i=4*maxres+1,4*maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,rpsiav,rpsisd,rpsimi,rpsima,xtot)
      write (*,6110) 'RANGE(psi)',rpsiav,rpsisd,rpsimi,rpsima,nok
      if (mode .eq. 2) then
        write (line,6110) 'RANGE(psi)',rpsiav,rpsisd,rpsimi,rpsima,nok
        call xps_legend (line)
      end if
c
      nok = 0
      do i=1,ndone
        if (buffb2(i) .ge. 0.0) then
          nok = nok + 1
          buffb3(nok) = buffb2(i)
        else
          buffb2(i) = 0.0
        end if
      end do
      call xstats (buffb3,nok,cphiav,cphisd,cphimi,cphima,xtot)
      write (*,6120) 'CV(phi)   ',cphiav,cphisd,cphimi,cphima,nok
      if (mode .eq. 2) then
        write (line,6110) 'CV(phi)',cphiav,cphisd,cphimi,cphima,nok
        call xps_legend (line)
      end if
      qmax = max (qmax,cphima)
c
      nok = 0
      do i=maxres+1,maxres+ndone
        if (buffb2(i) .ge. 0.0) then
          nok = nok + 1
          buffb3(nok) = buffb2(i)
        else
          buffb2(i) = 0.0
        end if
      end do
      call xstats (buffb3,nok,cpsiav,cpsisd,cpsimi,cpsima,xtot)
      write (*,6120) 'CV(psi)   ',cpsiav,cpsisd,cpsimi,cpsima,nok
      if (mode .eq. 2) then
        write (line,6110) 'CV(psi)',cpsiav,cpsisd,cpsimi,cpsima,nok
        call xps_legend (line)
      end if
      qmax = max (qmax,cpsima)
c
c ... done if multi Rama
c
      if (mode .eq. 2) then
c
        write (line,6200) 'Mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        call pretty (line)
        call xps_legend (line)
c
        call xps_close ()
        close (iunit)
        call prompt (' PostScript file written')
        return
      end if
c
c ... otherwise, open plot file
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
 6240 format (9f8.3)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Residue number'
c
      if (mode .eq. 1) then
        write (iunit,6200) 'YLABEL',
     +      'SIGMA(phi) and SIGMA(psi)'
        write (iunit,6200) 'REMARK',
     +      'Plot of SIGMA(phi) (solid blue) and SIGMA(psi) ',
     +      '(dotted red) as a function of residue nr'
      else if (mode .eq. 3) then
        write (iunit,6200) 'YLABEL',
     +      'CV(phi) and CV(psi)'
        write (iunit,6200) 'REMARK',
     +      'Plot of CV(phi) (solid blue) and CV(psi) ',
     +      '(dotted red) as a function of residue nr'
      end if
      write (iunit,6200) 'REMARK',
     +    'Values are for mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
      if (mode .eq. 1) then
        write (iunit,6213) 'XYVIEW',
     +    float(ir1-1),float(ir2+1),0.0,pmax+1.0
      else if (mode .eq. 3) then
        write (iunit,6213) 'XYVIEW',
     +    float(ir1-1),float(ir2+1),0.0,1.0
ccc     +    float(ir1-1),float(ir2+1),0.0,qmax+0.01
      end if
c
      write (iunit,6210) 'NPOINT',ndone
      write (iunit,6210) 'COLOUR',4
c
      write (iunit,6200) 'XVALUE','*'
      write (iunit,6220) (nint(buffb1(i)),i=1,ndone)
      write (iunit,6215) 'REMARK',
     +    'SIGMA(phi) ave, sdv, min, max ',
     +    sphiav,sphisd,sphimi,sphima
      write (iunit,6215) 'REMARK',
     +    'RANGE(phi) ave, sdv, min, max ',
     +    rphiav,rphisd,rphimi,rphima
      write (iunit,6215) 'REMARK',
     +    'CV(phi) ave, sdv, min, max ',
     +    cphiav,cphisd,cphimi,cphima
c
      if (mode .eq. 1) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_SIGPHI R ',
     +    ndone,' (9f8.2)'
      else if (mode .eq. 3) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_CVPHI R ',
     +    ndone,' (9f8.3)'
      end if
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
      write (iunit,6200) 'YVALUE','*'
c
      if (mode .eq. 1) then
        write (iunit,6230) (buffb1(maxres+i),i=1,ndone)
      else if (mode .eq. 3) then
        write (iunit,6240) (buffb2(i),i=1,ndone)
      end if
c
      write (iunit,6200) 'MORE  '
      write (iunit,6210) 'COLOUR',1
c
      write (iunit,6215) 'REMARK','SIGMA(psi) ave, sdv, min, max ',
     +  spsiav,spsisd,spsimi,spsima
      write (iunit,6215) 'REMARK','RANGE(psi) ave, sdv, min, max ',
     +  rpsiav,rpsisd,rpsimi,rpsima
      write (iunit,6215) 'REMARK','CV(psi) ave, sdv, min, max ',
     +  cpsiav,cpsisd,cpsimi,cpsima
c
      if (mode .eq. 1) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_SIGPSI R ',
     +    ndone,' (9f8.2)'
      else if (mode .eq. 3) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_CVPSI R ',
     +    ndone,' (9f8.3)'
      end if
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      write (iunit,6200) 'YVALUE','*'
c
      if (mode .eq. 1) then
        write (iunit,6230) (buffb1(3*maxres+i),i=1,ndone)
      else if (mode .eq. 3) then
        write (iunit,6240) (buffb2(maxres+i),i=1,ndone)
      end if
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      ierr = 0
      return
c
 6000 format (1x,a3,i6,' | ',2(4f7.1,:,' | '),2i5,' | ',2f7.3)
 6010 format (1x,a9,3x,2(4a7,3x),2a5,3x,2a7)
c
      end
c
c
c
      subroutine mulbfs (imol,refchn,filnam,cutoff,ierr)
c
      include 'lsqman.incl'
c
      real phi(0:maxchn)
      real xave,xsdv,xmin,xmax,xtot,pmax,cutoff
c
      integer imol,ierr,i,j,i1,i2,ir1,ir2
      integer nphi,iptr,jptr,ichn,ndone,leng1,nlist
c
      logical xinter
c
      character refchn*(1),filnam*(*),jchn*1,kchn*1
      character line*256
c
code ...
c
      ierr = 0
c
      do i=1,nchain(imol)
        if (chname(i,imol) .eq. refchn) then
          ichn = i
          goto 10
        end if
      end do
      call errcon ('Reference chain not found')
      ierr = -1
      return
c
   10 continue
      call textut (' Reference chain :',chname(ichn,imol))
      kchn = chname(ichn,imol)
      i1 = chnptr(1,ichn,imol)
      i2 = chnptr(2,ichn,imol)
      ir1 = iresid(i1,imol)
      ir2 = ir1
      do i=i1,i2
        ir1 = min (ir1,iresid(i,imol))
        ir2 = max (ir2,iresid(i,imol))
      end do
      write (*,'(1x,a,i5,a,i5)') 'Residue range : ',ir1,' - ',ir2
      call textut (' Central atom type :',atypes(1))
c
      if (atypes(1) .eq. 'ALL ' .or. atypes(1) .eq. 'NONH' .or.
     +    atypes(1) .eq. 'SIDE' .or. atypes(1) .eq. 'TRAC') then
        call errcon ('Invalid atom type !')
        ierr = -1
        return
      end if
c
      call fvalut (' Cut-off for printing :',1,cutoff)
c
      iptr = 1
      jptr = 1
c
      ndone = 0
      nlist = 0
c
c ... loop over the residues
c
      do i=ir1,ir2
c
        nphi = 0
        phi(1) = 0.0
c
        do j=1,nchain(imol)
          jchn = chname(j,imol)
          call getptr (imol,kchn,i,atypes(1),
     +                 imol,jchn,i,iptr,jptr,ierr)
          if (ierr .ne. 0) goto 6996
c
          nphi = nphi + 1
          phi (nphi) = batom(jptr,imol)
c
 6996     continue
c
        end do
c
        xave = phi(1)
        xsdv = 0.0
        xmin = phi(1)
        xmax = phi(1)
c
        if (nphi .gt. 1) then
          call xstats (phi(1),nphi,xave,xsdv,xmin,xmax,xtot)
        end if
c
        if (nphi .gt. 0) then
c
          ndone = ndone + 1
          buffb1 (ndone)          = float(i)
          buffb1 (maxres+ndone)   = xsdv
          buffb1 (2*maxres+ndone) = abs (xmax-xmin)
c
          if (buffb1(2*maxres+ndone) .ge. cutoff) then
            write (*,6000) resnam(iptr,imol),i,xave,xsdv,xmin,xmax,nphi
            nlist = nlist + 1
          end if
c
        end if
c
      end do
c
      call jvalut (' Nr of residues found :',1,ndone)
      call jvalut (' Nr of residues shown :',1,nlist)
      if (ndone .lt. 3) return
c
      call xstats (buffb1(1),ndone,xave,xsdv,xmin,xmax,xtot)
      ir1 = nint(xmin)
      ir2 = nint(xmax)
c
 6110 format (1x,a,' Ave, Sdv, Min, Max : ',4f8.2)
c
      pmax = 10.0
c
      call xstats (buffb1(maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
      write (*,6110) 'SIGMA(B)',xave,xsdv,xmin,xmax
      pmax = max (pmax,xmax)
      call xstats (buffb1(2*maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
      write (*,6110) 'RANGE(B)',xave,xsdv,xmin,xmax
      pmax = max (pmax,xmax)
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Residue number'
c
      write (iunit,6200) 'YLABEL',
     +    'SIGMA(B) and RANGE(B)'
      write (iunit,6200) 'REMARK',
     +    'Plot of SIGMA(B) (solid blue) and RANGE(B) ',
     +    '(dotted red) as a function of residue nr'
      write (iunit,6200) 'REMARK',
     +    'Atom type used ',atypes(1)
      write (iunit,6200) 'REMARK',
     +    'Values are for mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
      write (iunit,6213) 'XYVIEW',
     +    float(ir1-1),float(ir2+1),0.0,pmax+1.0
c
      write (iunit,6210) 'NPOINT',ndone
      write (iunit,6210) 'COLOUR',4
c
      write (iunit,6200) 'XVALUE','*'
      write (iunit,6220) (nint(buffb1(i)),i=1,ndone)
c
      call xstats (buffb1(maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
      write (iunit,6215) 'REMARK','SIGMA(B) ave, sdv, min, max ',
     +  xave,xsdv,xmin,xmax
c
      write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_SIGB R ',
     +    ndone,' (9f8.2)'
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
      write (iunit,6200) 'YVALUE','*'
c
      write (iunit,6230) (buffb1(maxres+i),i=1,ndone)
c
      write (iunit,6200) 'MORE  '
      write (iunit,6210) 'COLOUR',1
c
      call xstats (buffb1(2*maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
      write (iunit,6215) 'REMARK','RANGE(B) ave, sdv, min, max ',
     +  xave,xsdv,xmin,xmax
c
      write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_RANGEB R ',
     +  ndone,' (9f8.2)'
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      write (iunit,6200) 'YVALUE','*'
c
      write (iunit,6230) (buffb1(2*maxres+i),i=1,ndone)
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      ierr = 0
      return
c
 6000 format (1x,a3,1x,i5,' | ',1(4f8.1,:,' | '),2i3)
c
      end
c
c
c
      subroutine pstini (lunhp,psfile,how)
c
      implicit none
c
      integer lunhp,i,j
c
      real phi,psi
c
      character labx*40,laby*40,psfile*(*),how*1
c
code ...
c
      call xps_init ()
      call xps_open (lunhp,psfile,'LSQMAN')
c
      if (how .eq. 'P') then
        call xps_polar (0.0,360.0)
c
      else
        call xps_scale (0.,360.0,0.0,360.0)
c
        call xps_move (   0.,   0.)
        call xps_draw ( 360.,   0.)
        call xps_draw ( 360., 360.)
        call xps_draw (   0., 360.)
        call xps_draw (   0.,   0.)
      end if
c
      call xps_dash ()
      call xps_dash ()
c
c ---	chi1=60,180,-60
c
      if (how .eq. 'P') then
        do i=1,3
          phi = 60. + (i-1)*120.
          call xps_move (phi,0.0)
          do j=0,360,10
            call xps_draw (phi,float(j))
          end do
        end do
      else
        do 300 i=1,3
          phi = 60. + (i-1)*120.
          call xps_move (phi,  0.)
          call xps_draw (phi,359.9)
300     continue
      end if
c
c ---	chi2=60,180,-60
c
      do 310 i=1,3
        psi = 60. + (i-1)*120.
        call xps_move (  0.,psi)
        call xps_draw (359.9,psi)
310   continue
c
      call xps_solid ()
c
c ---	text
c
      labx = 'CHI1 mapped to [0,360>'
      laby = 'CHI2 mapped to [0,360>'
      call xps_label (labx,laby)
c
      return
      end
c
c
c
      subroutine mulsid (imol,refchn,filnam,imode,how,cutoff,ierr)
c
      include 'lsqman.incl'
c
      real phi(0:maxchn),psi(0:maxchn)
      real tangle
      real xave,xsdv,xmin,xmax,xtot,yave,ysdv,ymin,ymax,ytot,pmax
      real psize,cosphi,sinphi,cospsi,sinpsi,cvphi,cvpsi,qmax
      real sphiav,sphisd,sphimi,sphima,spsiav,spsisd,spsimi,spsima
      real rphiav,rphisd,rphimi,rphima,rpsiav,rpsisd,rpsimi,rpsima
      real cphiav,cphisd,cphimi,cphima,cpsiav,cpsisd,cpsimi,cpsima
      real cutoff
c
      integer imol,ierr,i,j,i1,i2,ir1,ir2,leng1
      integer nphi,npsi,iptr,jptr,ichn,ndone,mode,imode,nok,nlist
      integer inptr,jnptr,icptr,jcptr,ic1ptr,jc1ptr,in1ptr,jn1ptr
c
      logical xinter
c
      character refchn*(1),filnam*(*),jchn*1,kchn*1,how*1
      character line*256
c
code ...
c
      ierr = 0
c
      psize = 1.25
c
c ... mode = 1 -> SIGMA(chi1), SIGMA(chi2) plot
c            2 -> multiple chi1/chi2 plot
c            3 -> circular variance chi1/chi2 plot
c
      mode = imode
c
      do i=1,nchain(imol)
        if (chname(i,imol) .eq. refchn) then
          ichn = i
          goto 10
        end if
      end do
      call errcon ('Reference chain not found')
      ierr = -1
      return
c
   10 continue
      call textut (' Reference chain :',chname(ichn,imol))
      kchn = chname(ichn,imol)
      i1 = chnptr(1,ichn,imol)
      i2 = chnptr(2,ichn,imol)
      ir1 = iresid(i1,imol)
      ir2 = ir1
      do i=i1,i2
        ir1 = min (ir1,iresid(i,imol))
        ir2 = max (ir2,iresid(i,imol))
      end do
      write (*,'(1x,a,i5,a,i5)') 'Residue range : ',ir1,' - ',ir2
c
      call fvalut (' Cut-off for printing :',1,cutoff)
c
      iptr = 1
      jptr = 1
c
      ndone = 0
      nlist = 0
c
c ... open PostScript file for multi distribution plot
c
      if (mode .eq. 2) then
        call pstini (iunit,filnam,how)
      end if
c
c ... loop over the residues
c
      write (*,6010) 'Residue','<c1>','S(c1)','Min','Max',
     +  '<c2>','S(c2)','Min','Max','#c1','#c2',
     +  'V(c1)','V(c2)'
c
      do i=ir1,ir2
c
        nphi = 0
        npsi = 0
        phi(1) = 0.0
        psi(1) = 0.0
        cosphi = 0.0
        sinphi = 0.0
        cospsi = 0.0
        sinpsi = 0.0
c
        call getptr (imol,kchn,i,' CA ',
     +               imol,kchn,i,iptr,jptr,ierr)
        if (ierr .ne. 0) goto 6896
c
        do j=1,nchain(imol)
          jchn = chname(j,imol)
          call getptr (imol,kchn,i,' CA ',
     +                 imol,jchn,i,iptr,jptr,ierr)
          if (ierr .ne. 0) goto 6996
c
          call getptr (imol,kchn,i,' N  ',
     +                 imol,jchn,i,inptr,jnptr,ierr)
          if (ierr .ne. 0) goto 6996
          call getptr (imol,kchn,i,' CB ',
     +                 imol,jchn,i,icptr,jcptr,ierr)
          if (ierr .ne. 0) goto 6996
c
c ... CHI1 ?
c
          call getptr (imol,kchn,i,' CG ',
     +                 imol,jchn,i,ic1ptr,jc1ptr,ierr)
          if (ierr .eq. 0) goto 2150
c
          call getptr (imol,kchn,i,' SG ',
     +                 imol,jchn,i,ic1ptr,jc1ptr,ierr)
          if (ierr .eq. 0) goto 2150
c
          call getptr (imol,kchn,i,' OG ',
     +                 imol,jchn,i,ic1ptr,jc1ptr,ierr)
          if (ierr .eq. 0) goto 2150
c
          call getptr (imol,kchn,i,' CG1',
     +                 imol,jchn,i,ic1ptr,jc1ptr,ierr)
          if (ierr .eq. 0) goto 2150
c
          call getptr (imol,kchn,i,' OG1',
     +                 imol,jchn,i,ic1ptr,jc1ptr,ierr)
          if (ierr .eq. 0) goto 2150
c
          goto 6996
c
 2150     continue
          nphi = nphi + 1
          phi (nphi) = tangle (jnptr,jptr,jcptr,jc1ptr,
     +                         atmxyz(1,1,imol))
c
          cosphi = cosphi + cos (degtor*phi(nphi))
          sinphi = sinphi + sin (degtor*phi(nphi))
c
          if (nphi .eq. 1) then
            phi (0) = tangle (inptr,iptr,icptr,ic1ptr,
     +                        atmxyz(1,1,imol))
          end if
ccc              print *,i,(phi(k),k=0,nphi)
          call closan (phi(nphi),phi(0))
c
c ... CHI2 ?
c
          call getptr (imol,kchn,i,' CD ',
     +                 imol,jchn,i,in1ptr,jn1ptr,ierr)
          if (ierr .eq. 0) goto 2160
c
          call getptr (imol,kchn,i,' SD ',
     +                 imol,jchn,i,in1ptr,jn1ptr,ierr)
          if (ierr .eq. 0) goto 2160
c
          call getptr (imol,kchn,i,' CD1',
     +                 imol,jchn,i,in1ptr,jn1ptr,ierr)
          if (ierr .eq. 0) goto 2160
c
          call getptr (imol,kchn,i,' OD1',
     +                 imol,jchn,i,in1ptr,jn1ptr,ierr)
          if (ierr .eq. 0) goto 2160
c
          call getptr (imol,kchn,i,' ND1',
     +                 imol,jchn,i,in1ptr,jn1ptr,ierr)
          if (ierr .eq. 0) goto 2160
c
          goto 6996
c
 2160     continue
c
          npsi = npsi + 1
          psi (npsi) = tangle (jptr,jcptr,jc1ptr,jn1ptr,
     +                         atmxyz(1,1,imol))
c
          cospsi = cospsi + cos (degtor*psi(npsi))
          sinpsi = sinpsi + sin (degtor*psi(npsi))
c
          if (npsi .eq. 1) then
            psi (0) = tangle (iptr,icptr,ic1ptr,in1ptr,
     +                        atmxyz(1,1,imol))
          end if
          call closan (psi(npsi),psi(0))
c
 6996     continue
c
        end do
c
        xave = phi(1)
        yave = psi(1)
        xsdv = 0.0
        ysdv = 0.0
        xmin = phi(1)
        ymin = psi(1)
        xmax = phi(1)
        ymax = psi(1)
        cvphi = -1.0
        cvpsi = -1.0
c
        if (nphi .gt. 1) then
c
          cosphi = cosphi / float(nphi)
          sinphi = sinphi / float(nphi)
          xave = rtodeg * atan2 (sinphi,cosphi)
c
ccc          call xstats (phi(1),nphi,xave,xsdv,xmin,xmax,xtot)
 6910     continue
          cosphi = 0.0
          sinphi = 0.0
          do j=1,nphi
            call closan (phi(j),xave)
            cosphi = cosphi + cos (degtor*phi(j))
            sinphi = sinphi + sin (degtor*phi(j))
          end do
          call xstats (phi(1),nphi,xave,xsdv,xmin,xmax,xtot)
          cosphi = cosphi / float(nphi)
          sinphi = sinphi / float(nphi)
          xave = rtodeg * atan2 (sinphi,cosphi)
          call cirvar (phi(1),nphi,cvphi)
c
          if ( abs(xmax-xave) .gt. 180.0) goto 6910
          if ( abs(xmin-xave) .gt. 180.0) goto 6910
        end if
c
        if (npsi .gt. 1) then
c
          cospsi = cospsi / float(npsi)
          sinpsi = sinpsi / float(npsi)
          yave = rtodeg * atan2 (sinpsi,cospsi)
c
ccc          call xstats (psi(1),npsi,yave,ysdv,ymin,ymax,ytot)
 6920     continue
          cospsi = 0.0
          sinpsi = 0.0
          do j=1,npsi
            call closan (psi(j),yave)
            cospsi = cospsi + cos (degtor*psi(j))
            sinpsi = sinpsi + sin (degtor*psi(j))
          end do
          call xstats (psi(1),npsi,yave,ysdv,ymin,ymax,ytot)
          cospsi = cospsi / float(npsi)
          sinpsi = sinpsi / float(npsi)
          yave = rtodeg * atan2 (sinpsi,cospsi)
          call cirvar (psi(1),npsi,cvpsi)
c
          if ( abs(ymax-yave) .gt. 180.0) goto 6920
          if ( abs(ymin-yave) .gt. 180.0) goto 6920
        end if
c
        if ( (nphi+npsi) .gt. 0) then
c
          ndone = ndone + 1
          buffb1 (ndone)          = float(i)
          buffb1 (maxres+ndone)   = xsdv
          buffb1 (2*maxres+ndone) = abs (xmax-xmin)
          if (buffb1 (2*maxres+ndone) .gt. 360.0) then
            buffb1 (2*maxres+ndone) = buffb1 (2*maxres+ndone) - 360.0
          end if
          buffb1 (3*maxres+ndone) = ysdv
          buffb1 (4*maxres+ndone) = abs (ymax-ymin)
          if (buffb1 (4*maxres+ndone) .gt. 360.0) then
            buffb1 (4*maxres+ndone) = buffb1 (4*maxres+ndone) - 360.0
          end if
          buffb2 (ndone)          = cvphi
          buffb2 (maxres+ndone)   = cvpsi
c
          if (mode .eq. 1 .or. mode .eq. 2) then
            if (buffb1(2*maxres+ndone) .ge. cutoff .or.
     +          buffb1(4*maxres+ndone) .ge. cutoff) then
              write (*,6000) resnam(iptr,imol),i,xave,xsdv,xmin,xmax,
     +          yave,ysdv,ymin,ymax,nphi,npsi,cvphi,cvpsi
              nlist = nlist + 1
            end if
          else
            if (cvphi .ge. cutoff .or.
     +          cvpsi .ge. cutoff) then
              write (*,6000) resnam(iptr,imol),i,xave,xsdv,xmin,xmax,
     +          yave,ysdv,ymin,ymax,nphi,npsi,cvphi,cvpsi
              nlist = nlist + 1
            end if
          end if
c
c ... add point to multi chi1/chi2 plot
c
          if (mode .eq. 2) then
            if (nphi.gt.0 .and. npsi.gt.0 .and. nphi.eq.npsi) then
              call fix360 (xave)
              call fix360 (yave)
              if (how .eq.  'C') then
                call mmfold (xave,yave,nphi,phi(1),psi(1),2)
              end if
              do j=1,nphi
                call fix360 (phi(j))
                call fix360 (psi(j))
                call xps_move (phi(j)-psize,psi(j))
                call xps_draw (phi(j)+psize,psi(j))
                call xps_move (phi(j),psi(j)-psize)
                call xps_draw (phi(j),psi(j)+psize)
                if (how .ne. 'C') then
                  call xps_move (xave,yave)
                  call xps_draw (phi(j),psi(j))
                end if
              end do
            end if
          end if
c
        else
c
          ndone = ndone + 1
          buffb1 (ndone)          = float(i)
          buffb1 (maxres+ndone)   = xsdv
          buffb1 (2*maxres+ndone) = abs (xmax-xmin)
          buffb1 (3*maxres+ndone) = ysdv
          buffb1 (4*maxres+ndone) = abs (ymax-ymin)
          buffb2 (ndone)          = cvphi
          buffb2 (maxres+ndone)   = cvpsi
c
        end if
c
 6896   continue
c
      end do
c
      call jvalut (' Nr of residues found :',1,ndone)
      call jvalut (' Nr of residues shown :',1,nlist)
      if (ndone .lt. 3) return
c
      call xstats (buffb1(1),ndone,xave,xsdv,xmin,xmax,xtot)
      ir1 = nint(xmin)
      ir2 = nint(xmax)
c
 6110 format (1x,a,' Ave, Sdv, Min, Max, # : ',4f8.2,i8)
 6120 format (1x,a,' Ave, Sdv, Min, Max, # : ',4f8.3,i8)
c
      pmax = 90.0
c
      nok = 0
      do i=maxres+1,maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,sphiav,sphisd,sphimi,sphima,xtot)
      write (*,6110) 'SIGMA(chi1)',sphiav,sphisd,sphimi,sphima,nok
      if (mode .eq. 2) then
        write (line,6110) 'SIGMA(chi1)',sphiav,sphisd,sphimi,sphima,nok
        call xps_legend (line)
      end if
      pmax = max (pmax,sphima)
c
      nok = 0
      do i=2*maxres+1,2*maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,rphiav,rphisd,rphimi,rphima,xtot)
      write (*,6110) 'RANGE(chi1)',rphiav,rphisd,rphimi,rphima,nok
      if (mode .eq. 2) then
        write (line,6110) 'RANGE(chi1)',rphiav,rphisd,rphimi,rphima,nok
        call xps_legend (line)
      end if
c
      nok = 0
      do i=3*maxres+1,3*maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,spsiav,spsisd,spsimi,spsima,xtot)
      write (*,6110) 'SIGMA(chi2)',spsiav,spsisd,spsimi,spsima,nok
      if (mode .eq. 2) then
        write (line,6110) 'SIGMA(chi2)',spsiav,spsisd,spsimi,spsima,nok
        call xps_legend (line)
      end if
      pmax = max (pmax,spsima)
c
      nok = 0
      do i=4*maxres+1,4*maxres+ndone
        if (abs(buffb1(i)) .gt. 0.001) then
          nok = nok + 1
          buffb3(nok) = buffb1(i)
        end if
      end do
      call xstats (buffb3,nok,rpsiav,rpsisd,rpsimi,rpsima,xtot)
      write (*,6110) 'RANGE(chi2)',rpsiav,rpsisd,rpsimi,rpsima,nok
      if (mode .eq. 2) then
        write (line,6110) 'RANGE(chi2)',rpsiav,rpsisd,rpsimi,rpsima,nok
        call xps_legend (line)
      end if
c
      nok = 0
      do i=1,ndone
        if (buffb2(i) .ge. 0.0) then
          nok = nok + 1
          buffb3(nok) = buffb2(i)
        else
          buffb2(i) = 0.0
        end if
      end do
      call xstats (buffb3,nok,cphiav,cphisd,cphimi,cphima,xtot)
      write (*,6120) 'CV(chi1)   ',cphiav,cphisd,cphimi,cphima,nok
      if (mode .eq. 2) then
        write (line,6110) 'CV(chi1)',cphiav,cphisd,cphimi,cphima,nok
        call xps_legend (line)
      end if
      qmax = max (qmax,cphima)
c
      nok = 0
      do i=maxres+1,maxres+ndone
        if (buffb2(i) .ge. 0.0) then
          nok = nok + 1
          buffb3(nok) = buffb2(i)
        else
          buffb2(i) = 0.0
        end if
      end do
      call xstats (buffb3,nok,cpsiav,cpsisd,cpsimi,cpsima,xtot)
      write (*,6120) 'CV(chi2)   ',cpsiav,cpsisd,cpsimi,cpsima,nok
      if (mode .eq. 2) then
        write (line,6110) 'CV(chi2)',cpsiav,cpsisd,cpsimi,cpsima,nok
        call xps_legend (line)
      end if
      qmax = max (qmax,cpsima)
c
c ... done if multi chi1/2 plot
c
      if (mode .eq. 2) then
c
        write (line,6200) 'Mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        call pretty (line)
        call xps_legend (line)
c
        call xps_close ()
        close (iunit)
        call prompt (' PostScript file written')
        return
      end if
c
c ... otherwise, open plot file
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
 6240 format (9f8.3)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Residue number'
c
      if (mode .eq. 1) then
        write (iunit,6200) 'YLABEL',
     +    'SIGMA(chi1) and SIGMA(chi2)'
        write (iunit,6200) 'REMARK',
     +    'Plot of SIGMA(chi1) (solid blue) and SIGMA(chi2) ',
     +    '(dotted red) as a function of residue nr'
      else if (mode .eq. 3) then
        write (iunit,6200) 'YLABEL',
     +      'CV(chi1) and CV(chi2)'
        write (iunit,6200) 'REMARK',
     +      'Plot of CV(chi1) (solid blue) and CV(chi2) ',
     +      '(dotted red) as a function of residue nr'
      end if
      write (iunit,6200) 'REMARK',
     +    'Values are for mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
      if (mode .eq. 1) then
        write (iunit,6213) 'XYVIEW',
     +    float(ir1-1),float(ir2+1),0.0,pmax+1.0
      else if (mode .eq. 3) then
        write (iunit,6213) 'XYVIEW',
     +    float(ir1-1),float(ir2+1),0.0,1.0
ccc     +    float(ir1-1),float(ir2+1),0.0,qmax+0.01
      end if
c
      write (iunit,6210) 'NPOINT',ndone
      write (iunit,6210) 'COLOUR',4
c
      write (iunit,6200) 'XVALUE','*'
      write (iunit,6220) (nint(buffb1(i)),i=1,ndone)
c
      write (iunit,6215) 'REMARK',
     +    'SIGMA(chi1) ave, sdv, min, max ',
     +    sphiav,sphisd,sphimi,sphima
      write (iunit,6215) 'REMARK',
     +    'RANGE(chi1) ave, sdv, min, max ',
     +    rphiav,rphisd,rphimi,rphima
      write (iunit,6215) 'REMARK',
     +    'CV(chi1) ave, sdv, min, max ',
     +    cphiav,cphisd,cphimi,cphima
c
      if (mode .eq. 1) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_SIGCHI1 R ',
     +    ndone,' (9f8.2)'
      else if (mode .eq. 3) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_CVCHI1 R ',
     +    ndone,' (9f8.3)'
      end if
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
      write (iunit,6200) 'YVALUE','*'
c
      if (mode .eq. 1) then
        write (iunit,6230) (buffb1(maxres+i),i=1,ndone)
      else if (mode .eq. 3) then
        write (iunit,6240) (buffb2(i),i=1,ndone)
      end if
c
      write (iunit,6200) 'MORE  '
      write (iunit,6210) 'COLOUR',1
c
      write (iunit,6215) 'REMARK','SIGMA(chi2) ave, sdv, min, max ',
     +  spsiav,spsisd,spsimi,spsima
      write (iunit,6215) 'REMARK','RANGE(chi2) ave, sdv, min, max ',
     +  rpsiav,rpsisd,rpsimi,rpsima
      write (iunit,6215) 'REMARK','CV(chi2) ave, sdv, min, max ',
     +  cpsiav,cpsisd,cpsimi,cpsima
c
      if (mode .eq. 1) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_SIGCHI2 R ',
     +    ndone,' (9f8.2)'
      else if (mode .eq. 3) then
        write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_CVCHI2 R ',
     +    ndone,' (9f8.3)'
      end if
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      write (iunit,6200) 'YVALUE','*'
c
      if (mode .eq. 1) then
        write (iunit,6230) (buffb1(3*maxres+i),i=1,ndone)
      else if (mode .eq. 3) then
        write (iunit,6240) (buffb2(maxres+i),i=1,ndone)
      end if
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      ierr = 0
      return
c
 6000 format (1x,a3,i6,' | ',2(4f7.1,:,' | '),2i5,' | ',2f7.3)
 6010 format (1x,a9,3x,2(4a7,3x),2a5,3x,2a7)
c
      end
c
c
c
      subroutine nomenu (imol)
c
      include 'lsqman.incl'
c
      integer maxnat
      parameter (maxnat=5)
c
      character*1 quote,star
      parameter (quote = '''', star='*')
c
      integer imol,i,itype,nn,nr
c
      character typ3lc(maxnat)*3,dna3lc(maxnat)*3
c
      data typ3lc /'  A','  G','  C','  T','  U'/
      data dna3lc /' DA',' DG',' DC',' DT','  U'/
c
code ...
c
      nn = 0
      nr = 0
      do i=1,natoms(imol)
        call nuctyp (resnam(i,imol),itype)
        if (itype .gt. 0) then
          if (resnam(i,imol) .ne. typ3lc(itype) .and.
     +        resnam(i,imol) .ne. dna3lc(itype)) then
            resnam(i,imol) = typ3lc(itype)
            nr = nr + 1
          end if
          if (atmnam(i,imol)(4:4) .eq. quote) then
            atmnam(i,imol)(4:4) = star
            nn = nn + 1
          end if
        end if
      end do
c
      call jvalut (' Atoms with changed residue type :',1,nr)
      call jvalut (' Atoms with changed atom type    :',1,nn)
c
      return
      end
c
c
c
      subroutine nomen (imol)
c
      include 'lsqman.incl'
c
      integer nchk(5),nerr(5)
      integer imol,i,j,k,oldres,nres,ires
c
      logical lerror
c
code ...
c
      oldres = -99999
      nres = 0
c
      call jvalut (' Nr of atoms    :',1,natoms(imol))
      if (natoms(imol) .lt. 5) then
        call errcon ('Not enough atoms')
        return
      end if
c
      do i=1,natoms(imol)
        ires = iresid (i,imol)
        if (oldres .ne. ires) then
          nres = nres + 1
          oldres = ires
          buffb1 (nres) = float(i)
        end if
      end do
      buffb1 (nres+1) = float(natoms(imol) + 1)
c
      call jvalut (' Nr of residues :',1,nres)
      if (nres .lt. 1) return
c
      do i=1,5
        nchk (i) = 0
        nerr (i) = 0
      end do
c
      write (*,*)
      do i=1,nres
        j = nint(buffb1(i))
        k = nint(buffb1(i+1)-1)
        if (resnam(j,imol) .eq. 'PHE') then
          call chkres (resnam(j,imol),j,k,atmxyz(1,1,imol),
     +      atmnam(1,imol),achain(j,imol),
     +      iresid(j,imol),natoms(imol),lerror)
          nchk(1) = nchk(1) + 1
          if (lerror) nerr(1) = nerr(1) + 1
        else if (resnam(j,imol) .eq. 'TYR') then
          call chkres (resnam(j,imol),j,k,atmxyz(1,1,imol),
     +      atmnam(1,imol),achain(j,imol),
     +      iresid(j,imol),natoms(imol),lerror)
          nchk(2) = nchk(2) + 1
          if (lerror) nerr(2) = nerr(2) + 1
        else if (resnam(j,imol) .eq. 'ASP') then
          call chkres (resnam(j,imol),j,k,atmxyz(1,1,imol),
     +      atmnam(1,imol),achain(j,imol),
     +      iresid(j,imol),natoms(imol),lerror)
          nchk(3) = nchk(3) + 1
          if (lerror) nerr(3) = nerr(3) + 1
        else if (resnam(j,imol) .eq. 'GLU') then
          call chkres (resnam(j,imol),j,k,atmxyz(1,1,imol),
     +      atmnam(1,imol),achain(j,imol),
     +      iresid(j,imol),natoms(imol),lerror)
          nchk(4) = nchk(4) + 1
          if (lerror) nerr(4) = nerr(4) + 1
        else if (resnam(j,imol) .eq. 'ARG') then
          call chkres (resnam(j,imol),j,k,atmxyz(1,1,imol),
     +      atmnam(1,imol),achain(j,imol),
     +      iresid(j,imol),natoms(imol),lerror)
          nchk(5) = nchk(5) + 1
          if (lerror) nerr(5) = nerr(5) + 1
        end if
      end do
c
 6100 format (' # of ',a3,' checked : ',i5,' # errors : ',i5)
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
        call prompt (
     +    ' WARNING - any attached hydrogens NOT renamed')
      else
        call prompt (' No problem, mon !')
      end if
c
      return
      end
c
c
c
      subroutine chkres (resnam,i1,i2,allxyz,atmnam,chain,
     +                   id,natoms,lerror)
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
      logical lerror
c
      character atmnam(natoms)*4,resnam*3,chain*1
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
 6000 format (' Error in ',a3,1x,a1,1x,i6,' ...')
c
      if (resnam .eq. 'PHE' .or. resnam .eq. 'TYR') then
        if (min(ica,icb,icg,icd1,icd2) .le. 0) then
          call errcon ('Missing atom(s) !')
          return
        end if
        x1 = tangle(ica,icb,icg,icd1,allxyz)
        x2 = tangle(ica,icb,icg,icd2,allxyz)
        if (abs(x1) .gt. abs(x2)) then
            write (*,6000) resnam,chain,id
            atmnam (icd1) = ' CD2'
            atmnam (icd2) = ' CD1'
            atmnam (ice1) = ' CE2'
            atmnam (ice2) = ' CE1'
cc            call prompt (' Swapped CD1/2 and CE1/2')
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
            write (*,6000) resnam,chain,id
            atmnam (iod1) = ' OD2'
            atmnam (iod2) = ' OD1'
cc            call prompt (' Swapped OD1/2')
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
            write (*,6000) resnam,chain,id
            atmnam (ioe1) = ' OE2'
            atmnam (ioe2) = ' OE1'
cc            call prompt (' Swapped OE1/2')
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
            write (*,6000) resnam,chain,id
            atmnam (inh1) = ' NH2'
            atmnam (inh2) = ' NH1'
cc            call prompt (' Swapped NH1/2')
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
      subroutine mmfold (xave,yave,nphi,phi,psi,type)
c
c ... draw lines from centroid to individual points, but
c     use periodicity to avoid very long lines across the plot
c
c ... xave, yave = centroid angles
c     nphi = nr of individual angles
c     phi, psi = array of angles
c     type = 1 for [-180,+180> angles, 2 for [0,360> angles
c
      implicit none
c
      integer nphi,type
c
      real phi(nphi),psi(nphi),xave,yave,xlo,ylo,xhi,yhi
c
      integer i
c
code ...
c
      if (nphi .le. 1) return
c
      if (type .eq. 1) then
        xlo = -180.0
        xhi =  180.0
      else
        xlo =   0.0
        xhi = 360.0
      end if
c
      ylo = xlo
      yhi = xhi
c
      if (nphi .eq. 2) then
        call plout (phi(1),psi(1),phi(2),psi(2),xlo,xhi,ylo,yhi,type)
      else
        do i=1,nphi
          call plout (xave,yave,phi(i),psi(i),xlo,xhi,ylo,yhi,type)
        end do
      end if
c
      return
      end
c
c
c
      subroutine plout (x1,y1,x2,y2,xlo,xhi,ylo,yhi,type)
c
c ... draw connecting lines
c
      implicit none
c
      real x1,y1,x2,y2,xlo,xhi,ylo,yhi,p1,q1,p2,q2
c
      integer type
c
      logical ldone
c
code ...
c
      p1 = x1
      q1 = y1
      p2 = x2
      q2 = y2
c
ccc      write (*,'(1x,a,4(1x,f8.1))') 'PLOUT ',x1,y1,x2,y2
c
c ... get equivalent of (x1,y1) inside proper range -> (p1,q1)
c
      if (type .eq. 1) then
        call fixang (p1)
        call fixang (q1)
      else
        call fix360 (p1)
        call fix360 (q1)
      end if
c
c ... get equivalent of (x2,y2) nearest (p1,q1) -> (p2,q2)
c
      call closan (p2,p1)
      call closan (q2,q1)
      call shortest (p1,q1,p2,q2)
c
c ... first line (p1,q1) in direction of (p2,q2)
c
      call hidden (p1,q1,p2,q2,xlo,xhi,ylo,yhi,ldone)
c
      if (ldone) return
c
      p1 = x1
      q1 = y1
      p2 = x2
      q2 = y2
c
      if (type .eq. 1) then
        call fixang (p2)
        call fixang (q2)
      else
        call fix360 (p2)
        call fix360 (q2)
      end if
c
      call closan (p1,p2)
      call closan (q1,q2)
      call shortest (p2,q2,p1,q1)
c
c ... second line (p2,q2) in direction of (p1,q1)
c
      call hidden (p2,q2,p1,q1,xlo,xhi,ylo,yhi,ldone)
c
      return
      end
c
c
c
      subroutine hidden (x1,y1,x2,y2,xlo,xhi,ylo,yhi,ldone)
c
c ... draw a line but don't exceed the frame
c
      implicit none
c
      real small
      parameter (small=1.0E-6)
c
      real x1,y1,x2,y2,xlo,xhi,ylo,yhi
      real dx,dy,yx,xy,dmin,d,xp,yp
c
      logical ldone
c
code ...
c
ccc      write (*,'(1x,a,4(1x,f8.1))') 'HIDDEN',x1,y1,x2,y2
c
c ... is line direct ?
c
      if (x1 .ge. xlo .and. x1 .le. xhi) then
        if (x2 .ge. xlo .and. x2 .le. xhi) then
          if (y1 .ge. ylo .and. y1 .le. yhi) then
            if (y2 .ge. ylo .and. y2 .le. yhi) then
              ldone = .true.
              call xps_move (x1,y1)
              call xps_draw (x2,y2)
              return
            end if
          end if
        end if
      end if
c
      xp = x2
      yp = y2
      dmin = (x1-x2)**2 + (y1-y2)**2
c
c ... if the line is 'direct' return TRUE for LDONE
c
      ldone = .true.
c
c ... get direction vector of line       
c
      dx = (x2 - x1)
      dy = (y2 - y1)
c
      if (abs(dx) .lt. small) goto 100
c
c ... intersection with the line X=XLO
c
      yx = y1 + dy * (xlo-x1) / dx
      if (yx .ge. ylo .and. yx .le. yhi) then
        d = (x1-xlo)**2 + (y1-yx)**2
        if (d .le. dmin) then
          xp = xlo
          yp = yx
          dmin = d
          ldone = .false.
        end if
      end if
c
c ... intersection with the line X=XHI
c
      yx = y1 + dy * (xhi-x1) / dx
      if (yx .ge. ylo .and. yx .le. yhi) then
        d = (x1-xhi)**2 + (y1-yx)**2
        if (d .le. dmin) then
          xp = xhi
          yp = yx
          dmin = d
          ldone = .false.
        end if
      end if
c
  100 continue
      if (abs(dy) .lt. small) goto 200
c
c ... intersection with the line Y=YLO
c
      xy = x1 + dx * (ylo-y1) / dy
      if (xy .ge. xlo .and. xy .le. xhi) then
        d = (x1-xy)**2 + (y1-ylo)**2
        if (d .le. dmin) then
          xp = xy
          yp = ylo
          dmin = d
          ldone = .false.
        end if
      end if
c
c ... intersection with the line Y=YHI
c
      xy = x1 + dx * (yhi-y1) / dy
      if (xy .ge. xlo .and. xy .le. xhi) then
        d = (x1-xy)**2 + (y1-yhi)**2
        if (d .le. dmin) then
          xp = xy
          yp = yhi
          dmin = d
          ldone = .false.
        end if
      end if
c
c ... draw the line piece
c
  200 continue
      call xps_move (x1,y1)
      call xps_draw (xp,yp)
c
ccc      write (*,'(1x,a,4(1x,f8.1))') '   OUT',x1,y1,xp,yp
c
      return
      end
c
c
c
      subroutine shortest (p1,q1,p2,q2)
c
c ... find angle pair (p2,q2) which has the shortest "distance"
c     to (p1,q1)
c
      implicit none
c
      real p1,q1,p2,q2,pmin,qmin,dmin,pi,pj,d
c
      integer i,j
c
code ...
c
      dmin = (p1-p2)**2 + (q1-q2)**2
      pmin = p2
      qmin = q2
c
ccc      write (*,'(1x,a10,5f8.1)') 'SHORTEST ',p1,q1,p2,q2,dmin
c
      do i=-1,1,1
        pi = p2 + float(i)*360.0
        do j=-1,1,1
          pj = q2 + float(j)*360.0
          d = (p1-pi)**2 + (q1-pj)**2
          if (d .lt. dmin) then
            dmin = d
            pmin = pi
            qmin = pj
          end if
        end do
      end do
c
      p2 = pmin
      q2 = qmin
c
ccc      write (*,'(1x,a10,5f8.1)') 'AFTER ',p1,q1,p2,q2,dmin
c
      return
      end
c
c
c
      subroutine subavb (imol)
c
      include 'lsqman.incl'
c
      real ave
c
      integer imol,nh,ich,j
c
      logical lhydro
c
code ...
c
      call textut (' Subtract average chain B for :',name(imol))
c
 6000 format (' Chain ',a1,' # non-H atoms = ',i6,' <B> = ',
     +  f6.2,' A**2')
c
      do ich=1,nchain(imol)
        ave = 0.0
        nh = 0
        do j=chnptr(1,ich,imol),chnptr(2,ich,imol)
          if (.not. lhydro(atmnam(j,imol))) then
            ave = ave + batom(j,imol)
            nh = nh + 1
          end if
        end do
        ave = ave / float (nh)
        write (*,6000) chname(ich,imol),nh,ave
        do j=chnptr(1,ich,imol),chnptr(2,ich,imol)
          if (.not. lhydro(atmnam(j,imol))) then
            batom(j,imol) = batom(j,imol) - ave
          end if
        end do
      end do
c
      return
      end
c
c
c
      subroutine mcrho (nat,x1,x2,rms,rho,usm,nor,dsq)
c
c ... calculate Maiorov-Crippen RHO and some other statistics
c
      implicit none
c
      real x1(3,*),x2(3,*)
      real c1(3),c2(3)
      real rms,rho,r1,r2,usm,cn,dsq,nor,qqq
c
      integer nat,i,j
c
code ...
c
      rho = 999.9999
      usm = 999.9999
      nor = 999.9999
c
      if (nat .le. 0) return
c
c ... normalised RMSD (100); see: Carugo & Pongor, Prot Sci 10,
c     1470-1473 (2001)
c
      if (nat .gt. 14) then
        nor = rms / (1.0 + alog(sqrt(0.01*float(nat))))
      end if
c
c ... calc centres of gravity
c
      do i=1,3
        c1(i)=0.0
        c2(i)=0.0
      end do
c
      do i=1,nat
        do j=1,3
          c1(j)=c1(j)+x1(j,i)
          c2(j)=c2(j)+x2(j,i)
        end do
      end do
c
      do i=1,3
        c1(i)=c1(i)/float(nat)
        c2(i)=c2(i)/float(nat)
      end do
c
c ... calc radii of gyration
c
      r1 = 0.0
      r2 = 0.0
      do i=1,nat
        do j=1,3
          r1 = r1 + (x1(j,i)-c1(j))**2
          r2 = r2 + (x2(j,i)-c2(j))**2
        end do
      end do
      r1 = sqrt ( r1 / float(nat) )
      r2 = sqrt ( r2 / float(nat) )
c
c ... calc rho
c
      qqq = 2.0*r1*r1 + 2.0*r2*r2 - rms*rms
      if (qqq .ge. 0.0) then
        rho = 2.0 * rms / sqrt (2.0*r1*r1 + 2.0*r2*r2 - rms*rms)
      else
        call errcon ('MCRHO - Cannot calculate Crippen rho')
      end if
c
c ... calculate Universal Similarity Measure:
c     <D>**2 = r1**2 + r2**2 - 2*c*r1*r2
c     c = 0.42 - 0.05*(N-1)*EXP(-(N-1)/4.7) + 0.63*EXP(-(N-1)/37)
c     RRMSD = USM = rmsd / <D>
c
      cn = 0.42 - 0.05 * float(nat-1) * exp (- float(nat-1)/4.7) +
     +     0.63 * exp (- float(nat-1)/37.0)
c
      dsq = r1**2 + r2**2 - 2*cn*r1*r2
      if (dsq .gt. 0.0) then
        dsq = sqrt (dsq)
        usm = rms / dsq
      else
        call errcon ('MCRHO - Cannot calculate RRMSD')
      end if
c
c      print *,'r1,r2,rms,rho,cn,dsq,usm = ',
c     +         r1,r2,rms,rho,cn,dsq,usm
c
      return
      end
c
c
c
      subroutine hisdis (imol,jmol,bin)
c
      include 'lsqman.incl'
c
      integer maxbin
      parameter (maxbin = 25)
c
      real rtx(12),x,cw2,dmin,rms
      real xave,xsdv,xmin,xmax,xtot
      real q1,q2,q3,q4,bin
c
      integer imol,jmol,i,j,nmat,nmin,nmax
c
code ...
c
      call textut (' Distance histogram of mol 1 :',name(imol))
      call textut ('                   and mol 2 :',name(jmol))
c
c ... apply current operator to mol 2
c
      do i=1,12
        rtx (i) = rtlsq (i,imol,jmol)
      end do
      call fvalut (' Applying current operator to mol 2 :',12,rtx)
      call vecrtv (atmxyz(1,1,jmol),buffi,natoms(jmol),rtx(1),rtx(10))
c
c ... match atoms of mol 2
c
      cw2 = cuthis * cuthis
      nmat = 0
      rms = 0.0
      do i=1,natoms(jmol)
        dmin = cw2 + 0.001
        nmin = -1
        do j=1,natoms(imol)
          x = (buffi(1,i)-atmxyz(1,j,imol))**2 +
     +        (buffi(2,i)-atmxyz(2,j,imol))**2 +
     +        (buffi(3,i)-atmxyz(3,j,imol))**2
          if (x .lt. dmin) then
            dmin = x
            nmin = j
          end if
        end do
        if (nmin .gt. 0) then
          nmat = nmat + 1
          rms = rms + dmin
          x = sqrt (dmin)
          buffb3 (nmat) = x
        else
          write (*,6000) atmnam(i,jmol),resnam(i,jmol),
     +      achain(i,jmol),iresid(i,jmol)
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of matched atoms   :',1,nmat)
      call jvalut (' Nr of atoms in mol 2  :',1,natoms(jmol))
      x = 100.0 * float(nmat) / float(natoms(jmol))
      call fvalut (' % matched atoms mol 2 :',1,x)
      if (nmat .lt. 3) then
        call errcon ('Too few matched atoms')
        return
      end if
c
      rms = sqrt (rms / float(nmat))
      call fvalut (' RMS distance (A) :',1,rms)
c
      write (*,*)
      call xstats (buffb3,nmat,xave,xsdv,xmin,xmax,xtot)
      write (*,6010) 'Matching distances',xave,xsdv,xmin,xmax
c
 6000 format (' NOT MATCHED :',a4,1x,a3,'-',a1,i4)
c
 6010 format (1x,A,' :'/
     +  ' Average :',f8.2,'   St. dev. :',f8.2/
     +  ' Minimum :',f8.2,'   Maximum  :',f8.2)
c
 6020 format (1x,i8,' in [',f6.2,'-',f6.2,'> = ',f6.2,
     +  '% (Cumul ',f6.2,' %)')
c
      nmax = 1 + int (cuthis/bin)
      if (nmax .gt. maxatm) nmax = maxatm
      do i=1,nmax
        buffb1(i) = 0.0
      end do
      do i=1,nmat
        j = 1 + int(buffb3(i)/bin)
        if (j .gt. nmax) j = nmax
        buffb1(j) = buffb1(j) + 1.0
      end do
c
      write (*,*)
      q4 = 0.0
      do i=1,nmax
        if (buffb1(i) .gt. 0.5) then
          q1 = float(i-1)*bin
          q2 = q1 + bin
          q3 = 100.0 * buffb1(i) / float(nmat)
          q4 = q4 + q3
          write (*,6020) nint(buffb1(i)),q1,q2,q3,q4
        end if
      end do
c
      return
      end
c
c
c
      subroutine fixatm (imol,irange,jmol,jrange,
     +                   cacut,mode,how,what,delta,ierr)
c
      include 'lsqman.incl'
c
      integer maxopt
      parameter (maxopt = 25)
c
      real cacut,d1,d2,distce,t1,t2,t3,tangle,delta
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j
      integer idum,jdum,kdum,ldum,nuse,iptr,jptr,ires,jres,kres
      integer nbskip,leng1,i1,i2,i3,i4,j1,j2,j3,j4,i5,j5,i6,j6
c
      logical lstrict,lseq,lrmsd
c
      character irange*(*),jrange*(*),ichn*1,jchn*1,how*1,what*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20,mode*1
c
code ...
c
      ierr = 0
c
      write (*,6000) name(imol)(1:leng1(name(imol))),
     +    irange(1:leng1(irange)),
     +    name(jmol)(1:leng1(name(jmol))),
     +    jrange(1:leng1(jrange))
c
 6000 format (' Reference atoms ',a,1x,a/
     +        ' Fix atoms for   ',a,1x,a)
c
      if (mode(1:1) .eq. 'S') then
        lstrict = .true.
        call prompt (' Only fix Asp/Glu/Arg/Phe/Tyr')
      else
        lstrict = .false.
        call prompt (' Fix Asp/Glu/Arg/Phe/Tyr AND Asn/Gln/His')
      end if
c
      if (how(1:1) .eq. 'S') then
        lseq = .true.
        call prompt (' Use sequential residues (1:1 correspondence)')
      else
        lseq = .false.
        call prompt (' Find nearest residue (CA-CA)')
        call fvalut (' Cut-off distance :',1,cacut)
      end if
c
      if (what(1:1) .eq. 'R') then
        lrmsd = .true.
        call prompt (' Minimise RMSD')
      else
        lrmsd = .false.
        call prompt (' Minimise torsion-angle differences')
      end if
c
      delta = abs(delta)
      call fvalut (' Minimum improvement :',1,delta)
c
c ... go do it
c
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range1')
        goto 9900
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range2')
        goto 9900
      end if
c
      if (nopt1 .ne. nopt2) then
        call errcon ('Different nr of zones')
        goto 9900
      end if
c
c ... apply current operator to mol 2
c
      if ( lrmsd .or. (.not. lseq) ) then
        call prompt (' Applying current operator to Mol 2 ...')
        call anancs (1,rtlsq(1,imol,jmol),.true.,ierr)
      end if
      call vecrtv (atmxyz(1,1,jmol),buffk(1,1),natoms(jmol),
     +  rtlsq(1,imol,jmol),rtlsq(10,imol,jmol))
c
      nuse = 0
      nbskip = 0
      iptr = 0
      jptr = 0
c
 6100 format (' Fix sidechain of ',a3,'-',a1,'-',i4,' (',f8.2,
     +  ' versus ',f8.2,')')
c
      do i=1,nopt1
c
        write (*,*)
        call jvalut (' Zone :',1,i)
c
        ichn = optpa1(i)(1:1)
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range1 : '//optpa1(i))
ccc          goto 9900
        end if
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
          call errcon ('Not a zone in range1 : '//optpa1(i))
          goto 9900
        end if
ccc      print *,optpa1(i)(2:j-1)
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
ccc      print *,optpa1(i)(j+1:)
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon ('Invalid zone in range1 : '//optpa1(i))
          goto 9900
        end if
c
        jchn = optpa2(i)(1:1)
        if (index('1234567890',jchn) .gt. 0) then
          jchn = ' '
ccc        else if (jchn .lt. 'A' .or. jchn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range2 : '//optpa2(i))
ccc          goto 9900
        end if
        j = index (optpa2(i),'-')
        if (j .le. 0) j = index (optpa2(i),':')
        if (j .gt. 0) then
          call errcon (
     +      'Zone ignored in range2 : '//optpa2(i))
          optpa2(i) = optpa2(i)(1:j-1)
        end if
ccc      print *,optpa2(i)(2:)
        if (jchn .eq. ' ') then
          call str2i (optpa2(i)(1:),kdum,ierr)
        else
          call str2i (optpa2(i)(2:),kdum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        ldum = kdum + (jdum - idum)
c
        do ires=idum,jdum
c
c ... check if residue exists in mol 1
c
          call getptr (imol,ichn,ires,' CA ',
     +                 -1,jchn,jres,iptr,jptr,ierr)
          if (ierr .ne. 0) goto 6996
c
c ... check if interesting type
c
          if (resnam(iptr,imol) .eq. 'ASP' .or.
     +        resnam(iptr,imol) .eq. 'GLU' .or.
     +        resnam(iptr,imol) .eq. 'ARG' .or.
     +        resnam(iptr,imol) .eq. 'PHE' .or.
     +        resnam(iptr,imol) .eq. 'TYR') goto 6699
c
          if (.not. lstrict) then
            if (resnam(iptr,imol) .eq. 'ASN' .or.
     +          resnam(iptr,imol) .eq. 'GLN' .or.
     +          resnam(iptr,imol) .eq. 'HIS') goto 6699
          end if
c
          goto 6996
c
c ... okay
c
 6699     continue
c
c ... if sequential, look for residue in mol 2
c
          if (lseq) then
            jres = kdum + (ires - idum)
            call getptr (imol,ichn,ires,' CA ',
     +                   jmol,jchn,jres,iptr,jptr,ierr)
            if (ierr .ne. 0) goto 6996
          else
c
c ... otherwise, find nearest residue in mol 2
c
            d2 = cacut + 0.001
            jptr = -1
            jres = -1
            do kres = kdum,ldum
              call getptr (-1,ichn,ires,' CA ',
     +                     jmol,jchn,kres,i1,j1,ierr)
              if (ierr. eq. 0) then
                d1 = distce (atmxyz(1,iptr,imol),buffk(1,j1))
                if (d1 .lt. d2) then
                  d2 = d1
                  jptr = j1
                  jres = kres
                else
c
c ... if we already have an acceptable hit, and we now
c     get very far away, skip the rest of the molecule
c     (dangerous !)
c
                  if (jptr .gt. 0 .and. d1 .gt. (5.0*cacut))
     +              goto 6669
                end if
              end if
 6999         continue
            end do
c
c ... found one ?
c
 6669       continue
            if (jptr .le. 0) goto 6996
c
          end if
c
c ... check if interesting type
c
          if (resnam(jptr,jmol) .eq. 'ASP' .or.
     +        resnam(jptr,jmol) .eq. 'GLU' .or.
     +        resnam(jptr,jmol) .eq. 'ARG' .or.
     +        resnam(jptr,jmol) .eq. 'PHE' .or.
     +        resnam(jptr,jmol) .eq. 'TYR') goto 6698
c
          if (.not. lstrict) then
            if (resnam(jptr,jmol) .eq. 'ASN' .or.
     +          resnam(jptr,jmol) .eq. 'GLN' .or.
     +          resnam(jptr,jmol) .eq. 'HIS') goto 6698
          end if
c
          goto 6996
c
 6698     continue
c
ccc          print *
ccc          print *,resnam(iptr,imol),ichn,ires
ccc          print *,resnam(jptr,jmol),jchn,jres
c
c ... now we have IPTR = CA in residue in mol 1
c             and JPTR = CA in residue in mol 2
c
c ... ASP ... ASP ... ASP ... ASP ... ASP ... ASP ... ASP ... ASP ... ASP 
c
          if (resnam(iptr,imol) .eq. 'ASP' .and.
     +        resnam(jptr,jmol) .eq. 'ASP') then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' OD1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' OD2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1))
              d1 = sqrt(d1/2.0)
              d2 = sqrt(d2/2.0)
            else
              call getptr (imol,ichn,ires,' CB ',
     +                     jmol,jchn,jres,i3,j3,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CG ',
     +                     jmol,jchn,jres,i4,j4,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (iptr,i3,i4,i1,atmxyz(1,1,imol))
              t2 = tangle (jptr,j3,j4,j1,atmxyz(1,1,jmol))
              t3 = tangle (jptr,j3,j4,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' OD2'
              atmnam (j2,jmol) = ' OD1'
              write (*,6100) 'ASP',jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
c ... GLU ... GLU ... GLU ... GLU ... GLU ... GLU ... GLU ... GLU ... GLU
c
          else if (resnam(iptr,imol) .eq. 'GLU' .and.
     +             resnam(jptr,jmol) .eq. 'GLU') then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' OE1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' OE2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1))
              d1 = sqrt(d1/2.0)
              d2 = sqrt(d2/2.0)
            else
              call getptr (imol,ichn,ires,' CB ',
     +                     jmol,jchn,jres,i3,j3,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CG ',
     +                     jmol,jchn,jres,i4,j4,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CD ',
     +                     jmol,jchn,jres,i5,j5,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (i3,i4,i5,i1,atmxyz(1,1,imol))
              t2 = tangle (j3,j4,j5,j1,atmxyz(1,1,jmol))
              t3 = tangle (j3,j4,j5,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' OE2'
              atmnam (j2,jmol) = ' OE1'
              write (*,6100) 'GLU',jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
c ... ARG ... ARG ... ARG ... ARG ... ARG ... ARG ... ARG ... ARG ... ARG
c
          else if (resnam(iptr,imol) .eq. 'ARG' .and.
     +             resnam(jptr,jmol) .eq. 'ARG') then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' NH1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' NH2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1))
              d1 = sqrt(d1/2.0)
              d2 = sqrt(d2/2.0)
            else
              call getptr (imol,ichn,ires,' CD ',
     +                     jmol,jchn,jres,i3,j3,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' NE ',
     +                     jmol,jchn,jres,i4,j4,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CZ ',
     +                     jmol,jchn,jres,i5,j5,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (i3,i4,i5,i1,atmxyz(1,1,imol))
              t2 = tangle (j3,j4,j5,j1,atmxyz(1,1,jmol))
              t3 = tangle (j3,j4,j5,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' NH2'
              atmnam (j2,jmol) = ' NH1'
              write (*,6100) 'ARG',jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
c ... PHE ... PHE ... PHE ... PHE ... PHE ... PHE ... PHE ... PHE ... PHE
c
c ... TYR ... TYR ... TYR ... TYR ... TYR ... TYR ... TYR ... TYR ... TYR
c
          else if ( (resnam(iptr,imol) .eq. 'PHE' .or.
     +               resnam(iptr,imol) .eq. 'TYR') .and.
     +              (resnam(jptr,jmol) .eq. 'PHE' .or.
     +               resnam(jptr,jmol) .eq. 'TYR') ) then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' CD1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' CD2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' CE1',
     +                   jmol,jchn,jres,i3,j3,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' CE2',
     +                   jmol,jchn,jres,i4,j4,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i3,imol),buffk(1,j3)) +
     +             distce (atmxyz(1,i4,imol),buffk(1,j4))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i3,imol),buffk(1,j4)) +
     +             distce (atmxyz(1,i4,imol),buffk(1,j3))
              d1 = sqrt(d1/4.0)
              d2 = sqrt(d2/4.0)
            else
              call getptr (imol,ichn,ires,' CB ',
     +                     jmol,jchn,jres,i5,j5,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CG ',
     +                     jmol,jchn,jres,i6,j6,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (iptr,i5,i6,i1,atmxyz(1,1,imol))
              t2 = tangle (jptr,j5,j6,j1,atmxyz(1,1,jmol))
              t3 = tangle (jptr,j5,j6,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' CD2'
              atmnam (j2,jmol) = ' CD1'
              atmnam (j3,jmol) = ' CE2'
              atmnam (j4,jmol) = ' CE1'
              write (*,6100) resnam(iptr,imol),jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
          end if
c
          if (lstrict) goto 6996
c
c ... check ASN, GLN, HIS
c
c ... ASN ... ASN ... ASN ... ASN ... ASN ... ASN ... ASN ... ASN ... ASN
c
          if (resnam(iptr,imol) .eq. 'ASN' .and.
     +        resnam(jptr,jmol) .eq. 'ASN') then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' OD1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' ND2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1))
              d1 = sqrt(d1/2.0)
              d2 = sqrt(d2/2.0)
            else
              call getptr (imol,ichn,ires,' CB ',
     +                     jmol,jchn,jres,i3,j3,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CG ',
     +                     jmol,jchn,jres,i4,j4,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (iptr,i3,i4,i1,atmxyz(1,1,imol))
              t2 = tangle (jptr,j3,j4,j1,atmxyz(1,1,jmol))
              t3 = tangle (jptr,j3,j4,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' ND2'
              atmnam (j2,jmol) = ' OD1'
              write (*,6100) 'ASN',jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
c ... GLN ... GLN ... GLN ... GLN ... GLN ... GLN ... GLN ... GLN ... GLN
c
          else if (resnam(iptr,imol) .eq. 'GLN' .and.
     +             resnam(jptr,jmol) .eq. 'GLN') then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' OE1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' NE2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1))
              d1 = sqrt(d1/2.0)
              d2 = sqrt(d2/2.0)
            else
              call getptr (imol,ichn,ires,' CB ',
     +                     jmol,jchn,jres,i3,j3,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CG ',
     +                     jmol,jchn,jres,i4,j4,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CD ',
     +                     jmol,jchn,jres,i5,j5,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (i3,i4,i5,i1,atmxyz(1,1,imol))
              t2 = tangle (j3,j4,j5,j1,atmxyz(1,1,jmol))
              t3 = tangle (j3,j4,j5,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' NE2'
              atmnam (j2,jmol) = ' OE1'
              write (*,6100) 'GLN',jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
c ... HIS ... HIS ... HIS ... HIS ... HIS ... HIS ... HIS ... HIS ... HIS
c
          else if ( resnam(iptr,imol) .eq. 'HIS' .and.
     +              resnam(jptr,jmol) .eq. 'HIS' ) then
            nuse = nuse + 1
            call getptr (imol,ichn,ires,' ND1',
     +                   jmol,jchn,jres,i1,j1,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' CD2',
     +                   jmol,jchn,jres,i2,j2,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' CE1',
     +                   jmol,jchn,jres,i3,j3,ierr)
            if (ierr .ne. 0) goto 6996
            call getptr (imol,ichn,ires,' NE2',
     +                   jmol,jchn,jres,i4,j4,ierr)
            if (ierr .ne. 0) goto 6996
c
            if (lrmsd) then
              d1 = distce (atmxyz(1,i1,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i3,imol),buffk(1,j3)) +
     +             distce (atmxyz(1,i4,imol),buffk(1,j4))
              d2 = distce (atmxyz(1,i1,imol),buffk(1,j2)) +
     +             distce (atmxyz(1,i2,imol),buffk(1,j1)) +
     +             distce (atmxyz(1,i3,imol),buffk(1,j4)) +
     +             distce (atmxyz(1,i4,imol),buffk(1,j3))
              d1 = sqrt(d1/4.0)
              d2 = sqrt(d2/4.0)
            else
              call getptr (imol,ichn,ires,' CB ',
     +                     jmol,jchn,jres,i5,j5,ierr)
              if (ierr .ne. 0) goto 6996
              call getptr (imol,ichn,ires,' CG ',
     +                     jmol,jchn,jres,i6,j6,ierr)
              if (ierr .ne. 0) goto 6996
              t1 = tangle (iptr,i5,i6,i1,atmxyz(1,1,imol))
              t2 = tangle (jptr,j5,j6,j1,atmxyz(1,1,jmol))
              t3 = tangle (jptr,j5,j6,j2,atmxyz(1,1,jmol))
              d1 = abs(t1-t2)
              call fix360 (d1)
              d2 = abs(t1-t3)
              call fix360 (d2)
            end if
c
            if ((d2+delta) .lt. d1) then
              atmnam (j1,jmol) = ' CD2'
              atmnam (j2,jmol) = ' ND1'
              atmnam (j3,jmol) = ' NE2'
              atmnam (j4,jmol) = ' CE1'
              write (*,6100) resnam(iptr,imol),jchn,jres,d1,d2
              nbskip = nbskip + 1
            end if
            goto 6996
c
          end if
c
 6996     continue
c
        end do
c
      end do
c
      write (*,*)
      call jvalut (' Residues checked :',1,nuse)
      call jvalut (' Residues fixed   :',1,nbskip)
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      ierr = -1
c
      return
      end
c
c=====================================================
c
      subroutine fixadf (ang1,ang2,dif)
c
c ... fix difference of two angles in range <-180.0,+180.0]
c
      implicit none
c
      real ang1,ang2,dif
c
code ...
c
      dif = ang1 - ang2
c
  210 if (dif .le. -180.0) then
        ang2 = ang2 - 360.0
        dif = dif + 360.0
        goto 210
      end if
c
  220 if (dif .gt. 180.0) then
        ang2 = ang2 + 360.0
        dif = dif - 360.0
        goto 220
      end if
c
      return
      end
c
c
c
      subroutine initrt (imol,jmol)
c
c ... reset operators etc. bringing JMOL on top of IMOL
c
      include 'lsqman.incl'
c
      integer imol,jmol,i
c
code ...
c
      do i=1,12
        rtlsq(i,imol,jmol) = 0.0
      end do
      rtlsq (1,imol,jmol) = 1.0
      rtlsq (5,imol,jmol) = 1.0
      rtlsq (9,imol,jmol) = 1.0
c
      call inista (imol,jmol)
c
      return
      end
c
c
c
      subroutine inista (imol,jmol)
c
c ... reset statistics etc. bringing JMOL on top of IMOL
c
      include 'lsqman.incl'
c
      integer imol,jmol
c
code ...
c
      nmatch (imol,jmol) = 0
      last (imol,jmol) = 'none'
      rmsd (imol,jmol) = 999.99999
      rmsb (imol,jmol) = 999.99999
      corb (imol,jmol) = 999.99999
      simind (imol,jmol) = 999.99999
      matchi (imol,jmol) = 999.99999
      cripp  (imol,jmol) = 999.99999
      rrmsd  (imol,jmol) = 999.99999
      normsd (imol,jmol) = 999.99999
      rmsdna (imol,jmol) = 999.99999
      sas1 (imol,jmol) = 999.99999
      sas2 (imol,jmol) = 999.99999
      sas3 (imol,jmol) = 999.99999
      sas4 (imol,jmol) = 999.99999
c
      return
      end
c
c
c
      subroutine cirvar (angles,nang,cv)
c
c ... CIRVAR - calculate circular variance
c
      implicit none
c
      real twopi,pi,degtor,rtodeg,circle
      parameter (twopi=6.2831853071796)
      parameter (pi=0.50*twopi)
      parameter (circle=360.0)
      parameter (degtor=twopi/circle)
      parameter (rtodeg=circle/twopi)
c
      integer nang
c
      real angles(nang)
      real cv,sumc,sums
c
      integer i
c
code ...
c
      cv = -1.0
      if (nang .lt. 1) return
      if (nang .eq. 1) then
        cv = 0.0
        return
      end if
c
      sumc = 0.0
      sums = 0.0
      do i=1,nang
        sumc = sumc + cos(degtor*angles(i))
        sums = sums + sin(degtor*angles(i))
      end do
c
      cv = 1.0 - (sqrt(sumc*sumc+sums*sums)/float(nang))
c
      return
      end
c
c
c
      subroutine vrmlca (imol,chains,mxdist,caname)
c
      include 'lsqman.incl'
c
      integer imol,i,k,nca,ntot
c
      real mxdist,dist
c
      character chains*1,caname*4
c
code ...
c
      nca = 0
      ntot = 0
      k = -1
c
      do i=1,natoms(imol)
c
        if (chains .ne. '*') then
          if (achain(i,imol) .ne. chains) goto 10
        end if
c
        if (atmnam(i,imol) .eq. caname) then
c
          if (nca .gt. 0) then
ccc      print *,i,k,dist(i,k,atmxyz(1,1,imol))
            if (dist(i,k,atmxyz(1,1,imol)) .gt. mxdist) then
              call xvrml_polyline (nca,buffb1)
ccc      print *,' POLYLINE ',nca
              nca = 0
              k = -1
            end if
          end if
c
          nca = nca + 1
          ntot = ntot + 1
          buffb1(3*(nca-1)+1) = atmxyz(1,i,imol)
          buffb1(3*(nca-1)+2) = atmxyz(2,i,imol)
          buffb1(3*(nca-1)+3) = atmxyz(3,i,imol)
          k = i
        end if
   10   continue
      end do
c
      if (nca .gt. 0) then
        call xvrml_polyline (nca,buffb1)
ccc      print *,' POLYLINE ',nca
      end if
c
      call jvalut (' Nr of central atoms written :',1,ntot)
c
      return
      end
c
c
c
      subroutine lsqrms (n,x,y,rms)
c
      implicit none
c
      integer n,i
c
      real x(3,n),y(3,n)
      real rms
c
code ...
c
      rms = 999.99
      if (n .lt. 1) return
c
      rms = 0.0
      do i=1,n
        rms = rms + (x(1,i)-y(1,i))**2 + (x(2,i)-y(2,i))**2 +
     +        (x(3,i)-y(3,i))**2
      end do
c
      rms = sqrt ( rms / float(n) )
c
      return
      end
c
c
c
      subroutine casp2 (imol,ichn,kmol,kchn,castrt,castop,castep,how)
c
c ... CASP2 (...) - compare model to target using 1:1 correspondence
c
      include 'lsqman.incl'
c
      real rtx(12)
      real castrt,castop,castep,cutoff,rms,dummy
      real distce
c
      integer imol,ichn,kmol,kchn,i1,i2,i,ierr
      integer iptr,kptr,k,nuse,nnow,nresj
c
      character iichn*1,kkchn*1,how*1
c
code ...
c
      if (castop .lt. castrt) then
        castep = -1.0 * abs(castep)
      end if
      if (castep .eq. 0.0) castep = (castop-castrt)/5.0
c
      write (*,6000) name(imol),chname(ichn,imol),
     +  name(kmol),chname(kchn,kmol),
     +  castrt,castop,castep
c
 6000 format (/' CASP3 Model-Target RMSD assessment'/
     +  ' TARGET mol, chain = ',a,' (',a1,')'/
     +  ' MODEL  mol, chain = ',a,' (',a1,')'/
     +  ' Cut-off start (A) : ',f8.2/
     +  ' Cut-off stop  (A) : ',f8.2/
     +  ' Cut-off step  (A) : ',f8.2/)
c
      call prompt (' Getting pointers to atoms ...')
c
      i1 = iresid ( chnptr(1,ichn,imol) , imol )
      i2 = iresid ( chnptr(2,ichn,imol) , imol )
c
      do i=chnptr(1,ichn,imol),chnptr(2,ichn,imol)
        i2 = max (i2,iresid(i,imol))
      end do
c
ccc      print *,' residues ',i1,i2
c
      iichn = chname(ichn,imol)
      kkchn = chname(kchn,kmol)
      nresj = 0
      iptr = 0
      kptr = 0
c
      do i=i1,i2
        call getptr (imol,iichn,i,' CA ',kmol,kkchn,i,
     +               iptr,kptr,ierr)
ccc        print *,i,ierr,iptr,kptr
        if (ierr .ne. 0) goto 100
        if (iptr .gt. 0 .and. kptr .gt. 0) then
          nresj = nresj + 1
          ptri (nresj) = iptr
          ptrj (nresj) = kptr
          do k=1,3
            buffi(k,nresj) = atmxyz(k,iptr,imol)
            buffj(k,nresj) = atmxyz(k,kptr,kmol)
          end do
        end if
  100   continue
      end do
c
      call jvalut (' Residues in common   :',1,nresj)
c
      if (nresj .lt. 3) then
        call errcon ('Fewer than 3 residues matched')
        goto 200
      end if
c
      if (how .eq. 'E') then
        call prompt (' Calculating explicit operator')
        call lsqgjk (buffi,buffj,nresj,rms,rtx,ierr)
        call fvalut (' RMSD all CA atoms (A) :',1,rms)
      else
        call prompt (' Copying current operator')
        do i=1,12
          rtx (i) = rtlsq (i,imol,kmol)
        end do
      end if
c
      do cutoff=castrt,castop,castep
c
        nnow = 0
c
        write (*,*)
        call fvalut (' Cut-off now (A) :',1,cutoff)
c
c ... apply current operator
c
  300   continue
        nuse = nnow
        nnow = 0
        call vecrtv (buffj,buffk,nresj,rtx(1),rtx(10))
c
        do i=1,nresj
          dummy = distce (buffi(1,i),buffk(1,i))
          if (dummy .le. cutoff) then
            nnow = nnow + 1
            do k=1,3
              buffl (k,nnow) = buffi (k,i)
              buffm (k,nnow) = buffj (k,i)
            end do
          end if
        end do
c
        if (nnow .lt. 3) then
          call errcon ('Fewer than 3 residues matched')
          goto 200
        end if
c
        if (nnow .gt. nuse) then
          call jvalut (' Nr of matched residues :',1,nnow)
          call lsqgjk (buffl,buffm,nnow,rms,rtx,ierr)
          call fvalut (' RMSD (A) :',1,rms)
          goto 300
        end if
c
      end do
c
  200 continue
c
      return
      end
c
c
c
      subroutine casp3 (imol,ichn,jmol,jchn,kmol,kchn,
     +                  judis,juphi,juchi)
c
c ... CASP3 (...) - judge model (KMOL) versus target (IMOL)
c                   by reference to parent (JMOL)
c
      include 'lsqman.incl'
c
      integer maxins
      parameter (maxins=100)
c
      real rmsglo(maxins),rmsloc(maxins)
      real judis,juphi,juchi,xdum,dis,carmsd,inrmsd,pcloca
      real phirms,mphirms,pclophi,d1,d2,q3,q4,q5,q6,q7,q8
      real cordis,q1,q2,corall,psirms,mpsirms,pclopsi
      real chirms,mchirms,pclochi,pclo2,wavglo,wavloc,pnin
      real distce
c
      integer numins(maxins)
      integer imol,ichn,jmol,jchn,kmol,kchn,nresi,nresj,nresk
      integer ntot,ica,jca,kca,nca,i,j,k,nin,ncloca,nphi,mphi
      integer nclophi,jok,kok,npsi,mpsi,nclopsi,nchi,mchi,nclochi
      integer ncaclo2,nca2,ninsert,ninres,kk
c
      logical linsert
c
      character line*128,tres*10,pres*10,mres*10
c
code ...
c
      write (*,6000) name(imol),chname(ichn,imol),
     +  name(jmol),chname(jchn,jmol),name(kmol),chname(kchn,kmol),
     +  judis,juphi,juchi
c
 6000 format (/' CASP3 Judgment'/
     +  ' TARGET mol, chain = ',a,' (',a1,')'/
     +  ' PARENT mol, chain = ',a,' (',a1,')'/
     +  ' MODEL  mol, chain = ',a,' (',a1,')'/
     +  ' CA-CA   cut-off : ',f8.2/
     +  ' Phi-Psi cut-off : ',f8.2/
     +  ' Chi1-2  cut-off : ',f8.2/)
c
      call csetup (imol,ichn,nresi,ptri,buffi)
      call csetup (jmol,jchn,nresj,ptrj,buffj)
      call csetup (kmol,kchn,nresk,ptrij,buffk)
      if (nresi .lt. 10 .or. nresj .lt. 10 .or. nresk .lt. 10) then
        call errcon (' Not enough residues !')
        return
      end if
c
      ntot = 0
      nca = 0
      nin = 0
      ncloca = 0
      ncaclo2 = 0
      nca2 = 0
      carmsd = 0.0
      inrmsd = 0.0
      phirms = 0.0
      mphirms = 0.0
      nphi = 0
      mphi = 0
      nclophi = 0
      psirms = 0.0
      mpsirms = 0.0
      npsi = 0
      mpsi = 0
      nclopsi = 0
      chirms = 0.0
      mchirms = 0.0
      nchi = 0
      mchi = 0
      nclochi = 0
c
      write (*,6110) '..PARENT..','dist (A)','..TARGET..',
     +  'dist (A)','..MODEL...'
 6100 format (1x,a10,2(' <- ',f6.2,' -> ',a10))
 6110 format (/1x,a10,2(3x,a8,3x,a10))
c
      linsert = .false.
      ninsert = 0
      ninres = 0
c
      do i=1,nresi
        ntot = ntot + 1
        ica = ptri (i)
        write (tres,'(2a,i6)') resnam(ica,imol),'-',iresid(ica,imol)
        call remspa (tres)
c
        kca = -1
        kok = -1
        do k=1,nresk
          if (iresid(ptrij(k),kmol) .eq. iresid(ica,imol)) then
            kca = ptrij (k)
            kok = k
            goto 10
          end if
        end do
        write (line,'(3a,i6,a)')
     +    'Residue ',resnam(ica,imol),' - ',
     +    iresid(ica,imol),' not in MODEL ???'
        call pretty (line)
        call errcon (line)
        goto 100
c
   10   continue
        write (mres,'(2a,i6)') resnam(kca,kmol),'-',iresid(kca,kmol)
        call remspa (mres)
c
        jca = -1
        jok = -1
        dis = judis
        do j=1,nresj
          xdum = distce(atmxyz(1,ica,imol),atmxyz(1,ptrj(j),jmol))
          if (xdum .le. dis) then
            jca = ptrj (j)
            dis = xdum
            jok = j
          end if
        end do
c
c ... found in parent ?
c
        if (jca .gt. 0) then
          write (pres,'(2a,i6)') resnam(jca,jmol),'-',
     +                           iresid(jca,jmol)
          call remspa (pres)
c
          if (linsert) then
            call anains (ninsert,ninres,buffl,buffm,
     +        rmsglo(ninsert),rmsloc(ninsert))
            numins (ninsert) = ninres
            linsert = .false.
          end if
        else
          pres = 'NO EQUIV.'
          dis = 99999.99
c
          if (.not. linsert) then
            linsert = .true.
            ninsert = ninsert + 1
            ninres = 1
            do kk=1,3
              buffl(kk,ninres) = atmxyz (kk,ica,imol)
              buffm(kk,ninres) = atmxyz (kk,kca,kmol)
            end do
          else
            ninres = ninres + 1
            do kk=1,3
              buffl(kk,ninres) = atmxyz (kk,ica,imol)
              buffm(kk,ninres) = atmxyz (kk,kca,kmol)
            end do
          end if
c
        end if
c
        nin = nin + 1
        xdum = distce (atmxyz(1,ica,imol),atmxyz(1,kca,kmol))
        inrmsd = inrmsd + xdum*xdum
        buffb3 (nin) = batom(ica,imol)
        buffb4 (nin) = xdum
        if (jca .le. 0) then
          nca= nca + 1
          carmsd = carmsd + xdum*xdum
          if (xdum .le. judis) ncloca = ncloca + 1
          buffb1 (nca) = batom(ica,imol)
          buffb2 (nca) = xdum
        else
          nca2 = nca2 + 1
          if (xdum .le. dis) ncaclo2 = ncaclo2 + 1
        end if
c
        write (*,6100) pres,dis,tres,xdum,mres
c
c ... PHI ?
c
        if (buffi(1,i) .le. -199.9) goto 20
        if (buffk(1,kok) .le. -199.9) goto 20
        call fixadf (buffi(1,i),buffk(1,kok),d2)
        d1 = 10.0*juphi
        if (jok .gt. 0) then
          if (buffj(1,jok) .gt. -199.9) then
            call fixadf (buffi(1,i),buffj(1,jok),d1)
          end if
        end if
c
        nphi = nphi + 1
        phirms = phirms + d2*d2
        if (abs(d1) .gt. juphi) then
          mphi = mphi + 1
          mphirms = mphirms + d2*d2
          if (abs(d2) .le. juphi) nclophi = nclophi + 1
        end if
c
   20   continue
c
c ... PSI ?
c
        if (buffi(2,i) .le. -199.9) goto 30
        if (buffk(2,kok) .le. -199.9) goto 30
        call fixadf (buffi(2,i),buffk(2,kok),d2)
        d1 = 10.0*juphi
        if (jok .gt. 0) then
          if (buffj(2,jok) .gt. -199.9) then
            call fixadf (buffi(2,i),buffj(2,jok),d1)
          end if
        end if
c
        npsi = npsi + 1
        psirms = psirms + d2*d2
        if (abs(d1) .gt. juphi) then
          mpsi = mpsi + 1
          mpsirms = mpsirms + d2*d2
          if (abs(d2) .le. juphi) nclopsi = nclopsi + 1
        end if
c
   30   continue
c
c ... CHI-1 ?
c
        if (buffi(3,i) .le. -199.9) goto 40
        if (buffk(3,kok) .le. -199.9) goto 40
        call fixadf (buffi(3,i),buffk(3,kok),d2)
        d1 = 10.0*juchi
        if (jok .gt. 0) then
          if (buffj(3,jok) .gt. -199.9) then
            call fixadf (buffi(3,i),buffj(3,jok),d1)
          end if
        end if
c
        nchi = nchi + 1
        chirms = chirms + d2*d2
        if (abs(d1) .gt. juchi) then
          mchi = mchi + 1
          mchirms = mchirms + d2*d2
          if (abs(d2) .le. juchi) nclochi = nclochi + 1
        end if
c
   40   continue
c
  100   continue
      end do
c
      if (linsert) then
        call anains (ninsert,ninres,buffl,buffm,
     +    rmsglo(ninsert),rmsloc(ninsert))
        numins (ninsert) = ninres
        linsert = .false.
      end if
c
      pnin = 100.0*float(nin)/float(ntot)
      if (nin .gt. 0) then
        inrmsd = sqrt (inrmsd / float(nin))
      end if
c
      if (nca .gt. 0) then
        carmsd = sqrt (carmsd / float(nca))
      end if
c
      if (nphi .gt. 0) then
        phirms = sqrt (phirms / float(nphi))
      end if
c
      if (npsi .gt. 0) then
        psirms = sqrt (psirms / float(npsi))
      end if
c
      if (nchi .gt. 0) then
        chirms = sqrt (chirms / float(nchi))
      end if
c
      if (mphi .gt. 0) then
        mphirms = sqrt (mphirms / float(mphi))
      end if
c
      if (mpsi .gt. 0) then
        mpsirms = sqrt (mpsirms / float(mpsi))
      end if
c
      if (mchi .gt. 0) then
        mchirms = sqrt (mchirms / float(mchi))
      end if
c
      pcloca = 0.0
      if (nca.gt.0) pcloca = 100.0 * float(ncloca) / float(nca)
      pclo2 = 0.0
      if (nca2.gt.0) pclo2 = 100.0 * float(ncaclo2) / float(nca2)
      pclophi = 0.0
      if (mphi.gt.0) pclophi = 100.0 * float(nclophi) / float(mphi)
      pclopsi = 0.0
      if (mpsi.gt.0) pclopsi = 100.0 * float(nclopsi) / float(mpsi)
      pclochi = 0.0
      if (mchi.gt.0) pclochi = 100.0 * float(nclochi) / float(mchi)
c
      call xystat (buffb1,buffb2,nca,q1,q2,cordis,
     +  q3,q4,q5,q6,q7,q8)
c
      call xystat (buffb3,buffb4,nin,q1,q2,corall,
     +  q3,q4,q5,q6,q7,q8)
c
      write (*,6900) ntot,nin,pnin,inrmsd,corall,nphi,phirms,
     +               npsi,psirms,nchi,chirms,nca2,ncaclo2,pclo2,
     +  nca,carmsd,cordis,ncloca,pcloca,
     +  mphi,mphirms,nclophi,pclophi,
     +  mpsi,mpsirms,nclopsi,pclopsi,
     +  mchi,mchirms,nclochi,pclochi
c
 6900 format (/' RESULTS :'/
     +  ' Residues included in TARGET  : ',i8//
     +  ' Residues included in MODEL   : ',i8/
     +  '                   Percentage : ',f8.2/
     +  '            CA RMSD for these : ',f8.2/
     +  '      Corr. with CA B-factors : ',f8.2//
     +  ' Residues with valid PHI      : ',i8/
     +  '         RMS(D-PHI) for these : ',f8.2//
     +  ' Residues with valid PSI      : ',i8/
     +  '         RMS(D-PSI) for these : ',f8.2//
     +  ' Residues with valid CHI-1    : ',i8/
     +  '       RMS(D-CHI-1) for these : ',f8.2//
     +  ' Residues close to PARENT     : ',i8/
     +  '  Closer in MODEL than PARENT : ',i8/
     +  '                   Percentage : ',f8.2//
     +  ' Residues not close to PARENT : ',i8/
     +  '            CA RMSD for these : ',f8.2/
     +  '      Corr. with CA B-factors : ',f8.2/
     +  '        Number close in MODEL : ',i8/
     +  '                   Percentage : ',f8.2//
     +  ' Not close in distance & PHI  : ',i8/
     +  '         RMS(D-PHI) for these : ',f8.2/
     +  '        Number close in MODEL : ',i8/
     +  '                   Percentage : ',f8.2//
     +  ' Not close in distance & PSI  : ',i8/
     +  '         RMS(D-PSI) for these : ',f8.2/
     +  '        Number close in MODEL : ',i8/
     +  '                   Percentage : ',f8.2//
     +  ' Not close in distance & CHI1 : ',i8/
     +  '       RMS(D-CHI-1) for these : ',f8.2/
     +  '        Number close in MODEL : ',i8/
     +  '                   Percentage : ',f8.2/
     + )
c
      call jvalut (' Number of insertions :',1,ninsert)
c
      if (ninsert .gt. 0) then
        wavglo = 0.0
        wavloc = 0.0
        j = 0
        k = 0
        write (*,6910) '#','Nres','RMSD(global)','RMSD(local)'
        do i=1,ninsert
          write (*,6920) i,numins(i),rmsglo(i),rmsloc(i)
          wavglo = wavglo + float(numins(i))*rmsglo(i)
          j = j + numins(i)
          if (rmsloc(i) .ge. 0.0) then
            wavloc = wavloc + float(numins(i))*rmsloc(i)
            k = k + numins(i)
          end if
        end do
        if (j .gt. 0) then
          wavglo = wavglo / float(j)
          call fvalut (
     +      ' Weighted average global insert RMSD :',1,wavglo)
        end if
        if (k .gt. 0) then
          wavloc = wavloc / float(k)
          call fvalut (
     +      ' Weighted average local  insert RMSD :',1,wavloc)
        end if
      end if
      write (*,*)
c
 6910 format (2x,a2,3x,a4,2a15)
 6920 format (1x,i3,1x,i6,2(5x,f6.2,4x))
c
      return
      end
c
c
c
      subroutine csetup (imol,ichn,nresi,irptr,rbuf)
c
      include 'lsqman.incl'
c
      integer imol,ichn,nresi,i,i1,i2,jmol,jres
      integer iptr,jptr,ierr,nphi,npsi,nchi,in,ic,icb,icg
      integer ic1,in1
      integer irptr (*)
c
      real rbuf(3,*)
      real tangle,dist
c
      character iichn*1,jchn*1
c
code ...
c
      i1 = iresid ( chnptr(1,ichn,imol) , imol )
      i2 = iresid ( chnptr(2,ichn,imol) , imol )
c
      do i=chnptr(1,ichn,imol),chnptr(2,ichn,imol)
        i2 = max (i2,iresid(i,imol))
      end do
c
ccc      print *,' residues ',i1,i2
c
      iichn = chname(ichn,imol)
      jmol = -1
      jchn = '?'
      jres = -1
      nresi = 0
      nphi = 0
      npsi = 0
      nchi = 0
c
      do i=i1,i2
        call getptr (imol,iichn,i,' CA ',jmol,jchn,jres,
     +               iptr,jptr,ierr)
        if (ierr .ne. 0) goto 100
        nresi = nresi + 1
        irptr (nresi) = iptr
        rbuf (1,nresi) = -999.99
        rbuf (2,nresi) = -999.99
        rbuf (3,nresi) = -999.99
        call getptr (imol,iichn,i,' N  ',jmol,jchn,jres,
     +               in,jptr,ierr)
        if (ierr .ne. 0) goto 100
        call getptr (imol,iichn,i,' C  ',jmol,jchn,jres,
     +               ic,jptr,ierr)
        if (ierr .ne. 0) goto 100
c
c ... PHI ?
c
        call getptr (imol,iichn,i-1,' C  ',jmol,jchn,jres,
     +               ic1,jptr,ierr)
        if (ierr .eq. 0) then
          if (dist(ic1,in,atmxyz(1,1,imol)) .le. 2.0) then
            nphi = nphi + 1
            rbuf (1,nresi) = tangle (ic1,in,iptr,ic,
     +                               atmxyz(1,1,imol))
          end if
        end if
c
c ... PSI ?
c
        call getptr (imol,iichn,i+1,' N  ',jmol,jchn,jres,
     +               in1,jptr,ierr)
        if (ierr .eq. 0) then
          if (dist(ic,in1,atmxyz(1,1,imol)) .le. 2.0) then
            npsi = npsi + 1
            rbuf (2,nresi) = tangle (in,iptr,ic,in1,
     +                               atmxyz(1,1,imol))
          end if
        end if
c
c ... CHI1 ?
c
        call getptr (imol,iichn,i,' CB ',jmol,jchn,jres,
     +               icb,jptr,ierr)
        if (ierr .ne. 0) goto 100
c
        call getptr (imol,iichn,i,' CG ',jmol,jchn,jres,
     +               icg,jptr,ierr)
        if (ierr .eq. 0) goto 10
        call getptr (imol,iichn,i,' SG ',jmol,jchn,jres,
     +               icg,jptr,ierr)
        if (ierr .eq. 0) goto 10
        call getptr (imol,iichn,i,' OG ',jmol,jchn,jres,
     +               icg,jptr,ierr)
        if (ierr .eq. 0) goto 10
        call getptr (imol,iichn,i,' CG1',jmol,jchn,jres,
     +               icg,jptr,ierr)
        if (ierr .eq. 0) goto 10
        call getptr (imol,iichn,i,' OG1',jmol,jchn,jres,
     +               icg,jptr,ierr)
        if (ierr .eq. 0) goto 10
c
        goto 100
c
   10   continue
        nchi = nchi + 1
        rbuf (3,nresi) = tangle (in,iptr,icb,icg,
     +                           atmxyz(1,1,imol))
c
  100   continue
      end do
c
      write (*,6000) name(imol),chname(ichn,imol),
     +  nresi,nphi,npsi,nchi
c
 6000 format (' Molecule/Chain ',a,' / ',a/
     +  ' Amino acids  : ',i8/
     +  ' Nr of PHIs   : ',i8/
     +  ' Nr of PSIs   : ',i8/
     +  ' Nr of CHI-1s : ',i8)
c
      return
      end
c
c
c
      subroutine anains (ni,nres,xyz1,xyz2,rmsglo,rmsloc)
c
      implicit none
c
      integer ni,nres,ierr
c
      real xyz1(3,nres),xyz2(3,nres),rtdum(12)
      real rmsglo,rmsloc
c
code ...
c
      rmsloc = -1.0
      call lsqrms (nres,xyz1,xyz2,rmsglo)
      if (nres .ge. 3) then
        call lsqgjk (xyz1,xyz2,nres,rmsloc,rtdum,ierr)
        if (ierr .ne. 0) rmsloc = -1.0
      end if
c
      write (*,6000) ni,nres,rmsglo,rmsloc
 6000 format (' ==> INSERT # ',i3,' : ',i3,' res; RMSD (global) ',
     +  f6.2,' (local) ',f6.2)
c
      return
      end
c
c
c
      subroutine morphm (imol,irange,jmol,jrange,nstep,fbase,
     +                   mode,oid,krange,tarcut,ierr)
c
c ... MORPHM - morph transition between two conformational states
c
c ... buffers used: buffi,buffj,buffk,buffl,buffm,buffb1,buffb2,
c                   buffb3,ptri,ptrij,sptrij
c
      include 'lsqman.incl'
c
      integer maxopt
      parameter (maxopt = 25)
c
      real fact,xave,xsdv,xmin,xmax,xtot,xrmsd,tarcut
      real distce,dist
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j,k,kk
      integer idum,jdum,kdum,ldum,nuse,iptr,jptr,ires,jres
      integer leng1,nstep,junit,kunit,jlim,kptr,iii,ioff,ica
c
      logical xinter,lhydro,mainch,lmdist,lzone
c
      character irange*(*),jrange*(*),fbase*(*),ichn*1,jchn*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20,key*20
      character filenm*256,line*128,mode*1,odraw*10,oid*1
      character krange*(*),objnam*10
c
code ...
c
      ierr = -1
c
      if (mode .eq. 'C') then
        call prompt ('0Morph in Cartesian coordinate space')
        call prompt (
     +    ' NOTE - Molecules must have been superimposed !!!')
      else
        mode = 'I'
        call prompt ('0Morph in internal coordinate space')
      end if
c
      call asciut (' Atom type(s) :',natype,atypes)
c
      if (mode .eq. 'I') then
        if ( atypes(1).eq.'ALL ' .or.
     +       atypes(1).eq.'NONH' .or.
     +       atypes(1).eq.'SIDE' ) then
          call errcon ('Invalid central atom type')
          return
        end if
        if (atypes(1) .eq. 'TRAC') then
          call prompt (' Morph CA and side-chain atoms')
        else
          call prompt (' Only first atom type used in morphing !')
        end if
      end if
c
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range1')
        return
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range2')
        return
      end if
c
      if (nopt1 .ne. nopt2) then
        call errcon ('Different nr of zones')
        return
      end if
c
      nuse = 0
      iptr = 0
      jptr = 0
      xrmsd = 0.0
      lmdist = .false.
c
c ... for internal coordinate morphs, only use one atom per residue
c
      jlim = natype
      if (mode .eq. 'I') jlim = 1
c
c ... central atom flag
c
      ica = 1
c
c ... if CA plus side chains, generate three dummy anchor atoms
c
      ioff = 0
      if (mode .eq. 'I' .and. atypes(1) .eq. 'TRAC') then
        ioff = 3
c
        nuse = nuse + 1
        ptri (nuse) = -1
        sptrij (nuse) = 1
        do k=1,3
          buffi(k,nuse) = 3.8
          buffj(k,nuse) = 3.8
        end do
c
        nuse = nuse + 1
        ptri (nuse) = -1
        sptrij (nuse) = 1
        buffi(1,nuse) = 7.6
        buffj(1,nuse) = 7.6
        do k=2,3
          buffi(k,nuse) = 3.8
          buffj(k,nuse) = 3.8
        end do
c
        nuse = nuse + 1
        ptri (nuse) = -1
        sptrij (nuse) = 1
        buffi(1,nuse) = 7.6
        buffj(1,nuse) = 7.6
        buffi(2,nuse) = 7.6
        buffj(2,nuse) = 7.6
        buffi(3,nuse) = 3.8
        buffj(3,nuse) = 3.8
c
c        call fvalut (' Dummy I 1 :',3,buffi(1,1))
c        call fvalut (' Dummy I 2 :',3,buffi(1,2))
c        call fvalut (' Dummy I 3 :',3,buffi(1,3))
c        call fvalut (' Dummy J 1 :',3,buffj(1,1))
c        call fvalut (' Dummy J 2 :',3,buffj(1,2))
c        call fvalut (' Dummy J 3 :',3,buffj(1,3))
c
      end if       
c
      do i=1,nopt1
c
        ichn = optpa1(i)(1:1)
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon (
ccc     +      'Invalid chain id in range1 : '//optpa1(i))
ccc          return
        end if
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
c
c ... 971111 - allow e.g. "A78" to mean "A78-78"
c
           k = length(optpa1(i))
           if (ichn .eq. ' ') then
             optpa1(i) = optpa1(i)(1:k) // '-' // optpa1(i)(1:k)
           else
             optpa1(i) = optpa1(i)(1:k) // '-' // optpa1(i)(2:k)
           end if
           j = k + 1
ccc          call errcon ('Not a zone in range1 : '//optpa1(i))
ccc          return
        end if
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon (
     +      'Invalid zone in range1 : '//optpa1(i))
          return
        end if
c
        jchn = optpa2(i)(1:1)
        if (index('1234567890',jchn) .gt. 0) then
          jchn = ' '
ccc        else if (jchn .lt. 'A' .or. jchn .gt. 'Z') then
ccc          call errcon (
ccc     +      'Invalid chain id in range2 : '//optpa2(i))
ccc          return
        end if
        j = index (optpa2(i),'-')
        if (j .le. 0) j = index (optpa2(i),':')
        if (j .gt. 0) then
          call errcon (
     +      'Zone ignored in range2 : '//optpa2(i))
          optpa2(i) = optpa2(i)(1:j-1)
        end if
        if (jchn .eq. ' ') then
          call str2i (optpa2(i)(1:),kdum,ierr)
        else
          call str2i (optpa2(i)(2:),kdum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        ldum = kdum + (jdum - idum)
c
        do ires=idum,jdum
          jres = kdum + (ires - idum)
          do j=1,jlim
            call getptr (imol,ichn,ires,atypes(j),
     +                   jmol,jchn,jres,iptr,jptr,ierr)
            if (ierr .eq. 0) then
              nuse = nuse + 1
              ptri (nuse) = iptr
              sptrij (nuse) = ica - 1
              if (atmnam(iptr,imol) .eq. ' CA ') sptrij(nuse) = ica
              do k=1,3
                buffi(k,nuse) = atmxyz(k,iptr,imol)
                buffj(k,nuse) = atmxyz(k,jptr,jmol)
              end do
              xrmsd = xrmsd + distce(buffi(1,nuse),buffj(1,nuse))**2
c
c ... get others if ALL or NONH or SIDE or TRAC
c
              if ( natype .eq. 1 .and.
     +            (atypes(1).eq.'ALL ' .or.
     +             atypes(1).eq.'NONH' .or.
     +             atypes(1).eq.'SIDE' .or.
     +             atypes(1).eq.'TRAC') ) then
                do k=iptr+1,natoms(imol)
                  if (iresid(k,imol).ne.iresid(iptr,imol) .or.
     +                achain(k,imol).ne.achain(iptr,imol)) then
                    kptr = k - 1
                    goto 6511
                  else if (k .eq. natoms(imol)) then
                    kptr = k
                    goto 6511
                  end if
                end do
                kptr = iptr
 6511           continue
c
                if (kptr .gt. iptr) then
                  do k=iptr+1,kptr
c
c ... 960729 - find same atom in other molecule
c
                    call getptr (imol,ichn,ires,atmnam(k,imol),
     +                           jmol,jchn,jres,iii,jptr,ierr)
c
                    if (ierr .ne. 0) goto 6436
c
                      if (
     +   (atypes(1) .eq. 'ALL ')
     +      .or.
     +   (atypes(1) .eq. 'NONH' .and. (.not. lhydro(atmnam(k,imol))))
     +      .or.
     +   (atypes(1) .eq. 'SIDE' .and. (.not. mainch(atmnam(k,imol)))
     +                          .and. (.not. lhydro(atmnam(k,imol))))
     +      .or.
     +   (atypes(1) .eq. 'TRAC' .and. (.not. lhydro(atmnam(k,imol)))
     +       .and. (atmnam(k,imol) .eq. ' CA ' .or.
     +              (.not. mainch(atmnam(k,imol))) ))
     +      ) then
                nuse = nuse + 1
                ptri (nuse) = k
                sptrij (nuse) = ica - 1
                if (atmnam(k,imol) .eq. ' CA ') sptrij(nuse) = ica
                do kk=1,3
                  buffi(kk,nuse) = atmxyz(kk,k,imol)
                  buffj(kk,nuse) = atmxyz(kk,jptr,jmol)
                end do
                xrmsd = xrmsd +
     +            distce(buffi(1,nuse),buffj(1,nuse))**2
c
                      end if
 6436               continue
                  end do
                end if
              end if
c
            end if
          end do
        end do
c
      end do
c
      call jvalut (' Nr of atoms matched :',1,nuse-ioff)
      if (nuse .lt. ioff+3) then
        call errcon ('Fewer than 3 atoms; cannot morph')
        return
      end if
      if (nuse .gt. maxatm/3) then
        call errcon ('Too many atoms; cannot morph')
        return
      end if
c
      xrmsd = sqrt (xrmsd / float(nuse-ioff))
      call fvalut (' RMSD matching atoms :',1,xrmsd)
      if (mode .eq. 'C' .and. xrmsd .gt. 3.0) then
        if (xrmsd .le. 10.0) then
          call prompt ('0WARNING - RMSD high (> 3 A) !!!')
        else
          call prompt ('0WARNING - RMSD very high (> 10 A)!!!')
        end if
        call prompt (
     +    ' Are you sure the molecules have been superimposed ???')
      end if
c
      if (mode .eq. 'I') then
c
c ... convert to internal coordinates
c
        if (atypes(1) .ne. 'TRAC') then
          call cart2c (1,nuse,buffi,buffk,ptrij)
          call cart2c (2,nuse,buffj,buffl,ptrij)
          do i=1,nuse
            sptrij (i) = ica
          end do
        else
c
c ... now get more reasonable dummy atoms
c
c        call fvalut (' Dummy I 1 :',3,buffi(1,1))
c        call fvalut (' Dummy I 2 :',3,buffi(1,2))
c        call fvalut (' Dummy I 3 :',3,buffi(1,3))
c        call fvalut (' Dummy J 1 :',3,buffj(1,1))
c        call fvalut (' Dummy J 2 :',3,buffj(1,2))
c        call fvalut (' Dummy J 3 :',3,buffj(1,3))
c
          buffi (3,3) = buffi (3,4) - 0.5
          buffi (2,3) = buffi (2,4) - 0.5
          buffi (1,3) = buffi (1,4) - 0.5
c
          buffi (3,2) = buffi (3,3)
          buffi (2,2) = buffi (2,3) - 0.5
          buffi (1,2) = buffi (1,3) - 0.5
c
          buffi (3,1) = buffi (3,2)
          buffi (2,1) = buffi (2,2)
          buffi (1,1) = buffi (1,2) - 0.5
c
          buffj (3,3) = buffj (3,4) - 0.5
          buffj (2,3) = buffj (2,4) - 0.5
          buffj (1,3) = buffj (1,4) - 0.5
c
          buffj (3,2) = buffj (3,3)
          buffj (2,2) = buffj (2,3) - 0.5
          buffj (1,2) = buffj (1,3) - 0.5
c
          buffj (3,1) = buffj (3,2)
          buffj (2,1) = buffj (2,2)
          buffj (1,1) = buffj (1,2) - 0.5
c
c        call fvalut (' Dummy I 1 :',3,buffi(1,1))
c        call fvalut (' Dummy I 2 :',3,buffi(1,2))
c        call fvalut (' Dummy I 3 :',3,buffi(1,3))
c        call fvalut (' Dummy J 1 :',3,buffj(1,1))
c        call fvalut (' Dummy J 2 :',3,buffj(1,2))
c        call fvalut (' Dummy J 3 :',3,buffj(1,3))
c
c ... get topology
c
          call trac2c (1,nuse,buffi,sptrij,buffk,ptrij,ica)
          call trac2c (2,nuse,buffj,sptrij,buffl,ptrij,ica)
ccc          call xprint (25,buffi,buffk,ptrij,'MOL 1')
ccc          call xprint (25,buffj,buffl,ptrij,'MOL 2')
        end if
c
c ... map torsion angles properly
c
        do i=1,nuse
          if ( (buffk(3,i)-buffl(3,i)) .lt. -180.0) then
            buffl (3,i) = buffl (3,i) - 360.0
          else if ( (buffk(3,i)-buffl(3,i)) .gt. 180.0) then
            buffl (3,i) = buffl (3,i) + 360.0
          end if
          buffb1 (i) = abs(buffk(3,i)-buffl(3,i))
        end do
ccc        call xprint (25,buffj,buffl,ptrij,'MOL 1 180')
c
      else
c
        do i=1,nuse
          buffb1 (i) = distce (buffi(1,i),buffj(1,i))
        end do
c
      end if
c
      call xstats (buffb1(ioff+1),nuse-ioff,xave,xsdv,xmin,xmax,xtot)
      if (mode .eq. 'I') then
        call prompt (
     +    '0B-factors replaced by torsion angle differences:')
        call fvalut (' Average :',1,xave)
        call fvalut (' St.dev. :',1,xsdv)
        call fvalut (' Minimum :',1,xmin)
        call fvalut (' Maximum :',1,xmax)
        do i=1,nuse
          buffb3 (i) = buffb1 (i)
        end do
c
        if (atypes(1) .eq. 'TRAC') then
c
c ... use connectivity PTRIJ to figure out how to colour the
c     chemical bonds; first side chains, then CA trace
c
          do i=1,nuse
            buffb2 (i) = 0.0
          end do
          do i=1,nuse
            j = 1 + (i-1)*3
c
c ... skip if any antecedent undefined
c
            if (ptrij(j) .le. 0 .or. ptrij(j+1) .le. 0 .or.
     +          ptrij(j+2) .le. 0) goto 1200
c
c ... don't mix CA and SC torsions
c
              if (sptrij(ptrij(j+1)) .eq. ica .and.
     +            sptrij(ptrij(j+2)) .eq. ica .and.
     +            sptrij(ptrij(j))   .eq. ica .and.
     +            sptrij(i)          .ne. ica) goto 1200
c
c ... chemically connected CA-CA-CA-CA type
c
              if (sptrij(ptrij(j+1)) .eq. ica .and.
     +            sptrij(ptrij(j+2)) .eq. ica .and.
     +            sptrij(ptrij(j))   .eq. ica .and.
     +            sptrij(i)          .eq. ica) then
                if (dist(i,ptrij(j),buffi) .le. 4.5 .and.
     +              dist(ptrij(j),ptrij(j+1),buffi) .le. 4.5 .and.
     +              dist(ptrij(j+1),ptrij(j+2),buffi) .le. 4.5) then
                  buffb2 (ptrij(j)) =
     +                max(buffb2(ptrij(j)),buffb1(i))
                  buffb2 (ptrij(j+1)) =
     +                max(buffb2(ptrij(j+1)),buffb1(i))
                end if
                goto 1200
              end if
c
c ... chemically connected X-Y-SC-SC type (X-Y may be CA-CA)
c
              if (sptrij(ptrij(j))   .ne. ica .and.
     +            sptrij(i)          .ne. ica) then
                if (dist(i,ptrij(j),buffi) .le. 2.0 .and.
     +              dist(ptrij(j),ptrij(j+1),buffi) .le. 2.0) then
                  buffb2 (ptrij(j)) =
     +                max(buffb2(ptrij(j)),buffb1(i))
                  if (sptrij(ptrij(j+1))   .ne. -1) then
                    buffb2 (ptrij(j+1)) =
     +                max(buffb2(ptrij(j+1)),buffb1(i))
                  end if
                end if
                goto 1200
              end if
c
 1200       continue
          end do
          do i=1,nuse
            buffb1 (i) = buffb2 (i)
          end do
c
c ... re-map Bs such that actual torsion bonds are coloured
c
        else
          do i=1,nuse-2
            buffb1(i) = buffb1(i+2)
          end do
          buffb1 (nuse-1) = 0.0
          buffb1 (nuse) = 0.0
          do i=1,nuse-2
            if (buffb1(i+1) .gt. buffb1(i)) buffb1(i)=buffb1(i+1)
          end do
        end if
c
c ... need to fix large central-atom torsions ?
c
        if (xmax .gt. tarcut) then
          call prompt (
     +      '0WARNING - Large torsion angle change(s) !')
          call fixcat (nuse,buffb3,buffi,buffj,ptrij,
     +                 sptrij,ica,tarcut,k,kk)
          call jvalut (' Nr fixed  :',1,kk)
          call jvalut (' Unfixable :',1,k)
          if (kk .gt. 0) then
c
c ... change CA_MXDST entry !
c
            lmdist = .true.
            call cart2c (2,nuse,buffi,buffk,ptrij)
            call cart2c (2,nuse,buffj,buffl,ptrij)
            do i=1,nuse
              if ( (buffk(3,i)-buffl(3,i)) .lt. -180.0) then
                buffl (3,i) = buffl (3,i) - 360.0
              else if ( (buffk(3,i)-buffl(3,i)) .gt. 180.0) then
                buffl (3,i) = buffl (3,i) + 360.0
              end if
            end do
            call xstats (buffb3(ioff+1),nuse-ioff,
     +                   xave,xsdv,xmin,xmax,xtot)
            call prompt (
     +        '0Updated statistics for torsion angle differences:')
            call fvalut (' Average :',1,xave)
            call fvalut (' St.dev. :',1,xsdv)
            call fvalut (' Minimum :',1,xmin)
            call fvalut (' Maximum :',1,xmax)
          end if
        end if
c
      else
        call prompt (
     +    '0B-factors replaced by Cartesian distances:')
        call fvalut (' Average :',1,xave)
        call fvalut (' St.dev. :',1,xsdv)
        call fvalut (' Minimum :',1,xmin)
        call fvalut (' Maximum :',1,xmax)
      end if
c
      call prompt ('0Morphing ...')
c
c ... create LSQMAN macro
c
      junit = iunit + 1
      write (filenm,'(a,a)') fbase(1:leng1(fbase)),'.lsqmac'
      call remspa (filenm)
      call textut (' Creating LSQMAN macro :',filenm)
      call xopxua (junit,filenm,xinter(),ierr)
      if (ierr .ne. 0) goto 9900
      call stamp (line)
      write (junit,4712,err=9700) '!     ',line(1:leng1(line))
      if (mode .eq. 'I' .and. atypes(1) .eq. 'TRAC') then
        write (junit,'(9a)') 'atom_type CA'
      end if
c
c ... create first O macro
c
      kunit = iunit + 2
      write (filenm,'(a,a)') fbase(1:leng1(fbase)),'_read.omac'
      call remspa (filenm)
      call textut (' Creating first O macro :',filenm)
      call xopxua (kunit,filenm,xinter(),ierr)
      if (ierr .ne. 0) goto 9900
      call stamp (line)
      write (kunit,4712,err=9600) '!     ',line(1:leng1(line))
      if (mode .eq. 'I' .and. atypes(1) .eq. 'TRAC') then
        write (kunit,'(a)') 'conn o.dat'
        odraw = 'zo ; end'
        lzone = .true.
      else
        odraw = 'ca ; end'
        lzone = .false.
      end if
c
c ... generate the intermediates
c
      do i=1,nstep
c
        call jvalut (' Morph step number :',1,i)
        write (filenm,'(a,a,i8,a)') fbase(1:leng1(fbase)),
     +    '_',i,'.pdb'
        call remspa (filenm)
        call xopxua (iunit,filenm,xinter(),ierr)
        if (ierr .ne. 0) goto 9900
        call stamp (line)
        write (iunit,4712,err=9800) 'REMARK',line(1:leng1(line))
        write (line,*) 'LSQMAN MORPH STEP ',i,' OF ',nstep
        call pretty (line)
        write (iunit,4712,err=9800) 'REMARK',line(1:leng1(line))
c
        write (key,'(i8,a)') i,oid
        call remspa (key)
c
c ... LSQMAN macro
c
        if (i .eq. 1) then
          write (junit,'(9a)') 'read 1m ',filenm(1:leng1(filenm))
        else
          write (junit,'(9a)') 'read ',key(1:leng1(key)),' ',
     +      filenm(1:leng1(filenm))
          write (junit,'(9a)') 'explicit 1m ',
     +      krange(1:leng1(krange)),' ',key(1:leng1(key)),' ',
     +      krange(1:leng1(krange))
          write (junit,'(9a)') 'set dist 3.0'
          write (junit,'(9a)') 'improve 1m ',
     +      krange(1:leng1(krange)),' ',key(1:leng1(key)),' *'
          write (junit,'(9a)') 'set dist 2.0'
          write (junit,'(9a)') 'improve 1m ',
     +      krange(1:leng1(krange)),' ',key(1:leng1(key)),' *'
          write (junit,'(9a)') 'set dist 1.0'
          write (junit,'(9a)') 'improve 1m ',
     +      krange(1:leng1(krange)),' ',key(1:leng1(key)),' *'
          write (junit,'(9a)') 'apply 1m ',key(1:leng1(key))
          write (junit,'(9a)') 'write ',key(1:leng1(key)),' ',
     +      filenm(1:leng1(filenm))
          write (junit,'(9a)') 'delete ',key(1:leng1(key))
        end if
        if (i .eq. nstep) then
          write (junit,'(9a)') 'delete 1m'
          write (junit,'(9a)') '! RESET MATCHING DISTANCE CUT-OFF !!!'
          write (junit,'(9a)') '! FOR EXAMPLE: set dist 3.8'
          write (junit,'(9a)') '! OR: set reset'
          if (mode .eq. 'I' .and. atypes(1) .eq. 'TRAC') then
            write (junit,'(9a)') 'atom_type TRACE'
          end if
        end if
c
c ... O macro
c
c        write (kunit,'(9a)') 'sam_atom_in ',filenm(1:leng1(filenm)),
c     +    ' ',key(1:leng1(key))
ccc     +    ' ',key(1:leng1(key)),' pdb'
c
        write (kunit,'(9a)') 'read stereo_chem.odb'
        write (kunit,'(9a)') 'pdb_read ',filenm(1:leng1(filenm)),
     +    ' ',key(1:leng1(key)),' y y'
        write (kunit,'(9a)') 'mol ',key(1:leng1(key))
        if (lmdist) then
          write (kunit,'(9a)') 'db_set_dat ',key(1:leng1(key)),
     +      '_MOLECULE_CA_MXDST 1 1 15.0'
        end if
        write (kunit,'(9a)') 'paint_ramp atom_b ; blue red ',odraw
c
c ... morph coordinates
c
        fact = float (i-1) / float (nstep-1)
        if (mode .eq. 'I') then
          do j=1,nuse
            do k=1,3
              buffm(k,j)=buffk(k,j)+fact*(buffl(k,j)-buffk(k,j))
            end do
          end do
        else
          do j=1,nuse
            do k=1,3
              buffk(k,j)=buffi(k,j)+fact*(buffj(k,j)-buffi(k,j))
            end do
          end do
        end if
c
c ... convert back to Cartesian coordinates
c
        if (mode .eq. 'I') then
          call c2cart (nuse,buffi,buffm,ptrij)
ccc          call xprint (25,buffi,buffm,ptrij,'MOL 1 180')
        end if
c
c ... write the PDB file
c
 4711 format (a6,a)
 4712 format (a6,1x,a)
 4713 format (i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2,6x,a4)
c
        key = 'ATOM  '
c
        if (mode .eq. 'I') then
c
          do k=ioff+1,nuse
            kk = ptri(k)
            write (line,4713,err=9800) k-ioff,
     +        atmnam(kk,imol),resnam(kk,imol),achain(kk,imol),
     +        iresid(kk,imol),buffi(1,k),buffi(2,k),
     +        buffi(3,k),qatom(kk,imol),buffb1(k),
     +        axplor(kk,imol)
            write (iunit,4711,err=9800) key,line(1:leng1(line))
          end do
c
        else
c
          do k=ioff+1,nuse
            kk = ptri(k)
            write (line,4713,err=9800) k-ioff,
     +        atmnam(kk,imol),resnam(kk,imol),achain(kk,imol),
     +        iresid(kk,imol),buffk(1,k),buffk(2,k),
     +        buffk(3,k),qatom(kk,imol),buffb1(k),
     +        axplor(kk,imol)
            write (iunit,4711,err=9800) key,line(1:leng1(line))
          end do
c
        end if
c
        write (iunit,'(a)',err=9800) 'END'
        close (iunit)
c
      end do
c
c ... centre_zone on first molecule in O
c
      write (key,'(i8,a)') 1,oid
      call remspa (key)
      kk = ptri(ioff+1)
      write (line,'(a1,i4)') achain(kk,imol),iresid(kk,imol)
      call remspa (line)
      kk = ptri(nuse)
      write (line(10:),'(a1,i4)') achain(kk,imol),iresid(kk,imol)
      call remspa (line(10:))
      call upcase (line)
      call pretty (line)
      write (kunit,'(9a)') 'ce_zo ',key(1:leng1(key)),' ',
     +  line(1:leng1(line))
c
      write (kunit,'(9a)') 'print Morphed by LSQMAN - (C) ',
     +  '1998-2005 - G.J. Kleywegt - Uppsala'
      write (kunit,'(9a)') 'message Morphed by LSQMAN - (C) ',
     +  '1998-2005 - G.J. Kleywegt - Uppsala'
      close (kunit)
c
c ... create second O macro
c
      write (filenm,'(a,a)') fbase(1:leng1(fbase)),'_morph.omac'
      call remspa (filenm)
      call textut (' Creating second O macro :',filenm)
      call xopxua (kunit,filenm,xinter(),ierr)
      if (ierr .ne. 0) goto 9900
      call stamp (line)
      write (kunit,4712,err=9600) '!     ',line(1:leng1(line))
c
      do i=1,nstep
        write (key,'(i8,a)') i,oid
        call remspa (key)
        write (kunit,'(9a)') 'mol ',key(1:leng1(key)),' ',odraw
        j = i - 1
        if (j .lt. 1) j = nstep
        write (key,'(i8,a)') j,oid
        call remspa (key)
        write (kunit,'(9a)') 'del_obj ',key(1:leng1(key)),' ;'
      end do
      write (kunit,*)
c
      do i=nstep,1,-1
        write (key,'(i8,a)') i,oid
        call remspa (key)
        write (kunit,'(9a)') 'mol ',key(1:leng1(key)),' ',odraw
        j = i + 1
        if (j .gt. nstep) j = 1
        write (key,'(i8,a)') j,oid
        call remspa (key)
        write (kunit,'(9a)') 'del_obj ',key(1:leng1(key)),' ;'
      end do
      close (kunit)
c
c ... create third O macro
c
      write (filenm,'(a,a)') fbase(1:leng1(fbase)),'_plot.omac'
      call remspa (filenm)
      call textut (' Creating third O macro :',filenm)
      call xopxua (kunit,filenm,xinter(),ierr)
      if (ierr .ne. 0) goto 9900
      call stamp (line)
      write (kunit,'(9a)') '! ',line(1:leng1(line))
      write (kunit,'(9a)')
     +  '! Select viewpoint, zoom, etc. before running this macro'
      write (kunit,'(9a)')
     +  '!'
c
      write (kunit,'(9a)')
     +  'plot_setup ',fbase(1:leng1(fbase)),'_morph.plt 1 no no'
      if (lzone) then
        write (kunit,'(9a)')
     +    'sketch_setup stick smooth 0.1 8'
      else
        write (kunit,'(9a)')
     +    'sketch_setup stick smooth 0.3 8'
      end if
      write (kunit,'(9a)')
     +  'sketch_setup sphere smooth 0'
      write (kunit,'(9a)')
     +  'db_create .cpk_radii 110 r'
      if (lzone) then
        write (kunit,'(9a)')
     +    'db_set_dat .cpk_radii ; 0.2'
      else
        write (kunit,'(9a)')
     +    'db_set_dat .cpk_radii ; 0.3'
      end if
      write (kunit,'(9a)')
     +  '!'
      write (kunit,'(9a)')
     +  'plot_on'
c
      do i=1,nstep
        write (kunit,'(9a)') '!'
        write (key,'(i8,a)') i,oid
        call remspa (key)
        if (i .gt. 99) then
          write (objnam,'(a2,i8)') 'FM',i
        else
          write (objnam,'(a2,i2.2)') 'FM',i
        end if
        call remspa (objnam)
        write (kunit,'(9a)') 'mol ',key(1:leng1(key)),
     +    ' obj ',objnam(1:leng1(objnam)),' ',odraw
        write (kunit,'(9a)')
     +    'sketch_stick ',objnam(1:leng1(objnam))
        write (kunit,'(9a)')
     +    'sketch_cpk ',objnam(1:leng1(objnam))
        write (kunit,'(9a)')
     +    'del_obj ',objnam(1:leng1(objnam)),
     +    ' S_',objnam(1:leng1(objnam)),
     +    ' CPK_',objnam(1:leng1(objnam)),' ;'
      end do
c
      write (kunit,'(9a)')
     +  '!'
      write (kunit,'(9a)')
     +  'plot_off'
c
      write (kunit,'(9a)')
     +  '!'
      write (kunit,'(9a)') 'print Morphed by LSQMAN - (C) ',
     +  '1998-2005 - G.J. Kleywegt - Uppsala'
      write (kunit,'(9a)') 'message Morphed by LSQMAN - (C) ',
     +  '1998-2005 - G.J. Kleywegt - Uppsala'
      close (kunit)
c
      close (junit)
      call prompt ('0Morphing done ... now execute (@) the new')
      write (filenm,'(a,a)') fbase(1:leng1(fbase)),'.lsqmac'
      call remspa (filenm)
      call textut (' LSQMAN macro called :',filenm)
      call prompt (' to superimpose the morphed models ...')
c
      ierr = 0
      return
c
c ... errors
c
 9600 continue
      call errcon ('While writing O macro')
      goto 9900
c
 9700 continue
      call errcon ('While writing LSQMAN macro')
      goto 9900
c
 9800 continue
      call errcon ('While writing PDB file')
ccc      goto 9900
c
 9900 continue
      close (iunit)
      close (junit)
      close (kunit)
      return
c
      end
c
c
c
      subroutine cart2c (mode,nuse,xyz,int,ref)
c
c ... CART2C - convert Cartesian into internal coordinates, assuming
c              one continuous chain of atoms (e.g., CA atoms)
c
      implicit none
c
      integer nuse,i,mode
      integer ref(3,*)
c
      real xyz(3,*),int(3,*)
      real dist,angle,tangle
c
code ...
c
      if (mode .eq. 2) goto 500
c
      if (nuse .lt. 1) return
c
      int (1,1) = 0.0
      ref (1,1) = -1
      int (2,1) = 0.0
      ref (2,1) = -1
      int (3,1) = 0.0
      ref (3,1) = -1
      if (nuse .lt. 2) return
c
      int (1,2) = dist (1,2,xyz)
      ref (1,2) = 1
      int (2,2) = 0.0
      ref (2,2) = -1
      int (3,2) = 0.0
      ref (3,2) = -1
      if (nuse .lt. 3) return
c
      int (1,3) = dist (2,3,xyz)
      ref (1,3) = 2
      int (2,3) = angle (1,2,3,xyz)
      ref (2,3) = 1
      int (3,3) = 0.0
      ref (3,3) = -1
      if (nuse .lt. 4) return
c
      do i=4,nuse
        int (1,i) = dist(i-1,i,xyz)
        ref (1,i) = i-1
        int (2,i) = angle(i-2,i-1,i,xyz)
        ref (2,i) = i-2
        int (3,i) = tangle(i-3,i-2,i-1,i,xyz)
        ref (3,i) = i-3
      end do
c
      return
c
c ... subsequent molecule
c
  500 continue
c
      do i=1,nuse
        if (ref(1,i) .le. 0) then
          int (1,i) = 0.0
        else
          int (1,i) = dist(i,ref(1,i),xyz)
          if (ref(2,i) .le. 0) then
            int (2,i) = 0.0
          else
            int (2,i) = angle(i,ref(1,i),ref(2,i),xyz)
            if (ref(3,i) .le. 0) then
              int (3,i) = 0.0
            else
              int (3,i) = tangle(i,ref(1,i),ref(2,i),ref(3,i),xyz)
            end if
          end if
        end if
      end do
c
c
      return
      end
c
c
c
      subroutine trac2c (mode,nuse,xyz,typ,int,ref,ica)
c
c ... TRAC2C - convert Cartesian into internal coordinates, assuming
c              one continuous chain of atoms (e.g., CA atoms; TYP=1)
c              PLUS side chains (TYP=0)
c     NOTE THAT THE FIRST THREE ATOMS MUST BE (DUMMY) CA ATOMS !!!
c     MODE = 1 -> first molecule; define connectivity
c            2 -> subsequent molecule; use previous connectivity
c
      implicit none
c
      integer nuse,i,mode,iprev,icnt,jcnt,j,ica
      integer typ(*),ref(3,*)
c
      real xyz(3,*),int(3,*)
      real x,xmin
      real dist,angle,tangle
c
code ...
c
      if (mode .eq. 2) goto 500
c
      if (nuse .lt. 3) return
c
      iprev = -1
      icnt = 0
c
      int (1,1) = 0.0
      ref (1,1) = -1
      int (2,1) = 0.0
      ref (2,1) = -1
      int (3,1) = 0.0
      ref (3,1) = -1
      icnt = icnt + 1
      typ (1) = typ(1)-2
c
      int (1,2) = dist (1,2,xyz)
      ref (1,2) = 1
      int (2,2) = 0.0
      ref (2,2) = -1
      int (3,2) = 0.0
      ref (3,2) = -1
      icnt = icnt + 1
      typ (2) = typ(2)-2
c
      int (1,3) = dist (2,3,xyz)
      ref (1,3) = 2
      int (2,3) = angle (1,2,3,xyz)
      ref (2,3) = 1
      int (3,3) = 0.0
      ref (3,3) = -1
      icnt = icnt + 1
      typ (3) = typ(3)-2
c
      ica = ica - 2
c
      if (nuse .lt. 4) return
c
c ... first define the central atoms (CA)
c
      iprev = 3
      do i=4,nuse
        if (typ(i) .eq. 1) then
          icnt = icnt + 1
          int (1,i) = dist (i,iprev,xyz)
          ref (1,i) = iprev
          int (2,i) = angle(ref(1,iprev),iprev,i,xyz)
          ref (2,i) = ref(1,iprev)
          int (3,i) = tangle(ref(2,iprev),ref(1,iprev),
     +                       iprev,i,xyz)
          ref (3,i) = ref(2,iprev)
          iprev = i
          typ (i) = typ(i)-2
        end if
      end do
      call jvalut (' Nr of CA atoms :',1,icnt-3)
c
c ... now do the side chain atoms close to CAs
c
      icnt = 0
   30 continue
      jcnt = 0
      do i=1,nuse
        if (typ(i) .lt. 0) goto 100
c
c ... find nearest bonded atom (< 2 A)
c
        do j=4,nuse
          if (typ(j) .lt. 0 .and. ref(1,j) .gt. 0 .and.
     +        ref(2,j) .gt. 0) then
            x = dist(i,j,xyz)
            if (x .lt. 2.0) then
              iprev = j
              goto 20
            end if
          end if
        end do
        goto 100
c
   20   continue
        jcnt = jcnt + 1
        icnt = icnt + 1
        int (1,i) = dist (i,iprev,xyz)
        ref (1,i) = iprev
        int (2,i) = angle(ref(1,iprev),iprev,i,xyz)
        ref (2,i) = ref(1,iprev)
        int (3,i) = tangle(ref(2,iprev),ref(1,iprev),
     +                     iprev,i,xyz)
        ref (3,i) = ref(2,iprev)
        typ (i) = typ(i)-2
c
  100   continue
      end do
      if (jcnt .gt. 0) goto 30
c
c ... now do the rest
c
  110 continue
      jcnt = 0
      do i=1,nuse
        if (typ(i) .lt. 0) goto 200
        jcnt = jcnt + 1
c
c ... find nearest atom
c
        xmin = 9999.99
        iprev = -1
        do j=4,nuse
          if (typ(j) .lt. 0 .and. ref(1,j) .gt. 0 .and.
     +        ref(2,j) .gt. 0) then
            x = dist(i,j,xyz)
            if (x .lt. xmin .or. iprev .le. 0) then
              iprev = j
              xmin = x
            end if
          end if
        end do
        if (iprev .le. 0) then
          do j=1,3
            if (typ(j) .lt. 0 .and. ref(1,j) .gt. 0 .and.
     +          ref(2,j) .gt. 0) then
              x = dist(i,j,xyz)
              if (x .lt. xmin .or. iprev .le. 0) then
                iprev = j
                xmin = x
              end if
            end if
          end do
          if (iprev .le. 0) then
            call errcon ('BUG ??? Tell Gerard !')
            goto 200
          end if
        end if
c
        icnt = icnt + 1
        int (1,i) = dist (i,iprev,xyz)
        ref (1,i) = iprev
        int (2,i) = angle(ref(1,iprev),iprev,i,xyz)
        ref (2,i) = ref(1,iprev)
        int (3,i) = tangle(ref(2,iprev),ref(1,iprev),
     +                     iprev,i,xyz)
        ref (3,i) = ref(2,iprev)
        typ (i) = typ(i)-2
c
  200   continue
      end do
      if (jcnt .gt. 0) goto 30
c
      call jvalut (' Nr of side-chain atoms :',1,icnt)
c
      do i=4,nuse
        if (ref(1,i).eq.ref(2,i) .or.
     +      ref(2,i).eq.ref(3,i) .or.
     +      ref(3,i).eq.ref(1,i)) then
          call errcon ('BUG ! Inconsistent topology - tell Gerard !')
          call ivalut (' ATOM :',1,i)
          call ivalut (' REF  :',3,ref(1,i))
          call fvalut (' INT  :',3,int(1,i))
        end if
      end do
c
      return
c
c ... subsequent molecule
c
  500 continue
c
      do i=1,nuse
        if (ref(1,i) .le. 0) then
          int (1,i) = 0.0
        else
          int (1,i) = dist(i,ref(1,i),xyz)
          if (ref(2,i) .le. 0) then
            int (2,i) = 0.0
          else
            int (2,i) = angle(i,ref(1,i),ref(2,i),xyz)
            if (ref(3,i) .le. 0) then
              int (3,i) = 0.0
            else
              int (3,i) = tangle(i,ref(1,i),ref(2,i),ref(3,i),xyz)
            end if
          end if
        end if
      end do
c
      return
      end
c
c
c
      subroutine xprint (n,x1,x2,ptr,line)
c
      implicit none
c
      integer n,i,k
      integer ptr(3,*)
c
      real x1(3,*),x2(3,*)
c
      character line*(*)
c
code ...
c
      call textut (' LIST :',line)
      do i=1,n
        write (*,6000) i,(x1(k,i),k=1,3),(x2(k,i),k=1,3),
     +    (ptr(k,i),k=1,3)
      end do
c
 6000 format (1x,i3,3f8.3,1x,3f8.3,1x,3i6)
c
      return
      end
c
c
c
      subroutine fixcat (nuse,range,xyz,xyz2,ref,typ,
     +                   ica,tarcut,nno,nok)
c
      implicit none
c
      integer nuse,i,j,nok,nno,imin,jmin,kmin,k,l,ica
      integer ref(3,*),typ(*)
c
      real xyz(3,*),xyz2(3,*),range(*)
      real tangle,x1,x2,r,rmin,tarcut
c
code ...
c
      nok = 0
      nno = 0
c
      if (nuse .lt. 10) return
c
      do i=6,nuse-6
        if (typ(i) .ne. ica) goto 100
        if (range(i) .le. tarcut) goto 100
        rmin = range(i)
        imin = -1
        do j=i-2,4,-1
c
          if (typ(j) .ne. ica) goto 200
          if (ref(1,j) .le. 0) goto 200
          if (typ(ref(1,j)) .ne. ica) goto 200
          if (ref(2,j) .le. 0) goto 200
          if (typ(ref(2,j)) .ne. ica) goto 200
c
          x1 = tangle (i,j,ref(1,j),ref(2,j),xyz)
          x2 = tangle (i,j,ref(1,j),ref(2,j),xyz2)
          if ( (x1-x2) .lt. -180.0) then
            x2 = x2 - 360.0
          else if ( (x1-x2) .gt. 180.0) then
            x2 = x2 + 360.0
          end if
          r = abs(x1-x2)
          if (r .lt. rmin) then
            rmin = r
            imin = j
            if (rmin .lt. tarcut) goto 210
          end if
  200     continue
        end do
  210   continue
        if (imin .gt. 0 .and. rmin .lt. tarcut) then
          write (*,6000) 1,rmin,range(i)
c
c          print *,i,typ(i),imin,typ(imin),ref(1,imin),ref(2,imin)
c
          nok = nok + 1
          ref (1,i) = imin
          ref (2,i) = ref(1,imin)
          ref (3,i) = ref(2,imin)
          range (i) = rmin
          goto 100
        end if
        rmin = range(i)
        imin = -1
        do j=i-1,1,-1
          if (typ(j) .ne. ica) goto 500
          do k=i-1,1,-1
            if (typ(k) .ne. ica) goto 300
            if (k .eq. j) goto 300
            do l=i-1,1,-1
              if (typ(l) .ne. ica) goto 400
              if (l .eq. j) goto 400
              if (l .eq. k) goto 400
              x1 = tangle (i,j,k,l,xyz)
              x2 = tangle (i,j,k,l,xyz2)
              if ( (x1-x2) .lt. -180.0) then
                x2 = x2 - 360.0
              else if ( (x1-x2) .gt. 180.0) then
                x2 = x2 + 360.0
              end if
              r = abs(x1-x2)
              if (r .lt. rmin) then
                rmin = r
                imin = j
                jmin = k
                kmin = l
                if (rmin .lt. tarcut) goto 510
              end if
  400         continue
            end do
  300       continue
          end do
  500     continue
        end do
  510   continue
        if (imin .gt. 0 .and. rmin .lt. tarcut) then
          write (*,6000) 2,rmin,range(i)
c
c          print *,i,typ(i),imin,typ(imin),jmin,kmin
c
          nok = nok + 1
          ref (1,i) = imin
          ref (2,i) = jmin
          ref (3,i) = kmin
          range (i) = rmin
        else
          nno = nno + 1
        end if
  100   continue
      end do
c
 6000 format (' Fix type ',i1,' : ',f6.1,' degrees instead of ',f6.1)
c
      return
      end
c
c
c
      subroutine plotqd (imol,irange,jmol,jrange,filnam,ierr,
     +                   bufsiz,bufmat)
c
      include 'lsqman.incl'
c
      integer bufsiz
      real bufmat(bufsiz)
c
      integer maxopt
      parameter (maxopt = 25)
c
      real xave,xsdv,xtot,xrms,xhav,xmin,xmax,q1,q2,distce,qdmax
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j,leng1
      integer idum,jdum,kdum,ldum,nuse,iptr,jptr,ires,jres
      integer iii,ndone
c
      logical xinter
c
      character irange*(*),jrange*(*),ichn*1,jchn*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20
      character filnam*(*)
      character line*256
c
code ...
c
      ierr = 0
c
      write (*,6000) name(imol)(1:leng1(name(imol))),
     +    irange(1:leng1(irange)),
     +    name(jmol)(1:leng1(name(jmol))),
     +    jrange(1:leng1(jrange)),atypes(1)
c
 6000 format (' Difference-distance matrix plot'/
     +        ' For ',a,1x,a/
     +        ' And ',a,1x,a/
     +        ' Atom type |',a4,'|')
c
      if (atypes(1) .eq. 'ALL ' .or. atypes(1) .eq. 'NONH' .or.
     +    atypes(1) .eq. 'SIDE' .or. atypes(1) .eq. 'TRAC') then
        call errcon ('Invalid atom type !')
        ierr = -1
        return
      end if
c
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range1')
        goto 9900
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range2')
        goto 9900
      end if
c
      if (nopt1 .ne. nopt2) then
        call errcon ('Different nr of zones')
        goto 9900
      end if
c
      nuse = 0
      iptr = 0
      jptr = 0
c
      do i=1,nopt1
c
        ichn = optpa1(i)(1:1)
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range1 : '//optpa1(i))
ccc          goto 9900
        end if
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
          call errcon ('Not a zone in range1 : '//optpa1(i))
          goto 9900
        end if
ccc      print *,optpa1(i)(2:j-1)
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
ccc      print *,optpa1(i)(j+1:)
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon ('Invalid zone in range1 : '//optpa1(i))
          goto 9900
        end if
c
        jchn = optpa2(i)(1:1)
        if (index('1234567890',jchn) .gt. 0) then
          jchn = ' '
ccc        else if (jchn .lt. 'A' .or. jchn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range2 : '//optpa2(i))
ccc          goto 9900
        end if
        j = index (optpa2(i),'-')
        if (j .le. 0) j = index (optpa2(i),':')
        if (j .gt. 0) then
          call errcon (
     +      'Zone ignored in range2 : '//optpa2(i))
          optpa2(i) = optpa2(i)(1:j-1)
        end if
ccc      print *,optpa2(i)(2:)
        if (jchn .eq. ' ') then
          call str2i (optpa2(i)(1:),kdum,ierr)
        else
          call str2i (optpa2(i)(2:),kdum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        ldum = kdum + (jdum - idum)
c
        do ires=idum,jdum
            jres = kdum + (ires - idum)
ccc      print *,imol,ichn,ires,atypes(j)
ccc      print *,jmol,jchn,jres
            call getptr (imol,ichn,ires,atypes(1),
     +                   jmol,jchn,jres,iptr,jptr,ierr)
ccc      print *,ierr,iptr,jptr
            if (ierr .ne. 0) goto 6996
c
            nuse = nuse + 1
            do iii=1,3
              buffi (iii,nuse) = atmxyz (iii,iptr,imol)
              buffj (iii,nuse) = atmxyz (iii,jptr,jmol)
            end do
c
 6996       continue
c
        end do
c
      end do
c
      call jvalut (' Nr of atoms found :',1,nuse)
      if (nuse .lt. 3) then
        call errcon ('Fewer than 3 atoms; cannot plot')
        goto 9900
      end if
c
c ... store matrix in BUFMAT -> check if not too big
c
      if ( (nuse*nuse) .gt. bufsiz ) then
        call errcon ('2D Plot matrix too big for buffer !')
        call jvalut (' Requested points:',1,(nuse*nuse))
        call jvalut (' Max buffer size :',1,bufsiz)
        goto 9900
      end if
c
      ndone = 0
      do i=1,nuse
        do j=1,nuse
          q1 = distce (buffi(1,i),buffi(1,j))
          q2 = distce (buffj(1,i),buffj(1,j))
          ndone = ndone + 1
          bufmat (ndone) = q1 - q2
        end do
      end do
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Atom number'
      write (iunit,6200) 'YLABEL','Atom number'
      write (iunit,6200) 'REMARK','Difference-distance plot'
      write (iunit,6200) 'REMARK',
     +    'Values are for mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
      write (iunit,6200) 'REMARK',
     +    '       and for mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
c
      call xstats (bufmat,nuse*nuse,xave,xsdv,xmin,xmax,xtot)
      call xstat2 (bufmat,nuse*nuse,xrms,xhav)
      write (iunit,6215) 'REMARK','Ave, sdv, min, max=',
     +  xave,xsdv,xmin,xmax
      write (iunit,6215) 'REMARK','RMS, harmonic ave=',xrms,xhav
c
      call jvalut (' Nr of distances in matrix :',1,(nuse*nuse))
      call fvalut (' Average diff-dist (A):',1,xave)
      call fvalut (' St.dev. diff-dist (A):',1,xsdv)
      call fvalut (' Minimum diff-dist (A):',1,xmin)
      call fvalut (' Maximum diff-dist (A):',1,xmax)
      call fvalut (' RMS     diff-dist (A):',1,xrms)
      call fvalut (' Harm.av diff-dist (A):',1,xhav)
c
      qdmax = max (abs(xmin),abs(xmax))
      call fvalut (' Max absolute DDM element (A) :',1,qdmax)
c
      write (iunit,6200) 'NLEVEL','8'
      write (iunit,6200) 'LEVELS'
      write (iunit,'(a)') '-5.0 -3.0 -2.0 -1.0 1.0 2.0 3.0 5.0'
      write (iunit,6200) 'COLOUR'
      write (iunit,'(a)') '1 5 2 4 4 2 5 1'
      write (iunit,6210) 'XPOINT',nuse
      write (iunit,6210) 'YPOINT',nuse
      write (iunit,6213) 'XLIMIT',1.0,float(nuse)
      write (iunit,6213) 'YLIMIT',1.0,float(nuse)
      write (iunit,6200) 'ZVALUE','(13f6.1)'
      write (iunit,'(13f6.1)') (bufmat(i),i=1,nuse*nuse)
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine perrot (iomat,ampli,lprint)
c
c ... generate random rotation matrix
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real iomat(3,3),rotmat(3,3),dx,dy,dz,det3,sum,detr
      real rx(3,3),ry(3,3),rz(3,3),dum(3,3),dd,ampli
c
      integer i,j
c
      logical lprint
c
code ...
c
 6536 continue
c
      do i=1,3
        do j=1,3
          rotmat (i,j) = 0.0
          rx (i,j) = 0.0
          ry (i,j) = 0.0
          rz (i,j) = 0.0
          dum (i,j) = 0.0
        end do
        rotmat (i,i) = 1.0
        rx (i,i) = 1.0
        ry (i,i) = 1.0
        rz (i,i) = 1.0
        dum (i,i) = 1.0
      end do
c
      dd = max(abs(ampli),0.001)
      call gkrand (dx,-dd,dd,0)
      call gkrand (dy,-dd,dd,0)
      call gkrand (dz,-dd,dd,0)
c
      write (*,'(a,3f8.2)')
     +  ' Rotations around X,Y,Z (deg) :',dx,dy,dz
c
      dx = degtor*dx
      dy = degtor*dy
      dz = degtor*dz
c
      rx (2,2) = cos(dx)
      rx (3,3) = rx (2,2)
      rx (3,2) = sin(dx)
      rx (2,3) = -rx (3,2)
      dx = det3 (rx)
      if (lprint) write (*,'(a12,3(3f13.7,:,/,12x))')
     +   ' X Matrix : ',((rx(i,j),j=1,3),i=1,3)
      if (lprint) call rvalut (
     +  ' Determinant of X rotation matrix :',1,dx)
c
      ry (1,1) = cos(dy)
      ry (3,3) = ry (1,1)
      ry (3,1) = sin(dy)
      ry (1,3) = -ry (3,1)
      dy = det3 (ry)
      if (lprint) write (*,'(a12,3(3f13.7,:,/,12x))')
     +   ' Y Matrix : ',((ry(i,j),j=1,3),i=1,3)
      if (lprint) call rvalut (
     +  ' Determinant of Y rotation matrix :',1,dy)
c
      rz (2,2) = cos(dz)
      rz (1,1) = rz (2,2)
      rz (1,2) = sin(dz)
      rz (2,1) = -rz (1,2)
      dz = det3 (rz)
      if (lprint) write (*,'(a12,3(3f13.7,:,/,12x))')
     +   ' Z Matrix : ',((rz(i,j),j=1,3),i=1,3)
      if (lprint) call rvalut (
     +  ' Determinant of Z rotation matrix :',1,dz)
c
      call mulmat (rx,ry,dum)
      call mulmat (dum,rz,rotmat)
c
      detr = det3 (rotmat)
      if (lprint) write (*,'(a12,3(3f13.7,:,/,12x))')
     +   '   Matrix : ',((rotmat(i,j),j=1,3),i=1,3)
      if (lprint) call rvalut (
     +  ' Determinant of random rotation matrix :',1,detr)
c
      if (abs(detr-1.0) .gt. 1.0e-5) goto 6536
      do i=1,3
        sum = rotmat(i,1)**2 + rotmat(i,2)**2 + rotmat(i,3)**2
        if (abs(sum - 1.0) .gt. 1.0e-5 ) goto 6536
        sum = rotmat(1,i)**2 + rotmat(2,i)**2 + rotmat(3,i)**2
        if (abs(sum - 1.0) .gt. 1.0e-5 ) goto 6536
      end do
c
c ... apply
c
      call mulmat (iomat,rotmat,dum)
      do i=1,3
        do j=1,3
          iomat (i,j) = dum(i,j)
        end do
      end do
c
      return
      end
c
c
c
      subroutine copyrt (rt1,rt2)
c
      implicit none
c
      real rt1(12),rt2(12)
c
      integer i
c
code ...
c
      do i=1,12
        rt1(i) = rt2(i)
      end do
c
      return
      end
c
c
c
      subroutine convol (filnam,imol,ic,jmol,jc,n1,n2,ncolen,
     +                   x1,x2,xmat,maxsiz)
c
      include 'lsqman.incl'
c
      integer imol,ic,jmol,jc,n1,n2,maxsiz,ierr,leng1,ncolen
      integer i,j,np,nd,nlen
c
      real x1(3,n1),x2(3,n2),xmat(maxsiz),rt(12)
      real xrmsd,ave,sdv,xmin,xmax,xtot
c
      logical xinter
c
      character filnam*(*),line*128
c
code ...
c
      nlen = 2*ncolen + 1
      np = (n1-nlen+1)*(n2-nlen+1)
      if (np .gt. maxsiz) then
        call errcon ('2D Plot matrix too big for buffer !')
        call jvalut (' Requested points:',1,np)
        call jvalut (' Max buffer size :',1,maxsiz)
        return
      end if
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
      nd = 0
      do j=ncolen+1,n2-ncolen
        do i=ncolen+1,n1-ncolen
          nd = nd + 1
          call lsqgjk (x1(1,i-ncolen),x2(1,j-ncolen),nlen,
     +                 xrmsd,rt,ierr)
          if (ierr .ne. 0) xrmsd = 99.99
          xmat (nd) = xrmsd
        end do
      end do
c
      call jvalut (' Nr of RMSDs   :',1,nd)
      call xstats (xmat(1),nd,ave,sdv,xmin,xmax,xtot)
      call fvalut (' Average RMSD  :',1,ave)
      call fvalut (' St. deviation :',1,sdv)
      call fvalut (' Minimum RMSD  :',1,xmin)
      call fvalut (' Maximum RMSD  :',1,xmax)
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Residue in mol 1'
      write (iunit,6200) 'YLABEL','Residue in mol 2'
      write (iunit,6200) 'REMARK','Convolution plot'
      write (iunit,6217) 'REMARK',
     +  'Fragment length (one-sided, total) ',ncolen,nlen
      write (iunit,6215) 'REMARK','Average RMSD (A) ',ave
      write (iunit,6215) 'REMARK','St. dev. (A) ',sdv
      write (iunit,6200) 'REMARK',
     +    'Values are for mol 1 = ',name(imol)
     +    (1:leng1(name(imol))),' chain ',chname(ic,imol),
     +    ' = file ',file(imol)(1:leng1(file(imol)))
      write (iunit,6200) 'REMARK',
     +    '       and for mol 2 = ',name(jmol)
     +    (1:leng1(name(jmol))),' chain ',chname(jc,jmol),
     +    ' = file ',file(jmol)(1:leng1(file(jmol)))
c
      write (iunit,6200) 'NLEVEL','6'
      write (iunit,6200) 'LEVELS'
c
      if (ncolen .le. 7) then
        write (iunit,'(a)') '0.01 0.2 0.4 0.6 0.8 1.0'
      else if (ncolen .le. 12) then
        write (iunit,'(a)') '0.01 0.3 0.6 0.9 1.2 1.5'
      else
        write (iunit,'(a)') '0.01 0.4 0.8 1.2 1.6 2.0'
      end if
c
      write (iunit,6200) 'COLOUR'
      write (iunit,'(a)') '1 1 5 2 6 4'
      write (iunit,6210) 'XPOINT',n1-nlen+1
      write (iunit,6210) 'YPOINT',n2-nlen+1
      write (iunit,6213) 'XLIMIT',float(1+ncolen),float(n1-ncolen)
      write (iunit,6213) 'YLIMIT',float(1+ncolen),float(n2-ncolen)
      write (iunit,6200) 'ZVALUE','*'
      write (iunit,'(10f7.2)') (xmat(i),i=1,nd)
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      return
      end
c
c
c
      subroutine delchi (imol,irange,jmol,jrange,
     +                   dum1,ierr)
c
      include 'lsqman.incl'
c
      integer maxopt
      parameter (maxopt = 25)
c
      real tangle,dum1,outcut
      real c1i,c1j,d1,c2i,c2j,d2,c3i,c3j,d3,c4i,c4j,d4
c
      integer imol,jmol,length,ierr,nopt1,nopt2,i,j,leng1
      integer idum,jdum,kdum,ldum,nuse,iptr,jptr,ires,jres
      integer inptr,jnptr,jeptr,izptr,jzptr
      integer nbad,ibptr,jbptr,igptr,jgptr,idptr,jdptr,ieptr
c
      logicallprint,lchi1,lchi2,lchi3,lchi4
c
      character irange*(*),jrange*(*),ichn*1,jchn*1
      character optpa1(maxopt)*20,optpa2(maxopt)*20
      character a1*10,a2*10,a3*10,a4*10,chr*1
c
code ...
c
      ierr = 0
c
      call prompt (' Delta-CHI table')
      outcut = 10.0
      if (dum1 .gt. 0.0 .and. dum1 .le. 180.0) outcut = dum1
c
      write (*,6000) name(imol)(1:leng1(name(imol))),
     +    irange(1:leng1(irange)),
     +    name(jmol)(1:leng1(name(jmol))),
     +    jrange(1:leng1(jrange)),outcut
c
 6000 format (' Compare ',a,1x,a/
     +        ' And     ',a,1x,a/
     +        ' Outlier list cut-off  : ',f8.2)
c
      do i=1,length(irange)
        if (irange(i:i) .eq. '"') irange(i:i) = ' '
      end do
c
      call extrop (irange,nopt1,maxopt,optpa1,ierr)
      if (nopt1 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range1')
        goto 9900
      end if
c
      do i=1,length(jrange)
        if (jrange(i:i) .eq. '"') jrange(i:i) = ' '
      end do
c
      call extrop (jrange,nopt2,maxopt,optpa2,ierr)
      if (nopt2 .lt. 1 .or. ierr .ne. 0) then
        call errcon ('In range2')
        goto 9900
      end if
c
      if (nopt1 .ne. nopt2) then
        call errcon ('Different nr of zones')
        goto 9900
      end if
c
      nuse = 0
      nbad = 0
      iptr = 0
      jptr = 0
c
      do i=1,nopt1
c
        ichn = optpa1(i)(1:1)
        if (index('1234567890',ichn) .gt. 0) then
          ichn = ' '
ccc        else if (ichn .lt. 'A' .or. ichn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range1 : '//optpa1(i))
ccc          goto 9900
        end if
        j = index (optpa1(i),'-')
        if (j .le. 1) j = index (optpa1(i),':')
        if (j .le. 1) then
          call errcon ('Not a zone in range1 : '//optpa1(i))
          goto 9900
        end if
ccc      print *,optpa1(i)(2:j-1)
        if (ichn .eq. ' ') then
          call str2i (optpa1(i)(1:j-1),idum,ierr)
        else
          call str2i (optpa1(i)(2:j-1),idum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
ccc      print *,optpa1(i)(j+1:)
        if (optpa1(i)(j+1:j+1) .eq. ichn)
     +    optpa1(i)(j+1:j+1) = ' '
        call str2i (optpa1(i)(j+1:),jdum,ierr)
        if (ierr .ne. 0) goto 9900
        if (idum .le. 0 .or. jdum .le. 0 .or. idum .gt. jdum) then
          call errcon ('Invalid zone in range1 : '//optpa1(i))
          goto 9900
        end if
c
        jchn = optpa2(i)(1:1)
        if (index('1234567890',jchn) .gt. 0) then
          jchn = ' '
ccc        else if (jchn .lt. 'A' .or. jchn .gt. 'Z') then
ccc          call errcon ('Invalid chain id in range2 : '//optpa2(i))
ccc          goto 9900
        end if
        j = index (optpa2(i),'-')
        if (j .le. 0) j = index (optpa2(i),':')
        if (j .gt. 0) then
          call errcon (
     +      'Zone ignored in range2 : '//optpa2(i))
          optpa2(i) = optpa2(i)(1:j-1)
        end if
ccc      print *,optpa2(i)(2:)
        if (jchn .eq. ' ') then
          call str2i (optpa2(i)(1:),kdum,ierr)
        else
          call str2i (optpa2(i)(2:),kdum,ierr)
        end if
        if (ierr .ne. 0) goto 9900
        ldum = kdum + (jdum - idum)
c
        do ires=idum,jdum
c
          lprint = .false.
          lchi1 = .false.
          lchi2 = .false.
          lchi3 = .false.
          lchi4 = .false.
c
          jres = kdum + (ires - idum)
          call getptr (imol,ichn,ires,' CA ',
     +                 jmol,jchn,jres,iptr,jptr,ierr)
          if (ierr .ne. 0) goto 6196
c
          call getptr (imol,ichn,ires,' N  ',
     +                 jmol,jchn,jres,inptr,jnptr,ierr)
          if (ierr .ne. 0) goto 6196
c
          call getptr (imol,ichn,ires,' CB ',
     +                 jmol,jchn,jres,ibptr,jbptr,ierr)
          if (ierr .ne. 0) goto 6196
c
c ... CHI1 defined ?
c
          call getptr (imol,ichn,ires,' CG ',
     +                 jmol,jchn,jres,igptr,jgptr,ierr)
          if (ierr .eq. 0) goto 6110
c
          call getptr (imol,ichn,ires,' OG ',
     +                 jmol,jchn,jres,igptr,jgptr,ierr)
          if (ierr .eq. 0) goto 6110
c
          call getptr (imol,ichn,ires,' SG ',
     +                 jmol,jchn,jres,igptr,jgptr,ierr)
          if (ierr .eq. 0) goto 6110
c
          call getptr (imol,ichn,ires,' CG1',
     +                 jmol,jchn,jres,igptr,jgptr,ierr)
          if (ierr .eq. 0) goto 6110
c
          call getptr (imol,ichn,ires,' OG1',
     +                 jmol,jchn,jres,igptr,jgptr,ierr)
          if (ierr .eq. 0) goto 6110
c
          goto 6196
c
 6110     continue
          lchi1 = .true.
c
c ... CHI2 defined ?
c
          call getptr (imol,ichn,ires,' CD ',
     +                 jmol,jchn,jres,idptr,jdptr,ierr)
          if (ierr .eq. 0) goto 6120
c
          call getptr (imol,ichn,ires,' SD ',
     +                 jmol,jchn,jres,idptr,jdptr,ierr)
          if (ierr .eq. 0) goto 6120
c
          call getptr (imol,ichn,ires,' CD1',
     +                 jmol,jchn,jres,idptr,jdptr,ierr)
          if (ierr .eq. 0) goto 6120
c
          call getptr (imol,ichn,ires,' OD1',
     +                 jmol,jchn,jres,idptr,jdptr,ierr)
          if (ierr .eq. 0) goto 6120
c
          call getptr (imol,ichn,ires,' ND1',
     +                 jmol,jchn,jres,idptr,jdptr,ierr)
          if (ierr .eq. 0) goto 6120
c
          goto 6200
c
 6120     continue
          lchi2 = .true.
c
c ... CHI3 defined ?
c
          call getptr (imol,ichn,ires,' CE ',
     +                 jmol,jchn,jres,ieptr,jeptr,ierr)
          if (ierr .eq. 0) goto 6130
c
          call getptr (imol,ichn,ires,' CE1',
     +                 jmol,jchn,jres,ieptr,jeptr,ierr)
          if (ierr .eq. 0) goto 6130
c
          call getptr (imol,ichn,ires,' NE ',
     +                 jmol,jchn,jres,ieptr,jeptr,ierr)
          if (ierr .eq. 0) goto 6130
c
          call getptr (imol,ichn,ires,' NE1',
     +                 jmol,jchn,jres,ieptr,jeptr,ierr)
          if (ierr .eq. 0) goto 6130
c
          call getptr (imol,ichn,ires,' OE1',
     +                 jmol,jchn,jres,ieptr,jeptr,ierr)
          if (ierr .eq. 0) goto 6130
c
          goto 6200
c
 6130     continue
          lchi3 = .true.
c
c ... CHI4 defined ?
c
          call getptr (imol,ichn,ires,' CZ ',
     +                 jmol,jchn,jres,izptr,jzptr,ierr)
          if (ierr .eq. 0) goto 6140
c
          call getptr (imol,ichn,ires,' NZ ',
     +                 jmol,jchn,jres,izptr,jzptr,ierr)
          if (ierr .eq. 0) goto 6140
c
          goto 6200
c
 6140     continue
          lchi4 = .true.
c
 6200     continue
          a1 = ' '
          a2 = ' '
          a3 = ' '
          a4 = ' '
c
          c1i = tangle (inptr,iptr,ibptr,igptr,atmxyz(1,1,imol))
          c1j = tangle (jnptr,jptr,jbptr,jgptr,atmxyz(1,1,jmol))
          call fixadf (c1i,c1j,d1)
          chr = ' '
          if (abs(d1) .ge. outcut) then
            lprint = .true.
            chr = '*'
          end if
          write (a1,'(f8.2,1x,a1)') d1,chr
          if (.not. lchi2) goto 6300
c
          c2i = tangle (iptr,ibptr,igptr,idptr,atmxyz(1,1,imol))
          c2j = tangle (jptr,jbptr,jgptr,jdptr,atmxyz(1,1,jmol))
          call fixadf (c2i,c2j,d2)
          chr = ' '
          if (abs(d2) .ge. outcut) then
            lprint = .true.
            chr = '*'
          end if
          write (a2,'(f8.2,1x,a1)') d2,chr
          if (.not. lchi3) goto 6300
c
          c3i = tangle (ibptr,igptr,idptr,ieptr,atmxyz(1,1,imol))
          c3j = tangle (jbptr,jgptr,jdptr,jeptr,atmxyz(1,1,jmol))
          call fixadf (c3i,c3j,d3)
          chr = ' '
          if (abs(d3) .ge. outcut) then
            lprint = .true.
            chr = '*'
          end if
          write (a3,'(f8.2,1x,a1)') d3,chr
          if (.not. lchi4) goto 6300
c
          c4i = tangle (igptr,idptr,ieptr,izptr,atmxyz(1,1,imol))
          c4j = tangle (jgptr,jdptr,jeptr,jzptr,atmxyz(1,1,jmol))
          call fixadf (c4i,c4j,d4)
          chr = ' '
          if (abs(d4) .ge. outcut) then
            lprint = .true.
            chr = '*'
          end if
          write (a4,'(f8.2,1x,a1)') d4,chr
c
 6300     continue
          nuse = nuse + 1
          if (.not. lprint) goto 6196
          nbad = nbad + 1
          if (nbad .eq. 1) write (*,6204) 'Residue 1',
     +      'Residue 2','CHI-1','CHI-2','CHI-3','CHI-4'
c
          write (*,6206)
     +        resnam(iptr,imol),
     +        achain(iptr,imol),
     +        iresid(iptr,imol),
     +        resnam(jptr,jmol),
     +        achain(jptr,jmol),
     +        iresid(jptr,jmol),
     +        a1,a2,a3,a4
c
 6204 format (/1x,a9,3x,a9,3x,4(a10,2x))
 6206 format (1x,a3,1x,a1,i4,' - ',
     +        a3,1x,a1,i4,' = ',4(a10,2x))
c
 6196     continue
c
        end do
c
      end do
c
      call jvalut (' Nr of residues :',1,nuse)
      call jvalut (' Nr of outliers :',1,nbad)
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      ierr = -1
c
      return
      end
c
c
c
      subroutine nwlocs (imol,ic,jmol,jc,n1,n2,ncolen,
     +                   x1,x2,xmat,gappen,fmat,pmat,
     +                   seq1,seq2,rtij)
c
      include 'lsqman.incl'
c
      integer maxg11
      parameter (maxg11=2*maxres)
c
      integer imol,ic,jmol,jc,n1,n2,ierr,ncolen,nma
      integer i,j,nlen,n1n,inow,jnow,n1d,nsi,k,k1,k2
      integer pmat(0:n1,0:n2)
c
      real x1(3,n1),x2(3,n2),xmat(n1,n2),fmat(0:n1,0:n2),rt(12)
      real rmsent(0:maxg11),rtij(12)
      real gappen,r1,r2,r3,xdum,xrmsd,qqq
c
      character*1 seq1(n1),seq2(n2),l1(maxg11),l2(maxg11),l3(maxg11)
c
code ...
c
      nlen = 2*ncolen + 1
c
      write (*,*)
      call prompt (' Calculating local RMSD matrix ...')
      do i=1,n1
        do j=1,n2
          xmat (i,j) = 0.0
        end do
      end do
c
      do j=ncolen+1,n2-ncolen
        do i=ncolen+1,n1-ncolen
          call lsqgjk (x1(1,i-ncolen),x2(1,j-ncolen),nlen,
     +                 xrmsd,rt,ierr)
          if (ierr .ne. 0) xrmsd = 99.99
          xmat (i,j) = xrmsd
        end do
      end do
c
      write (*,*)
      call prompt (' Executing Needleman-Wunsch ...')
      do i=0,n1
        do j=0,n2
          fmat (i,j) = 0.0
          pmat (i,j) = -1
        end do
      end do
c
      do i=0,n1
        fmat (i,0) = -float(i)*gappen
        pmat (i,0) = 0
      end do
      do j=0,n2
        fmat (0,j) = -float(j)*gappen
        pmat (0,j) = 0
      end do
c
      do j=1,n2
        do i=1,n1
          r1 = fmat(i-1,j-1) - xmat(i,j)
          r2 = fmat(i-1,j) - gappen
          r3 = fmat(i,j-1) - gappen
          if (r1 .ge. r2 .and. r1 .ge. r3) then
            fmat (i,j) = r1
            pmat (i,j) = 3
          else if (r2 .ge. r3) then
            fmat (i,j) = r2
            pmat (i,j) = 2
          else
            fmat (i,j) = r3
            pmat (i,j) = 1
          end if
        end do
      end do
c
      n1n = maxg11
      inow = n1
      jnow = n2
      n1d = 0
      nsi = 0
      k1 = n1
      k2 = n2
c
  100 continue
      if (pmat(inow,jnow) .eq. 3) then
        l1(n1n) = seq1(inow)
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
        rmsent(n1n) = xmat(inow,jnow)
c
        x1(1,k1) = x1(1,inow)
        x1(2,k1) = x1(2,inow)
        x1(3,k1) = x1(3,inow)
        k1 = k1 - 1
        x2(1,k2) = x2(1,jnow)
        x2(2,k2) = x2(2,jnow)
        x2(3,k2) = x2(3,jnow)
        k2 = k2 - 1
c
        if (l1(n1n) .eq. l2(n1n) .and.
     +      l1(n1n) .ne. '?') then
          n1d = n1d + 1
          nsi = nsi + 1
          l3(n1n) = '|'
c        else if (cmpmat(iptr(inow),jptr(jnow)) .gt. 0) then
c          nsi = nsi + 1
c          l3(n1n) = '+'
        end if
        inow = inow - 1
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 2) then
        l1(n1n) = seq1(inow)
        l2(n1n) = '-'
        l3(n1n) = ' '
        rmsent(n1n) = 999.999
        inow = inow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 1) then
        l1(n1n) = '-'
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
        rmsent(n1n) = 999.999
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      end if
c
c      print *,' INOW, JNOW ',inow,jnow
c
      if (inow .eq. 0 .and. jnow .gt. 0) then
        do j=jnow,1,-1
          l1(n1n) = '-'
          l2(n1n) = seq2(j)
          l3(n1n) = ' '
          rmsent(n1n) = 999.999
          n1n = n1n - 1
        end do
      else if (inow .gt. 0 .and. jnow .eq. 0) then
        do i=inow,1,-1
          l1(n1n) = seq1(i)
          l2(n1n) = '-'
          l3(n1n) = ' '
          rmsent(n1n) = 999.999
          n1n = n1n - 1
        end do
      end if
c
      n1n = n1n + 1
c
      write (*,*)
      k = 0
      do i=n1n,maxg11
        k = k + 1
        write (*,6010) k,l1(i),l3(i),l2(i),rmsent(i)
      end do
c
 6010 format (1x,i6,1x,3(a1,1x),' RMSD = ',f8.2,' A')
 6020 format (1x,a10,1x,60a1)
c
      write (*,*)
      do i=n1n,maxg11,60
        j = min (i+59,maxg11)
        write (*,*)
        write (*,6020) 'Sequence 1 ',(l1(k),k=i,j)
        write (*,6020) '   |=ID    ',(l3(k),k=i,j)
        write (*,6020) 'Sequence 2 ',(l2(k),k=i,j)
      end do
c
      write (*,*)
      call fvalut (' Gap penalty         :',1,gappen)
      call rvalut (' Raw alignment score :',1,fmat(n1,n2))
      call ivalut (' Length sequence 1   :',1,n1)
      call ivalut (' Length sequence 2   :',1,n2)
      i = maxg11 - n1n + 1
      call ivalut (' Alignment length    :',1,i)
      call ivalut (' Nr of identities    :',1,n1d)
      xdum = 100.0 * float(n1d) / float(min(n1,n2))
      call fvalut (' Perc identities     :',1,xdum)
c      call ivalut (' Nr of similarities  :',1,nsi)
c      xdum = 100.0 * float(nsi) / float(min(n1,n2))
c      call fvalut (' Perc similarities   :',1,xdum)
c
      k1 = k1 + 1
      k2 = k2 + 1
      i = (n1-k1+1)
      call ivalut (' Nr of matched res   :',1,i)
c
      if (i .lt. 3) then
        call errcon ('Not enough residues matched (<3) !')
        return
      end if
c
      call lsqgjk (x1(1,k1),x2(1,k2),i,
     +             xrmsd,rtij,ierr)
      call fvalut (' RMSD (A) for those  :',1,xrmsd)
      call fvalut (' Operator stored     :',12,rtij)
c
c ... update statistics
c
      nma = i
      simind (imol,jmol) = xrmsd * float(min(n1,n2)) / float(nma)
      matchi (imol,jmol) = (1.0+float(nma)) /
     +       ( (1.0+rmswgt*xrmsd) * (1.0+float(min (n1,n2))) )
      call mcrho (nma,x1(1,k1),x2(1,k2),xrmsd,
     +            cripp(imol,jmol),rrmsd(imol,jmol),
     +            normsd(imol,jmol),xdum)
c
      nmatch (imol,jmol) = nma
      rmsd (imol,jmol) = xrmsd
      corb (imol,jmol) = 999.99999
      rmsb (imol,jmol) = 999.99999
      rmsdna (imol,jmol) = xrmsd / float(nma)
      qqq = 100.0/float(nma)
      sas1 (imol,jmol) = xrmsd * qqq
      sas2 (imol,jmol) = sas1 (imol,jmol) * qqq
      sas3 (imol,jmol) = sas2 (imol,jmol) * qqq
      sas4 (imol,jmol) = sas3 (imol,jmol) * qqq
c
      write (*,6030) nma,xrmsd,simind(imol,jmol),matchi(imol,jmol),
     +  cripp(imol,jmol),xdum,rrmsd(imol,jmol),normsd(imol,jmol),
     +  rmsdna(imol,jmol),sas1(imol,jmol),sas2(imol,jmol),
     +  sas3(imol,jmol),sas4(imol,jmol)
c
 6030 format (/
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' SI = RMS * Nmin / Nmatch             = ',f12.5/
     +  ' MI = (1+Nmatch)/{(1+W*RMS)*(1+Nmin)} = ',f12.5/
     +  ' CR = Maiorov-Crippen RHO (0-2)       = ',f12.5/
     +  ' Estimated RMSD for 2 random proteins = ',f10.3,' A'/
     +  ' RR = Relative RMSD                   = ',f12.5/
     +  ' NR = Normalised RMSD (100)           = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A'/
     +  ' SAS(1) = cRMS * (100/Nmatch)         = ',f10.3,' A'/
     +  ' SAS(2) = cRMS * (100/Nmatch)^2       = ',f10.3,' A'/
     +  ' SAS(3) = cRMS * (100/Nmatch)^3       = ',f10.3,' A'/
     +  ' SAS(4) = cRMS * (100/Nmatch)^4       = ',f10.3,' A')
c
      return
      end
c
c
c
      subroutine nwglob (imol,ic,jmol,jc,n1,n2,nwcut,
     +                   x1,x2,xmat,fmat,pmat,
     +                   seq1,seq2,rtij,filnam,iremin,ireshi,
     +                   chnami,chnamj,iresi,iresj)
c
      include 'lsqman.incl'
c
      integer maxg11
      parameter (maxg11=2*maxres)
c
      integer imol,ic,jmol,jc,n1,n2,nali,nmres,alilen,length
      integer i,j,n1n,inow,jnow,n1d,nsi,k,k1,k2,k3,k4,ierr,ngap
      integer iremin,ireshi,nres,i1,i2,m,ltarg
      integer pmat(0:n1,0:n2)
      integer iresi(n1),iresj(n2),lresi(maxg11),lresj(maxg11)
c
      real*8 dz,dprob
c
      real x1(3,n1),x2(3,n2),xmat(n1,n2),fmat(0:n1,0:n2)
      real rmsent(0:maxg11),rtij(12),rt(12)
      real nwcut,gappen,r1,r2,r3,xdum,xrmsd,z,sstr,ssum,xdum2
      real tmscor,d0
c
      logical lgap,logfil,xinter
c
      character*1 seq1(n1),seq2(n2),l1(maxg11),l2(maxg11),l3(maxg11)
      character res1*(maxg11),res2*(maxg11)
      character chnami*1,chnamj*1,number*10,filnam*(*),nami*7,namj*7
c
code ...
c
      logfil = (length(filnam) .gt. 0 .and.
     +          filnam .ne. '-')
      if (logfil) then
        close (iunit)
        call xopxua (iunit,filnam,xinter(),ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening log file')
          logfil = .false.
        end if
      end if
c
c ... gap penalty = half of square of cut-off distance
c                   this means that no two residues will be
c                   aligned if they are more than the cut-off
c                   distance apart, since it then is cheaper
c                   to introduce two gaps instead
c
      gappen = 0.5*nwcut*nwcut
c
c ... Ltarget = target alignment length (Zhang & Skolnick)
c               needed to calculate the TM-score; to make the
c               score symmetric we define it as the shortest
c               of the two sequences
c
      ltarg = min (n1,n2)
      d0 = ( 1.24 * (float(ltarg)-15.0)**0.333 ) - 1.8
c
c ... first apply operator to coordinates of mol 2
c
      write (*,*)
      call fvalut (' Applying current operator to mol 2 :',12,rtij)
c      call fvalut (' BEFORE :',3,x2(1,1))
      call vecrtv (x2,x2,n2,rtij(1),rtij(10))
c      call fvalut (' AFTER  :',3,x2(1,1))
c
      write (*,*)
      call prompt (' Calculating superposition-distance matrix ...')
c
      do i=1,n1
        do j=1,n2
          xmat (i,j) = 0.0
        end do
      end do
c
c      xdum = nwcut*nwcut
      do j=1,n2
        do i=1,n1
          xrmsd = (x1(1,i)-x2(1,j))**2 +
     +            (x1(2,i)-x2(2,j))**2 +
     +            (x1(3,i)-x2(3,j))**2
          xmat (i,j) = - xrmsd
c          if (xrmsd .le. xdum) then
c            xmat (i,j) = xdum - xrmsd
c            print *,'I/J ',i,j,xrmsd,xdum
c          end if
        end do
      end do
c
      write (*,*)
      call prompt (' Executing Needleman-Wunsch ...')
      do i=0,n1
        do j=0,n2
          fmat (i,j) = 0.0
          pmat (i,j) = -1
        end do
      end do
c
      do i=0,n1
        fmat (i,0) = -float(i)*gappen
        pmat (i,0) = 0
      end do
      do j=0,n2
        fmat (0,j) = -float(j)*gappen
        pmat (0,j) = 0
      end do
c
      do j=1,n2
        do i=1,n1
          r1 = fmat(i-1,j-1) + xmat(i,j)
          r2 = fmat(i-1,j) - gappen
          r3 = fmat(i,j-1) - gappen
          if (r1 .ge. r2 .and. r1 .ge. r3) then
            fmat (i,j) = r1
            pmat (i,j) = 3
          else if (r2 .ge. r3) then
            fmat (i,j) = r2
            pmat (i,j) = 2
          else
            fmat (i,j) = r3
            pmat (i,j) = 1
          end if
        end do
      end do
c
      n1n = maxg11
      inow = n1
      jnow = n2
      n1d = 0
      nsi = 0
      k1 = n1
      k2 = n2
      sstr = 0.0
      tmscor = 0.0
c
  100 continue
      if (pmat(inow,jnow) .eq. 3) then
c
        l1(n1n) = seq1(inow)
        l2(n1n) = seq2(jnow)
        l3(n1n) = '.'
        lresi(n1n) = iresi(inow)
        lresj(n1n) = iresj(jnow)
        rmsent(n1n) = sqrt (- xmat(inow,jnow))
c
        sstr = sstr + (1.0 / (1.0 + (-xmat(inow,jnow)/25.0)))
c
        tmscor = tmscor + (1.0 / (1.0 + (-xmat(inow,jnow)/d0**2)))
c
        x1(1,k1) = x1(1,inow)
        x1(2,k1) = x1(2,inow)
        x1(3,k1) = x1(3,inow)
        k1 = k1 - 1
        x2(1,k2) = x2(1,jnow)
        x2(2,k2) = x2(2,jnow)
        x2(3,k2) = x2(3,jnow)
        k2 = k2 - 1
c
        if (l1(n1n) .eq. l2(n1n) .and.
     +      l1(n1n) .ne. '?') then
          n1d = n1d + 1
          nsi = nsi + 1
          l3(n1n) = '|'
c        else if (cmpmat(iptr(inow),jptr(jnow)) .gt. 0) then
c          nsi = nsi + 1
c          l3(n1n) = '+'
        end if
        inow = inow - 1
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 2) then
        l1(n1n) = seq1(inow)
        l2(n1n) = '-'
        l3(n1n) = ' '
        lresi(n1n) = iresi(inow)
        lresj(n1n) = -1
        rmsent(n1n) = -999.999
        inow = inow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 1) then
        l1(n1n) = '-'
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
        lresi(n1n) = -1
        lresj(n1n) = iresj(jnow)
        rmsent(n1n) = -999.999
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      end if
c
c      print *,' INOW, JNOW ',inow,jnow
c
      if (inow .eq. 0 .and. jnow .gt. 0) then
        do j=jnow,1,-1
          l1(n1n) = '-'
          l2(n1n) = seq2(j)
          l3(n1n) = ' '
          lresi(n1n) = -1
          lresj(n1n) = iresj(j)
          rmsent(n1n) = -999.999
          n1n = n1n - 1
        end do
      else if (inow .gt. 0 .and. jnow .eq. 0) then
        do i=inow,1,-1
          l1(n1n) = seq1(i)
          l2(n1n) = '-'
          l3(n1n) = ' '
          lresi(n1n) = iresi(i)
          lresj(n1n) = -1
          rmsent(n1n) = -999.999
          n1n = n1n - 1
        end do
      end if
c
      n1n = n1n + 1
c
      write (*,*)
      k = 0
      nali = 0
      if (logfil) write (iunit,*)
      do i=n1n,maxg11
        k = k + 1
        if (rmsent(i) .ge. 0.0) then
          write (number,'(f8.2,1x,a1)') rmsent(i),'A'
          nali = nali + 1
          buffb3 (nali) = rmsent(i)
        else
          number = '       -'
        end if
        if (lresi(i) .gt. 0) then
          write (nami,'(a1,a1,i4,a1)') '[',chnami,lresi(i),']'
        else
          nami = '   -   '
        end if
        if (lresj(i) .gt. 0) then
          write (namj,'(a1,a1,i4,a1)') '[',chnamj,lresj(i),']'
        else
          namj = '   -   '
        end if
c
ccc        write (*,6010) k,l1(i),l3(i),l2(i),number
ccc        write (*,6015) k,nami,l1(i),l3(i),l2(i),namj,number
c
        if (logfil) then
          write (iunit,6015) k,nami,l1(i),l3(i),l2(i),
     +      namj,number
        else
          write (*,6015) k,nami,l1(i),l3(i),l2(i),namj,number
        end if
      end do
c
      if (logfil) call prompt (
     +   ' Detailed structural alignment suppressed - see log file')
c
ccc 6010 format (1x,i6,1x,3(a1,1x),' DIST = ',a10)
 6015 format (1x,i6,1x,a7,1x,3(a1,1x),a8,' ... DIST = ',a10)
 6020 format (1x,a10,1x,60a1)
c
      do i=n1n,maxg11,60
        j = min (i+59,maxg11)
        write (*,*)
        write (*,6020) 'Sequence 1 ',(l1(k),k=i,j)
        write (*,6020) '.=ALI |=ID ',(l3(k),k=i,j)
        write (*,6020) 'Sequence 2 ',(l2(k),k=i,j)
        if (logfil) then
          write (iunit,*)
          write (iunit,6020) 'Sequence 1 ',(l1(k),k=i,j)
          write (iunit,6020) '.=ALI |=ID ',(l3(k),k=i,j)
          write (iunit,6020) 'Sequence 2 ',(l2(k),k=i,j)
        end if
      end do
c
      write (*,*)
      call distat (nali,buffb3)
c
      alilen = maxg11 - n1n + 1
c
      write (*,6100) gappen,fmat(n1,n2),n1,n2,alilen,n1d
      if (logfil) write (iunit,6100) gappen,fmat(n1,n2),n1,n2,
     +  alilen,n1d
c
 6100 format (/' Gap penalty            : ',f12.3/
     +         ' Raw alignment score    : ',1p,e12.4,0p/
     +         ' L1 = Length sequence 1 : ',i12/
     +         ' L2 = Length sequence 2 : ',i12/
     +         ' Alignment length       : ',i12/
     +         ' NI = Nr of identities  : ',i12)
c
c      call ivalut (' Nr of similarities  :',1,nsi)
c      xdum = 100.0 * float(nsi) / float(min(n1,n2))
c      call fvalut (' Perc similarities   :',1,xdum)
c
 6110 format (' L3 = Nr of matched res : ',i12/
     +        ' RMSD for those (A)     : ',f12.3/
     +        ' (Note: RMSD calculated after superimposing matched'/
     +        '  residues; corresponding RT operator not stored)'/
     +        ' ID = NI/min(L1,L2) (%) : ',f12.2/
     +        ' ID = NI/L3 (%)         : ',f12.2)
   
 6120 format (' L3 = Nr of matched res : ',i12/
     +        ' ID = NI/min(L1,L2) (%) : ',f12.2/
     +        ' ID = NI/L3 (%)         : ',f12.2/
     +        ' Fewer than 3 residues - not enough for superposition'/
     +        ' Cannot calculate Levitt-Gerstein statistics !'/)
c
      k1 = k1 + 1
      k2 = k2 + 1
      nmres = (n1-k1+1)
      if (n1d .gt. 0) then
        xdum  = 100.0 * float(n1d) / float(min(n1,n2))
        xdum2 = 100.0 * float(n1d) / float(nmres)
      else
        xdum = 0.0
        xdum2 = 0.0
      end if
c
      if (nmres .ge. 3) then
        call lsqgjk (x1(1,k1),x2(1,k2),nmres,
     +               xrmsd,rt,ierr)
        write (*,6110) nmres,xrmsd,xdum,xdum2
        if (logfil) write (iunit,6110) nmres,xrmsd,xdum,xdum2
ccc        call fvalut (' Operator NOT stored :',12,rt)
      else
        write (*,6120) nmres,xdum,xdum2
        if (logfil) write (iunit,6120) nmres,xdum,xdum2
        goto 9999
      end if
c
c ... Levitt-Gerstein statistics from M Levitt & M Gerstein,
c     PNAS 95, 5913-5920 (1998)
c
c ... what is Ngap ? sum of # gaps in mol 1 & mol 2, excl. termini
c
      ngap = 0
      lgap = .false.
      do i=n1n,maxg11
        if (.not. lgap) then
          if (l1(i) .eq. '-' .or. l2(i) .eq. '-') then
            ngap = ngap + 1
          end if
        end if
        lgap = (l1(i) .eq. '-' .or. l2(i) .eq. '-')
      end do
      if (l1(n1n) .eq. '-' .or.
     +    l2(n1n) .eq. '-') ngap = ngap - 1
      if (l1(maxg11) .eq. '-' .or.
     +    l2(maxg11) .eq. '-') ngap = ngap - 1
      if (ngap .lt. 0) ngap = 0
c
      ssum = sstr
      sstr = 20.0 * (ssum - 0.5 * float(ngap))
      xdum = alog(float(nmres))
      if (nmres .lt. 120) then
        z = (sstr - (18.4*xdum*xdum -4.50*xdum +2.64)) /
     +      (21.4*xdum -37.5)
      else
        z = (sstr - (171.7*xdum - 419.2)) /
     +      (21.4*alog(120.0) -37.5)
      end if
c
ccc      prob = 1.0 - exp(-exp(-z))
c
c ... evaluate P-value in double precision
c
      dz = z
      dprob = 1.0D0 - dexp ( -dexp (-dz))
c
      write (*,6130) ngap,sstr,dz,dprob
      if (logfil) write (iunit,6130) ngap,sstr,dz,dprob
 6130 format (/' Levitt-Gerstein statistics:'/
     +         ' Nr of gaps       : ',i12/
     +         ' Similarity score : ',1p,e12.4/
     +         ' Z-score          : ',e12.4/
     +         ' P (z > Z)        : ',e12.4/
     +  ' P (z > Z) is the probability of matching any two'/
     +  ' random structures and finding a Z-score z which'/
     +  ' is greater than the Z-score Z of the current pair.')
c
c ... TM-score
c
      tmscor = tmscor / float(ltarg)
c
      write (*,6133) ltarg,d0,tmscor
      if (logfil) write (iunit,6133) ltarg,d0,tmscor
 6133 format (/' TM-score statistics:'/
     +         ' Ltarget          : ',i12/
     +         ' d0 (Ltarget) (A) : ',f12.3/
     +         ' TM-score         : ',f12.3)
c
      write (*,*)
c
c ... check for register shifts ?
c
      if (iremin .le. 0) goto 9999
c
      res1 = ' '
      res2 = ' '
      do i=n1n,maxg11
        j = i - n1n + 1
        res1(j:j) = l1(i)
        if (res1(j:j) .eq. '?') res1(j:j) = 'X'
        res2(j:j) = l2(i)
        if (res2(j:j) .eq. '?') res2(j:j) = 'X'
        lresi(j) = lresi(i)
        lresj(j) = lresj(i)
      end do
c
      nres = maxg11 - n1n + 1
      i = 0
c
c ... loop over residues in sequence alignment
c
 2128 continue
      i = i + 1
c
c ... end of sequence reached ?
c
      if (i .gt. nres) goto 2133
c
c ... skip if sequence contains unaligned residues
c
      if (index(res1(i:i+iremin-1),'-') .gt. 0) goto 2128
c
c ... try all shifts
c
      do j=i-ireshi,i+ireshi
        if (j .ge. 1 .and. (j+iremin-1) .le. nres .and.
     +      j .ne. i) then
c
c ... check if short sequences are identical (and no unaligned
c     residues in the other sequence)
c
          if (res1(i:i+iremin-1) .eq. res2(j:j+iremin-1) .and.
     +        index(res2(j:j+iremin-1),'-') .le. 0 .and.
     +        res1(i:i+iremin-1) .ne. res1(j:j+iremin-1)) then
c
ccc            write (*,*)
c
c ... find out how much of the sequences is out of register
c
c ... K2 = number of residues involved in shift
c
            k2 = 0
            do k1=0,nres
              if ( (i+k1) .gt. nres) goto 2138
              if ( (j+k1) .gt. nres) goto 2138
              if (res1(i+k1:i+k1) .eq. '-') goto 2138
              if (res2(j+k1:j+k1) .eq. '-') goto 2138
              if (res1(i+k1:i+k1) .ne. res2(j+k1:j+k1)) goto 2138
              k2 = k1 + 1
            end do
c
c ... K1 = shift (residues); value can be positive or negative
c     K3 = residues misaligned in sequence 1
c     K4 = residues misaligned in sequence 2
c
 2138       continue
            k1 = j - i
            k3 = 0
            k4 = 0
ccc            print *,'K1, K2, K3, K4 = ',k1,k2,k3,k4
            do k=i,i+k2-1
              if (res2(k:k) .ne. '-') k3 = k3 + 1
            end do
            do k=j,j+k2-1
              if (res1(k:k) .ne. '-') k4 = k4 + 1
            end do
c
ccc            print *,'K1, K2, K3, K4 = ',k1,k2,k3,k4
c
c ... check if sufficient residues involved
c
            if (k3 .lt. iremin .and. k4 .lt. iremin) then
ccc              print *,'SKIP !'
ccc              i = i + k2 - 1
              goto 2128
            end if
c
            i1 = 0
            i2 = 0
            if (k1 .gt. 0) then
              do k=i,i+k2-1
                do m=k,k+k1-1
                  if (res2(m:m) .ne. '-') i1 = i1 + 1
                end do
              end do
              do k=j,j+k2-1
                do m=k,k-k1+1,-1
                  if (res1(m:m) .ne. '-') i2 = i2 + 1
                end do
              end do
              i2 = -i2
            else
              do k=i,i+k2-1
                do m=k,k+k1+1,-1
                  if (res2(m:m) .ne. '-') i1 = i1 + 1
                end do
              end do
              i1 = -i1
              do k=j,j+k2-1
                do m=k,k-k1-1
                  if (res1(m:m) .ne. '-') i2 = i2 + 1
                end do
              end do
            end if
c
c ... calculated average per-residue shifts
c
            xdum  = float (i1) / float (k2)
            xdum2 = float (i2) / float (k2)
c
c ... print some information
c
            write (*,6140) k2,k3,(k2-k3),k4,(k2-k4),k1,xdum,xdum2
            if (logfil) write (iunit,6140)
     +        k2,k3,(k2-k3),k4,(k2-k4),k1,xdum,xdum2
c
 6140 format (/' Possible register shift !'/
     +         ' Nr of residues involved   : ',i12/
     +         ' Misaligned in sequence 1  : ',i12/
     +         ' Not aligned in sequence 1 : ',i12/
     +         ' Misaligned in sequence 2  : ',i12/
     +         ' Not aligned in sequence 2 : ',i12/
     +         ' Shift (incl. gaps !)      : ',i12/
     +         ' Average real shift seq 1  : ',f12.1/
     +         ' Average real shift seq 2  : ',f12.1)
c
c ... write residue numbers involved
c
            write (*,6145) 1,res1(i:i),chnami,lresi(i),
     +        res1(i+k2-1:i+k2-1),chnami,lresi(i+k2-1)
            write (*,6145) 2,res2(j:j),chnamj,lresj(j),
     +        res2(j+k2-1:j+k2-1),chnamj,lresj(j+k2-1)
c
            if (logfil) then
              write (iunit,6145) 1,res1(i:i),chnami,lresi(i),
     +          res1(i+k2-1:i+k2-1),chnami,lresi(i+k2-1)
              write (iunit,6145) 2,res2(j:j),chnamj,lresj(j),
     +          res2(j+k2-1:j+k2-1),chnamj,lresj(j+k2-1)
            end if
c
 6145 format (' Affected in seq ',i1,' : ',a1,' [',a1,i4,'] - ',
     +  a1,' [',a1,i4,']')
c
c ... show residues involved in upper case, rest in lower case
c
            i1 = max (1, min (i-5,j-5))
            i2 = min (nres, max(i+k2+5,j+k2+5))
c
            if (i1 .lt. i-1) call locase (res1(i1:i-1))
            if (i1 .lt. j-1) call locase (res2(i1:j-1))
            if (i2 .gt. i+k2) call locase (res1(i+k2:i2))
            if (i2 .gt. j+k2) call locase (res2(j+k2:i2))
c
            call textut (' Seq 1 :',res1(i1:i2))
            if (logfil) write (iunit,'(2a)') ' Seq 1 : ',res1(i1:i2)
            call textut (' Seq 2 :',res2(i1:i2))
            if (logfil) write (iunit,'(2a)') ' Seq 2 : ',res2(i1:i2)
c
c ... restore everything to upper case
c
            call upcase (res1(i1:i2))
            call upcase (res2(i1:i2))
c
c ... don't consider this part of the sequence anymore
c
            i = i + k2 - 1
c
            goto 2128
c
          end if
        end if
      end do
c
      goto 2128
c
c ... done with register-shift detection
c
 2133 continue
c
      write (*,*)
c
 9999 continue
c
      if (logfil) close (iunit)
c
      return
      end
c
c
c
      subroutine nwseq (imol,ic,jmol,jc,n1,n2,
     +                  x1,x2,cmpmat,maxaat,gappen,
     +                  fmat,pmat,seq1,seq2,iseq1,iseq2,rtij)
c
      include 'lsqman.incl'
c
      integer maxg11
      parameter (maxg11=2*maxres)
c
      integer imol,ic,jmol,jc,n1,n2,ierr,maxaat
      integer i,j,n1n,inow,jnow,n1d,nsi,k,k1,k2,nma
      integer pmat(0:n1,0:n2),iseq1(n1),iseq2(n2)
c
      real x1(3,n1),x2(3,n2),fmat(0:n1,0:n2)
      real rtij(12),cmpmat(maxaat,maxaat)
      real gappen,r1,r2,r3,xdum,xrmsd,xave,xsdv,xmin,xmax,xtot
      real qqq
c
      character*1 seq1(n1),seq2(n2),l1(maxg11),l2(maxg11),l3(maxg11)
c
code ...
c
      i = maxaat*maxaat
      call xstats (cmpmat,i,xave,xsdv,xmin,xmax,xtot)
      xmin = -2.0 * abs(xmin)
c
      write (*,*)
      call prompt (' Executing Needleman-Wunsch ...')
      do i=0,n1
        do j=0,n2
          fmat (i,j) = 0.0
          pmat (i,j) = -1
        end do
      end do
c
      do i=0,n1
        fmat (i,0) = -float(i)*gappen
        pmat (i,0) = 0
      end do
      do j=0,n2
        fmat (0,j) = -float(j)*gappen
        pmat (0,j) = 0
      end do
c
      do j=1,n2
        do i=1,n1
          if (iseq1(i) .gt. 0 .and. iseq2(j) .gt. 0) then
            r1 = fmat(i-1,j-1) + cmpmat(iseq1(i),iseq2(j))
          else
            r1 = fmat(i-1,j-1) + xmin
          end if
          r2 = fmat(i-1,j) - gappen
          r3 = fmat(i,j-1) - gappen
          if (r1 .ge. r2 .and. r1 .ge. r3) then
            fmat (i,j) = r1
            pmat (i,j) = 3
          else if (r2 .ge. r3) then
            fmat (i,j) = r2
            pmat (i,j) = 2
          else
            fmat (i,j) = r3
            pmat (i,j) = 1
          end if
        end do
      end do
c
      n1n = maxg11
      inow = n1
      jnow = n2
      n1d = 0
      nsi = 0
      k1 = n1
      k2 = n2
c
  100 continue
      if (pmat(inow,jnow) .eq. 3) then
        l1(n1n) = seq1(inow)
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
c
        x1(1,k1) = x1(1,inow)
        x1(2,k1) = x1(2,inow)
        x1(3,k1) = x1(3,inow)
        k1 = k1 - 1
        x2(1,k2) = x2(1,jnow)
        x2(2,k2) = x2(2,jnow)
        x2(3,k2) = x2(3,jnow)
        k2 = k2 - 1
c
        if (l1(n1n) .eq. l2(n1n) .and.
     +      l1(n1n) .ne. '?') then
          n1d = n1d + 1
          nsi = nsi + 1
          l3(n1n) = '|'
        else if (iseq1(inow) .gt. 0 .and. iseq2(jnow) .gt. 0) then
          if (cmpmat(iseq1(inow),iseq2(jnow)) .gt. 0) then
            nsi = nsi + 1
            l3(n1n) = '+'
          end if
        end if
        inow = inow - 1
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 2) then
        l1(n1n) = seq1(inow)
        l2(n1n) = '-'
        l3(n1n) = ' '
        inow = inow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 1) then
        l1(n1n) = '-'
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      end if
c
c      print *,' INOW, JNOW ',inow,jnow
c
      if (inow .eq. 0 .and. jnow .gt. 0) then
        do j=jnow,1,-1
          l1(n1n) = '-'
          l2(n1n) = seq2(j)
          l3(n1n) = ' '
          n1n = n1n - 1
        end do
      else if (inow .gt. 0 .and. jnow .eq. 0) then
        do i=inow,1,-1
          l1(n1n) = seq1(i)
          l2(n1n) = '-'
          l3(n1n) = ' '
          n1n = n1n - 1
        end do
      end if
c
      n1n = n1n + 1
c
 6010 format (1x,i6,1x,3(a1,1x),' RMSD = ',f8.2,' A')
 6020 format (1x,a10,1x,60a1)
c
      write (*,*)
      do i=n1n,maxg11,60
        j = min (i+59,maxg11)
        write (*,*)
        write (*,6020) 'Sequence 1 ',(l1(k),k=i,j)
        write (*,6020) '   |=ID    ',(l3(k),k=i,j)
        write (*,6020) 'Sequence 2 ',(l2(k),k=i,j)
      end do
c
      write (*,*)
      call fvalut (' Gap penalty         :',1,gappen)
      call rvalut (' Raw alignment score :',1,fmat(n1,n2))
      call ivalut (' Length sequence 1   :',1,n1)
      call ivalut (' Length sequence 2   :',1,n2)
      i = maxg11 - n1n + 1
      call ivalut (' Alignment length    :',1,i)
      call ivalut (' Nr of identities    :',1,n1d)
      xdum = 100.0 * float(n1d) / float(min(n1,n2))
      call fvalut (' Perc identities     :',1,xdum)
      call ivalut (' Nr of similarities  :',1,nsi)
      xdum = 100.0 * float(nsi) / float(min(n1,n2))
      call fvalut (' Perc similarities   :',1,xdum)
c
      k1 = k1 + 1
      k2 = k2 + 1
      i = (n1-k1+1)
      call ivalut (' Nr of matched res   :',1,i)
      if (i .ge. 3) then
        call lsqgjk (x1(1,k1),x2(1,k2),i,
     +               xrmsd,rtij,ierr)
        call fvalut (' RMSD (A) for those  :',1,xrmsd)
        call fvalut (' Operator stored     :',12,rtij)
      else
        call errcon ('Fewer than 3 residues aligned')
        return
      end if
c
c ... update statistics
c
      nma = i
      simind (imol,jmol) = xrmsd * float(min(n1,n2)) / float(nma)
      matchi (imol,jmol) = (1.0+float(nma)) /
     +       ( (1.0+rmswgt*xrmsd) * (1.0+float(min (n1,n2))) )
      call mcrho (nma,x1(1,k1),x2(1,k2),xrmsd,
     +            cripp(imol,jmol),rrmsd(imol,jmol),
     +            normsd(imol,jmol),xdum)
c
      nmatch (imol,jmol) = nma
      rmsd (imol,jmol) = xrmsd
      corb (imol,jmol) = 999.99999
      rmsb (imol,jmol) = 999.99999
      rmsdna (imol,jmol) = xrmsd / float(nma)
      qqq = 100.0/float(nma)
      sas1 (imol,jmol) = xrmsd * qqq
      sas2 (imol,jmol) = sas1 (imol,jmol) * qqq
      sas3 (imol,jmol) = sas2 (imol,jmol) * qqq
      sas4 (imol,jmol) = sas3 (imol,jmol) * qqq
c
      write (*,6030) nma,xrmsd,simind(imol,jmol),matchi(imol,jmol),
     +  cripp(imol,jmol),xdum,rrmsd(imol,jmol),normsd(imol,jmol),
     +  rmsdna(imol,jmol),sas1(imol,jmol),sas2(imol,jmol),
     +  sas3(imol,jmol),sas4(imol,jmol)
c
 6030 format (/
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' SI = RMS * Nmin / Nmatch             = ',f12.5/
     +  ' MI = (1+Nmatch)/{(1+W*RMS)*(1+Nmin)} = ',f12.5/
     +  ' CR = Maiorov-Crippen RHO (0-2)       = ',f12.5/
     +  ' Estimated RMSD for 2 random proteins = ',f10.3,' A'/
     +  ' RR = Relative RMSD                   = ',f12.5/
     +  ' NR = Normalised RMSD (100)           = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A'/
     +  ' SAS(1) = cRMS * (100/Nmatch)         = ',f10.3,' A'/
     +  ' SAS(2) = cRMS * (100/Nmatch)^2       = ',f10.3,' A'/
     +  ' SAS(3) = cRMS * (100/Nmatch)^3       = ',f10.3,' A'/
     +  ' SAS(4) = cRMS * (100/Nmatch)^4       = ',f10.3,' A')
c
      return
      end
c
c
c
      subroutine mulrms (imol,refchn,filnam,cdplot,dmax,cutoff,
     +                   ierr,prognm)
c
      include 'lsqman.incl'
c
      character*12 prognm
c
      integer mxband
      parameter (mxband = 1 + ((maxchn*maxchn) / 2))
c
      real cacoor(3,maxchn)
ccc      real dumdis(10)
      real xave,xsdv,xmin,xmax,xtot,pmax,cutoff,xrms,dis2
      real rxmin,rxmax,rymin,rymax,gxmin,gxmax,gymin,gymax
      real dmax,xoff,xdum,grey,yoff,ymin,ymax,x1,x2,ptsize,qmax
c
      integer iband(maxchn,maxchn),ibptr(maxchn)
      integer imol,ierr,i,j,k,i1,i2,ir1,ir2,nband
      integer nphi,iptr,jptr,ichn,ndone,leng1,nlist,ndis
c
      logical xinter
c
      character refchn*(1),filnam*(*),cdplot*(*),jchn*1,kchn*1
      character line*256,bandnm(0:mxband)*5
c
code ...
c
      ierr = 0
c
      dmax = max (0.1, min (100.0, dmax))
      qmax = 0.0
c
      do i=1,nchain(imol)
        if (chname(i,imol) .eq. refchn) then
          ichn = i
          goto 10
        end if
      end do
      call errcon ('Reference chain not found')
      ierr = -1
      return
c
   10 continue
      call textut (' Reference chain :',chname(ichn,imol))
      kchn = chname(ichn,imol)
      i1 = chnptr(1,ichn,imol)
      i2 = chnptr(2,ichn,imol)
      ir1 = iresid(i1,imol)
      ir2 = ir1
      do i=i1,i2
        ir1 = min (ir1,iresid(i,imol))
        ir2 = max (ir2,iresid(i,imol))
      end do
      write (*,'(1x,a,i5,a,i5)') 'Residue range : ',ir1,' - ',ir2
      call textut (' Central atom type :',atypes(1))
c
      if (atypes(1) .eq. 'ALL ' .or. atypes(1) .eq. 'NONH' .or.
     +    atypes(1) .eq. 'SIDE' .or. atypes(1) .eq. 'TRAC') then
        call errcon ('Invalid atom type !')
        ierr = -1
        return
      end if
c
      call fvalut (' Max dist on grey-scale (A) :',1,dmax)
      call fvalut (' Cut-off for printing   (A) :',1,cutoff)
c
      iptr = 1
      jptr = 1
c
      ndone = 0
      nlist = 0
c
      i = 0
      bandnm (0) = 'mRMSD'
      do k=1,nchain(imol)-1
        do j=k+1,nchain(imol)
          i = i + 1
          iband (k,j) = i
          bandnm(i) = chname(k,imol)//'-'//chname(j,imol)
        end do
      end do
      nband = i
      gxmin = float(ir1)
      gxmax = float(ir2+1)
      gymin = 0.0
      gymax = float (nband+1)
      call xps_open (iunit,cdplot,prognm)
      call xps_inquire (rxmin,rxmax,rymin,rymax)
      call xps_scale (gxmin,gxmax,gymin,gymax)
      call xps_stroke ()
      call xps_label ('Residue nr','Comparison')
c
c ... paint whole plot area pink so missing residues can easily
c     be detected
c
      call xps_filled_box (gxmin,gxmax,gymin,gymax,
     +      1.0,0.3,0.3)
c
c ... loop over the residues
c
      do i=ir1,ir2
c
        nphi = 0
        xoff = float(i)
c
        do j=1,nchain(imol)
          jchn = chname(j,imol)
          call getptr (imol,jchn,i,atypes(1),
     +                 -1,jchn,k,iptr,jptr,ierr)
          if (ierr .ne. 0) goto 6996
c
          nphi = nphi + 1
          ibptr (nphi) = j
          cacoor (1,nphi) = atmxyz (1,iptr,imol)
          cacoor (2,nphi) = atmxyz (2,iptr,imol)
          cacoor (3,nphi) = atmxyz (3,iptr,imol)
c
 6996     continue
c
        end do
c
ccc        print *,'I, NPHI ',i,nphi
c
        xrms = 0.0
        ndis = 0
c
        if (nphi .gt. 1) then
          do k=1,nphi-1
            do j=k+1,nphi
              ndis = ndis + 1
              dis2 = (cacoor(1,k)-cacoor(1,j))**2 +
     +               (cacoor(2,k)-cacoor(2,j))**2 +
     +               (cacoor(3,k)-cacoor(3,j))**2
              xrms = xrms + dis2
ccc              if (ndis .le. 10) dumdis(ndis) = sqrt(dis2)
              x1 = sqrt(dis2)
              xdum = max (0.0, min (dmax, x1))
              qmax = max (qmax, x1)
              grey = 1.0 - (xdum/dmax)
              yoff = float(iband(ibptr(k),ibptr(j)))
              call xps_filled_box (xoff,xoff+1.0,yoff,yoff+1.0,
     +          grey,grey,grey)
            end do
          end do
          xrms = sqrt ( xrms / float(ndis) )
          yoff = 0.0
          xdum = max (0.0, min (dmax, xrms))
          grey = 1.0 - (xdum/dmax)
          call xps_filled_box (xoff,xoff+1.0,yoff,yoff+1.0,
     +      grey,grey,grey)
          call xps_colour (1)
        end if
c
ccc        print *,'I, NPHI, NDIS ',i,nphi,ndis,xrms
c
        if (nphi .gt. 1) then
c
          ndone = ndone + 1
          buffb1 (ndone)          = float(i)
          buffb1 (maxres+ndone)   = xrms
c
          if (xrms .ge. cutoff) then
            write (*,6000) resnam(iptr,imol),i,xrms,nphi,ndis
            nlist = nlist + 1
ccc            call fvalut (' DISTANCES :',min(10,ndis),dumdis)
          end if
c
        end if
c
      end do
c
      call xps_colour (1)
      do i=1,nband
        yoff = float(i)
        call xps_move (gxmin,yoff)
        call xps_draw (gxmax,yoff)
      end do
c
      call xps_colour (0)
      call xps_move (gxmin,gymin)
      call xps_draw (gxmax,gymin)
      call xps_draw (gxmax,gymax)
      call xps_draw (gxmin,gymax)
      call xps_draw (gxmin,gymin)
c
      call xps_scale (rxmin,rxmax,rymin,rymax)
      xmin = 350.0
      xmax = 500.0
      ymin = 200.0
      ymax = 225.0
      do i=1,30
        x1 = xmin + float(i-1)*5.0
        x2 = x1 + 5.0
        grey = 1.0 - (float(i-1)/30.0)
        call xps_filled_box (x1,x2,ymin,ymax,grey,grey,grey)
      end do
c
      call xps_colour (1)
      call xps_move (xmin,ymin)
      call xps_draw (xmax,ymin)
      call xps_draw (xmax,ymax)
      call xps_draw (xmin,ymax)
      call xps_draw (xmin,ymin)
c
      call xps_colour (0)
      write (line,'(f10.2,a)') 0.0,' A'
      call pretty (line)
      call xps_text (xmin,ymin-16.0,12.0,line)
      write (line,'(f10.2,a)') dmax,' A'
      call pretty (line)
      call xps_text (xmin+125.0,ymin-16.0,12.0,line)
c
      ptsize = int((700.0-250.0)/float(nband+1))
      ptsize = min (ptsize,14.0)
      call fvalut (' Point size :',1,ptsize)
      call xps_scale (rxmin,rxmax,gymin,gymax)
      do i=0,nband
        k = leng1(bandnm(i))
        call xps_text (525.0,float(i)+0.1,ptsize,bandnm(i)(1:k))
      end do
      call xps_scale (gxmin,gxmax,gymin,gymax)
c
      line = '"CD plot" for mol ' // name(imol) // ' = file '
     +       // file(imol)
      call pretty (line)
      call xps_legend (line)
c
      write (line,*) 'Residue range ',ir1,' - ',ir2
      call pretty (line)
      call xps_legend (line)
c
      call fvalut (' Max observed distance (A) :',1,qmax)
      if (dmax .lt. qmax) then
        write (line,'(a,f10.2,a,f10.2,a)')
     +    'NOTE: DMAX_BLACK = ',dmax,
     +    ' A is smaller than max observed distance = ',
     +    qmax,' A !'
        call pretty (line)
        call xps_legend (line)
        line = ' ' // line
        call prompt (line)
      end if
c
      if (ndone .gt. 2) then
        call xstats (buffb1(maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
        write (line,6110) 'Multi-RMS dist (A)',xave,xsdv,xmin,xmax
        call pretty (line)
        call xps_legend (line)
      end if
c
      line = 'Any pink areas in the plot indicate residues' //
     +       ' missing in one or more chains'
      call pretty (line)
      call xps_legend (line)
c
      line = 'Reference for CD plots: T.A. Jones & G.J. Kleywegt' //
     +       ' (1999). CASP3 comparative'
      call pretty (line)
      call xps_legend (line)
c
      line = 'modelling evaluation. Proteins: Struct. Funct.' //
     +       ' Genet. Suppl. 3, 30-46.'
      call pretty (line)
      call xps_legend (line)
c
      call xps_close ()
      call prompt (' PostScript file written')
c
 6000 format (1x,a3,1x,i5,' | ',1(f8.2,:,' A | '),2i5)
c
      call jvalut (' Nr of residues found :',1,ndone)
      call jvalut (' Nr of residues shown :',1,nlist)
      if (ndone .lt. 3) return
c
      call xstats (buffb1(1),ndone,xave,xsdv,xmin,xmax,xtot)
      ir1 = nint(xmin)
      ir2 = nint(xmax)
c
 6110 format (1x,a,' Ave, Sdv, Min, Max : ',4f8.2)
c
      call xstats (buffb1(maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
      write (*,6110) 'Multi-RMS dist (A)',xave,xsdv,xmin,xmax
      pmax = max (0.1, 1.1*xmax)
c
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6200 format (a6,1x,8a)
 6210 format (a6,1x,8i8)
 6213 format (a6,1x,8f8.2)
 6215 format (a6,1x,a,8f8.2)
 6217 format (a6,1x,a,8i8)
 6220 format (12i6)
 6230 format (9f8.2)
c
      call stamp (line)
      write (iunit,6200) 'REMARK',line(1:leng1(line))
      write (iunit,6200) 'XLABEL','Residue number'
c
      write (iunit,6200) 'YLABEL',
     +    'Multi-RMS (distances between all unique pairs)'
      write (iunit,6200) 'REMARK',
     +    'Plot of Multi-RMS (distances between all unique ',
     +    'pairs) as a function of residue nr'
      write (iunit,6200) 'REMARK',
     +    'Atom type used ',atypes(1)
      write (iunit,6200) 'REMARK',
     +    'Values are for mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
      write (iunit,6213) 'XYVIEW',
     +    float(ir1-1),float(ir2+1),0.0,pmax
c
      write (iunit,6210) 'NPOINT',ndone
      write (iunit,6210) 'COLOUR',4
c
      write (iunit,6200) 'XVALUE','*'
      write (iunit,6220) (nint(buffb1(i)),i=1,ndone)
c
      call xstats (buffb1(maxres+1),ndone,xave,xsdv,xmin,xmax,xtot)
      write (iunit,6215) 'REMARK','Multi-RMS ave, sdv, min, max ',
     +  xave,xsdv,xmin,xmax
c
      write (line,'(a,i6,a)') '! MOLNAM_RESIDUE_MULTI_RMS R ',
     +    ndone,' (9f8.2)'
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
      write (iunit,6200) 'YVALUE','*'
c
      write (iunit,6230) (buffb1(maxres+i),i=1,ndone)
c
      write (iunit,6200) 'END   '
c
      close (iunit)
      call prompt (' Plot file written')
c
      ierr = 0
      return
c
      end
c
c
c
      subroutine nwimpr (imol,ic,jmol,jc,n1,n2,nwcut,
     +                   x1,x2,x3,x4,xmat,fmat,pmat,
     +                   seq1,seq2,rtij,mode,verbose)
c
      include 'lsqman.incl'
c
      integer maxg11
      parameter (maxg11=2*maxres)
c
      integer imol,ic,jmol,jc,n1,n2,ngjk
      integer i,j,n1n,inow,jnow,n1d,nsi,k,k1,k2,ierr
      integer pmat(0:n1,0:n2)
c
      real x1(3,n1),x2(3,n2),x3(3,n2),x4(3,n1)
      real xmat(n1,n2),fmat(0:n1,0:n2)
      real rmsent(0:maxg11),rtij(12),rt(12)
      real nwcut,gappen,r1,r2,r3,xdum,xrmsd,dsq,qqq
c
      logical lprint
c
      character*1 seq1(n1),seq2(n2),l1(maxg11),l2(maxg11),l3(maxg11)
      character number*10,mode*10,verbose*(*)
c
code ...
c
      lprint = (verbose(1:1) .eq. 'Y')
c
c ... first apply operator to coordinates of mol 2
c
      call vecrtv (x2,x3,n2,rtij(1),rtij(10))
c
      do i=1,n1
        do j=1,n2
          xmat (i,j) = 0.0
        end do
      end do
c     
      if (mode(1:1) .eq. 'D') then
        call prompt (' Calculating distance matrix ...')
        gappen = 0.5*nwcut
        do j=1,n2
          do i=1,n1
            xrmsd = (x1(1,i)-x3(1,j))**2 +
     +              (x1(2,i)-x3(2,j))**2 +
     +              (x1(3,i)-x3(3,j))**2
            xmat (i,j) = - sqrt(xrmsd)
          end do
        end do
      else
        call prompt (' Calculating squared distance matrix ...')
        gappen = 0.5*nwcut*nwcut
        do j=1,n2
          do i=1,n1
            xrmsd = (x1(1,i)-x3(1,j))**2 +
     +              (x1(2,i)-x3(2,j))**2 +
     +              (x1(3,i)-x3(3,j))**2
            xmat (i,j) = - xrmsd
          end do
        end do
      end if
c
      write (*,*)
      call prompt (' Executing Needleman-Wunsch ...')
      do i=0,n1
        do j=0,n2
          fmat (i,j) = 0.0
          pmat (i,j) = -1
        end do
      end do
c
      do i=0,n1
        fmat (i,0) = -float(i)*gappen
        pmat (i,0) = 0
      end do
      do j=0,n2
        fmat (0,j) = -float(j)*gappen
        pmat (0,j) = 0
      end do
c
      do j=1,n2
        do i=1,n1
          r1 = fmat(i-1,j-1) + xmat(i,j)
          r2 = fmat(i-1,j) - gappen
          r3 = fmat(i,j-1) - gappen
          if (r1 .ge. r2 .and. r1 .ge. r3) then
            fmat (i,j) = r1
            pmat (i,j) = 3
          else if (r2 .ge. r3) then
            fmat (i,j) = r2
            pmat (i,j) = 2
          else
            fmat (i,j) = r3
            pmat (i,j) = 1
          end if
        end do
      end do
c
      n1n = maxg11
      inow = n1
      jnow = n2
      n1d = 0
      nsi = 0
      k1 = n1
      k2 = n2
      ngjk = 0
c
  100 continue
      if (pmat(inow,jnow) .eq. 3) then
        ngjk = ngjk + 1
        l1(n1n) = seq1(inow)
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
        if (mode(1:1) .eq. 'D') then
          rmsent(n1n) = - xmat(inow,jnow)
        else
          rmsent(n1n) = sqrt (- xmat(inow,jnow))
        end if
c
        x4(1,k1) = x1(1,inow)
        x4(2,k1) = x1(2,inow)
        x4(3,k1) = x1(3,inow)
        k1 = k1 - 1
        x3(1,k2) = x2(1,jnow)
        x3(2,k2) = x2(2,jnow)
        x3(3,k2) = x2(3,jnow)
        k2 = k2 - 1
c
        if (l1(n1n) .eq. l2(n1n) .and.
     +      l1(n1n) .ne. '?') then
          n1d = n1d + 1
          nsi = nsi + 1
          l3(n1n) = '|'
c        else if (cmpmat(iptr(inow),jptr(jnow)) .gt. 0) then
c          nsi = nsi + 1
c          l3(n1n) = '+'
        end if
        inow = inow - 1
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 2) then
        l1(n1n) = seq1(inow)
        l2(n1n) = '-'
        l3(n1n) = ' '
        rmsent(n1n) = -999.999
        inow = inow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      else if (pmat(inow,jnow) .eq. 1) then
        l1(n1n) = '-'
        l2(n1n) = seq2(jnow)
        l3(n1n) = ' '
        rmsent(n1n) = -999.999
        jnow = jnow - 1
        n1n = n1n - 1
        if (inow .gt. 0 .and. jnow .gt. 0) goto 100
      end if
c
c      print *,' INOW, JNOW ',inow,jnow
c
      if (inow .eq. 0 .and. jnow .gt. 0) then
        do j=jnow,1,-1
          l1(n1n) = '-'
          l2(n1n) = seq2(j)
          l3(n1n) = ' '
          rmsent(n1n) = -999.999
          n1n = n1n - 1
        end do
      else if (inow .gt. 0 .and. jnow .eq. 0) then
        do i=inow,1,-1
          l1(n1n) = seq1(i)
          l2(n1n) = '-'
          l3(n1n) = ' '
          rmsent(n1n) = -999.999
          n1n = n1n - 1
        end do
      end if
c
      n1n = n1n + 1
c
      if (lprint) then
        write (*,*)
        k = 0
        do i=n1n,maxg11
          k = k + 1
          if (rmsent(i) .ge. 0.0) then
            write (number,'(f8.2,1x,a1)') rmsent(i),'A'
          else
            number = '       -'
          end if
          write (*,6010) k,l1(i),l3(i),l2(i),number
        end do
      end if
c
 6010 format (1x,i6,1x,3(a1,1x),' DIST = ',a10)
 6020 format (1x,a10,1x,60a1)
c
      if (lprint) then
        write (*,*)
        do i=n1n,maxg11,60
          j = min (i+59,maxg11)
          write (*,*)
          write (*,6020) 'Sequence 1 ',(l1(k),k=i,j)
          write (*,6020) '   |=ID    ',(l3(k),k=i,j)
          write (*,6020) 'Sequence 2 ',(l2(k),k=i,j)
        end do
      end if
c
      write (*,*)
      call fvalut (' Gap penalty         :',1,gappen)
      call rvalut (' Raw alignment score :',1,fmat(n1,n2))
      call ivalut (' Length sequence 1   :',1,n1)
      call ivalut (' Length sequence 2   :',1,n2)
      i = maxg11 - n1n + 1
      call ivalut (' Alignment length    :',1,i)
      call ivalut (' Nr of identities    :',1,n1d)
      xdum = 100.0 * float(n1d) / float(min(n1,n2))
      call fvalut (' Perc identities     :',1,xdum)
c      call ivalut (' Nr of similarities  :',1,nsi)
c      xdum = 100.0 * float(nsi) / float(min(n1,n2))
c      call fvalut (' Perc similarities   :',1,xdum)
c
      k1 = k1 + 1
      k2 = k2 + 1
      i = (n1-k1+1)
      call ivalut (' Nr of matched res   :',1,i)
ccc      call ivalut (' NGJK :',1,ngjk)
      if (i .ge. 3) then
        call lsqgjk (x4(1,k1),x3(1,k2),i,xrmsd,rt,ierr)
        call fvalut (' RMSD for those (A)  :',1,xrmsd)
c
        simind (imol,jmol) = xrmsd * float(min(n1,n2)) / float(i)
        matchi (imol,jmol) = (1.0+float(i)) /
     +       ( (1.0+rmswgt*xrmsd) * (1.0+float(min (n1,n2))) )
        call mcrho (i,x4(1,k1),x3(1,k2),xrmsd,
     +              cripp(imol,jmol),rrmsd(imol,jmol),
     +              normsd(imol,jmol),dsq)
c
        nmatch (imol,jmol) = i
        rmsd (imol,jmol) = xrmsd
        corb (imol,jmol) = 999.99999
        rmsb (imol,jmol) = 999.99999
        rmsdna (imol,jmol) = xrmsd / float(i)
        qqq = 100.0/float(i)
        sas1 (imol,jmol) = xrmsd * qqq
        sas2 (imol,jmol) = sas1 (imol,jmol) * qqq
        sas3 (imol,jmol) = sas2 (imol,jmol) * qqq
        sas4 (imol,jmol) = sas3 (imol,jmol) * qqq
c
        write (*,6030) i,xrmsd,simind(imol,jmol),matchi(imol,jmol),
     +                cripp(imol,jmol),dsq,rrmsd(imol,jmol),
     +                normsd(imol,jmol),rmsdna(imol,jmol),
     +                sas1(imol,jmol),sas2(imol,jmol),
     +                sas3(imol,jmol),sas4(imol,jmol),
     +                rt(1),rt(4),rt(7),rt(2),rt(5),rt(8),
     +                rt(3),rt(6),rt(9),rt(10),rt(11),rt(12)
c
        do j=1,12
          rtij(j) = rt(j)
        end do
      end if
c
 6030 format (/
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' SI = RMS * Nmin / Nmatch             = ',f12.5/
     +  ' MI = (1+Nmatch)/{(1+W*RMS)*(1+Nmin)} = ',f12.5/
     +  ' CR = Maiorov-Crippen RHO (0-2)       = ',f12.5/
     +  ' Estimated RMSD for 2 random proteins = ',f10.3,' A'/
     +  ' RR = Relative RMSD                   = ',f12.5/
     +  ' NR = Normalised RMSD (100)           = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A'/
     +  ' SAS(1) = cRMS * (100/Nmatch)         = ',f10.3,' A'/
     +  ' SAS(2) = cRMS * (100/Nmatch)^2       = ',f10.3,' A'/
     +  ' SAS(3) = cRMS * (100/Nmatch)^3       = ',f10.3,' A'/
     +  ' SAS(4) = cRMS * (100/Nmatch)^4       = ',f10.3,' A'/
     +  ' Rotation     : ',3f12.8/16x,3f12.8/16x,3f12.8/
     +  ' Translation  : ',3f12.4)
c
      return
      end
c
c
c
      subroutine nwsoap (n1,n2,filnam,
     +                   x1,x2,xmat,fmat,pmat,
     +                   rtij,verbose)
c
      include 'lsqman.incl'
c
      integer i,j,inow,jnow,k,ierr,n1,n2
      integer pmat(0:n1,0:n2)
c
      real x1(3,n1),x2(3,n2),fmat(0:n1,0:n2),xmat(1:n1,1:n2)
      real rtij(12)
      real r1,r2,area3,sum,delta,test
c
      logical xinter,lprint
c
      character filnam*(*),verbose*(*)
c
code ...
c
      lprint = (verbose(1:1) .eq. 'Y')
c
      close (iunit)
      call xopxua (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening ODL file')
        return
      end if
c
      write (iunit,'(a)') 'begin soap'
      write (iunit,'(a)') 'mode solid'
      write (iunit,'(a)') 'colour sky_blue'
c
c ... first apply operator to coordinates of mol 2
c
      write (*,*)
      call fvalut (' Applying current operator to mol 2 :',12,rtij)
c      call fvalut (' BEFORE :',3,x2(1,1))
      call vecrtv (x2,x2,n2,rtij(1),rtij(10))
c      call fvalut (' AFTER  :',3,x2(1,1))
c
      write (*,*)
      call prompt (' Calculating triangle-area matrix ...')
c
c      do i=1,n1
c        do j=1,n2
c          xmat (i,j) = 0.0
c        end do
c      end do
c
c      do j=2,n2
c        do i=2,n1
c          xmat(i,j) = area3 (x1(1,i-1),x1(1,i),x2(1,j))
c          xmat(j,i) = area3 (x1(1,i),x2(1,j-1),x2(1,j))
c          if (i .eq. j) xmat(i,i) = 9999.99
c        end do
c      end do
c
      write (*,*)
      call prompt (' Executing Needleman-Wunsch ...')
      do i=0,n1
        do j=0,n2
          fmat (i,j) = 0.0
          pmat (i,j) = -1
        end do
      end do
c
      do i=0,n1
        fmat (i,0) = 0.0
        pmat (i,0) = 1
        fmat (i,1) = 0.0
        if (i .gt. 1) fmat (i,1) = fmat (i-1,1) +
     +    area3 (x1(1,i-1),x1(1,i),x2(1,1))
        pmat (i,1) = 1
      end do
      do j=0,n2
        fmat (0,j) = 0.0
        pmat (0,j) = 2
        fmat (1,j) = 0.0
        if (j .gt. 1) fmat (1,j) = fmat (1,j-1) +
     +    area3 (x1(1,1),x2(1,j-1),x2(1,j))
        pmat (1,j) = 2
      end do
c
      do j=2,n2
        do i=2,n1
c          if (i .ne. j) then
c            r1 = fmat(i-1,j) + xmat(i,j)
c            r2 = fmat(i,j-1) + xmat(j,i)          
c          else
            r1 = fmat(i-1,j) + area3 (x1(1,i-1),x1(1,i),x2(1,j))
            r2 = fmat(i,j-1) + area3 (x1(1,i),x2(1,j-1),x2(1,j))          
c          end if
          if (r1 .le. r2) then
            fmat (i,j) = r1
            pmat (i,j) = 1
          else
            fmat (i,j) = r2
            pmat (i,j) = 2
          end if
        end do
      end do
c
      inow = n1
      jnow = n2
      sum = fmat(inow,jnow)
      test = 0.0
c
  100 continue
c
ccc      if (inow .gt. 1 .and. jnow .gt. 1) print *,inow,jnow,
ccc     +  pmat(inow,jnow),
ccc     +  area3 (x1(1,inow-1),x1(1,inow),x2(1,jnow)),
ccc     +  area3 (x1(1,inow),x2(1,jnow-1),x2(1,jnow))
c
      if (pmat(inow,jnow) .eq. 1) then
        write (iunit,'(a)') ' poly 3'
        write (iunit,'(3x,3f10.3)') (x1(k,inow-1),k=1,3)
        write (iunit,'(3x,3f10.3)') (x1(k,inow),k=1,3)
        write (iunit,'(3x,3f10.3)') (x2(k,jnow),k=1,3)
        delta = sum - fmat(inow-1,jnow)
        sum = fmat(inow-1,jnow)
        if (lprint) write (*,6010) inow-1,inow,jnow,delta
        test = test + delta
ccc        print *,'TYPE 1 -> ',inow-1,inow,jnow,delta
        inow = inow - 1
        if (inow .gt. 1 .or. jnow .gt. 1) goto 100
      else if (pmat(inow,jnow) .eq. 2) then
        write (iunit,'(a)') ' poly 3'
        write (iunit,'(3x,3f10.3)') (x1(k,inow),k=1,3)
        write (iunit,'(3x,3f10.3)') (x2(k,jnow-1),k=1,3)
        write (iunit,'(3x,3f10.3)') (x2(k,jnow),k=1,3)
        delta = sum - fmat(inow,jnow-1)
        sum = fmat(inow,jnow-1)
        test = test + delta
        if (lprint) write (*,6020) inow,jnow-1,jnow,delta
ccc        print *,'TYPE 2 -> ',inow,jnow-1,jnow,delta
        jnow = jnow - 1
        if (inow .gt. 1 .or. jnow .gt. 1) goto 100
      end if
c
ccc      print *,'SUM (zero ?)   ',sum
ccc      print *,'TEST (total ?) ',test
c
 6010 format (' Type 1 -> I-1,I,J = ',3i6,' Area = ',f10.2,' A2')
 6020 format (' Type 2 -> I,J-1,J = ',3i6,' Area = ',f10.2,' A2')
c
      write (iunit,'(a)') 'end_object'
c
      close (iunit)
c
      call prompt (' ODL file written')
      call rvalut (' Total area (A2) :',1,fmat(n1,n2))
c
      return
      end
c
c
c
      subroutine distat (nd,xd)
c
      implicit none
c
      integer nd
c
      real xd(nd)
      real x1,x2,x3,x4,x5
c
      integer i,j,k
c
code ...
c
      if (nd .lt. 3) return
c
      call xstats (xd,nd,x1,x2,x3,x4,x5)
c
      write (*,6000) nd,x1,x2,(x2*x2),x3,x4,x4-x3,x5
c
 6000 format (' Analysis of distance distribution:'/
     +  ' Number of distances                    : ',i10/
     +  ' Average (A)                            : ',f10.2/
     +  ' Standard deviation (A)                 : ',f10.2/
     +  ' Variance (A**2)                        : ',f10.2/
     +  ' Minimum (A)                            : ',f10.2/
     +  ' Maximum (A)                            : ',f10.2/
     +  ' Range (A)                              : ',f10.2/
     +  ' Sum (A)                                : ',f10.2)
c
      call xstat2 (xd,nd,x1,x2)
c
ccc      call fvalut (' XD array:',nd,xd)
c
      write (*,6010) x1,x2
c
 6010 format (
     +  ' Root-mean-square (A)                   : ',f10.2/
     +  ' Harmonic average (A)                   : ',f10.2)
c
      if (nd .lt. 10) return
c
c ... heap-sort
c
      call hsortr (nd,xd)
c
c ... median
c
      i = nd/2
      if ( (2*i) .eq. nd ) then
        x1 = 0.5 * (xd(i) + xd(i+1))
      else
        x1 = xd (i)
      end if
c
      i = nd/4
      x2 = 0.5 * (xd(i) + xd(i+1))
c
      j = (3*nd)/4
      x3 = 0.5 * (xd(j) + xd(j+1))
c
      x4 = 0.25 * (x2 + x1 + x1 + x3)
c
      x5 = 0.0
      do k=i+1,j
        x5 = x5 + xd (k)
      end do
      x5 = x5 / float (j-i)
c
      write (*,6020) x1,x2,x3,x3-x2,x4,x5
c
 6020 format (
     +  ' Median (A)                             : ',f10.2/
     +  ' 25th Percentile (A)                    : ',f10.2/
     +  ' 75th Percentile (A)                    : ',f10.2/
     +  ' Semi-interquartile range (A)           : ',f10.2/
     +  ' Trimean (A)                            : ',f10.2/
     +  ' 50% Trimmed mean (A)                   : ',f10.2)
c
      i = nd/10
      x1 = 0.5 * (xd(i) + xd(i+1))
c
      j = (9*nd)/10
      x2 = 0.5 * (xd(j) + xd(j+1))
c
      x3 = 0.0
      do k=i+1,j
        x3 = x3 + xd (k)
      end do
      x3 = x3 / float (j-i)
c
      write (*,6030) x1,x2,x3
c
 6030 format (
     +  ' 10th Percentile (A)                    : ',f10.2/
     +  ' 90th Percentile (A)                    : ',f10.2/
     +  ' 20% Trimmed mean (A)                   : ',f10.2)
c
      return
      end
c
c
c
      subroutine extali (imol,ic,jmol,jc,n1,n2,
     +                  x1,x2,alifil,
     +                  seq1,seq2,rtij)
c
      include 'lsqman.incl'
c
      integer maxg11
      parameter (maxg11=2*maxres)
c
      integer imol,ic,jmol,jc,n1,n2,ierr,ni1,ni2,nma
      integer i,j,n1d,k,k1,k2,length
c
      real x1(3,n1),x2(3,n2)
      real rtij(12)
      real xdum,xrmsd,qqq
c
      logical xinter
c
      character*1 seq1(*),seq2(*),l1(maxg11),l2(maxg11),l3(maxg11)
      character line*128,alifil*(*)
c
code ...
c
      close (iunit)
      call xopxua (iunit,alifil,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening alignment file')
        return
      end if
c
      write (*,*)
      call prompt (' Reading external alignment ...')
c
      k = 1
      n1d = 0
c
   10 continue
      read (iunit,'(a)',err=9000,end=9000) line
      if (line(1:1) .ne. '>') goto 10
c
      call textut (' Found sequence :',line)
      read (iunit,'(a)',err=9000,end=9000) line
      call textut (' Title :',line)
      i = 0
   11 continue
      read (iunit,'(a)',err=9000,end=9000) line
      if (line(1:1) .eq. '*') goto 15
      call subchr (line,'X','?',j)
      call subchr (line,'.','-',j)
      do j=1,length(line)
        i = i + 1
        if (k .eq. 1) then
          l1 (i) = line(j:j)
        else
          l2 (i) = line(j:j)
          l3 (i) = ' '
          if (l1(i) .eq. l2(i) .and. l1(i) .ne. '-' .and.
     +        l1(i) .ne. '?') then
            l3(i) = '|'
            n1d = n1d + 1
          end if
        end if
      end do
      goto 11
c
   15 continue
      call jvalut (' Length (incl. gaps) :',1,i)
      if (k .eq. 1) then
        k1 = i
        k = 2
        goto 10
      end if
c
      k2 = i
      close (iunit)
c
      if (k1 .ne. k2) then
        call errcon ('Sequences have different lengths')
        return
      end if
c
 6010 format (1x,i6,1x,3(a1,1x),' RMSD = ',f8.2,' A')
 6020 format (1x,a10,1x,60a1)
c
      write (*,*)
      call prompt (' Sequences deduced from PDB files:')
      write (*,*)
      do i=1,n1,60
        j = min (i+59,n1)
        write (*,6020) 'Sequence 1 ',(seq1(k),k=i,j)
      end do
      write (*,*)
      do i=1,n2,60
        j = min (i+59,n2)
        write (*,6020) 'Sequence 2 ',(seq2(k),k=i,j)
      end do
c
      write (*,*)
      call prompt (' Sequences read from alignment file:')
      do i=1,k1,60
        j = min (i+59,k1)
        write (*,*)
        write (*,6020) 'Sequence 1 ',(l1(k),k=i,j)
        write (*,6020) '   |=ID    ',(l3(k),k=i,j)
        write (*,6020) 'Sequence 2 ',(l2(k),k=i,j)
      end do
c
c ... check integrity of the sequences
c
      ni1 = 0
      nma = 0
      ni2 = 0
c
      write (*,*)
      call prompt (' Checking integrity of sequences ...')
c
      do i=1,k1
c
        if (l1(i) .ne. '-') then
          ni1 = ni1 + 1
          if (seq1(ni1) .ne. l1(i)) goto 8000
        end if
c
        if (l2(i) .ne. '-') then
          ni2 = ni2 + 1
          if (seq2(ni2) .ne. l2(i)) goto 8100
        end if
c
        if (l1(i) .ne. '-' .and. l2(i) .ne. '-') then
          nma = nma + 1
          x1(1,nma) = x1(1,ni1)
          x1(2,nma) = x1(2,ni1)
          x1(3,nma) = x1(3,ni1)
          x2(1,nma) = x2(1,ni2)
          x2(2,nma) = x2(2,ni2)
          x2(3,nma) = x2(3,ni2)
        end if
c
      end do
c
      call prompt (' PDB and alignment-file sequences identical !')
      call jvalut (' Nr of aligned residues :',1,nma)
c
      write (*,*)
      call ivalut (' Length sequence 1   :',1,n1)
      call ivalut (' Length sequence 2   :',1,n2)
      call ivalut (' Alignment length    :',1,k1)
      call ivalut (' Nr of identities    :',1,n1d)
      xdum = 100.0 * float(n1d) / float(min(n1,n2))
      call fvalut (' Perc identities     :',1,xdum)
c
      if (nma .lt. 3) then
        call errcon ('Too few aligned residues (< 3) !')
        return
      end if
c
      write (*,*)
      call lsqgjk (x1,x2,nma,xrmsd,rtij,ierr)
      call jvalut (' Nr of aligned residues :',1,nma)
      call fvalut (' RMSD (A) for those     :',1,xrmsd)
      call fvalut (' Operator stored        :',12,rtij)
c
c ... update statistics
c
      simind (imol,jmol) = xrmsd * float(min(n1,n2)) / float(nma)
      matchi (imol,jmol) = (1.0+float(nma)) /
     +       ( (1.0+rmswgt*xrmsd) * (1.0+float(min (n1,n2))) )
      call mcrho (nma,x1,x2,xrmsd,
     +            cripp(imol,jmol),rrmsd(imol,jmol),
     +            normsd(imol,jmol),xdum)
c
      nmatch (imol,jmol) = nma
      rmsd (imol,jmol) = xrmsd
      corb (imol,jmol) = 999.99999
      rmsb (imol,jmol) = 999.99999
      rmsdna (imol,jmol) = xrmsd / float(nma)
      qqq = 100.0/float(nma)
      sas1 (imol,jmol) = xrmsd * qqq
      sas2 (imol,jmol) = sas1 (imol,jmol) * qqq
      sas3 (imol,jmol) = sas2 (imol,jmol) * qqq
      sas4 (imol,jmol) = sas3 (imol,jmol) * qqq
c
      write (*,6030) nma,xrmsd,simind(imol,jmol),matchi(imol,jmol),
     +               cripp(imol,jmol),xdum,rrmsd(imol,jmol),
     +               normsd(imol,jmol),rmsdna(imol,jmol),
     +               sas1(imol,jmol),sas2(imol,jmol),
     +               sas3(imol,jmol),sas4(imol,jmol)
c
 6030 format (/
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' SI = RMS * Nmin / Nmatch             = ',f12.5/
     +  ' MI = (1+Nmatch)/{(1+W*RMS)*(1+Nmin)} = ',f12.5/
     +  ' CR = Maiorov-Crippen RHO (0-2)       = ',f12.5/
     +  ' Estimated RMSD for 2 random proteins = ',f10.3,' A'/
     +  ' RR = Relative RMSD                   = ',f12.5/
     +  ' NR = Normalised RMSD (100)           = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A'/
     +  ' SAS(1) = cRMS * (100/Nmatch)         = ',f10.3,' A'/
     +  ' SAS(2) = cRMS * (100/Nmatch)^2       = ',f10.3,' A'/
     +  ' SAS(3) = cRMS * (100/Nmatch)^3       = ',f10.3,' A'/
     +  ' SAS(4) = cRMS * (100/Nmatch)^4       = ',f10.3,' A')
c
      return
c
c ... sequence 1 inconsistency
c
 8000 continue
      call errcon ('Inconsistency for sequence 1')
      write (*,6800) ni1,seq1(ni1),i,l1(i)
 6800 format (
     +  ' Sequential position in PDB sequence = ',i6,' type = ',a1/
     +  ' Sequential position in alignment    = ',i6,' type = ',a1)
c
      return
c
c ... sequence 2 inconsistency
c
 8100 continue
      call errcon ('Inconsistency for sequence 2')
      write (*,6810) ni2,seq2(ni2),i,l2(i)
 6810 format (
     +  ' Sequential position in PDB sequence = ',i6,' type = ',a1/
     +  ' Sequential position in alignment    = ',i6,' type = ',a1)
c
      return
c
c ... read error
c
 9000 continue
      call errcon ('While reading alignment file')
      close (iunit)
c
      return
      end
