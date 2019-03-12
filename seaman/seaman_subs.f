c
c ... subroutines for seaman.f
c
      subroutine readdb (iunit,ierr)
c
      include 'seaman.incl'
c
      integer iunit,ierr,nl,i,j
c
      character line*256
c
code ...
c
      ntype = 0
      nl = 0
      ierr = 0
c
   10 continue
      read (iunit,'(a)',end=100,err=9999) line
      nl = nl + 1
c
      if (line(1:1) .eq. '!') goto 10
      call upcase (line)
c
      if (line(1:4) .eq. 'RESI') then
        if (ntype .eq. maxtyp) then
          call errcon ('Too many residue types')
          ierr = -1
          return
        end if
        ntype = ntype + 1
        if (ntype .gt. 1) then
          if (natype(ntype-1) .le. 0) ntype=ntype-1
        end if
        read (line(5:),*,err=9998) typnam(ntype)
        natype(ntype) = 0
        goto 10
      end if
c
      if (line(1:4) .eq. 'ATOM') then
        if (ntype .lt. 1) then
          call errcon ('Atom card before Residue card')
          ierr = -1
          return
        end if
        natype(ntype) = natype(ntype) + 1
        i = natype(ntype)
        if (i .gt. maxatp) then
          call errcon ('Too many atom types')
          ierr = -1
          return
        end if
        read (line,6000,err=9997) typanm (i,ntype),
     +    (typxyz(j,i,ntype),j=1,3)
        goto 10
      end if
c
c ... on END, sort N,CA,C,CB,O,rest
c
      if (line(1:4) .eq. 'END ') then
        if (typanm(1,ntype) .ne. ' N  ') then
          do i=2,natype(ntype)
            if (typanm(i,ntype) .eq. ' N  ') then
              call typswp (ntype,1,i)
              goto 20
            end if
          end do
          goto 9010
        end if
   20   continue
        if (typanm(2,ntype) .ne. ' CA ') then
          do i=3,natype(ntype)
            if (typanm(i,ntype) .eq. ' CA ') then
              call typswp (ntype,2,i)
              goto 22
            end if
          end do
          goto 9010
        end if
   22   continue
        if (typanm(3,ntype) .ne. ' C  ') then
          do i=4,natype(ntype)
            if (typanm(i,ntype) .eq. ' C  ') then
              call typswp (ntype,3,i)
              goto 24
            end if
          end do
          goto 9010
        end if
   24   continue
        if (typanm(4,ntype) .ne. ' O  ') then
          do i=4,natype(ntype)
            if (typanm(i,ntype) .eq. ' O  ') then
              call typswp (ntype,4,i)
              goto 26
            end if
          end do
          goto 9010
        end if
   26   continue
        goto 10
      end if
c
      call errcon ('Unrecognised line')
      call textut (' >',line)
c
      goto 10
c
 6000 format (12x,a4,14x,3f8.3)
c
cATOM      1  N   TRP A  27       2.123   1.324   0.720  1.00 20.00   6
c1234567890123456789012345678901234567890123456789012345678901234567890
c
  100 continue
      call jvalut (' Lines read    :',1,nl)
      call jvalut (' Residue types :',1,ntype)
c
  101 continue
      close (iunit)
      return
c
 9997 call errcon ('While reading ATOM card')
      ierr = -1
      goto 101
c
 9998 call errcon ('While reading RESI card')
      ierr = -1
      goto 101
c
 9999 call errcon ('While reading from file')
      ierr = -1
      goto 101
c
 9010 call errcon ('Missing main-chain atom type')
      ierr = -1
      goto 101
c
      end
c
c
c
      subroutine typswp (it,i1,i2)
c
      include 'seaman.incl'
c
      integer it,i1,i2
c
      character name*4
c
code ...
c
      name = typanm(i1,it)
      typanm(i1,it) = typanm(i2,it)
      typanm(i2,it) = name
c
      call rswap (typxyz(1,i1,it),typxyz(1,i2,it))
      call rswap (typxyz(2,i1,it),typxyz(2,i2,it))
      call rswap (typxyz(3,i1,it),typxyz(3,i2,it))
c
      return
      end
c
c
c
      subroutine putpdb (iunit,ierr)
c
      include 'seaman.incl'
c
      real ca2fa(12)
c
      integer iunit,ierr,kk,length,i,leng1
c
      character key*6,line*256
c
code ...
c
      ierr = 0
c
 4711 format (a6,10a)
 4713 format (i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2,6x,a4)
c
      call stamp (line)
      write (iunit,4711,err=997) 'REMARK',' ',line(1:leng1(line))
c
c ... write cell etc. if known
c
      if (lcell) then
        call prompt (' Including CRYST1 etc. cards')
c
        write (iunit,'(a6,3f9.3,3f7.2,1x,a11,i4)',err=997) 
     +    'CRYST1',(cell(i),i=1,6)
c
        write (iunit,'(a)',err=997)
     +    'ORIGX1      1.000000  0.000000  0.000000        0.00000'
        write (iunit,'(a)',err=997)
     +    'ORIGX2      0.000000  1.000000  0.000000        0.00000'
        write (iunit,'(a)',err=997)
     +    'ORIGX3      0.000000  0.000000  1.000000        0.00000'
c
        do i=1,12
          ca2fa (i) = 0.0
        end do
c
        call orthog (cell,ca2fa,1)
c
        write (iunit,'(a6,4x,3f10.6,f15.5)',err=997)
     +    'SCALE1',(ca2fa(i),i=1,10,3)
        write (iunit,'(a6,4x,3f10.6,f15.5)',err=997)
     +    'SCALE2',(ca2fa(i),i=2,11,3)
        write (iunit,'(a6,4x,3f10.6,f15.5)',err=997)
     +    'SCALE3',(ca2fa(i),i=3,12,3)
c
      else
        call prompt (' Cell unknown; no CRYST1 etc. cards written')
      end if
c
      key = 'ATOM  '
c
      do kk=1,natoms
c
        write (line,4713,err=5372) kk,
     +    atmnam(kk),resnam(kk),achain(kk),
     +    iresid(kk),atmxyz(1,kk),atmxyz(2,kk),
     +    atmxyz(3,kk),qatom(kk),batom(kk)
        write (iunit,4711,err=997) key,line(1:leng1(line))
c
 1938   continue
      end do
c
      write (iunit,'(a)',err=997) 'END'
      call jvalut (' Number of atoms written :',1,natoms)
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
      subroutine getpdb (iunit,ierr)
c
      include 'seaman.incl'
c
      integer iunit,ierr,nl,i,j,ichain,jchain,iprev,nhet
      integer length,nstrip,nunk,kk,leng1
c
      logical nmr,lhydro
c
      character key*6,line*256,prev*1,chain*1
c
code ...
c
      nl = 0
      ierr = 0
      ichain = ichar('A') - 1
      jchain = 0
      prev = '?'
      nmr = .false.
      iprev = 999999
      nhet = 0
      nunk = 0
      nstrip = 0
      natoms = 0
      lcell = .false.
c
c ... read the PDB file
c
 4711 format (a6,a)
 4713 format (i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2,6x,a4)
c
cATOM      1  N   TRP A  27       2.123   1.324   0.720  1.00 20.00   6
c1234567890123456789012345678901234567890123456789012345678901234567890
c
   10 continue
      read (iunit,4711,end=21,err=997) key,line
      nl = nl + 1
      call upcase (key)
c
      if (key .eq. 'MODEL ') then
        if (.not. nmr) call prompt (' Multiple NMR models')
        nmr = .true.
        ichain = ichain + 1
        if (ichain .gt. ichar('Z')) then
          ichain = ichain - 1
          call errcon ('Too many MODELs; max = 26; rest skipped')
          goto 21
        end if
        chain = char(ichain)
        prev = chain
        write (*,'(1x,a,i3,a,a1,a)')
     +      'NMR model ',ichain - ichar('A') + 1,
     +      ' becomes chain |',chain,'|'
c
        jchain = jchain + 1
        chnptr (1,jchain) = natoms + 1
        chname (jchain) = chain
        goto 10
      end if
c
c ... read cell from CRYST1 card if present
c
      if (key .eq. 'CRYST1') then
        read (line,*) (cell(i),i=1,6)
        call fvalut (' Cell :',6,cell(1))
        lcell = .true.
        goto 10
      end if
c
c ... skip hetero atoms
c
      if (key .eq. 'HETATM') then
        nhet = nhet + 1
        goto 10
      end if
c
      if (key .ne. 'ATOM  ') goto 10
c
c ... XPLOR terminal oxygen
c
      if (line(7:10) .eq. ' OT1') line(7:10) = ' O  '
c
c ... skip hydrogens
c
      if (lhydro(line(7:10))) then
        nstrip = nstrip + 1
        goto 10
      end if
c
      do i=1,ntype
        if (line(12:14) .eq. typnam(i)) then
          do j=1,natype(i)
            if (line(7:10) .eq. typanm(j,i)) goto 20
          end do
        end if
      end do
      nunk = nunk + 1
      write (*,4711) ' UNK  ',line(1:leng1(line))
      goto 10
c
   20 continue
      natoms = natoms + 1
c
      if (natoms .gt. maxatm) then
        call errcon ('Too many atoms -- rest skipped')
        call jvalut (' Max nr of atoms =',1,maxatm)
        natoms = natoms - 1
        goto 21
      end if
c
      kk = natoms
      read (line,4713,err=5372) i,
     +  atmnam(kk),resnam(kk),achain(kk),
     +  iresid(kk),atmxyz(1,kk),atmxyz(2,kk),
     +  atmxyz(3,kk),qatom(kk),batom(kk)
c
      if ((.not. nmr) .and. achain(kk) .ne. prev) then
        ichain = ichain + 1
        if (ichain .gt. ichar('Z')) then
          ichain = ichain - 1
          call errcon ('Too many CHAINs; max = 26; rest skipped')
          goto 21
        end if
        chain = char(ichain)
        prev = achain (kk)
        write (*,'(1x,a,a1,a,a1,a)')
     +    'Old chain |',prev,'| becomes chain |',chain,'|'
        jchain = jchain + 1
        chnptr (1,jchain) = natoms
        chname (jchain) = chain
      end if
c
      achain (kk) = chain
c
      goto 10
c
c ... read from file error
c
  997 continue
      close (iunit)
      call errcon ('While reading PDB file')
      call jvalut (' Line nr :',1,nl)
      call jvalut (' Atom nr :',1,natoms)
      call textut (' Line :',line)
      ierr = -1
      return
c
c ... read from string error
c
 5372 continue
      close (iunit)
      call errcon ('While reading PDB file')
      call jvalut (' Line nr :',1,nl)
      call jvalut (' Atom nr :',1,natoms)
      call textut (' Line :',line)
      ierr = -1
      return
c
   21 continue
      write (*,*)
      ichain = ichain - ichar('A') + 1
      call jvalut (' Nr of lines read from file  :',1,nl)
      call jvalut (' Nr of atoms in molecule     :',1,natoms)
      call jvalut (' Nr of chains or models      :',1,ichain)
      call jvalut (' Nr of hydrogens stripped    :',1,nstrip)
      call jvalut (' Nr of HETATMs stripped      :',1,nhet)
      call jvalut (' Nr of non-AA atoms stripped :',1,nunk)
c
      chnptr (2,ichain) = natoms
      if (ichain .gt. 1) then
        do i=1,ichain-1
          chnptr (2,i) = chnptr (1,i+1) - 1
        end do
      end if
      nchain = ichain
c
      close (iunit)
      write (*,*)
c
      return
      end
c
c
c
      subroutine yasspa ()
c
      include 'seaman.incl'
c
      real alpha(3,5),beta(3,5),temp(3,5),rt(12),rmsd
c
      integer i,j,i1,i2,j1,j2,ierr
c
      logical wassec
c
      data alpha /1.456,0.389,3.482,    1.670,1.639,-0.126,
     +            1.840,-1.957,-1.468,  -1.271,-2.969,0.489,
     +            -3.183,0.105,-0.711/
c
      data beta /-0.737,3.837,-4.763,   0.927,1.784,-2.056,
     +           -0.843,-1.207,-0.492,  0.030,-2.870,2.832,
     +           -1.072,-6.256,4.133/
c
code ...
c
      call prompt (' Running YASSPA ...')
c
      do i=1,nres
        struct (i) = ' '
      end do
c
      do i=1,nres
        j1 = max (1,i-2)
        j2 = min (nres,i+2)
        if ( (j2-j1) .ne. 4 ) goto 10
        i1 = captr(j1)
        i2 = captr(j2)
        if (i1 .le. 0) goto 10
        if (i2 .le. 0) goto 10
        if (achain(i1) .ne. achain(i2)) goto 10
        if ( (iresid(i2)-iresid(i1)) .ne. 4 ) goto 10
c
c ... get CA coordinates
c
        do j=j1,j2
          i1=captr(j)
          if (i1 .le. 0) goto 10
          do i2=1,3
            temp(i2,j-j1+1) = atmxyz(i2,i1)
          end do
        end do
c
c ... do LSQ superpositioning
c
        call lsqgjk (alpha,temp,5,rmsd,rt,ierr)
c
ccc      call fvalut (' ALPHA :',1,rmsd)
c
        if (ierr .ne. 0) goto 10
        if (rmsd .le. alpyas) then
          struct (i-1) = 'ALPHA '
          struct (i)   = 'ALPHA '
          struct (i+1) = 'ALPHA '
          goto 10
        end if
c
        call lsqgjk (beta,temp,5,rmsd,rt,ierr)
c
ccc      call fvalut (' BETA  :',1,rmsd)
c
        if (ierr .ne. 0) goto 10
        if (rmsd .le. betyas) then
          struct (i-1) = 'BETA  '
          struct (i)   = 'BETA  '
          struct (i+1) = 'BETA  '
          goto 10
        end if
c
   10   continue
c
      end do
c
      struct (0) = ' '
      struct (nres+1) = ' '
c
      do i=1,nres
c
        if (struct(i).eq.' ') then
          if (struct(i-1).eq.'ALPHA ' .and.
     +        struct(i+1).eq.'ALPHA ') then
            struct (i) = 'ALPHA '
          else if (struct(i-1).eq.'BETA  ' .and.
     +             struct(i+1).eq.'BETA  ') then
            struct (i) = 'BETA  '
          end if
        end if
      end do
c
      do i=1,nres
        if (struct(i).ne.' ') then
          if (struct(i-1).eq.' ' .and. struct(i+1).eq.' ') then
            struct(i) = ' '
cc          else if (struct(i+1).eq.struct(i)) then
cc            if (struct(i-1).eq.' ' .and. struct(i+2).eq.' ') then
cc              struct(i) = ' '
cc              struct(i+1) = ' '
cc            end if
cc          else if (struct(i+1).ne.struct(i)) then
cc            struct(i) = ' '
cc            struct(i+1) = ' '
          end if
        end if
      end do
c
      i1 = 0
      i2 = 0
      j1 = 0
      j2 = 0
      wassec = .false.
      i = 0
c
   20 continue
      i = i + 1
      if (i .gt. nres) goto 30
c
      if (struct(i) .eq. 'ALPHA ') then
        i1 = i1 + 1
        wassec = .true.
        goto 20
      else if (struct(i) .eq. 'BETA  ') then
        i2 = i2 + 1
        wassec = .true.
        goto 20
      end if
c
      if (wassec) then
        wassec = .false.
        if ( i .lt. nres) then
          if (struct(i+1) .ne. ' ') then
            if (conect(i-1) .and. conect(i)) then
              struct(i) = 'TURN  '
              j1 = j1 + 1
              goto 20
            end if
          end if
        end if
c
        if ( (i+1) .lt. nres) then
          if (struct(i+2) .ne. ' ') then
            if (struct(i+1).eq.' ' .and.
     +          conect(i-1) .and. conect(i) .and.
     +          conect(i+1)) then
              struct(i) = 'TURN  '
              struct(i+1) = 'TURN  '
              i = i + 1
              j1 = j1 + 2
              goto 20
            end if
          end if
        end if
c
        if ( (i+2) .lt. nres) then
          if (struct(i+3) .ne. ' ') then
            if (struct(i+1).eq.' ' .and.
     +          struct(i+2).eq.' ' .and.
     +          conect(i-1) .and. conect(i) .and.
     +          conect(i+1) .and. conect(i+2)) then
              struct(i) = 'TURN  '
              struct(i+1) = 'TURN  '
              struct(i+2) = 'TURN  '
              i = i + 2
              j1 = j1 + 3
              goto 20
            end if
          end if
        end if
c
        if ( (i+3) .lt. nres) then
          if (struct(i+4) .ne. ' ') then
            if (struct(i+1).eq.' ' .and.
     +          struct(i+2).eq.' ' .and.
     +          struct(i+3).eq.' ' .and.
     +          conect(i-1) .and. conect(i) .and.
     +          conect(i+1) .and. conect(i+2) .and.
     +          conect(i+3)) then
              struct(i) = 'TURN  '
              struct(i+1) = 'TURN  '
              struct(i+2) = 'TURN  '
              struct(i+3) = 'TURN  '
              i = i + 3
              j1 = j1 + 4
              goto 20
            end if
          end if
        end if
      end if
c
      j2 = j2 + 1
      struct (i) = 'LOOP  '
      goto 20
c
   30 continue
      call jvalut (' Nr of ALPHA :',1,i1)
      call jvalut (' Nr of BETA  :',1,i2)
      call jvalut (' Nr of TURN  :',1,j1)
      call jvalut (' Nr of LOOP  :',1,j2)
c
      return
      end
c
c
c
      subroutine surfer ()
c
      include 'seaman.incl'
c
      real cog(3),rcog(3)
      real dist,x
c
      integer i,j,k,l,nn,ica,cosang
c
      logical mainch
c
code ...
c
      if (natoms .ge. 0) return
c
      call prompt (' Looking for surface residues ...')
c
      do i=1,nchain
        cog (1) = 0.0
        cog (2) = 0.0
        cog (3) = 0.0
        do j=chnptr(1,i),chnptr(2,i)
          do k=1,3
            cog(k)=cog(k)+atmxyz(k,j)
          end do
        end do
        do k=1,3
          cog(k)=cog(k)/float(chnptr(2,i)-chnptr(1,i)+1)
        end do
c
        do j=1,nres
c
c ... rsidue in chain ?
c
          if (resptr(1,j) .lt. chnptr(1,i)) goto 10
          if (resptr(2,j) .gt. chnptr(2,i)) goto 10
c
          surfac (j) = .false.
c
          ica= captr(j)
          if (ica .le. 0) goto 10
c
c ... get centre of gravity of the side chain atoms
c
          do k=1,3
            rcog(k) = 0.0
          end do
          do k=resptr(1,j),resptr(2,j)
            if (.not. mainch(atmnam(k))) then
              nn=nn+1
              do l=1,3
                rcog(l)=rcog(l)+atmxyz(l,k)
              end do
            end if
          end do
c
c .. if gly, skip
c
          if (nn .le. 0) goto 10
          do k=1,3
            rcog(k)=rcog(k)/float(nn)
          end do
c
c ... if sidechain points towards CoG of chain, assume buried
c
          cosang = 0.0
          do k=1,3
            cosang = cosang + (rcog(k)-atmxyz(k,ica)) *
     +               (cog(k)-atmxyz(k,ica))
          end do
          if (cosang .le. 0.0) goto 10
c
c ... if not, count neigbouring atoms
c
          nn = 0
          do k=resptr(1,j),resptr(2,j)
            if (.not. mainch(atmnam(k))) then
              do l=chnptr(1,i),chnptr(2,i)
                if (iresid(l) .ne. iresid(k)) then
                  x = dist(l,k,atmxyz)
                  if (x .le. nbrdst) nn=nn+1
                end if
              end do
              if (nn .ge. 3) goto 10
            end if
          end do
          if (nn .ge. 3) goto 10
c
          surfac (j) = .true.
c
   10     continue
c
        end do
c
      end do
c
      return
      end
c
c
c
      subroutine bookkp (ierr)
c
      include 'seaman.incl'
c
      real x
      real dist
c
      integer i,j,nowres,ierr,i1,i2
c
      character nownam*3,nowchn*1
c
code ...
c
      call prompt (' Doing some book-keeping ...')
c
      ierr = 0
      nres = 0
      nowres = -93723
      nownam = '?@#'
      nowchn = '?'
c
c      print *,natoms,' ATOMS'
c
c ... find residues
c
      do i=1,natoms
        if (iresid(i).ne.nowres .or. resnam(i).ne.nownam .or.
     +      achain(i).ne.nowchn) then
          nres = nres + 1
          if (nres .gt. maxres) then
            call errcon ('Too many residues')
            ierr = -1
            return
          end if
          resptr (1,nres) = i
          nowres = iresid(i)
          nownam = resnam(i)
          nowchn = achain(i)
          write (namres(nres),'(a1,i5)') nowchn,nowres
          call remspa (namres(nres))
          call upcase (namres(nres))
c
c      print *,nres,i,' |',namres(nres),'|',iresid(i)
c
        end if
      end do
c
c ... get residue pointers
c
      resptr (2,nres) = natoms
      do i=1,nres-1
        resptr(2,i) = resptr(1,i+1)-1
      end do
c
c ... get CA pointers
c
      do i=1,nres
        captr(i) = -1
        do j=resptr(1,i),resptr(2,i)
          if (atmnam(j) .eq. ' CA ') then
            captr(i) = j
            goto 10
          end if
        end do
   10   continue
      end do
c
c ... get connectivity
c
      conect (0) = .false.
      conect (nres) = .false.
      do i=1,nres-1
        i1=captr(i)
        i2=captr(i+1)
        conect (i) = .false.
        if (i1 .le. 0) goto 20
        if (i2 .le. 0) goto 20
        if (achain(i1) .ne. achain(i2)) goto 20
        if ((iresid(i1)+1) .ne. iresid(i2)) goto 20
        x = dist (i1,i2,atmxyz)
        if (x .gt. cacamx) goto 20
        conect (i) = .true.
   20   continue
      end do
c
c      do i=1,nres
c        write (*,6912) i,namres(i),resptr(1,i),resptr(2,i),
c     +    captr(i),conect(i)
c      end do
c
 6912 format (1x,i6,1x,a6,3(1x,i8),1x,l1)
c
      return
      end
c
c
c
      subroutine bfadel (blo,bhi)
c
      include 'seaman.incl'
c
      real blo,bhi
c
      integer i,nn
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      call rlohi (blo,bhi)
      do i=1,natoms
        delatm(i)=.false.
        if (batom(i).lt.blo .or. batom(i).gt.bhi) then
          delatm(i)=.true.
          nn=nn+1
        end if
      end do
c
      call deleta (nn,.true.)
c
      return
      end
c
c
c
      subroutine occdel (qlo,qhi)
c
      include 'seaman.incl'
c
      real qlo,qhi
c
      integer i,nn
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      call rlohi (qlo,qhi)
      do i=1,natoms
        delatm(i)=.false.
        if (qatom(i).lt.qlo .or. qatom(i).gt.qhi) then
          delatm(i)=.true.
          nn=nn+1
        end if
      end do
c
      call deleta (nn,.true.)
c
      return
      end
c
c
c
      subroutine zondel ()
c
      include 'seaman.incl'
c
      integer i,nn
c
code ...
c
      nn = 0
      do i=1,natoms
        delatm(i)=.false.
        if (select(i)) then
          delatm(i)=.true.
          nn=nn+1
        end if
      end do
c
      call deleta (nn,.true.)
c
      return
      end
c
c
c
      subroutine isodel ()
c
c ... delete atoms which are not connected to at least one
c     other atom in the residue
c
      include 'seaman.incl'
c
      real dist
c
      integer i,nn,j,k
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      do i=1,nres
        do j=resptr(1,i),resptr(2,i)
          do k=resptr(1,i),resptr(2,i)
            if (j.eq.k) goto 20
            if (dist(j,k,atmxyz).le.bondis) goto 10
   20       continue
          end do
          delatm(i)=.true.
          nn=nn+1
   10     continue
        end do
      end do          
c
      call deleta (nn,.true.)
c
      return
      end
c
c
c
      subroutine loodel ()
c
      include 'seaman.incl'
c
      integer i,nn,j
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      do i=1,nres
        if (struct(i) .eq. 'LOOP  ') then
          do j=resptr(1,i),resptr(2,i)
            delatm(j)=.true.
            nn=nn+1
          end do
        end if          
      end do          
c
      call deleta (nn,.true.)
c
      return
      end
c
c
c
      subroutine turdel ()
c
      include 'seaman.incl'
c
      integer i,nn,j
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      do i=1,nres
        if (struct(i) .eq. 'TURN  ') then
          do j=resptr(1,i),resptr(2,i)
            delatm(j)=.true.
            nn=nn+1
          end do
        end if          
      end do          
c
      call deleta (nn,.true.)
c
      return
      end
c
c
c
      subroutine deleta (nn,askem)
c
      include 'seaman.incl'
c
      integer i,j,nn,na,k,ierr
c
      logical xinter,isitok,askem
c
      character now*1
c
code ...
c
      do i=1,natoms
        if (.not. delatm(i)) then
          if (atmnam(i) .eq. 'XXXX') then
            delatm(i) = .true.
            nn = nn + 1
          end if
        end if
      end do
c
      if (nn .le. 0) then
        call prompt (' NO atoms marked for deletion')
        return
      end if
c
      if (askem) then
        if (xinter()) then
          call jvalut (' Nr of atoms marked for deletion :',1,nn)
          if (.not. isitok(' Delete them ?')) return
        end if
      else
        call jvalut (' Nr of atoms marked for deletion :',1,nn)
      end if
c
      j=0
      do i=1,nres
        na = 0
        do k=resptr(1,i),resptr(2,i)
          if (.not. delatm(k)) na = na + 1
        end do
        if (na .gt. 0) then
          j = j + 1
          if (i.ne.j) then
            struct(j)=struct(i)
            surfac(j)=surfac(i)
          end if
        end if
      end do
c
      j=0
      do i=1,natoms
        if (delatm(i)) goto 10
        j = j + 1
        if (j .ne. i) then
          atmxyz(1,j) = atmxyz(1,i)
          atmxyz(2,j) = atmxyz(2,i)
          atmxyz(3,j) = atmxyz(3,i)
          qatom(j) = qatom(i)
          batom(j) = batom(i)
          iresid(j) = iresid(i)
          atmnam(j) = atmnam(i)
          resnam(j) = resnam(i)
          achain(j) = achain(i)
        end if
   10   continue
      end do
c
      natoms = j
      call jvalut (' Nr of atoms now :',1,natoms)
c
      nchain = 0
      now = '?'
      do i=1,natoms
        if (achain(i) .ne. now) then
          nchain = nchain + 1
          chname(nchain)=achain(i)
          now = achain(i)
          chnptr(1,nchain)=i
        end if
      end do
c
      chnptr(2,nchain)=natoms
      do i=1,nchain-1
        chnptr(2,i)=chnptr(1,i+1)-1
      end do
c
      call bookkp (ierr)
c
      return
      end
c
c
c
      subroutine selecm (zone,chain,nn)
c
      include 'seaman.incl'
c
      integer zone(2),i,j,nn
c
      character chain*1
c
code ...
c
      call prompt (' Select zone (0 0 = all residues)')
      call jvalin (' Zone ?',2,zone)
      if (zone(1).ne.0 .and. zone(2).ne.0) then
        call ilohi (zone(1),zone(2))
      end if
c
      if (nchain .gt. 1) then
        call prompt (' Select chain (* = all chains)')
        call textin (' Chain ?',chain)
        call upcase (chain)
      else
        chain = chname(1)
      end if
c
      if (chain .eq. '*') then
        if (zone(1).eq.0 .and. zone(2).eq.0) then
          do i=1,natoms
            select(i)=.true.
          end do
          nn = natoms
          goto 20
        else if (zone(1).eq.0) then
          do i=1,natoms
            select(i) = (iresid(i).le.zone(2))
          end do
        else if (zone(2).eq.0) then
          do i=1,natoms
            select(i) = (iresid(i).ge.zone(1))
          end do
        else
          do i=1,natoms
            select(i) = (iresid(i).ge.zone(1) .and.
     +                   iresid(i).le.zone(2))
          end do
        end if
c
      else
        do i=1,natoms
          select(i)=.false.
        end do
        do j=1,nchain
          if (chname(j).eq.chain) then
            if (zone(1).eq.0 .and. zone(2).eq.0) then
              do i=chnptr(1,j),chnptr(2,j)
                select(i)=.true.
              end do
              nn = chnptr(2,j)-chnptr(1,j)+1
              goto 20
            else if (zone(1).eq.0) then
              do i=chnptr(1,j),chnptr(2,j)
                select(i) = (iresid(i).le.zone(2))
              end do
            else if (zone(2).eq.0) then
              do i=chnptr(1,j),chnptr(2,j)
                select(i) = (iresid(i).ge.zone(1))
              end do
            else
              do i=chnptr(1,j),chnptr(2,j)
                select(i) = (iresid(i).ge.zone(1) .and.
     +                       iresid(i).le.zone(2))
              end do
            end if
            goto 10
          end if
        end do
      end if
c
   10 continue
      nn = 0
      do i=1,natoms
        if (select(i)) nn=nn+1
      end do
c
   20 continue
      call jvalut (' Nr of selected atoms :',1,nn)
      if (nn .le. 0) call errcon ('No atoms selected')
c
      return
      end
c
c
c
      subroutine polgly ()
c
      include 'seaman.incl'
c
      integer i,nn
c
      logical mainch
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      do i=1,natoms
        if (select(i)) then
          if (resnam(i) .ne. 'GLY') then
            if (.not. mainch(atmnam(i))) then
              delatm(i) = .true.
              nn = nn + 1
            else
              resnam(i) = 'GLY'
            end if
          end if
        end if
      end do
c
      call deleta (nn,.false.)
c
      return
      end
c
c
c
      subroutine polala ()
c
      include 'seaman.incl'
c
      integer i,nn
c
      logical mainch
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      do i=1,natoms
        if (select(i)) then
          if (resnam(i) .ne. 'GLY' .and.
     +        resnam(i) .ne. 'ALA') then
            if ( (.not. mainch(atmnam(i))) .and.
     +           (atmnam(i).ne.' CB ') ) then
              delatm(i) = .true.
              nn = nn + 1
            else
              resnam(i) = 'ALA'
            end if
          end if
        end if
      end do
c
      call deleta (nn,.false.)
c
      return
      end
c
c
c
      subroutine polser ()
c
      include 'seaman.incl'
c
      integer i,nn
c
      logical mainch
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
      do i=1,natoms
        if (select(i)) then
          if (resnam(i) .ne. 'GLY' .and.
     +        resnam(i) .ne. 'ALA' .and.
     +        resnam(i) .ne. 'SER') then
            resnam (i) = 'SER'
            if (mainch(atmnam(i))) goto 10
            if (atmnam(i) .eq. ' CB ') goto 10
            if (atmnam(i) .eq. ' CG ' .or.
     +          atmnam(i) .eq. ' OG ' .or.
     +          atmnam(i) .eq. ' SG ' .or.
     +          atmnam(i) .eq. ' CG1' .or.
     +          atmnam(i) .eq. ' OG1') then
              atmnam (i) = ' OG '
              goto 10
            end if
c
            delatm(i) = .true.
            nn = nn + 1
c
   10       continue
          end if
        end if
      end do
c
      call deleta (nn,.false.)
c
      return
      end
c
c
c
      subroutine selecr (type,nn)
c
      include 'seaman.incl'
c
      integer i,nn
c
      character type*3
c
code ...
c
      call textin (' Residue type ?',type)
      call upcase (type)
c
      do i=1,natoms
        select (i) = (resnam(i) .eq. type)
      end do
c
      nn = 0
      do i=1,natoms
        if (select(i)) nn=nn+1
      end do
c
      call jvalut (' Nr of selected atoms :',1,nn)
      if (nn .le. 0) call errcon ('No atoms selected')
c
      return
      end
c
c
c
      subroutine selec1 (zone,chain,nn)
c
      include 'seaman.incl'
c
      integer i,j,nn,zone
c
      character chain*1
c
code ...
c
      nn = 0
      call jvalin (' Residue number ?',1,zone)
      if (zone.le.0) then
        call errcon ('Invalid residue number')
        return
      end if
c
      if (nchain .gt. 1) then
        call prompt (' Select chain (* = all chains)')
        call textin (' Chain ?',chain)
        call upcase (chain)
      else
        chain = chname(1)
      end if
c
      if (chain .eq. '*') then
        do i=1,nres
          j = resptr(1,i)
          selres (i) = (iresid(j).eq.zone)
        end do
      else
        do i=1,nres
          selres(i)=.false.
        end do
        do i=1,nres
          j = resptr(1,i)
          selres (i) = (iresid(j).eq.zone .and.
     +                  achain(j).eq.chain)
        end do
      end if
c
   10 continue
      nn = 0
      do i=1,nres
        if (selres(i)) nn=nn+1
      end do
c
   20 continue
      call jvalut (' Nr of selected residues :',1,nn)
      if (nn .le. 0) call errcon ('No residues selected')
c
      return
      end
c
c
c
      subroutine minima (ires,newtyp)
c
      include 'seaman.incl'
c
      integer i,nn,ires,ityp,jtyp,n1,n2
c
      character newtyp*3,type*3
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
c
      if (newtyp .eq. 'GLY' .or. newtyp .eq. 'ALA' .or.
     +    newtyp .eq. 'SER' .or. newtyp .eq. 'PRO' .or.
     +    newtyp .eq. 'TRP' .or. newtyp .eq. 'ILE' .or.
     +    newtyp .eq. 'MET' .or. newtyp .eq. 'ARG' .or.
     +    newtyp .eq. 'LEU' .or. newtyp .eq. 'TYR') then
        call errcon (' Cannot change to G/A/S/P/W/I/L/Y/M/R')
        return
      end if
c
      ityp = 0
      jtyp = 0
      type = resnam(resptr(1,ires))
      do i=1,ntype
        if (type.eq.typnam(i)) ityp=i
        if (newtyp.eq.typnam(i)) jtyp=i
      end do
c
      if (ityp .le. 0 .or. jtyp .le. 0) then
        call errcon (' Residue type not in library')
        return
      else if (ityp .eq. jtyp) then
        call errcon (' Identical residue types; no change')
        return
      end if
c
      n1 = resptr(2,ires)-resptr(1,ires)+1
      n2 = natype(jtyp)
c
      if (n2 .gt. n1) then
        call prompt (' Sorry - not a minimalist substitution')
        return
      end if
c
      if (newtyp .eq. 'CYS') then
        call tocys (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'THR') then
        call tothr (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'VAL') then
        call toval (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'LYS') then
        call tolys (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'ASN') then
        call toasx (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'ASP') then
        call toasx (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'GLU') then
        call toglx (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'GLN') then
        call toglx (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'HIS') then
        call tohis (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else if (newtyp .eq. 'PHE') then
        call tophe (ires,nn)
        if (nn .gt. 0) call deleta (nn,.false.)
      else
        call prompt (' Sorry - substitution not available')
      end if
c
      return
      end
c
c
c
      subroutine insert (ires,nnew,ierr)
c
      include 'seaman.incl'
c
      integer ires,nnew,ierr,i,j,jlo,jhi
c
code ...
c
      ierr = -1
      if (nnew .le. 0) then
        call prompt (' No atoms to insert')
        return
      end if
c
      ierr = -2
      jhi = natoms + nnew
      if (jhi .gt. maxatm) then
        call errcon (' No room for more atoms')
        return
      end if
c
      ierr = 0
      if (ires .eq. nres) then
        natoms = natoms + nnew
        i = resptr(1,ires)
        do j=resptr(2,ires)+1,natoms
          atmxyz(1,j) = 0.0
          atmxyz(2,j) = 0.0
          atmxyz(3,j) = 0.0
          qatom(j) = defq
          batom(j) = defb
          iresid(j) = iresid(i)
          atmnam(j) = 'XXXX'
          resnam(j) = resnam(i)
          achain(j) = achain(i)
        end do
        resptr (2,ires) = natoms
        return
      end if
c
      jlo = resptr(1,ires+1) + nnew
      jhi = natoms + nnew
      if (jhi .gt. maxatm) then
        ierr = -3
        call errcon (' No room for more atoms')
        return
      end if
c
      do i=ires+1,nres
        resptr(1,i)=resptr(1,i)+nnew
        resptr(2,i)=resptr(2,i)+nnew
      end do
c
cc      print *,resptr(1,ires),resptr(2,ires),nnew
cc      print *,jhi,jlo
c
      do j=jhi,jlo,-1
        i = j - nnew
        atmxyz(1,j) = atmxyz(1,i)
        atmxyz(2,j) = atmxyz(2,i)
        atmxyz(3,j) = atmxyz(3,i)
        qatom(j) = qatom(i)
        batom(j) = batom(i)
        iresid(j) = iresid(i)
        atmnam(j) = atmnam(i)
        resnam(j) = resnam(i)
        achain(j) = achain(i)
      end do
c
      i = resptr(1,ires)
      do j=resptr(2,ires)+1,resptr(2,ires)+nnew
        atmxyz(1,j) = 0.0
        atmxyz(2,j) = 0.0
        atmxyz(3,j) = 0.0
        qatom(j) = defq
        batom(j) = defb
        iresid(j) = iresid(i)
        atmnam(j) = 'XXXX'
        resnam(j) = resnam(i)
        achain(j) = achain(i)
      end do
      resptr (2,ires) = resptr(2,ires)+nnew
c
      natoms = jhi
c
      return
      end
c
c
c
      subroutine mutate (ires,newtyp)
c
      include 'seaman.incl'
c
      real myxyz(3,4),rt(12)
      real rmsd
c
      integer i,nn,ierr,ires,ityp,jtyp,n1,n2,j,nuse
c
      logical allok(3)
      logical mainch
c
      character newtyp*3,type*3
c
code ...
c
      do i=1,natoms
        delatm(i) = .false.
      end do
      nn = 0
c
      if (newtyp .eq. 'GLY') then
        call errcon (' Cannot change to Gly')
        return
      end if
c
      ityp = 0
      jtyp = 0
      type = resnam(resptr(1,ires))
ccc      print *,type,' ',newtyp
      do i=1,ntype
        if (type.eq.typnam(i)) ityp=i
        if (newtyp.eq.typnam(i)) jtyp=i
ccc        print *,typnam(i)
      end do
c
      if (ityp .le. 0 .or. jtyp .le. 0) then
        call errcon (' Residue type not in library')
        return
      else if (ityp .eq. jtyp) then
        call prompt (' Identical residue types; replace by rotamer')
      end if
c
      n1 = resptr(2,ires)-resptr(1,ires)+1
      n2 = natype(jtyp)
c
      ierr = 0
      if (n2 .gt. n1) call insert (ires,(n2-n1),ierr)
      if (ierr .ne. 0) then
        call prompt (' Unable to insert atoms - oops')
        return
      end if
c
      do i=1,3
        allok(i)=.false.
      end do
c
      do i=resptr(1,ires),resptr(2,ires)
        if (atmnam(i) .eq. ' N  ') then
          j=1
        else if (atmnam(i) .eq. ' CA ') then
          j=2
        else if (atmnam(i) .eq. ' C  ') then
          j=3
        else
          goto 10
        end if
        allok(j)=.true.
        myxyz(1,j)=atmxyz(1,i)
        myxyz(2,j)=atmxyz(2,i)
        myxyz(3,j)=atmxyz(3,i)
   10   continue
      end do
c
      if (.not. (allok(1).and.allok(2).and.allok(3))) then
        call errcon (' Missing main-chain atoms')
        return
      end if
      nuse = 3
c
cc      print *,'Do LSQ'
c
      ierr = 0
      call lsqgjk (myxyz,typxyz(1,1,jtyp),nuse,rmsd,rt,ierr)
      call fvalut (' RMSD (A) :',1,rmsd)
      if (ierr .ne. 0) then
        call prompt (' Unable to do LSQ - oops')
        return
      end if
c
      call vecrtv (typxyz(1,1,jtyp),typxyz(1,1,jtyp),
     +             natype(jtyp),rt(1),rt(10))
c
      do i=resptr(1,ires),resptr(1,ires)+n2-1
        resnam (i) = newtyp
      end do
c
      j = 0
      do i=resptr(1,ires),resptr(1,ires)+n2-1
        if (.not. mainch(atmnam(i))) then
   19     continue
          j = j + 1
          if (mainch(typanm(j,jtyp))) goto 19
          qatom (i) = defq
          batom (i) = defb
          atmnam (i) = typanm(j,jtyp)
          atmxyz (1,i) = typxyz(1,j,jtyp)
          atmxyz (2,i) = typxyz(2,j,jtyp)
          atmxyz (3,i) = typxyz(3,j,jtyp)
        end if
      end do
c
      nn = 0
      if (n1 .gt. n2) then
        do i=resptr(1,ires)+n2,resptr(1,ires)+n1-1
          delatm (i) = .true.
          nn = nn + 1
        end do
        call deleta (nn,.false.)
c
      else
        call deleta (nn,.false.)
      end if
c
ccc      call bookkp (ierr)
c
      return
      end
c
c
c
      subroutine tocys (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).eq.'GLY' .or. resnam(i1).eq.'ALA') then
        call prompt (' Sorry, only [~GLY|~ALA] -> CYS')
        return
      end if
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = 'CYS'
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG ' .or.
     +      atmnam(i) .eq. ' OG ' .or.
     +      atmnam(i) .eq. ' SG ' .or.
     +      atmnam(i) .eq. ' CG1' .or.
     +      atmnam(i) .eq. ' OG1') then
          atmnam (i) = ' SG '
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine tothr (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'ILE' .and. resnam(i1).ne.'VAL') then
        call prompt (' Sorry, only [ILE|VAL] -> THR')
        return
      end if
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = 'THR'
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG2') goto 10
        if (atmnam(i) .eq. ' CG1') then
          atmnam (i) = ' OG1'
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine toval (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'ILE' .and. resnam(i1).ne.'THR') then
        call prompt (' Sorry, only [ILE|THR] -> VAL')
        return
      end if
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = 'VAL'
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG1') goto 10
        if (atmnam(i) .eq. ' CG2') goto 10
        if (atmnam(i) .eq. ' OG1') then
          atmnam (i) = ' CG1'
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine tolys (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'ARG') then
        call prompt (' Sorry, only [ARG] -> LYS')
        return
      end if
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = 'LYS'
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG ') goto 10
        if (atmnam(i) .eq. ' CD ') goto 10
        if (atmnam(i) .eq. ' NE ') then
          atmnam (i) = ' CE '
          goto 10
        end if
        if (atmnam(i) .eq. ' CZ ') then
          atmnam (i) = ' NZ '
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine toasx (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
      character new*3
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'ASP' .and. resnam(i1).ne.'ASN') then
        call prompt (' Sorry, only ASP <-> ASN')
        return
      end if
c
      new = 'ASP'
      if (resnam(i1).eq.'ASP') new='ASN'
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = new
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG ') goto 10
        if (atmnam(i) .eq. ' OD1') goto 10
        if (atmnam(i) .eq. ' OD2') then
          atmnam (i) = ' ND2'
          goto 10
        end if
        if (atmnam(i) .eq. ' ND2') then
          atmnam (i) = ' OD2'
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine toglx (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
      character new*3
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'GLU' .and. resnam(i1).ne.'GLN') then
        call prompt (' Sorry, only GLU <-> GLN')
        return
      end if
c
      new = 'GLU'
      if (resnam(i1).eq.'GLU') new='GLN'
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = new
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG ') goto 10
        if (atmnam(i) .eq. ' CD ') goto 10
        if (atmnam(i) .eq. ' OE1') goto 10
        if (atmnam(i) .eq. ' OE2') then
          atmnam (i) = ' NE2'
          goto 10
        end if
        if (atmnam(i) .eq. ' NE2') then
          atmnam (i) = ' OE2'
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine tohis (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'TRP') then
        call prompt (' Sorry, only TRP -> HIS')
        return
      end if
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = 'HIS'
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .eq. ' CB ') goto 10
        if (atmnam(i) .eq. ' CG ') goto 10
        if (atmnam(i) .eq. ' CD2') goto 10
        if (atmnam(i) .eq. ' CD1') then
          atmnam (i) = ' ND1'
          goto 10
        end if
        if (atmnam(i) .eq. ' NE1') then
          atmnam (i) = ' CE1'
          goto 10
        end if
        if (atmnam(i) .eq. ' CE2') then
          atmnam (i) = ' NE2'
          goto 10
        end if
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine tophe (ires,nn)
c
      include 'seaman.incl'
c
      integer i,nn,ires,i1
c
      logical mainch
c
code ...
c
      nn = 0
      i1=resptr(1,ires)
      if (resnam(i1).ne.'TYR') then
        call prompt (' Sorry, only TYR -> PHE')
        return
      end if
c
      do i=resptr(1,ires),resptr(2,ires)
        resnam (i) = 'PHE'
        if (mainch(atmnam(i))) goto 10
        if (atmnam(i) .ne. ' OH ') goto 10
c
        delatm(i) = .true.
        nn = nn + 1
c
   10   continue
      end do
c
      return
      end
c
c
c
      subroutine psea (factor)
c
      include 'seaman.incl'
c
      real alpha(2,5),beta(2,5),temp(3,5)
      real a,t,d1,d2,d3,dist,angle,tangle,factor,f,t1
c
      integer i,j,i1,i2,j1,j2,ierr
c
      logical crital(5),critbe(5)
      logical lalpha(maxres),lalext(maxres)
      logical lbeta(maxres),lbeext(maxres)
      logical wassec
c
      data alpha /92.0, 15.0,   50.0, 25.0,
     +             5.5,  0.5,    5.15, 0.4,
     +             6.1,  0.35/
c
      data beta  /125.0, 20.0,   -165.0, 55.0,
     +              6.75, 0.75,     9.9,  1.0,
     +             12.4,  1.1/
c
code ...
c
      call prompt (' Running P-SEA ...')
      call prompt (' Labesse et al., CABIOS 13, 291 (1997)')
      f = max (0.1, min (10.0, factor))
      call fvalut (' Sensitivity :',1,f)
c
      do i=1,nres
c
        do j=1,5
          crital(j) = .false.
          critbe(j) = .false.
        end do
        struct (i) = ' '
        j1 = max (1,i-1)
        j2 = min (nres,i+3)
        if ( (j2-j1) .ne. 4 ) goto 10
        i1 = captr(j1)
        i2 = captr(j2)
        if (i1 .le. 0) goto 10
        if (i2 .le. 0) goto 10
        if (achain(i1) .ne. achain(i2)) goto 10
        if ( (iresid(i2)-iresid(i1)) .ne. 4 ) goto 10
c
c ... get CA coordinates
c
        do j=j1,j2
          i1=captr(j)
          if (i1 .le. 0) goto 10
          do i2=1,3
            temp(i2,j-j1+1) = atmxyz(i2,i1)
          end do
        end do
c
        a = angle (1,2,3,temp)
        if (abs(a-alpha(1,1)) .le. f*alpha(2,1))
     +    crital(1) = .true.
        if (abs(a-beta(1,1)) .le. f*beta(2,1))
     +    critbe(1) = .true.
c
        t = tangle (1,2,3,4,temp)
        call fixdif (t,alpha(1,2),t1)
        if (abs(t1) .le. f*alpha(2,2))
     +    crital(2) = .true.
        call fixdif (t,beta(1,2),t1)
        if (abs(t1) .le. f*beta(2,2))
     +    critbe(2) = .true.
c
        d1 = dist (1,3,temp)
        if (abs(d1-alpha(1,3)) .le. f*alpha(2,3))
     +    crital(3) = .true.
        if (abs(d1-beta(1,3)) .le. f*beta(2,3))
     +    critbe(3) = .true.
c
        d2 = dist (1,4,temp)
        if (abs(d2-alpha(1,4)) .le. f*alpha(2,4))
     +    crital(4) = .true.
        if (abs(d2-beta(1,4)) .le. f*beta(2,4))
     +    critbe(4) = .true.
c
        d3 = dist (1,5,temp)
        if (abs(d3-alpha(1,5)) .le. f*alpha(2,5))
     +    crital(5) = .true.
        if (abs(d3-beta(1,5)) .le. f*beta(2,5))
     +    critbe(5) = .true.
c
        lalpha (i) = ( (crital(4) .and. crital(5)) .or.
     +                 (crital(1) .and. crital(2)) )
        lalext (i) = ( crital(4) .or. crital(1) .or. lalpha(i))
c
        lbeta (i) = ( (critbe(3).and.critbe(4).and.critbe(5)) .or.
     +                (critbe(1) .and. critbe(2)) )
        lbeext (i) = (critbe(4) .or. lbeta(i))
c
c        write (*,6010) iresid(captr(i)),a,t,d1,d2,d3
c        print *,' AB ',lalpha(i),lalext(i),lbeta(i),lbeext(i)
c 6010 format (' Residue ',i4,' Angles ',2f8.1,' D123 ',3f6.1)
c
   10   continue
c
      end do
c
c ... find helices
c
  534 continue
      j = 0
      do i=1,nres
        if (i .le. (nres-4)) then
        if (lalpha(i) .and. lalpha(i+1) .and. lalpha(i+2) .and.
     +      lalpha(i+3) .and. lalpha(i+4) .and.
     +      struct(i) .eq. ' ' .and. struct(i+1) .eq. ' ' .and.
     +      struct(i+2) .eq. ' ' .and. struct(i+3) .eq. ' ' .and.
     +      struct(i+4) .eq. ' ' ) then
          struct (i)   = 'ALPHA '
          struct (i+1) = 'ALPHA '
          struct (i+2) = 'ALPHA '
          struct (i+3) = 'ALPHA '
          struct (i+4) = 'ALPHA '
          j = j + 5
        end if
        end if
c        if (lalpha(i) .and. lalpha(i+1) .and. lalext(i+2) .and.
c     +      lalpha(i+3) .and. lalpha(i+4) .and.
c     +      struct(i) .eq. ' ' .and. struct(i+1) .eq. ' ' .and.
c     +      struct(i+2) .eq. ' ' .and. struct(i+3) .eq. ' ' .and.
c     +      struct(i+4) .eq. ' ' ) then
c          struct (i)   = 'ALPHA '
c          struct (i+1) = 'ALPHA '
c          struct (i+2) = 'ALPHA '
c          struct (i+3) = 'ALPHA '
c          struct (i+4) = 'ALPHA '
c          j = j + 5
c        end if
        if (struct(i) .eq. ' ' .and. struct(i+1).eq. 'ALPHA '
     +      .and. lalext(i)) then
          struct (i) = 'ALPHA '
          j = j + 1
        end if
        if (struct(i) .eq. ' ' .and. struct(i-1).eq. 'ALPHA '
     +      .and. lalext(i)) then
          struct (i) = 'ALPHA '
          j = j + 1
        end if
      end do
      if (j .gt. 0) goto 534
c
c ... find strands
c
  536 continue
      j = 0
      do i=1,nres
        if (i .le. (nres-2)) then
        if (lbeta(i) .and. lbeta(i+1) .and. lbeta(i+2) .and.
     +      struct(i) .eq. ' ' .and. struct(i+1) .eq. ' ' .and.
     +      struct(i+2) .eq. ' ' ) then
          struct (i)   = 'BETA  '
          struct (i+1) = 'BETA  '
          struct (i+2) = 'BETA  '
          j = j + 3
        end if
        end if
        if (struct(i) .eq. ' ' .and. struct(i+1).eq. 'BETA  '
     +      .and. lbeext(i)) then
          struct (i) = 'BETA  '
          j = j + 1
        end if
        if (struct(i) .eq. ' ' .and. struct(i-1).eq. 'BETA  '
     +      .and. lbeext(i)) then
          struct (i) = 'BETA  '
          j = j + 1
        end if
      end do
      if (j .gt. 0) goto 536
c
ccc      if (j .le. 0) goto 30
c
      struct (0) = ' '
      struct (nres+1) = ' '
      struct (nres+2) = ' '
      do i=1,nres
c
        if (struct(i).eq.' ') then
          if (struct(i-1).eq.'ALPHA ' .and.
     +        struct(i+1).eq.'ALPHA ') then
            struct (i) = 'ALPHA '
          else if (struct(i-1).eq.'BETA  ' .and.
     +             struct(i+1).eq.'BETA  ') then
            struct (i) = 'BETA  '
          end if
        end if
        if (struct(i) .eq. 'ALPHA ') then
          if (struct(i-1) .ne. 'ALPHA ') then
            if (struct(i+1) .ne. 'ALPHA ' .or.
     +          struct(i+2) .ne. 'ALPHA ') then
              struct (i) = ' '
            end if
          end if
        end if
        if (struct(i) .eq. 'BETA  ') then
          if (struct(i-1) .ne. 'BETA  ') then
            if (struct(i+1) .ne. 'BETA  ' .or.
     +          struct(i+2) .ne. 'BETA  ') then
              struct (i) = ' '
            end if
          end if
        end if
      end do
c
      do i=1,nres
        if (struct(i).ne.' ') then
          if (struct(i-1).eq.' ' .and. struct(i+1).eq.' ') then
            struct(i) = ' '
cc          else if (struct(i+1).eq.struct(i)) then
cc            if (struct(i-1).eq.' ' .and. struct(i+2).eq.' ') then
cc              struct(i) = ' '
cc              struct(i+1) = ' '
cc            end if
cc          else if (struct(i+1).ne.struct(i)) then
cc            struct(i) = ' '
cc            struct(i+1) = ' '
          end if
        end if
      end do
c
      i1 = 0
      i2 = 0
      j1 = 0
      j2 = 0
      wassec = .false.
      i = 0
c
   20 continue
      i = i + 1
      if (i .gt. nres) goto 30
c
      if (struct(i) .eq. 'ALPHA ') then
        i1 = i1 + 1
        wassec = .true.
        goto 20
      else if (struct(i) .eq. 'BETA  ') then
        i2 = i2 + 1
        wassec = .true.
        goto 20
      end if
c
      if (wassec) then
        wassec = .false.
        if ( i .lt. nres) then
          if (struct(i+1) .ne. ' ') then
            if (conect(i-1) .and. conect(i)) then
              struct(i) = 'TURN  '
              j1 = j1 + 1
              goto 20
            end if
          end if
        end if
c
        if ( (i+1) .lt. nres) then
          if (struct(i+2) .ne. ' ') then
            if (struct(i+1).eq.' ' .and.
     +          conect(i-1) .and. conect(i) .and.
     +          conect(i+1)) then
              struct(i) = 'TURN  '
              struct(i+1) = 'TURN  '
              i = i + 1
              j1 = j1 + 2
              goto 20
            end if
          end if
        end if
c
        if ( (i+2) .lt. nres) then
          if (struct(i+3) .ne. ' ') then
            if (struct(i+1).eq.' ' .and.
     +          struct(i+2).eq.' ' .and.
     +          conect(i-1) .and. conect(i) .and.
     +          conect(i+1) .and. conect(i+2)) then
              struct(i) = 'TURN  '
              struct(i+1) = 'TURN  '
              struct(i+2) = 'TURN  '
              i = i + 2
              j1 = j1 + 3
              goto 20
            end if
          end if
        end if
c
        if ( (i+3) .lt. nres) then
          if (struct(i+4) .ne. ' ') then
            if (struct(i+1).eq.' ' .and.
     +          struct(i+2).eq.' ' .and.
     +          struct(i+3).eq.' ' .and.
     +          conect(i-1) .and. conect(i) .and.
     +          conect(i+1) .and. conect(i+2) .and.
     +          conect(i+3)) then
              struct(i) = 'TURN  '
              struct(i+1) = 'TURN  '
              struct(i+2) = 'TURN  '
              struct(i+3) = 'TURN  '
              i = i + 3
              j1 = j1 + 4
              goto 20
            end if
          end if
        end if
      end if
c
      j2 = j2 + 1
      struct (i) = 'LOOP  '
      goto 20
c
   30 continue
      call jvalut (' Nr of ALPHA :',1,i1)
      call jvalut (' Nr of BETA  :',1,i2)
      call jvalut (' Nr of TURN  :',1,j1)
      call jvalut (' Nr of LOOP  :',1,j2)
c
      return
      end
