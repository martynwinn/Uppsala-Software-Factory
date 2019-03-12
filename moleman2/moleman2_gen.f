c
c ... mole2_subs.f - SPECIFIC subroutines for MOLEMAN2
c
c ... i.e., subroutines that *DO* include 'moleman2.incl'
c
c
      subroutine readdb (iunit,ierr)
c
      include 'moleman2.incl'
c
      integer maxopt
      parameter (maxopt=50)
c
      integer nline,iunit,ierr,ii,jj,i,j,k,nopt,kk,length
c
      logical lmore,lmch
c
      character line*256,optpar(maxopt)*80,quote*1
c
      parameter (quote='''')
c
code ...
c
      nline  = 0
      nalias = 0
      nmrtyp = 0
      nmatyp = 0
      nlring = 0
      ierr = 0
c
      libcat(1)='PROT'
      libcat(2)='NUCL'
      libcat(3)='WATE'
      libcat(4)='META'
      libcat(5)='INOR'
      libcat(6)='CARB'
      libcat(7)='ORGA'
      libcat(8)='HETE'
c
 6000 format (A,A,A,A,A)
c
   10 continue
      read (iunit,6000,err=8000,end=9000) line
      nline = nline + 1
      if (line(1:1) .eq. '!') goto 10
      call upcase (line)
      if (line(1:3) .eq. 'REM') then
        call textut (' >',line(4:))
        goto 10
      end if
      if (line(1:3) .eq. 'RES') goto 20
      if (line(1:3) .eq. 'RNG') goto 120
      call errcon ('Invalid keyword (expected !, REM, RES, RNG)')
      call jvalut (' Line number    :',1,nline)
      call textut (' Offending line :',line)
      ierr = -1
      return
c
c ... new residue type
c
   20 continue
      if (nmrtyp .ge. mxrtyp) then
        call errcon ('Too many residue definitions; skipping rest')
        call jvalut (' Max allowed    :',1,mxrtyp)
        call jvalut (' Line number    :',1,nline)
        call textut (' Offending line :',line)
        goto 9000
      end if
      nmrtyp = nmrtyp + 1
      ii = nmrtyp
      nmrptr (1,ii) = -1
      nmrptr (2,ii) = -1
      nmaptr (1,ii) = -1
      nmaptr (2,ii) = -1
      lrname (ii) = line(5:7)
      lrdesc (ii) = line(8:)
      call pretty (lrdesc(ii))
      lrtype (ii) = ihete
      libolc (ii) = '?'
      jj = nmatyp
      kk = nalias
c
   30 continue
      read (iunit,6000,err=8000,end=9000) line
      nline = nline + 1
      if (line(1:1) .eq. '!') goto 30
      call upcase (line)
      if (line(1:3) .eq. 'REM') then
        call textut (' >',line(4:))
        goto 30
      end if
c
      if (line(1:3) .eq. 'TYP') then
        do i=1,nrlcat
          if (line(5:8) .eq. libcat(i)) lrtype (ii) = i
        end do
c
      else if (line(1:3) .eq. 'OLC') then
        libolc (ii) = line(5:5)
c
      else if (line(1:3) .eq. 'AKA') then
        call extrop (line(5:),nopt,maxopt,optpar,ierr)
        if (ierr .ne. 0) then
          call errcon ('While parsing AKA line')
          call jvalut (' Line number    :',1,nline)
          call textut (' Offending line :',line)
          ierr = -1
          return
        end if
        if ( (nalias+nopt) .gt. mxrtal) then
          call errcon ('Too many AKA aliases; rest skipped')
          call jvalut (' Max allowed    :',1,mxrtal)
          call jvalut (' Line number    :',1,nline)
          call textut (' Offending line :',line)
        else
          do i=1,nopt
            nalias = nalias + 1
            lalias (nalias) = optpar(i)(1:3)
          end do
        end if
c
      else if (line(1:3) .eq. 'MCH' .or.
     +         line(1:3) .eq. 'SCH') then
c
        lmore = .false.
        lmch  = (line(1:3) .eq. 'MCH')
        line  = line (4:)
   40   continue
        j = length (line)
        if (line(j:j) .eq. '-') then
          lmore = .true.
          line(j:) = ' '
        end if
   50   continue
        j = index (line,quote)
        if (j .le. 0) goto 60
        k = j + index (line(j+1:),quote)
        if (k .le. 0) k = j + 5
        if (nmatyp .ge. mxrtat) then
          call errcon ('Too many atom types')
          call jvalut (' Max allowed    :',1,mxrtat)
          call jvalut (' Line number    :',1,nline)
          call textut (' Offending line :',line)
          ierr = -1
          return
        end if
        nmatyp = nmatyp + 1
        lratom (nmatyp) = line(j+1:j+4)
        ismain (nmatyp) = lmch
ccc          call textut (' >',line)
ccc          print *,j,k,' |',lratom (nmatyp),'| '
        line = line(k+1:)
        goto 50
c
   60   continue
c
        if (lmore) then
          lmore = .false.
          read (iunit,6000,err=8000,end=9000) line
          nline = nline + 1
          goto 40
        end if
c
      else if (line(1:3) .eq. 'END') then
c
        if (nmatyp .eq. jj) then
          call errcon ('Residue type has no atoms')
          call textut (' Type :',lrname(ii))
          call textut (' Name :',lrdesc(ii))
          nmrtyp = nmrtyp - 1
          goto 10
        end if
c
        nmrptr (1,ii) = jj + 1
        nmrptr (2,ii) = nmatyp
c
        if (nalias .ne. kk) then
          nmaptr (1,ii) = kk + 1
          nmaptr (2,ii) = nalias
        end if
c
        goto 10
c
      else
c
        call errcon ('Invalid line in residue definition')
        call jvalut (' Line number    :',1,nline)
        call textut (' Offending line :',line)
        ierr = -1
        return
c
      end if
c
      goto 30
c
c ... RING DEFINITION
c
  120 continue
      if (nlring .ge. mxring) then
        call errcon ('Too many ring definitions')
        call jvalut (' Max allowed    :',1,mxring)
        call jvalut (' Line number    :',1,nline)
        call textut (' Offending line :',line)
        ierr = -1
        return
      end if
      nlring = nlring + 1
      rngres (nlring) = line(5:7)
      read (line(9:10),'(i2)',err=8000,end=8000) rngnat (nlring)
      if (rngnat(nlring) .gt. 10) then
        call errcon ('Max 10 ring atoms; rest ignored !')
        rngnat(nlring) = 10
      end if
      read (line(12:),'(10(a4,1x))')
     +  (rngatm(j,nlring),j=1,rngnat(nlring))
c
      goto 10
c
 8000 continue
      call errcon ('While reading file')
      ierr = -1
      return
c
 9000 continue
      write (*,*)
      call jvalut (' Lines read       :',1,nline)
      call jvalut (' Residue types    :',1,nmrtyp)
      call jvalut (' Atom types       :',1,nmatyp)
      call jvalut (' Aliases          :',1,nalias)
      call jvalut (' Ring definitions :',1,nlring)
c
      if (nmrtyp .lt. 10) then
        call errcon ('Fewer than 10 defined residue types')
        ierr = -1
        return
      end if
c
      call prompt ('0First and last residue types:')
      call tellib (lrname(1),i,.true.)
      call tellib (lrname(nmrtyp),i,.true.)
c
c ... check integrity
c
      call prompt ('0Check integrity:')
      lmore = .false.
      do i=1,nmrtyp
        do j=1,nmrtyp
          if (i.ne.j .and. lrname(i) .eq. lrname(j)) then
            lmore = .true.
            write (*,6200) i,lrname(i),j,lrname(j)
          end if
          if (nmaptr(1,j) .gt. 0 .and. i.ne.j) then
            do k=nmaptr(1,j),nmaptr(2,j)
              if (lalias(k) .eq. lrname(i)) then
                lmore = .true.
                write (*,6200) i,lrname(i),j,lalias(k)
              end if
            end do
          end if
        end do
      end do
      if (lmore) then
        call errcon ('Non-unique residue names/aliases')
      else
        call prompt ('There are no name or alias conflicts')
      end if
c
 6200 format (' WARNING - name or alias conflict: ',i5,' = ',
     +  a3,' and ',i5,' = ',a3)
c
      call prompt ('0Count types:')
      do i=1,nrlcat+1
        icnt(i) = 0
      end do
      do i=1,nmrtyp
        icnt(lrtype(i)) = icnt(lrtype(i)) + 1
      end do
      call jvalut (' Nr of amino acid residue types :',
     +  1,icnt(iprot))
      call jvalut (' Nr of nucleic acid types       :',
     +  1,icnt(inucl))
      call jvalut (' Nr of water types              :',
     +  1,icnt(iwate))
      call jvalut (' Nr of metal types              :',
     +  1,icnt(imeta))
      call jvalut (' Nr of inorganic types          :',
     +  1,icnt(iinor))
      call jvalut (' Nr of carbohydrate types       :',
     +  1,icnt(icarb))
      call jvalut (' Nr of organic compound types   :',
     +  1,icnt(iorga))
      call jvalut (' Nr of other compound types     :',
     +  1,icnt(ihete))
c
      if (icnt(iprot) .lt. 20) call errstp (
     +  'READDB - Not enough amino acid types')
      if (icnt(iwate) .lt. 1) call errstp (
     +  'READDB - No water type defined')
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine tellib (myname,ires,lprint)
c
      include 'moleman2.incl'
c
      integer i,j,k,ires,leng1
c
      logical lprint
c
      character myname*(*),name*3
c
code ...
c
      name = myname
      call upcase (name)
c
c ... first: check normal names
c
      do i=1,nmrtyp
        if (name .eq. lrname(i)) then
          ires = i
          if (lprint) then
            write (*,6100) i,lrname(i),libolc(i),
     +       libcat(lrtype(i)),
     +       lrdesc(i)(1:leng1(lrdesc(i)))
            if (nmaptr(1,i) .gt. 0) then
              write (*,6110) (lalias(j),j=nmaptr(1,i),nmaptr(2,i))
            end if
            write (*,6120)
     +        (lratom(j),ismain(j),j=nmrptr(1,i),nmrptr(2,i))
          end if
          return
        end if
      end do
c
c ... second: check aliases
c
      do i=1,nmrtyp
        if (nmaptr(1,i) .gt. 0) then
          do k=nmaptr(1,i),nmaptr(2,i)
            if (name .eq. lalias(k)) then
              ires = i
              if (lprint) then
                write (*,6100) i,lrname(i),libolc(i),
     +           libcat(lrtype(i)),
     +           lrdesc(i)(1:leng1(lrdesc(i)))
                if (nmaptr(1,i) .gt. 0) then
                  write (*,6110) (lalias(j),j=nmaptr(1,i),nmaptr(2,i))
                end if
                write (*,6120)
     +            (lratom(j),ismain(j),j=nmrptr(1,i),nmrptr(2,i))
              end if
              return
            end if
          end do
        end if
      end do
c
      ires = -1
      if (lprint) call textut (' Residue type not found :',name)
c
 6100 format (/' Residue # ',i5,' = ',a3,1x,a1,' (',a4,') = ',a)
 6110 format (10(' A.k.a. ',10(1x,a3),:,/))
 6120 format (10(' Atoms  ',6('  |',a4,'| (',l1,')',:),:,/))
c
      return
      end
c
c
c
      subroutine print_atom (kk)
c
      include 'moleman2.incl'
c
      integer kk,j,leng1
c
code ...
c
 6000 format (' ATOM ',i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a)
      write (*,6000) atomnr(kk),atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    (xyz(j,kk),j=1,3),qatom(kk),batom(kk),
     +    inote(kk)(1:leng1(inote(kk)))
c
      return
      end
c
c
c
      subroutine print_res (i,n)
c
c ... print N blank lines and then the name of the residue with pointer I
c
      include 'moleman2.incl'
c
      integer i,kk,n,j,length
c
      character mystr*80,ustr*(*)
c
code ...
c
      mystr = 'RESIDUE'
      goto 100
c
      entry string_res (i,n,ustr)
c
      mystr = ustr
c
  100 continue
c
      kk = atmptr(1,i)
      if (n .gt. 0) then
        do j=1,min(100,n)
          write (*,*)
        end do
      end if
c
 6000 format (1x,a,1x,a1,a3,1x,a1,i4,a1,1x,a4)
      write (*,6000) mystr(1:length(mystr)),altloc(kk),resnam(kk),
     +  achain(kk),iresid(kk),insert(kk),inote(kk)(7:10)
c
      return
      end
c
c
c
      subroutine do_select (how,what,which)
c
      include 'moleman2.incl'
c
      integer i,j
c
      logical lhydro
c
      character*(*) how,what,which
      character line*132
c
code ...
c
ccc      call textut (' How   :',how)
ccc      call textut (' What  :',what)
ccc      call textut (' Which :',which)
c
      if (how(1:2) .eq. 'AL') then
        call prompt (' Select ALL atoms')
        do i=1,natoms
          select(i) = .true.
        end do
        nselect = natoms
        selstr = 'ALL |'
c
      else if (how(1:1) .eq. '?') then
        goto 1000
c
      else if (how(1:2) .eq. 'NO') then
        call prompt (' Select NO atoms')
        do i=1,natoms
          select(i) = .false.
        end do
        nselect = 0
        selstr = 'NONE |'
c
      else if (how(1:2) .eq. 'HY') then
        call prompt (' Select HYDROGEN atoms')
        nselect = 0
        do i=1,natoms
          select(i) = lhydro(atmnam(i))
          if (select(i)) nselect = nselect + 1
        end do
        selstr = 'HYDROGEN |'
c
      else if (how(1:2) .eq. 'EX') then
        call prompt (' Select NON-HYDROGEN atoms')
        nselect = 0
        do i=1,natoms
          select(i) = (.not. lhydro(atmnam(i)))
          if (select(i)) nselect = nselect + 1
        end do
        selstr = 'NON-HYDROGEN |'
c
      else if (how(1:2) .eq. 'NE') then
        call prompt (' NEGATE atom selection')
        nselect = 0
        do i=1,natoms
          select(i) = (.not. select(i))
          if (select(i)) nselect = nselect + 1
        end do
        call appstr (selstr,' | NEGATE')
c
      else if (how(1:2) .eq. 'AN' .or.
     +         how(1:2) .eq. 'OR' .or.
     +         how(1:2) .eq. 'BU') then
c
        if (how(1:2) .eq. 'AN') then
          call prompt (' AND atom selection')
        else if (how(1:2) .eq. 'BU') then
          call prompt (' BUtnot atom selection')
        else if (how(1:2) .eq. 'OR') then
          call prompt (' OR atom selection')
        end if
        call textut (' With atoms for which :',what)
        call textut (' Equals :',which)
c
        if (what(1:2) .eq. 'CH') then
          what = 'CHain'
          do i=1,natoms
            lbuf (i) = (achain(i) .eq. which(1:1))
          end do
        else if (what(1:2) .eq. 'SE') then
          what = 'SEgid'
          do i=1,natoms
            lbuf (i) = (inote(i)(7:10) .eq. which(1:4))
          end do
        else if (what(1:2) .eq. 'RE') then
          what = 'REsidue'
          do i=1,natoms
            lbuf (i) = (resnam(i) .eq. which(1:3))
          end do
        else if (what(1:2) .eq. 'AT') then
          what = 'ATom'
          do i=1,natoms
            lbuf (i) = (atmnam(i) .eq. which(1:4))
          end do
        else if (what(1:2) .eq. 'AL') then
          what = 'ALtloc'
          do i=1,natoms
            lbuf (i) = (altloc(i) .eq. which(1:1))
          end do
        else if (what(1:2) .eq. 'AN') then
          what = 'ANisou'
          if (which(1:1) .eq. 'F') then
            do i=1,natoms
              lbuf (i) = (.not. laniso(i))
            end do
          else if (which(1:1) .eq. 'T') then
            do i=1,natoms
              lbuf (i) = (laniso(i))
            end do
          else
            call errcon ('ANiso must be True or False')
          end if
        else if (what(1:2) .eq. 'TY') then
          what = 'TYpe'
          j = -1
          do i=1,nrlcat+1
            if (libcat(i) .eq. which(1:4)) j=i
          end do
          if (j .le. 0) then
            call errcon ('Invalid type')
            call asciut (' Must be one of :',nrlcat+1,libcat)
            goto 1000
          end if
          which = libcat (j)
          do i=1,natoms
            lbuf (i) = (restyp(resptr(i)) .eq. j)
          end do
        else if (what(1:2) .eq. 'ST') then
          what = 'STruc_sec'
          j = -1
          if (which(1:2) .eq. 'LO') j=0
          if (which(1:2) .eq. 'TU') j=0
          if (which(1:2) .eq. 'AL') j=1
          if (which(1:2) .eq. 'BE') j=2
          if (which(1:2) .eq. 'LE') j=3
          if (which(1:2) .eq. 'NO') j=-1
          if (j .le. 0) then
            call errcon ('Invalid secondary structure type')
            call prompt (' Must be one of : LOop, TUrn, ALpha,')
            call prompt ('     BEta, LEft-handed, NOn-protein')
            goto 1000
          end if
          which = ssenam(j)
          do i=1,natoms
            lbuf (i) = (sstype(resptr(i)) .eq. j)
          end do
        else if (what(1:2) .eq. 'CL') then
          what = 'CLass'
          if (which(1:1) .eq. 'M') then
            which = 'Main'
            do i=1,natoms
              lbuf (i) = lmain (i)
            end do
          else if (which(1:1) .eq. 'S') then
            which = 'Side'
            do i=1,natoms
              lbuf (i) = (.not. lmain (i))
            end do
          else
            call errcon ('Invalid class')
            call prompt (' Must be Main or Side chain')
            goto 1000
          end if
        else
          call errcon ('Invalid AND/OR/BUTNOT option')
          goto 1000
        end if
c
        if (how(1:2) .eq. 'AN') then
          nselect = 0
          do i=1,natoms
            select(i) = (select(i) .and. lbuf(i))
            if (select(i)) nselect = nselect + 1
          end do
          line = ' AND '//what(1:6)//' = '//which(1:6)//' | '
          call pretty (line(2:))
          call appstr (selstr,line)
        else if (how(1:2) .eq. 'OR') then
          nselect = 0
          do i=1,natoms
            select(i) = (select(i) .or. lbuf(i))
            if (select(i)) nselect = nselect + 1
          end do
          line = ' OR '//what(1:6)//' = '//which(1:6)//' | '
          call pretty (line(2:))
          call appstr (selstr,line)
        else if (how(1:2) .eq. 'BU') then
          nselect = 0
          do i=1,natoms
            select(i) = (select(i) .and. (.not.lbuf(i)))
            if (select(i)) nselect = nselect + 1
          end do
          line = ' BUTNOT '//what(1:6)//' = '//which(1:6)//' | '
          call pretty (line(2:))
          call appstr (selstr,line)
        end if
c
      else
        call errcon ('Invalid SElection command')
      end if
c
 1000 continue
ccc      print *,'okay'
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine do_numeric_select (andor,what,lo,hi)
c
      include 'moleman2.incl'
c
      real rlo,rhi
c
      integer i,ilo,ihi,ierr,length,leng1,j
c
      character*(*) andor,what,lo,hi
      character line*132
c
code ...
c
      call remspa (andor)
      call remspa (what)
      call remspa (lo)
      call remspa (hi)
c
      if (andor(1:1) .eq. 'A') then
        line = ' AND '
      else if (andor(1:1) .eq. 'B') then
        line = ' BUTNOT '
      else
        line = ' OR '
      end if
c
      do i=1,natoms
        lbuf (i) = .false.
      end do
c
      if (what(1:1) .eq. 'R') then
        call str2i (lo,ilo,ierr)
        if (ierr .ne. 0) return
        call str2i (hi,ihi,ierr)
        if (ierr .ne. 0) return
        call ilohi (ilo,ihi)
        call appstr (line,' Residue_nr ')
        write (line(length(line)+1:),*) ilo,ihi
        do i=1,natoms
          lbuf (i) = (iresid(i) .ge. ilo .and.
     +                iresid(i) .le. ihi)
        end do
c
      else if (what(1:1) .eq. 'S') then
        call str2i (lo,ilo,ierr)
        if (ierr .ne. 0) return
        call str2i (hi,ihi,ierr)
        if (ierr .ne. 0) return
        call ilohi (ilo,ihi)
        call appstr (line,' Sec_struc ')
        write (line(length(line)+1:),*) ilo,ihi
        do i=1,natoms
          j = resptr (i)
          lbuf (i) = (sstype(j) .ge. ilo .and.
     +                sstype(j) .le. ihi)
        end do
c
      else if (what(1:1) .eq. 'E') then
        call str2i (lo,ilo,ierr)
        if (ierr .ne. 0) return
        call str2i (hi,ihi,ierr)
        if (ierr .ne. 0) return
        call ilohi (ilo,ihi)
        call appstr (line,' Element ')
        write (line(length(line)+1:),*) ilo,ihi
        do i=1,natoms
          lbuf (i) = (element(i) .ge. ilo .and.
     +                element(i) .le. ihi)
        end do
c
      else if (what(1:1) .eq. 'B') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' B-factor ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (batom(i) .ge. rlo .and.
     +                batom(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'O') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Occupancy ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (qatom(i) .ge. rlo .and.
     +                qatom(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'M') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Mass ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (atmass(i) .ge. rlo .and.
     +                atmass(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'C') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Cov_radius ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (cvbrad(i) .ge. rlo .and.
     +                cvbrad(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'X') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' X-coordinate ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (xyz(1,i) .ge. rlo .and.
     +                xyz(1,i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'Y') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Y-coordinate ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (xyz(2,i) .ge. rlo .and.
     +                xyz(2,i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'Z') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Z-coordinate ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (xyz(3,i) .ge. rlo .and.
     +                xyz(3,i) .le. rhi)
        end do
c
      else
        call errcon ('Invalid numeric property selected')
        return
      end if
c
      call pretty (line(2:))
      call textut (' Select Numeric :',line)
      line = line(1:leng1(line))//' |'
      call appstr (selstr,line)
c
      if (andor(1:1) .eq. 'A') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      else if (andor(1:1) .eq. 'B') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. (.not. lbuf(i)) )
          if (select(i)) nselect = nselect + 1
        end do
      else
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .or. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      end if
c
 1000 continue
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine do_point_select (andor,mx,my,mz,lo,hi)
c
      include 'moleman2.incl'
c
      real rlo,rhi,xxx,rxyz(3)
      real distce
c
      integer i,ierr,leng1
c
      character*(*) andor,mx,my,mz,lo,hi
      character line*132
c
code ...
c
      call remspa (andor)
      call remspa (mx)
      call remspa (my)
      call remspa (mz)
      call remspa (lo)
      call remspa (hi)
c
      if (andor(1:1) .eq. 'A') then
        line = ' AND '
      else if (andor(1:1) .eq. 'B') then
        line = ' BUTNOT '
      else
        line = ' OR '
      end if
c
      call appstr (line,' POint_distance ')
      call pretty (line(2:))
      line = line(1:leng1(line))//' |'
      call appstr (selstr,line)
c
      do i=1,natoms
        lbuf (i) = .false.
      end do
c
      call str2r (mx,rxyz(1),ierr)
      if (ierr .ne. 0) return
      call str2r (my,rxyz(2),ierr)
      if (ierr .ne. 0) return
      call str2r (mz,rxyz(3),ierr)
      if (ierr .ne. 0) return
      call str2r (lo,rlo,ierr)
      if (ierr .ne. 0) return
      call str2r (hi,rhi,ierr)
      if (ierr .ne. 0) return
      call rlohi (rlo,rhi)
c
      call fvalut (' Point (A) :',3,rxyz)
      call fvalut (' Minimum distance (A) :',1,rlo)
      call fvalut (' Maximum distance (A) :',1,rhi)
c
      do i=1,natoms
        xxx = distce (xyz(1,i),rxyz)
        lbuf (i) = (xxx .ge. rlo .and. xxx .le. rhi)
      end do
c
      if (andor(1:1) .eq. 'A') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      else if (andor(1:1) .eq. 'B') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. (.not. lbuf(i)) )
          if (select(i)) nselect = nselect + 1
        end do
      else
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .or. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      end if
c
 1000 continue
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine select_near (lo,hi)
c
      include 'moleman2.incl'
c
      real dist,qq,rlo,rhi
c
      integer i,j,ierr,leng1
c
      character*(*) lo,hi
      character line*132
c
code ...
c
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      if (nselect .lt. 1) then
        call errcon ('No atoms selected')
        return
      end if
c
      call remspa (lo)
      call remspa (hi)
c
      do i=1,natoms
        lbuf (i) = .false.
      end do
c
      call str2r (lo,rlo,ierr)
      if (ierr .ne. 0) return
      call str2r (hi,rhi,ierr)
      if (ierr .ne. 0) return
      call rlohi (rlo,rhi)
c
      line = ' DIstance '
      write (line(leng1(line)+1:),'(2f6.2)') rlo,rhi
      line = line(1:leng1(line))//' |'
      call appstr (selstr,line)
c
      do i=1,natoms
        if (select(i)) then
          do j=1,natoms
            if (i .eq. j) goto 1000
            if (select(j)) goto 1000
            if (abs(xyz(1,i)-xyz(1,j)) .gt. rhi) goto 1000
            if (abs(xyz(2,i)-xyz(2,j)) .gt. rhi) goto 1000
            if (abs(xyz(3,i)-xyz(3,j)) .gt. rhi) goto 1000
            qq = dist (i,j,xyz)
            if (qq .ge. rlo .and. qq .le. rhi) then
              lbuf (j) = .true.
            end if
 1000       continue
          end do
        end if
      end do
c
      nselect = 0
      do i=1,natoms
        select(i) = (select(i) .or. lbuf(i))
        if (select(i)) nselect = nselect + 1
      end do
c
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine select_by_residue ()
c
      include 'moleman2.incl'
c
      integer i,j,k,nok
c
code ...
c
      nok = 0
      call prompt (' Select by residue')
      do i=1,nres
        do j=atmptr(1,i),atmptr(2,i)
          if (select(j)) then
            do k=atmptr(1,i),atmptr(2,i)
              select (k) = .true.
              nok = nok + 1
            end do
            goto 10
          end if
        end do
   10   continue
      end do
c
      call appstr (selstr,' BY_residue | ')
c
      nselect = nok
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine list_select (which)
c
      include 'moleman2.incl'
c
      integer i,j,nok
c
      character which*(*)
c
code ...
c
      call upcase (which)
      call remspa (which)
c
      if (which(1:1) .ne. 'A') then
        call prompt (' List first selected atom of every residue')
        nok = 0
        do i=1,nres
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              call print_atom (j)
              nok = nok + 1
              goto 10
            end if
          end do
   10     continue
        end do
        call jvalut (' Nr of residues listed :',1,nok)
c
      else
        call prompt (' List all selected atoms')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            call print_atom (i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms listed :',1,nok)
      end if
c
      return
      end
c
c
c
      subroutine pdb_chem_charge (lprint,lcol)
c
      include 'moleman2.incl'
c
      real mass,radius
c
      integer i,j,ne,nz,np,nn
c
      logical lhydro,lprint,lcol
c
      character try*2,chem*2,charge*2,line*80
c
code ...
c
      if (lprint) call prompt (
     +  ' Deriving chemical name and charge ...')
c
      nz = 0
      nn = 0
      np = 0
      ne = 0
c
      do i=1,natoms
c
c ... is it hydrogen ?
c
        if (lhydro(atmnam(i))) then
          chem = ' H'
          call elinfo (chem,line,j,mass,radius,.false.)
          goto 50
        end if
c
c ... try first two characters of atom name
c
        try = atmnam(i)(1:2)
        call elinfo (try,line,j,mass,radius,.false.)
        if (j .ge. 1) then
          chem = try
          goto 50
        end if
c
c ... try space + second char of atom name
c
        try = ' '//atmnam(i)(2:2)
        call elinfo (try,line,j,mass,radius,.false.)
        if (j .ge. 1) then
          chem = try
          goto 50
        end if
c
c ... try space + first char of atom name
c
        try = ' '//atmnam(i)(1:1)
        call elinfo (try,line,j,mass,radius,.false.)
        if (j .ge. 1) then
          chem = try
          goto 50
        end if
c
c ... unknown element
c
        ne = ne + 1
        call errcon ('Unknown chemical element; assume Carbon')
        call print_atom (i)
        chem = ' C'
        call elinfo (chem,line,j,mass,radius,.false.)
c
c ... try to read charge from atom name
c
   50   continue
c
        element (i) = j
        atmass  (i) = mass
        cvbrad  (i) = 1.0
        if (radius .gt. 0.0) cvbrad (i) = radius
c
        if (lcol) then
          if (element(i) .eq. 1) then
            call xvrml_encode_rgb (1.0,0.0,1.0,colour(i))
          else if (element(i) .eq. 6) then
            call xvrml_encode_rgb (1.0,1.0,0.0,colour(i))
          else if (element(i) .eq. 7) then
            call xvrml_encode_rgb (0.0,0.0,1.0,colour(i))
          else if (element(i) .eq. 8) then
            call xvrml_encode_rgb (1.0,0.0,0.0,colour(i))
          else if (element(i) .eq. 16) then
            call xvrml_encode_rgb (0.0,1.0,0.0,colour(i))
          else if (element(i) .eq. 15) then
            call xvrml_encode_rgb (1.0,0.5,0.0,colour(i))
          else
            call xvrml_col_index (20+3*element(i),colour(i))
          end if
        end if
c
        if (atmnam(i)(3:4) .eq. '++') then
          charge = '+2'
        else if (atmnam(i)(3:4) .eq. '--') then
          charge = '-2'
        else if (atmnam(i)(3:3) .eq. '+') then
          charge = atmnam(i)(3:4)
        else if (atmnam(i)(3:3) .eq. '-') then
          charge = atmnam(i)(3:4)
        else if (atmnam(i)(4:4) .eq. '+') then
          charge = '+'//atmnam(i)(3:3)
        else if (atmnam(i)(4:4) .eq. '-') then
          charge = '-'//atmnam(i)(3:3)
        else
          charge = ' 0'
        end if
c
        if (charge(1:1) .eq. '+' .or. charge(1:1) .eq. '-') then
          if (charge(2:2) .eq. ' ') then
            charge (2:2) = '1'
          else if (charge(2:2) .lt. '1' .or.
     +             charge(2:2) .gt. '9') then
            charge = ' 0'
          end if
        end if
c
        if (lprint) then
          inote (i)(11:12) = chem
          inote (i)(13:14) = charge
        end if
c
        if (charge .eq. ' 0') then
          nz = nz + 1
        else if (charge(1:1) .eq. '+') then
          np = np + 1
        else if (charge(1:1) .eq. '-') then
          nn = nn + 1
        end if
c
      end do
c
      if (lprint) then
        call jvalut (' Nr of atoms processed    :',1,natoms)
        call jvalut (' Unknown chemical element :',1,ne)
        call jvalut (' Nr of positive atoms     :',1,np)
        call jvalut (' Nr of negative atoms     :',1,nn)
        call jvalut (' Nr uncharged or unknown  :',1,nz)
      end if
c
      return
      end
c
c
c
      subroutine pdb_number (first)
c
      include 'moleman2.incl'
c
      integer i,j,nok,ierr,i1,inow
c
      logical select_res
c
      character*(*) first
c
code ...
c
      call str2i (first,i1,ierr)
      if (ierr .ne. 0) return
c
      call jvalut (' Renumber selected residues starting at :',1,i1)
      nok = 0
      inow = i1
      do i=1,nres
        if (select_res(i)) then
          do j=atmptr(1,i),atmptr(2,i)
            iresid(j) = inow
          end do
          nok = nok + 1
          inow = inow + 1
        end if
      end do
      call jvalut (' Nr of last changed residue :',1,(inow-1))
      call jvalut (' Nr of residues changed :',1,nok)
c
      return
      end
c
c
c
      subroutine pdb_name (which,old,new)
c
      include 'moleman2.incl'
c
      integer i,j,nok,length
c
      logical select_res
c
      character*(*) which,old,new
c
code ...
c
      call upcase (which)
      call upcase (old)
      call upcase (new)
      call remspa (which)
c
      if (length(which) .lt. 0) then
        call errcon ('Invalid selection')
        return
      else if (length(old) .lt. 0) then
        call errcon ('Invalid value for old name')
        return
      else if (length(new) .lt. 0) then
        call errcon ('Invalid value for new name')
        return
      end if
c
      if (which(1:1) .ne. 'R') which(1:1) = 'A'
c
 6000 format (' Replace atom name |',a4,'| by |',a4,'|')
 6010 format (' Replace residue name |',a3,'| by |',a3,'|')
c
      if (which(1:1) .eq. 'A') then
        write (*,6000) old,new
        nok = 0
        do i=1,natoms
          if (select(i)) then
            if (atmnam(i) .eq. old(1:4)) then
              nok = nok + 1
              atmnam (i) = new(1:4)
            end if
          end if
        end do
        call jvalut (' Nr of atom names changed :',1,nok)
      else
        write (*,6010) old,new
        nok = 0
        do i=1,nres
          if (resnam(atmptr(1,i)) .eq. old(1:3)) then
            if (select_res(i)) then
              do j=atmptr(1,i),atmptr(2,i)
                resnam(j) = new(1:3)
              end do
              nok = nok + 1
            end if
          end if
        end do
        call jvalut (' Nr of residues changed :',1,nok)
      end if
c
      return
      end
c
c
c
      subroutine pdb_remark (what,which)
c
      include 'moleman2.incl'
c
      integer i,ii,ierr,leng1
c
      character what*(*),which*(*)
c
code ...
c
      call upcase (what)
      call remspa (what)
      if (what(1:1) .ne. 'D' .and. what(1:1) .ne. 'A') what = 'L'
c
 6000 format (1x,i5,': ',a)
c
      if (what(1:1) .eq. 'L') then
        call prompt (' List REMARK records')
        call jvalut (' Nr of REMARK records :',1,nrem)
        if (nrem .gt. 0) then
          do i=1,nrem
            write (*,6000) i,remark(i)(1:leng1(remark(i)))
          end do
        end if
c
      else if (what(1:1) .eq. 'D') then
        if (which(1:1) .eq. '*') then
          call prompt (' Delete all REMARK records')
          nrem = 0
        else
          call str2i (which,ii,ierr)
          if (ierr .ne. 0) return
          if (ii .le. 0 .or. ii .gt. nrem) then
            call errcon ('Invalid REMARK number')
            return
          end if
          call jvalut (' Delete REMARK record nr:',1,ii)
          if (ii .eq. nrem) then
            nrem = nrem - 1
          else
            do i=ii,nrem-1
              remark(i) = remark(i+1)
            end do
            nrem = nrem - 1
          end if
        end if
c
      else if (what(1:1) .eq. 'A') then
        call textut (' Add REMARK record :',which)
        if (nrem .lt. maxcom) then
          nrem = nrem + 1
          remark (nrem) = 'REMARK '//which(1:leng1(which))
          write (*,6000) nrem,remark(nrem)(1:leng1(remark(nrem)))
        else
          call errcon ('Too many REMARK records')
          call jvalut (' Maximum :',1,maxcom)
        end if
      end if
c
      return
      end
c
c
c
      subroutine pdb_farout ()
c
      include 'moleman2.incl'
c
      integer i,i1,i2,nah,nbs,ityp,ifirst,nn
c
      logical busy
c
      character line*128
c
code ...
c
      call prompt (
     +  ' Generating quick-n-dirty HELIX and SHEET records ...')
c
      call pdb_remark ('A',
     +  'YASSPA quick-n-dirty HELIX and SHEET records')
c
      nah = 0
      nbs = 0
c
      busy = .false.
      ityp = 0
c
      do i=1,nres
c
c ... start of a new SSE ?
c
        if (.not. busy .and. sstype(i) .gt. 0) then
          ityp = sstype (i)
          busy = .true.
          ifirst = i
          goto 50
        end if
c
c ... end of an SSE ?
c
        if (busy .and.
     +      (sstype(i+1) .ne. ityp .or. i .eq. nres)) then
          nn = i - ifirst + 1
          i1 = atmptr (1,ifirst)
          i2 = atmptr (1,i)
          if (nn .ge. 3) then
            if (ityp .eq. 1) then
              nah = nah + 1
              write (line,6000) 'HELIX ',nah,nah,
     +          resnam(i1),achain(i1),iresid(i1),insert(i1),
     +          resnam(i2),achain(i2),iresid(i2),insert(i2),
     +          1
            else
              nbs = nbs + 1
              write (line,6100) 'SHEET ',nbs,nbs,1,
     +          resnam(i1),achain(i1),iresid(i1),insert(i1),
     +          resnam(i2),achain(i2),iresid(i2),insert(i2),
     +          0
            end if
c
            call textut (' Record :',line)
            if (nother .lt. maxcom) then
              nother = nother + 1
              other (nother) = line
            else
              call errcon ('Too many records')
              call jvalut (' Maximum :',1,maxcom)
            end if
c
          end if
          busy = .false.
          ityp = 0
          goto 50
        end if
c
   50   continue
      end do
c
 6000 format (a6,1x,i3,1x,i3,2(1x,a3,1x,a1,1x,i4,a1),i2)
 6100 format (a6,1x,i3,1x,i3,i2,2(1x,a3,1x,a1,i4,a1),i2)
c
cATOM    230  CA  ARG    29      23.498  24.524  11.798  1.00 18.43      1CBS 452
cATOM     67  CA  MET A   9      -2.009  19.060  -6.222  1.00 33.52      1CBR 370
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
cHELIX    1   1 PHE     15  LEU     22  1
cHELIX    1   1 GLY A   41  MET A   43  5                                1FSS 170
cSHEET    6   6 1 GLU    70  GLU    72  0) 
cSHEET    1   A 3 LEU A   7  THR A  10  0                                1FSS 194
cSHEET    2   A 3 GLY A  13  MET A  16 -1  N  VAL A  15   O  VAL A   8   1FSS 195
cSHEET    3   A 3 VAL A  57  ASN A  59  1  N  TRP A  58   O  LYS A  14   1FSS 196
c
      call jvalut (' Nr of HELIX records :',1,nah)
      call jvalut (' Nr of SHEET records :',1,nbs)
c
      return
      end
c
c
c
      subroutine pdb_seqres ()
c
      include 'moleman2.incl'
c
      integer maxchn
      parameter (maxchn = 100)
c
      integer i,i1,i2,i3,i4,j,k,ndone
c
      logical rsdone(maxres)
c
      character seqres(maxres)*3
      character line*128,chdone*(maxchn)
c
code ...
c
      call prompt (
     +  ' Generating quick-n-dirty SEQRES records ...')
c
      call pdb_remark ('A',
     +  'Quick-n-dirty SEQRES records')
c
      if (nres .lt. 3) then
        call errcon ('Fewer than 3 residues')
        return
      end if
c
      chdone = ' '
      ndone = 0
      do i=1,nres
        rsdone (i) = .true.
        if (restyp (i) .eq. iprot .or.
     +      restyp (i) .eq. inucl) then
          rsdone (i) = .false.
          i1 = atmptr (1,i)
          if (ndone .lt. 1) then
            ndone = 1
            chdone = achain (i1)
          else
            if (index(chdone,achain(i1)) .le. 0) then
              if (ndone .lt. maxchn) then
                ndone = ndone + 1
                chdone(ndone:ndone) = achain(i1)
              else
                call errcon ('Too many chain IDs - rest skipped')
              end if
            end if
          end if
        end if
      end do
c
      if (ndone .lt. 1) then
        call errcon ('No chain IDs found ???')
        return
      end if
c
      call textut (' Chain IDs found :',chdone)
c
      do j=1,ndone
c
        i2 = 0
        do i=1,nres
          if (.not. rsdone(i)) then
            i1 = atmptr (1,i)
            if (achain(i1) .eq. chdone(j:j)) then
              i2 = i2 + 1
              seqres (i2) = resres (i)
            end if
          end if
        end do
c
        write (*,*)
        call textut (' Chain ID :',chdone(j:j))
        call jvalut (' Nr of residues :',1,i2)
c
        if (i2 .gt. 0) then
c
          i4 = 0
          do i=1,i2,13
            i3 = min(i2,i+12)
            i4 = i4 + 1
            write (line,6000) i4,chdone(j:j),i2,(seqres(k),k=i,i3)
c
ccc            write (*,'(1x,a)') line(1:length(line))
c
            if (nother .lt. maxcom) then
              nother = nother + 1
              other (nother) = line
            else
              call errcon ('Too many records')
              call jvalut (' Maximum :',1,maxcom)
              return
            end if
          end do
          call jvalut (' Nr of SEQRES records :',1,i4)
c
        end if
c
      end do
c
 6000 format ('SEQRES',i4,1x,a1,i5,1x,13(1x,a3))
c
      return
      end
c
c
c
      subroutine pdb_ssbond (mywhat)
c
      include 'moleman2.incl'
c
      real dist,dd
c
      integer i,nok,j,ii,jj,ncys,leng1
c
      character mywhat*(*),what*10
c
code ...
c
      what = mywhat
      call upcase (what)
      call remspa (what)
      if (what(1:1) .ne. 'D' .and. what(1:1) .ne. 'G') what = 'L'
c
      if (what(1:1) .eq. 'L') then
        call prompt (' List SSBOND records')
        nok = 0
        if (nother .gt. 0) then
          do i=1,nother
            if (other(i)(1:6) .eq. 'SSBOND') then
              write (*,'(1x,a)') other(i)(1:leng1(other(i)))
              nok = nok + 1
            end if
          end do
        end if
        call jvalut (' Nr of SSBOND records listed :',1,nok)
c
      else if (what(1:1) .eq. 'D') then
        call prompt (' Delete SSBOND records')
        nok = 0
        if (nother .gt. 0) then
          j = 0
          do i=1,nother
            if (other(i)(1:6) .eq. 'SSBOND') then
              nok = nok + 1
            else
              j = j + 1
              if (i .ne. j) other(j) = other(i)
            end if
          end do
          nother = j
        end if
        call jvalut (' Nr of SSBOND records deleted :',1,nok)
c
      else if (what(1:1) .eq. 'G') then
        call prompt (' Generate SSBOND records')
        nok = 0
        call find_type (natoms,resnam,atmnam,'CYS',' SG ',
     +    ncys,ibuf,.true.)
        if (ncys .lt. 2) then
          call prompt (' No disulfide links')
          return
        end if
c
        call fvalut (' Max SG-SG distance for link :',1,mxcyss)
        nok = 0
        do i=1,ncys-1
          ii = ibuf(i)
          do j=i+1,ncys
            jj = ibuf(j)
            dd = dist (ii,jj,xyz)
            if (dd .le. mxcyss) then
              if (nother .eq. maxcom) then
                call errcon ('No room for more records')
                call jvalut (' Maximum :',1,maxcom)
                return
              end if
              nok = nok + 1
              nother = nother + 1
              write (other(nother),
     + '(a6,1x,i3,1x,a3,a2,1x,i4,a1,3x,a3,a2,1x,i4,a1,4x,a,f6.2,a)')
     +            'SSBOND',nok,'CYS',achain(ii),iresid(ii),
     +            insert(ii),'CYS',achain(jj),iresid(jj),
     +            insert(jj),'S-S = ',dd,' A'
              write (*,'(1x,a)') other(nother)(1:leng1(other(nother)))
            end if
          end do
        end do
c
        call jvalut (' Nr of SSBOND records generated :',1,nok)
c
      end if
c
      return
      end
c
c
c
      subroutine pdb_hetero (what)
c
      include 'moleman2.incl'
c
      integer i,j,n1,n2
c
      character what*(*)
c
code ...
c
      call upcase (what)
      call remspa (what)
      if (what(1:1) .ne. 'A') what = 'D'
c
      if (what(1:1) .eq. 'A') then
        do i=1,natoms
          lhet (i) = .false.
        end do
        call prompt (' All atoms set to type ATOM')
      else if (what(1:1) .eq. 'D') then
        n1 = 0
        n2 = 0
        call prompt (' Deducing ATOM/HETATM types ...')
        do i=1,natoms
          j = restyp(resptr(i))
          if (j .eq. iprot .or. j .eq. inucl) then
            lhet (i) = .false.
            n1 = n1 + 1
          else
            lhet (i) = .true.
            n2 = n2 + 1
          end if
        end do
        call jvalut (' Nr set to ATOM   :',1,n1)
        call jvalut (' Nr set to HETATM :',1,n2)
      end if
c
      return
      end
c
c
c
      subroutine mc_check (iunit,file,what,which)
c
      include 'moleman2.incl'
c
      real dist,tangle,omega,plan,phi,psi,dx,dy,ptsize,chir
      real lefthh(4)
c
      integer coregn(37,37)
      integer i,nok,j,iphi,ipsi,nend,nerr,ngly,nout,nyes,ndaa
      integer ncis,nraar,nbent,iunit,mode,length,ierr,kk,leng1
      integer nrleft,i1left,inleft
c
      logical lstart,lend,lomega,lcis,lraar,lplan,lbent,lpsi,lphi
      logical lpos,lrama,lp,label,ldaa,lchir
c
      character*(*) file,what,which
      character line*256,mychn*1
c
      data lefthh /30.0,130.0,-50.0,100.0/
c
code ...
c
      call defcor (coregn)
c
      mode = 0
      if (length(file) .gt. 0) then
        call upcase (what)
        call remspa (what)
        if (what(1:1) .eq. 'T') then
          mode = 1
        else if (what(1:1) .eq. 'L') then
          mode = 2
          label = .true.
        else
          mode = 2
          label = .false.
        end if
      end if
c
      mychn = which(1:1)
      call upcase (mychn)
c
      if (mode .eq. 1) then
c
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
c
      else if (mode .eq. 2) then
c
        dx = 2.0
        dy = 2.0
        ptsize = 10.0
c
        call psrini (iunit,file,prognm)
c
      end if
c
      do i=1,5*nres
        ibuf (i) = -1
      end do
c
      do i=1,2*nres
        rbuf (i) = -999.99
      end do
c
      write (*,*)
      call textut (' Chain ID to check (* = all) :',mychn)
c
      nok = 0
      do i=1,nres
        lbuf (i) = (restyp(i).eq.iprot)
        if (mychn .ne. '*') then
          lbuf(i) = ( lbuf(i) .and.
     +                (achain(atmptr(1,i)) .eq. mychn) )
        end if
        if (lbuf(i)) then
          nok = nok + 1
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' N  ') then
              ibuf (i) = j
            else if (atmnam(j) .eq. ' CA ') then
              ibuf (nres+i) = j
            else if (atmnam(j) .eq. ' C  ') then
              ibuf (2*nres+i) = j
            else if (atmnam(j) .eq. ' O  ') then
              ibuf (3*nres+i) = j
            else if (atmnam(j) .eq. ' CB ') then
              ibuf (4*nres+i) = j
            end if
          end do
        end if
      end do
c
      call jvalut (' Nr of residues to check :',1,nok)
      if (nok .lt. 3) then
        call errcon ('I don''t call this a protein ...')
        call xps_delete ()
        call prompt (' PostScript file empty and deleted')
        return
      end if
c
      nend = 0
      nerr = 0
      ngly = 0
      ndaa = 0
      nout = 0
      nyes = 0
      ncis = 0
      nraar = 0
      nbent = 0
c
      nrleft = 0
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,4(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,4(1x,f10.2))
 6010 format (/
     +  ' Listing of Phi, Psi, Omega and Planarity (pseudo) torsion'/
     +  ' angles for PDB file ',a,' chain ',a1/
     + /' PHI = C(i-1) - N(i) - CA(i) - C(i) - usually negative'/
     +  ' PSI = N(i) - CA(i) - C(i) - N(i+1)'/
     +  ' OMEGA = CA(i) - C(i) - N(i+1) - CA (i+1) - cis=0, trans=180'/
     +  ' PLANARITY = C(i) - CA(i) - N(i+1) - O(i) - planar < 5'/
     + /' An entry of "-999.9" means that the torsion angle could'/
     +  ' not be calculated (for terminal residues and residues with'/
     +  ' missing atoms).'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil)),mychn
        write (iunit,6100,err=9000)
     +    'Residue','Phi','Psi','Omega','Planarity'
        write (iunit,6100,err=9000)
     +    '-------','---','---','-----','---------'
      end if
c
      do i=1,nres
c
        if (.not. lbuf(i)) then
          nrleft = 0
          goto 1000
        end if
c
        if (captr(i) .le. 0) then
          call print_res (i,1)
          call prompt (' ERROR - missing CA atom')
          nerr = nerr + 1
          nrleft = 0
          goto 1000
        end if
c
        lstart = .false.
        lend = .false.
        lomega = .false.
        lcis = .false.
        lraar = .false.
        lchir = .false.
        ldaa = .false.
        lbent = .false.
        lphi = .false.
        lpsi = .false.
        lpos = .false.
        lrama = .false.
c
        omega = -999.9
        plan = -999.9
        phi = -999.9
        psi = -999.9
c
        if (i .eq. 1) then
          lstart = .true.
        else if (captr(i-1) .gt. 0) then
          lstart =  (dist(captr(i-1),captr(i),xyz) .gt. mxcaca)
        else
          lstart = .true.
        end if
c
        if (i .eq. nres) then
          lend = .true.
        else  if (captr(i+1) .gt. 0) then
          lend =  (dist(captr(i+1),captr(i),xyz) .gt. mxcaca)
        else
          lend = .true.
        end if
c
        if (.not. lend) then
c
c ... OMEGA = CA(i) - C(i) - N(i+1) - CA (i+1)
c
          lomega = (ibuf(nres+i) .gt. 0 .and.
     +              ibuf(2*nres+i) .gt. 0 .and.
     +              ibuf(i+1) .gt. 0 .and.
     +              ibuf(nres+i+1) .gt. 0)
          if (lomega) then
            omega = tangle (ibuf(nres+i),ibuf(2*nres+i),
     +                ibuf(i+1),ibuf(nres+i+1),xyz)
            lcis = (abs(omega) .le. 30)
            lraar = (abs(omega) .gt. 30 .and.
     +               abs(omega) .lt. 150)
          end if
c
c ... PLANARITY = C(i) - CA(i) - N(i+1) - O(i)
c
          lplan = (ibuf(2*nres+i) .gt. 0 .and.
     +             ibuf(nres+i) .gt. 0 .and.
     +             ibuf(i+1) .gt. 0 .and.
     +             ibuf(3*nres+i) .gt. 0)
          if (lplan) then
            plan = tangle (ibuf(2*nres+i),ibuf(nres+i),
     +               ibuf(i+1),ibuf(3*nres+i),xyz)
            lbent = (abs(plan) .gt. 5.0)
          end if
c
c ... PSI = N(i) - CA(i) - C(i) - N(i+1)
c
          lpsi = (ibuf(i) .gt. 0 .and.
     +            ibuf(nres+i) .gt. 0 .and.
     +            ibuf(2*nres+i) .gt. 0 .and.
     +            ibuf(i+1) .gt. 0)
          if (lpsi) then
            psi = tangle (ibuf(i),ibuf(nres+i),
     +              ibuf(2*nres+i),ibuf(i+1),xyz)
          end if
        end if
c
c ... PHI = C(i-1) - N(i) - CA(i) - C(i)
c
        if (.not. lstart) then
          lphi = (ibuf(2*nres+i-1) .gt. 0 .and.
     +            ibuf(i) .gt. 0 .and.
     +            ibuf(nres+i) .gt. 0 .and.
     +            ibuf(2*nres+i) .gt. 0)
          if (lphi) then
            phi = tangle (ibuf(2*nres+i-1),ibuf(i),
     +              ibuf(nres+i),ibuf(2*nres+i),xyz)
            if (resnam(atmptr(1,i)) .ne. 'GLY') then
              lpos = (phi .ge. 0.0)
            end if
          end if
        end if
c
c ... CHIR = CA - N - C - CB
c
        lchir = (ibuf(nres+i) .gt. 0 .and.
     +           ibuf(i) .gt. 0 .and.
     +           ibuf(2*nres+i) .gt. 0 .and.
     +           ibuf(4*nres+i) .gt. 0)
        if (lchir) then
          chir = tangle (ibuf(nres+i),ibuf(i),
     +             ibuf(2*nres+i),ibuf(4*nres+i),xyz)
          ldaa = (chir .lt. 0.0)
        end if
c
        if (lphi .and. lpsi) then
          rbuf (i) = phi
          rbuf (nres+i) = psi
          if (resnam(atmptr(1,i)) .ne. 'GLY') then
            if (.not. ldaa) then
              iphi = int ( (180.0 + phi) / 10.0 ) + 1
              ipsi = int ( (180.0 + psi) / 10.0 ) + 1
              lrama = (coregn(iphi,ipsi) .eq. 1)
            else
              iphi = int ( (180.0 - phi) / 10.0 ) + 1
              ipsi = int ( (180.0 - psi) / 10.0 ) + 1
              lrama = (coregn(iphi,ipsi) .eq. 1)
              iphi = int ( (180.0 + phi) / 10.0 ) + 1
              ipsi = int ( (180.0 + psi) / 10.0 ) + 1
            end if
          end if
        end if
c
        lp = .false.
c
        if (lstart) then
          call print_res(i,2)
          lp = .true.
          call prompt (' <<<<< Start of new chain')
          if (mode .eq. 1) write (iunit,*,err=9000)
        end if
c
        if (lend) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' >>>>> End of chain')
        end if
c
        if (.not. lend) then
          if (.not. lomega) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate OMEGA')
          else
            if (lcis) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,f8.1,a,a3)')
     +          ' Cis-peptide; omega = ',omega,'; next residue is ',
     +          resnam(atmptr(1,i+1))
              call prompt (line)
            end if
            if (lraar) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - strange omega = ',omega
              call prompt (line)
            end if
          end if
          if (.not. lplan) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate planarity')
          else
            if (lbent) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - non-planar peptide; planarity = ',plan
              call prompt (line)
            end if
          end if
          if (.not. lpsi) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate PSI')
          end if
        end if
c
        if (.not. lstart) then
          if (.not. lphi) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate PHI')
          else
            if (resnam(atmptr(1,i)) .ne. 'GLY') then
              if (lpos) then
                if (.not. lp) call print_res(i,1)
                lp = .true.
                write (line,'(a,f8.1)')
     +            ' Warning - positive PHI = ',phi
                call prompt (line)
              end if
            end if
          end if
        end if
c
        if (lchir) then
          if (ldaa) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            write (line,'(a,f8.1)')
     +        ' Warning - D-amino acid; improper = ',chir
            call prompt (line)
          end if
        end if
c
        if (lphi .and. lpsi) then
          if (resnam(atmptr(1,i)) .ne. 'GLY') then
c
c ... plus for all non-glycines
c
            if (mode .eq. 2) then
              if (.not. ldaa) then
                call xps_colour (0)
                call xps_symbol (1,phi-dx,phi+dx,psi-dy,psi+dy)
              else
                call xps_colour (1)
                call xps_symbol (4,phi-dx,phi+dx,psi-dy,psi+dy)
              end if
            end if
            if (.not. lrama) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,2f8.1)')
     +          ' Warning - unusual PHI,PSI combination = ',phi,psi
              call prompt (line)
c
c ... add cross for outliers and optional label with residue name etc.
c
              if (mode .eq. 2) then
                call xps_symbol (2,phi-dx,phi+dx,psi-dy,psi+dy)
                if (label) then
                  kk = atmptr(1,i)
                  write (line,6000) resnam(kk),achain(kk),iresid(kk),
     +              insert(kk),inote(kk)(7:10)
                  call remspa(line)
                  call xps_colour (4)
                  call xps_text (phi+dx+0.5,psi-dy-0.5,ptsize,line)
                end if
              end if
            end if
          else
c
c ... square for glycines
c
            if (mode .eq. 2) then
              call xps_colour (0)
              call xps_symbol (0,phi-dx,phi+dx,psi-dy,psi+dy)
            end if
          end if
        end if
c
        if (lstart .or. lend) then
          nend = nend + 1
        else if (.not. (lphi .and. lpsi)) then
          nerr = nerr + 1
        else if (resnam(atmptr(1,i)) .eq. 'GLY') then
          ngly = ngly + 1
        else if (.not. lrama) then
          nout = nout + 1
        else
          nyes = nyes + 1
        end if
        if (lcis) ncis = ncis + 1
        if (lraar) nraar = nraar + 1
        if (lbent) nbent = nbent + 1
        if (ldaa) ndaa = ndaa + 1
c
        if (mode .eq. 1) then
          kk = atmptr(1,i)
          write (iunit,6200,err=9000)
     +      resnam(kk),achain(kk),iresid(kk),
     +      insert(kk),inote(kk)(7:10),phi,psi,omega,plan
        end if
c
c ... check for left-handed helices
c
        if (lphi .and. lpsi) then
          if (phi .ge. lefthh(1) .and. phi .le. lefthh(2) .and.
     +        psi .ge. lefthh(3) .and. psi .le. lefthh(4)) then
            nrleft = nrleft + 1
            if (nrleft .eq. 1) i1left = i
            inleft = i
            goto 1000
          end if
        end if
c
        if (nrleft .ge. 4) then
          write (*,*)
          call jvalut (
     +      ' Left-handed helix - nr of residues :',1,nrleft)
          call string_res (i1left,0,'Left-handed helix first :')
          call string_res (inleft,0,'Left-handed helix last  :')
        end if
        nrleft = 0
c
 1000   continue
c
      end do
c
      write (*,*)
      call textut (' Chain ID checked (* = all)   :',mychn)
      call jvalut (' Total nr of residues checked :',1,nok)
      call jvalut (' D-amino acids                :',1,ndaa)
      call jvalut (' Cis-peptide bonds            :',1,ncis)
      call jvalut (' Unusual OMEGA values         :',1,nraar)
      call jvalut (' Non-planar peptide bonds     :',1,nbent)
      call jvalut (' Start/end residues           :',1,nend)
      call jvalut (' Problems (missing atoms ?)   :',1,nerr)
      call jvalut (' Glycine residues             :',1,ngly)
      call jvalut (' Remaining residues in Ramach :',1,(nout+nyes))
      call jvalut ('   In core regions            :',1,nyes)
      call jvalut ('   Outliers                   :',1,nout)
      if ( (nout+nyes) .gt. 0) then
        call fvalut ('   Outlier percentage         :',1,
     +    100.0*float(nout)/float(nout+nyes))
      end if
      call prompt (' An average <= 2.0 A model has ~0-5% outliers')
      call prompt (' See: Kleywegt, G.J. and Jones, T.A. (1996).')
      call prompt ('      Structure 4, 1395-1400.')
      if (ndaa .gt. 0) then
        call prompt (' For D-amino acids, -Phi and -Psi were used')
      end if
c
      if (mode .eq. 1) then
c
        write (iunit,*)
        close (iunit)
c
      else if (mode .eq. 2) then
c
        line = ' PDB file : '//pdbfil
        call xps_legend (line)
c
        write (line,*)
     +    ' Glycines (open squares):',ngly,
     +    ' ; Start/end residues :',nend
        call xps_legend (line)
c
        write (line,*)
     +    ' D-amino acids :',ndaa,
     +    ' ; Residues with missing atoms :',nerr
        call xps_legend (line)
c
        write (line,*) ' Residues in Ramachandran plot checked :',
     +    (nout+nyes),' out of ',nok
        call xps_legend (line)
c
        write (line,*) ' In core regions (plus signs):',nyes,
     +    ' ; Outliers (asterisks):',nout
        call xps_legend (line)
c
        if ( (nout+nyes) .gt. 0) then
          write (line,'(a,1x,f8.1)')
     +      ' Percentage outliers:',
     +      100.0*float(nout)/float(nout+nyes)
          call xps_legend (line)
        end if
c
        line = ' An average <= 2.0 A model has ~0-5% outliers'
        call xps_legend (line)
c
        line = ' See: Kleywegt, G.J. and Jones, T.A. (1996). '//
     +         'Structure 4, 1395-1400.'
        call xps_legend (line)
c
        if (ndaa .gt. 0) then
          line = ' For D-amino acids, -Phi and -Psi were used !'
          call xps_legend (line)
        end if
c
        if ((nout+nyes) .gt. 0) then
          call xps_close ()
          call prompt (' PostScript file created')
        else
          call xps_delete ()
          call prompt (' PostScript file empty and deleted')
        end if
c
      end if
c
      return
c
 9000 continue
      call errcon ('While writing text file')
c
      return
      end
c
c
c
      subroutine sc_check (iunit,file)
c
      include 'moleman2.incl'
c
      real chi(5)
      real tangle,chir,cideal,ctoler,xbad
c
      integer iptr(9)
      integer i,nok,j,iunit,mode,length,ierr,k,nchi,kk
      integer ndaa,nbad,nswap,leng1
c
      logical llaa,lbadc,lp,lswap
c
      character*(*) file
      character line*256
c
code ...
c
      cideal = 34.0
      ctoler = tortol
c
      mode = 0
      if (length(file) .gt. 0) then
        mode = 1
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
      end if
c
      call fvalut (' Ideal chirality improper :',1,cideal)
      call fvalut (' Ideal flatness improper  :',1,0.0)
      call fvalut (' Tolerance                :',1,ctoler)
      write (*,*)
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,6(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,6(1x,f10.2))
 6010 format (/
     +  ' Listing of Chi1-5 and Chirality (pseudo) torsion'/
     +  ' angles for PDB file ',a/
     + /' CHIRAL = CA    - N     - C     - CB (+34=L-aa; -34=D-aa)'/
     +  ' CHI-1  = N     - CA    - CB    - ?G(1)'/
     +  ' CHI-2  = CA    - CB    - ?G(1) - ?D(1)'/
     +  ' CHI-3  = CB    - ?G(1) - ?D(1) - ?E(1)'/
     +  ' CHI-4  = ?G(1) - ?D(1) - ?E(1) - ?Z(1)'/
     +  ' CHI-5  = ?D(1) - ?E(1) - ?Z(1) - ?H(1)'/
     + /' An entry of "-999.9" means that the torsion angle could'/
     +  ' not be calculated (e.g., for residues with missing atoms).'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil))
        write (iunit,6100,err=9000)
     +    'Residue','CA-chir','Chi-1','Chi-2','Chi-3','Chi-4','Chi-5'
        write (iunit,6100,err=9000)
     +    '-------','-------','-----','-----','-----','-----','-----'
      end if
c
      nok = 0
      ndaa = 0
      nbad = 0
      nswap = 0
c
      do i=1,nres
        if (restyp(i) .ne. iprot) goto 100
        if (resnam(atmptr(1,i)) .eq. 'GLY') goto 100
c
        do k=1,9
          iptr(k) = -1
        end do
c
        do k=1,5
          chi(k) = -999.9
        end do
c
        llaa = .true.
        lbadc = .false.
        lswap = .false.
c
        do j=atmptr(1,i),atmptr(2,i)
          if (atmnam(j) .eq. ' C  ') then
            iptr(1) = j
          else if (atmnam(j) .eq. ' N  ') then
            iptr(2) = j
          else if (atmnam(j) .eq. ' CA ') then
            iptr(3) = j
          else if (atmnam(j) .eq. ' CB ') then
            iptr(4) = j
          else if (atmnam(j)(3:4) .eq. 'G ' .or.
     +             atmnam(j)(3:4) .eq. 'G1') then
            iptr(5) = j
          else if (atmnam(j)(3:4) .eq. 'D ' .or.
     +             atmnam(j)(3:4) .eq. 'D1') then
            iptr(6) = j
          else if (atmnam(j)(3:4) .eq. 'E ' .or.
     +             atmnam(j)(3:4) .eq. 'E1') then
            iptr(7) = j
          else if (atmnam(j)(3:4) .eq. 'Z ' .or.
     +             atmnam(j)(3:4) .eq. 'Z1') then
            iptr(8) = j
          else if (atmnam(j)(3:4) .eq. 'H ' .or.
     +             atmnam(j)(3:4) .eq. 'H1') then
            iptr(9) = j
          end if
        end do
c
c ... CA chirality
c
        if (iptr(1).le.0 .or. iptr(2).le.0 .or.
     +      iptr(3).le.0 .or. iptr(4).le.0) then
          call print_res (i,1)
          call prompt (' ERROR - missing atom(s)')
          goto 100
        end if
c
        chir = tangle (iptr(3),iptr(2),iptr(1),iptr(4),xyz)
        if (chir .lt. 0) llaa = .false.
        if (llaa) then
          lbadc = ( abs(chir - cideal) .gt. ctoler)
        else
          lbadc = ( abs(-chir - cideal) .gt. ctoler)
        end if
c
        nchi = 0
        if (iptr(5) .le. 0) goto 90
        chi(1) = tangle (iptr(2),iptr(3),iptr(4),iptr(5),xyz)
        nchi = nchi + 1
c
        if (iptr(6) .le. 0) goto 90
        chi(2) = tangle (iptr(3),iptr(4),iptr(5),iptr(6),xyz)
        nchi = nchi + 1
c
        if (resnam(atmptr(1,i)) .eq. 'ASP' .or.
     +      resnam(atmptr(1,i)) .eq. 'PHE' .or.
     +      resnam(atmptr(1,i)) .eq. 'TYR') then
          lswap = (abs(chi(2)) .gt. 90.0)
          xbad = chi(2)
        end if
c
        if (iptr(7) .le. 0) goto 90
        chi(3) = tangle (iptr(4),iptr(5),iptr(6),iptr(7),xyz)
        nchi = nchi + 1
c
        if (resnam(atmptr(1,i)) .eq. 'GLU') then
          lswap = (abs(chi(3)) .gt. 90.0)
          xbad = chi(3)
        end if
c
        if (iptr(8) .le. 0) goto 90
        chi(4) = tangle (iptr(5),iptr(6),iptr(7),iptr(8),xyz)
        nchi = nchi + 1
c
        if (iptr(9) .le. 0) goto 90
        chi(5) = tangle (iptr(6),iptr(7),iptr(8),iptr(9),xyz)
        nchi = nchi + 1
c
        if (resnam(atmptr(1,i)) .eq. 'ARG') then
          lswap = (abs(chi(5)) .gt. 90.0)
          xbad = chi(5)
        end if
c
   90   continue
        nok = nok + 1
        lp = .false.
c
        if (.not. llaa) then
          if (.not. lp) call print_res (i,1)
          lp = .true.
          call prompt (' Warning - D-amino acid !')
          ndaa = ndaa + 1
        end if
c
        if (lbadc) then
          if (.not. lp) call print_res (i,1)
          lp = .true.
          write (line,'(a,f8.1)')
     +      ' Warning - poor CA chirality improper = ',chir
          call prompt (line)
          nbad = nbad + 1
        end if
c
        if (lswap) then
          if (.not. lp) call print_res (i,1)
          lp = .true.
          write (line,'(a,f8.1)')
     +      ' Warning - sidechain 1/2 atom names incorrect; Chi = ',
     +      xbad
          call prompt (line)
          nswap = nswap + 1
        end if
c
c ... residue-specific checks can go here
c
c ... Ile - CB chirality
c
        if (resnam(atmptr(1,i)) .eq. 'ILE') then
          call get_geom (i,' CB ',' CG1',' CG2',' CA ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' ERROR - wrong CB chirality; improper = ',chir
              call prompt (line)
            end if
          end if
c
c ... Leu, Val, Thr - CG / CB tetrahedral
c
        else if (resnam(atmptr(1,i)) .eq. 'LEU') then
          call get_geom (i,' CG ',' CD2',' CD1',' CB ',chir,ierr)
ccc      print *,'LEU ',chir,ierr
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG tetrahedral; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'VAL') then
          call get_geom (i,' CB ',' CG2',' CG1',' CA ',chir,ierr)
ccc      print *,'VAL ',chir,ierr
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CB tetrahedral; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'THR') then
          call get_geom (i,' CB ',' OG1',' CG2',' CA ',chir,ierr)
ccc      print *,'THR ',chir,ierr
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CB tetrahedral; improper = ',chir
              call prompt (line)
            end if
          end if
        end if
c
c ... sp2-carbon flatness
c
        if (resnam(atmptr(1,i)) .eq. 'ARG') then
          call get_geom (i,' CZ ',' NH1',' NH2',' NE ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CZ flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'ASN') then
          call get_geom (i,' CG ',' OD1',' ND2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'ASP') then
          call get_geom (i,' CG ',' OD1',' OD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'GLN') then
          call get_geom (i,' CD ',' OE1',' NE2',' CG ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CD flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'GLU') then
          call get_geom (i,' CD ',' OE1',' OE2',' CG ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CD flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'HIS') then
          call get_geom (i,' CG ',' ND1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'PHE') then
          call get_geom (i,' CG ',' CD1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'TYR') then
          call get_geom (i,' CG ',' CD1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'TRP') then
          call get_geom (i,' CG ',' CD1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        end if
c
c ... ring flatness here ?
c
c
c ... favourable chi1-4 torsion angles here ?
c
c
c ... favourable chi1-2 rotamers here ?
c
c
c ... add to text file
c
        if (mode .eq. 1) then
          kk = atmptr(1,i)
          if (nchi .le. 0) then
            write (iunit,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),chir
          else
            if (resnam(atmptr(1,i)) .eq. 'TYR' .or.
     +          resnam(atmptr(1,i)) .eq. 'PHE' .or.
     +          resnam(atmptr(1,i)) .eq. 'HIS' .or.
     +          resnam(atmptr(1,i)) .eq. 'TRP') then
              nchi = min (2, nchi)
            end if
            write (iunit,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),chir,(chi(j),j=1,nchi)
          end if
        end if
c
  100   continue
      end do
c
 8000 continue
      write (*,*)
      call jvalut (' Nr of residues checked    :',1,nok)
      call jvalut (' Nr of D-amino acids       :',1,ndaa)
      call jvalut (' Nr with poor CA chirality :',1,nbad)
      call jvalut (' Nr wrong 1/2 atom names   :',1,nswap)
c
      if (mode .eq. 1) close (iunit)
c
      return
c
 9000 continue
      call errcon ('While writing file')
      return
c
      end
c
c
c
      subroutine get_geom (ires,a1,a2,a3,a4,val,ierr)
c
      include 'moleman2.incl'
c
      real dist,angle,tangle,val
c
      integer ip(4),ires,ierr,ndo,nok,i,j
c
      character*4 a1,a2,a3,a4,ats(4)
c
code ...
c
      ierr = -1
      if (a1 .eq. ' ' .or. a2 .eq. ' ') return
c
      ats(1) = a1
      ats(2) = a2
      ats(3) = a3
      ats(4) = a4
c
      ndo = 4
      do i=1,4
        if (ats(i) .eq. ' ') then
          ndo = i-1
          goto 10
        end if
      end do
c
   10 continue
      if (ndo .lt. 2) return
c
      nok = 0
      do i=1,ndo
        ip(i) = -1
        do j=atmptr(1,ires),atmptr(2,ires)
          if (atmnam(j) .eq. ats(i)) then
            ip(i) = j
            nok = nok + 1
            goto 20
          end if
        end do
   20   continue
      end do
c
      if (nok .ne. ndo) return
c
      if (ndo .eq. 2) then
        val = dist (ip(1),ip(2),xyz)
      else if (ndo .eq. 3) then
        val = angle (ip(1),ip(2),ip(3),xyz)
      else if (ndo .eq. 4) then
        val = tangle (ip(1),ip(2),ip(3),ip(4),xyz)
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine ca_check (iunit,file,what,which)
c
      include 'moleman2.incl'
c
      real dist,angle,tangle,dx,dy,ptsize,phi,psi,cadist,cators
      real cangle,ave,sdv,xmin,xmax,xtot,x,y
c
      integer caqual(0:60,-60:60)
      integer lstat(maxres),dcnts(5),atcnts(0:3)
      integer i,nok,j,iunit,mode,length,ierr,kk,i1,natok,i2,ii
      integer ncheck,niso,leng1
c
      logical lp,label,lend,lgly,lbadd,lbadat,lstart
c
      character*(*) file,what,which
      character line*256,labx*40,laby*40,mychn*1
c
code ...
c
      call initca (caqual(0,-60))
c
c ... bug fix
c
      do i=0,60
        caqual (i,-60) = max (caqual (i,-59),caqual (i,-60),
     +                        caqual (i,60))
      end do
c
      do i=1,5
        dcnts(i) = 0
      end do
c
      do i=0,3
        atcnts(i) = 0
      end do
c
      mode = 0
      if (length(file) .gt. 0) then
        call upcase (what)
        call remspa (what)
        if (what(1:1) .eq. 'T') then
          mode = 1
        else if (what(1:1) .eq. 'L') then
          mode = 2
          label = .true.
        else
          mode = 2
          label = .false.
        end if
      end if
c
      mychn = which(1:1)
      call upcase (mychn)
c
      if (mode .eq. 1) then
c
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
c
      else if (mode .eq. 2) then
c
        dx = 1.0
        dy = 2.0
        ptsize = 10.0
c
        call xps_init ()
        call xps_open (iunit,file,'MOLEMAN2')
        call xps_scale (0.0,180.0,-180.0,180.0)
        call xps_stroke ()
        call xps_ps_comment ('Set up CA-Ramachandran plot')
c
        do i=0,60
          do j=-60,60
            phi = float(i)*3.0
            psi = float(j)*3.0
            if (caqual(i,j) .eq. 1) then
              call xps_dark_box (phi,phi+3.0,psi,psi+3.0)
            else if (caqual(i,j) .eq. 2) then
              call xps_grey_box (phi,phi+3.0,psi,psi+3.0)
            else if (caqual(i,j) .eq. 3) then
              call xps_light_box (phi,phi+3.0,psi,psi+3.0)
            end if
          end do
        end do
c
        call xps_move (  0.,-180.)
        call xps_draw (180.,-180.)
        call xps_draw (180., 180.)
        call xps_draw (  0., 180.)
        call xps_draw (  0.,-180.)
c
c ---	phi axis ticks
c
        do 300 i=1,5
          phi = i*30.
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
          call xps_move (  0.,psi)
          call xps_draw (  3.,psi)
          call xps_move (180.,psi)
          call xps_draw (177.,psi)
310     continue
c
        labx = 'CA-CA-CA angle mapped to [0,180>'
        laby = 'CA-CA-CA-CA dihedral mapped to [-180,180>'
        call xps_label (labx,laby)
        call xps_ps_comment ('CA-Ramachandran initialisation done')
c
      end if
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,4(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,4(1x,f10.2))
 6010 format (/
     +  ' Listing for PDB file ',a,' chain ',a1/
     + /' CA(i)   - CA(i+1) distance (2.9=cis, 3.8=trans)'/
     +  ' CA(i-1) - CA(i) - CA(i+1) angle'/
     +  ' CA(i-1) - CA(i) - CA(i+1) - CA(i+2) pseudo-torsion angle'/
     + /' An entry of "-999.9" means that the value could not be'/
     +  ' calculated (for near-terminal residues)'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil)),mychn
        write (iunit,6100,err=9000)
     +    'Residue','Distance','Angle','Torsion'
        write (iunit,6100,err=9000)
     +    '-------','--------','-----','-------'
      end if
c
      nok = 0
      natok = 0
      ncheck = 0
      niso = 0
c
      do i=1,nres
c
        lstat (i) = 0
c
        if (mychn .ne. '*') then
          if (achain(atmptr(1,i)) .ne. mychn) then
            lstat (i) = -10
            goto 1000
          end if
        end if
c
        if (captr(i) .le. 0) then
          lstat (i) = -10
          goto 1000
        end if
c
        lstart = .false.
        lend = .false.
        lgly = (resnam(captr(i)) .eq. 'GLY')
        lbadd = .false.
        lbadat = .false.
c
        cadist = -999.9
        cangle = -999.9
        cators = -999.9
c
        if (i .eq. 1) then
          lstart = .true.
        else if (captr(i-1) .gt. 0) then
          lstart = (dist(captr(i-1),captr(i),xyz) .gt. mxcaca)
        else
          lstart = .true.
        end if
c
        if (i .eq. nres) then
          lend = .true.
        else if (captr(i+1) .gt. 0) then
          lend = (dist(captr(i+1),captr(i),xyz) .gt. mxcaca)
        else
          lend = .true.
        end if
c
        if (lend) then
c
          lstat (i) = -20
c
        else
c
          nok = nok + 1
          cadist = dist(captr(i+1),captr(i),xyz)
          rbuf (nok) = cadist
          if (cadist .le. 2.80) then
            dcnts(1)=dcnts(1)+1
            lbadd = .true.
          else if (cadist .le. 3.00) then
            dcnts(2)=dcnts(2)+1
          else if (cadist .le. 3.70) then
            dcnts(3)=dcnts(3)+1
            lbadd = .true.
          else if (cadist .le. 3.90) then
            dcnts(4)=dcnts(4)+1
          else
            dcnts(5)=dcnts(5)+1
            lbadd = .true.
          end if
c
          if (lstart) goto 100
c
          if (captr(i+2) .le. 0) goto 100
          if (dist(captr(i+1),captr(i+2),xyz) .gt. mxcaca) goto 100
c
          cangle = angle (captr(i-1),captr(i),captr(i+1),xyz)
c
          if (captr(i+3) .le. 0) goto 100
          if (dist(captr(i+2),captr(i+3),xyz) .gt. mxcaca) goto 100
c
          natok = natok + 1
          cators = tangle (captr(i-1),captr(i),captr(i+1),
     +                     captr(i+2),xyz)
c
          if (.not. lgly) then
            ncheck = ncheck + 1
            i1 = max (0, min (60, nint (cangle/3.0) ) )
            i2 = max (-60, min (60, nint (cators/3.0) ) )
            ii = caqual (i1,i2)
            atcnts(ii) = atcnts(ii) + 1
            lbadat = (ii .eq. 0 .and. (.not. lgly))
          end if
c
          if (mode .eq. 2) then
            if (lgly) then
              call xps_colour (0)
              call xps_symbol (0,cangle-dx,cangle+dx,
     +          cators-dy,cators+dy)
            else
              if (lbadat) then
                call xps_colour (4)
                call xps_symbol (1,cangle-dx,cangle+dx,
     +            cators-dy,cators+dy)
                call xps_symbol (2,cangle-dx,cangle+dx,
     +            cators-dy,cators+dy)
                if (label) then
                  kk = atmptr(1,i)
                  write (line,6000) resnam(kk),achain(kk),iresid(kk),
     +              insert(kk),inote(kk)(7:10)
                  call remspa(line)
                  call xps_text (cangle+dx+0.5,cators-dy-0.5,
     +                           ptsize,line)
                end if
              else
                call xps_colour (0)
                call xps_symbol (1,cangle-dx,cangle+dx,
     +            cators-dy,cators+dy)
              end if
            end if
          end if
c
        end if
c
  100   continue
c
        lp = .false.
c
        if (lstart) then
          call print_res(i,2)
          lp = .true.
          call prompt (' <<<<< Start of new chain')
          if (mode .eq. 1) write (iunit,*,err=9000)
        end if
c
        if (lstart .and. lend) then
          call prompt (' Warning - isolated CA atom !!!')
          niso = niso + 1
        end if
c
        if (lend) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' >>>>> End of chain')
        end if
c
        if (lbadd) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          write (line,'(a,f8.2)') ' Warning - unusual CA-CA distance ',
     +      cadist
          call prompt (line)
        end if
c
        if (lbadat) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' Warning - unusual CA geometry !')
          lstat (i) = 1
        end if
c
        if (mode .eq. 1) then
          kk = atmptr(1,i)
          write (iunit,6200,err=9000)
     +      resnam(kk),achain(kk),iresid(kk),
     +      insert(kk),inote(kk)(7:10),cadist,cangle,cators
        end if
c
 1000   continue
      end do
c
      write (*,*)
      call textut (' Chain ID checked (* = all) :',mychn)
      call jvalut (' Nr of residues checked     :',1,nres)
      call jvalut (' Nr of isolated CAs         :',1,niso)
      call jvalut (' Nr with CA-CA distance     :',1,nok)
      call jvalut (' Nr with CA ang/torsion     :',1,natok)
      call jvalut ('   Ditto, non-Gly           :',1,ncheck)
c
      if (nok .gt. 2) then
        call prompt ('0CA-CA distances:')
        call xstats (rbuf,nok,ave,sdv,xmin,xmax,xtot)
        write (*,7000) nok,ave,sdv,xmin,xmax
c
        x = 100.0*float(dcnts(1))/float(nok)
        y = abs(x-0.008)/0.079
        write (*,7010) 'Short (<= 2.8 A) ',
     +      dcnts(1),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(2))/float(nok)
        y = abs (x-0.240)/0.458
        write (*,7010) 'CIS   (<= 3.0 A) ',
     +      dcnts(2),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(3))/float(nok)
        y = abs (x-1.517)/3.492
        write (*,7010) 'Poor  (<= 3.7 A) ',
     +      dcnts(3),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(4))/float(nok)
        y = abs (x-96.818)/7.044
        write (*,7010) 'TRANS (<= 3.9 A) ',
     +      dcnts(4),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(5))/float(nok)
        y = abs (x-1.416)/4.004
        write (*,7010) 'Long  (>  3.9 A) ',
     +      dcnts(5),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
      end if
c
      if (ncheck .gt. 2) then
        call prompt ('0CA geometry:')
c
        x=100.0*float(atcnts(3))/float(ncheck)
        y=abs(x-72.8)/8.9
        write (line,5010) 'CORE',atcnts(3),x,y
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        if (y .gt. 3.0) then
          line = ' >>> WARNING - 3 SIGMA deviant !'
          call prompt (line)
        end if
c
        write (line,5000) 'Additional',atcnts(2),
     +    100.0*float(atcnts(2))/float(ncheck)
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
c
        write (line,5000) 'Generous',atcnts(1),
     +    100.0*float(atcnts(1))/float(ncheck)
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
c
        x=100.0*float(atcnts(0))/float(ncheck)
        y=abs(x-3.1)/2.2
        write (line,5010) 'DISALLOWED',atcnts(0),x,y
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        if (y .gt. 3.0) then
          line = ' >>> WARNING - 3 SIGMA deviant !'
          call prompt (line)
        end if
c
        line = ' For <= 2.0 A structures :'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Core       -  7.1 % area - 72.8 % (8.9) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Additional -  5.2 % area - 12.8 % (4.0) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Generous   - 15.0 % area - 11.3 % (4.6) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Disallowed - 72.6 % area -  3.1 % (2.2) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
c
      end if
c
c ... check for stretches of poor residues
c
      if (atcnts(0) .ge. 2) then
c
        i = 0
 3241   continue
        i = i + 1
        if (i .ge. nres) goto 3242
c
c ... 5 out of 8 sequential ?
c
          if (i .le. (nres-7)) then
            j = lstat(i) + lstat(i+1) + lstat(i+2) +
     +          lstat(i+3) + lstat(i+4) + lstat(i+5) +
     +          lstat(i+6) + lstat(i+7)
            if (j .ge. 5) then
              write (*,*)
              call prompt (
     +  ' WARNING !!!  At least five of eight sequential residues')
              call prompt (
     +  '              have poor CA geometry, starting at:')
              call print_res (i,0)
              i = i + 7
              goto 3241
            end if
          end if
c
c ... 3 out of 5 sequential ?
c
          if (i .le. (nres-4)) then
            j = lstat(i) + lstat(i+1) + lstat(i+2) +
     +          lstat(i+3) + lstat(i+4)
            if (j .ge. 3) then
              write (*,*)
              call prompt (
     +  ' WARNING !!!  At least three of five sequential residues')
              call prompt (
     +  '              have poor CA geometry, starting at:')
              call print_res (i,0)
              i = i + 4
              goto 3241
            end if
          end if
c
c ... 2 out of 3 sequential ?
c
          if (i .le. (nres-2)) then
            j = lstat(i) + lstat(i+1) + lstat(i+2)
            if (j .ge. 2) then
              write (*,*)
              call prompt (
     +  ' WARNING !!!  At least two of three sequential residues')
              call prompt (
     +  '              have poor CA geometry, starting at:')
              call print_res (i,0)
              i = i + 2
              goto 3241
            end if
          end if
c
        goto 3241
c
 3242   continue
c
      end if
c
 5000 format (1x,a10,' res : ',i6,' = ',f8.2,' %')
 5010 format (1x,a10,' res : ',i6,' = ',f8.2,' % = ',f4.1,
     +  ' SIGMA from mean')
c
 7000 format (1x,i7,' CA-CA distances'/
     +  ' Average CA-CA distance = ',f8.3,' Sigma = ',f8.3/
     +  ' Minimum CA-CA distance = ',f8.3,' Maxim = ',f8.3)
 7010 format (1x,a,' CA-CA dists : ',i6,' res = ',f8.2,' % = ',
     +  f4.1,' SIGMA from mean')
c
 8000 continue
c
      if (mode .eq. 1) then
        close (iunit)
      else if (mode .eq. 2) then
        if (natok .gt. 0) then
          call xps_close ()
          call prompt (' PostScript file created')
        else
          call xps_delete ()
          call prompt (' PostScript file empty and deleted')
        end if
      end if
c
      return
c
 9000 continue
      call errcon ('While writing file')
      return
c
      end
c
c
c
      subroutine mole_stats ()
c
      include 'moleman2.incl'
c
      real xn,yn,zn,ave,sdv,xmin,xmax,rog,xtot,xrms,xhave
      real xdim,ydim,zdim
c
      integer i,nh,ns,na
c
      logical lhydro
c
code ...
c
      call jvalut (' Nr of atoms    :',1,natoms)
      call jvalut (' Nr of residues :',1,nres)
      if (natoms .le. 0) return
c
      write (*,*)
      call jvalut (' Nr of amino acid residues :',
     +  1,icnt(iprot))
      call jvalut (' Nr of nucleic acids       :',
     +  1,icnt(inucl))
      call jvalut (' Nr of waters              :',
     +  1,icnt(iwate))
      call jvalut (' Nr of metals              :',
     +  1,icnt(imeta))
      call jvalut (' Nr of inorganics          :',
     +  1,icnt(iinor))
      call jvalut (' Nr of carbohydrates       :',
     +  1,icnt(icarb))
      call jvalut (' Nr of organic compounds   :',
     +  1,icnt(iorga))
      call jvalut (' Nr of other compounds     :',
     +  1,icnt(ihete))
c
 8702 format (' ',a10,6a11)
 8704 format (' ',a10,6f11.3)
c
      nh = 0
      ns = 0
      na = 0
      do i=1,natoms
        if (select(i)) then
          ns = ns + 1
          rbuf (ns) = xyz(1,i)
          rbuf (natoms+ns) = xyz(2,i)
          rbuf (2*natoms+ns) = xyz(3,i)
          rbuf (3*natoms+ns) = batom(i)
          rbuf (4*natoms+ns) = qatom(i)
          if (lhydro(atmnam(i))) nh = nh + 1
          if (laniso(i)) na = na + 1
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of selected atoms :',1,ns)
      call jvalut ('      Ditto, hydrogen :',1,nh)
      call jvalut ('      Ditto, ANISOU   :',1,na)
      if (ns .le. 0) return
c
      write (*,*)
      write (*,8702) 'Item','Average','St.Dev','Min','Max',
     +  'RMS','Harm.ave.'
      write (*,8702) '----','-------','------','---','---',
     +  '---','---------'
c
      call xstats (rbuf(1),ns,ave,sdv,xmin,xmax,xtot)
      write (*,8704) 'X-coord',ave,sdv,xmin,xmax
      xn = ave
      xdim = xmax - xmin
      call xstats (rbuf(natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      write (*,8704) 'Y-coord',ave,sdv,xmin,xmax
      yn = ave
      ydim = xmax - xmin
      call xstats (rbuf(2*natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      write (*,8704) 'Z-coord',ave,sdv,xmin,xmax
      zn = ave
      zdim = xmax - xmin
c
      call xstats (rbuf(3*natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      call xstat2 (rbuf(3*natoms+1),ns,xrms,xhave)
      write (*,8704) 'B-factor',ave,sdv,xmin,xmax,xrms,xhave
      if (xmin .lt. 2.0) then
        call prompt (' Warning - there are B-factors < 2.0 A**2 !')
      end if
      if (xmax .gt. 100.0) then
        call prompt (' Warning - there are B-factors > 100 A**2 !')
      end if
c
      call xstats (rbuf(4*natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      call xstat2 (rbuf(4*natoms+1),ns,xrms,xhave)
      write (*,8704) 'Occpncy',ave,sdv,xmin,xmax,xrms,xhave
      if (xmin .lt. 0.0) then
        call prompt (' Warning - there are occupancies < 0 !')
      end if
      if (xmax .gt. 1.0) then
        call prompt (' Warning - there are occupancies > 1 !')
      end if
c
      rog = 0.0
      do i=1,ns
        rog = rog + (rbuf(i)-xn)**2 +
     +              (rbuf(natoms+i)-yn)**2 +
     +              (rbuf(2*natoms+i)-zn)**2
      end do
      rog = max (0.0, rog / float(ns))
      rog = sqrt (rog)
      write (*,6010) rog
 6010 format (/' The radius of gyration is ',f6.1,' A')
c
      write (*,6000) xdim,ydim,zdim
 6000 format (/' Range of X, Y, and Z coordinates: ',
     +  f6.1,' A * ',f6.1,' A * ',f6.1,' A'/
     +  ' If you have used XYz ALign_inertia_axes, these numbers'/
     +  ' give you an indication of the dimensions of the selected'/
     +  ' molecule (or set of atoms).')
c
c ... more (centre-of-mass; sum of masses ?)
c
      return
      end
c
c
c
      subroutine geom_multi (mytyp)
c
      include 'moleman2.incl'
c
      real dist,angle,tangle,xave,xsdv,xmin,xmax,xtot,cosave,sinave
c
      integer iptr(maxcpy,maxapr)
      integer i,j,nok,ne,k,l,m,ilib,nlib,nn
c
      logical lbond(maxapr,maxapr)
      logical ldum,select_res
c
      character mytyp*3,text*11
c
code ...
c
      call textut (' Multiple copy geometry for :',mytyp)
      ilib = 0
      do i=1,nmrtyp
        if (lrname(i) .eq. mytyp) then
          ilib = i
          goto 10
        end if
      end do
   10 continue
      if (ilib .le. 0) then
        call errcon ('Residue type not in library')
        return
      end if
      nlib = nmrptr(2,ilib) - nmrptr(1,ilib) + 1
      call jvalut (' Nr of atoms :',1,nlib)
      call asciut (' Atoms :',nlib,lratom(nmrptr(1,ilib)))
      if (nlib .gt. maxapr) then
        call errcon ('Residue type has too many atoms')
        call jvalut (' Maximum :',1,maxapr)
        return
      end if
      if (nlib .le. 1) then
        call errcon ('Residue type has fewer than two atoms')
        return
      end if
c
      call prompt (' Looking for selected residues ...')
      nok = 0
      do i=1,nres
        if (resnam(atmptr(1,i)) .eq. mytyp) then
          if (select_res(i)) then
            nok = nok + 1
            if (nok .gt. maxcpy) then
              call errcon ('Too many copies')
              call jvalut (' Maximum :',1,maxcpy)
              nok = maxcpy
              goto 40
            end if
            goto 20
          end if
          goto 30
c
   20     continue
          call print_res (i,0)
          do j=1,nlib
            iptr (nok,j) = 0
ccc         print *,' ... ',lratom(nmrptr(1,ilib)+j-1)
            do k=atmptr(1,i),atmptr(2,i)
              if (lratom(nmrptr(1,ilib)+j-1) .eq.
     +            atmnam(k)) then
                iptr (nok,j) = k
ccc             print *,'     ',atmnam(k)
              end if
            end do
          end do
c
   30     continue
        end if
      end do
   40 continue
c
      call jvalut (' Nr of copies found :',1,nok)
      if (nok .lt. 2) then
        call errcon ('Fewer than two copies found')
        return
      end if
c
c ... two atoms are bonded if they are in bonding distance
c     in at least one of the copies
c
      ne = 0
      do i=1,nlib-1
        do j=i+1,nlib
          ldum = .false.
          do k=1,nok
            if (iptr(k,i) .gt. 0 .and.
     +          iptr(k,j) .gt. 0) then
              if (dist(iptr(k,i),iptr(k,j),xyz) .le. mxbond) then
                ldum = .true.
                ne = ne + 1
                goto 50
              end if
            end if
          end do
   50     continue
          lbond (i,j) = ldum
          lbond (j,i) = ldum
        end do
      end do
      call jvalut (' Nr of bonds :',1,ne)
      if (ne .lt. 1) return
c
      call fvalut (' Bond distance range large if >',1,largeb)
      call fvalut (' Bond angle    range large if >',1,largea)
c
c ... do all bonds
c
      write (*,6000) mxbond
      do i=1,nlib-1
        do j=i+1,nlib
          if (.not. lbond(i,j)) goto 60
          nn = 0
          do k=1,nok
            if (iptr(k,i) .gt. 0 .and.
     +          iptr(k,j) .gt. 0) then
              nn = nn + 1
              rbuf (nn) = dist(iptr(k,i),iptr(k,j),xyz)
            end if
          end do
          if (nn .le. 0) then
            xave = 0.0
            xsdv = 0.0
            xmin = 0.0
            xmax = 0.0
          else if (nn .eq. 1) then
            xave = rbuf (1)
            xsdv = 0.0
            xmin = rbuf (1)
            xmax = rbuf (1)
          else
            call xstats (rbuf(1),nn,xave,xsdv,xmin,xmax,xtot)
          end if
          text = ' '
          if ( (xmax-xmin) .gt. largeb) text='Large range'
          write (*,6010) lratom(nmrptr(1,ilib)+i-1),
     +      lratom(nmrptr(1,ilib)+j-1),nn,xave,xsdv,xmin,xmax,text
   60     continue
        end do
      end do
c
 6000 format (/' Bonded distances with cut-off : ',f8.3,' A'/
     +         ' ==========================================')
 6010 format (1x,a4,' - ',a4,' # ',i3,' Ave, Sdv, Min, Max ',
     +  4f8.3,1x,a)
c
      if (ne .lt. 2) return
c
c ... do all angles
c
      write (*,7000)
      do i=1,nlib
        do j=1,nlib-1
          if (.not. lbond(i,j)) goto 70
          do k=j+1,nlib
            if (.not. lbond(i,k)) goto 80
            nn = 0
            do l=1,nok
              if (iptr(l,i) .gt. 0 .and.
     +            iptr(l,j) .gt. 0 .and.
     +            iptr(l,k) .gt. 0) then
                nn = nn + 1
                rbuf (nn) = angle (iptr(l,j),iptr(l,i),iptr(l,k),xyz)
                rbuf (maxatm+nn) = dist (iptr(l,j),iptr(l,k),xyz)
              end if
            end do
c
            if (nn .le. 0) then
              xave = 0.0
              xsdv = 0.0
              xmin = 0.0
              xmax = 0.0
            else if (nn .eq. 1) then
              xave = rbuf (1)
              xsdv = 0.0
              xmin = rbuf (1)
              xmax = rbuf (1)
            else
              call xstats (rbuf(1),nn,xave,xsdv,xmin,xmax,xtot)
            end if
            text = ' '
            if ( (xmax-xmin) .gt. largea) text='Large range'
            write (*,7010) lratom(nmrptr(1,ilib)+j-1),
     +        lratom(nmrptr(1,ilib)+i-1),lratom(nmrptr(1,ilib)+k-1),
     +        nn,xave,xsdv,xmin,xmax,text
c
            if (nn .le. 0) then
              xave = 0.0
              xsdv = 0.0
              xmin = 0.0
              xmax = 0.0
            else if (nn .eq. 1) then
              xave = rbuf (maxatm+1)
              xsdv = 0.0
              xmin = rbuf (maxatm+1)
              xmax = rbuf (maxatm+1)
            else
              call xstats (rbuf(maxatm+1),nn,xave,xsdv,xmin,xmax,xtot)
            end if
            write (*,7020) nn,xave,xsdv,xmin,xmax
c
   80       continue
          end do
   70     continue
        end do
      end do
c
 7000 format (/' Angles and 1-3 angle distances'/
     +         ' ==============================')
 7010 format (1x,a4,' - ',a4,' - ',a4,' Angle     : ',i3,4f8.2,1x,a)
 7020 format (19x,' 1-3 Dist  : ',i3,4f8.3,1x,a)
c
      if (ne .lt. 3) return
c
c ... do all dihedrals
c
      write (*,8000)
      do i=1,nlib-1
        do j=i+1,nlib
          if (.not. lbond(i,j)) goto 100
          do k=1,nlib
            if (k.eq.i .or. k.eq.j) goto 110
            if (.not. lbond(k,i)) goto 110
            do l=1,nlib
              if (l.eq.i .or. l.eq.j .or. l.eq.k) goto 120
              if (.not. lbond(l,j)) goto 120
              nn = 0
              do m=1,nok
                if (iptr(m,i) .gt. 0 .and.
     +              iptr(m,j) .gt. 0 .and.
     +              iptr(m,k) .gt. 0 .and.
     +              iptr(m,l) .gt. 0) then
                  nn = nn + 1
                  rbuf (nn) = tangle (iptr(m,k),iptr(m,i),
     +                        iptr(m,j),iptr(m,l),xyz)
c
ccc                  if (rbuf(nn) .le. -140.) rbuf(nn)=rbuf(nn)+360.
ccc                  if (rbuf(nn) .gt.  240.) rbuf(nn)=rbuf(nn)-360.
c
c  (The proper way to average angles is to calculate ATAN2 ( <SIN>, <COS> )
c
                  rbuf (maxatm+nn) = dist (iptr(m,k),iptr(m,l),xyz)
                end if
              end do
c
              if (nn .le. 0) then
                xave = 0.0
                xsdv = 0.0
                xmin = 0.0
                xmax = 0.0
              else if (nn .eq. 1) then
                xave = rbuf (1)
                xsdv = 0.0
                xmin = rbuf (1)
                xmax = rbuf (1)
              else
c
                call xstats (rbuf(1),nn,xave,xsdv,xmin,xmax,xtot)
c
                cosave = 0.0
                sinave = 0.0
                do m=1,nn
                  cosave = cosave + cos (degtor*rbuf(m))
                  sinave = sinave + sin (degtor*rbuf(m))
                end do
                cosave = cosave / float (nn)
                sinave = sinave / float (nn)
                xave = rtodeg * atan2 (sinave,cosave)
              end if
c
ccc              text = ' '
ccc              if ( (xmax-xmin) .gt. larged) text='Large range'
c
              write (*,8010) lratom(nmrptr(1,ilib)+k-1),
     +          lratom(nmrptr(1,ilib)+i-1),lratom(nmrptr(1,ilib)+j-1),
     +          lratom(nmrptr(1,ilib)+l-1),nn,xave,xsdv,xmin,xmax
c
              if (nn .le. 0) then
                xave = 0.0
                xsdv = 0.0
                xmin = 0.0
                xmax = 0.0
              else if (nn .eq. 1) then
                xave = rbuf (maxatm+1)
                xsdv = 0.0
                xmin = rbuf (maxatm+1)
                xmax = rbuf (maxatm+1)
              else
                call xstats (rbuf(maxatm+1),nn,xave,xsdv,xmin,xmax,xtot)
              end if
              write (*,8020) nn,xave,xsdv,xmin,xmax
c
  120         continue
            end do
  110       continue
          end do
  100     continue
        end do
      end do
c
 8000 format (/' Dihedrals and 1-4 torsion distances'/
     +         ' ===================================')
 8010 format (1x,a4,' - ',a4,' - ',a4,' - ',a4,' Dihedral : ',
     +  i3,4f8.2,1x,a)
 8020 format (26x,' 1-4 Dist : ',i3,4f8.3,1x,a)
c
      return
      end
c
c
c
      subroutine geom_select ()
c
      include 'moleman2.incl'
c
      real dist,angle,tangle
c
      integer i,j,nok,maxp,np,nd,ne,ii,jj,k,kk,l,ll
c
code ...
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          ibuf (nok) = i
        end if
      end do
      call jvalut (' Nr of selected atoms :',1,nok)
c
      maxp = maxbuf
      np = int (sqrt(float(maxp)))
      if (nok .gt. np) then
        call errcon ('Too many atoms selected')
        call jvalut (' Maximum allowed :',1,np)
        return
      end if
c
      nd = 0
      do i=1,nok-1
        ii = ibuf(i)
        do j=i+1,nok
          jj = ibuf(j)
          nd = (i-1)*nok + j
          rbuf (nd) = dist (ii,jj,xyz)
          ne = (j-1)*nok + i
          rbuf (ne) = rbuf (nd)
        end do
      end do
c
c ... bonded distances
c
      write (*,6000) mxbond
      ne = 0
      do i=1,nok-1
        ii = ibuf(i)
        do j=i+1,nok
          jj = ibuf(j)
          nd = (i-1)*nok + j
          if (rbuf(nd) .le. mxbond) then
            ne = ne + 1
            write (*,6010) atmnam(ii),achain(ii),iresid(ii),insert(ii),
     +        atmnam(jj),achain(jj),iresid(jj),insert(jj),rbuf(nd)
          end if
        end do
      end do
      call jvalut (' Nr of bonded distances :',1,ne)
 6000 format (/' Bonded distances with cut-off : ',f8.3,' A'/
     +         ' ==========================================')
 6010 format (1x,(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' = ',f8.3,' A')
c
      if (ne .lt. 2) return
c
c ... angles
c
      write (*,7000)
      ne = 0
      do i=1,nok
        ii = ibuf (i)
        do j=1,nok-1
          jj = ibuf (j)
          nd = (i-1)*nok + j
          if (rbuf(nd) .le. mxbond .and. (ii.ne.jj)) then
            do k=j+1,nok
              kk = ibuf (k)
              if (kk .ne. ii .and. kk .ne. jj) then
                nd = (i-1)*nok + k
                if (rbuf(nd) .le. mxbond) then
                  ne = ne + 1
                  write (*,7010)
     +        atmnam(jj),achain(jj),iresid(jj),insert(jj),
     +        atmnam(ii),achain(ii),iresid(ii),insert(ii),
     +        atmnam(kk),achain(kk),iresid(kk),insert(kk),
     +        angle(jj,ii,kk,xyz),dist(jj,kk,xyz)
                end if
              end if
            end do
          end if
        end do
      end do
      call jvalut (' Nr of angles :',1,ne)
 7000 format (/' Angles and 1-3 angle distances'/
     +         ' ==============================')
 7010 format (1x,(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' - ',(a4,' [',a1,i4,a1,']'),' = ',f8.3,' deg = ',f8.3,' A')
c
c ... dihedrals
c
      write (*,8000)
      ne = 0
      do i=1,nok-1
        ii = ibuf (i)
        do j=i+1,nok
          jj = ibuf (j)
          nd = (i-1)*nok + j
          if (rbuf(nd) .le. mxbond .and. (ii.ne.jj)) then
            do k=1,nok
              kk = ibuf (k)
              if (kk .ne. ii .and. kk .ne. jj) then
                nd = (i-1)*nok + k
                if (rbuf(nd) .le. mxbond) then
                  do l=1,nok
                    ll = ibuf (l)
                    if (ll .ne. ii .and. ll .ne. jj .and.
     +                  ll .ne. kk) then
                      nd = (j-1)*nok + l
                      if (rbuf(nd) .le. mxbond) then
                        ne = ne + 1
                        write (*,8010)
     +        atmnam(kk),achain(kk),iresid(kk),insert(kk),
     +        atmnam(ii),achain(ii),iresid(ii),insert(ii),
     +        atmnam(jj),achain(jj),iresid(jj),insert(jj),
     +        atmnam(ll),achain(ll),iresid(ll),insert(ll),
     +        tangle(kk,ii,jj,ll,xyz),dist(kk,ll,xyz)
                      end if
                    end if
                  end do
                end if
              end if
            end do
          end if
        end do
      end do
      call jvalut (' Nr of dihedrals :',1,ne)
 8000 format (/' Dihedrals and 1-4 torsion distances'/
     +         ' ===================================')
 8010 format (1x,(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' - ',(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' = ',f8.3,' deg = ',f8.3,' A')
c
      return
      end
c
c
c
      logical function select_res (ires)
c
      include 'moleman2.incl'
c
      integer i,ires
c
code ...
c
      select_res = .false.
c
      if (ires .gt. 0 .and. ires .le. nres) then
        do i=atmptr(1,ires),atmptr(2,ires)
          if (select(i)) then
            select_res = .true.
            return
          end if
        end do
      end if
c
      return
      end
c
c
c
      subroutine na_duarte (iunit,file,what,which)
c
      include 'moleman2.incl'
c
      real dist,tangle,dx,dy,ptsize,eta,theta
c
      integer i,nok,j,nend,nerr,nyes,isym,icol
      integer iunit,mode,length,ierr,kk,leng1
      integer na,nc,ng,nt,nu,nx,inu
c
      logical lstart,lend,leta,ltheta,lp,label
c
      character*(*) file,what,which
      character typ1lc(6)*1
      character line*256,mychn*1
c
      data typ1lc /'A','G','C','T','U','X'/
c
code ...
c
      mode = 0
      if (length(file) .gt. 0) then
        call upcase (what)
        call remspa (what)
        if (what(1:1) .eq. 'T') then
          mode = 1
        else if (what(1:1) .eq. 'A') then
          mode = 2
          label = .true.
        else
          mode = 2
          label = .false.
        end if
      end if
c
      mychn = which(1:1)
      call upcase (mychn)
c
      if (mode .eq. 1) then
c
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
c
      else if (mode .eq. 2) then
c
        dx = 2.0
        dy = 2.0
        ptsize = 10.0
c
        call psnini (iunit,file,prognm)
c
      end if
c
      do i=1,2*nres
        ibuf (i) = -1
      end do
c
      do i=1,2*nres
        rbuf (i) = -999.99
      end do
c
      write (*,*)
      call textut (' Chain ID to check (* = all) :',mychn)
c
      nok = 0
      do i=1,nres
        lbuf (i) = (restyp(i).eq.inucl)
        if (mychn .ne. '*') then
          lbuf(i) = ( lbuf(i) .and.
     +                (achain(atmptr(1,i)) .eq. mychn) )
        end if
        if (lbuf(i)) then
          nok = nok + 1
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' P  ') then
              ibuf (i) = j
            else if (atmnam(j) .eq. ' C4*') then
              ibuf (nres+i) = j
            else if (atmnam(j) .eq. ' C4''') then
              ibuf (nres+i) = j
            end if
          end do
        end if
      end do
c
      call jvalut (' Nr of residues to check :',1,nok)
      if (nok .lt. 3) then
        call errcon ('Fewer than 3 nucleic acid residues ...')
        call xps_delete ()
        call prompt (' PostScript file empty and deleted')
        return
      end if
c
      nend = 0
      nerr = 0
      nyes = 0
c
      na = 0
      nc = 0
      ng = 0
      nt = 0
      nu = 0
      nx = 0
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,4(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,2(1x,f10.2),' [',a1,']')
 6010 format (/
     +  ' Listing of Eta and Theta pseudo-torsion'/
     +  ' angles for PDB file ',a,' chain ',a1/
     + /' ETA   = C4*(i-1) - P(i) - C4*(i) - P(i+1)'/
     +  ' THETA = P(i) - C4*(i) - P(i+1) - C4*(i+1)'/
     + /' An entry of "-999.9" means that the pseudo-torsion could'/
     +  ' not be calculated (for terminal residues and residues with'/
     +  ' missing atoms).'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil)),mychn
        write (iunit,6100,err=9000)
     +    'Residue','Eta','Theta'
        write (iunit,6100,err=9000)
     +    '-------','---','-----'
      end if
c
      do i=1,nres
c
        if (.not. lbuf(i)) goto 1000
c
        if (ibuf(i) .le. 0) goto 1000
c
        lstart = .false.
        lend = .false.
        leta = .false.
        ltheta = .false.
c
        eta = -999.9
        theta = -999.9
c
        if (i .eq. 1) then
          lstart = .true.
        else if (ibuf(i-1) .gt. 0) then
          lstart =  (dist(ibuf(i-1),ibuf(i),xyz) .gt. mxpp)
        else
          lstart = .true.
        end if
c
        if (i .eq. nres) then
          lend = .true.
        else  if (ibuf(i+1) .gt. 0) then
          lend =  (dist(ibuf(i+1),ibuf(i),xyz) .gt. mxpp)
        else
          lend = .true.
        end if
c
        if (.not. lend) then
c
c ... THETA = P(i) - C4*(i) - P(i+1) - C4*(i+1)
c
          ltheta = (ibuf(i) .gt. 0 .and.
     +              ibuf(nres+i) .gt. 0 .and.
     +              ibuf(i+1) .gt. 0 .and.
     +              ibuf(nres+i+1) .gt. 0)
          if (ltheta) then
            theta = tangle (ibuf(i),ibuf(nres+i),
     +              ibuf(i+1),ibuf(nres+i+1),xyz)
            call fix360 (theta)
          end if
c
c ... ETA   = C4*(i-1) - P(i) - C4*(i) - P(i+1)
c
          if (.not. lstart) then
            leta = (ibuf(nres+i-1) .gt. 0 .and.
     +              ibuf(i) .gt. 0 .and.
     +              ibuf(nres+i) .gt. 0 .and.
     +              ibuf(i+1) .gt. 0)
            if (leta) then
              eta = tangle (ibuf(nres+i-1),ibuf(i),
     +              ibuf(nres+i),ibuf(i+1),xyz)
              call fix360 (eta)
            end if
          end if
c
        end if
c
        if (ltheta .and. leta) then
          rbuf (i) = eta
          rbuf (nres+i) = theta
        end if
c
        lp = .false.
c
        if (lstart) then
          call print_res(i,2)
          lp = .true.
          call prompt (' <<<<< Start of new chain')
          if (mode .eq. 1) write (iunit,*,err=9000)
        end if
c
        if (lend) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' >>>>> End of chain')
        end if
c
        if (.not. lend) then
          if (.not. ltheta) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate THETA')
          end if
          if (.not. lstart) then
            if (.not. leta) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              call prompt (' ERROR - could not calculate ETA')
            end if
          end if
        end if
c
        if (lstart .or. lend) then
          nend = nend + 1
        else if (.not. (leta .and. ltheta)) then
          nerr = nerr + 1
        else
          nyes = nyes + 1
        end if
c
        if (ltheta .and. leta) then
c
          call nuctyp (resnam(atmptr(1,i)),inu)
c
          if (inu .eq. 1) then
            isym = 0
            icol = 0
            na = na + 1
          else if (inu .eq. 3) then
            isym = 1
            icol = 1
            nc = nc + 1
          else if (inu .eq. 2) then
            isym = 2
            icol = 6
            ng = ng + 1
          else if (inu .eq. 4) then
            isym = 4
            icol = 4
            nt = nt + 1
          else if (inu .eq. 5) then
            isym = 5
            icol = 5
            nu = nu + 1
          else
            isym = 3
            icol = 2
            nx = nx + 1
            inu = 6
          end if
c
          kk = atmptr(1,i)
          write (*,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),eta,theta,typ1lc(inu)
          if (mode .eq. 0) goto 1000
c
          if (mode .eq. 1) then
            write (iunit,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),eta,theta,typ1lc(inu)
            goto 1000
          end if
c
          if (.not. label) then
            call xps_colour (0)
            call xps_symbol (1,eta-dx,eta+dx,theta-dy,theta+dy)
          else
            call xps_colour (icol)
            call xps_symbol (isym,eta-dx,eta+dx,theta-dy,theta+dy)
          end if
        end if
c
 1000   continue
      end do
c
      write (*,*)
      call textut (' Chain ID checked (* = all)   :',mychn)
      call jvalut (' Total nr of residues checked :',1,nok)
      call jvalut (' Start/end residues           :',1,nend)
      call jvalut (' Problems (missing atoms ?)   :',1,nerr)
      call jvalut (' Remaining residues in plot   :',1,nyes)
      call jvalut (' ... Number of A :',1,na)
      call jvalut (' ... Number of C :',1,nc)
      call jvalut (' ... Number of G :',1,ng)
      call jvalut (' ... Number of T :',1,nt)
      call jvalut (' ... Number of U :',1,nu)
      call jvalut (' ... Others      :',1,nx)
c
      if (mode .eq. 2) then
        line = ' PDB file : '//pdbfil
        call xps_legend (line)
c
        write (line,*) ' Residues in plot :',nyes
        call xps_legend (line)
c
        if (label) then
          line = ' A=black box, C=red plus, G=cyan cross,'//
     +      ' T=blue diamond, U=magenta triangle, other=green Z'
          call xps_legend (line)
        end if
c
        line = ' Shaded area has ETA in [150,190]'//
     +    ' and/or THETA in [190,260]'
        call xps_legend (line)
c
        line = ' See: CM Duarte & AM Pyle (1998). '//
     +         'J. Mol. Biol. 284, 1465-1478.'
        call xps_legend (line)
c
        if (mode .eq. 1) then
          close (iunit)
        else if (mode .eq. 2) then
          if (nyes .gt. 0) then
            call xps_close ()
            call prompt (' PostScript file created')
          else
            call xps_delete ()
            call prompt (' PostScript file empty and deleted')
          end if
        end if
      end if
c
      return
c
 9000 continue
      call errcon ('While writing text file')
c
      return
      end
c
c
c
      subroutine pdb_sanity_check ()
c
      include 'moleman2.incl'
c
c ... max nr of allowed alternative conformations for any atom
c
      integer maxcon
      parameter (maxcon=20)
c
      real sumocc(maxapr)
      real distce,dd
c
      integer basptr(maxapr),numcon(maxapr),conptr(maxcon,maxapr)
      integer i,j,k,m,n,nat,nuniq
c
code ...
c
      do i=1,nres
        nat = atmptr(2,i) - atmptr(1,i) + 1
        if (nat .gt. maxapr) then
          write (*,*)
          call errcon ('Too many atoms in residue')
          call jvalut (' Maximum allowed :',1,maxapr)
          call print_res (i,1)
          goto 900
        end if
        do j=1,nat
          basptr (j) = 0
          numcon (j) = 0
          sumocc (j) = 0.0
        end do
        nuniq = 1
        basptr (1) = atmptr(1,i)
        numcon (1) = 1
        conptr (1,1) = atmptr(1,i)
        sumocc (1) = qatom (atmptr(1,i))
        do j=atmptr(1,i)+1,atmptr(2,i)
          do k=1,nuniq
            if (atmnam(j) .eq. atmnam(basptr(k))) then
              numcon (k) = numcon (k) + 1
              conptr (numcon(k),k) = j
              sumocc (k) = sumocc (k) + qatom (j)
              goto 800
            end if
          end do
          nuniq = nuniq + 1
          basptr (nuniq) = j
          numcon (nuniq) = 1
          conptr (1,nuniq) = j
          sumocc (nuniq) = qatom (j)
  800     continue
        end do
c
        do k=1,nuniq
c
c ... check if occupancies sum to 1.0 (not neccessary for atoms on symmetry
c     axes etc.)
c
          j = nint (100*sumocc(k))
          if (j .ne. 100) then
            write (*,6000) numcon(k),sumocc(k)
            do m=1,numcon(k)
              call print_atom (conptr(m,k))
            end do
          end if
c
c ... if only one alt. loc., then normally the flag should be blank
c
          if (numcon(k) .eq. 1) then
            if (altloc(basptr(k)) .ne. ' ') then
              write (*,6010) altloc(basptr(k))
              call print_atom (basptr(k))
            end if
          end if
c
c ... if more than 3 alt. loc., issue a warning
c
          if (numcon(k) .gt. 3) then
            write (*,6015) numcon(k)
            do m=1,numcon(k)
              call print_atom (conptr(m,k))
            end do
          end if
c
c ... if more than one alt. loc., none of the flags should be blank
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)
              if (altloc(conptr(m,k)) .eq. ' ') then
                write (*,6020) numcon(k)
                call print_atom (conptr(m,k))
              end if
            end do
          end if
c
c ... if more than one alt. loc., none of the flags should be duplicated
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)-1
              do n=m+1,numcon(k)
                if (altloc(conptr(m,k)) .eq. altloc(conptr(n,k))) then
                  write (*,6030) altloc(conptr(m,k))
                  call print_atom (conptr(m,k))
                  call print_atom (conptr(n,k))
                end if
              end do
            end do
          end if
c
c ... if more than one alt. loc., the flags ought to be alphabetic
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)-1
              do n=m+1,numcon(k)
                if (ichar(altloc(conptr(n,k))) .lt.
     +              ichar(altloc(conptr(n,k)))) then
                  write (*,6040)
                  call print_atom (conptr(m,k))
                  call print_atom (conptr(n,k))
                end if
              end do
            end do
          end if
c
c ... if more than one alt. loc., they should not be too close in space
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)-1
              do n=m+1,numcon(k)
                dd = distce (xyz(1,conptr(m,k)),xyz(1,conptr(n,k)))
                if (dd .le. 0.05) then
                  write (*,6050) dd
                else if (dd .le. 0.2) then
                  write (*,6052) dd
                else if (dd .le. 1.0) then
                  write (*,6054) dd
                else
                  goto 700
                end if
                call print_atom (conptr(m,k))
                call print_atom (conptr(n,k))
  700           continue
              end do
            end do
          end if
c
        end do
c
  900   continue
      end do
c
 6000 format (/
     +  ' WARNING - OCCUPANCIES DO NOT SUM TO 1.00'/
     +  '           for the following atom, the occupancies of the ',i3/
     +  '           alternate locations add up to ',f5.2,' instead of'/
     +  '           1.00 (this can be okay if you are sure that the'/
     +  '           atom has partial occupancy, or if it lies in a'/
     +  '           special position, such as on a twofold axis):')
c
 6010 format (/
     +  ' WARNING - SINGLE LOCATION WITH NON-BLANK LABEL'/
     +  '           the following atom has only one location, but its'/
     +  '           alternate location label is "',a1,'" instead of'/
     +  '           blank (this can be okay in some cases, e.g. for'/
     +  '           waters that only interact with one of a set of'/
     +  '           alternative conformations of an amino acid'/
     +  '           residue):')
c
 6015 format (/
     +  ' WARNING - MORE THAN 3 ALTERNATE LOCATIONS'/
     +  '           the following atom has ',i3,' alternate locations'/
     +  '           (you may want to verify that these are supported'/
     +  '           by the electron density):')
c
 6020 format (/
     +  ' ERROR   - ALTERNATE LOCATION WITH BLANK LABEL'/
     +  '           the following atom occupies one of ',i3/
     +  '           alternate locations and should therefore have a'/
     +  '           non-blank label):')
c
 6030 format (/
     +  ' ERROR   - ALTERNATE LOCATIONS WITH IDENTICAL LABELS'/
     +  '           the following atoms have identical alternate'/
     +  '           location labels ("',a1,'"):')
c
 6040 format (/
     +  ' NOTE    - LABELS NOT IN ALPHABETICAL ORDER'/
     +  '           the following two atoms have alternate location'/
     +  '           labels that are not in alphabetical order (this'/
     +  '           is not strictly necessary either, but may help'/
     +  '           avoid problems with programs that assume that'/
     +  '           labels will be called A, B, C, etc.):')
c
 6050 format (/
     +  ' ERROR   - ALTERNATE LOCATIONS IN IDENTICAL POSITIONS'/
     +  '           the following two atoms have alternate locations'/
     +  '           but their distance is only ',f5.2,' A, suggesting'/
     +  '           that they are in identical positions in space and'/
     +  '           should be merged into a single location:')
c
 6052 format (/
     +  ' WARNING - ALTERNATE LOCATIONS IN ALMOST IDENTICAL POSITIONS'/
     +  '           the following two atoms have alternate locations'/
     +  '           but their distance is only ',f5.2,' A, suggesting'/
     +  '           that they are in almost identical positions in'/
     +  '           space and could be merged into a single location:')
c
 6054 format (/
     +  ' NOTE    - ALTERNATE LOCATIONS IN SIMILAR POSITIONS'/
     +  '           the following two atoms have alternate locations'/
     +  '           and their distance is ',f5.2,' A, suggesting that'/
     +  '           they are in similar positions in space and could'/
     +  '           perhaps be merged into a single location:')
c
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      return
      end
c
c
c
      subroutine rel_cont_order (rcocut)
c
      include 'moleman2.incl'
c
      real rcocut,rcosqr,dd,rco
c
      integer i,j,k,l,nrok,ncon,sumsep
c
      logical select_res
c
code ...
c
      call prompt (
     +  ' Calculate relative contact order of selected atoms')
      call fvalut (' Cut-off contact distance (A) :',1,rcocut)
c
      rcosqr = rcocut * rcocut
      nrok = 0
      ncon = 0
      sumsep = 0
c
      do i=1,nres-1
        if (select_res(i)) then
          nrok = nrok + 1
          do j=i+1,nres
            if (select_res(j)) then
              do k=atmptr(1,i),atmptr(2,i)
                if (.not. select(k)) goto 300
                do l=atmptr(1,j),atmptr(2,j)
                  if (.not. select(l)) goto 100
                  dd = (xyz(1,k)-xyz(1,l))**2
                  if (dd .gt. rcosqr) goto 100
                  dd = dd + (xyz(2,k)-xyz(2,l))**2
                  if (dd .gt. rcosqr) goto 100
                  dd = dd + (xyz(3,k)-xyz(3,l))**2
                  if (dd .gt. rcosqr) goto 100
c
                  ncon = ncon + 1
c
c ... assume sequence separation is |res_nr(atom_l)-res_nr(atom_k)|
c
                  sumsep = sumsep + abs(iresid(l) - iresid(k))
c
                  goto 200
c
  100             continue
                end do
  300           continue
              end do
  200         continue
            end if
          end do
        end if
      end do
c
      if (select_res(nres)) nrok = nrok + 1
      call jvalut (' Nr of selected residues :',1,nrok)
      call jvalut (' Nr of contacting pairs  :',1,ncon)
      call jvalut (' Sum of separations      :',1,sumsep)
c
      rco = float(sumsep)/(float(nrok)*float(ncon))
      call fvalut (' Relative contact order  :',1,rco)
c
      return
      end
