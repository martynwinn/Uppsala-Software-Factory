      program sod
c
c ... SOD === Sequences to "O" Datablocks
c
c ... Gerard Kleywegt @ 930917
c
c ... Modified 930920
c ... Modified 940412
c
      include 'sod.incl'
c
c ... VARIABLES
c
      integer nopt,length,nline,ierr,inx,iseq,i,ilo,ihi
      integer j,nn1,nn2,nn3,nn4,leng1
c
      character line*(maxchr),optpar(maxopt)*(maxchr)
      character key*4,dumstr*(maxseq)
c
code ...
c
      call gkinit (prognm,vers)
c
      call jvalut (' Max nr of sequences ......... :',1,maxseq)
      call jvalut (' Max nr of residues .......... :',1,maxres)
      call jvalut (' Max nr of residue types ..... :',1,maxtyp)
      call jvalut (' Max length of input lines ... :',1,maxchr)
c
      nline = 0
c
      task = 'mult'
      form = 'mega'
      molnam = 'A'
      prefix = ' '
      libfil = ' '
      outfil = 'sod.odb'
      keep (1) = 'all'
      keep (2) = ' '
      keep (3) = ' '
      ofirst = '1'
      refer = '1'
      replac = 'all'
      insert = 'yes'
      delmut = 'yes'
c
      nseq = 0
      iseq = 0
      imax = 0
c
      do i=1,maxseq
        write (seqnam(i),'(a,i3)') 'sequence_',i
        call remspa (seqnam(i))
        lenseq (i) = 0
        seqnce (i) = ' '
      end do
c
      write (*,*)
      call prompt (' Reading input file ...')
c
c ... main loop
c
   10 continue
      read (*,6000,err=9990,end=9999) line
      nline = nline + 1
c
      if (length(line) .lt. 1) goto 10
      if (line(1:1) .eq. '!') goto 10
c
      do i=1,maxopt
        optpar(i) = ' '
      end do
      call extrop (line,nopt,maxopt,optpar,ierr)
c
      if (nopt .lt. 1 .or. ierr .ne. 0) then
        call errcon ('Not a valid option')
        call textut (' Line >',line)
        goto 10
      end if
c
      key = optpar(1)(1:4)
      call upcase (key)
c
 6000 format (a)
 6100 format (1x,a4,' > ',a)
 6200 format (' > ',a)
 6300 format (1x,a4,1x,i3,' > ',a)
 6400 format (1x,a4,' > ',8(a,1x))
c
      if (key .eq. 'TASK') then
        task = optpar(2)
        call upcase (task)
        write (*,6100) optpar(1),task
c
      else if (key .eq. 'FORM') then
        form = optpar(2)
        call upcase (form)
        write (*,6100) optpar(1),form
        fmtstr = optpar(3)
        call upcase (fmtstr)
        call remspa (fmtstr)
c
      else if (key .eq. 'REMA') then
        write (*,6200) line(1:leng1(line))
c
      else if (key .eq. 'MOLN') then
        molnam = optpar(2)
        call upcase (molnam)
        write (*,6100) optpar(1),molnam
c
      else if (key .eq. 'PREF') then
        prefix = optpar(2)
        call upcase (prefix)
        write (*,6100) optpar(1),prefix
c
      else if (key .eq. 'OFIR') then
        ofirst = optpar(2)
        call upcase (ofirst)
        write (*,6100) optpar(1),ofirst
c
      else if (key .eq. 'REFE') then
        refer = optpar(2)
        call upcase (refer)
        write (*,6100) optpar(1),refer
c
      else if (key .eq. 'LIBF') then
        libfil = optpar(2)
        write (*,6100) optpar(1),libfil
c
      else if (key .eq. 'REPL') then
        replac = optpar(2)
        write (*,6100) optpar(1),replac
c
      else if (key .eq. 'INSE') then
        insert = optpar(2)
        write (*,6100) optpar(1),insert
c
      else if (key .eq. 'DELE') then
        delmut = optpar(2)
        write (*,6100) optpar(1),delmut
c
      else if (key .eq. 'OUTF') then
        outfil = optpar(2)
        write (*,6100) optpar(1),outfil
c
      else if (key .eq. 'NAME') then
        call str2i (optpar(2),inx,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading NAME record')
          goto 10
        end if
        if (inx .le. 0 .or. inx .gt. maxseq) then
          call errcon ('Invalid index on NAME record')
          goto 10
        end if
        seqnam (inx) = optpar(3)
        write (*,6300) optpar(1),inx,
     +    seqnam(inx)(1:leng1(seqnam(inx)))
c
      else if (key .eq. 'KEEP') then
        do i=1,3
          keep(i) = optpar(i+1)
          call upcase (keep(i))
        end do
        write (*,6400) optpar(1),(keep(i),i=1,3)
c
      else if (key .eq. 'SEQU') then
c
        if (form .eq. 'MEGA') then
c
c ... MEGAlign format
c
          nseq = 0
          iseq = 0
          call prompt (' Reading sequences (MEGA format) ...')
c
  100     continue
          read (*,6000,err=8990,end=8999) line
          nline = nline + 1
          if (line(1:1) .eq. '!') goto 100
          if (length(line) .le. 0) then
            if (nseq .eq. 0 .and. iseq .eq. 0) then
              goto 100
            else if (nseq .eq. 0) then
              nseq = iseq
              call jvalut (' Nr of sequences detected :',1,nseq)
              iseq = 0
              goto 200
            else
              if (iseq .eq. 0) goto 100
              if (iseq .ne .nseq) then
                call errcon ('Inconsistent number of sequences')
                call jvalut (' Expected :',1,nseq)
                call jvalut (' Found    :',1,iseq)
                call errstp ('Error in input file')
              end if
              iseq = 0
            end if
  200       continue
            imax = -1
            do i=1,nseq
              imax = max (imax, length(seqnce(i)))
            end do
            do i=1,nseq
              lenseq (i) = imax
            end do
            call jvalut (' Nr of residues so far :',1,imax)
            if (imax .ge. maxres) then
              call errcon ('Max nr of residues read !!!')
            end if
            goto 100
          end if
c
          iseq = iseq + 1
          seqnce (iseq) (lenseq(iseq)+1:) = line
c
          goto 100
c
        else if (form .eq. 'EMBL') then
c
c ... EMBL format
c
          nseq = 0
          iseq = 0
          imax = 0
          call prompt (' Reading sequences (EMBL format) ...')
c
  300     continue
          read (*,6000,err=8990,end=8999) line
          nline = nline + 1
          if (line(1:1) .eq. '!') goto 300
          if (length(line) .le. 0) goto 300
c
          iseq = iseq + 1
          dumstr (1:1) = line (15:15)
          dumstr (2:)  = line (52:)
          nseq = length(dumstr)
          imax = max (imax, nseq)
          do i=1,nseq
            if (dumstr (i:i) .eq. '.' ) dumstr (i:i) = delete
            seqnce (i) (iseq:iseq) = dumstr (i:i)
          end do
          nseq = imax
          goto 300
c
        else if (form .eq. 'PIR ') then
c
c ... PIR format
c
          nseq = 0
          iseq = 0
          imax = 0
          call prompt (' Reading sequences (PIR format) ...')
c
  400     continue
          read (*,6000,err=8990,end=8999) line
          nline = nline + 1
          if (line(1:1) .eq. '!') goto 400
          if (length(line) .le. 0) goto 400
c
          read (line,*,err=8990) iseq
          if (iseq .eq. 1) nseq = nseq + 1
          dumseq = line(10:)
          call remspa (dumseq)
          seqnce (nseq) (iseq:) = dumseq
          imax = max (imax, length(seqnce(nseq)))
c
          goto 400
c
        else if (form .eq. 'EXPL') then
c
c ... EXPLicit format
c
          nseq = 0
          iseq = 0
          imax = 0
          call prompt (' Reading sequences (EXPLicit format) ...')
          call textut (' Using format >',fmtstr)
c
  500     continue
          read (*,6000,err=8990,end=8999) line
          nline = nline + 1
          if (line(1:1) .eq. '!') goto 500
          if (length(line) .le. 0) then
            iseq = 0
            goto 500
          end if
          call upcase (line)
c
          if (iseq .eq. 0) nseq = nseq + 1
          ihi = length(line)
          dumseq = ' '
          read (line,fmt=fmtstr,err=510,end=510)
     +      (dumseq(i:i),i=1,ihi)
  510     continue
c
          do i=1,ihi
            if ( (dumseq(i:i) .ge. 'A' .and. dumseq(i:i) .le. 'Z') .or.
     +           (dumseq(i:i) .eq. ' ') .or.
     +           (dumseq(i:i) .eq. delete) .or.
     +           (dumseq(i:i) .eq. nothing) ) then
              iseq = iseq + 1
              seqnce (nseq) (iseq:iseq) = dumseq(i:i)
            end if
          end do
          imax = max (imax, iseq)
c
          goto 500
c
c ... invalid format
c
        else
          call errcon ('Invalid format type')
          call prompt (' Implemented FORMats : MEGA EMBL PIR EXPL')
          call errstp ('Error in input file')
        end if
c
        goto 8999
c
c ... read error
c
 8990   continue
        call errstp ('While reading sequences ...')
c
c ... end-of-sequences
c
 8999   continue
        call jvalut (' Nr of sequences read :',1,nseq)
        imax = -1
        do i=1,nseq
          imax = max (imax, length(seqnce(i)))
        end do
        do i=1,nseq
          lenseq (i) = imax
        end do
        call jvalut (' Nr of residues total :',1,imax)
c
c ... invalid option
c
      else
        call errcon ('Invalid keyword !')
        call textut (' >',line)
      end if
c
      goto 10
c
c ... read error
c
 9990 continue
      call errstp ('While reading input file')
c
c ... end-of-file
c
 9999 continue
      write (*,*)
      call jvalut (' Nr of lines read     :',1,nline)
      write (*,*)
      call textut (' Task .................. :',task)
      call textut (' Format ................ :',form)
      call textut (' Molecule name in O .... :',molnam)
      call textut (' Residue prefix in O ... :',prefix)
      call textut (' First residue in O .... :',ofirst)
      call textut (' Library file name  .... :',libfil)
      call textut (' Output file name  ..... :',outfil)
      call asciut (' Keep mode ............. :',3,keep)
      call textut (' Reference sequence .... :',refer)
      call textut (' Mutate_replace mode ... :',replac)
      call textut (' Mutate_insert mode .... :',insert)
      call textut (' Mutate_delete mode .... :',delmut)
c
      write (*,*)
      call jvalut (' Nr of sequences read :',1,nseq)
      call jvalut (' Nr of residues read  :',1,imax)
c
      if (nseq .lt. 1) call errstp ('No sequences read')
      if (imax .lt. 1) call errstp ('Sequence(s) empty')
c
      call str2i (refer,nref,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading reference sequence number')
        call errstp ('Error in REFE record')
      end if
      if (nref .lt. 1 .or. nref .gt. nseq) then
        call errcon ('Invalid reference sequence number')
        call errstp ('Error in REFE record')
      end if
c
      write (*,*)
      if (keep(1) .eq. 'ALL ') then
c
        call prompt (' Keeping ALL residues')
c
        if (task .eq. 'HOMO') then
          imax = max(length(seqnce(1)),length(seqnce(2)))
        else
          imax = length (seqnce(nref))
        end if
c
        do i=1,nseq
          seqnce(i) = seqnce(i)(1:imax)
          lenseq(i) = imax
        end do
c
      else if (keep(1) .eq. 'RANG') then
c
        call prompt (' Keeping residues by RANGE')
        call str2i (keep(2),ilo,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading first residue number')
          call errstp ('Error in KEEP record')
        end if
        call str2i (keep(3),ihi,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading second residue number')
          call errstp ('Error in KEEP record')
        end if
        call ilohi (ilo,ihi)
        call jvalut (' First residue to keep :',1,ilo)
        call jvalut (' Last  residue to keep :',1,ihi)
        if (ilo .lt. 1 .or. ihi .gt. imax .or. ilo .ge. ihi) then
          call errcon ('Invalid range of residues')
          call errstp ('Error in KEEP record')
        end if
        imax = ihi - ilo + 1
        do i=1,nseq
          dumseq = seqnce(i)(ilo:ihi)
          seqnce(i) = dumseq
          lenseq(i) = imax
        end do
c
      else if (keep(1) .eq. 'MARK') then
c
        call prompt (' Keeping residues by MARKER')
        call remspa(keep(2))
        mark = keep(2)(1:1)
        call textut (' Marker :',mark)
        ilo = 0
        ihi = 0
        iseq = 0
        call prompt (' Searching for marker')
        do i=1,nseq
          inx = index (seqnce(i),mark)
          if (inx .gt. 0) then
            call jvalut (' Found first marker in sequence :',1,i)
            ilo = inx
            inx = index (seqnce(i)(ilo+1:),mark)
            if (inx .gt. 0) then
              call prompt (' Found second marker as well')
              ihi = ilo + inx
              iseq = i
              goto 1100
            else
              call errcon ('Second marker not found')
            endif
          end if
        end do
        call errcon (' Marker not found; keep ALL residues !')
        goto 1110
c
 1100   continue
        inx = 0
        imax = ihi - ilo + 1
        if (nref .eq. iseq) then
          call errcon ('Reference sequence IS dummy marker')
          call errstp ('Error in REFE record')
        end if
c
        do i=1,nseq
          if (i .ne. iseq) then
            inx = inx + 1
            seqnce(inx) = seqnce(i)(ilo:ihi)
            seqnam(inx) = seqnam(i)
            lenseq(inx) = imax
            if (nref .eq. i) nref = inx
          end if
        end do
        nseq = inx
        call jvalut (' Nr of sequences left   :',1,nseq)
        call jvalut (' Index of first residue :',1,ilo)
        call jvalut (' Index of last  residue :',1,ihi)
        call jvalut (' Number of residues now :',1,imax)
        if (task .eq. 'HOMO') then
          call jvalut (' Reference sequence nr  :',1,1)
          call jvalut (' Comparison sequence nr :',1,2)
          nref = 1
        else
          call jvalut (' Reference sequence nr  :',1,nref)
        end if
c
 1110   continue
c
      else
        call errcon ('Invalid KEEP option; keep ALL residues')
c
        if (task .eq. 'HOMO') then
          imax = max(length(seqnce(1)),length(seqnce(2)))
        else
          imax = length (seqnce(nref))
        end if
c
        do i=1,nseq
          seqnce(i) = seqnce(i)(1:imax)
          lenseq(i) = imax
        end do
      end if
c
c ... replace non-existent residues at the termini by '+'
c     count nr of deletions & print some info
c
      write (*,*)
      write (*,6900) 'Sequence name    ','!N-term','!C-term',
     +  'Ndel','Nres'
      write (*,6900) '=============    ','=======','=======',
     +  '====','===='
c
 6900 format (1x,a20,4(1x,a8))
 6910 format (1x,a20,4(1x,i8))
c
      do i=1,nseq
c
        nn1 = 0
        do j=1,lenseq(i)
          if (seqnce(i)(j:j) .eq. ' ') then
            seqnce(i)(j:j) = nothing
          else
            nn1 = j - 1
            goto 9910
          end if
        end do
c
 9910   continue
        do j=lenseq(i),1,-1
          if (seqnce(i)(j:j) .eq. ' ') then
            seqnce(i)(j:j) = nothing
          else
            nn2 = lenseq(i) - j
            goto 9920
          end if
        end do
c
 9920   continue
        do j=1,lenseq(i)
          if (seqnce(i)(j:j) .eq. ' ') then
            seqnce(i)(j:j) = delete
          end if
        end do
c
        nn3 = 0
        do j=1,lenseq(i)
          if (seqnce(i)(j:j) .eq. delete) nn3 = nn3 + 1
        end do
c
        nn4 = lenseq(i) - nn1 - nn2 - nn3
        write (*,6910) seqnam(i),nn1,nn2,nn3,nn4
c
      end do
c
c ... if task = HOMOlogy model, then skip this
c
      if (task .ne. 'HOMO') then
c
        write (*,*)
        call prompt (' Removing deletions in reference sequence ...')
 9030   continue
        inx = index (seqnce(nref), delete)
ccc        print *,' INX = ',inx
        if (inx .le. 0) goto 9040
        do i=1,nseq
c
c ... keep track of insertions in other sequences
c
          if (i .ne. nref) then
            if (seqnce(i)(inx:inx) .ne. delete) then
              if (inx .gt. 1) then
                call locase (seqnce(i)(inx-1:inx-1))
              end if
              call upcase (seqnce(i)(inx+1:inx+1))
            end if
          end if
c
          if (inx .gt. 1) then
            dumseq = seqnce(i)(1:inx-1) // seqnce(i)(inx+1:)
          else
            dumseq = seqnce(i)(inx+1:)
          end if
          seqnce (i) = dumseq
        end do
        goto 9030
c
      end if
c
 9040 continue
      imax = -1
      do i=1,nseq
        imax = max (imax, length(seqnce(i)))
      end do
      do i=1,nseq
        lenseq (i) = imax
      end do
      call jvalut (' Nr of residues left :',1,imax)
c
      if (task .eq. 'MULT') then
        call domult ()
      else if (task .eq. 'INIT') then
        call doinit ()
      else if (task .eq. 'PAIR') then
        call dopair ()
      else if (task .eq. 'HOMO') then
        call dohomo ()
      else
        call errcon ('Invalid TASK !')
        call prompt (' Implemented TASKs : MULT INIT PAIR HOMO')
      end if
c
      call gkquit
c
      end
c
c
c
      subroutine domult ()
c
c ... domult - analyse multiple aligned sequences
c
      include 'sod.incl'
c
      real cons(maxres),entropy(maxres)
      real l102,xdum
c
      integer ndiff(maxres),counts(0:26)
      integer i,j,nn1,inx,iunit,nfirst,ierr,leng1
      integer length,ip,np,k,ia
c
      character temp1*(maxres),temp2*(maxres),refres*1
      character resnam*6,odbnam*40,line*(maxchr)
      character form1*10,form2*10,form3*10
c
      data iunit /11/
      data form1 / '(13f6.2)' /
      data form2 / '(26i3)' /
      data form3 / '(13f6.3)' /
c
code ...
c
      write (*,*)
      call prompt (' Starting task MULT ...')
c
      if (nseq .lt. 2) call errstp ('Need at least TWO sequences')
c
c ... convert sequences to upper case
c
      do i=1,nseq
        call upcase (seqnce(i))
      end do
c
      call str2i (ofirst,nfirst,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading first O residue number')
        call errstp ('Error in OFIR record')
      end if
c
      call xopxua (iunit,outfil,.false.,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening output ODB file')
      end if
c
      call stamp (line)
      write (iunit,'(a,a)',err=9990) '! ',line(1:leng1(line))
      write (iunit,'(a)',err=9990) '!'
c
c ... do it all
c
      np = 50
      ip = np
      l102 = log10(2.0)
      ia = ichar ('A')
c
      odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_POSSIBLE'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! list of possible types for each residue'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'c',imax,'(a)'
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      do i=1,imax
c
        if (ip .eq. np) then
          write (*,6110)
          ip = 1
        else
          ip = ip + 1
        end if
c
        do j=1,nseq
          temp1(j:j) = seqnce(j)(i:i)
        end do
c
        refres = seqnce(nref)(i:i)
        temp2 = refres
        nn1 = 0
c
        do j=1,nseq
          inx = index (temp2, temp1(j:j))
          if (inx .le. 0)
     +      temp2 = temp2(1:leng1(temp2)) // temp1(j:j)
          if (temp1(j:j) .eq. refres) nn1 = nn1 + 1
        end do
        ndiff(i) = length(temp2)
        cons(i) = 100.0 * float(nn1) / float(nseq)
c
c ... calculate entropy
c
        do j=0,26
          counts (j) = 0
        end do
        do j=1,nseq
          k = ichar (temp1(j:j)) - ia + 1
          if (k .lt. 0 .or. k. gt. 26) k = 0
          counts (k) = counts (k) + 1
        end do
        entropy (i) = 0.0
ccc        print *,' COUNTS ',counts
        do j=0,26
          if (counts(j) .gt. 0) then
            xdum = float(counts(j))/float(nseq)
            entropy (i) = entropy (i) - (xdum * log10(xdum) / l102)
          end if
        end do
ccc        print *,' ENTROPY ',entropy(i)
c
        write (resnam,'(a1,i5)') prefix,(i+nfirst-1)
        call remspa (resnam)
        write (*,6100) i,resnam,refres,cons(i),entropy(i),ndiff(i),
     +    temp2(1:leng1(temp2))
c
        write (iunit,'(a)',err=9990) temp2(1:leng1(temp2))
c
      end do
      write (*,*)
c
      odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_ALIGNED'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! list of possible types for each residue'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'c',imax,'(a)'
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      do i=1,imax
        write (iunit,'(999a1)') (seqnce(j)(i:i),j=1,nseq)
      end do
c
      odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_CONSERVED'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! list of percentage sequence conservation'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'r',imax,form1
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (iunit,fmt=form1,err=9990) (cons(i),i=1,imax)
c
      odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_VARIATION'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! nr of possible types for each residue'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'i',imax,form2
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (iunit,fmt=form2,err=9990) (ndiff(i),i=1,imax)
c
      odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_ENTROPY'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! entropy for each residue (alignment position)'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'r',imax,form3
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (iunit,fmt=form3,err=9990) (entropy(i),i=1,imax)
c
      odbnam = '.ID_SOD'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! replaces your .id_template'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! reset with: copy_db .id_template .id_old'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,i8)') odbnam,'t',5,40
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (iunit,'(a)',err=9990) '%RESNAM'
      write (iunit,'(a)',err=9990) 'residue_conserved'
      write (iunit,'(a)',err=9990) 'residue_variation'
      write (iunit,'(a)',err=9990) 'residue_entropy'
      write (iunit,'(a)',err=9990) 'residue_possible'
c
      odbnam = '@' // molnam(1:leng1(molnam)) // '_SOD'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! macro ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! does the work for you'
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (line,'(a,1x,a,1x,i8,1x,i8)') odbnam,'t',17,80
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 1
c
      write (line,*) 'mol ',molnam,' delete cons vari grad ;'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 2
c
      write (line,*) 'paint_ramp RESIDUE_CONSERVED ; red blue'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 3
c
      write (line,*) 'object cons ca ; end'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 4
c
      write (line,*) 'paint_ramp RESIDUE_VARIATION ; blue red'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 5
c
      write (line,*) 'object vari ca ; end'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 16
c
      write (line,*) 'paint_ramp RESIDUE_ENTROPY ; blue red'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 17
c
      write (line,*) 'object entr ca ; end'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 6
c
      write (line,*) 'paint_prop RESIDUE_CONSERVED > -1.0 red'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 7
c
      write (line,*) 'paint_prop RESIDUE_CONSERVED > 20.0 orange'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 8
c
      write (line,*) 'paint_prop RESIDUE_CONSERVED > 40.0 green'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 9
c
      write (line,*) 'paint_prop RESIDUE_CONSERVED > 60.0 steel_blue'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 10
c
      write (line,*) 'paint_prop RESIDUE_CONSERVED > 80.0 blue'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 11
c
      write (line,*) 'object grad ca ; end'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 12
c
      write (temp1,'(a1,i5)') prefix,nfirst
      call remspa (temp1)
      write (temp2,'(a1,i5)') prefix,(nfirst+imax-1)
      call remspa (temp2)
      write (line,*) 'centre_zone ',molnam(1:leng1(molnam)),
     +  ' ',temp1(1:leng1(temp1)),' ',temp2(1:leng1(temp2))
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 13
c
      write (line,*) 'copy_db .id_old .id_template'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 14
c
      write (line,*) 'copy_db .id_template .id_sod'
      call pretty (line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
c line 15
c
      write (iunit,'(a)',err=9990) 'bell message Done'
c
c ... footer
c
      write (iunit,'(a)',err=9990) '!'
      write (iunit,'(a)',err=9990) '! File read OK'
      write (iunit,'(a)',err=9990) '!'
c
 6100 format (1x,i6,1x,a6,1x,a1,1x,f6.2,1x,f6.3,1x,i3,1x,a)
 6110 format (/1x,' Index',1x,'ResNam',1x,'S',1x,'% cons',
     +  1x,'Entrop',1x,' Nr',1x,'Possible'/
     +         1x,'======',1x,'======',1x,'=',1x,'======',
     +  1x,'======',1x,' ==',1x,'========')
c
      close (iunit)
c
      return
c
 9990 continue
      call errstp ('While writing ODB output file')
c
      end
c
c
c
      subroutine doinit ()
c
c ... doinit - generate residue_type datablock for sam_init_db
c
      include 'sod.incl'
c
      integer i,inx,iunit,ierr,leng1
c
      character dumtyp(maxres)*3,dummy*(maxtyp)
      character odbnam*40,line*(maxchr)
c
      data iunit /11/
c
code ...
c
      write (*,*)
      call prompt (' Starting task INIT ...')
c
c ... convert sequences to upper case
c
      do i=1,nseq
        call upcase (seqnce(i))
      end do
c
      call redlib (iunit)
c
      call xopxua (iunit,outfil,.false.,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening output ODB file')
      end if
c
      call stamp (line)
      write (iunit,'(a,a)',err=9990) '! ',line(1:leng1(line))
c
      do i=1,ntyp
        dummy (i:i) = typone(i)
      end do
c
      do i=1,imax
        inx = index (dummy,seqnce(nref)(i:i))
        if (inx .le. 0) then
          call textut (' Unknown residue :',seqnce(nref)(i:i))
          call errstp ('Unknown residue type !')
        end if
        dumtyp(i)=typnam(inx)
      end do
c
      odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_TYPE'
      call upcase (odbnam)
      call remspa (odbnam)
      write (iunit,'(a)',err=9990) '!'
      write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
      write (iunit,'(a)',err=9990) line(1:leng1(line))
      write (line,*) '! residue types for mol ',
     +  molnam(1:leng1(molnam))
c
      write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'c',imax,'(12a)'
      call pretty (line)
      call upcase (line)
      call textut (' Writing datablock :',line)
      write (iunit,'(a)',err=9990) line(1:leng1(line))
c
      write (iunit,'(12(a3,3x))',err=9990)
     +  (dumtyp(i),i=1,imax)
c
c ... footer
c
      write (iunit,'(a)',err=9990) '!'
      write (iunit,'(a)',err=9990) '! File read OK'
      write (iunit,'(a)',err=9990) '!'
c
      close (iunit)
c
      return
c
 9990 continue
      call errstp ('While writing ODB output file')
c
      end
c
c
c
      subroutine redlib (iunit)
c
c ... redlib - read library file
c
      include 'sod.incl'
c
      integer iunit,nopt,nline,i,ierr,length,leng1
c
      character line*(maxchr),optpar(maxopt)*(maxchr),key*4
c
code ...
c
      write (*,*)
      call textut (' Opening library file >',libfil)
c
      call xopxoa (iunit,libfil,.false.,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening library file')
      end if
      call prompt (' Reading library file ...')
c
      ntyp = 0
      nline = 0
c
  100 continue
      read (iunit,'(a)',err=800,end=900) line
      nline = nline + 1
      if (length(line) .lt. 1) goto 100
      if (line(1:1) .eq. '!') goto 100
c
      do i=1,maxopt
        optpar(i) = ' '
      end do
      call extrop (line,nopt,maxopt,optpar,ierr)
c
      if (nopt .lt. 1 .or. ierr .ne. 0) then
        call errcon ('Not a valid option')
        call textut (' Line >',line)
        goto 100
      end if
c
      key = optpar(1)(1:4)
      call upcase (key)
c
      if (key .eq. 'RESI') then
        ntyp = ntyp + 1
        if (ntyp .le. maxtyp) then
          typone (ntyp) = optpar(2)
          call upcase (typone (ntyp))
          typnam (ntyp) = optpar(3)
          call upcase (typnam (ntyp))
          typtxt (ntyp) = optpar(4)
        else if (ntyp .eq. (maxtyp +1) ) then
          call errcon ('Too many residue types in library')
        end if
      else if (key .eq. 'END ') then
        goto 900
      else
        call errcon ('Invalid keyword in library file')
        call textut (' >',line)
      end if
c
      goto 100
c
  800 continue
      call errstp ('While reading library file')
c
  900 continue
      close (iunit)
      call jvalut (' Nr of lines read    :',1,nline)
      call jvalut (' Nr of residue types :',1,ntyp)
c
      write (*,*)
      write (*,6000) ' Nr','Codes','Comments'
      write (*,6000) ' ==','=====','========'
c
 6000 format (1x,a3,1x,a5,1x,a)
 6010 format (1x,i3,1x,a1,1x,a3,1x,a)
c
      do i=1,ntyp
        write (*,6010) i,typone(i),typnam(i),
     +    typtxt(i)(1:leng1(typtxt(i)))
      end do
c
      write (*,*)
c
      return
      end
c
c
c
      subroutine dopair ()
c
c ... dopair - pair-wise comparison of sequences
c
      include 'sod.incl'
c
      integer i,j,iunit,ierr,leng1
      integer ncode(maxres)
c
      character odbnam*40,line*(maxchr),form1*10
c
      data iunit /11/
      data form1 / '(35i2)' /
c
code ...
c
      write (*,*)
      call prompt (' Starting task PAIR ...')
c
      if (nseq .lt. 2) call errstp ('Need at least TWO sequences')
c
      call xopxua (iunit,outfil,.false.,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening output ODB file')
      end if
c
      call stamp (line)
      write (iunit,'(a,a)',err=9990) '! ',line(1:leng1(line))
      write (iunit,'(a)',err=9990) '!'
c
      write (iunit,'(a)',err=9990) '! pair-wise sequence comparisons'
      write (iunit,'(a)',err=9990) '! codes (use with paint_case) :'
      write (iunit,'(a)',err=9990) '! 0 = identical residues'
      write (iunit,'(a)',err=9990) '! 1 = mutation'
      write (iunit,'(a)',err=9990) '! 2 = insertion in other sequence'
      write (iunit,'(a)',err=9990) '! 3 = deletion in other sequence'
      write (iunit,'(a)',err=9990) '! 4 = outside other sequence'
c
c ... do it
c
      do i=1,nseq
c
        if (i .ne. nref) then
c
          odbnam = molnam(1:leng1(molnam)) // '_RESIDUE_VS_' //
     +             seqnam(i)(1:leng1(seqnam(i)))
          call upcase (odbnam)
          call remspa (odbnam)
          write (iunit,'(a)',err=9990) '!'
          write (line,*) '! datablock ',odbnam(1:leng1(odbnam))
          write (iunit,'(a)',err=9990) line(1:leng1(line))
          write (line,*) '! pair-wise comparison with ',
     +      seqnam(i)(1:leng1(seqnam(i)))
          write (iunit,'(a)',err=9990) line(1:leng1(line))
c
          write (line,'(a,1x,a,1x,i8,1x,a)') odbnam,'i',imax,form1
          call pretty (line)
          call upcase (line)
          call textut (' Writing datablock :',line)
          write (iunit,'(a)',err=9990) line(1:leng1(line))
c
          do j=1,imax
            if (seqnce(i)(j:j) .eq. seqnce(nref)(j:j) ) then
              ncode (j) = 0
            else if (seqnce(i)(j:j) .ge. 'a' .and.
     +               seqnce(i)(j:j) .le. 'z' ) then
              ncode (j) = 2
            else if (seqnce(i)(j:j) .eq. delete) then
              ncode (j) = 3
            else if (seqnce(i)(j:j) .eq. nothing) then
              ncode (j) = 4
            else
              ncode (j) = 1
            end if
          end do
c
          write (iunit,fmt=form1,err=9990) (ncode(j),j=1,imax)
c
        end if
c
      end do
c
c ... footer
c
      write (iunit,'(a)',err=9990) '!'
      write (iunit,'(a)',err=9990) '! File read OK'
      write (iunit,'(a)',err=9990) '!'
c
      close (iunit)
c
      return
c
 9990 continue
      call errstp ('While writing ODB output file')
c
      end
c
c
c
      subroutine dohomo ()
c
c ... dohomo - build an homology model
c
      include 'sod.incl'
c
      integer i,iunit,ierr,inx,irep,jrep,inow,iins,idel
      integer length,nfirst,leng1
c
      logical dodel,doins
c
      character dummy*(maxtyp),newres*6,insres*6
      character line*(maxchr),form1*10,resnam*6
c
      data iunit /11/
      data form1 / '(35i2)' /
c
code ...
c
      write (*,*)
      call prompt (' Starting task HOMO ...')
c
      if (nseq .lt. 2) call errstp ('Need at least TWO sequences')
c
c ... convert sequences to upper case
c
      do i=1,nseq
        call upcase (seqnce(i))
      end do
c
      call redlib (iunit)
c
      do i=1,ntyp
        dummy (i:i) = typone(i)
      end do
c
      call upcase (replac)
      irep = 0
      if (replac.eq.'ALL') then
        call prompt (' All differing residues will be mutated')
      else
        irep = 0
        do i=1,ntyp
          if (replac .eq. typnam(i)) irep = i
        end do
        call textut (' Replace differing residues by :',replac)
        if (irep .le. 0) then
          call errstp (' Replace residue type unknown')
        end if
      end if
c
      call upcase (insert)
      doins =  (insert .ne. 'NO ')
      if (doins) then
        call prompt (' Will insert residues')
      else
        call prompt (' Will not insert residues')
      end if
c
      call upcase (delmut)
      dodel =  (delmut .ne. 'NO ')
      if (doins) then
        call prompt (' Will delete residues')
      else
        call prompt (' Will not delete residues')
      end if
c
      call str2i (ofirst,nfirst,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading first O residue number')
        call errstp ('Error in OFIR record')
      end if
c
      call xopxua (iunit,outfil,.false.,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening O macro file')
      end if
c
      call stamp (line)
      write (iunit,'(a,a)',err=9990) '! ',line(1:leng1(line))
      write (iunit,'(a)',err=9990) '!'
c
      idel = 0
      iins = 0
      jrep = 0
      inow = 0
      insres = ' '
c
      call textut (' Sequence 1 :',seqnam(1))
      call textut (' Sequence 1 :',seqnce(1))
      call textut (' Sequence 2 :',seqnam(2))
      call textut (' Sequence 2 :',seqnce(2))
c
      do i=1,imax
c
        if (seqnce(2)(i:i).ne.delete .and.
     +      seqnce(2)(i:i).ne.nothing) inow = inow + 1
c
        write (resnam,'(a1,i5)') prefix,(inow+nfirst-1)
        call remspa (resnam)
c
        write (newres,'(a1,i5)') 'X',(iins+1)
        call remspa (newres)
c
c ... residue deleted in our new sequence
c
        if ((seqnce(1)(i:i).eq.delete .or.
     +       seqnce(1)(i:i).eq.nothing) .and.
     +      seqnce(2)(i:i).ne.delete .and.
     +      seqnce(2)(i:i).ne.nothing) then
          write (line,*) 'mutate_delete ',molnam,' ',resnam,' ;'
          call pretty (line)
          if (dodel) then
            write (iunit,'(a)',err=9990) line(1:leng1(line))
            call textut (' DO :',line)
          else
            write (iunit,'(a,a)',err=9990) '! ',line(1:leng1(line))
            call textut (' !  :',line)
          end if
          idel = idel + 1
c
c ... residue inserted in our new sequence
c
        else if ((seqnce(2)(i:i).eq.delete .or.
     +            seqnce(2)(i:i).eq.nothing) .and.
     +           seqnce(1)(i:i).ne.delete .and.
     +           seqnce(1)(i:i).ne.nothing) then
          inx = index (dummy,seqnce(1)(i:i))
          if (inx .le. 0) then
            call textut (' Unknown residue :',seqnce(1)(i:i))
            call errstp ('Unknown residue type !')
          end if
c
          write (line,*) 'mutate_insert ',molnam,' ',insres,
     +      ' ',newres,' ',typnam(inx),' ;'
          if (length(insres) .eq. 0) then
            call errcon ('Insert N-terminus yourself !')
            write (iunit,'(a)',err=9990) 'bell message NOTE'
            write (iunit,'(a)',err=9990)
     +        'print ... Insert N-terminus yourself !!!'
            write (line,*) 'mutate_insert ',molnam,' ',resnam,
     +        ' ',newres,' ',typnam(inx),' ;'
          end if
          call pretty (line)
c
          if (doins .and. length(insres) .gt. 0) then
            write (iunit,'(a)',err=9990) line(1:leng1(line))
            call textut (' DO :',line)
            insres = newres
          else
            write (iunit,'(a,a)',err=9990) '! ',line(1:leng1(line))
            call textut (' !  :',line)
            if (length(insres) .gt. 0) insres = newres
          end if
          iins = iins + 1
c
c ... residues differ
c
        else if (seqnce(1)(i:i).ne.seqnce(2)(i:i)) then
          if (irep .lt. 1) then
            inx = index (dummy,seqnce(1)(i:i))
            if (inx .le. 0) then
              call textut (' Unknown residue :',seqnce(1)(i:i))
              call errstp ('Unknown residue type !')
            end if
            write (line,*) 'mutate_replace ',molnam,' ',resnam,
     +      ' ',typnam(inx),' ;'
            call pretty (line)
            write (iunit,'(a)',err=9990) line(1:leng1(line))
            call textut (' DO :',line)
          else
            write (line,*) 'mutate_replace ',molnam,' ',resnam,
     +      ' ',typnam(irep),' ;'
            call pretty (line)
            write (iunit,'(a)',err=9990) line(1:leng1(line))
            call textut (' DO :',line)
          end if
          jrep = jrep + 1
          insres = resnam
c
        else
          if (seqnce(2)(i:i).ne.delete .and.
     +        seqnce(2)(i:i).ne.nothing .and.
     +        seqnce(1)(i:i).ne.delete .and.
     +        seqnce(1)(i:i).ne.nothing) insres = resnam
        end if
      end do
c
c ... footer
c
      write (iunit,'(a)',err=9990) '!'
      write (iunit,'(a)',err=9990) 'bell message Done'
      write (iunit,'(a)',err=9990) '!'
c
      write (*,*)
      call ivalut (' Nr of deletions    :',1,idel)
      call ivalut (' Nr of insertions   :',1,iins)
      call ivalut (' Nr of replacements :',1,jrep)
c
      close (iunit)
c
      return
c
 9990 continue
      call errstp ('While writing O macro file')
c
      end
