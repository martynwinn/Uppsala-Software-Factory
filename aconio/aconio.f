      program aconio
c
c ... ACONIO - convert PDB files with alternative conformations
c              from and to O-type multiple PDB file (also account
c              for ANISOU cards)
c
c ... Gerard J Kleywegt @ 970818
c
c ... f77 -o ACONIO aconio.f ../gklib/6d_kleylib; strip ACONIO
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'ACONIO', vers = '081202/1.0.3')
c
c ... MAXATM = max nr of atoms per residue
c     MAXALT = max nr of alt conformations
c     MAXAAT = max nr of atoms for with alt. conf.
c     MAXEXT = max nr of "extra" records kept (CRYST, ORIGX, SCALE)
c
      integer maxatm,maxalt,maxaat,maxext
      parameter (maxatm=1000, maxalt=5, maxaat=1000000, maxext=10)
c
      real cutdis,cutocc
c
      integer numalt(maxalt)
      integer ioin,iunit,junit,ierr,nrifil,natom,i,j,k
      integer nalt,nextra,numres,natres,natlin,length,kunit
      integer naltres,naltat,naltcon
c
      logical linter,xinter,lmain,isitok,lxplor
c
      character atline(maxatm)*80,anisou(maxatm)*80
      character altatm(maxaat)*80,altani(maxaat)*80
      character infile(maxalt+1)*128,utfile(maxalt+1)*128
      character line*80,extras(maxext)*80,nowres*10,nowatm*80
      character insert*(maxalt+1),segid(maxalt+1)*4
c
code ...
c
      call gkinit (prognm,vers)
c
      linter = xinter()
      iunit = 11
      junit = 31
      kunit = 41
      cutdis = -1.0
      cutocc = 0.99
      segid (1) = ' AC1'
      segid (2) = ' AC2'
      do i=3,6
        segid (i) = ' '
      end do
c
      call jvalut (' Max nr of atoms per residue     :',1,maxatm)
      call jvalut (' Max nr of alt. confns.          :',1,maxalt)
      call jvalut (' Max nr of atoms with alt. conf. :',1,maxaat)
c
      write (*,*)
      write (*,*) 'Supported operations:'
      write (*,*) '(1) PDB, SHELX, CCP4, TNT ---> O'
      write (*,*) '(2) X-PLOR, CNS ---> O'
      write (*,*) '(3) O ---> PDB, SHELX, CCP4, TNT'
      write (*,*) '(4) O ---> X-PLOR, CNS'
      write (*,*) '(5) Re-associate ANISOU cards (quick-n-dirty)'
      write (*,*)
      ioin = 1
      call jvalin (' Operation (1,2,3,4,5) ?',1,ioin)
      if (ioin .lt. 1 .or. ioin .gt. 5) then
        call errstp ('Invalid operation')
      end if
c
      if (ioin .eq. 1. or. ioin. eq. 2) then
        write (*,*)
        infile(1) = 'm1.pdb'
        call textin (' Input PDB file name ?',infile(1))
        call xopxoa (iunit,infile(1),linter,ierr)
        if (ierr .ne. 0) then
          call errstp ('While opening input PDB file')
        end if
        nrifil = 1
        if (ioin .eq. 2) then
          write (*,*)
          call textin (' SEGID of *main* conformation ?',segid(1))
          do i=2,maxalt
            call textin (' Alt. SEGID (<CR> to end) ?',segid(i))
            if (segid(i) .eq. ' ') then
              naltcon = i - 1
              goto 11
            end if
          end do
          naltcon = maxalt
   11     continue
        end if
      else if (ioin .eq. 3 .or. ioin .eq. 4) then
        write (*,*)
        nrifil = 0
   10   continue
        infile (nrifil+1) = ' '
        call textin (' Input PDB file from O (<CR> to end) ?',
     +    infile(nrifil+1))
        if (infile(nrifil+1) .eq. ' ') goto 20
        call xopxoa (iunit+nrifil,infile(nrifil+1),linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('Could not open file; try again')
          goto 10
        end if
        nrifil = nrifil + 1
        if (nrifil .lt. maxalt) goto 10
   20   continue
        call jvalut (' Nr of PDB files from O :',1,nrifil)
        if (nrifil .lt. 1) then
          call errstp ('No O files specified')
        end if
c
        write (*,*)
        insert = 'ABCDEFGHIJKLMNOPQ'
        if (ioin .eq. 3) then
          call textin (' Alt. conf. flags ?',insert)
          call upcase (insert)
          call remspa (insert)
        else
          do i=1,nrifil
            if (i .le. 9) then
              write (segid(i),'(1x,a2,i1)') 'AC',i
            else
              write (segid(i),'(a2,i2)') 'AC',i
            end if
            write (line,'(a,i2,a)') ' SEGID # ',i,' ?'
            call textin (line,segid(i))
          end do
        end if
c
        call prompt (
     +    '0Alternative conformations for all main-chain')
        call prompt (
     +    ' atoms can be ignored if you like.')
        lmain = isitok (' Ignore alt. main-chain conformations ?')
c
        call prompt (
     +    '0An alternative conformation for an atom is')
        call prompt (
     +    ' ignored if the occupancy of its main-')
        call prompt (
     +    ' conformation partner is high.')
        call fvalin (' Cut-off occupancy ?',1,cutocc)
        cutocc = max ( 0.0, min (cutocc,9.99) )
c
        call prompt (
     +    '0An alternative conformation for an atom is')
        call prompt (
     +    ' ignored if it lies very close to its main-')
        call prompt (
     +    ' conformation partner.')
        call fvalin (' Cut-off distance ?',1,cutdis)
        cutdis = max ( -1.0 , cutdis)
c
      else if (ioin .eq. 5) then
c
        write (*,*)
        infile(1) = 'm1.pdb'
        call textin (' Current PDB file name ?',infile(1))
        call xopxoa (iunit,infile(1),linter,ierr)
        if (ierr .ne. 0) then
          call errstp ('While opening current PDB file')
        end if
c
        write (*,*)
        infile(2) = 'm2.pdb'
        call textin (' Old PDB file name with ANISOU ?',infile(2))
        call xopxoa (junit,infile(2),linter,ierr)
        if (ierr .ne. 0) then
          call errstp ('While opening old PDB file')
        end if
c
        write (*,*)
        utfile(1) = 'm3.pdb'
        call textin (' Output PDB file name ?',utfile(1))
        call xopxua (kunit,utfile(1),linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('Could not open output file')
c
c ... changed 'goto 10'
c
          goto 9100
        end if
c
      end if
c
      if (ioin .eq. 3. or. ioin. eq. 4) then
        write (*,*)
        utfile(1) = 'm2.pdb'
        call textin (' Output PDB file name ?',utfile(1))
        call xopxua (junit,utfile(1),linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('Could not open output file')
c
c ... changed 'goto 10'
c
          goto 9100
        end if
      else if (ioin .eq. 1 .or. ioin .eq. 2) then
        write (*,*)
        utfile(1) = 'm2'
        call textin (' Base name of output PDB files ?',utfile)
        nalt = maxalt + 1
        if (ioin.eq.2) nalt = naltcon
        write (*,*)
        do i=2,nalt
          write (line,'(2a,i2,a)') utfile(1)(1:length(utfile(1))),
     +      '_',i,'.pdb'
          call remspa (line)
          utfile (i) = line
          call textut (' Open PDB file :',utfile(i))
          call xopxua (junit+i-1,utfile(i),linter,ierr)
          if (ierr .ne. 0) then
            call errcon ('Could not open output file')
c
c ... changed 'goto 10'
c
            goto 9100
          end if
        end do
        utfile(1) = utfile(1)(1:length(utfile(1))) // '_1.pdb'
        call textut (' Open PDB file :',utfile(1))
        call xopxua (junit,utfile(1),linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('Could not open output file')
c
c ... changed 'goto 10'
c
          goto 9100
        end if
      end if
c
c ... read PDB file(s)
c
      natom = 0
      nalt = 0
      do i=1,maxalt
        numalt (i) = 0
      end do
      nextra = 0
      nowres = '??????????'
      numres = 0
      natres = 0
      natlin = 0
      write (*,*)
c
cATOM      1  N  AGLU     1       4.127  26.179  -7.903  0.49 57.53           N
cANISOU    1  N  AGLU     1     9336   7394   4591      4   2737   2771       N
cATOM      2  N  BGLU     1       3.535  25.488 -12.889  0.51 54.52           N
cANISOU    2  N  BGLU     1     8406   5015   6783   -887   3093    161       N
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      if (ioin .eq. 1) then
c
  100   continue
        read (iunit,'(a)',end=110,err=9000) line
  101   continue
        if (line(1:5) .eq. 'CRYST' .or.
     +      line(1:5) .eq. 'ORIGX' .or.
     +      line(1:5) .eq. 'SCALE') then
          if (nextra .lt. maxext) then
            nextra = nextra + 1
            extras (nextra) = line
          end if
          goto 100
        end if
        if (line(1:6) .ne. 'ATOM  ' .and.
     +      line(1:6) .ne. 'HETATM') goto 100
        natom = natom + 1
        if (nowres .ne. line(18:27)) then
          if (numres .gt. 0) then
            if (numres .eq. 1) then
              do i=1,nextra
                do j=1,maxalt+1
                  write (junit+j-1,'(a)',err=9100)
     +              extras(i)(1:length(extras(i)))
                end do
              end do
            end if
            call put_o_res (junit,natres,atline,anisou)
          end if
          nowres = line(18:27)
          numres = numres + 1
          natres = 0
          natlin = 0
        end if
        natres = natres + 1
        natlin = natlin + 1
        atline(natres) = line
        anisou(natres) = ' '
        read (iunit,'(a)',end=110,err=9000) line
        if (line(1:6) .eq. 'ANISOU') then
          anisou(natres) = line
          goto 100
        else
          goto 101
        end if
c
  110   continue
        call put_o_res (junit,natres,atline,anisou)
c
c ... delete unused alt. conf. files
c
        write (*,*)
        do i=1,maxalt+1
          write (junit+i-1,'(a)',err=9100) 'END'
          close (junit+i-1)
          call xopxoa (junit,utfile(i),linter,ierr)
  120     continue
          read (junit,'(a)',end=130) line
          if (line(1:6) .eq. 'ATOM  ') goto 140
          if (line(1:6) .eq. 'HETATM') goto 140
          goto 120
  130     continue
          call textut (' Delete unused file :',utfile(i))
          close (unit=junit,status='DELETE')
  140     continue
        end do
c
      else if (ioin .eq. 2) then
c
ccc        call prompt (' Not yet implemented')
c
cATOM    149  O   ARG    19      20.399  26.123  11.069  1.00 15.00          
cATOM    150  N   TYR    20      20.356  26.091  13.275  1.00 15.00          
cATOM    151  CA  TYR    20      19.127  25.254  13.216  1.00 15.00          
cATOM    152  CB  TYR    20      19.291  23.851  13.854  0.33 15.00      AC1 
cATOM    153  CG  TYR    20      20.680  23.156  13.762  0.33 15.00      AC1 
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
c ... collect all residues in segid(1) [main conf. of alt. confs.]
c
        naltres = 0
  200   continue
        read (iunit,'(a)',end=210,err=9000) line
        if (line(1:6) .ne. 'ATOM  ') goto 200
        if (line(72:75) .ne. segid(1)) goto 200
        if (naltres .gt. 0) then
          do i=1,naltres
            if (anisou(i)(1:9) .eq. line(18:26)) goto 200
          end do
        end if
        naltres = naltres + 1
        anisou (naltres) = line(18:26)
        call textut (' Residue with alt. conf. :',line(18:26))
        goto 200
c
  210   continue
        rewind (iunit)
        write (*,*)
        call jvalut (' Nr of residues with alt. confs. :',1,naltres)
        if (naltres .lt. 1) call errstp ('No alt. confs. ???')
c
c ... collect all atoms of the residues with alt. confs.
c
        nalt = 0
        naltat = 0
  220   continue
        read (iunit,'(a)',end=230,err=9000) line
        if (line(1:6) .ne. 'ATOM  ') then
          write (junit,'(a)',err=9100) line(1:length(line))
          goto 220
        end if
        do i=1,naltres
          if (anisou(i)(1:9) .eq. line(18:26)) goto 222
        end do
        write (junit,'(a)',err=9100) line(1:length(line))
        goto 220
c
  222   continue
        do i=2,naltcon
          if (line(72:75) .eq. segid(i)) then
            nalt = nalt + 1
            altatm (nalt) = line
            goto 220
          end if
        end do
        if (line(72:75) .eq. segid(1)) then
          write (junit,'(a)',err=9100) line(1:length(line))
          goto 220
        end if
c
c ... part of residue but not in alt. conf.
c
        naltat = naltat + 1
        altani (naltat) = line
        write (junit,'(a)',err=9100) line(1:length(line))
        goto 220
c
  230   continue
ccc        call jvalut (' Nr non-alt :',1,naltat)
ccc        call jvalut (' Nr alt :',1,nalt)
c
c ... now write alt. confs.
c
        write (*,*)
        do i=2,naltcon
c
          call textut (' Writing file for SEGID :',segid(i))
          j = 1
c
  240     continue
          if (j .gt. nalt) goto 250
          if (altatm(j)(72:75) .eq. segid(i)) then
            nowres = altatm (j)(18:26)
            do k=j,nalt
              if (altatm(k)(18:26) .ne. nowres(1:9) .or.
     +            altatm(k)(72:75) .ne. segid(i)) then
                j = k
                goto 245
              end if
              write (junit+i-1,'(a)',err=9100)
     +          altatm(k)(1:length(altatm(k)))
            end do
            j = nalt + 1
c
  245       continue
            if (naltat .gt. 0) then
              do k=1,naltat
                if (altani(k)(18:26) .eq. nowres(1:9)) then
                  write (junit+i-1,'(a)',err=9100)
     +              altani(k)(1:length(altani(k)))
                end if
              end do
            end if
          else
            j = j + 1
          end if
c
          goto 240
c
  250     continue
          write (junit+i-1,'(a)',err=9100) 'END   '
c
        end do
c
      else if (ioin .eq. 3 .or. ioin .eq. 4) then
c
        lxplor = (ioin .eq. 4)
c
c ... read alt. confns. first
c
        nalt = 0
        do i=2,nrifil
  300     continue
          read (iunit+i-1,'(a)',end=310,err=9000) line
  301     continue
          if (line(1:6) .ne. 'ATOM  ' .and.
     +        line(1:6) .ne. 'HETATM') goto 300
ccc          if (line(17:17) .eq. ' ') goto 300
c
          line (17:17) = insert(i:i)
          nalt = nalt + 1
          altatm(nalt) = line
          altani(nalt) = ' '
          read (iunit+i-1,'(a)',end=310,err=9000) line
          if (line(1:6) .eq. 'ANISOU') then
            line(17:17) = insert(i:i)
            altani(nalt) = line
            goto 300
          else
            goto 301
          end if
c
  310     continue
          call jvalut (' File #   :',1,i)
          call jvalut (' # AC now :',1,nalt)
        end do
c
c ... now process main file
c
  320   continue
        read (iunit,'(a)',end=330,err=9000) line
  321   continue
        if (line(1:6) .eq. 'END   ') goto 330
        if (line(1:6) .ne. 'ATOM  ' .and.
     +      line(1:6) .ne. 'HETATM') then
          write (junit,'(a)',err=9100) line(1:length(line))
          goto 320
        end if
c
        if (line(17:17) .ne. ' ') line(17:17) = insert(1:1)
ccc        write (junit,'(a)',err=9100) line(1:length(line))
        nowatm = line
        read (iunit,'(a)',end=322,err=9000) line
        if (line(1:6) .eq. 'ANISOU') then
ccc          if (line(17:17) .ne. ' ') line(17:17) = insert(1:1)
ccc          write (junit,'(a)',err=9100) line(1:length(line))
          call put_mates (junit,nowatm,nalt,altatm,altani,
     +      line,insert(1:1),cutdis,cutocc,lmain,lxplor,segid(1))
          goto 320
        else
          line = ' '
          call put_mates (junit,nowatm,nalt,altatm,altani,
     +      line,insert(1:1),cutdis,cutocc,lmain,lxplor,segid(1))
          backspace (iunit)
          goto 320
        end if
c
  322   continue
        line = ' '
        call put_mates (junit,nowatm,nalt,altatm,altani,
     +    line,insert(1:1),cutdis,cutocc,lmain,lxplor,segid(1))
c
  330   continue
c
c ... for X-PLOR/CNS, append alt. confs.
c
        if (lxplor) then
          write (*,*)
          do i=2,nrifil
            write (junit,'(a)',err=9100) 'REMARK'
            call textut (' SEGID :',segid(i))
            k = 0
            do j=1,nalt
              if (altatm(j)(17:17) .eq. insert(i:i) .and.
     +            altatm(j)(72:75) .eq. 'OKAY') then
                altatm (j) (17:17) = ' '
                altatm (j) (72:75) = segid(i)
                write (junit,'(a)',err=9100)
     +            altatm(j)(1:length(altatm(j)))
                altatm (j) = ' '
                k = k + 1
              end if
            end do
            call jvalut (' Nr of atoms :',1,k)
          end do
        end if
c
        write (junit,'(a)',err=9100) 'END   '
c
c ... check that no alt. conf. atoms are left
c
        j = 0
        do i=1,nalt
          if (altatm(i) .ne. ' ') then
            call textut (' Orphan ?',altatm(i))
            j = j + 1
          end if
        end do
        write (*,*)
        call jvalut (' Nr of "orphans" :',1,j)
        if (j .gt. 0) call errstp ('Help ! Orphans found !')
c
      else if (ioin .eq. 4) then
c
        call prompt (' Not yet implemented')
c
      else if (ioin .eq. 5) then
c
cATOM      1  N  AGLU     1       4.127  26.179  -7.903  0.49 57.53           N
cATOM      2  N  BGLU     1       3.535  25.488 -12.889  0.51 54.52           N
c...
cATOM      1  N  AGLU     1       4.127  26.179  -7.903  0.49 57.53           N  
cANISOU    1  N  AGLU     1     9336   7394   4591      4   2737   2771       N  
cATOM      2  N  BGLU     1       3.535  25.488 -12.889  0.51 54.52           N  
cANISOU    2  N  BGLU     1     8406   5015   6783   -887   3093    161       N  
c1234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7
c
        natom = 0
        call prompt (' Working ...')
c
  500   continue
        read (iunit,'(A)',end=590,err=9000) line
        if (line(1:6) .ne. 'ATOM  ' .and.
     +      line(1:6) .ne. 'HETATM') then
          write (kunit,'(A)',err=9100) line(1:length(line))
          goto 500
        end if
c
c ... found an ATOM or HETATM card
c
        nowatm = line (13:27)
        write (kunit,'(A)',err=9100) line(1:length(line))
c
  510   continue
        read (junit,'(A)',end=520,err=9000) line
        if (line(1:6) .ne. 'ANISOU') goto 510
        if (line(13:27) .ne. nowatm(1:15)) goto 510
c
        write (kunit,'(A)',err=9100) line(1:length(line))
        natom = natom + 1
c
        goto 500
c
  520   continue
        call errcon ('ANISOU card not found !')
        call textut (' For atom :',nowatm)
        call prompt (' Order/names of atoms differ ???')
        call errstp ('ANISOU card not found !')
c
  590   continue
        call jvalut (' Nr of atoms processed :',1,natom)
c
      end if
c
      goto 9900
c
c ... read error
c
 9000 continue
      call errstp ('While reading PDB file')
c
c ... write error
c
 9100 continue
      call errstp ('While writing PDB file')
c
c ... all done
c
 9900 continue
      write (*,*)
      call prompt (' All done !')
c
      call gkquit ()
      end
c
c
c
      subroutine put_o_res (junit,natres,atline,anisou)
c
      implicit none
c
      integer natres,junit,i,length,j,ialt
c
      character atline(natres)*80,anisou(natres)*80
      character insert*10,altnow*1,myline*128
c
code ...
c
cATOM      1  N  AGLU     1       4.127  26.179  -7.903  0.49 57.53           N
cANISOU    1  N  AGLU     1     9336   7394   4591      4   2737   2771       N
cATOM      2  N  BGLU     1       3.535  25.488 -12.889  0.51 54.52           N
cANISOU    2  N  BGLU     1     8406   5015   6783   -887   3093    161       N
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c
c ... (1) write atoms with blank insert to first file
c
      call textut (' Residue :',atline(1)(13:27))
      do i=1,natres
        if (atline(i)(17:17) .eq. ' ') then
          write (junit,'(a)',err=9100)
     +      atline(i)(1:length(atline(i)))
          if (anisou(i) .ne. ' ') then
            write (junit,'(a)',err=9100)
     +        anisou(i)(1:length(anisou(i)))
          end if
        end if
      end do
c
c ... (2) write atoms with first non-blank insert to first
c         file as well; if there aren't any, return
c
      insert = ' '
      do i=1,natres
        if (atline(i)(17:17) .ne. ' ') then
          insert (1:1) = atline(i)(17:17)
          goto 100
        end if
      end do
c
      return
c
  100 continue
      call textut (' Add insert :',insert(1:1))
      do i=1,natres
        if (atline(i)(17:17) .eq. insert(1:1)) then
          myline = atline(i)
          myline (17:17) = ' '
          write (junit,'(a)',err=9100)
     +      myline(1:length(myline))
          if (anisou(i) .ne. ' ') then
            myline = anisou(i)
            myline (17:17) = ' '
            write (junit,'(a)',err=9100)
     +        myline(1:length(myline))
          end if
        end if
      end do
c
c ... (3) look for other insert characters
c
      ialt = 0
      altnow = ' '
c
  200 continue
c
      do i=1,natres
        if (atline(i)(17:17) .ne. ' ') then
          if (index(insert,atline(i)(17:17)) .le. 0) then
            ialt = ialt + 1
            altnow = atline(i)(17:17)
            insert (length(insert)+1:) = atline(i)(17:17)
            goto 210
          end if
        end if
      end do
c
      return
c
  210 continue
      call textut (' Alt. conf. :',altnow)
      do i=1,natres
        if (atline(i)(17:17) .eq. altnow) then
          myline = atline(i)
          myline (17:17) = ' '
          write (junit+ialt,'(a)',err=9100)
     +      myline(1:length(myline))
          if (anisou(i) .ne. ' ') then
            myline = anisou(i)
            myline (17:17) = ' '
            write (junit+ialt,'(a)',err=9100)
     +        myline(1:length(myline))
          end if
        end if
      end do
c
c ... copy other atoms without alt. conf.
c
      do i=1,natres
        if (atline(i)(17:17) .eq. ' ') then
          do j=1,natres
            if (atline(j)(17:17) .eq. altnow) then
              if (atline(j)(13:16) .eq. atline(i)(13:16)) goto 220
            end if
          end do
          write (junit+ialt,'(a)',err=9100)
     +      atline(i)(1:length(atline(i)))
          if (anisou(i) .ne. ' ') then
            write (junit+ialt,'(a)',err=9100)
     +        anisou(i)(1:length(anisou(i)))
          end if
  220     continue
        end if
      end do
c
      goto 200
c
c ... write error
c
 9100 continue
      call errstp ('PUT_O_RES - While writing PDB file')
c
      end
c
c
c
      subroutine put_mates (junit,nowatm,nalt,altatm,altani,
     +                      anisou,ichar,cutdis,cutocc,lmain,
     +                      lxplor,segid1)
c
      implicit none
c
      real x1(3),x2(3),cutdis,distce,dis,cutocc,occ
c
      integer nalt,i,junit,length
c
      logical lmain,mainch,ldone,lxplor
c
      character altatm(nalt)*80,altani(nalt)*80
      character nowatm*80,anisou*80,ichar*1,segid1*4
c
code ...
c
cATOM      1  N  AGLU     1       4.127  26.179  -7.903  0.49 57.53           N
cANISOU    1  N  AGLU     1     9336   7394   4591      4   2737   2771       N
cATOM      2  N  BGLU     1       3.535  25.488 -12.889  0.51 54.52           N
cANISOU    2  N  BGLU     1     8406   5015   6783   -887   3093    161       N
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
ccc      if (nowatm(17:17) .eq. ' ') return
c
      if (nalt .le. 0) then
        write (junit,'(a)',err=9100) nowatm(1:length(nowatm))
        if (anisou .ne. ' ') then
          write (junit,'(a)',err=9100) anisou(1:length(anisou))
        end if
        return
      end if
c
      ldone = .false.
c
      do i=1,nalt
        if (altatm(i)(18:27) .eq. nowatm(18:27)) then
          if (altatm(i)(13:16) .eq. nowatm(13:16)) then
            write (*,*)
            call textut (' AC :',nowatm(1:66))
            call textut ('    :',altatm(i)(1:66))
c
c ... main-chain check (if required)
c
            if (lmain .and. mainch(nowatm(13:16))) then
              call prompt (' Merged; main-chain atom')
c              write (junit,'(a)',err=9100)
c     +          nowatm(1:length(nowatm))
c              if (anisou .ne. ' ') then
c                write (junit,'(a)',err=9100)
c     +            anisou(1:length(anisou))
c              end if
              altatm(i) = ' '
              goto 100
            end if
c
c ... occupancy check
c
            read (nowatm(56:60),'(f5.2)',err=9100) occ
            if (occ .ge. cutocc) then
              call fvalut (' Merged; main occupancy :',1,occ)
c              write (junit,'(a)',err=9100)
c     +          nowatm(1:length(nowatm))
c              if (anisou .ne. ' ') then
c                write (junit,'(a)',err=9100)
c     +            anisou(1:length(anisou))
c              end if
              altatm(i) = ' '
              goto 100
            end if
c
c ... distance check
c
            read (altatm(i)(31:54),'(3f8.3)',err=9100) x1
            read (nowatm(31:54),'(3f8.3)',err=9100) x2
            dis = distce(x1,x2)
            if (dis .le. cutdis) then
              call fvalut (' Merged; distance :',1,dis)
c              write (junit,'(a)',err=9100)
c     +          nowatm(1:length(nowatm))
c              if (anisou .ne. ' ') then
c                write (junit,'(a)',err=9100)
c     +            anisou(1:length(anisou))
c              end if
              altatm(i) = ' '
              goto 100
            end if
c
c ... accept as alternative conformation
c
            write (*,6045) occ,dis
 6045 format (' Alt. conf.: main occupancy ',f5.2,
     +        ' and main-alt atom distance ',f5.2,' A')
c
ccc            call fvalut (' Alt. conf.; dist :',1,dis)
c
            if (.not. ldone) then
              if (lxplor) then
                nowatm (17:17) = ' '
                nowatm (72:75) = segid1
              else
                nowatm (17:17) = ichar
              end if
              write (junit,'(a)',err=9100) nowatm(1:length(nowatm))
              if (anisou .ne. ' ') then
                anisou (17:17) = ichar
                write (junit,'(a)',err=9100) anisou(1:length(anisou))
              end if
              ldone = .true.
            end if
c
            if (lxplor) then
              altatm (i)(72:75) = 'OKAY'
            else
              write (junit,'(a)',err=9100)
     +          altatm(i)(1:length(altatm(i)))
              altatm (i) = ' '
              if (altani(i) .ne. ' ') then
                write (junit,'(a)',err=9100)
     +            altani(i)(1:length(altani(i)))
                altani (i) = ' '
              end if
            end if
c
            goto 100
c
          end if
        end if
  100   continue
      end do
c
      if (.not. ldone) then
        write (junit,'(a)',err=9100) nowatm(1:length(nowatm))
        if (anisou .ne. ' ') then
          write (junit,'(a)',err=9100) anisou(1:length(anisou))
        end if
      end if
c
      return
c
c ... write error
c
 9100 continue
      call errstp ('PUT_MATES - While writing PDB file')
c
      end
