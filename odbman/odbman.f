      program odbman
c
c ... ODBMAN - MANipulation of O DataBlocks
c
c ... Gerard Kleywegt @ 931115
c
c ... 0.1 @ 931115 - first version
c
      include 'odbman.incl'
c
      integer maxopt,maxhis
c
      parameter (maxopt = 25)
      parameter (maxhis = maxopt)
c
c ... length of a big line to hold 2*MAXODB ODBNAMs + TABs
c
      integer linbig
      parameter (linbig=2*maxodb*26)
c
      real buff1(maxsiz),buff2(maxsiz)
      real rhis(maxhis),total,user,sys,ave,sdv,xmin,xmax,xtot
      real xdum,rmsd,shap,corr,rf1,rf2,ymin,ymax,ydum,x6,x7,x8,x9
      real rdummy(maxsiz),cormin,cormax,xave,xsdv,yave,ysdv
      real rxmin,rxmax,rymin,rymax,d,zd,probd,rs,probrs
c
      integer idummy(maxsiz),nhis(maxhis+1)
      integer i,length,nempty,nopt,ierr,iunit,j,iptr,idum,jdum
      integer whichm,numhis,irf,iri,irc,k,i1,i2,i3,jptr,jmax
      integer munit,leng1,nwind,nw1,nw2,k1,k2,nhit,l
c
      logical xinter,unsave,ldone,linit,lfirst,lecho,ltable
c
      character biglin*(linbig),tab*1
      character line*256,optpar(maxopt)*256,reply*1,inimac*128
      character pro*12,prev*80,deffmt(4)*20,deftyp*4
c
      data deftyp / 'RICT' /
      data deffmt / '(1P,5E14.4)', '(6I12)', '(14A6)', '72' /
c
      equivalence (rdummy(1),idummy(1))
c
code ...
c
      call gkinit (prognm,vers)
c
c ... initialise history
c
      call dohist ('*INIT*',ldone)
c
      lecho = .false.
c
      tab = char(9)
c
      call jvalut (' Max nr of O data blocks (odb) :',1,maxodb)
      call jvalut (' Total nr of O data blocks     :',1,maxtdb)
      call jvalut (' Max nr of elements per odb    :',1,maxsiz)
      call jvalut (' Max length of text odbs       :',1,maxlen)
c
c ... define some symbols
c
      write (*,*)
      nopt = 3
      optpar(1) = '&'
      optpar(2) = ' '
      optpar(3) = ' '
c
c      optpar(2) = 'PROGRAM'
c      optpar(3) = prognm
c      call dosymb (nopt,optpar,ldone)
c      optpar(2) = 'VERSION'
c      optpar(3) = vers
c      call dosymb (nopt,optpar,ldone)
c
      optpar(2) = 'START_TIME'
      call gkdate (optpar(3)(1:24))
      call dosymb (nopt,optpar,ldone)
      optpar(2) = 'USERNAME'
      call gkuser (optpar(3))
      call dosymb (nopt,optpar,ldone)
c
c ... foreplay
c
      do i=1,maxopt
        optpar (i) = ' '
      end do
      iunit = 1
      nempty = 0
c
c ... user input unit (5=interactive; other=macro)
c
      munit = 5
c
      linit = .false.
c
      linter = xinter()
      if (linter) then
        pro='$ODBMAN > '
      else
        pro=' ODBMAN > '
      end if
c
      do i=1,maxtdb
        odbptr (i) = -1
        odblen (i) = 0
        odbuse (i) = .false.
        odbcha (i) = .false.
        odbsel (i) = .false.
        odbtyp (i) = 'c'
        odbfmt (i) = '*'
        odbnam (i) = 'mol_residue_junk'
        odbcom (i) = 'not used'
      end do
c
c ... formats
c
 6000 format (/
     +  ' ODBMAN options :'//
     +  ' ? (list options)                    ',
     +              ' ! (comment)'/
     +  ' QUit                                ',
     +              ' $ shell_command'/
     +  ' & symbol value                      ',
     +              ' & ? (list symbols)'/
     +  ' @ macro_file                        ',
     +              ' ECho on_off'/
     +  ' # parameter(s) (command history)    ',
     +              ' '/
     + /' REad filename                       ',
     +              ' WRite odb* filename'/
     +  ' DElete odb*                         ',
     +              ' TYpe odb*'/
     +  ' FOrmat RIC_odb* new_format          ',
     +              ' NAme odb new_name'/
     +  ' CReate name type nr_elements format ',
     +              ' DUplicate old_odb new_odb'/
     + /' LIst odb*                           ',
     +              ' STats RI_odb*'/
     +  ' HIsto RI_odb* x1 x2 x3 [...]        ',
     +              ' BOxcar R_odb* window_size '/
     +  ' SImilarity RI_odb1 RI_odb2* [min] [max] [how]'/
     +  ' ALl_correlations [min] [max] [how]'/
     +  ' NOn_parametric RI_odb1 RI_odb2* [min] [max] [how]'/
     + /' PLot_file RI_odb file [y_lo y_hi label_x label-y]'/
     +  ' SCatter RI_odb1 RI_odb2 file [label_x label_y]'/
     +  ' 3D_odl_file RI_odb1 RI_odb2 RI_odb3 file [colour]'/
     +  ' CGraph_file RI_odb* file'/
     + /' CHaracter_function CT_odb function  ',
     +              ' INteger_func I_odb function value'/
     +  ' FLoat_func R_odb function value     ',
     +              ' '/
     + /' SEt ALl odb value                   ',
     +              ' SEt MAny odb first last value'/
     +  ' SEt INdiv odb first last [values]   ',
     +              ' SEt ONe odb index value'/
     +  ' SEt IF odb other_odb operator value ',
     +              ' SEt KEep odb first last [step]      '/
     +  ' COmbine RI_odb1 RI_odb2 operator    ',
     +              ' '/
     + /' EXtract FIeld odb type file head_skip line_skip',
     +      ' field_nr'/
     +  ' EXtract FOrmat odb type file head_skip line_skip',
     +      ' format'/
     +  ' EXtract MUltiple odb type file head_skip format'/
     +  ' EXtract PRocheck mol_name file'/
     + /)
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6000)
c
      call jvalut (' Max nr of O data blocks (odb) :',1,maxodb)
      call jvalut (' Total nr of O data blocks     :',1,maxtdb)
      call jvalut (' Max nr of elements per odb    :',1,maxsiz)
      call jvalut (' Max length of text odbs       :',1,maxlen)
      write (*,*)
c
c --- MAIN EVENT LOOP
c
   10 continue
      call gkdcpu (total,user,sys)
      if (total .ge. 1.0)
     +  write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
c ... execute initialisation macro ?
c
      if (.not. linit) then
        linit = .true.
        call gknval ('GKODBMAN',inimac,ierr)
        if (ierr .eq. 0) then
          call textut (' Execute initialisation macro :',inimac)
          line = '@ ' // inimac
          call getcom ('*INIT*',line,munit,lecho)
          if (line(1:1) .eq. '@') goto 10
        else
          call getcom (pro,line,munit,lecho)
        end if
      else
c
c ... get next command line from terminal or macro
c
        call getcom (pro,line,munit,lecho)
c
      end if
c
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (nopt .lt. 1 .or. ierr .ne. 0) then
        call textut (' ERROR - No valid option :',line)
        goto 100
      end if
c
c ... handle symbols
c
      call dosymb (nopt,optpar,ldone)
      if (ldone) goto 10
c
      if (optpar(1)(1:1) .eq. '?') goto 100
c
      call upcase (optpar(1))
c
c ... QUIT
c
      if (optpar(1)(1:2) .eq. 'QU') then
c
        goto 9000
c
c ... ECHO
c
      else if (optpar(1)(1:2) .eq. 'EC') then
c
        if (nopt .lt. 2) then
          optpar (2) = '?'
          call textin (' Option (ON/OFf/?) ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:1) .eq. '?') then
          if (lecho) then
            call prompt (' Echo is ON')
          else
            call prompt (' Echo is OFF')
          end if
        else if (optpar(2)(1:2) .eq. 'ON') then
          lecho = .true.
        else if (optpar(2)(1:2) .eq. 'OF') then
          lecho = .false.
        else
          call errcon ('Invalid parameter')
        end if
c
        goto 10
c
c ... EXTRACT
c
      else if (optpar(1)(1:2) .eq. 'EX') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'PRocheck'
          call textin (' Option ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'FI' .or.
     +      optpar(2)(1:2) .eq. 'FO' .or.
     +      optpar(2)(1:2) .eq. 'MU') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'junk'
            call textin (' Name for new odb ?',optpar(3))
          end if
          call upcase (optpar(3))
          optpar (3) = optpar(3)(1:25)
c
          if (nopt .lt. 4) then
            optpar (4) = 'R'
            call textin (' Type (R/I) ?',optpar(4))
          end if
          call upcase (optpar(4))
          optpar (4) = optpar(4)(1:1)
          jdum = index(deftyp,optpar(4)(1:1))
          if (jdum .le. 0 .or.
     +        (optpar(4)(1:1) .ne. 'R' .and.
     +         optpar(4)(1:1) .ne. 'I') ) then
            call errcon ('Invalid datablock type')
            goto 10
          end if
c
c ... allocate space (if any left)
c
          call allocm (optpar(3),optpar(4),maxsiz,deffmt(jdum),
     +                 j,iptr,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            optpar (5) = 'junk.out'
            call textin (' File to scan ?',optpar(5))
          end if
c
          if (nopt .lt. 6) then
            optpar (6) = '0'
            call textin (' Lines to skip at top of file ?',
     +        optpar(6))
          end if
          call str2i (optpar(6),i1,ierr)
          if (ierr .ne. 0) goto 10
c
        end if
c
        if (optpar(2)(1:2) .eq. 'FI') then
c
          if (nopt .lt. 7) then
            optpar (7) = '0'
            call textin (' Lines to skip after each line ?',
     +        optpar(7))
          end if
          call str2i (optpar(7),i2,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 8) then
            optpar (8) = '1'
            call textin (' Field to read ?',optpar(8))
          end if
          call str2i (optpar(8),i3,ierr)
          if (ierr .ne. 0) goto 10
c
          call exfiel (iunit,optpar(5),rdummy,idummy,
     +                 optpar(4)(1:1),idum,i1,i2,i3,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .le. 0) then
            call errcon ('No values read')
            goto 10
          end if
          if (idum .gt. maxsiz) then
            call errcon ('Too many values')
            goto 10
          end if
c
c ... copy to datablock etc.
c
          odbuse (iptr) = .true.
          odbptr (iptr) = j
          odbcha (iptr) = .true.
          odbsel (iptr) = .false.
          odblen (iptr) = idum
          odbtyp (iptr) = optpar (4) (1:1)
          odbfmt (iptr) = deffmt (jdum)
          odbnam (iptr) = optpar (3)
          odbcom (iptr) = 'Extracted by field from '//optpar(5)
c
          if (odbtyp (iptr) .eq. 'R') then
            do i=1,idum
              rodb (i,j) = rdummy (i)
            end do
          else if (odbtyp (iptr) .eq. 'I') then
            do i=1,idum
              iodb (i,j) = idummy (i)
            end do
          end if
c
          call textut (' Extracted new odb :',odbnam(iptr))
c
        else if (optpar(2)(1:2) .eq. 'FO') then
c
          if (nopt .lt. 7) then
            optpar (7) = '0'
            call textin (' Lines to skip after each line ?',
     +        optpar(7))
          end if
          call str2i (optpar(7),i2,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 8) then
            optpar (8) = '(1x,f8.3)'
            call textin (' Format to read from line ?',optpar(8))
          end if
c
          call exform (iunit,optpar(5),rdummy,idummy,
     +                 optpar(4)(1:1),idum,i1,i2,optpar(8),ierr)
          if (ierr .ne. 0) goto 10
          if (idum .le. 0) then
            call errcon ('No values read')
            goto 10
          end if
          if (idum .gt. maxsiz) then
            call errcon ('Too many values')
            goto 10
          end if
c
c ... copy to datablock etc.
c
          odbuse (iptr) = .true.
          odbptr (iptr) = j
          odbcha (iptr) = .true.
          odbsel (iptr) = .false.
          odblen (iptr) = idum
          odbtyp (iptr) = optpar (4) (1:1)
          odbfmt (iptr) = deffmt (jdum)
          odbnam (iptr) = optpar (3)
          odbcom (iptr) = 'Extracted by field from '//optpar(5)
c
          if (odbtyp (iptr) .eq. 'R') then
            do i=1,idum
              rodb (i,j) = rdummy (i)
            end do
          else if (odbtyp (iptr) .eq. 'I') then
            do i=1,idum
              iodb (i,j) = idummy (i)
            end do
          end if
c
          call textut (' Extracted new odb :',odbnam(iptr))
c
        else if (optpar(2)(1:2) .eq. 'MU') then
c
          if (nopt .lt. 7) then
            optpar (7) = '(1x,5f8.3)'
            call textin (' Format to read from file ?',optpar(7))
          end if
c
          call exmult (iunit,optpar(5),rdummy,idummy,maxsiz,
     +                 optpar(4)(1:1),idum,i1,optpar(7),ierr)
          if (ierr .ne. 0) goto 10
          if (idum .le. 0) then
            call errcon ('No values read')
            goto 10
          end if
          if (idum .gt. maxsiz) then
            call errcon ('Too many values')
            goto 10
          end if
c
c ... copy to datablock etc.
c
          odbuse (iptr) = .true.
          odbptr (iptr) = j
          odbcha (iptr) = .true.
          odbsel (iptr) = .false.
          odblen (iptr) = idum
          odbtyp (iptr) = optpar (4) (1:1)
          odbfmt (iptr) = deffmt (jdum)
          odbnam (iptr) = optpar (3)
          odbcom (iptr) = 'Extracted by field from '//optpar(5)
c
          if (odbtyp (iptr) .eq. 'R') then
            do i=1,idum
              rodb (i,j) = rdummy (i)
            end do
          else if (odbtyp (iptr) .eq. 'I') then
            do i=1,idum
              iodb (i,j) = idummy (i)
            end do
          end if
c
          call textut (' Extracted new odb :',odbnam(iptr))
c
        else if (optpar(2)(1:2) .eq. 'PR') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'pro'
            call textin (' Molecule name (6 chars) ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
          optpar(3) = optpar(3)(1:6)
c
          if (nopt .lt. 4) then
            optpar (4) = 'procheck.out'
            call textin (' PROCHECK output file ?',optpar(4))
          end if
c
          call exproc (iunit,optpar(3),optpar(4),ierr)
          if (ierr .ne. 0) goto 10
c
        else
          call errcon ('Invalid EXtract command')
        end if
c
c ... SET
c
      else if (optpar(1)(1:2) .eq. 'SE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'ALl'
          call textin (' Option ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(3))
        end if
        i = whichm (optpar(3),ierr)
        idum = i
        if (ierr .ne. 0) goto 10
        prev = optpar (3)
        iptr = odbptr (i)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
c ... ALL
c
        if (optpar(2)(1:2) .eq. 'AL') then
c
          if (nopt .lt. 4) then
            optpar (4) = '0'
            call textin (' Value ?',optpar(4))
          end if
c
          call setodb (i,iptr,1,odblen(i),1,optpar(4),ierr)
          if (ierr .ne. 0) goto 10
c
          odbcha (idum) = .true.
c
c ... MANY
c
        else if (optpar(2)(1:2) .eq. 'MA') then
c
          if (nopt .lt. 4) then
            optpar (4) = '1'
            call textin (' First entry to set ?',optpar(4))
          end if
          call str2i (optpar(4),i1,ierr)
          if (ierr .ne. 0) goto 10
          i1 = max(1,i1)
c
          if (nopt .lt. 5) then
            write (optpar(5),*) odblen(i)
            call pretty (optpar(5))
            call textin (' Last entry to set ?',optpar(5))
          end if
          call str2i (optpar(5),i2,ierr)
          if (ierr .ne. 0) goto 10
          i2 = max(i1,min(i2,odblen(i)))
c
          if (nopt .lt. 6) then
            optpar (6) = '0'
            call textin (' Value ?',optpar(6))
          end if
c
          call setodb (i,iptr,i1,i2,1,optpar(6),ierr)
          if (ierr .ne. 0) goto 10
c
          odbcha (idum) = .true.
c
c ... ONE
c
        else if (optpar(2)(1:2) .eq. 'ON') then
c
          if (nopt .lt. 4) then
            optpar (4) = '1'
            call textin (' Entry to set ?',optpar(4))
          end if
          call str2i (optpar(4),i1,ierr)
          if (ierr .ne. 0) goto 10
          i1 = max(1,min(i1,odblen(i)))
c
          if (nopt .lt. 5) then
            optpar (5) = '0'
            call textin (' Value ?',optpar(5))
          end if
c
          call setodb (i,iptr,i1,i1,1,optpar(5),ierr)
          if (ierr .ne. 0) goto 10
c
          odbcha (idum) = .true.
c
c ... INDIVIDUAL
c
        else if (optpar(2)(1:2) .eq. 'IN') then
c
          if (nopt .lt. 4) then
            optpar (4) = '1'
            call textin (' First entry to set ?',optpar(4))
          end if
          call str2i (optpar(4),i1,ierr)
          if (ierr .ne. 0) goto 10
          i1 = max(1,i1)
c
          if (nopt .lt. 5) then
            write (optpar(5),*) odblen(i)
            call pretty (optpar(5))
            call textin (' Last entry to set ?',optpar(5))
          end if
          call str2i (optpar(5),i2,ierr)
          if (ierr .ne. 0) goto 10
          i2 = max(i1,min(i2,odblen(i)))
c
          call setind (i,iptr,i1,i2,ierr)
          if (ierr .ne. 0) goto 10
c
          odbcha (idum) = .true.
c
c ... IF
c
        else if (optpar(2)(1:2) .eq. 'IF') then
c
          if (nopt .lt. 4) then
            optpar (4) = '0'
            call textin (' Value ?',optpar(4))
          end if
c
          if (nopt .lt. 5) then
            optpar(5) = prev
            call textin (' Comparison ODB ?',optpar(5))
          end if
          j = whichm (optpar(5),ierr)
          if (ierr .ne. 0) goto 10
          jptr = odbptr (j)
c
          if (j .lt. 1) then
            call errcon ('ODB may not be a wildcard')
            goto 10
          end if
c
          if (odblen(j) .lt. odblen(i)) then
            call errcon ('Comparison ODB has too few elements')
            goto 10
          end if
c
          if (nopt .lt. 6) then
            optpar(6) = '='
            call textin (' Operator (=|<|>|#) ?',optpar(6))
          end if
          if (index('=<>#',optpar(6)(1:1)) .le. 0) then
            call errcon ('Invalid operator')
            goto 10
          end if
c
          if (nopt .lt. 7) then
            optpar (7) = '0'
            call textin (' Value ?',optpar(7))
          end if
c
          call setif (i,iptr,optpar(4),j,jptr,optpar(6),
     +                optpar(7),ierr)
          if (ierr .ne. 0) goto 10
c
          odbcha (idum) = .true.
c
c ... KEEP
c
        else if (optpar(2)(1:2) .eq. 'KE') then
c
          if (nopt .lt. 4) then
            optpar (4) = '1'
            call textin (' First entry to keep ?',optpar(4))
          end if
          call str2i (optpar(4),i1,ierr)
          if (ierr .ne. 0) goto 10
          i1 = max(1,i1)
c
          if (nopt .lt. 5) then
            write (optpar(5),*) odblen(i)
            call pretty (optpar(5))
            call textin (' Last entry to keep ?',optpar(5))
          end if
          call str2i (optpar(5),i2,ierr)
          if (ierr .ne. 0) goto 10
          i2 = max(i1,min(i2,odblen(i)))
c
          if (nopt .lt. 6) then
            optpar (6) = '1'
            call textin (' Step ?',optpar(6))
          end if
          call str2i (optpar(6),i3,ierr)
          if (ierr .ne. 0) goto 10
          i3 = max(1,i3)
c
          call setkep (i,iptr,i1,i2,i3,ierr)
          if (ierr .ne. 0) goto 10
c
          odbcha (idum) = .true.
c
c ... invalid SET command
c
        else
          call errcon ('Invalid SET command')
        end if
c
c ... READ
c
      else if (optpar(1)(1:2) .eq. 'RE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'junk.odb'
          call textin (' File name ?',optpar(2))
        end if
c
        call impodb (iunit,optpar(2),ierr)
c
        if (ierr .ne. 0) goto 10
c
        call prompt (' Datablock read okay')
c
c ... WRITE
c
      else if (optpar(1)(1:2) .eq. 'WR') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'junk.odb'
          call textin (' File name ?',optpar(3))
        end if
c
        call expodb (iunit,optpar(3),ierr)
c
        if (ierr .ne. 0) goto 10
c
        call prompt (' Datablock write okay')
c
c ... DUPLICATE
c
      else if (optpar(1)(1:2) .eq. 'DU') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
        idum = odbptr (i)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = 'junk'
          call textin (' Name for new odb ?',optpar(3))
        end if
        call upcase (optpar(3))
        optpar (3) = optpar(3)(1:25)
c
c ... allocate space (if any left)
c
        call allocm (optpar(3),odbtyp(i),odblen(i),odbfmt(i),
     +               j,iptr,ierr)
        if (ierr .ne. 0) goto 10
c
c ... okay
c
        odbuse (iptr) = .true.
        odbptr (iptr) = j
        odbcha (iptr) = .true.
        odbsel (iptr) = .false.
        odblen (iptr) = odblen(i)
        odbtyp (iptr) = odbtyp(i)
        odbfmt (iptr) = odbfmt(i)
        odbnam (iptr) = optpar(3)
        odbcom (iptr) = 'Duplicate of '//odbnam(i)
c
        if (odbtyp (iptr) .eq. 'R') then
          do i=1,odblen(i)
            rodb (i,j) = rodb(i,idum)
          end do
        else if (odbtyp (iptr) .eq. 'I') then
          do i=1,odblen(i)
            iodb (i,j) = iodb(i,idum)
          end do
        else if (odbtyp (iptr) .eq. 'C') then
          do i=1,odblen(i)
            codb (i,j) = codb(i,idum)
          end do
        else if (odbtyp (iptr) .eq. 'T') then
          do i=1,odblen(i)
            todb (i,j) = todb(i,idum)
          end do
        end if
c
        call textut (' Duplicate odb :',odbnam(iptr))
c
c ... CREATE
c
      else if (optpar(1)(1:2) .eq. 'CR') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'junk'
          call textin (' Name for new odb ?',optpar(2))
        end if
        call upcase (optpar(2))
        optpar (2) = optpar(2)(1:25)
c
        if (nopt .lt. 3) then
          optpar (3) = 'C'
          call textin (' Type (R/I/C/T) ?',optpar(3))
        end if
        call upcase (optpar(3))
        optpar (3) = optpar(3)(1:1)
        jdum = index(deftyp,optpar(3)(1:1))
        if (jdum .le. 0) then
          call errcon ('Invalid datablock type')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '1'
          call textin (' Number of elements ?',optpar(4))
        end if
        call str2i (optpar(4),idum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 5) then
          optpar (5) = deffmt(jdum)
          call textin (' Format ?',optpar(5))
        end if
        call upcase (optpar(5))
        call remspa (optpar(5))
c
        irf = max (index(optpar(5),'F'),
     +             index(optpar(5),'E'),
     +             index(optpar(5),'G'))
        iri = index(optpar(5),'I')
        irc = index(optpar(5),'A')
c
        if (optpar(3)(1:1) .eq. 'R' .and. irf .le. 0) then
          call errcon ('Unsuitable format for Real')
          optpar (5) = deffmt(jdum)
        else if (optpar(3)(1:1) .eq. 'I' .and. iri .le. 0) then
          call errcon ('Unsuitable format for Integer')
          optpar (5) = deffmt(jdum)
        else if (optpar(3)(1:1) .eq. 'C' .and. irc .le. 0) then
          call errcon ('Unsuitable format for Character')
          optpar (5) = deffmt(jdum)
        end if
c
c ... allocate space (if any left)
c
        call allocm (optpar(2),optpar(3),idum,optpar(5),
     +               j,iptr,ierr)
        if (ierr .ne. 0) goto 10
c
c ... okay
c
        odbuse (iptr) = .true.
        odbptr (iptr) = j
        odbcha (iptr) = .true.
        odbsel (iptr) = .false.
        odblen (iptr) = idum
        odbtyp (iptr) = optpar (3) (1:1)
        odbfmt (iptr) = optpar (5)
        odbnam (iptr) = optpar (2)
        odbcom (iptr) = 'Created from scratch'
c
        if (odbtyp (iptr) .eq. 'R') then
          do i=1,idum
            rodb (i,j) = 0.0
          end do
        else if (odbtyp (iptr) .eq. 'I') then
          do i=1,idum
            iodb (i,j) = 0
          end do
        else if (odbtyp (iptr) .eq. 'C') then
          do i=1,idum
            codb (i,j) = ' '
          end do
        else if (odbtyp (iptr) .eq. 'T') then
          do i=1,idum
            todb (i,j) = ' '
          end do
        end if
c
        call textut (' Initialised new odb :',odbnam(iptr))
c
c ... TYPE
c
      else if (optpar(1)(1:2) .eq. 'TY') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        call expodb (6,optpar(3),ierr)
c
c ... LIST
c
      else if (optpar(1)(1:2) .eq. 'LI') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        write (*,6100) 'Datablock name                ',
     +    'Type  Nr Alt? Format'
        do i=1,maxtdb
          if (odbsel(i) .and. odbuse(i)) then
            write (*,6110)
     +        odbnam(i),
     +        odbtyp(i),odblen(i),odbcha(i),
     +        odbfmt(i)(1:leng1(odbfmt(i)))
          end if
        end do
c
c ... STATS
c
      else if (optpar(1)(1:2) .eq. 'ST') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            write (*,*)
            call textut (' Stats :',odbnam(i))
            call xstats (rodb(1,iptr),odblen(i),
     +                   ave,sdv,xmin,xmax,xtot)
            call ivalut (' Nr of elements :',1,odblen(i))
            call rvalut (' Minimum  :',1,xmin)
            call rvalut (' Maximum  :',1,xmax)
            call rvalut (' Sum      :',1,xtot)
            call rvalut (' Average  :',1,ave)
            call rvalut (' St-devn  :',1,sdv)
            call rvalut (' Variance :',1,(sdv*sdv))
          end if
        end do
c
        do i=maxodb+1,2*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            write (*,*)
            call textut (' Stats :',odbnam(i))
            do j=1,odblen(i)
              rodb(j,0) = float(iodb(j,iptr))
            end do
            call xstats (rodb(1,0),odblen(i),
     +                   ave,sdv,xmin,xmax,xtot)
            call ivalut (' Nr of elements :',1,odblen(i))
            call rvalut (' Minimum  :',1,xmin)
            call rvalut (' Maximum  :',1,xmax)
            call rvalut (' Sum      :',1,xtot)
            call rvalut (' Average  :',1,ave)
            call rvalut (' St-devn  :',1,sdv)
            call rvalut (' Variance :',1,(sdv*sdv))
          end if
        end do
c
c ... FORMAT
c
      else if (optpar(1)(1:2) .eq. 'FO') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          optpar (3) = '???'
          do i=1,3*maxodb
            if (odbuse(i) .and. odbsel(i)) then
              optpar(3) = odbfmt(i)
              goto 6901
            end if
          end do
 6901     continue
          call textin (' New format ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) then
          call errcon ('No format provided')
          goto 10
        end if
c
        call upcase (optpar(3))
        call remspa (optpar(3))
        irf = max (index(optpar(3),'F'),
     +             index(optpar(3),'E'),
     +             index(optpar(3),'G'))
        iri = index(optpar(3),'I')
        irc = index(optpar(3),'A')
c
        if (irf .gt. 0 .and. iri .eq. 0 .and.
     +      irc .eq. 0) then
          call textut (' New Real format :',optpar(3))
          do i=1,maxodb
            if (odbuse(i) .and. odbsel(i)) then
              call textut (' Processing :',odbnam(i))
              odbfmt (i) = optpar(3)
              odbcha (i) = .true.
            end if
          end do
        else if (irf .eq. 0 .and. iri .gt. 0 .and.
     +           irc .eq. 0) then
          call textut (' New Integer format :',optpar(3))
          do i=maxodb+1,2*maxodb
            if (odbuse(i) .and. odbsel(i)) then
              call textut (' Processing :',odbnam(i))
              odbfmt (i) = optpar(3)
              odbcha (i) = .true.
            end if
          end do
        else if (irf .eq. 0 .and. iri .eq. 0 .and.
     +           irc .gt. 0) then
          call textut (' New Character format :',optpar(3))
          do i=2*maxodb+1,3*maxodb
            if (odbuse(i) .and. odbsel(i)) then
              call textut (' Processing :',odbnam(i))
              odbfmt (i) = optpar(3)
              odbcha (i) = .true.
            end if
          end do
        else
          call errcon ('Invalid format')
          call textut (' Format :',optpar(3))
        end if
c
c ... CGRAPH_FILE
c
      else if (optpar(1)(1:2) .eq. 'CG') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'odbman.cg'
          call textin (' CricketGraph file ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) then
          call errcon ('No filename provided')
          goto 10
        end if
c
        jmax = -1
        biglin = ' '
        k = 0
        lfirst = .false.
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            jmax = max (jmax, odblen(i))
            if (.not. lfirst) then
              biglin = odbnam(i)
              lfirst = .true.
            else
              call appstr (biglin,(tab//odbnam(i)))
            end if
            k = k + 1
          end if
          if (odbsel(i+maxodb) .and. odbuse(i+maxodb)) then
            jmax = max (jmax, odblen(i+maxodb))
            if (.not. lfirst) then
              biglin = odbnam(i+maxodb)
              lfirst = .true.
            else
              call appstr (biglin,(tab//odbnam(i+maxodb)))
            end if
            k = k + 1
          end if
        end do
c
        call jvalut (' Nr of datablocks selected :',1,k)
        call jvalut (' Max nr of elements        :',1,jmax)
        if (k .lt. 1) then
          call errcon ('No datablocks selected !')
          goto 10
        end if
c
        call xopxua (iunit,optpar(3),linter,ierr)
        if (ierr .ne. 0) goto 10
c
ccc        write (iunit,'(a1)') '*'
c
        call pretty (biglin)
        write (iunit,'(a)') biglin(1:leng1(biglin))
c
        do k=1,jmax
          lfirst = .false.
          do j=1,maxodb
            i = j
            if (odbsel(i) .and. odbuse(i)) then
              iptr = odbptr(i)
              if (odblen(i) .ge. k) then
                if (.not. lfirst) then
                  write (biglin,'(f15.5)') rodb(k,iptr)
                  lfirst = .true.
                else
                  write (line,'(a1,f15.5)') tab,rodb(k,iptr)
                  call appstr (biglin,line)
                end if
              else
                if (.not. lfirst) then
                  biglin = ' '
                  lfirst = .true.
                else
                  call appstr (biglin,tab)
                end if
              end if
            end if
c
            i = maxodb + j
            if (odbsel(i) .and. odbuse(i)) then
              iptr = odbptr(i)
              if (odblen(i) .ge. k) then
                if (.not. lfirst) then
                  write (biglin,'(f15.5)') float(iodb(k,iptr))
                  lfirst = .true.
                else
                  write (line,'(a1,f15.5)') tab,float(iodb(k,iptr))
                  call appstr (biglin,line)
                end if
              else
                if (.not. lfirst) then
                  biglin = ' '
                  lfirst = .true.
                else
                  call appstr (biglin,tab)
                end if
              end if
            end if
          end do
          call pretty (biglin)
          write (iunit,'(a)') biglin(1:leng1(biglin))
        end do
c
        close (iunit)
        call prompt (' CricketGraph file written')
c
c ... HISTO
c
      else if (optpar(1)(1:2) .eq. 'HI') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        numhis = 0
        do i=1,3
          if (nopt .lt. (i+2)) then
            write (optpar(i+2),*) float(i)
            call textin (' Histogram value ?',optpar(i+2))
          end if
          call str2r (optpar(i+2),xdum,ierr)
          if (ierr .ne. 0) goto 10
          numhis = numhis + 1
          rhis (numhis) = xdum
        end do
c
        if (nopt .gt. 5) then
          do i=6,nopt
            call str2r (optpar(i),xdum,ierr)
            if (ierr .ne. 0) goto 10
            numhis = numhis + 1
            rhis (numhis) = xdum
          end do
        end if
c
        call qsortg (rhis,numhis)
        call rvalut (' Histrogram limits :',numhis,rhis)
        write (*,*)
c
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Histogram :',odbnam(i))
            call histo (odblen(i),rodb(1,iptr),
     +                  numhis,rhis,nhis)
          end if
        end do
c
        do i=maxodb+1,2*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Histogram :',odbnam(i))
            do j=1,odblen(i)
              rodb(j,0) = float(iodb(j,iptr))
            end do
            call histo (odblen(i),rodb(1,0),
     +                  numhis,rhis,nhis)
          end if
        end do
c
c ... BOXCAR AVERAGING
c
      else if (optpar(1)(1:2) .eq. 'BO') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          write (optpar(3),*) '7'
          call textin (' Window size ?',optpar(3))
        end if
        call str2i (optpar(3),nwind,ierr)
        if (ierr .ne. 0) goto 10
        nw1 = (nwind-1)/2
        nw2 = nwind - nw1 - 1
        if (nw1 .le. 0 .or. nw2 .le. 0) then
          call errcon ('Window size too small')
          goto 10
        end if
c
ccc        print *,'NW ',nwind,nw1,nw2
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Box-car average :',odbnam(i))
            do j=1,odblen(i)
              k1 = max (1,j-nw1)
              k2 = min (odblen(i),j+nw2)
              xdum = 1.0 / float(k2-k1+1)
              buffer (j) = 0.0
ccc        print *,j,k1,k2
              do k=k1,k2
                buffer (j) = buffer (j) + rodb(k,iptr)
              end do
              buffer (j) = xdum * buffer (j)
            end do
            do j=1,odblen(i)
              rodb (j,iptr) = buffer (j)
            end do
            odbcha (i) = .true.
          end if
        end do
c
c ... 3D_ODL_FILE
c
      else if (optpar(1)(1:2) .eq. '3D') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which "X" odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(i) .ne. 'R' .and.
     +      odbtyp(i) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which "Y" odb (NO wildcards) ?',optpar(3))
        end if
        j = whichm (optpar(3),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(3)
c
        if (j .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(j) .ne. 'R' .and.
     +      odbtyp(j) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (odblen(i) .ne. odblen(j)) then
          call errcon ('Datablocks must have same size')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Which "Z" odb (NO wildcards) ?',optpar(4))
        end if
        k = whichm (optpar(4),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(4)
c
        if (k .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(k) .ne. 'R' .and.
     +      odbtyp(k) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (odblen(i) .ne. odblen(k)) then
          call errcon ('Datablocks must have same size')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = odbnam(i)(1:leng1(odbnam(i))) // 
     +      '_' // odbnam(j)(1:leng1(odbnam(j))) // 
     +      '_' // odbnam(k)(1:leng1(odbnam(k))) // 
     +      '.odl'
          call locase(optpar(5))
          call textin (' Plot file ?',optpar(5))
        end if
        if (length(optpar(5)) .lt. 1) then
          call errcon ('No filename provided')
          goto 10
        end if
c
        if (nopt .lt. 6) optpar (6) = 'cyan'
c
        iptr = odbptr (i)
        if (odbtyp(i) .eq. 'R') then
          do l=1,odblen(i)
            buffer (l) = rodb(l,iptr)
          end do
        else
          do l=1,odblen(i)
            buffer (l) = iodb(l,iptr)
          end do
        end if
c
        iptr = odbptr (j)
        if (odbtyp(j) .eq. 'R') then
          do l=1,odblen(j)
            rodb(l,0) = rodb(l,iptr)
          end do
        else
          do l=1,odblen(j)
            rodb(l,0) = iodb(l,iptr)
          end do
        end if
c
        iptr = odbptr (k)
        if (odbtyp(k) .eq. 'R') then
          do l=1,odblen(k)
            buff2(l) = rodb(l,iptr)
          end do
        else
          do l=1,odblen(k)
            buff2(l) = iodb(l,iptr)
          end do
        end if
c
        call xopxua (iunit,optpar(5),linter,ierr)
        if (ierr .ne. 0) goto 10
c
        call odl3d (iunit,odblen(i),buffer,rodb(1,0),buff2,
     +              optpar(6),ierr)
        if (ierr .ne. 0) goto 10
c
c ... COMBINE
c
      else if (optpar(1)(1:2) .eq. 'CO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(i) .ne. 'R' .and.
     +      odbtyp(i) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(3))
        end if
        j = whichm (optpar(3),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(3)
c
        if (j .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(j) .ne. 'R' .and.
     +      odbtyp(j) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (odblen(i) .ne. odblen(j)) then
          call errcon ('Datablocks must have same size')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '+'
          call textin (' Operator (+|-|*|/|min|max) ?',optpar(4))
        end if
        call remspa (optpar(4))
        call locase (optpar(4))
        if (length(optpar(4)) .lt. 1) then
          call errcon ('No operator provided')
          goto 10
        end if
c
        if (index('+|-|*|/',optpar(4)(1:1)) .gt. 0 .or.
     +      optpar(4)(1:3) .eq. 'min' .or.
     +      optpar(4)(1:3) .eq. 'max') goto 4253
c
        call errcon ('Not a valid operator !')
        goto 10
c
 4253   continue
c
        iptr = odbptr (i)
        if (odbtyp(i) .eq. 'R') then
          do k=1,odblen(i)
            buffer (k) = rodb(k,iptr)
          end do
        else
          do k=1,odblen(i)
            buffer (k) = iodb(k,iptr)
          end do
        end if
c
        iptr = odbptr (j)
        if (odbtyp(j) .eq. 'R') then
          do k=1,odblen(j)
            rodb(k,0) = rodb(k,iptr)
          end do
        else
          do k=1,odblen(j)
            rodb(k,0) = iodb(k,iptr)
          end do
        end if
c
        call textut (' Datablock 1 :',odbnam(i))
        call textut (' Datablock 2 :',odbnam(j))
c
        if (optpar(4)(1:1) .eq. '+') then
          call prompt (' DB1 = DB1 + DB2')
          do k=1,odblen(i)
            buffer (k) = buffer (k) + rodb (k,0)
          end do
        else if (optpar(4)(1:1) .eq. '-') then
          call prompt (' DB1 = DB1 - DB2')
          do k=1,odblen(i)
            buffer (k) = buffer (k) - rodb (k,0)
          end do
        else if (optpar(4)(1:1) .eq. '*') then
          call prompt (' DB1 = DB1 * DB2')
          do k=1,odblen(i)
            buffer (k) = buffer (k) * rodb (k,0)
          end do
        else if (optpar(4)(1:1) .eq. '/') then
          call prompt (' DB1 = DB1 / DB2')
          do k=1,odblen(i)
            if (rodb(k,0) .eq. 0.0) then
              call errcon (' Attempt to divide by zero !')
              call jvalut (' Zero element in DB2 is nr :',1,k)
              goto 10
            end if
            buffer (k) = buffer (k) / rodb (k,0)
          end do
        else if (optpar(4)(1:3) .eq. 'min') then
          call prompt (' DB1 = MIN (DB1, DB2)')
          do k=1,odblen(i)
            buffer (k) = min (buffer (k), rodb (k,0))
          end do
        else if (optpar(4)(1:3) .eq. 'max') then
          call prompt (' DB1 = MAX (DB1, DB2)')
          do k=1,odblen(i)
            buffer (k) = max (buffer (k), rodb (k,0))
          end do
        end if
c
        iptr = odbptr (i)
        if (odbtyp(i) .eq. 'R') then
ccc        print *,'RRR'
          do k=1,odblen(i)
            rodb(k,iptr) = buffer (k)
          end do
        else
ccc        print *,'III'
          do k=1,odblen(i)
            iodb(k,iptr) = nint (buffer (k))
          end do
        end if
c
        odbcha (i) = .true.
c
c ... SCATTER_PLOT
c
      else if (optpar(1)(1:2) .eq. 'SC') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(i) .ne. 'R' .and.
     +      odbtyp(i) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(3))
        end if
        j = whichm (optpar(3),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(3)
c
        if (j .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(j) .ne. 'R' .and.
     +      odbtyp(j) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (odblen(i) .ne. odblen(j)) then
          call errcon ('Datablocks must have same size')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = odbnam(j)(1:leng1(odbnam(j))) // 
     +      '_vs_' // odbnam(i)(1:leng1(odbnam(i))) // 
     +      '.plt'
          call locase(optpar(4))
          call textin (' Plot file ?',optpar(4))
        end if
        if (length(optpar(4)) .lt. 1) then
          call errcon ('No filename provided')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'Value of ' //
     +                 odbnam(i)(1:leng1(odbnam(i)))
        end if
        call textut (' X label :',optpar(5))
c
        if (nopt .lt. 6) then
          optpar (6) = 'Value of ' //
     +                 odbnam(j)(1:leng1(odbnam(j)))
        end if
        call textut (' Y label :',optpar(6))
c
        iptr = odbptr (i)
        if (odbtyp(i) .eq. 'R') then
          do k=1,odblen(i)
            buffer (k) = rodb(k,iptr)
          end do
          call xstats (buffer,odblen(i),
     +                 xave,xsdv,xmin,xmax,xtot)
        else
          do k=1,odblen(i)
            buffer (k) = iodb(k,iptr)
          end do
          call xstats (buffer,odblen(i),
     +                 xave,xsdv,xmin,xmax,xtot)
        end if
        rxmin = xmin
        rxmax = xmax
        xdum = 0.05 * (xmax - xmin)
        xmin = xmin - xdum
        xmax = xmax + xdum
        call rvalut (' X lower :',1,xmin)
        call rvalut (' X upper :',1,xmax)
c
        iptr = odbptr (j)
        if (odbtyp(j) .eq. 'R') then
          do k=1,odblen(j)
            rodb(k,0) = rodb(k,iptr)
          end do
          call xstats (rodb(1,0),odblen(j),
     +                 yave,ysdv,ymin,ymax,xtot)
        else
          do k=1,odblen(j)
            rodb(k,0) = iodb(k,iptr)
          end do
          call xstats (rodb(1,0),odblen(j),
     +                 yave,ysdv,ymin,ymax,xtot)
        end if
        rymin = ymin
        rymax = ymax
        ydum = 0.05 * (ymax - ymin)
        ymin = ymin - ydum
        ymax = ymax + ydum
        call rvalut (' Y lower :',1,ymin)
        call rvalut (' Y upper :',1,ymax)
c
        call xystat (buffer,rodb(1,0),min(odblen(i),odblen(j)),
     +               rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
c
        call xopxua (iunit,optpar(4),linter,ierr)
        if (ierr .ne. 0) goto 10
c
c
        call stamp (line)
        write (iunit,'(9a)',err=6908) 'REMARK ',
     +    line(1:leng1(line))
        write (iunit,'(9a)',err=6908) 'REMARK ',
     +    'Scatter plot of datablock (X) ',
     +     odbnam(i)(1:leng1(odbnam(i)))
        write (iunit,'(9a)',err=6908) 'REMARK ',
     +    'And datablock (Y) ',
     +     odbnam(j)(1:leng1(odbnam(j)))
c
        write (iunit,'(a,4(1x,1pe15.4))',err=6908)
     +    'REMARK X ave,sdv,min,max = ',
     +    xave,xsdv,rxmin,rxmax
c
        write (iunit,'(a,4(1x,1pe15.4))',err=6908)
     +    'REMARK Y ave,sdv,min,max = ',
     +    yave,ysdv,rymin,rymax
c
        write (iunit,'(a,1x,f8.3,1x,1pe15.4)',err=6908)
     +    'REMARK X-Y corr. coeff., rmsd = ',
     +    corr,rmsd
c
        write (iunit,'(a)',err=6908) 'LINFIT'
        write (iunit,'(a,i8)',err=6908) 'NPOINT ',odblen(i)
        write (iunit,'(9a)',err=6908) 'XLABEL ',
     +    optpar(5)(1:leng1(optpar(5)))
        write (iunit,'(9a)',err=6908) 'YLABEL ',
     +    optpar(6)(1:leng1(optpar(6)))
        write (iunit,'(a)',err=6908) 'COLOUR 4'
        write (iunit,'(a,4(1x,1p,e12.4))',err=6908)
     +    'XYVIEW ',xmin,xmax,ymin,ymax
        write (iunit,'(a,4(1x,1p,e12.4))',err=6908)
     +    'BOXPAR ',rxmin-0.5*xdum,xdum,rxmax+0.5*xdum
        write (iunit,'(a)',err=6908) 'XVALUE *'
        write (iunit,'(1p,6(1x,e12.4))',err=6908)
     +    (buffer(k),k=1,odblen(i))
        write (iunit,'(a)',err=6908) 'YVALUE *'
        write (iunit,'(1p,6(1x,e12.4))',err=6908)
     +    (rodb(k,0),k=1,odblen(i))
        write (iunit,'(a)',err=6908) 'END'
c
        call prompt (' O2D plot file written')
        goto 6906
c
 6908   continue
        call errcon ('While writing file')
c
 6906   continue
        close (iunit)
c
c ... PLOT_FILE
c
      else if (optpar(1)(1:2) .eq. 'PL') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(i) .ne. 'R' .and.
     +      odbtyp(i) .ne. 'I') then
          call errcon ('ODB must be R or I')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = odbnam(i)(1:leng1(odbnam(i))) // '.plt'
          call locase(optpar(3))
          call textin (' Plot file ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) then
          call errcon ('No filename provided')
          goto 10
        end if
c
        iptr = odbptr (i)
        if (odbtyp(i) .eq. 'R') then
          do j=1,odblen(i)
            buffer (j) = rodb(j,iptr)
          end do
          call xstats (buffer,odblen(i),
     +                 xave,xsdv,xmin,xmax,xtot)
        else
          do j=1,odblen(i)
            buffer (j) = iodb(j,iptr)
          end do
          call xstats (buffer,odblen(i),
     +                 xave,xsdv,xmin,xmax,xtot)
        end if
        rxmin = xmin
        rxmax = xmax
        xdum = 0.05 * (xmax - xmin)
c
        if (nopt .lt. 4) then
          xmin = xmin - xdum
          write (optpar(4),*) xmin
        end if
        call str2r (optpar(4),xmin,ierr)
        if (ierr .ne. 0) goto 10
        call rvalut (' Y lower :',1,xmin)
c
        if (nopt .lt. 5) then
          xmax = xmax + xdum
          write (optpar(5),*) xmax
        end if
        call str2r (optpar(5),xmax,ierr)
        if (ierr .ne. 0) goto 10
        call rvalut (' Y upper :',1,xmax)
c
        if (nopt .lt. 6) then
          optpar (6) = 'Index in datablock'
        end if
        call textut (' X label :',optpar(6))
c
        if (nopt .lt. 7) then
          optpar (7) = 'Value of ' //
     +                 odbnam(i)(1:leng1(odbnam(i)))
        end if
        call textut (' Y label :',optpar(7))
c
        call xopxua (iunit,optpar(3),linter,ierr)
        if (ierr .ne. 0) goto 10
c
        call stamp (line)
        write (iunit,'(9a)',err=6902) 'REMARK ',
     +    line(1:leng1(line))
        write (iunit,'(9a)',err=6908) 'REMARK ',
     +    'Plot of datablock ',
     +     odbnam(i)(1:leng1(odbnam(i)))
        write (iunit,'(a,4(1x,1pe15.4))',err=6908)
     +    'REMARK Datablock ave,sdv,min,max = ',
     +    xave,xsdv,rxmin,rxmax
        write (iunit,'(a)',err=6902) 'LINFIT'
        write (iunit,'(a,i8)',err=6902) 'NPOINT ',odblen(i)
        write (iunit,'(9a)',err=6902) 'XLABEL ',
     +    optpar(6)(1:leng1(optpar(6)))
        write (iunit,'(9a)',err=6902) 'YLABEL ',
     +    optpar(7)(1:leng1(optpar(7)))
        write (iunit,'(a)',err=6902) 'COLOUR 4'
        write (iunit,'(a)',err=6902) 'XLIMIT 1.0 1.0'
        write (iunit,'(a,4(1x,1p,e12.4))',err=6902)
     +    'XYVIEW ',0.0,float(odblen(i)+1),xmin,xmax
        write (iunit,'(a)',err=6902) 'YVALUE *'
c
        write (iunit,'(1p,6(1x,e12.4))',err=6902)
     +    (buffer(j),j=1,odblen(i))
        write (iunit,'(a)',err=6902) 'END'
c
        call prompt (' O2D plot file written')
        goto 6904
c
 6902   continue
        call errcon ('While writing file')
c
 6904   continue
        close (iunit)
c
c ... NAME
c
      else if (optpar(1)(1:2) .eq. 'NA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('ODB may not be a wildcard')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = odbnam(i)
          call textin (' New name ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) then
          call errcon ('No name provided')
          goto 10
        end if
c
        call textut (' Old name :',odbnam (i))
        odbnam (i) = optpar(3)
        call upcase (odbnam(i))
        call remspa (odbnam(i))
        call textut (' New name :',odbnam (i))
        odbcha (i) = .true.
c
c ... SIMILARITY
c
      else if (optpar(1)(1:2) .eq. 'SI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('First ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(i) .ne. 'R' .and.
     +      odbtyp(i) .ne. 'I') then
          call errcon ('ODB must be numerical')
          goto 10
        end if
c
        iptr = odbptr(i)
        if (odbtyp(i) .eq. 'R') then
          do j=1,odblen(i)
            buffer (j) = rodb(j,iptr)
          end do
        else
          do j=1,odblen(i)
            buffer (j) = float (iodb(j,iptr))
          end do
        end if
c
        idum = odblen(i)
        jdum = i
c
        if (nopt .lt. 3) then
          optpar (3) = '*'
          call textin (' Which odb ?',optpar(3))
        end if
        i = whichm (optpar(3),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(3)
c
        cormin = -1.1
        if (nopt .ge. 4) then
          call str2r (optpar(4),cormin,ierr)
          if (ierr .ne. 0 .or. cormin .gt. 1.0) then
            cormin = -1.1
          end if
          call fvalut (' Minimum correlation coefficient :',1,cormin)
        end if
c
        cormax = 1.1
        if (nopt .ge. 5) then
          call str2r (optpar(5),cormax,ierr)
          if (ierr .ne. 0 .or. cormax .lt. cormin) then
            cormax = min (cormin,1.1)
          end if
          call fvalut (' Maximum correlation coefficient :',1,cormax)
        end if
c
        if (nopt .ge. 6) then
          call upcase (optpar(6))
          call remspa (optpar(6))
        else
          optpar(6) = 'Table'
        end if
        ltable = .true.
        if (optpar(6)(1:1) .eq. 'L') ltable = .false.
c
        nhit = 0
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            if (idum .eq. odblen(i) .and. i.ne.jdum) then
              iptr = odbptr (i)
              call xystat (buffer,rodb(1,iptr),idum,
     +                     rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
              if (corr .ge. cormin .and. corr .le. cormax) then
                nhit = nhit + 1
                if (ltable) then
                  if (nhit .eq. 1) write (*,6740) 'Hit #','#Elem',
     +              'Corr-C','RMS Diff',
     +              '--- Datablock 1 ---','--- Datablock 2 ---'
                  write (*,6741) nhit,idum,corr,rmsd,odbnam(jdum),
     +                           odbnam(i)
                else
                  write (*,*)
                  call textut (' Comparing :',odbnam(jdum))
                  call textut (' And       :',odbnam(i))
                  call ivalut (' Nr of elements   :',1,odblen(i))
                  call fvalut (' Correl Coeff     :',1,corr)
                  call rvalut (' RMS difference:',1,rmsd)
                  call fvalut (' Shape Similarity :',1,shap)
                  call rvalut (' R-factor (1)  :',1,rf1)
                  call fvalut (' Ditto, 2 scaled  :',1,x6)
                  call fvalut (' Scale for 2      :',1,x8)
                  call rvalut (' R-factor (2)  :',1,rf2)
                  call fvalut (' Ditto, 1 scaled  :',1,x7)
                  call fvalut (' Scale for 1      :',1,x9)
                end if
              end if
            end if
          end if
        end do
c
        do i=maxodb+1,2*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            if (idum .eq. odblen(i) .and. i.ne.jdum) then
              iptr = odbptr (i)
              do j=1,odblen(i)
                rodb(j,0) = float(iodb(j,iptr))
              end do
              call xystat (buffer,rodb(1,0),idum,
     +                     rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
              if (corr .ge. cormin .and. corr .le. cormax) then
                nhit = nhit + 1
                if (ltable) then
                  if (nhit .eq. 1) write (*,6740) 'Hit #','#Elem',
     +              'Corr-C','RMS Diff',
     +              '--- Datablock 1 ---','--- Datablock 2 ---'
                  write (*,6741) nhit,idum,corr,rmsd,odbnam(jdum),
     +                           odbnam(i)
                else
                  write (*,*)
                  call textut (' Comparing :',odbnam(jdum))
                  call textut (' And       :',odbnam(i))
                  call ivalut (' Nr of elements   :',1,odblen(i))
                  call fvalut (' Correl Coeff     :',1,corr)
                  call rvalut (' RMS difference   :',1,rmsd)
                  call fvalut (' Shape Similarity :',1,shap)
                  call rvalut (' R-factor (1)     :',1,rf1)
                  call fvalut (' Ditto, 2 scaled  :',1,x6)
                  call fvalut (' Scale for 2      :',1,x8)
                  call rvalut (' R-factor (2)     :',1,rf2)
                  call fvalut (' Ditto, 1 scaled  :',1,x7)
                  call fvalut (' Scale for 1      :',1,x9)
                end if
              end if
            end if
          end if
        end do
c
        call jvalut (' Nr of hits :',1,nhit)
c
c ... NON-PARAMETRIC
c
      else if (optpar(1)(1:2) .eq. 'NO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which odb (NO wildcards) ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (i .lt. 1) then
          call errcon ('First ODB may not be a wildcard')
          goto 10
        end if
c
        if (odbtyp(i) .ne. 'R' .and.
     +      odbtyp(i) .ne. 'I') then
          call errcon ('ODB must be numerical')
          goto 10
        end if
c
        iptr = odbptr(i)
        if (odbtyp(i) .eq. 'R') then
          do j=1,odblen(i)
            buffer (j) = rodb(j,iptr)
          end do
        else
          do j=1,odblen(i)
            buffer (j) = float (iodb(j,iptr))
          end do
        end if
c
        idum = odblen(i)
        jdum = i
c
        if (nopt .lt. 3) then
          optpar (3) = '*'
          call textin (' Which odb ?',optpar(3))
        end if
        i = whichm (optpar(3),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(3)
c
        cormin = -1.1
        if (nopt .ge. 4) then
          call str2r (optpar(4),cormin,ierr)
          if (ierr .ne. 0 .or. cormin .gt. 1.0) then
            cormin = -1.1
          end if
          call fvalut (' Minimum correlation coefficient :',1,cormin)
        end if
c
        cormax = 1.1
        if (nopt .ge. 5) then
          call str2r (optpar(5),cormax,ierr)
          if (ierr .ne. 0 .or. cormax .lt. cormin) then
            cormax = min (cormin,1.1)
          end if
          call fvalut (' Maximum correlation coefficient :',1,cormax)
        end if
c
        if (nopt .ge. 6) then
          call upcase (optpar(6))
          call remspa (optpar(6))
        else
          optpar(6) = 'Table'
        end if
        ltable = .true.
        if (optpar(6)(1:1) .eq. 'L') ltable = .false.
c
        nhit = 0
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            if (idum .eq. odblen(i) .and. i.ne.jdum) then
              iptr = odbptr (i)
              call nonpar (buffer,rodb(1,iptr),idum,
     +                     buff1,buff2,d,zd,probd,rs,probrs)
              if (rs .ge. cormin .and. rs .le. cormax) then
                nhit = nhit + 1
                if (ltable) then
                  if (nhit .eq. 1) write (*,6750) 'Hit #','#Elem',
     +              'Rank-CC','Prob','Z(D)','Prob',
     +              '--- Datablock 1 ---','--- Datablock 2 ---'
                  write (*,6751) nhit,idum,rs,probrs,zd,probd,
     +                           odbnam(jdum),odbnam(i)
                else
                  write (*,*)
                  call textut (' Comparing :',odbnam(jdum))
                  call textut (' And       :',odbnam(i))
                  call ivalut (' Nr of elements    :',1,odblen(i))
                  call fvalut (' Rank corr. coeff. :',1,rs)
                  call rvalut (' Probability    :',1,probrs)
                  call rvalut (' D-value        :',1,d)
                  call fvalut (' Z-score (D)       :',1,zd)
                  call rvalut (' Probability    :',1,probd)
                end if
              end if
            end if
          end if
        end do
c
        do i=maxodb+1,2*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            if (idum .eq. odblen(i) .and. i.ne.jdum) then
              iptr = odbptr (i)
              do j=1,odblen(i)
                rodb(j,0) = float(iodb(j,iptr))
              end do
              call nonpar (buffer,rodb(1,0),idum,
     +                     buff1,buff2,d,zd,probd,rs,probrs)
              if (rs .ge. cormin .and. rs .le. cormax) then
                nhit = nhit + 1
                if (ltable) then
                  if (nhit .eq. 1) write (*,6750) 'Hit #','#Elem',
     +              'Rank-CC','Prob','Z(D)','Prob',
     +              '--- Datablock 1 ---','--- Datablock 2 ---'
                  write (*,6751) nhit,idum,rs,probrs,zd,probd,
     +                           odbnam(jdum),odbnam(i)
                else
                  write (*,*)
                  call textut (' Comparing :',odbnam(jdum))
                  call textut (' And       :',odbnam(i))
                  call ivalut (' Nr of elements    :',1,odblen(i))
                  call fvalut (' Rank corr. coeff. :',1,rs)
                  call rvalut (' Probability    :',1,probrs)
                  call rvalut (' D-value        :',1,d)
                  call fvalut (' Z-score (D)       :',1,zd)
                  call rvalut (' Probability    :',1,probd)
                end if
              end if
            end if
          end if
        end do
c
        call jvalut (' Nr of hits :',1,nhit)
c
c ... ALL_CORRELATIONS
c
      else if (optpar(1)(1:2) .eq. 'AL') then
c
        cormin = -1.1
        if (nopt .ge. 2) then
          call str2r (optpar(2),cormin,ierr)
          if (ierr .ne. 0 .or. cormin .gt. 1.0) then
            cormin = -1.1
          end if
          call fvalut (' Minimum correlation coefficient :',1,cormin)
        end if
c
        cormax = 1.1
        if (nopt .ge. 3) then
          call str2r (optpar(3),cormax,ierr)
          if (ierr .ne. 0 .or. cormax .lt. cormin) then
            cormax = min (cormin,1.1)
          end if
          call fvalut (' Maximum correlation coefficient :',1,cormax)
        end if
c
        if (nopt .ge. 4) then
          call upcase (optpar(4))
          call remspa (optpar(4))
        else
          optpar(4) = 'Table'
        end if
        ltable = .true.
        if (optpar(4)(1:1) .eq. 'L') ltable = .false.
c
        nhit = 0
c
        do i=1,2*maxodb-1
c
          if (.not. odbuse(i)) goto 6206
c
          iptr = odbptr(i)
          if (odbtyp(i) .eq. 'R') then
            do j=1,odblen(i)
              buffer (j) = rodb(j,iptr)
            end do
          else
            do j=1,odblen(i)
              buffer (j) = float (iodb(j,iptr))
            end do
          end if
c
          idum = odblen(i)
          jdum = i
c
          do j=i+1,2*maxodb
c
            if (.not. odbuse(j)) goto 6205
            if (odblen(j) .ne. idum) goto 6205
c
            jptr = odbptr (j)
            if (odbtyp(j) .eq. 'R') then
              call xystat (buffer,rodb(1,jptr),idum,
     +                     rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
            else
              do k=1,odblen(i)
                rodb(k,0) = float(iodb(k,jptr))
              end do
              call xystat (buffer,rodb(1,0),idum,
     +                     rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
            end if
            if (corr .ge. cormin .and. corr .le. cormax) then
              nhit = nhit + 1
              if (ltable) then
                if (nhit .eq. 1) write (*,6740) 'Hit #','#Elem',
     +            'Corr-C','RMS Diff',
     +            '--- Datablock 1 ---','--- Datablock 2 ---'
                write (*,6741) nhit,idum,corr,rmsd,odbnam(jdum),
     +                         odbnam(j)
              else
                write (*,*)
                call textut (' Comparing :',odbnam(jdum))
                call textut (' And       :',odbnam(j))
                call ivalut (' Nr of elements   :',1,odblen(i))
                call fvalut (' Correl Coeff     :',1,corr)
                call rvalut (' RMS difference   :',1,rmsd)
                call fvalut (' Shape Similarity :',1,shap)
              end if
            end if
 6205       continue
          end do
 6206     continue
        end do
c
        call jvalut (' Nr of hits :',1,nhit)
c
 6740 format (/2(1x,a6),1x,a6,1x,a12,2(1x,a25))
 6741 format (2(1x,i6),1x,f6.3,1x,1pe12.4,2(1x,a25))
c
 6750 format (/2(1x,a6),2(1x,a6,1x,a8),2(1x,a25))
 6751 format (2(1x,i6),(1x,f6.3,1x,1pe8.1,0p),
     +        (1x,f8.3,1x,1pe8.1,0p),2(1x,a25))
c
c ... FLOAT_FUNCTION
c
      else if (optpar(1)(1:2) .eq. 'FL') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          optpar (3) = '?'
          call textin (' Function ?',optpar(3))
        end if
        call remspa (optpar(3))
        call upcase (optpar(3))
        if (optpar(3)(1:1) .eq. '?') then
          write (*,'(99(1x,a/))')
     +      'Available floating point functions :',
     +      'LOwer_bound_reset',
     +      'UPper_bound_reset',
     +      'ADd',
     +      'MUltiply',
     +      'INteger_cutoff',
     +      'NEarest_integer',
     +      'SQuare_root',
     +      'POwer_raise'
          goto 10
        end if
c
        if (optpar(3)(1:2) .ne. 'UP' .and.
     +      optpar(3)(1:2) .ne. 'LO' .and.
     +      optpar(3)(1:2) .ne. 'AD' .and.
     +      optpar(3)(1:2) .ne. 'MU' .and.
     +      optpar(3)(1:2) .ne. 'IN' .and.
     +      optpar(3)(1:2) .ne. 'NE' .and.
     +      optpar(3)(1:2) .ne. 'SQ' .and.
     +      optpar(3)(1:2) .ne. 'PO') then
          call errcon ('Invalid floating point function')
          goto 10
        end if
        call textut (' Floating point function :',optpar(3))
c
        if (optpar(3)(1:2) .eq. 'UP' .or.
     +      optpar(3)(1:2) .eq. 'LO' .or.
     +      optpar(3)(1:2) .eq. 'AD' .or.
     +      optpar(3)(1:2) .eq. 'MU' .or.
     +      optpar(3)(1:2) .eq. 'PO') then
          if (nopt .lt. 4) then
            optpar(4) = '1.0'
            call textin (' Value ?',optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .ne. 0) goto 10
          call rvalut (' Value :',1,xdum)
        end if
c
        do i=1,maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Processing :',odbnam(i))
            if (optpar(3)(1:2) .eq. 'UP') then
              do j=1,odblen(i)
                rodb(j,iptr) = min (xdum,rodb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'LO') then
              do j=1,odblen(i)
                rodb(j,iptr) = max (xdum,rodb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'AD') then
              do j=1,odblen(i)
                rodb(j,iptr) = xdum + rodb(j,iptr)
              end do
            else if (optpar(3)(1:2) .eq. 'MU') then
              do j=1,odblen(i)
                rodb(j,iptr) = xdum * rodb(j,iptr)
              end do
            else if (optpar(3)(1:2) .eq. 'PO') then
              do j=1,odblen(i)
                rodb(j,iptr) = rodb(j,iptr) ** xdum
              end do
            else if (optpar(3)(1:2) .eq. 'IN') then
              do j=1,odblen(i)
                rodb(j,iptr) = int (rodb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'NE') then
              do j=1,odblen(i)
                rodb(j,iptr) = nint (rodb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'SQ') then
              do j=1,odblen(i)
                if (rodb(j,iptr) .ge. 0.0) then
                  rodb(j,iptr) = sqrt (rodb(j,iptr))
                else
                  call rvalut (
     +              ' ERROR --- Cannot take SQRT :',1,
     +              rodb(j,iptr))
                end if
              end do
            end if
            odbcha (i) = .true.
          end if
        end do
c
c ... INTEGER_FUNCTION
c
      else if (optpar(1)(1:2) .eq. 'IN') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          optpar (3) = '?'
          call textin (' Function ?',optpar(3))
        end if
        call remspa (optpar(3))
        call upcase (optpar(3))
        if (optpar(3)(1:1) .eq. '?') then
          write (*,'(99(1x,a/))')
     +      'Available integer functions :',
     +      'LOwer_bound_reset',
     +      'UPper_bound_reset',
     +      'ADd',
     +      'MUltiply'
          goto 10
        end if
c
        if (optpar(3)(1:2) .ne. 'UP' .and.
     +      optpar(3)(1:2) .ne. 'LO' .and.
     +      optpar(3)(1:2) .ne. 'AD' .and.
     +      optpar(3)(1:2) .ne. 'MU') then
          call errcon ('Invalid integer function')
          goto 10
        end if
        call textut (' Integer function :',optpar(3))
c
        if (nopt .lt. 4) then
          optpar(4) = '0'
          call textin (' Value ?',optpar(4))
        end if
        call str2i (optpar(4),idum,ierr)
        if (ierr .ne. 0) goto 10
        call ivalut (' Value :',1,idum)
c
        do i=maxodb+1,2*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Processing :',odbnam(i))
            if (optpar(3)(1:2) .eq. 'UP') then
              do j=1,odblen(i)
                iodb(j,iptr) = min (idum,iodb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'LO') then
              do j=1,odblen(i)
                iodb(j,iptr) = max (idum,iodb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'AD') then
              do j=1,odblen(i)
                iodb(j,iptr) = idum + iodb(j,iptr)
              end do
            else if (optpar(3)(1:2) .eq. 'MU') then
              do j=1,odblen(i)
                iodb(j,iptr) = idum * iodb(j,iptr)
              end do
            end if
            odbcha (i) = .true.
          end if
        end do
c
c ... CHARACTER_FUNCTION
c
      else if (optpar(1)(1:2) .eq. 'CH') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
        if (nopt .lt. 3) then
          optpar (3) = '?'
          call textin (' Function ?',optpar(3))
        end if
        call remspa (optpar(3))
        call upcase (optpar(3))
        if (optpar(3)(1:1) .eq. '?') then
          write (*,'(99(1x,a/))')
     +      'Available character functions :',
     +      'UPpercase',
     +      'LOwercase',
     +      'PRetty',
     +      'REmove_spaces'
          goto 10
        end if
c
        if (optpar(3)(1:2) .ne. 'UP' .and.
     +      optpar(3)(1:2) .ne. 'LO' .and.
     +      optpar(3)(1:2) .ne. 'PR' .and.
     +      optpar(3)(1:2) .ne. 'RE') then
          call errcon ('Invalid character function')
          goto 10
        end if
c
        call textut (' Character function :',optpar(3))
c
        do i=2*maxodb+1,3*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Processing :',odbnam(i))
            if (optpar(3)(1:2) .eq. 'UP') then
              do j=1,odblen(i)
                call upcase (codb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'LO') then
              do j=1,odblen(i)
                call locase (codb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'PR') then
              do j=1,odblen(i)
                call pretty (codb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'RE') then
              do j=1,odblen(i)
                call remspa (codb(j,iptr))
              end do
            end if
            odbcha (i) = .true.
          end if
        end do
c
        do i=3*maxodb+1,4*maxodb
          if (odbsel(i) .and. odbuse(i)) then
            iptr = odbptr(i)
            call textut (' Processing :',odbnam(i))
            if (optpar(3)(1:2) .eq. 'UP') then
              do j=1,odblen(i)
                call upcase (todb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'LO') then
              do j=1,odblen(i)
                call locase (todb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'PR') then
              do j=1,odblen(i)
                call pretty (todb(j,iptr))
              end do
            else if (optpar(3)(1:2) .eq. 'RE') then
              do j=1,odblen(i)
                call remspa (todb(j,iptr))
              end do
            end if
            odbcha (i) = .true.
          end if
        end do
c
c ... DELETE
c
      else if (optpar(1)(1:2) .eq. 'DE') then
c
        if (nopt .lt. 2) then
          optpar (2) = '*'
          call textin (' Which odb ?',optpar(2))
        end if
        i = whichm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar(2)
c
c ... if some wild-card used, ask confirmation if interactive
c
        if (i .lt. 1 .and. linter) then
          reply = 'N'
          call textut (' Delete :',optpar(2))
          call textin (' Sure (Y/N) ?',reply)
          call upcase (reply)
          if (reply .ne. 'Y') goto 10
        end if
c
        do i=1,maxtdb
          if (odbsel(i) .and. odbuse(i)) then
            call textut (' Deleting :',odbnam(i))
            odbuse (i) = .false.
            odbsel (i) = .false.
            odbnam (i) = '!@#$%^&*?'
          end if
        end do
c
 6100 format (/1x,a25,a)
 6110 format (1x,a25,1x,a1,1x,i6,1x,l1,1x,a)
c
c ... INVALID OPTION
c
      else
c
        call errcon ('Invalid option')
        call textut (' ==>',line)
c
      end if
c
      goto 10
c
c --- END OF MAIN EVENT LOOP
c
c ... check if all changes saved
c
 9000 continue
c
      write (*,*)
      unsave = .false.
      do i=1,maxtdb
        if (odbuse(i)) then
          if (odbcha(i)) then
              call textut (' There are unsaved changes to :',
     +          odbnam(i))
            unsave = .true.
          end if
        end if
      end do
c
c ... if unsaved changes and interactive, ask if user is sure
c
      if (unsave .and. linter) then
        reply = 'N'
        call textin (
     +    ' Do you really want to quit (Y/N) ?',reply)
        call upcase (reply)
        if (reply .ne. 'Y') goto 10
      end if
c
      call gkquit
c
      end
c
c
c
      subroutine allocm (name,type,ndata,fmt,
     +                   j,iptr,ierr)
c
c ... allocate a datablock
c
c .... J = returned = index in arrays of this type of datablock
c      IPTR =  ,,   = index in arrays of ALL datablocks
c
      include 'odbman.incl'
c
      integer length,j,iptr,i,l1,l2,ierr,ndata
c
      character name*(*),type*1,fmt*(*)
c
code ...
c
      ierr = -1
      j = -1
      iptr = -1
c
c ... check if datablock not too big
c
      if (ndata .gt. maxsiz) then
        call errcon ('Datablock too big')
        call jvalut (' Requested size :',1,ndata)
        call jvalut (' Maximum size   :',1,maxsiz)
        return
      end if
c
c ... check if name is unique
c
      l1 = length(name)
      do i=1,maxtdb
        if (odbuse(i)) then
          l2 = length(odbnam(i))
          if (l1 .eq. l2) then
            if (name(1:l1) .eq. odbnam(i)(1:l1)) then
              call errcon ('Datablock has non-unique name')
              call textut (' Name :',name)
              return
            end if
          end if
        end if
      end do
c
c ... allocate it if there's space left
c
      if (type.eq.'R') then
        do i=1,maxodb
          if (.not. odbuse(i)) then
            j = i
            iptr = i
            ierr = 0
            return
          end if
        end do
      else if (type.eq.'I') then
        do i=maxodb+1,2*maxodb
          if (.not. odbuse(i)) then
            j = i - maxodb
            iptr = i
            ierr = 0
            return
          end if
        end do
      else if (type.eq.'C') then
        do i=2*maxodb+1,3*maxodb
          if (.not. odbuse(i)) then
            j = i - 2*maxodb
            iptr = i
            ierr = 0
            return
          end if
        end do
      else if (type.eq.'T') then
        do i=3*maxodb+1,4*maxodb
          if (.not. odbuse(i)) then
            j = i - 3*maxodb
            iptr = i
            ierr = 0
            return
          end if
        end do
      else
        call errcon ('Invalid datablock type')
        call textut (' Type (should be I,R,C,T) :',type)
        return
      end if
c
      call errcon ('No more space for datablock')
      call textut (' Type :',type)
c
      return
      end
c
c
c
      integer function whichm (nam,ierr)
c
c ... which datablock does the name "nam" correspond to ?
c
c ... the following are special names:
c
c     * = all datablocks
c     # = all INTEGER datablocks
c     % = all REAL datablocks
c     #% = %# = all NUMERICAL datablocks
c     $ = all CHARACTER*6 datablocks
c     @ = all TEXT datablocks
c     $@ = @$ = all CHARACTER datablocks
c     *XYZ*, *XYZ, XYZ* = any matching datablocks
c     ? = list all these (returns error)
c
      include 'odbman.incl'
c
      integer i,ll,l2,l3,length,ierr
c
      character nam*(*)
c
code ...
c
      ierr = -1
      whichm = -999
      ll = length(nam)
      if (ll .le. 0) return
c
      ierr = 0
      do i=1,maxtdb
        odbsel (i) = .false.
      end do
c
c ... is it a special symbol ?
c
      if (ll .eq. 1) then
c
        if (nam(1:1) .eq. '?') then
          write (*,'(99(1x,a/))')
     +      'ODB selector may be:',
     +      'datablock name',
     +      '* for all datablocks',
     +      '# for all integer datablocks',
     +      '% for all real datablocks',
     +      '#% or %# for all numerical datablocks',
     +      '$ for all character datablocks',
     +      '@ for all text datablocks',
     +      '$@ or @$ for all character and text datablocks',
     +      '*XYZ* for all datablocks containing "XYZ"',
     +      '*XYZ  for all datablocks ending in "XYZ"',
     +      'XYZ*  for all datablocks beginning with "XYZ"',
     +      '? for this list'
          ierr = -2
          return
c
        else if (nam(1:1) .eq. '*') then
          do i=1,maxtdb
            odbsel (i) = odbuse (i)
          end do
          whichm = all
          return
        else if (nam(1:1) .eq. '%') then
          do i=1,maxodb
            odbsel (i) = odbuse (i)
          end do
          whichm = allr
          return
        else if (nam(1:1) .eq. '#') then
          do i=maxodb+1,2*maxodb
            odbsel (i) = odbuse (i)
          end do
          whichm = alli
          return
        else if (nam(1:1) .eq. '$') then
          do i=2*maxodb+1,3*maxodb
            odbsel (i) = odbuse (i)
          end do
          whichm = allc
          return
        else if (nam(1:1) .eq. '@') then
          do i=3*maxodb+1,4*maxodb
            odbsel (i) = odbuse (i)
          end do
          whichm = allt
          return
        end if
      end if
c
      if (ll .eq. 2) then
        if (nam(1:2) .eq. '#%' .or.
     +      nam(1:2) .eq. '%#') then
          do i=1,2*maxodb
            odbsel (i) = odbuse (i)
          end do
          whichm = alln
          return
        else if (nam(1:2) .eq. '$@' .or.
     +          nam(1:2) .eq. '@$') then
          do i=2*maxodb+1,4*maxodb
            odbsel (i) = odbuse (i)
          end do
          whichm = alls
          return
        end if
      end if
c
      call upcase (nam)
      whichm = -998
c
c ... *string*
c
      if (nam(1:1) .eq. '*' .and. nam(ll:ll) .eq. '*') then
        do i=1,maxtdb
          if (odbuse(i)) then
            if (index(odbnam(i),nam(2:ll-1)) .gt. 0) then
              odbsel (i) = .true.
            end if
          end if
        end do
        whichm = many
        return
c
c ... *string
c
      else if (nam(1:1) .eq. '*') then
        do i=1,maxtdb
          if (odbuse(i)) then
            l2 = length(odbnam(i))
            l3 = l2 - ll + 2
            if (odbnam(i)(l3:l2) .eq. nam(2:ll)) then
              odbsel (i) = .true.
            end if
          end if
        end do
        whichm = many
        return
c
c ... string*
c
      else if (nam(ll:ll) .eq. '*') then
        do i=1,maxtdb
          if (odbuse(i)) then
            l3 = ll - 1
            if (odbnam(i)(1:l3) .eq. nam(1:l3)) then
              odbsel (i) = .true.
            end if
          end if
        end do
        whichm = many
        return
      end if
c
c ... it must be a full datablock name, then
c
      do i=1,maxtdb
        if (odbuse(i)) then
          l2 = length(odbnam(i))
          if (ll .eq. l2) then
            if (nam(1:ll) .eq. odbnam(i)(1:ll)) then
              odbsel (i) = .true.
              whichm = i
              return
            end if
          end if
        end if
      end do
c
      ierr = -1
      whichm = -997
c
      call errcon ('Unrecognised datablock name')
      call textut (' Name :',nam)
c
      return
      end
c
c
c
      subroutine impodb (iunit,file,ierr)
c
      include 'odbman.incl'
c
      integer iunit,ierr,nlin,nodb,length,nopt,ndata,j,iptr,i
      integer lmax,leng1
c
      character line*256,pars(4)*80,file*(*)
c
code ...
c
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening ODB file')
        return
      end if
c
      nlin = 0
      nodb = 0
      ierr = 0
c
c ... read loop
c
   10 continue
      read (iunit,'(a)',end=9000,err=9900) line
      nlin = nlin + 1
c
      if (length(line) .lt. 1) then
        write (*,*)
        goto 10
      end if
c
      if (line(1:1) .eq. '!') then
        write (*,'(1x,a)') line(1:leng1(line))
        goto 10
      end if
c
c ... if here, this must be the first line of a new
c     datablock
c
      call upcase (line)
      call extrop (line,nopt,4,pars,ierr)
      if (ierr .ne. 0) then
        call errcon ('While parsing datablock header')
        call textut (' Offending line :',line)
        call prompt (' Skipping rest of file !')
        goto 9000
      end if
c
      call str2i (pars(3),ndata,ierr)
      if (ierr .ne. 0) then
        call errcon ('While parsing datablock header')
        call textut (' Offending line :',line)
        call prompt (' Skipping rest of file !')
        goto 9000
      end if
c
      call allocm (pars(1),pars(2)(1:1),ndata,pars(4),
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('While allocating datablock space')
        call textut (' At line :',line)
        call prompt (' Skipping rest of file !')
        goto 9000
      end if
c
c ... okay, datablock allocated; now read it
c
      if (pars(2)(1:1) .eq. 'R') then
        read (iunit,pars(4),end=9900,err=9900)
     +    (rodb(i,j),i=1,ndata)
      else if (pars(2)(1:1) .eq. 'I') then
        read (iunit,pars(4),end=9900,err=9900)
     +    (iodb(i,j),i=1,ndata)
      else if (pars(2)(1:1) .eq. 'C') then
        read (iunit,pars(4),end=9900,err=9900)
     +    (codb(i,j),i=1,ndata)
      else
        lmax = 1
        do i=1,ndata
          read (iunit,'(a)',end=9900,err=9900)
     +      todb(i,j)
          lmax = max(lmax, length(todb(i,j)))
        end do
        write (pars(4),*) lmax
        call remspa (pars(4))
      end if
c
c ... okay
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .false.
      odbsel (iptr) = .false.
      odblen (iptr) = ndata
      odbtyp (iptr) = pars (2) (1:1)
      odbfmt (iptr) = pars (4)
      odbnam (iptr) = pars (1)
      odbcom (iptr) = 'Read from ' // file
c
      nodb = nodb + 1
      call pretty (line)
      call textut (' Read :',line)
c
      goto 10
c
c ... end of file (or no more allocation space)
c
 9000 continue
      close (iunit)
      call ivalut (' Nr of datablocks read :',1,nodb)
      call ivalut (' Nr of header/comments :',1,nlin)
c
      ierr = 0
      return
c
c ... error
c
 9900 continue
      call errcon ('While reading file')
      goto 9000
c
      end
c
c
c
      subroutine expodb (iunit,file,ierr)
c
      include 'odbman.incl'
c
      integer iunit,ierr,nodb,j,k,i,iptr,leng1
c
      character line*256,file*(*)
c
code ...
c
      if (iunit .ne. 6) then
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening ODB file')
          return
        end if
      end if
c
      nodb = 0
      ierr = 0
c
      if (iunit .ne. 6) then
        call stamp(line)
        line = '! ' // line
        write (iunit,'(a)',err=9900) line(1:leng1(line))
      end if
c
      do i=1,maxtdb
        if (odbuse(i) .and. odbsel(i)) then
          iptr = odbptr(i)
          write (line,*)
     +      odbnam(i)(1:leng1(odbnam(i))),
     +      ' ',odbtyp(i),' ',odblen(i),' ',
     +      odbfmt(i)
          call pretty (line)
          write (iunit,'(a)',err=9900) line(1:leng1(line))
          if (odbtyp(i) .eq. 'R') then
            write (iunit,odbfmt(i),err=9900)
     +        (rodb(j,iptr),j=1,odblen(i))
          else if (odbtyp(i) .eq. 'I') then
            write (iunit,odbfmt(i),err=9900)
     +        (iodb(j,iptr),j=1,odblen(i))
          else if (odbtyp(i) .eq. 'C') then
            write (iunit,odbfmt(i),err=9900)
     +        (codb(j,iptr),j=1,odblen(i))
          else
            call str2i (odbfmt(i),k,ierr)
            if (ierr .ne. 0) k = 72
            do j=1,odblen(i)
              write (iunit,'(a)',err=9900) todb(j,iptr)(1:k)
            end do
          end if
c
          if (iunit .ne. 6) call textut (' Written :',line)
          nodb = nodb + 1
          odbcha (i) = .false.
c
        end if
      end do
c
      ierr = 0
c
 9000 continue
      if (iunit .ne. 6) close (iunit)
      call ivalut (' Nr of datablocks written :',1,nodb)
      return
c
c ... write error
c
 9900 continue
      call errcon ('While writing file')
      ierr = -1
      goto 9000
c
      end
c
c
c
      subroutine setodb (iptr,i,first,last,step,valstr,ierr)
c
      include 'odbman.incl'
c
      real xdum
c
      integer i,iptr,first,last,step,ierr,idum,k
c
      character valstr*(*)
c
code ...
c
      if (odbtyp (iptr) .eq. 'R') then
        call str2r (valstr,xdum,ierr)
        if (ierr .ne. 0) return
        do k=first,last,step
          rodb (k,i) = xdum
        end do    
      else if (odbtyp (iptr) .eq. 'I') then
        call str2i (valstr,idum,ierr)
        if (ierr .ne. 0) return
        do k=first,last,step
          iodb (k,i) = idum
        end do
      else if (odbtyp (iptr) .eq. 'C') then
        do k=first,last,step
          codb (k,i) = valstr(1:6)
        end do
      else if (odbtyp (iptr) .eq. 'T') then
        do k=first,last,step
          todb (k,i) = valstr(1:maxlen)
        end do
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine setind (iptr,i,first,last,ierr)
c
      include 'odbman.incl'
c
      integer i,iptr,first,last,step,ierr,k,i1,i2,i3
c
      character line*80
c
code ...
c
      step = 4
c
      if (odbtyp (iptr) .eq. 'R') then
        do k=first,last,step
          i1 = k
          i2 = min (k+step-1,last)
          i3 = i2 - i1 + 1
          write (line,*) ' Entries ',i1,' - ',i2,' ?'
          call pretty (line)
          line = ' '//line
          call gvalin (line,i3,rodb(i1,i))
        end do
      else if (odbtyp (iptr) .eq. 'I') then
        do k=first,last,step
          i1 = k
          i2 = min (k+step-1,last)
          i3 = i2 - i1 + 1
          write (line,*) ' Entries ',i1,' - ',i2,' ?'
          call pretty (line)
          line = ' '//line
          call jvalin (line,i3,iodb(i1,i))
        end do
      else if (odbtyp (iptr) .eq. 'C') then
        do k=first,last
          write (line,*) ' Entry ',k,' ?'
          call pretty (line)
          line = ' '//line
          call textin (line,codb(k,i))
        end do
      else if (odbtyp (iptr) .eq. 'T') then
        do k=first,last
          write (line,*) ' Entry ',k,' ?'
          call pretty (line)
          line = ' '//line
          call textin (line,todb(k,i))
        end do
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine setif (i,iptr,setstr,j,jptr,oper,valstr,ierr)
c
      include 'odbman.incl'
c
      real xdum
c
      integer i,iptr,j,jptr,ierr,k,idum,nok
c
      logical eval
c
      character setstr*(*),valstr*(*),oper*(*)
c
code ...
c
ccc      print *,i,iptr
ccc      print *,setstr
ccc      print *,j,jptr
ccc      print *,oper
ccc      print *,valstr
c
      nok = 0
c
      if (odbtyp (i) .eq. 'R') then
        call str2r (setstr,xdum,ierr)
        if (ierr .ne. 0) return
        do k=1,odblen(i)
          if (eval(j,jptr,k,oper,valstr,ierr)) then
            rodb(k,iptr)=xdum
            nok = nok + 1
          end if
          if (ierr .ne. 0) return
        end do
      else if (odbtyp (i) .eq. 'I') then
        call str2i (setstr,idum,ierr)
        if (ierr .ne. 0) return
        do k=1,odblen(i)
          if (eval(j,jptr,k,oper,valstr,ierr)) then
            iodb(k,iptr)=idum
            nok = nok + 1
          end if
          if (ierr .ne. 0) return
        end do
      else if (odbtyp (i) .eq. 'C') then
        do k=1,odblen(i)
          if (eval(j,jptr,k,oper,valstr,ierr)) then
            codb(k,iptr)=setstr
            nok = nok + 1
          end if
          if (ierr .ne. 0) return
        end do
      else if (odbtyp (i) .eq. 'T') then
        do k=1,odblen(i)
          if (eval(j,jptr,k,oper,valstr,ierr)) then
            todb(k,iptr)=setstr
            nok = nok + 1
          end if
          if (ierr .ne. 0) return
        end do
      end if
c
      call jvalut (' Entries set :',1,nok)
      call jvalut (' Out of      :',1,odblen(i))
c
      ierr = 0
c
      return
      end
c
c
c
      logical function eval (j,jptr,k,oper,valstr,ierr)
c
      include 'odbman.incl'
c
      real xdum
c
      integer j,jptr,ierr,k,idum
c
      logical temp
c
      character valstr*(*),oper*(*)
c
code ...
c
      ierr = -1
      temp = .false.
c
      if (k .gt. odblen(j)) then
        call errcon ('EVAL -- out of bounds')
        goto 9999
      end if 
c
      if (odbtyp(j) .eq. 'R') then
        call str2r (valstr,xdum,ierr)
        if (ierr .ne. 0) goto 9999
        if (oper(1:1) .eq. '=') then
          temp = (xdum .eq. rodb(k,jptr))
        else if (oper(1:1) .eq. '>') then
          temp = (rodb(k,jptr) .gt. xdum)
        else if (oper(1:1) .eq. '<') then
          temp = (rodb(k,jptr) .lt. xdum)
        else if (oper(1:1) .eq. '#') then
          temp = (rodb(k,jptr) .ne. xdum)
        else
          goto 9999
        end if
        goto 9000
      else if (odbtyp(j) .eq. 'I') then
        call str2i (valstr,idum,ierr)
        if (ierr .ne. 0) goto 9999
        if (oper(1:1) .eq. '=') then
          temp = (idum .eq. iodb(k,jptr))
        else if (oper(1:1) .eq. '>') then
          temp = (iodb(k,jptr) .gt. idum)
        else if (oper(1:1) .eq. '<') then
          temp = (iodb(k,jptr) .lt. idum)
        else if (oper(1:1) .eq. '#') then
          temp = (iodb(k,jptr) .ne. idum)
        else
          goto 9999
        end if
        goto 9000
      else if (odbtyp(j) .eq. 'C') then
        if (oper(1:1) .eq. '=') then
          temp = (valstr(1:6) .eq. codb(k,jptr))
        else if (oper(1:1) .eq. '>') then
          temp = (codb(k,jptr) .gt. valstr(1:6))
        else if (oper(1:1) .eq. '<') then
          temp = (codb(k,jptr) .lt. valstr(1:6))
        else if (oper(1:1) .eq. '#') then
          temp = (codb(k,jptr) .ne. valstr(1:6))
        else
          goto 9999
        end if
        goto 9000
      else if (odbtyp(j) .eq. 'T') then
        if (oper(1:1) .eq. '=') then
          temp = (valstr(1:maxlen) .eq. todb(k,jptr))
        else if (oper(1:1) .eq. '>') then
          temp = (todb(k,jptr) .gt. valstr(1:maxlen))
        else if (oper(1:1) .eq. '<') then
          temp = (todb(k,jptr) .lt. valstr(1:maxlen))
        else if (oper(1:1) .eq. '#') then
          temp = (todb(k,jptr) .ne. valstr(1:maxlen))
        else
          goto 9999
        end if
        goto 9000
      end if
c
 9000 continue
      ierr = 0
c
 9999 continue
      eval = temp
c
      return
      end
c
c
c
      subroutine setkep (i,iptr,first,last,step,ierr)
c
      include 'odbman.incl'
c
      integer i,iptr,first,last,step,ierr,k,i1
c
code ...
c
      i1 = 0
      call jvalut (' Nr of entries before :',1,odblen(i))
c
      if (odbtyp (iptr) .eq. 'R') then
        do k=first,last,step
          i1 = i1 + 1
          rodb (i1,iptr) = rodb (k,iptr)
        end do
        odblen (i) = i1
      else if (odbtyp (iptr) .eq. 'I') then
        do k=first,last,step
          i1 = i1 + 1
          iodb (i1,iptr) = iodb (k,iptr)
        end do
        odblen (i) = i1
      else if (odbtyp (iptr) .eq. 'C') then
        do k=first,last,step
          i1 = i1 + 1
          codb (i1,iptr) = codb (k,iptr)
        end do
        odblen (i) = i1
      else if (odbtyp (iptr) .eq. 'T') then
        do k=first,last,step
          i1 = i1 + 1
          todb (i1,iptr) = todb (k,iptr)
        end do
        odblen (i) = i1
      end if
c
      call jvalut (' Nr of entries after  :',1,odblen(i))
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine exfiel (iunit,file,rdummy,idummy,type,
     +                   idum,i1,i2,i3,ierr)
c
      implicit none
c
      integer maxopt
      parameter (maxopt=50)
c
      real rdummy(*)
c
      integer idummy(*)
      integer i1,i2,i3,idum,iunit,ierr,nline,nopt,i,length
c
      logical xinter
c
      character file*(*),type*1
      character line*1024,optpar(maxopt)*80
c
code ...
c
      ierr = -1
      idum = 0
      nline = 0
c
      if (i1 .lt. 0 .or. i2 .lt. 0 .or.
     +    i3 .le. 0 .or. 
     +   (type.ne.'R' .and. type.ne.'I') ) then
        call errcon ('Invalid parameters')
        return
      end if
c
      call xopxoa (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) return
c
c ... skip header lines
c
      if (i1 .gt. 0) then
        do i=1,i1
          read (iunit,*,end=8010,err=8020)
          nline = nline + 1
        end do
      end if
c
   15 continue
      read (iunit,'(a)',end=25,err=8020) line
      nline = nline + 1
      if (length(line) .lt. 1) goto 8021
      call extrop (line,nopt,maxopt,optpar,ierr)
      if (ierr .ne. 0) goto 8022
      if (nopt .lt. i3) goto 8023
      idum = idum + 1
      if (type .eq. 'R') then
        call str2r (optpar(i3),rdummy(idum),ierr)
        if (ierr .ne. 0) goto 8025
      else 
        call str2i (optpar(i3),idummy(idum),ierr)
        if (ierr .ne. 0) goto 8024
      end if
c
c ... skip more lines if necessary
c
      if (i2 .gt. 0) then
        do i=1,i2
          read (iunit,*,end=25,err=8020)
          nline = nline + 1
        end do
      end if
c
      goto 15
c
c ... error traps
c
 8010 continue
      call errcon ('Unexpected end-of-file')
      call jvalut (' Near line :',1,nline)
      goto 25
c
 8020 continue
      call errcon ('Unexpected read error')
      call jvalut (' Near line :',1,nline)
      goto 25
c
 8021 continue
      call errcon ('Unexpected empty line')
      call jvalut (' Near line :',1,nline)
      goto 25
c
 8022 continue
      call errcon ('Unexpected error while parsing line')
      call jvalut (' Near line :',1,nline)
      call textut (' Line :',line)
      goto 25
c
 8023 continue
      call errcon ('Not enough items on line')
      call jvalut (' Near line :',1,nline)
      call textut (' Line :',line)
      goto 25
c
 8024 continue
      call errcon ('While converting item to integer')
      call jvalut (' Near line :',1,nline)
      call textut (' Line :',line)
      idum = idum - 1
      goto 25
c
 8025 continue
      call errcon ('While converting item to real')
      call jvalut (' Near line :',1,nline)
      call textut (' Line :',line)
      idum = idum - 1
      goto 25
c
c ... end of file reached
c
   25 continue
c
      call jvalut (' Number of lines read  :',1,nline)
      call jvalut (' Number of values read :',1,idum)
c
      close (iunit)
      ierr = 0
c
      return
      end
c
c
c
      subroutine exform (iunit,file,rdummy,idummy,type,
     +                   idum,i1,i2,myfmt,ierr)
c
      implicit none
c
      real rdummy(*)
c
      integer idummy(*)
      integer i1,i2,idum,iunit,ierr,nline,i,length
c
      logical xinter
c
      character file*(*),myfmt*(*),type*1
      character line*1024
c
code ...
c
      ierr = -1
      idum = 0
      nline = 0
c
      if (i1 .lt. 0 .or. i2 .lt. 0 .or.
     +    length(myfmt) .lt. 1 .or.
     +   (type.ne.'R' .and. type.ne.'I') ) then
        call errcon ('Invalid parameters')
        return
      end if
c
      call xopxoa (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) return
c
c ... skip header lines
c
      if (i1 .gt. 0) then
        do i=1,i1
          read (iunit,*,end=8010,err=8020)
          nline = nline + 1
        end do
      end if
c
   15 continue
      read (iunit,'(a)',end=25,err=8020) line
      nline = nline + 1
      if (length(line) .lt. 1) goto 8021
      idum = idum + 1
      if (type .eq. 'R') then
        read (line,myfmt,end=8022,err=8022) rdummy(idum)
      else
        read (line,myfmt,end=8022,err=8022) idummy(idum)
      end if
c
c ... skip more lines if necessary
c
      if (i2 .gt. 0) then
        do i=1,i2
          read (iunit,*,end=25,err=8020)
          nline = nline + 1
        end do
      end if
c
      goto 15
c
c ... error traps
c
 8010 continue
      call errcon ('Unexpected end-of-file')
      call jvalut (' Near line :',1,nline)
      goto 25
c
 8020 continue
      call errcon ('Unexpected read error')
      call jvalut (' Near line :',1,nline)
      goto 25
c
 8021 continue
      call errcon ('Unexpected empty line')
      call jvalut (' Near line :',1,nline)
      goto 25
c
 8022 continue
      call errcon ('Unexpected error reading item from line')
      call jvalut (' Near line :',1,nline)
      call textut (' Line :',line)
      idum = idum - 1
      goto 25
c
c ... end of file reached
c
   25 continue
c
      call jvalut (' Number of lines read  :',1,nline)
      call jvalut (' Number of values read :',1,idum)
c
      close (iunit)
      ierr = 0
c
      return
      end
c
c
c
      subroutine exmult (iunit,file,rdummy,idummy,nmax,type,
     +                   idum,i1,myfmt,ierr)
c
      implicit none
c
      integer nmax
c
      real rdummy(nmax)
c
      integer idummy(nmax)
      integer i1,idum,iunit,ierr,i,length
c
      logical xinter
c
      character file*(*),myfmt*(*),type*1
c
code ...
c
      ierr = -1
c
      if (i1 .lt. 0 .or.
     +    length(myfmt) .lt. 1 .or.
     +   (type.ne.'R' .and. type.ne.'I') ) then
        call errcon ('Invalid parameters')
        return
      end if
c
      call xopxoa (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) return
c
c ... skip header lines
c
      if (i1 .gt. 0) then
        do i=1,i1
          read (iunit,*,end=25,err=25)
        end do
      end if
c
      if (type .eq. 'R') then
        do i=1,nmax
          rdummy(i)=-620605.0
        end do
        read (iunit,myfmt,end=25,err=20)
     +    (rdummy(idum),idum=1,nmax)
      else
        do i=1,nmax
          idummy(i)=-620605
        end do
        read (iunit,myfmt,end=25,err=20)
     +    (idummy(idum),idum=1,nmax)
      end if
c
      call errcon (' Max nr of numbers read; rest skipped')
      goto 25
c
   20 continue
      call errcon ('While reading from file')
c
c ... end of file reached
c
   25 continue
c
      call prompt (' Counting number of values ...')
      idum = 0
      if (type .eq. 'R') then
        do i=nmax,1,-1
          if (rdummy(i) .ne. -620605.0) then
            idum = i
            goto 30
          end if
        end do
      else
        do i=nmax,1,-1
          if (idummy(i) .ne. -620605) then
            idum = i
            goto 30
          end if
        end do
      end if
c
   30 continue
      call jvalut (' Number of values read :',1,idum)
c
      close (iunit)
      ierr = 0
c
      return
      end
c
c
c
      subroutine exproc (iunit,molnam,file,ierr)
c
      include 'odbman.incl'
c
      real hbond(maxsiz)
c
      integer nrbadc(maxsiz)
      integer iunit,ierr,nlines,nres,i,j,iptr,ires
      integer length,leng1
c
      logical xinter
c
      character molnam*(*),file*(*)
      character resnam(maxsiz)*6,restyp(maxsiz)*6
      character kabsch(maxsiz)*6,ramach(maxsiz)*6
      character line*256,profil*80
c
code ...
c
      profil = file
c
      call xopxoa (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening file')
        return
      end if
c
      nlines = 0
      nres = 0
c
 6000 format (a)
 6100 format (1x,i6,4(1x,a6),1x,i3,1x,f5.1)
c
c ... read until " Residue     Kabsch Region"
c
   10 continue
      read (iunit,6000,err=9000,end=9000) line
      nlines = nlines + 1
      if (line(1:26) .ne.
     +    ' Residue     Kabsch Region') goto 10
c
c ... read 4 more lines
c
      do i=1,4
        read (iunit,6000,err=9000,end=9000) line
        nlines = nlines + 1
      end do
c
      call prompt (' Start of residue listing; extracting')
c
c ... now process residue lines
c
   20 continue
      read (iunit,6000,err=9000,end=9000) line
      nlines = nlines + 1
c
      if (line(1:13) .eq. '             ') goto 20
c
      if (length(line) .lt. 10) goto 20
c
      if (line(124:127) .eq. 'Page') then
        do i=1,8
          read (iunit,6000,err=9000,end=9000) line
          nlines = nlines + 1
        end do
        goto 20
      end if
c
      if (line(1:13) .eq. '-------------') goto 20
c
      if (line(1:15) .eq. 'Max deviations:') goto 100
c
c ... we have a residue !
c
      nres = nres + 1
c
c --> UPPERCASE
c
      call upcase (line)
c
      restyp (nres) = line (7:9)
      call remspa (restyp (nres))
      call upcase (restyp (nres))
c
      read (line(10:13),*) ires
      write (resnam(nres),'(a1,i5)')
     +  line(5:5),ires
      call remspa (resnam (nres))
      call upcase (resnam (nres))
c
      kabsch (nres) = line(17:17)
      call remspa (kabsch (nres))
c
      ramach (nres) = line(23:24)
      call remspa (ramach (nres))
c
      if (line(119:121) .eq. '  -') then
        nrbadc (nres) = 0
      else
        read (line(119:121),*) nrbadc (nres)
      end if
c
      if (line(102:105) .eq. '  - ') then
        hbond (nres) = 0.0
      else
        read (line(102:105),*) hbond (nres)
      end if
c
      write (*,6100) nres,restyp(nres),resnam(nres),kabsch(nres),
     +  ramach(nres),nrbadc(nres),hbond(nres)
c
      if (nres .eq. maxsiz) then
        call errcon ('Too many residues; rest skipped')
        goto 100
      end if
c
      goto 20
c
c ... done
c
  100 continue
      write (*,*)
      call prompt (' End of residue listing; skip rest of file')
      call jvalut (' Nr of lines read :',1,nlines)
      close (iunit)
c
      if (nres .lt. 1) then
        call errcon ('No residues found in file')
        return
      end if
c
c ... residue type
c
      file = molnam(1:leng1(molnam))//'_residue_type'
      call upcase (file)
c
      call allocm (file,'C',nres,'(1x,5a)',
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot allocate ODB; skipped')
        call textut (' ODB :',file)
        goto 200
      end if
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .true.
      odbsel (iptr) = .false.
      odblen (iptr) = nres
      odbtyp (iptr) = 'C'
      odbfmt (iptr) = '(1x,5a)'
      odbnam (iptr) = file
      odbcom (iptr) = 'Extracted from PROCHECK file '//profil
      do i=1,nres
        codb (i,j) = restyp(i)
      end do
c
c ... residue name
c
  200 continue
      file = molnam(1:leng1(molnam))//'_residue_name'
      call upcase (file)
c
      call allocm (file,'C',nres,'(1x,5a)',
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot allocate ODB; skipped')
        call textut (' ODB :',file)
        goto 300
      end if
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .true.
      odbsel (iptr) = .false.
      odblen (iptr) = nres
      odbtyp (iptr) = 'C'
      odbfmt (iptr) = '(1x,5a)'
      odbnam (iptr) = file
      odbcom (iptr) = 'Extracted from PROCHECK file '//profil
      do i=1,nres
        codb (i,j) = resnam(i)
      end do
c
c ... DSSP
c
  300 continue
      file = molnam(1:leng1(molnam))//'_residue_dssp'
      call upcase (file)
c
      call allocm (file,'C',nres,'(1x,5a)',
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot allocate ODB; skipped')
        call textut (' ODB :',file)
        goto 400
      end if
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .true.
      odbsel (iptr) = .false.
      odblen (iptr) = nres
      odbtyp (iptr) = 'C'
      odbfmt (iptr) = '(1x,5a)'
      odbnam (iptr) = file
      odbcom (iptr) = 'Extracted from PROCHECK file '//profil
      do i=1,nres
        codb (i,j) = kabsch(i)
      end do
c
c ... Ramachandran area
c
  400 continue
      file = molnam(1:leng1(molnam))//'_residue_rama'
      call upcase (file)
c
      call allocm (file,'C',nres,'(1x,5a)',
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot allocate ODB; skipped')
        call textut (' ODB :',file)
        goto 500
      end if
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .true.
      odbsel (iptr) = .false.
      odblen (iptr) = nres
      odbtyp (iptr) = 'C'
      odbfmt (iptr) = '(1x,5a)'
      odbnam (iptr) = file
      odbcom (iptr) = 'Extracted from PROCHECK file '//profil
      do i=1,nres
        codb (i,j) = ramach(i)
      end do
c
c ... Bad contacts
c
  500 continue
      file = molnam(1:leng1(molnam))//'_residue_badcon'
      call upcase (file)
c
      call allocm (file,'I',nres,'(1x,25i3)',
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot allocate ODB; skipped')
        call textut (' ODB :',file)
        goto 600
      end if
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .true.
      odbsel (iptr) = .false.
      odblen (iptr) = nres
      odbtyp (iptr) = 'I'
      odbfmt (iptr) = '(1x,25i3)'
      odbnam (iptr) = file
      odbcom (iptr) = 'Extracted from PROCHECK file '//profil
      do i=1,nres
        iodb (i,j) = nrbadc(i)
      end do
c
c ... H-bond energy
c
  600 continue
      file = molnam(1:leng1(molnam))//'_residue_hbond'
      call upcase (file)
c
      call allocm (file,'R',nres,'(1x,15f5.1)',
     +             j,iptr,ierr)
      if (ierr .ne. 0) then
        call errcon ('Cannot allocate ODB; skipped')
        call textut (' ODB :',file)
        goto 700
      end if
c
      odbuse (iptr) = .true.
      odbptr (iptr) = j
      odbcha (iptr) = .true.
      odbsel (iptr) = .false.
      odblen (iptr) = nres
      odbtyp (iptr) = 'R'
      odbfmt (iptr) = '(1x,15f5.1)'
      odbnam (iptr) = file
      odbcom (iptr) = 'Extracted from PROCHECK file '//profil
      do i=1,nres
        rodb (i,j) = hbond(i)
      end do
c
c ... done
c
  700 continue
      ierr = 0
      close (iunit)
      call prompt (' All done !')
c
      return
c
 9000 continue
      call errcon ('While reading PROCHECK file')
      close (iunit)
      return
c
      end
c
c
c
      subroutine odl3d (iunit,ndata,xdata,ydata,zdata,colour,ierr)
c
      implicit none
c
      real xdata(*),ydata(*),zdata(*)
      real ave,sdv,xmin,xmax,xtot,scale,radius
c
      integer iunit,ndata,ierr,i,leng1
c
      character line*128,colour*(*)
c
code...
c
      ierr = -1
c
c ... scale data to range 0-100
c
      call xstats (xdata,ndata,ave,sdv,xmin,xmax,xtot)
      scale = 100.0/(xmax-xmin)
      write (*,6100) 'X',xmin,xmax,scale
      do i=1,ndata
c        if (i.le.10) print *,'X ',i,xdata(i)
        xdata (i) = scale*(xdata(i)-xmin)
c        if (i.le.10) print *,'X ',i,xdata(i)
      end do
c
      call xstats (ydata,ndata,ave,sdv,xmin,xmax,xtot)
      scale = 100.0/(xmax-xmin)
      write (*,6100) 'Y',xmin,xmax,scale
      do i=1,ndata
c        if (i.le.10) print *,'Y ',i,ydata(i)
        ydata (i) = scale*(ydata(i)-xmin)
c        if (i.le.10) print *,'Y ',i,ydata(i)
      end do
c
      call xstats (zdata,ndata,ave,sdv,xmin,xmax,xtot)
      scale = 100.0/(xmax-xmin)
      write (*,6100) 'Z',xmin,xmax,scale
      do i=1,ndata
c        if (i.le.10) print *,'Z ',i,zdata(i)
        zdata (i) = scale*(zdata(i)-xmin)
c        if (i.le.10) print *,'Z ',i,zdata(i)
      end do
c
      write (iunit,6000,err=900) 'begin odb3d'
      write (iunit,6000,err=900)
     +  'colour ',colour(1:leng1(colour))
      write (iunit,6000,err=900)
     +  'text_colour ',colour(1:leng1(colour))
c
      write (iunit,6000,err=900)
     +  'move 0.0 0.0 0.0'
      write (iunit,6000,err=900)
     +  'line 100.0 0.0 0.0'
      write (iunit,6000,err=900)
     +  'text 100.0 0.0 0.0 X='
c
      write (iunit,6000,err=900)
     +  'move 0.0 0.0 0.0'
      write (iunit,6000,err=900)
     +  'line 0.0 100.0 0.0'
      write (iunit,6000,err=900)
     +  'text 0.0 100.0 0.0 Y='
c
      write (iunit,6000,err=900)
     +  'move 0.0 0.0 0.0'
      write (iunit,6000,err=900)
     +  'line 0.0 0.0 100.0'
      write (iunit,6000,err=900)
     +  'text 0.0 0.0 100.0 Z='
c
      write (iunit,6000,err=900) 'mode solid'
c
 6000 format (99a)
 6010 format ('sphere_xyz',3(1x,f8.3),1x,f8.1)
 6100 format (1x,a1,'-data range ',1p,2e12.4,' (scale by ',
     +        e12.4,')')
c
      radius = 0.5
c
      do i=1,ndata
        write (line,6010,err=900)
     +    xdata(i),ydata(i),zdata(i),radius
        call pretty (line)
        write (iunit,6000,err=900) line(1:leng1(line))
      end do
c
      write (iunit,6000,err=900) 'end_object'
      close (iunit)
c
      call prompt (' ODL file written')
      ierr = 0
c
      return
c
  900 continue
      call errcon ('While writing ODL file')
      close (iunit)
      ierr = -1
c
      return
      end


