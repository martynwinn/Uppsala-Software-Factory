      program cello
c
c === CELLO === play with your unit cell (parameters)
c
c Gerard Kleywegt @ 920714
c
c Version 921103
c
c ... note: for more conversions of rotation matrices, see:
c           R Diamond, Int Tables (new), vol B, ch 3.3
c
c ===========================================================================
c
      include 'cello.incl'
c
      real xyzf(3),xyzc(3),xyzs(3),xyzn(3),xyzt(3),xyzm(3)
      real det3
c
      integer length,iunit,ierror,i,j,k,i1,i2,is,ierr,leng1
c
      logical first
c
      character option*80,spgrpf*80,dcellf*80,dtype*1,ofmt*40
      character line*120,stampy*40
c
      data iunit  /99/
      data option /'redefine_cell_constants'/
      data spgrpf /'spgrp.o'/
      data dcellf /'cell.o'/
      data ofmt /'(3f15.7)'/
      data xyzf, xyzc /6*0.0/
c
code ...
c
      call gkinit (prognm,version)
c
      write (tty,*) '... Initialising ...'
c
c ... initialise program
c
      call lm_init
c
      ncsfil = ' '
c
c ... get relevant data (cell parameters, symm-ops)
c
cc      call get_data
c
      first = .true.
      option = 'L'
      goto 111
c
c ... MAIN LOOP
c
   99 continue
      if (first) first = .false.
      write (tty,'(/99(1x,a/:))')
     +  'Select an option :',
     +  '==================',
     +  ' ',
     +  '?                           Quit',
     +  'Cartesian_to_fractional     Fractional_to_Cartesian',
     +  'List_cell_constants         Redefine_cell_constants',
     +  'Spacegroup_datablock_for_O  Draw_cell_datablock_for_O',
     +  'Import_NCS_data             Export_NCS_data',
     +  'Apply_symmetry_operators'
c
  100 continue
      if (first) goto 99
      write (tty,*)
      call textin (' Option ?',option)
      if (option(1:1) .eq. '!')  goto 100
      if (option(1:1) .eq. '?')  goto  99
      if (length(option) .lt. 1) goto 100
      call upcase (option)
c
  111 continue
c
c Q   = QUIT
c
      if (option(1:1) .eq. 'Q') then
c
        goto 999
c
c C   = CARTESIAN_TO_FRACTIONAL
c
      else if (option(1:1) .eq. 'C') then
c
        call fvalin (' Cartesian coordinates  ?',3,xyzc)
c
        call ca2fra (xyzc,xyzf)
c
        call fvalut (' Fractional coordinates :',3,xyzf)
c
c F   = FRACTIONAL_TO_CARTESIAN
c
      else if (option(1:1) .eq. 'F') then
c
        call fvalin (' Fractional coordinates ?',3,xyzf)
c
        call fra2ca (xyzf,xyzc)
c
        call fvalut (' Cartesian coordinates  :',3,xyzc)
c
c A   = APPLY_SYMMETRY_OPERATORS
c
      else if (option(1:1) .eq. 'A') then
c
        call fvalin (' Fractional coordinates ?',3,xyzf)
c
        do i=1,nsym
c
          xyzt (1) = xyzf(1)*symmat(1,1,i) + xyzf(2)*symmat(1,2,i) + 
     +               xyzf(3)*symmat(1,3,i) + symmat(1,4,i)
          xyzt (2) = xyzf(1)*symmat(2,1,i) + xyzf(2)*symmat(2,2,i) + 
     +               xyzf(3)*symmat(2,3,i) + symmat(2,4,i)
          xyzt (3) = xyzf(1)*symmat(3,1,i) + xyzf(2)*symmat(3,2,i) + 
     +               xyzf(3)*symmat(3,3,i) + symmat(3,4,i)
c
          call fra2ca (xyzt,xyzs)
c
          do j=1,ncs
            xyzn (1) = xyzs(1)*ncsmat(1,1,j) + xyzs(2)*ncsmat(1,2,j) +
     +                 xyzs(3)*ncsmat(1,3,j) + ncsmat(1,4,j)
            xyzn (2) = xyzs(1)*ncsmat(2,1,j) + xyzs(2)*ncsmat(2,2,j) +
     +                 xyzs(3)*ncsmat(2,3,j) + ncsmat(2,4,j)
            xyzn (3) = xyzs(1)*ncsmat(3,1,j) + xyzs(2)*ncsmat(3,2,j) +
     +                 xyzs(3)*ncsmat(3,3,j) + ncsmat(3,4,j)
c
            call ca2fra (xyzn,xyzm)
c
            write (tty,'(1x,a,i2,a,i2,a,3f10.3)')
     +        'Fract coords SYM # ',i,' & NCS # ',j,' : ',
     +        (xyzm(k),k=1,3)
            write (tty,'(1x,a,i2,a,i2,a,3f10.3)')
     +        'Cart  coords SYM # ',i,' & NCS # ',j,' : ',
     +        (xyzn(k),k=1,3)
c
          end do
        end do
c
c R   = REDEFINE_CELL_CONSTANTS
c
      else if (option(1:1) .eq. 'R') then
c
        call get_data
c
c L   = LIST_CELL_CONSTANTS
c
      else if (option(1:1) .eq. 'L') then
c
        call fvalut (' Cell axes    (A) :',3,axes)
        call fvalut (' Angles (degrees) :',3,angles)
        call jvalut (' Nr SYM operators :',1,nsym)
        call jvalut (' Nr NCS operators :',1,ncs)
c
c        do is=1,nsym
c          write (*,*)
c          write (tty,'(a7,i2,a3,3(3f13.3,10x,f13.3,:,/12x))')
c     +     ' SYM # ',is,' : ',((symmat(i,j,is),j=1,4),i=1,3)
c          call anamat (symmat(1,1,is),sdet(is))
c        end do
c
        call anasgs (nsym,symmat,.true.,ierr)
        do is=1,nsym
          sdet (is) = det3 (symmat(1,1,is))
        end do
c
c        do is=1,ncs
c          write (*,*)
c          write (tty,'(a7,i2,a3,3(3f13.7,10x,f13.4,:,/12x))')
c     +     ' NCS # ',is,' : ',((ncsmat(i,j,is),j=1,4),i=1,3)
c          call anamat (ncsmat(1,1,is),ndet(is))
c          ndet (is) = det3(ncsmat(1,1,is))
c          if (is .gt. 1) call matana (ncsmat(1,1,is))
c        end do
c
        call anancs (ncs,ncsmat,.true.,ierr)
        do is=1,nsym
          ndet (is) = det3 (ncsmat(1,1,is))
        end do
c
c S   = SPACEGROUP_DATABLOCK_FOR_O
c
      else if (option(1:1) .eq. 'S') then
c
        call textin (' Filename ?',spgrpf)
        close (iunit)
        call xopxua (iunit,spgrpf,.true.,ierror)
        if (ierror .ne. 0) goto 100
c
        write (iunit,'(a,i4,a)') 
     +    '.space_group_operators  R ',12*nsym,'  (3f10.2)'
        write (iunit,'(3f10.2)')
     +    (((symmat(i,j,k),i=1,3),j=1,4),k=1,nsym)
        close (iunit)
c
c D   = DRAW_CELL_DATABLOCK_FOR_O
c
      else if (option(1:1) .eq. 'D') then
c
        call textin (' Filename ?',dcellf)
        close (iunit)
        call xopxua (iunit,dcellf,.true.,ierror)
        if (ierror .ne. 0) goto 100
c
        write (iunit,'(a)') 'CELL  T  24  40'
c
        write (iunit,'(a/a)') ' begin cell',' colour 16799999'
c
 6000 format (' t ',3f10.3,' ',a)
 6010 format (' ',a1,' ',3f10.3)
c
        call fra3ca (0.1,0.0,0.0,xyzs)
        write (iunit,6000) (xyzs(j),j=1,3),'X'
        call fra3ca (0.0,0.1,0.0,xyzs)
        write (iunit,6000) (xyzs(j),j=1,3),'Y'
        call fra3ca (0.0,0.0,0.1,xyzs)
        write (iunit,6000) (xyzs(j),j=1,3),'Z'
c
        call fra3ca (0.0, 0.0, 0.0, xyzs)
        write (iunit,6010) 'm', (xyzs(j),j=1,3)
        call fra3ca (1.0, 0.0, 0.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
        call fra3ca (1.0, 1.0, 0.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
        call fra3ca (0.0, 1.0, 0.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
        call fra3ca (0.0, 0.0, 0.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
c
        call fra3ca (0.0, 0.0, 1.0, xyzs)
        write (iunit,6010) 'm', (xyzs(j),j=1,3)
        call fra3ca (1.0, 0.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
        call fra3ca (1.0, 1.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
        call fra3ca (0.0, 1.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
        call fra3ca (0.0, 0.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
c
        call fra3ca (0.0, 0.0, 0.0, xyzs)
        write (iunit,6010) 'm', (xyzs(j),j=1,3)
        call fra3ca (0.0, 0.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
c
        call fra3ca (1.0, 0.0, 0.0, xyzs)
        write (iunit,6010) 'm', (xyzs(j),j=1,3)
        call fra3ca (1.0, 0.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
c
        call fra3ca (1.0, 1.0, 0.0, xyzs)
        write (iunit,6010) 'm', (xyzs(j),j=1,3)
        call fra3ca (1.0, 1.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
c
        call fra3ca (0.0, 1.0, 0.0, xyzs)
        write (iunit,6010) 'm', (xyzs(j),j=1,3)
        call fra3ca (0.0, 1.0, 1.0, xyzs)
        write (iunit,6010) 'l', (xyzs(j),j=1,3)
c
        write (iunit,'(a)') ' end_object'
c
        close (iunit)
c
c I   = IMPORT_NCS_DATA
c
      else if (option(1:1) .eq. 'I') then
c
 6100 format (3f13.7,10x,f13.4)
c
        call textin (
     +    ' X-PLOR-file, O-file or Manual input (X/O/M) ?',
     +    mode)
        call upcase (mode)
        if (index ('XOM',mode) .le. 0) mode = 'M'
c
        if (mode .eq. 'M') then
c
          call jvalin (' Nr of NCS operators ?',1,ncs)
          call chkdim (ncs,1,maxncs,'MAIN',
     +      'NCS - nr of NCS operators')
c
          do is=1,ncs
            if (is .eq. 1) then
              write (tty,*) 'Unit operator fixed'
            else
              write (tty,*) 'Enter matrix # ',is,' ...'
              call fvalin (' Operator ?',12,ncsmat(1,1,is))
            end if
            write (tty,6100) ((ncsmat(i,j,is),j=1,4),i=1,3)
            call textin (' Operator name ?',ncsnam(is))
          end do
c
        else if (mode .eq. 'O') then
c
          ncs = 0
          j = 4
cccc          write (tty,*) 'Unit operator fixed'
c
   23     continue
          if (j .ne. 4) ncs = ncs - 1
          if (ncs .ge. maxncs) goto 25
          ncsfil = ' '
          call textin (' O filename (RETURN to stop) ?',ncsfil)
c
          if (length(ncsfil) .lt. 1) goto 25
          close (iunit)
          call xopxoa (iunit,ncsfil,.true.,ierror)
          if (ierror .ne. 0) goto 100
c
   26     continue
          ncs = ncs + 1
          j = 0
c
   24     continue
          read (iunit,'(a)',END=23) line
          if (line(1:1) .eq. '!') goto 24
c
          call jvalut (' NCS operator    :',1,ncs)
c
          call pretty (line)
          call upcase (line)
c
          i1 = 1
          i2 = index (line,' ')
          ncsnam (ncs) = line(i1:i2-1)
          call upcase (ncsnam(ncs))
          line = line (i2+1:)
c
          i1 = 1
          i2 = index (line,' ')
          dtype = line(i1:i2-1)
          line = line (i2+1:)
c
          i1 = 1
          i2 = index (line,' ')
          read (line(i1:i2-1),*) i
          line = line (i2+1:)
c
          ofmt = line (1:leng1(line))
          call remspa (ofmt)
c
          call textut (' Data block name :',ncsnam(ncs))
          call textut (' Format          :',ofmt)
c
          read (iunit,ofmt) ((ncsmat(i,j,ncs),i=1,3),j=1,4)
          write (tty,6100) ((ncsmat(i,j,ncs),j=1,4),i=1,3)
          j = 4
c
          goto 26
c
   25     continue
          call jvalut (' Nr of NCS operators read :',1,ncs)
c
          close (iunit)
c
        else if (mode .eq. 'X') then
c
          ncs = 0
          call textin (' XPLOR filename ?',ncsfil)
          close (iunit)
          call xopxoa (iunit,ncsfil,.true.,ierror)
          if (ierror .ne. 0) goto 100
c
   33     continue
          read (iunit,'(a)',end=35) line
          i = length(line)
          if (i .lt. 1) goto 33
          call upcase (line)
c
          if (index (line,'XNCSREL') .le. 0) goto 33
c
c ... found a new XNCSREL line
c
          ncs = ncs + 1
          j = 0
c
   37     continue
          read (iunit,'(a)',end=35) line
          i = length(line)
          if (i .lt. 1) goto 37
c
          i1 = index (line,'(')
          i2 = index (line,')')
          if (i1 .lt. 1 .or. i2 .lt. 1 .or. i2 .le. i1) goto 37
c
cc          call textut (' Line :',line)
c
          j = j + 1
          if (j .le. 3) then
            read (line(i1+1:i2-1),*) (ncsmat(j,i,ncs),i=1,3)
          else if (j .eq. 4) then
            read (line(i1+1:i2-1),*) (ncsmat(i,4,ncs),i=1,3)
          end if
          if (j .eq. 4) goto 38
          goto 37
c
   38     continue
          call jvalut (' NCS operator :',1,ncs)
          write (tty,6100) ((ncsmat(i,j,ncs),j=1,4),i=1,3)
          write (ncsnam(ncs),'(a,i6)') '.LSQ_RT_',ncs
          call remspa (ncsnam(ncs))
          call textut (' Operator name :',ncsnam(ncs))
          goto 33
c
c ... end
c
   35     continue
          if (j .ne. 4) ncs = ncs - 1
          call jvalut (' Nr of NCS operators read :',1,ncs)
c
          close (iunit)
c
        end if
c
c E   = EXPORT_NCS_DATA
c
      else if (option(1:1) .eq. 'E') then
c
        if (ncs .le. 0) then
          call errcon ('There are no NCS operators in memory !')
          goto 100
        end if
c
        call textin (
     +    ' X-PLOR-file, O-file or Merely to the screen (X/O/M) ?',
     +    mode)
        call upcase (mode)
        if (index ('XOM',mode) .le. 0) mode = 'M'
c
        if (mode .eq. 'M') then
c
          do is=1,ncs
            call textut (' Operator :',ncsnam(is))
            write (tty,6100) ((ncsmat(i,j,is),j=1,4),i=1,3)
          end do
c
        else if (mode .eq. 'O') then
c
          ncsfil = ' '
          call textin (' O filename ?',ncsfil)
          close (iunit)
          call xopxua (iunit,ncsfil,.true.,ierror)
          if (ierror .ne. 0) goto 100
          call textin (' Format ?',ofmt)
          stampy = ' '
          call stamp (stampy)
          write (iunit,'(a,a)') '! ',stampy(1:leng1(stampy))
c
          do is=1,ncs
c
            call textut (' Operator   :',ncsnam(is))
c
            write (iunit,'(a)') '!'
            write (line,'(a,1x,a,1x,i4,1x,a)')
     +        ncsnam(is),'R',12,ofmt
            call pretty (line)
            write (iunit,'(a)') line(1:leng1(line))
c
            write (iunit,fmt=ofmt) ((ncsmat(i,j,is),i=1,3),j=1,4)
c
          end do
c
          close (iunit)
c
        else if (mode .eq. 'X') then
c
          ncsfil = ' '
          call textin (' XPLOR filename ?',ncsfil)
          close (iunit)
          call xopxua (iunit,ncsfil,.true.,ierror)
          if (ierror .ne. 0) goto 100
c
          stampy = ' '
          call gkdate (stampy)
c
          write (line,'(9a)') ' Remarks NCS File : ',ncsfil
          write (iunit,'(a)') line(1:leng1(line))
c
          write (line,'(9a)') ' { Created at : ',
     +      stampy(1:leng1(stampy)),' }'
          write (iunit,'(a)') line(1:leng1(line))
c
          write (line,'(9a)') ' { By program : ',prognm,' }'
          write (iunit,'(a)') line(1:leng1(line))
c
          write (line,'(9a)') ' { Version : ',version,' }'
          write (iunit,'(a)') line(1:leng1(line))
c
          write (line,'(a)')
     +      ' { Invoke strict non-crystallographic symmetry }'
          write (iunit,'(a)') line(1:leng1(line))
c
          write (iunit,*)
          write (iunit,*) ' ncs strict'
          write (iunit,*)
          write (iunit,*) ' { ==> Assuming identity skew matrix }'
          write (iunit,*)
          write (iunit,*) '  skew'
          write (iunit,*)
     +      '    matrix = ( 1.000000 0.000000 0.000000 )'
          write (iunit,*)
     +      '             ( 0.000000 1.000000 0.000000 )'
           write (iunit,*)
     +      '             ( 0.000000 0.000000 1.000000 )'
          write (iunit,*)
     +      '    translation = ( 0.0000 0.0000 0.0000 )'
          write (iunit,*) '  end'
          write (iunit,*)
c
          do j=1,ncs
            write (iunit,'(a,i3,a)') '   xncsrel { # ',j,' }'
            write (iunit,'(a,3f13.6,a)')
     +        '    matrix      = ( ',(ncsmat(1,i,j),i=1,3),' )'
            write (iunit,'(a,3f13.6,a)')
     +        '                  ( ',(ncsmat(2,i,j),i=1,3),' )'
            write (iunit,'(a,3f13.6,a)')
     +        '                  ( ',(ncsmat(3,i,j),i=1,3),' )'
            write (iunit,'(a,3f13.4,a)')
     +        '    translation = ( ',(ncsmat(i,4,j),i=1,3),' )'
            write (iunit,*) '  end'
            write (iunit,*)
          end do
c
          write (iunit,*)
          write (iunit,*) '  { xncsrel end }'
          write (iunit,*)
          write (iunit,*) '  ?'
          write (iunit,*)
          write (iunit,*) ' end'
c
          close (iunit)
c
        end if
c
c INVALID OPTION
c
      else
        call errcon ('Invalid option')
        call textut (' Option :',option)
      end if
c
      goto 100
c
c ... THE END
c
  999 continue
      call gkquit
c
      stop
      end
c
c ===========================================================================
c
      subroutine lm_init
c
      include 'cello.incl'
c
      integer is,i,j
c
c ... data statements to initialise
c
      data alpha,beta,gamma /3*90.0/, a1,a2,a3 /3*100.0/
c
code ...
c
c
      nsym = 1
      do is=1,maxsym
        sdet (is) = 1.0
        do i=1,3
          do j=1,4
            symmat (i,j,is) = 0.0
          end do
          symmat (i,i,is) = 1.0
        end do
      end do
c
      mode = 'M'
      ncs  = 1
      do is=1,maxncs
        ndet (is) = 1.0
        write (ncsnam(is),'(a,i5)') 'ncs_operator_',is
        call remspa (ncsnam(is))
        do i=1,3
          do j=1,4
            ncsmat (i,j,is) = 0.0
          end do
          ncsmat (i,i,is) = 1.0
        end do
      end do
      ncsnam (1) = 'unit_operator'
c
      write (tty,1200)
      write (tty,1000) 'Current version',version
      write (tty,1100) 'Max nr of symmetry operators',maxsym
      write (tty,1100) 'Max nr of NCS operators',maxncs
      write (tty,1200)
c
      return
c
 1000 format (1x,a30,' : ',a)
 1100 format (1x,a30,' : ',i8)
 1200 format (/' ***** ',4('CELLO ***** ')/)
c
      end
c
c ===========================================================================
c
      subroutine get_data
c
c ... get problem-related data
c
      include 'cello.incl'
c
      real arg,dom
c
      integer is,j,i,iunit,ierror,i1,i2,k,leng1
c
      character smode*1,ofile*80,line*120,ofmt*40
c
      data smode /'O'/, ofile /'c2.o'/, iunit /98/
c
      save smode,ofile
c
code ...
c
      call fvalin (' Axes (A) ?',3,axes)
      write (tty,6911) 'a,b,c :',a1,a2,a3
c
      call fvalin (' Angles (deg) ?',3,angles)
      write (tty,6911) 'alpha,beta,gamma :',alpha,beta,gamma
c
c ... compute some values for conversion FRACT <-> COORD
c
      CA = COS(ALPHA*DEGTOR)
      SA = SIN(ALPHA*DEGTOR)
      CB = COS(BETA *DEGTOR)
      SB = SIN(BETA *DEGTOR)
      CG = COS(GAMMA*DEGTOR)
      SG = SIN(GAMMA*DEGTOR)
      CARTYZ = (CA-CB*CG)/SG
      arg = SQRT(one-CA*CA-CB*CB-CG*CG+two*CA*CB*CG)
      CARTZZ = arg/SG
      volume = a1*a2*a3*arg
      astar = (a2*a3*sa)/volume
      bstar = (a1*a3*sb)/volume
      cstar = (a1*a2*sg)/volume
      car = (cb*cg-ca)/(sb*sg)
      cbr = (ca*cg-cb)/(sa*sg)
      cgr = (ca*cb-cg)/(sa*sb)
c
      dom = sqrt(one-car*car)
      a11 = one/a1
      a12 = -cg/(sg*a1)
      a13 = -(cg*sb*ca+cb*sg)/(sb*dom*sg*a1)
      a22 = one/(sg*a2)
      a23 = car/(dom*sg*a2)
      a33 = one/(sb*dom*a3)
c
 6911 format (1x,a25,3x,5f10.2)
c
  100 continue
      call textin (' Symmetry operators: Manual or O-file (M/O) ?',
     +  smode)
      call upcase (smode)
      if (smode .ne. 'O') smode = 'M'
c
      if (smode .eq. 'M') then
c
        call ivalin (' Nr of symmetry operators ?',1,nsym)
        write (tty,*) nsym
        call chkdim (nsym,1,maxsym,'GET_DATA',
     +    'NSYM - nr of symmetry operators')
c
        write (tty,*) 'Enter symmetry operators ...'
        do is=1,nsym
          if (is .eq. 1) then
            write (tty,*) ' Unit operator fixed'
          else
            write (tty,*) 'Enter matrix # ',is,' ...'
            call fvalin (' Operator ?',12,symmat(1,1,is))
          end if
          write (tty,'(3f8.2,10x,f8.2)') ((symmat(i,j,is),j=1,4),i=1,3)
        end do
c
      else
c
        nsym = 1
        call textin (' O file name ?',ofile)
        close (iunit)
        call xopxoa (iunit,ofile,.true.,ierror)
        if (ierror .ne. 0) goto 100
c
        read (iunit,'(a)') line
        call pretty (line)
        call upcase (line)
c
        i1 = 1
        i2 = index (line,' ')
        call textut (' Data block name :',line(i1:i2-1))
        line = line (i2+1:)
c
        i1 = 1
        i2 = index (line,' ')
        line = line (i2+1:)
c
        i1 = 1
        i2 = index (line,' ')
        read (line(i1:i2-1),*) nsym
        nsym = nsym / 12
        call jvalut (' Nr symmops      :',1,nsym)
        line = line (i2+1:)
c
        ofmt = line (1:leng1(line))
        call remspa (ofmt)
        call textut (' Format          :',ofmt)
c
        read (iunit,ofmt) (((symmat(i,j,k),i=1,3),j=1,4),k=1,nsym)
        write (tty,'(3f8.2,10x,f8.2)')
     +    (((symmat(i,j,k),i=1,3),j=1,4),k=1,nsym)
c
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine ca2fra (xyzc,xyzf)
c
c ... Cartesian to Fractional
c
      include 'cello.incl'
c
      real xyzc(3),xyzf(3)
c
c
code ...
c
      xyzf (1) = a11*xyzc(1) + a12*xyzc(2) + a13*xyzc(3)
      xyzf (2) =               a22*xyzc(2) + a23*xyzc(3)
      xyzf (3) =                             a33*xyzc(3)
c
      return
      end
c
c ===========================================================================
c
      subroutine fra2ca (xyzf,xyzc)
c
c ... Fractional to Cartesian
c
      include 'cello.incl'
c
      real xyzc(3),xyzf(3)
c
c
code ...
c
      xyzc (1) = a1*xyzf(1) + a2*xyzf(2)*cg + a3*xyzf(3)*cb
      xyzc (2) =              a2*xyzf(2)*sg + a3*xyzf(3)*cartyz
      xyzc (3) =                              a3*xyzf(3)*cartzz
c
      return
      end
c
c ===========================================================================
c
      subroutine fra3ca (xf,yf,zf,xyzc)
c
c ... Fractional to Cartesian
c
      include 'cello.incl'
c
      real xyzc(3),xf,yf,zf
c
c
code ...
c
      xyzc (1) = a1*xf + a2*yf*cg + a3*zf*cb
      xyzc (2) =         a2*yf*sg + a3*zf*cartyz
      xyzc (3) =                    a3*zf*cartzz
c
      return
      end
c
c ===========================================================================
c
      subroutine matana (q)
c
c ... fudged version specifically for CELLO
c
      implicit none
c
      real twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      real r(3,3),det,det3,theta,dtheta,v(3),vlen,p(3),e(3),q(3,3)
      real xsign,qqq
c
      integer i,j
c
code ...
c
      do i=1,3
        do j=1,3
          r(i,j)=q(j,i)
        end do
      end do
c
      det = det3(r)
c
c ... find POLAR angles from matrix
c
      qqq = max (-1.0, min (1.0, (r(1,1)+r(2,2)+r(3,3)-1.0)*0.5))
      theta = acos ( qqq )
      dtheta = theta*rtodeg
ccc      if (dtheta .gt. 180.0) dtheta = 360.0 - dtheta
c
c ... changed 940316
c
      xsign = 1.0
      if (dtheta .gt. 180.0) then
        dtheta = 360.0 - dtheta
        xsign = -1.0
      end if
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
        do i=1,3
          v(i) = v(i)*vlen
        end do
        write (*,'(a,3f10.4)')
     +    ' Direction cosines of rotation axis : ',
     +    (v(i),i=1,3)
      end if
c
c ... changed 940316
c
      det = sqrt( abs(1.0-v(3)*v(3)) )
      p(1) = atan2 (det,v(3)) * rtodeg
ccc      p(1) = acos (v(3)) * rtodeg
      qqq = max (-1.0, min (1.0, v(3)))
      p(1) = acos ( qqq ) * rtodeg
      p(2) = atan2 (v(2),v(1)) * rtodeg 
      p(3) = dtheta
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
      return
      end
