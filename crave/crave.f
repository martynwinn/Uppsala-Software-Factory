      program crave
c
c ... CRAVE - generate C-shell script to do RAVE multiple-crystal averaging
c
c     Usage: CRAVE < input_file
c
c     Gerard Kleywegt @ 960409/10
c
c     f77 -o CRAVE jiffy.f ../gklib/alpha_kleylib
c
c     f77 -Olimit 3000 -C -O -u -v -recursive -check underflow -check overflow -o CRAVE jiffy.f ../gklib/alpha_kleylib 
c
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'CRAVE', vers = '960415/0.03')
c
c ... MAXFRM = max nr of crystal forms
c     MAXNCS = max nr of NCS operators per crystal form
c     MAXSTR = max string length for filenames etc.
c     MAXKEY = max nr of general keywords
c
      integer maxfrm,maxncs,maxstr,maxkey
      parameter (maxfrm=25, maxncs=100, maxstr=256, maxkey=50)
c
      integer numncs(maxncs)
      integer iunit,isfall,irstats,ifft,imave,icomap,imapman
      integer iwork,iden,icycle,iscra,imask,ifile,maxi,jcopy
      integer i,jmtz,jmap,jlabf,jsigf,jgrid,jexte,jreso,josym
      integer jxxrt,maxj,j,ncopy,length,ncycle,mcycle,ierr,leng1
c
      logical xinter
c
      character key(maxkey)*4,value(maxkey)*(maxstr)
      character ckey(maxkey)*4,cvalue(maxfrm,maxkey)*(maxstr)
      character ncsfil(maxncs,maxfrm)*(maxstr)
      character*(maxstr) line
c
code ...
c
      call gkinit (prognm,vers)
c
      call jvalut (' Max nr of crystal forms    :',1,maxfrm)
      call jvalut (' Max nr of NCS-ops per form :',1,maxncs)
c
      iunit = 11
c
      isfall = 1
      key (1) = 'SFAL'
      value (1) = 'sfall'
c
      irstats = 2
      key (2) = 'RSTA'
      value (2) = 'rstats'
c
      ifft = 3
      key (3) = 'FFT '
      value (3) = 'fft'
c
      imave = 4
      key (4) = 'MAVE'
      value (4) = 'mave'
c
      icomap = 5
      key (5) = 'COMA'
      value (5) = 'comap'
c
      ifile = 6
      key (6) = 'FILE'
      value (6) = 'xxtal.csh'
c
      imapman = 7
      key (7) = 'MAPM'
      value (7) = 'mapman'
c
      iwork = 8
      key (8) = 'WORK'
      value (8) = '.'
c
      iden = 9
      key (9) = 'IDEN'
      value (9) = 'xxtal'
c
      icycle = 10
      key (10) = 'CYCL'
      value (10) = '10'
c
      iscra = 11
      key (11) = 'SCRA'
      value (11) = '$CCP4_SCR'
c
      imask = 12
      key (12) = 'MASK'
      value (12) = 'm1.mask'
c
c
c
      maxi = 12
      call asciut (' Main keywords :',maxi,key)
c
c
c
      jcopy = 1
      ckey (1) = 'COPY'
      do i=1,maxfrm
        write (cvalue(i,1),'(a4,i3)') 'form',i
        call remspa (cvalue(i,1))
      end do
c
      jmtz = 2
      ckey (2) = 'MTZ '
      do i=1,maxfrm
        cvalue(i,2) = '?'
      end do
c
      jmap = 3
      ckey (3) = 'MAP '
      do i=1,maxfrm
        cvalue(i,3) = '?'
      end do
c
      jlabf = 4
      ckey (4) = 'LABF'
      do i=1,maxfrm
        cvalue(i,4) = 'F'
      end do
c
      jsigf = 5
      ckey (5) = 'SIGF'
      do i=1,maxfrm
        cvalue(i,5) = 'SIGF'
      end do
c
      jgrid = 6
      ckey (6) = 'GRID'
      do i=1,maxfrm
        cvalue(i,6) = '100 100 100'
      end do
c
      jexte = 7
      ckey (7) = 'EXTE'
      do i=1,maxfrm
        cvalue(i,7) = '0 99 0 99 0 99'
      end do
c
      jreso = 8
      ckey (8) = 'RESO'
      do i=1,maxfrm
        cvalue(i,8) = '8.0 2.5'
      end do
c
      josym = 9
      ckey (9) = 'OSYM'
      do i=1,maxfrm
        cvalue(i,9) = 'p1.sym'
      end do
c
      jxxrt = 10
      ckey (10) = 'XXRT'
      do i=1,maxfrm
        cvalue(i,10) = 'rt_unit.o'
      end do
c
c
c
      maxj = 10
      ckey (maxj+1) = 'NCSO'
      call asciut (' Form keywords :',maxj+1,ckey)
c
c
c
      do i=1,maxncs
        do j=1,maxfrm
          ncsfil (i,j) = ' '
        end do
      end do
c
c
c
      ncopy = 0
c
   10 continue
      read (*,'(a)',end=1000) line
ccc      call textut (' Line >',line)
      if (line(1:1) .eq. '!') goto 10
c
      do i=1,maxi
        if (line(1:4) .eq. key(i)) then
          value (i) = line(6:)
          goto 10
        end if
      end do
c
      do i=1,maxj
        if (line(1:4) .eq. ckey(i)) then
          if (i .eq. jcopy) then
            ncopy = ncopy + 1
            call jvalut (' Crystal form :',1,ncopy)
            if (ncopy .gt. maxfrm) then
              call errstp ('Too many crystal forms; recompile')
            end if
          end if
          cvalue (ncopy,i) = line(6:)
          goto 10
        end if
      end do
c
      if (line(1:4) .eq. 'NCSO') then
        do i=1,maxncs
          if (length(ncsfil(i,ncopy)) .lt. 1) then
            ncsfil (i,ncopy) = line(6:)
            goto 10
          end if
        end do
      end if
c
      call errcon ('Unrecognised keyword')
      call textut (' >',line)
      goto 10
c
c
c
 1000 continue
      write (*,*)
      do i=1,maxi
        line = ' '//key(i)//' >'
        call textut (line,value(i))
      end do
      do i=1,ncopy
        numncs (i) = 0
        write (*,*)
        do j=1,maxj
          line = '   '//ckey(j)//' >'
          call textut (line,cvalue(i,j))
        end do
        do j=1,maxncs
          if (length(ncsfil(j,i)) .gt. 0) then
            call textut ('   NCSO >',ncsfil(j,i))
            numncs (i) = numncs (i) + 1
          end if
        end do
      end do
c
      write (*,*)
      call jvalut (' Nr of crystal forms    :',1,ncopy)
      if (ncopy .lt. 2) call errstp ('Fewer than 2 xtal forms')
c
      read (value(icycle),*) mcycle
      call jvalut (' Nr of averaging cycles :',1,mcycle)
      if (mcycle .lt. 1) call errstp ('No cycles to be executed')
c
c ... write the C-shell script
c
      call xopxua (iunit,value(ifile),xinter(),ierr)
      if (ierr .ne. 0) call errstp ('While opening file')      
c
 6000 format (20(a,1x))
c
      write (iunit,6000) '#!/bin/csh -f'
      call stamp (line)
      write (iunit,6000) '#',line(1:leng1(line))
      write (iunit,6000)
      write (iunit,6000) 'cd',value(iwork)(1:leng1(value(iwork)))
      write (iunit,6000)
c
      ncycle = 0
      write (*,*)
      call message (iunit,'CYCLE',value(iden),ncycle,line)
c
c ... MAVE 0
c
      call chkfil (iunit,value(imask))
      do i=1,ncopy
        call chkfil (iunit,cvalue(i,jmap))
        call chkfil (iunit,cvalue(i,jxxrt))
        call chkfil (iunit,cvalue(i,josym))
        do j=1,maxncs
          if (length(ncsfil(j,i)) .gt. 0) then
            call chkfil (iunit,ncsfil(j,i))
          end if
        end do
        call message (iunit,'average',cvalue(i,jcopy),ncycle,line)
        call getlogf (cvalue(i,jcopy),'average',ncycle,line)
        write (iunit,6000) value(imave)(1:leng1(value(imave))),
     +    '-b << EOF >&',line(1:leng1(line))
        write (iunit,6000) 'average'
        write (iunit,6000) cvalue(i,jmap)(1:leng1(cvalue(i,jmap)))
        write (iunit,6000) value(imask)(1:leng1(value(imask)))
        write (iunit,6000) cvalue(i,jxxrt)(1:leng1(cvalue(i,jxxrt)))
        write (iunit,6000) cvalue(i,josym)(1:leng1(cvalue(i,josym)))
        do j=1,maxncs
          if (length(ncsfil(j,i)) .gt. 0) then
            write (iunit,6000) ncsfil(j,i)(1:leng1(ncsfil(j,i)))
          end if
        end do
        write (iunit,6000)
        if (i .eq. 1) then
          line = cvalue(1,jmap)
        else
          call getmap (1,ncycle,1,value(iscra),cvalue(1,jcopy),line)
        end if
        write (iunit,6000) line(1:leng1(line))
        call getmap (1,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,6000) 'EOF'
        write (iunit,6000)
        call chkfil (iunit,line)
        call getlogf (cvalue(i,jcopy),'average',ncycle,line)
        call grep (iunit,'error',line)
        call grep (iunit,'Corr. coeff. for operator',line)
        call grep (iunit,'R-factor for operator',line)
      end do
c
c ... COMAP 0
c
      call message (iunit,'comap',value(iden),ncycle,line)
      call getlogf (value(iden),'comap',ncycle,line)
      write (iunit,6000) value(icomap)(1:leng1(value(icomap))),
     +  '-b << EOF >&',line(1:leng1(line))
      do i=1,ncopy
        call getmap (1,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,*) float(numncs(i))
      end do
      write (iunit,6000)
      call getmap (0,ncycle,i,value(iscra),value(iden),line)
      write (iunit,6000) line(1:leng1(line))
      write (iunit,6000) 'EOF'
      write (iunit,6000)
      call chkfil (iunit,line)
      call getlogf (value(iden),'comap',ncycle,line)
      call grep (iunit,'error',line)
      call grep (iunit,'Corr coeff       :',line)
c
c ... loop
c
      do ncycle=1,mcycle
c
        write (*,*)
        call message (iunit,'CYCLE',value(iden),ncycle,line)
c
c ... MAVE (expand) NCYCLE-1
c
      do i=1,ncopy
        call message (iunit,'expand',cvalue(i,jcopy),ncycle,line)
        call getlogf (cvalue(i,jcopy),'expand',ncycle,line)
        write (iunit,6000) value(imave)(1:leng1(value(imave))),
     +    '-b << EOF >&',line(1:leng1(line))
        write (iunit,6000) 'expand'
        call getmap (0,ncycle-1,i,value(iscra),value(iden),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,6000) value(imask)(1:leng1(value(imask)))
        write (iunit,6000) cvalue(i,jxxrt)(1:leng1(cvalue(i,jxxrt)))
        write (iunit,6000) cvalue(1,josym)(1:leng1(cvalue(1,josym)))
        write (iunit,6000) cvalue(i,josym)(1:leng1(cvalue(i,josym)))
        do j=1,maxncs
          if (length(ncsfil(j,i)) .gt. 0) then
            write (iunit,6000) ncsfil(j,i)(1:leng1(ncsfil(j,i)))
          end if
        end do
        write (iunit,6000)
        write (iunit,6000) cvalue(i,jmap)(1:leng1(cvalue(i,jmap)))
        call getmap (0,ncycle-1,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,6000) 'EOF'
        write (iunit,6000)
        call chkfil (iunit,line)
        call getlogf (cvalue(i,jcopy),'expand',ncycle,line)
        call grep (iunit,'error',line)
      end do
c
c ... SFALL/RSTATS/FFT/MAVE (average)
c
      do i=1,ncopy
c
c ... sfall
c
        call message (iunit,'sfall',cvalue(i,jcopy),ncycle,line)
        write (iunit,6000)
     +    value(isfall)(1:leng1(value(isfall))),' HKLIN',
     +    cvalue(i,jmtz)(1:leng1(cvalue(i,jmtz))),'\\'
        call getmap (0,ncycle-1,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000)
     +    '  MAPIN',line(1:leng1(line)),'\\'
        call getlogf (cvalue(i,jcopy),'sfall',ncycle,line)
        write (iunit,6000)
     +    '  HKLOUT temp1.mtz << EOF >',line(1:leng1(line))
        write (iunit,6000) 'TITLE calculate Fc and PHIc'
        write (iunit,6000) 'MODE SFCALC MAPIN HKLIN'
        write (iunit,6000) 'RESOLUTION',
     +    cvalue(i,jreso)(1:leng1(cvalue(i,jreso)))
        write (iunit,6000) 'SFSGRP 1'
        write (iunit,6000) 'BINS 40'
        write (iunit,6000) 'BADD 0'
        write (iunit,6000) 'NGAU 2'
        write (iunit,6000) 'FORM C H N O S P'
        write (iunit,6000) 'RSCB',
     +    cvalue(i,jreso)(1:leng1(cvalue(i,jreso)))
        write (iunit,6000) 'GRID',
     +    cvalue(i,jgrid)(1:leng1(cvalue(i,jgrid)))
        write (iunit,6000) 'LABIN',
     +    'FP='//cvalue(i,jlabf)(1:leng1(cvalue(i,jlabf))),
     +    'SIGFP='//cvalue(i,jsigf)(1:leng1(cvalue(i,jsigf)))
        write (iunit,6000) 'LABOUT FC=FCALC PHIC=PHICALC'
        write (iunit,6000) 'END'
        write (iunit,6000) 'EOF'
        write (iunit,6000)
        call chkfil (iunit,'temp1.mtz')
        call getlogf (cvalue(i,jcopy),'sfall',ncycle,line)
        call grep (iunit,'error',line)
        call grep (iunit,'Overall Reliability index is',line)
c
c ... rstats
c
        call message (iunit,'rstats',cvalue(i,jcopy),ncycle,line)
        write (iunit,6000)
     +    value(irstats)(1:leng1(value(irstats))),
     +    ' HKLIN temp1.mtz HKLOUT temp2.mtz \\'
        call getlogf (cvalue(i,jcopy),'rstats',ncycle,line)
        write (iunit,6000)
     +    '  << EOF >',line(1:leng1(line))
        write (iunit,6000) 'TITLE scale Fc and Fo'
        write (iunit,6000) 'PROCESS FCAL'
        write (iunit,6000) 'SCALE 1.0'
        write (iunit,6000) 'OUTPUT ASIN'
        write (iunit,6000) 'PRINT ALL'
        write (iunit,6000) 'LABIN',
     +    'FP='//cvalue(i,jlabf)(1:leng1(cvalue(i,jlabf))),
     +    'SIGFP='//cvalue(i,jsigf)(1:leng1(cvalue(i,jsigf))),
     +    'FC=FCALC PHIC=PHICALC'
        write (iunit,6000) 'RESOLUTION',
     +    cvalue(i,jreso)(1:leng1(cvalue(i,jreso)))
        write (iunit,6000) 'WIDTH_OF_BINS RTHETA=0.01 FBINR=500'
        write (iunit,6000) 'CYCLES 5'
        write (iunit,6000) 'LIST 1000000'
        write (iunit,6000) 'WEIGHTING_SCHEME NONE'
        write (iunit,6000) 'EOF'
        write (iunit,6000)
        call chkfil (iunit,'temp1.mtz')
        call getlogf (cvalue(i,jcopy),'rstats',ncycle,line)
        call grep (iunit,'error',line)
        call grep (iunit,'Range 4',line)
        call grep (iunit,'Overall Totals:',line)
c
c ... fft
c
        call message (iunit,'fft',cvalue(i,jcopy),ncycle,line)
        write (iunit,6000)
     +    value(ifft)(1:leng1(value(ifft))),
     +    ' HKLIN temp2.mtz \\'
        call getmap (0,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000)
     +    '  MAPOUT',line(1:leng1(line)),'\\'
        call getlogf (cvalue(i,jcopy),'fft',ncycle,line)
        write (iunit,6000)
     +    '  << EOF >',line(1:leng1(line))
        write (iunit,6000) 'TITLE calculate new map'
        write (iunit,6000) 'RESOLUTION',
     +    cvalue(i,jreso)(1:leng1(cvalue(i,jreso)))
        write (iunit,6000) 'SCALE F1 2.0 0.0'
        write (iunit,6000) 'SCALE F2 1.0 0.0'
        write (iunit,6000) 'FFTSYMMETRY 1'
        write (iunit,6000) 'GRID',
     +    cvalue(i,jgrid)(1:leng1(cvalue(i,jgrid)))
        write (iunit,6000) 'XYZLIMIT',
     +    cvalue(i,jexte)(1:leng1(cvalue(i,jexte)))
        write (iunit,6000) 'LABIN',
     +    'F1='//cvalue(i,jlabf)(1:leng1(cvalue(i,jlabf))),
     +    'SIG1='//cvalue(i,jsigf)(1:leng1(cvalue(i,jsigf))),
     +    'F2=FCALC',
     +    'SIG2='//cvalue(i,jsigf)(1:leng1(cvalue(i,jsigf))),
     +    'PHI=PHICALC'
        write (iunit,6000) 'RHOLIM 100.0'
        write (iunit,6000) 'EOF'
        write (iunit,6000)
        call getmap (0,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        call chkfil (iunit,line)
        call getlogf (cvalue(i,jcopy),'fft',ncycle,line)
        call grep (iunit,'error',line)
        call grep (iunit,'Rms deviation from mean density',line)
c
c ... mave (average)
c
        call message (iunit,'average',cvalue(i,jcopy),ncycle,line)
        call getlogf (cvalue(i,jcopy),'average',ncycle,line)
        write (iunit,6000) value(imave)(1:leng1(value(imave))),
     +    '-b << EOF >&',line(1:leng1(line))
        write (iunit,6000) 'average'
        call getmap (0,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,6000) value(imask)(1:leng1(value(imask)))
        write (iunit,6000) cvalue(i,jxxrt)(1:leng1(cvalue(i,jxxrt)))
        write (iunit,6000) cvalue(i,josym)(1:leng1(cvalue(i,josym)))
        do j=1,maxncs
          if (length(ncsfil(j,i)) .gt. 0) then
            write (iunit,6000) ncsfil(j,i)(1:leng1(ncsfil(j,i)))
          end if
        end do
        write (iunit,6000)
        write (iunit,6000) cvalue(1,jmap)(1:leng1(cvalue(1,jmap)))
        call getmap (1,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,6000) 'EOF'
        write (iunit,6000)
        call chkfil (iunit,line)
        write (iunit,6000) '\\rm temp1.mtz temp2.mtz'
        call getlogf (cvalue(i,jcopy),'average',ncycle,line)
        call grep (iunit,'error',line)
        call grep (iunit,'Corr. coeff. for operator',line)
        call grep (iunit,'R-factor for operator',line)
c
      end do
c
c ... comap
c
      call message (iunit,'comap',value(iden),ncycle,line)
      call getlogf (value(iden),'comap',ncycle,line)
      write (iunit,6000) value(icomap)(1:leng1(value(icomap))),
     +  '-b << EOF >&',line(1:leng1(line))
      do i=1,ncopy
        call getmap (1,ncycle,i,value(iscra),cvalue(i,jcopy),line)
        write (iunit,6000) line(1:leng1(line))
        write (iunit,*) float(numncs(i))
      end do
      write (iunit,6000)
      call getmap (0,ncycle,i,value(iscra),value(iden),line)
      write (iunit,6000) line(1:leng1(line))
      write (iunit,6000) 'EOF'
      write (iunit,6000)
      call chkfil (iunit,line)
      call getlogf (value(iden),'comap',ncycle,line)
      call grep (iunit,'error',line)
      call grep (iunit,'Corr coeff       :',line)
c
c ... end do NCYCLE
c
      end do
c
c ... MAPMAN
c
      write (*,*)
      call message (iunit,'mapman',value(iden),mcycle,line)
      call getlogf (value(iden),'mapman',mcycle,line)
      write (iunit,6000) value(imapman)(1:leng1(value(imapman))),
     +  '-b << EOF >&',line(1:leng1(line))
      call getmap (0,mcycle,i,value(iscra),value(iden),line)
      write (iunit,6000) 'read m1',line(1:leng1(line)),'ccp4'
      write (iunit,6000) 'mappage m1 final.omap'
      write (iunit,6000) 'EOF'
      write (iunit,6000)
      call chkfil (iunit,'final.omap')
c
      call message (iunit,'all done',value(iden),mcycle,line)
      write (iunit,6000) 'exit 0'
c
      call gkquit ()
c
      end
c
c
c
      subroutine getmap (mode,icycle,ixtal,scra,id,line)
c
      implicit none
c
      integer mode,icycle,ixtal,ll,length
c
      character*(*) scra,id,line
c
code ...
c
      line = scra
      ll = length(line)
      if (line(ll:ll) .ne. '/') then
        line(ll+1:ll+1) = '/'
        ll = ll + 1
      end if
c
      call appstr (line,id)
      ll = length(line)
      write (line(ll+1:),'(a1,i6)') '_',icycle
      ll = length(line)
      if (mode .eq. 1) then
        line(ll+1:ll+1) = 'x'
        ll = ll + 1
      end if
      line (ll+1:) = '.E'
c
      call remspa (line)
c
      return
      end
c
c
c
      subroutine getlogf (id,prog,icycle,line)
c
      implicit none
c
      integer icycle,ll,length
c
      character*(*) id,prog,line
c
code ...
c
      line = id
      call appstr (line,'_')
      call appstr (line,prog)
      call appstr (line,'_')
      ll = length(line)
      write (line(ll+1:),'(i6)') icycle
      call appstr (line,'.log')
c
      call remspa (line)
      call locase (line)
c
      return
      end
c
c
c
      subroutine message (iunit,task,id,icycle,line)
c
      implicit none
c
      integer iunit,icycle,ll,length,leng1
c
      character*(*) task,id,line
c
code ...
c
      line = ' '//task
      ll = length(line)
      line (ll+2:) = id
      ll = length(line)
      write (line (ll+2:),*) icycle
      call pretty (line)
      write (iunit,'(a5)') 'echo '
      write (iunit,'(a5,a)') 'echo ',line(1:leng1(line))
      write (iunit,'(a5)') 'echo '
      call textut (' TASK >',line)
c
      return
      end
c
c
c
      subroutine chkfil (iunit,file)
c
      implicit none
c
      integer iunit,leng1
c
      character*(*) file
c
code ...
c
      write (iunit,6000) 'if (! -e',file(1:leng1(file)),') then'
      write (iunit,6000) '  echo ERROR ... file',file(1:leng1(file)),
     +  'not found ... aborting'
      write (iunit,6000) '  exit -1'
      write (iunit,6000) 'endif'
c
 6000 format (20(a,1x))
c
      return
      end
c
c
c
      subroutine grep (iunit,text,file)
c
      implicit none
c
      integer iunit,leng1
c
      character*(*) file,text
c
code ...
c
      write (iunit,'(4a)') 'grep -i "',text(1:leng1(text)),'" ',
     +  file(1:leng1(file))
c
      return
      end
