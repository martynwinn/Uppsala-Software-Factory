      program dataman
c
c ... DATAMAN - manipulation of ASCII HKL-files
c
c ... Gerard Kleywegt @ 930319
c
c ... 0.1 @ 930319 - first version
c
      implicit none
c
      include 'dataman_dim.incl'
c
      character*12 prognm,vers
      parameter (prognm = 'DATAMAN', vers = '061208/6.4.2')
c
c      pointer (iaptr,fobs)
c      pointer (ibptr,sigfob)
c      pointer (icptr,reso)
c      pointer (idptr,hkl)
c      pointer (ieptr,morbit)
c      pointer (ifptr,rfree)
c      pointer (igptr,centri)
c      pointer (ihptr,buffer)
c
c      real fobs(1),sigfob(1),reso(1),buffer(1)
c
c      integer hkl(1),morbit(1),rfree(1),malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr,idptr,ieptr,ifptr,igptr,ihptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr,idptr,ieptr,ifptr,igptr,ihptr
      integer fmalloc
#endif
c
      integer nb,numhkl,numset,numbuf,minsiz,i1,i2
c
      logical lretry,ldone
c
c      character centri(1)*1
c
code ...
c
      call gainit (prognm,vers)
c
c ... initialise history
c
      call dohist ('*INIT*',ldone)
c
      numhkl = 200000
      minsiz = 10000
      numset = 4
c
      call extint ('SETSIZE',numhkl)
      numhkl = max ( numhkl , minsiz )
      call jvalut (' Allocate data sets of size :',1,numhkl)
c
      call extint ('NUMSETS',numset)
      numset = max ( min ( numset, maxgk9 ), 1 )
      call jvalut (' Max number of data sets    :',1,numset)
      write (*,*)
c
c ... WRDBYT accounts for 4 or 8 bytes per word
c
   10 continue
      nb = wrdbyt*numhkl*numset
      iaptr = fmalloc (nb)
      ibptr = fmalloc (nb)
      icptr = fmalloc (nb)
      ieptr = fmalloc (nb)
      ifptr = fmalloc (nb)
c
      nb = 3*wrdbyt*numhkl*numset
      idptr = fmalloc (nb)
c
c ... character*1 array
      nb = 1*numhkl*numset
      igptr = fmalloc (nb)
c
      numbuf = 5 * numhkl
      nb = wrdbyt*numbuf
      ihptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0 .or.
     +    idptr .eq. 0 .or. ieptr .eq. 0 .or. ifptr .eq. 0 .or.
     +    igptr .eq. 0 .or. ihptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0 .or.
     +    idptr .le. 0 .or. ieptr .le. 0 .or. ifptr .le. 0 .or.
     +    igptr .le. 0 .or. ihptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
c ... BUFFER and IBUFF were equivalenced in the old version of MAMA
c     without dynamic memory allocation (they should never both be used
c     at the same time); passing BUFFER for both has the same net effect
c     (i.e., using the same memory addresses for two different arrays)
c
      lretry = .false.
      i1 = numhkl
      i2 = numset
c
      call dodata (numset,numhkl,numbuf,
     +             %val(iaptr),%val(ibptr),%val(icptr),
     +             %val(idptr),%val(ieptr),%val(ifptr),
     +             %val(igptr),%val(ihptr),%val(ihptr),
     +             lretry,i1,i2)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
      call ffree (idptr)
      call ffree (ieptr)
      call ffree (ifptr)
      call ffree (igptr)
      call ffree (ihptr)
c
      if (lretry) then
        numhkl = i1
        numhkl = max ( numhkl , minsiz )
        numset = i2
        numset = max ( min ( numset, maxgk9 ), 1 )
        write (*,*)
        call jvalut (' Allocate data sets of size :',1,numhkl)
        call jvalut (' Max number of data sets    :',1,numset)
        write (*,*)
        goto 10
      end if
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine dodata (maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,
     +                   buffer,ibuff,lretry,ii1,ii2)
c
c ... DATAMAN - manipulation of ASCII HKL-files
c
c ... Gerard Kleywegt @ 930319
c
c ... 0.1 @ 930319 - first version
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer morbit(maxhkl,maxset),rfree(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      integer maxopt,maxhis,maxlau
c
      parameter (maxopt = 17)
      parameter (maxhis = maxopt - 3)
c
      parameter (maxlau = 15)
c
      real total,user,sys,xdum,xdum2,dx,rhis(maxhis)
      real rfperc,step,dpbin,resol,mw,cvol,x,reshi,vm
      real xminoi,xmanoi,sigmin,sigmax
c
      integer numpar,i,j,k,length,nempty,nopt,ierr,iset,iunit
      integer whichm,nhis(maxhis+1),numhis,jset,kset,laue,i123(3)
      integer i1,i2,i3,mathkl(3,3),inow,iseed,nfneg,nsneg,nasu
      integer nr,nres,nncs,nbins,idum,jdum,kdum,ldum,nstart
      integer munit,leng1,ii1,ii2
c
      logical xinter,linter,unsave
      logical ldone,linit,lretry,lecho
c
      character line*256,optpar(maxopt)*128,reply*1,lattic*1
      character parnam*40,partyp*1,parfmt*40,pro*12,prev*10
      character inimac*128
      character*75 lautxt(maxlau)
c
code ...
c
      lretry = .false.
      lecho = .false.
c
ccc      call gkinit (prognm,vers)
c
      call jvalut (' Max nr of data sets           :',1,maxset)
      call jvalut (' Max nr of reflections per set :',1,maxhkl)
      call jvalut (' Max nr of symmetry operators  :',1,maxsym)
ccc      call jvalut (' Max BYTES for reflection data :',1,totbyt)
      write (*,*)
c
      call gkrand (dx,0.0,0.0,-1)
      iseed = 620605
      rfperc = 10.0
      step = 0.0025
      dpbin = -15.0
      nbins = 15
      resol = 2.0
      reshi = 100.0
      nasu = 4
      lattic = 'P'
      nres = 100
      mw = 11200.0
      nncs = 1
c
c ... user input unit (5=interactive; other=macro)
c
      munit = 5
c
      linit = .false.
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
      nempty = 0
      iunit = 1
c
      linter = xinter()
      if (linter) then
        pro='$DATAMAN > '
      else
        pro=' DATAMAN > '
      end if
c
      do i=1,maxset
        numhkl (i) = 0
        nsymop (i) = 0
        name   (i) = '$#@%^&@^62'
        file   (i) = 'not_saved_yet'
        coment (i) = 'No comment'
        incore (i) = .false.
        change (i) = .false.
        do j=1,maxitm
          know (j,i) = .false.
        end do
        do j=1,3
          cell (j,i)   = 100.0
          cell (j+3,i) = 90.0
        end do
      end do
c
      lautxt(1) =
     +   'LAUE = 1, 1bar,   hkl:h>=0   0kl:k>=0   00l:l>=0'
      lautxt(2) =
     +   'LAUE = 2, 1bar,   hkl:k>=0   h0l:l>=0   h00:h>=0'
      lautxt(3) =
     +   'LAUE = 3, 1bar,   hkl:l>=0   hk0:h>=0   0k0:k>=0'
cxyz changed 940711 (added k>=0 to hk0)
      lautxt(4) =
     +   'LAUE = 4, 2/m,    hkl:k>=0, l>=0     hk0:h>=0, k>=0'
cxyz changed 940711 (added l>=0 to 0kl)
      lautxt(5) =
     +   'LAUE = 5, 2/m,    hkl:h>=0, l>=0     0kl:k>=0, l>=0'//
     +   '   (2-nd sett.)'
      lautxt(6) =
     +   'LAUE = 6, mmm,    hkl:h>=0, k>=0, l>=0'
      lautxt(7) =
     +   'LAUE = 7, 4/m,    hkl:h>=0, k>0, l>=0 with  k>=0 for h=0'
      lautxt(8) =
     +   'LAUE = 8, 4/mmm,  hkl:h>=0, h>=k>=0, l>=0'
      lautxt(9) =
     +   'LAUE = 9, 3bar,   hkl:h>=0, k<0, l>=0 including 00l'
      lautxt(10) =
     +   'LAUE = 10, 3bar,  hkl:h>=0, k>0  including  00l:l>0'
      lautxt(11) =
     +   'LAUE = 11, 3barm, hkl:h>=0, k>=0 with k<=h; if h=k l>=0'
      lautxt(12) =
     +   'LAUE = 12, 6/m,   hkl:h>=0, k>0, l>=0  with  k>=0 for h=0'
      lautxt(13) =
     +   'LAUE = 13, 6/mmm, hkl:h>=0, h>=k>=0, l>=0'
      lautxt(14) =
     +   'LAUE = 14, m3,    hkl:h>=0, k>=0, l>=0 with l>=h, k>=h'//
     +   ' for l=h, k>h if l>h'
      lautxt(15) =
     +   'LAUE = 15, m3m,   hkl:k>=l>=h>=0'
c
c ... formats
c
 6000 format (/
     +  ' DATAMAN options :'//
     +  ' ? (list options)                    ',
     +              ' ! (comment)'/
     +  ' QUit                                ',
     +              ' $ shell_command'/
     +  ' & symbol value                      ',
     +              ' & ? (list symbols)'/
     +  ' @ macro_file                        ',
     +              ' ZP_restart setsize numsets'/
     +  ' ECho on_off                         ',
     +              ' # parameter(s) (command history)'/
     +  ' DElete set                          ',
     +              ' RAmp_odl set file ramp_option'/
     +  ' REad_refl set file type [format]    ',
     +              ' APpend_refl set file type [format]'/
     +  ' WRite_ref set file type [format] [which_ATW] [which_ABC]'/
     + /' LIst set                            ',
     +              ' STats set'/
     +  ' HIsto set which x1 x2 x3 [...]      ',
     +              ' SHow_hkl set criterion operand value'/
     +  ' CEll set a b c al be ga             ',
     +              ' ANnotate set "text"'/
     +  ' SYmmop set o_file                   ',
     +              ' TYpe_hkl set start end step'/
     +  ' SPecial set hkl_type                ',
     +              ' RSym_hkl_khl set'/
     +  ' ABsences set [list_or_kill]         ',
     +              ' RInt set'/
     +  ' PArity_test set                     ',
     +              ' FIll_in set nbins'/
     + /' TWin_stats set                      ',
     +              ' GEmini set plotf1 plotf2'/
     +  ' LAue newset set laue_group          ',
     +              ' SOrt_hkl newset set hkl_order'/
     +  ' KIll_hkl set criterion operand value',
     +              ' PRod_plus set which prod plus'/
     +  ' ODd_kill set h_k_l                  ',
     +              ' EVen_kill set h_k_l'/
     +  ' CAlc set what                       ',
     +              ' TEmp_factor set value'/
     +  ' CHange_index set newh newk newl     ',
     +              ' ROgue_kill set h1 k1 l1 [...]'/
     +  ' NOise set nbins min% max%           ',
     +              ' DUplicate newset set'/
     +  ' HEmisphere newset set resolution    ',
     +              ' ASym_unit newset set resol laue_group'/
     +  ' SIgmas FAke set                     ',
     +              ' SIgmas LImit set lower upper'/
     +  ' SIgmas CEntric_vs_acentric set      ',
     +              ' PY_stats set plotfile'/
     +  ' MUltiplicity set                    ',
     +              ' '/
     + /' WIlson set1 set2 plotf1 plotf2 step ',
     +              ' DF newset set1 set2'/
     +  ' COmpare set1 set2                   ',
     +              ' MErge newset set1 set2 how'/
     +  ' YEates_stats set1 set2              ',
     +              ' '/
     + /' RFree INit seed                     ',
     +              ' RFree LIst set'/
     +  ' RFree GEnerate set %_or_#           ',
     +              ' RFree REset set'/
     +  ' RFree SHell set %_or_# nbins        ',
     +              ' RFree COmplete set nsets basename'/
     +  ' RFree GSheldrick set nth            ',
     +              ' RFree SPheres set %_or_# radius'/
     +  ' RFree TRansfer set old_set          ',
     +              ' RFree ADjust set new%'/
     +  ' RFree FIll_bins set target% nbins   ',
     +              ' RFree CUt_bins set target% nbins'/
     +  ' RFree BIn_list set nbins            ',
     +              ' RFree MUlti set nsets how'/
     +  ' RFree SUggest set                   ',
     +              ' '/
     + /' EStimate_unique set resol latt nasu ',
     +              ' EFfective_resolution set latt nasu'/
     +  ' GUess MW nres                       ',
     +              ' GUess NRes MW'/
     +  ' GUess VM set nres nasu nncs         ',
     +              ' GUess COmpl set res1 res2 latt nasu'/
     +  ' GUess RHo set latt nasu nncs nres  ',' '/
     + /' SCatter_plot set file hori vert     ',
     +              ' BIn_plot set file hori vert bin'/
     +  ' HKl_aniso_plot set file             ',
     +              ' DOuble_plot set1 set2 file hor ver bin'/
     +  ' EO_plot set file hkl                ',
     +              ' '/
     +  )
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6000)
c
      call jvalut (' Max nr of data sets           :',1,maxset)
      call jvalut (' Max nr of reflections per set :',1,maxhkl)
      call jvalut (' Max nr of symmetry operators  :',1,maxsym)
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
        call gknval ('GKDATAMAN',inimac,ierr)
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
c ... ZP_RESTART
c
      else if (optpar(1)(1:2) .eq. 'ZP') then
c
        if (nopt .lt. 2) then
          write (optpar(2),*) ii1
          call textin (' New SETSIZE ?',optpar(2))
        end if
        call str2i (optpar(2),idum,ierr)
        if (ierr .ne. 0) goto 10
        ii1 = idum
        if (ii1 .lt. 10) then
          call errcon ('Silly SETSIZE')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          write (optpar(3),*) ii2
          call textin (' New NUMSETS ?',optpar(3))
        end if
        call str2i (optpar(3),idum,ierr)
        if (ierr .ne. 0) goto 10
        ii2 = idum
        if (ii2 .lt. 1) then
          call errcon ('Silly NUMSETS')
          goto 10
        end if
c
        lretry = .true.
        return
c
c ... RAMP_ODL
c
      else if (optpar(1)(1:2) .eq. 'RA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'reciprocal_lattice.odl'
          call textin (' ODL file name ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'RES'
          call textin (
     +      ' Colour ramping criterion [FOB|SIG|F/S|RES|NONe] ?',
     +      optpar(4))
        end if
c
        call reclat (iset,iunit,optpar(3),optpar(4),ierr,
     +               maxset,maxhkl,maxbuf,
     +               fobs,sigfob,reso,hkl,buffer)
c
c ... READ
c
      else if (optpar(1)(1:2) .eq. 'RE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),iset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(iset)
c
        if (nopt .lt. 3) then
          optpar (3) = file(iset)
          call textin (' File name ?',optpar(3))
        end if
        file (iset) = optpar(3)
c
        if (nopt .lt. 4) then
          optpar (4) = '*'
          call textin (' File type ?',optpar(4))
        end if
c
        if (nopt .lt. 5) optpar (5) = '*'
c
        nstart = 0
c
        call datain (iunit,file(iset),optpar(4),optpar(5),maxhkl,
     +    nstart,hkl(1,1,iset),numhkl(iset),fobs(1,iset),
     +    sigfob(1,iset),rfree(1,iset),reso(1,iset),morbit(1,iset),
     +    centri(1,iset),know(kcell,iset),cell(1,iset),ierr)
c
        if (ierr .ne. 0) then
          name (iset) = '!@#$%^&*()'
          goto 10
        end if
c
        if (numhkl(iset) .le. 0) then
          call errcon ('No reflections read; deleting set')
          incore (iset) = .false.
          goto 10
        end if
c
        change (iset) = .false.
        incore (iset) = .true.
        coment (iset) = 'Read from '//file(iset)
        know (kdata,iset) = .true.
c
        call testrf (iset,maxset,maxhkl,rfree)
c
        goto 10
c
c ... APPEND
c
      else if (optpar(1)(1:2) .eq. 'AP') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = file(iset)
          call textin (' File name ?',optpar(3))
        end if
        file (iset) = optpar(3)
c
        if (nopt .lt. 4) then
          optpar (4) = '*'
          call textin (' File type ?',optpar(4))
        end if
c
        if (nopt .lt. 5) optpar (5) = '*'
c
        nstart = numhkl(iset)
c
        call datain (iunit,file(iset),optpar(4),optpar(5),maxhkl,
     +    nstart,hkl(1,1,iset),numhkl(iset),fobs(1,iset),
     +    sigfob(1,iset),rfree(1,iset),reso(1,iset),morbit(1,iset),
     +    centri(1,iset),know(kcell,iset),cell(1,iset),ierr)
c
c        if (ierr .ne. 0) then
c          name (iset) = '!@#$%^&*()'
c          goto 10
c        end if
c
c        if (numhkl(iset) .le. 0) then
c          call errcon ('No reflections read; deleting set')
c          incore (iset) = .false.
c          goto 10
c        end if
c
        change (iset) = .true.
        incore (iset) = .true.
        coment (iset) = 'Appended from '//file(iset)
        know (kdata,iset) = .true.
c
        call testrf (iset,maxset,maxhkl,rfree)
c
        goto 10
c
c ... WRITE
c
      else if (optpar(1)(1:2) .eq. 'WR') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = file(iset)
          call textin (' File name ?',optpar(3))
        end if
        file (iset) = optpar(3)
c
        if (nopt .lt. 4) then
          optpar (4) = '*'
          call textin (' File type ?',optpar(4))
        end if
c
        if (nopt .lt. 5) optpar (5) = '*'
c
        if (nopt .lt. 6) optpar (6) = 'ALL'
        call upcase (optpar(6))
c
        if (nopt .lt. 7) optpar (7) = 'BOTH'
        call upcase (optpar(7))
        if (optpar(7)(1:1) .ne. 'B' .and.
     +      (.not. know (kcent,iset))) then
          call errcon (' (A)centrics have NOT been deduced')
          goto 10
        end if
c
        call testrf (iset,maxset,maxhkl,rfree)
c
        call dataut (iunit,file(iset),optpar(4),optpar(5),maxhkl,
     +    hkl(1,1,iset),numhkl(iset),fobs(1,iset),sigfob(1,iset),
     +    rfree(1,iset),optpar(6),know(kcell,iset),cell(1,iset),
     +    optpar(7),centri(1,iset),ierr)
c
        if (ierr .ne. 0) goto 10
c
        change (iset) = .false.
c
        goto 10
c
c ... COMPARE
c
      else if (optpar(1)(1:2) .eq. 'CO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set 1 ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Set 2 ?',optpar(3))
        end if
        jset = whichm (optpar(3),maxset)
c
        if (jset .le. 0 .or. jset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(jset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (iset .eq. jset) then
          call errcon ('Data sets must be different')
          goto 10
        end if
c
c        call compar (iset,jset,
c     +               maxset,maxhkl,maxbuf,
c     +               fobs,hkl,buffer,ibuff)
c
        call compar (iset,jset,
     +               maxset,maxhkl,maxbuf,
     +               fobs,hkl,buffer,ibuff)
c
c ... YEATES_STATS
c
      else if (optpar(1)(1:2) .eq. 'YE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set 1 ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Set 2 ?',optpar(3))
        end if
        jset = whichm (optpar(3),maxset)
c
        if (jset .le. 0 .or. jset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(jset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (iset .eq. jset) then
          call errcon ('Data sets must be different')
          goto 10
        end if
c
        if (.not.know(kcent,iset)) then
          if (.not.know(kcent,jset)) then
            call errcon (' (A)centrics have NOT been deduced')
            call prompt (
     +        ' Calculate centrics for at least one set')
            goto 10
          else
            call prompt (' (A)centrics only known for set 2')
            call prompt (' Swapping set 1 and set 2')
            call iswap (iset,jset)
          end if
        end if
c
        call yeates (iset,jset,
     +               maxset,maxhkl,maxbuf,
     +               fobs,hkl,centri,ibuff)
c
c ... DF (SHELXS deltaF)
c
      else if (optpar(1)(1:2) .eq. 'DF') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),kset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(kset)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Set 1 ?',optpar(3))
        end if
        iset = whichm (optpar(3),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (3)
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Set 2 ?',optpar(4))
        end if
        jset = whichm (optpar(4),maxset)
c
        if (jset .le. 0 .or. jset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(jset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (iset .eq. jset) then
          call errcon ('Data sets must be different')
          goto 10
        end if
c
        call deltaf (kset,iset,jset,ierr,
     +               maxset,maxhkl,maxbuf,
     +               fobs,sigfob,reso,hkl,rfree,
     +               morbit,centri,ibuff)
c
        if (ierr .ne. 0) goto 10
c
        change (kset) = .true.
        incore (kset) = .true.
        coment (kset) = 'DeltaF from '//name(iset)//' & '//
     +                  name(jset)
        call pretty (coment (kset))
        know (kdata,kset) = .true.
        call cpknow (iset,kset)
        call testrf (kset,maxset,maxhkl,rfree)
c
c ... MERGE
c
      else if (optpar(1)(1:2) .eq. 'ME') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),kset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(kset)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Set 1 ?',optpar(3))
        end if
        iset = whichm (optpar(3),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (3)
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Set 2 ?',optpar(4))
        end if
        jset = whichm (optpar(4),maxset)
c
        if (jset .le. 0 .or. jset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(jset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (iset .eq. jset) then
          call errcon ('Data sets must be different')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'SIG'
          call textin (' Method [?|SIG|AVE|COM] ?',optpar(5))
        end if
        call upcase (optpar(5))
c
        if (optpar(5)(1:1) .eq. '?') then
          write (*,'(a/a,a/a,a/a,a)')
     +      ' Select one of:',
     +      ' SIG = sigma weighting: Fnew = (S2*F1+S1*F2), ',
     +        'Snew = 2*S1*S2/(S1+S2)',
     +      ' AVE = average: Fnew = (F1+F2)/2, ',
     +        'Snew = 1/2*SQRT(S1^2+S2^2)',
     +      ' COM = complement: new set = set1 + all data ',
     +        'from set2 not in set1'
          goto 10
        end if
c
        if (index ('SIG|AVE|COM',optpar(5)(1:3)) .le. 0) then
          call errcon ('Invalid method')
          call textut (' Method :',optpar(5))
          goto 10
        end if
c
        call merger (kset,iset,jset,optpar(5)(1:3),ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,buffer,
     +                   morbit,centri,ibuff)
c
        if (ierr .ne. 0) goto 10
c
        change (kset) = .true.
        incore (kset) = .true.
        coment (kset) = 'Merged from '//name(iset)//' & '//
     +                  name(jset) // ' Method ' // optpar(5)
        call pretty (coment (kset))
        know (kdata,kset) = .true.
        call cpknow (iset,kset)
        call testrf (kset,maxset,maxhkl,rfree)
        call prompt (' The new dataset is UNSORTED !')
c
c ... WILSON
c
      else if (optpar(1)(1:2) .eq. 'WI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set 1 ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Set 2 ?',optpar(3))
        end if
        jset = whichm (optpar(3),maxset)
c
        if (jset .le. 0 .or. jset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(jset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (iset .eq. jset) then
          call errcon ('Data sets must be different')
          goto 10
        end if
c
        if (.not.know(kreso,iset)) then
          call errcon ('Resolution of set 1 unknown')
          goto 10
        end if
c
        if (.not.know(kreso,jset)) then
          call errcon ('Resolution of set 2 unknown')
          goto 10
        end if
c
        if (.not.know(korbi,iset)) then
          call errcon ('Orbital multiplicity of set 1 unknown')
          goto 10
        end if
c
        if (.not.know(korbi,jset)) then
          call errcon ('Orbital multiplicity of set 2 unknown')
          goto 10
        end if
c
        line = 'wilson_'//name(iset)//'_'//name(jset)//'_'
        call remspa (line)
        call locase (line)
c
        if (nopt .lt. 4) then
          optpar (4) = line(1:leng1(line))//'1.plt'
          call textin (' Name of first plot file ?',optpar(4))
        end if
        if (length(optpar(4)) .lt. 1) then
          optpar (4) = line(1:leng1(line))//'1.plt'
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = line(1:leng1(line))//'2.plt'
          call textin (' Name of second plot file ?',optpar(5))
        end if
        if (length(optpar(5)) .lt. 1) then
          optpar (5) = line(1:leng1(line))//'2.plt'
        end if
c
        if (nopt .lt. 6) then
          write (optpar(6),*) step
          call remspa (optpar(6))
          call textin (' Step size ?',optpar(6))
        end if
        call str2r (optpar(6),step,ierr)
        if (ierr .ne. 0) goto 10
c
        call gerard (iset,jset,optpar(4),optpar(5),step,ierr,
     +               maxset,maxhkl,maxbuf,
     +               fobs,reso,morbit,buffer)
        if (ierr .eq. 0) change (jset) = .true.
c
c ... GEMINI
c
      else if (optpar(1)(1:2) .eq. 'GE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (.not.know(kcent,iset)) then
          call errcon (' (A)centrics have NOT been deduced')
          goto 10
        end if
c
        line = 'gemini_'//name(iset)//'_'
        call remspa (line)
c
        if (nopt .lt. 3) then
          optpar (3) = line(1:leng1(line))//'1.ps'
          call textin (' Name of first PostScript file ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) then
          optpar (3) = line(1:leng1(line))//'1.ps'
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = line(1:leng1(line))//'2.ps'
          call textin (' Name of second PostScript file ?',optpar(4))
        end if
        if (length(optpar(4)) .lt. 1) then
          optpar (4) = line(1:leng1(line))//'2.ps'
        end if
c
        call gemini (iset,optpar(3),optpar(4),
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,centri,buffer)
c
c ... PY_STATS
c
      else if (optpar(1)(1:2) .eq. 'PY') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (.not.know(kcent,iset)) then
          call errcon (' (A)centrics have NOT been deduced')
          goto 10
        end if
c
        line = 'local_intensity_plot_'//name(iset)//'.ps'
        call remspa (line)
c
        if (nopt .lt. 3) then
          optpar (3) = line
          call textin (' Name of PostScript file ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) then
          optpar (3) = line
        end if
c
        call locint (iset,optpar(3),
     +               maxset,maxhkl,maxhkl,
     +               fobs,hkl,centri,
     +               ibuff(1),ibuff(maxhkl+1),
     +               ibuff(2*maxhkl+1),buffer(3*maxhkl+1))
c
c ... MUltiplicity
c
      else if (optpar(1)(1:2) .eq. 'MU') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Multiplicity :',name(iset))
            call multip (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   hkl,ibuff)
          end if
        end do
c
c ... RSym_hkl_khl
c
      else if (optpar(1)(1:2) .eq. 'RS') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Rsym (hkl,khl) :',name(iset))
            call hklkhl (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,buffer,ibuff)
          end if
        end do
c
c ... RInt
c
      else if (optpar(1)(1:2) .eq. 'RI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Rint :',name(iset))
            call rint (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,buffer,ibuff)
          end if
        end do
c
c ... TWin_stats
c
      else if (optpar(1)(1:2) .eq. 'TW') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            if (know(kcent,iset)) then
              write (*,*)
              call textut (' Twin_stats :',name(iset))
              call twin (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,centri,buffer)
            else
              call errcon ('(A)centrics have NOT been deduced')
            end if
          end if
        end do
c
c ... ABsences
c
      else if (optpar(1)(1:2) .eq. 'AB') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'L'
        end if
        call upcase (optpar(3))
        if (optpar(3)(1:1) .ne. 'K') then
          optpar (3) = 'List'
        else
          optpar (3) = 'Kill'
        end if
        line = ' ' // optpar(3)(1:4) //
     +         ' systematic absences for :'
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (line,name(iset))
            if (know(ksymm,iset)) then
              write (*,*)
              call sysabs (iset,optpar(3),
     +                   maxset,maxhkl,maxbuf,fobs,sigfob,hkl,
     +                   rfree,reso,morbit,centri,ibuff)
            else
              call errcon ('I don''t know the symmops')
            end if
          end if
        end do
c
c ... PArity_test
c
      else if (optpar(1)(1:2) .eq. 'PA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Parity test for :',name(iset))
            call parity (iset,maxset,maxhkl,fobs,hkl)
          end if
        end do
c
c ... FIll_in
c
      else if (optpar(1)(1:2) .eq. 'FI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (.not.know(kreso,iset)) then
          call errcon (' Resolution has not been calculated')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          write (optpar(3),*) nbins
          call textin (' Nr of bins ?',optpar(3))
        end if
        call str2i (optpar(3),idum,ierr)
        if (ierr .ne. 0) goto 10
        if (idum .lt. 5) then
          call errcon ('Too few bins')
          goto 10
        end if
        nbins = idum
c
        call fillin (iset,maxset,maxhkl,fobs,sigfob,
     +               reso,buffer,nbins,ierr)
        if (ierr .ne. 0) goto 10
        call testrf (iset,maxset,maxhkl,rfree)
        change (iset) = .true.
c
c ... GUess
c
      else if (optpar(1)(1:2) .eq. 'GU') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'NResidues'
          call textin (' Guess what ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'NR') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) mw
            call pretty (optpar(3))
            call textin (' Mol Weight (Da) ?',optpar(3))
          end if
          call str2r (optpar(3),mw,ierr)
          if (ierr .ne. 0) goto 10
c
          nres = nint ( mw / 112.0 )
          call prompt (' Nres ~ MW / 112')
          call rvalut (' Mol Weight (Da) :',1,mw)
          call ivalut (' Nr of residues  ~',1,nres)
c
        else if (optpar(2)(1:2) .eq. 'MW') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) nres
            call pretty (optpar(3))
            call textin (' Nr of residues ?',optpar(3))
          end if
          call str2i (optpar(3),nres,ierr)
          if (ierr .ne. 0) goto 10
c
          mw = float(nres) * 112.0 
          call prompt (' MW ~ Nres * 112')
          call ivalut (' Nr of residues  ~',1,nres)
          call rvalut (' Mol Weight (Da) :',1,mw)
c
        else if (optpar(2)(1:2) .eq. 'VM') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (.not.know(kcell,iset)) then
            call errcon (' Cell has NOT been entered')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            write (optpar(4),*) nres
            call pretty (optpar(4))
            call textin (' Nr of residues ?',optpar(4))
          end if
          call str2i (optpar(4),nres,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            write (optpar(5),*) nasu
            call pretty (optpar(5))
            call textin (' Nr of asymm. units ?',optpar(5))
          end if
          call str2i (optpar(5),nasu,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            write (optpar(6),*) nncs
            call pretty (optpar(6))
            call textin (' Nr of NCS molecules ?',optpar(6))
          end if
          call str2i (optpar(6),nncs,ierr)
          if (ierr .ne. 0) goto 10
c
          call textut (' Set :',name(iset))
          call fvalut (' Cell constants    :',6,cell(1,iset))
          call ivalut (' Nr of residues    :',1,nres)
          call ivalut (' Asymm. units      :',1,nasu)
          call ivalut (' NCS molecules     :',1,nncs)
c
          call celvol (cell(1,iset),cvol)
          call rvalut (' Cell volume (A3)  :',1,cvol)
c
          x = 112.0 * float(nres) * float (nncs) * float (nasu)
          call prompt (' Assuming average residue mass = 112 Da')
          call prompt (' Mass ~ 112 * Nres * Nncs * Nasu')
          call rvalut (' Mass in cell (Da) ~',1,x)
          vm = cvol/x
          call prompt (' Vm = Volume / Mass')
          call fvalut (' Vm (A3/Da)        ~',1,vm)
c
          x = 140.0 * float (nres)
          call prompt (' Assuming average residue volume = 140 A3')
          call prompt (' Mol_vol ~ Nres * 140')
          call rvalut (' Mol volume (A3)   ~',1,x)
          x = 100.0 * x * float(nncs) * float (nasu) / cvol
          call fvalut (' Protein cntnt (%) ~',1,x)
          x = 100.0 - x
          call fvalut (' Solvent cntnt (%) ~',1,x)
          call fvalut (' 100%*(1-1.23/Vm)  ~',1,(100.0-(123.0/vm)))
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (.not.know(kcell,iset)) then
            call errcon (' Cell has NOT been entered')
            goto 10
          end if
c
          if (.not.know(kreso,iset)) then
            call errcon (' Resolution has NOT been calculated')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            write (optpar(4),*) reshi
            call pretty (optpar(4))
            call textin (' Resolution cut-off ?',optpar(4))
          end if
          call str2r (optpar(4),reshi,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            write (optpar(5),*) resol
            call pretty (optpar(5))
            call textin (' Resolution cut-off ?',optpar(5))
          end if
          call str2r (optpar(5),resol,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            optpar(6) = lattic
            call textin (' Lattice (P/R/C/F/I) ?',optpar(6))
          end if
          call upcase (optpar(6))
          if (index('PRCFI',optpar(6)(1:1)) .le. 0) goto 10
          lattic = optpar(6)(1:1)
c
          if (nopt .lt. 7) then
            write (optpar(7),*) nasu
            call pretty (optpar(7))
            call textin (' Nr of asymm. units ?',optpar(7))
          end if
          call str2i (optpar(7),nasu,ierr)
          if (ierr .ne. 0) goto 10
c
          call compl (iset,reshi,resol,lattic,nasu,
     +                   maxset,maxhkl,reso)
c
        else if (optpar(2)(1:2) .eq. 'RH') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (.not.know(kcell,iset)) then
            call errcon (' Cell has NOT been entered')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = lattic
            call textin (' Lattice (P/R/C/F/I) ?',optpar(4))
          end if
          call upcase (optpar(4))
          if (index('PRCFI',optpar(4)(1:1)) .le. 0) goto 10
          lattic = optpar(4)(1:1)
c
          if (nopt .lt. 5) then
            write (optpar(5),*) nasu
            call pretty (optpar(5))
            call textin (' Nr of asymm. units ?',optpar(5))
          end if
          call str2i (optpar(5),nasu,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            write (optpar(6),*) nncs
            call pretty (optpar(6))
            call textin (' Nr of NCS molecules ?',optpar(6))
          end if
          call str2i (optpar(6),nncs,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 7) then
            write (optpar(7),*) nres
            call pretty (optpar(7))
            call textin (' Nr of residues ?',optpar(7))
          end if
          call str2i (optpar(7),nres,ierr)
          if (ierr .ne. 0) goto 10
c
          call minres (iset,lattic,nasu,nncs,nres)
c
c ... room for more GUess commands
c
        else
          call errcon ('Invalid GUess command')
        end if
c
c ... EStimate_unique
c
      else if (optpar(1)(1:2) .eq. 'ES') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          write (optpar(3),*) resol
          call pretty (optpar(3))
          call textin (' Resolution (A) ?',optpar(3))
        end if
        call str2r (optpar(3),resol,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          optpar(4) = lattic
          call textin (' Lattice (P/R/C/F/I) ?',optpar(4))
        end if
        call upcase (optpar(4))
        if (index('PRCFI',optpar(4)(1:1)) .le. 0) goto 10
        lattic = optpar(4)(1:1)
c
        if (nopt .lt. 5) then
          write (optpar(5),*) nasu
          call pretty (optpar(5))
          call textin (' Nr of asymm. units ?',optpar(5))
        end if
        call str2i (optpar(5),nasu,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            if (know(kcell,iset)) then
              write (*,*)
              call textut (' Estimate unique reflections :',
     +          name(iset))
              call fvalut (' Unit cell axis lengths :',3,cell(1,iset))
              call fvalut (' Resolution limit (A)   :',1,resol)
              call asciut (' Lattice type           :',1,lattic)
              call ivalut (' Nr asymm. units/cell   :',1,nasu)
              call estuni (cell(1,iset),resol,lattic,nasu,nr)
              call jvalut (' Est. nr of reflections :',1,nr)
            else
              call errcon (' Unit cell has NOT been entered')
            end if
          end if
        end do
c
c ... EFfective_resolution
c
      else if (optpar(1)(1:2) .eq. 'EF') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = lattic
          call textin (' Lattice (P/R/C/F/I) ?',optpar(3))
        end if
        call upcase (optpar(3))
        if (index('PRCFI',optpar(3)(1:1)) .le. 0) goto 10
        lattic = optpar(3)(1:1)
c
        if (nopt .lt. 4) then
          write (optpar(4),*) nasu
          call pretty (optpar(4))
          call textin (' Nr of asymm. units ?',optpar(4))
        end if
        call str2i (optpar(4),nasu,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            if (know(kcell,iset)) then
              write (*,*)
              call textut (' Effective resolution :',
     +          name(iset))
              call fvalut (' Unit cell axis lengths :',3,cell(1,iset))
              call asciut (' Lattice type           :',1,lattic)
              call ivalut (' Nr asymm. units/cell   :',1,nasu)
              call effres (iset,lattic,nasu,
     +                   maxset,maxhkl,fobs,sigfob)
            else
              call errcon (' Unit cell has NOT been entered')
            end if
          end if
        end do
c
c ... DELETE
c
      else if (optpar(1)(1:2) .eq. 'DE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            if (change(iset)) then
              write (*,*)
              call textut (' There are unsaved changes to :',
     +          name(iset))
              if (linter) then
                reply = 'N'
                call textin (
     +            ' Really DELETE this set (Y/N) ?',reply)
                call upcase (reply)
                if (reply .ne. 'Y') goto 1210
              end if
            end if
c
            call textut (' Deleted :',name(iset))
c
            incore (iset) = .false.
            change (iset) = .false.
            coment (iset) = 'No comment'
            do i=1,maxitm
              know (i,iset) = .false.
            end do
            numhkl (iset) = 0
            nsymop (iset) = 0
            name (iset) = '&^%$%'
            file (iset) = 'not_saved_yet'
            do j=1,3
              cell (j,iset)   = 100.0
              cell (j+3,iset) = 90.0
            end do
c
 1210       continue
          end if
        end do
c
c ... SORT
c
      else if (optpar(1)(1:2) .eq. 'SO') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),kset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(kset)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which set ?',optpar(3))
        end if
        iset = whichm (optpar(3),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = 'LKH'
          call textin (' Sort order (fast-medium-slow) ?',optpar(4))
        end if
c
        call upcase (optpar(4))
        do i=1,3
          if (optpar(4)(i:i) .eq. 'H') then
            i123(i) = 1
          else if (optpar(4)(i:i) .eq. 'K') then
            i123(i) = 2
          else if (optpar(4)(i:i) .eq. 'L') then
            i123(i) = 3
          else
            call errcon ('Invalid sort order')
            call textut (' Input :',optpar(4))
            goto 10
          end if
        end do
c
        if ((i123(1)+i123(2)+i123(3)) .ne. 6) then
          call errcon ('In the sort order')
          call textut (' Input :',optpar(4))
          call ivalut (' Order :',3,i123)
          goto 10
        end if
c
        call textut (' Sort :',name(iset))
c
        call sortem (kset,iset,i123,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,
     +                   morbit,centri,ibuff)
        if (ierr .ne. 0) goto 10
c
        change (kset) = .true.
        incore (kset) = .true.
        coment (kset) = 'Sorted from '//name(iset)
        know (kdata,kset) = .true.
        call cpknow (iset,kset)
        call testrf (kset,maxset,maxhkl,rfree)
c
c ... DUPLICATE
c
      else if (optpar(1)(1:2) .eq. 'DU') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),kset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(kset)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which set ?',optpar(3))
        end if
        iset = whichm (optpar(3),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        call duplic (kset,iset,maxset,maxhkl,
     +               fobs,sigfob,reso,hkl,rfree,morbit,centri)
c
        change (kset) = .true.
        incore (kset) = .true.
        coment (kset) = 'Copied from '//name(iset)
        know (kdata,kset) = .true.
        call cpknow (iset,kset)
        call testrf (kset,maxset,maxhkl,rfree)
c
c ... HEMISPHERE/ASYM_UNIT
c
      else if (optpar(1)(1:2) .eq. 'HE' .or.
     +         optpar(1)(1:2) .eq. 'AS') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),kset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(kset)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which set ?',optpar(3))
        end if
        iset = whichm (optpar(3),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (.not. know(kcell,iset)) then
          call errcon ('Don''t know the cell parameters')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '2.0'
          call textin (' Resolution (A) ?',optpar(4))
        end if
        call str2r (optpar(4),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if ( optpar(1)(1:2) .eq. 'AS') then
          if (nopt .lt. 5) then
            optpar (5) = '3'
            call textin (' Laue group (? for list) ?',optpar(5))
          end if
c
          if (optpar(5)(1:1) .eq. '?') then
            do i=1,maxlau
              write (*,'(1x,a)') lautxt(i)(1:leng1(lautxt(i)))
            end do
            goto 10
          end if
c
          call str2i (optpar(5),idum,ierr)
          if (ierr .ne. 0) goto 10
c
          idum = max(1,min(idum,maxlau))
c
          optpar (6) = 'Asymmetric unit'
        else
          idum = 3
          optpar (6) = 'Hemisphere'
        end if
c
        call hemisp (kset,iset,maxset,maxhkl,xdum,idum,
     +               fobs,sigfob,reso,hkl,rfree,optpar(6),ierr)
c
        if (ierr .ne. 0) goto 10
c
        change (kset) = .true.
        incore (kset) = .true.
        coment (kset) = 'Hemisphere based on set '//name(iset)
        know (kdata,kset) = .true.
        know (kcell,kset) = .true.
        know (ksymm,kset) = .false.
        know (kreso,kset) = .true.
        know (kcent,kset) = .false.
        know (korbi,kset) = .false.
        know (kfree,kset) = .false.
        call testrf (kset,maxset,maxhkl,rfree)
c
c ... SIGMAS
c
      else if (optpar(1)(1:2) .eq. 'SI') then
c
        if (nopt .lt. 2) then
          optpar (2) = '?'
          call textin (' SIgmas option ?',optpar(2))
        end if
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'FA') then
c
c ... SIGMAS FAKE
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
c
          call textut (' Fake sigmas :',name(iset))
          call dosigs (1,iset,ierr,maxset,maxhkl,fobs,sigfob,
     +                 sigmin,sigmax)
          if (ierr .ne. 0) goto 10
c
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'CE') then
c
c ... SIGMAS CENTRIC_VS_ACENTRIC
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
c
          if (.not. know(kcent,iset)) then
            call errcon ('(A-)Centrics have not been deduced yet')
            goto 10
          end if
c
          call textut (' Sigmas centric vs acentric :',name(iset))
          call sigace (iset,maxset,maxhkl,fobs,sigfob,centri)
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
c ... SIGMAS LIMIT
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = '0.01'
            call textin (' Minimum sigma value ?',optpar(4))
          end if
          call str2r (optpar(4),sigmin,ierr)
          if (ierr .ne. 0) goto 10
          if (sigmin .le. 0.0) then
            call errcon ('Invalid minimum sigma value (<= 0.0)')
            goto 10
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = '1.0E30'
            call textin (' Maximum sigma value ?',optpar(5))
          end if
          call str2r (optpar(5),sigmax,ierr)
          if (ierr .ne. 0) goto 10
          if (sigmax .le. sigmin) then
            call errcon ('Invalid maximum sigma value (<= minimum)')
            goto 10
          end if
c
          call textut (' Limit sigmas :',name(iset))
          call dosigs (2,iset,ierr,maxset,maxhkl,fobs,sigfob,
     +                 sigmin,sigmax)
          if (ierr .ne. 0) goto 10
c
          change (iset) = .true.
c
c ... room for more options
c
        else
          call errcon ('Unknown SIgmas option')
        end if
c
c ... NOISE
c
      else if (optpar(1)(1:2) .eq. 'NO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = '15'
          call textin (' Number of bins ?',optpar(3))
        end if
        call str2i (optpar(3),nbins,ierr)
        if (ierr .ne. 0) goto 10
        if (nbins .lt. 3 .or. nbins .gt. 100) then
          call errcon ('Invalid number of bins (3-100)')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '2.5'
          call textin (' Minimum % noise ?',optpar(4))
        end if
        call str2r (optpar(4),xminoi,ierr)
        if (ierr .ne. 0) goto 10
        if (xminoi .lt. 0.0 .or. xminoi .gt. 100.0) then
          call errcon ('Invalid minimum noise level')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = '7.5'
          call textin (' Maximum % noise ?',optpar(5))
        end if
        call str2r (optpar(5),xmanoi,ierr)
        if (ierr .ne. 0) goto 10
        if (xmanoi .lt. 0.0 .or. xmanoi .gt. 100.0) then
          call errcon ('Invalid maximum noise level')
          goto 10
        end if
        call rlohi (xminoi,xmanoi)
c
        if (know(kreso,iset)) then
c
          call randno (iset,nbins,xminoi,xmanoi,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,reso,buffer,ibuff)
          if (ierr .ne. 0) goto 10
c
          change (iset) = .true.
c
        else
          call errcon ('I don''t know the resolution')
        end if
c
c ... LAUE
c
      else if (optpar(1)(1:2) .eq. 'LA') then
c
        if (nopt .lt. 2) then
          optpar (2) = 's99'
          call textin (' New set ?',optpar(2))
        end if
        call allocm (optpar(2),kset,ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = name(kset)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which set ?',optpar(3))
        end if
        iset = whichm (optpar(3),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '?'
          call textin (' Laue group (?=list) ?',optpar(4))
        end if
c
        if (optpar(4)(1:1) .eq. '?') then
          do i=1,maxlau
            write (*,'(1x,a)') lautxt(i)(1:leng1(lautxt(i)))
          end do
          goto 10
        end if
c
        call str2i (optpar(4),laue,ierr)
        if (ierr .ne. 0) goto 10
        if (laue .le. 0 .or. laue. gt. maxlau) then
          call errcon ('Invalid Laue group')
          goto 10
        end if
c
        call textut (' Laue old set :',name(iset))
        call textut ('      New set :',name(kset))
c
        if (know(ksymm,iset)) then
c
          call lauegr (kset,iset,laue,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,
     +                   morbit,centri)
          if (ierr .ne. 0) goto 10
c
          change (kset) = .true.
          incore (kset) = .true.
          coment (kset) = 'Laue from '//name(iset)
          know (kdata,kset) = .true.
          call cpknow (iset,kset)
          call testrf (kset,maxset,maxhkl,rfree)
c
        else
          call errcon ('I don''t know the symmops')
        end if
c
c ... CHANGE_INDEX
c
      else if (optpar(1)(1:2) .eq. 'CH') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = '+h'
          call textin (' New expression for h ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '+k'
          call textin (' New expression for k ?',optpar(4))
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = '+l'
          call textin (' New expression for l ?',optpar(5))
        end if
c
        call hklmat (optpar(3),optpar(4),optpar(5),
     +               mathkl,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Re-index :',name(iset))
            call reindx (numhkl(iset),hkl(1,1,iset),
     +                   mathkl,ierr)
            if (ierr .eq. 0) change (iset) = .true.
          end if
        end do
c
c ... RFREE
c
      else if (optpar(1)(1:2) .eq. 'RF') then
c
        if (nopt .lt. 2) then
          optpar(2) = 'LIST'
          call textin (' RFREE option ?',optpar(2))
          nopt = 2
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'IN') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) iseed
            call remspa (optpar(3))
            call textin (' Integer seed ?',optpar(3))
          end if
          call str2i (optpar(3),iseed,ierr)
          if (ierr .ne. 0) goto 10
c
          call gkrand (dx,0.0,0.0,iseed)
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          call selecm (optpar(3),ierr,maxset)
          if (ierr .ne. 0) goto 10
          prev = optpar (3)
c
          do iset=1,maxset
            if (select(iset)) then
              write (*,*)
              call textut (' Rfree :',name(iset))
              call testrf (iset,maxset,maxhkl,rfree)
            end if
          end do
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          call selecm (optpar(3),ierr,maxset)
          if (ierr .ne. 0) goto 10
          prev = optpar (3)
c
          do iset=1,maxset
            if (select(iset)) then
              write (*,*)
              call textut (' Rfree reset:',name(iset))
              do i=1,numhkl(iset)
                rfree(i,iset) = 0
              end do
              call testrf (iset,maxset,maxhkl,rfree)
              change (iset) = .true.
            end if
          end do
c
        else if (optpar(2)(1:2) .eq. 'AD') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = '10.0'
            call textin (' New % TEST reflections ?',
     +                   optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .lt. 0.0001 .or. xdum .gt. 99.999) then
            call errcon ('Invalid percentage')
            goto 10
          end if
c
          call rfadju (iset,maxset,maxhkl,rfree,xdum)
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'FI' .or.
     +           optpar(2)(1:2) .eq. 'CU' .or.
     +           optpar(2)(1:2) .eq. 'BI') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (.not. know(kreso,iset)) then
            call errcon ('Resolution not calculated yet')
            goto 10
          end if
c
          if (optpar(2)(1:2) .eq. 'BI') then
            if (nopt .lt. 4) then
              optpar (4) = '10'
              call textin (' Number of bins ?',
     +                     optpar(4))
            end if
            call str2i (optpar(4),idum,ierr)
            if (ierr .ne. 0) goto 10
            if (idum .lt. 3 .or. idum .gt. 99) then
              call errcon ('Invalid number of bins')
              goto 10
            end if
          else
            if (nopt .lt. 4) then
              optpar (4) = '10.0'
              call textin (' Target % TEST reflections ?',
     +                     optpar(4))
            end if
            call str2r (optpar(4),xdum,ierr)
            if (ierr .ne. 0) goto 10
            if (xdum .lt. 0.0001 .or. xdum .gt. 99.999) then
              call errcon ('Invalid percentage')
              goto 10
            end if
c
            if (nopt .lt. 5) then
              optpar (5) = '10'
              call textin (' Number of bins ?',
     +                     optpar(5))
            end if
            call str2i (optpar(5),idum,ierr)
            if (ierr .ne. 0) goto 10
            if (idum .lt. 3 .or. idum .gt. 99) then
              call errcon ('Invalid number of bins')
              goto 10
            end if
          end if
c
          call rfbins (iset,maxset,maxhkl,rfree,reso,buffer,
     +                 xdum,idum,optpar(2)(1:2))
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'GS') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = '10'
            call textin (' Flag every N-th hkl; value for N ?',
     +                   optpar(4))
          end if
          call str2i (optpar(4),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .lt. 5 .or. idum .gt. 250) then
            call errcon ('Invalid partitioning')
            goto 10
          end if
c
          do i=idum,numhkl(iset),idum
            rfree (i,iset) = 1
          end do
c
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'TR') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' TO which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = prev
            call textin (' FROM which set ?',optpar(4))
          end if
          jset = whichm (optpar(4),maxset)
c
          if (jset .le. 0 .or. jset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(jset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
c
          call textut (' Transferring TEST flags FROM :',name(jset))
          call testrf (jset,maxset,maxhkl,rfree)
          call textut (' TO :',name(iset))
          call testrf (iset,maxset,maxhkl,rfree)
c
          call prompt (' Encoding reflections of set 1 ...')
          do i=1,numhkl(iset)
            call packin (hkl(1,i,iset),hkl(2,i,iset),
     +        hkl(3,i,iset),0,ibuff(i))
            ibuff (maxhkl+i) = i
          end do
c
          call prompt (' Sorting reflections of set 1 ...')
          call shell (ibuff,ibuff(maxhkl+1),numhkl(iset))
c
          call prompt (' Transferring flags from 2 to 1 ...')
          do j=1,numhkl(jset)
            if (rfree(j,jset) .eq. 1) then
              call packin (hkl(1,j,jset),hkl(2,j,jset),
     +                     hkl(3,j,jset),0,idum)
c
              call bindex (idum,ibuff,ibuff(maxhkl+1),
     +          numhkl(iset),i)
              if (i .gt. 0) rfree (i,iset) = 1
            end if
          end do
c
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'SP') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = '10'
            call textin (' Percentage TEST data ?',
     +                   optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .lt. 0.001) then
            call errcon ('Invalid percentage')
            goto 10
          end if
          if (xdum .ge. 100.0) then
            xdum = 100.0 * xdum / float(numhkl(iset))
            call fvalut (' Converted to percentage :',1,xdum)
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = '1'
            call textin (' Reciprocal sphere radius ?',
     +                   optpar(5))
          end if
          call str2i (optpar(5),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .lt. 1 .or. idum .gt. 25) then
            call errcon ('Invalid radius')
            goto 10
          end if
c
          call prompt (' Encoding reflections ...')
          do i=1,numhkl(iset)
            call packin (hkl(1,i,iset),hkl(2,i,iset),
     +        hkl(3,i,iset),0,ibuff(i))
            ibuff (maxhkl+i) = i
          end do
c
          call prompt (' Sorting reflections ...')
          call shell (ibuff,ibuff(maxhkl+1),numhkl(iset))
c
          ldum = 0
          k = 0
 2359     continue
          call gkrand (dx,0.0,float(numhkl(iset)),0)
          j = 1 + int(dx)
          if (rfree(j,iset) .eq. 1) goto 2359
          ldum = ldum + 1
ccc          print *,j,hkl(1,j,iset),hkl(2,j,iset),hkl(3,j,iset)
          do i1=hkl(1,j,iset)-idum,hkl(1,j,iset)+idum
            do i2=hkl(2,j,iset)-idum,hkl(2,j,iset)+idum
              do i3=hkl(3,j,iset)-idum,hkl(3,j,iset)+idum
                dx = (float(i1-hkl(1,j,iset)))**2 +
     +               (float(i2-hkl(2,j,iset)))**2 +
     +               (float(i3-hkl(3,j,iset)))**2
                if (dx .le. float(idum*idum)) then
                  call packin (i1,i2,i3,0,jdum)
                  call bindex (jdum,ibuff,ibuff(maxhkl+1),
     +              numhkl(iset),kdum)
c
                  if (kdum .gt. 0) then
                    if (rfree(kdum,iset) .eq. 0) then
                      rfree (kdum,iset) = 1
                      k = k + 1
ccc              print *,' ADD ',i1,i2,i3
                    end if
                  end if
c
                end if
              end do
            end do
          end do
c
          dx = 100.0 * float(k) / float(numhkl(iset))
          if (dx .lt. xdum) goto 2359
c
          call jvalut (' Nr of TEST spheres :',1,ldum)
c
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = '10'
            call textin (' Number of partitionings ?',optpar(4))
          end if
          call str2i (optpar(4),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .lt. 5 .or. idum .gt. 250) then
            call errcon ('Invalid number of partitionings')
            goto 10
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = 'complete_xval_'
            call textin (' Basename for output files ?',optpar(5))
          end if
          call remspa (optpar(5))
c
c .. save current Rfree flags & generate test set index
c
          do i=1,numhkl(iset)
            ibuff (numhkl(iset)+i) = rfree (i,iset)
            call gkrand (dx,0.0,float(idum),0)
            ibuff (i) = 1 + int (dx)
          end do
c
          optpar (6) = 'RXPLOR'
          optpar (7) = '*'
c
          i1 = 0
          do i=1,idum
            k = 0
            do j=1,numhkl(iset)
              if (ibuff(j) .eq. i) then
                rfree (j,iset) = 1
                k = k + 1
              else
                rfree (j,iset) = 0
              end if
            end do
            i1 = i1 + k
            xdum = 100.0 * float(k) / float(numhkl(iset))
            write (line,'(a,i2,a)') optpar(5)(1:leng1(optpar(5))),
     +        i,'.rxplor'
            call remspa (line)
            write (*,2433) i,k,xdum,line(1:leng1(line))
c
            call dataut (iunit,line,optpar(6),optpar(7),maxhkl,
     +                   hkl(1,1,iset),numhkl(iset),fobs(1,iset),
     +                   sigfob(1,iset),rfree(1,iset),'ALL',
     +                   know(kcell,iset),cell(1,iset),
     +                   'BOTH',centri(1,iset),ierr)
c
          end do
c
          call jvalut (' Total nr of reflexions :',1,numhkl(iset))
          call jvalut (' Total TEST  reflexions :',1,i1)
c
c ... reset Rfree flags
c
          do i=1,numhkl(iset)
            rfree (i,iset) = ibuff (numhkl(iset)+i)
          end do
c
          call testrf (iset,maxset,maxhkl,rfree)
c
 2433 format (' Test set ',i2,' #hkl = ',i8,' = ',f8.2,' %'/
     +        ' ... file name = ',a)
c
        else if (optpar(2)(1:2) .eq. 'MU') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = '10'
            call textin (' Number of partitionings ?',optpar(4))
          end if
          call str2i (optpar(4),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .lt. 5 .or. idum .gt. 250) then
            call errcon ('Invalid number of partitionings')
            goto 10
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = 'R'
            call textin (' Random or Systematic (R|S) ?',optpar(5))
          end if
          call remspa (optpar(5))
          call upcase (optpar(5))
c
c .. generate test flags in range 1..NSETS
c
          do i=1,idum
            ibuff (i) = 0
          end do
c
          if (optpar(5)(1:1) .eq. 'S') then
            j = 0
            do i=1,numhkl(iset)
              j = j + 1
              rfree (i,iset) = j
              ibuff (j) = ibuff (j) + 1
              if (j .eq. idum) j = 0
            end do
          else
            do i=1,numhkl(iset)
              call gkrand (dx,0.0,float(idum),0)
              j = 1 + int(dx)
              rfree (i,iset) = j
              ibuff (j) = ibuff (j) + 1
            end do
          end if
c
          do i=1,idum
            xdum = 100.0 * float(ibuff(i)) / float(numhkl(iset))
            write (*,2437) i,ibuff(i),xdum
          end do
c
 2437 format (' Flag = ',i6,' for ',i8,' reflections (',f6.2,
     +        ' %)')
c
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'SU') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          write (*,*)
          call jvalut (' Total nr of reflections :',1,
     +      numhkl(iset))
c
          write (*,'(/99(1x,a,:,/))')
     +      'This command can help you decide how many, or',
     +      'what fraction of reflections should be set aside',
     +      'for cross-validation purposes. An estimate of the',
     +      'relative error in Rfree is provided, calculated as',
     +      '1 / SQRT (Nr HKLs).'
c
          write (*,'(/1x,a10,1x,a10,1x,a12)')
     +      '% HKLs','Nr HKLs','Rel. error'
c
          do i=1,15
            xdum = float(i) * float(numhkl(iset)) / 100.0
            j = nint (xdum)
            xdum2 = 1.0 / sqrt (float(j))
            write (*,'(1x,i10,1x,i10,1x,f10.3)') i,j,xdum2
          end do
c
          write (*,'(/1x,a10,1x,a10,1x,a12)')
     +      'Nr HKLs','% HKLs','Rel. error'
c
          do i=500,2500,100
            xdum = 100.0 * float(i) / float(numhkl(iset))
            xdum2 = 1.0 / sqrt (float(i))
            write (*,'(1x,i10,1x,f10.2,1x,f10.3)') i,xdum,xdum2
          end do
c
        else if (optpar(2)(1:2) .eq. 'GE') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            write (optpar(4),*) rfperc
            call remspa (optpar(4))
            call textin (' Percentage TEST data ?',optpar(4))
          end if
          call str2r (optpar(4),rfperc,ierr)
          if (ierr .ne. 0) goto 10
          if (rfperc .lt. 0.001) then
            call errcon ('Invalid percentage')
            goto 10
          end if
          if (rfperc .ge. 100.0) then
            rfperc = 100.0 * rfperc / float(numhkl(iset))
            call fvalut (' Converted to percentage :',1,rfperc)
          end if
c
          call textut (' Rfree generate:',name(iset))
          do i=1,numhkl(iset)
            rfree(i,iset) = 0
            call gkrand (dx,0.0,100.0,0)
            if (dx .le. rfperc) rfree(i,iset) = 1
          end do
          call testrf (iset,maxset,maxhkl,rfree)
          change (iset) = .true.
c
        else if (optpar(2)(1:2) .eq. 'SH') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Which set ?',optpar(3))
          end if
          iset = whichm (optpar(3),maxset)
c
          if (iset .le. 0 .or. iset .gt. maxset) then
            call errcon ('Invalid set selection')
            goto 10
          end if
c
          if (.not. incore(iset)) then
            call errcon ('Set not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            write (optpar(4),*) rfperc
            call remspa (optpar(4))
            call textin (' Percentage TEST data ?',optpar(4))
          end if
          call str2r (optpar(4),rfperc,ierr)
          if (ierr .ne. 0) goto 10
          if (rfperc .lt. 0.001) then
            call errcon ('Invalid percentage')
            goto 10
          end if
          if (rfperc .ge. 100.0) then
            rfperc = 100.0 * rfperc / float(numhkl(iset))
            call fvalut (' Converted to percentage :',1,rfperc)
          end if
c
          if (nopt .lt. 5) then
            write (optpar(5),*) nbins
            call textin (' Number of resolution bins ?',optpar(5))
          end if
          call str2i (optpar(5),nbins,ierr)
          if (ierr .ne. 0) goto 10
          if (nbins .lt. 3) then
            call errcon ('Must have more than two bins')
            goto 10
          end if
c
          call textut (' Rfree shell:',name(iset))
          if (know(kreso,iset)) then
            call rfshel (iset,rfperc,nbins,
     +                   maxset,maxhkl,maxbuf,
     +                   reso,hkl,rfree,buffer,ibuff)
            call testrf (iset,maxset,maxhkl,rfree)
            change (iset) = .true.
          else
            call errcon ('Resolution not calculated yet')
          end if
c
        else
          call errcon ('Invalid RFREE option')
          call textut (' Option :',optpar(2))
        end if
c
c ... STATS
c
      else if (optpar(1)(1:2) .eq. 'ST') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Stats :',name(iset))
            call datast (numhkl(iset),hkl(1,1,iset),
     +                   fobs(1,iset),sigfob(1,iset),
     +                   reso(1,iset),know(kreso,iset),
     +                   morbit(1,iset),know(korbi,iset),
     +                   buffer)
            call testrf (iset,maxset,maxhkl,rfree)
          end if
        end do
c
c ... SPECIAL
c
      else if (optpar(1)(1:2) .eq. 'SP') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = '0k0'
          call textin (' HKL-type (00l, h0l, etc.) ?',optpar(3))
        end if
c
        call upcase (optpar(3))
        call remspa (optpar(3))
        if ( index('HKL0',optpar(3)(1:1)) .le. 0 .or.
     +       index('HKL0',optpar(3)(2:2)) .le. 0 .or.
     +       index('HKL0',optpar(3)(3:3)) .le. 0) then
          call errcon ('Unrecognised character in HKL-type')
          goto 10
        end if
c
c        if ( (optpar(3)(1:1) .ne. 'H' .and.
c     +        optpar(3)(1:1) .ne. '0') .or.
c     +       (optpar(3)(2:2) .ne. 'K' .and.
c     +        optpar(3)(2:2) .ne. '0') .or.
c     +       (optpar(3)(3:3) .ne. 'L' .and.
c     +        optpar(3)(3:3) .ne. '0') ) then
c          call errcon ('Invalid HKL-type entered')
c          goto 10
c        end if
c
c        j = 0
c        do i=1,3
c          if (optpar(3)(i:i) .eq. '0') j=j+1
c        end do
c        if (j.eq.0) then
c          call errcon ('Cannot list ALL reflections')
c          goto 10
c        end if
c
        do iset=1,maxset
          if (select(iset)) then
c
            write (*,*)
            call textut (' Special  :',name(iset))
            call textut (' HKL-type :',optpar(3))
c
            call hklspe (iset,optpar(3),
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri)
c
          end if
        end do
c
c ... TYPE_HKL
c
      else if (optpar(1)(1:2) .eq. 'TY') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        i123(1) = 1
        i123(2) = 100
        i123(3) = 10
c
        do i=1,3
          if (nopt.lt.(i+2)) then
            write (optpar(i+2),*) i123(i)
            call textin (' Index :',optpar(i+2))
          end if
          call str2i (optpar(i+2),i123(i),ierr)
          if (ierr .ne. 0) goto 10
        end do
c
        do iset=1,maxset
          if (select(iset)) then
c
            write (*,*)
            call textut (' Type :',name(iset))
c
            i1 = max(1,i123(1))
            i2 = i123(2)
            if (i2 .le. 0) i2=numhkl(iset)
            i2 = max(i1,min(i2,numhkl(iset)))
            i3 = i123(3)
            if (i3 .eq. 0) i3 = i2-1
            if (i3 .lt. 0) i3=(i2-i1+1)/i3
            write (*,'(a,i8,a,i8,a,i8)')
     +        ' Type from ',i1,' to ',i2,' in steps of ',i3
c
            write (*,6901) 'Reflxn','   H','   K','   L',
     +        'Resoln','   Fobs   ',' Sigma(Fobs)',' F/Sig',
     +        'Test','A/C','Orb'
c
            do i=i1,i2,i3
              xdum = -1.0
              if (sigfob(i,1) .ne. 0.0) then
                xdum = fobs(i,1)/sigfob(i,1)
              end if
              write (*,6900) i,(hkl(j,i,1),j=1,3),reso(i,1),
     +          fobs(i,1),sigfob(i,1),xdum,rfree(i,1),
     +          centri(i,1),morbit(i,1)
ccc              write (*,6990) i,(hkl(j,i,iset),j=1,3),
ccc     +          fobs(i,iset),sigfob(i,iset),rfree(i,iset)
            end do
c
          end if
        end do
c
 6900 format (1x,i6,1x,3i4,1x,f6.2,1x,1p,2e12.4,1x,0p,f6.2,
     +  1x,i4,3x,a1,1x,i3)
 6901 format (1x,a6,1x,3a4,1x,a6,1x,2a12,1x,a6,1x,a4,1x,a3,1x,a3)
 6990 format (' # ',i8,' HKL ',3i6,' Fo, S(Fo) = ',1p,2e12.4,
     +  ' Test ',i1)
c
c ... LIST
c
      else if (optpar(1)(1:2) .eq. 'LI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
 6100 format (
     +  ' Number of reflections : ',i10/
     +  ' File name : ',a/
     +  ' Label     : ',a)
 6110 format (
     +  ' Cell      : ',6f10.3)
 6120 format (
     +  ' Nr of symmetry operators : ',i3)
c
        do iset=1,maxset
          if (select(iset)) then
c
            write (*,*)
            call textut (' List :',name(iset))
            write (*,6100) numhkl(iset),
     +        file(iset)(1:leng1(file(iset))),
     +        coment(iset)(1:leng1(coment(iset)))
c
            if (know(kcell,iset)) then
              write (*,6110) (cell(i,iset),i=1,6)
            else
              call prompt (' Cell constants not supplied')
            end if
c
            if (know(ksymm,iset)) then
              write (*,6120) nsymop(iset)
            else
              call prompt (' Symmetry operators not supplied')
            end if
c
            if (know(kreso,iset)) then
              call prompt (' Resolution has been calculated')
            else
              call prompt (' Resolution has NOT been calculated')
            end if
c
            if (know(kcent,iset)) then
              call prompt (' (A)centrics have been deduced')
            else
              call prompt (' (A)centrics have NOT been deduced')
            end if
c
            if (know(korbi,iset)) then
              call prompt (' Orbital multiplicities calculated')
            else
              call prompt (' Orbital multiplicities NOT calculated')
            end if
c
            if (know(kfree,iset)) then
              call prompt (' This is a cross-validation dataset')
            else
              call prompt (' This is NOT a cross-validation dataset')
            end if
c
            if (change(iset)) then
              call prompt (' There are UNSAVED changes')
            else
              call prompt (' There are no unsaved changes')
            end if
c
          end if
        end do
c
c ... SHOW_HKL
c
      else if (optpar(1)(1:2) .eq. 'SH') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'FOB'
          call textin (' Criterion [FOB|SIG|F/S|RES] ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = '>'
          call textin (' Operand [<|>| ?',optpar(4))
        end if
c
        if (nopt .lt. 5) then
          optpar(5) = '0.0'
          call textin (' Value ?',optpar(5))
        end if
        call str2r (optpar(5),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Show_hkl :',name(iset))
            call kilhkl ('S',iset,optpar(3),optpar(4),xdum,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
          end if
        end do
c
c ... ROGUE_KILL
c
      else if (optpar(1)(1:2) .eq. 'RO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = '0'
          call textin (' Index H of this rogue ?',optpar(3))
        end if
        call str2i (optpar(3),i,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          optpar(4) = '0'
          call textin (' Index K of this rogue ?',optpar(4))
        end if
        call str2i (optpar(4),j,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 5) then
          optpar(5) = '0'
          call textin (' Index L of this rogue ?',optpar(5))
        end if
        call str2i (optpar(5),k,ierr)
        if (ierr .ne. 0) goto 10
c
        inow = 5
c
 6262   continue
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Rogue_kill :',name(iset))
            call kilrog (iset,i,j,k,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
            if (ierr .eq. 0) change (iset) = .true.
          end if
        end do
c
        if (nopt .ge. (inow+3)) then
          call str2i (optpar(inow+1),i,ierr)
          call str2i (optpar(inow+2),j,ierr)
          call str2i (optpar(inow+3),k,ierr)
          if (ierr .ne. 0) goto 10
          inow = inow + 3
          goto 6262
        end if
c
        do iset=1,maxset
          if (select(iset) .and. change(iset)) then
            write (*,*)
            call textut (' Checking Rfree flags of :',name(iset))
            call testrf (iset,maxset,maxhkl,rfree)
          end if
        end do
c
c ... KILL_HKL
c
      else if (optpar(1)(1:2) .eq. 'KI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'FOB'
          call textin (' Criterion [FOB|SIG|F/S|RES] ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = '>'
          call textin (' Operand [<|>| ?',optpar(4))
        end if
c
        if (nopt .lt. 5) then
          optpar(5) = '0.0'
          call textin (' Value ?',optpar(5))
        end if
        call str2r (optpar(5),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Kill_hkl :',name(iset))
            call kilhkl ('K',iset,optpar(3),optpar(4),xdum,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
            if (ierr .eq. 0) change (iset) = .true.
            call testrf (iset,maxset,maxhkl,rfree)
          end if
        end do
c
c ... ODD_KILL / EVEN_KILL
c
      else if (optpar(1)(1:2) .eq. 'OD' .or.
     +         optpar(1)(1:2) .eq. 'EV') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'L'
          call textin (' Use H, K or L ?',optpar(3))
        end if
c
        call upcase (optpar(3))
        if (optpar(3)(1:1) .ne. 'H' .and.
     +      optpar(3)(1:1) .ne. 'K' .and.
     +      optpar(3)(1:1) .ne. 'L') then
          call errcon ('You must select H, K or L')
          goto 10
        end if
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            if (optpar(1)(1:2) .eq. 'OD') then
              call textut (' Odd kill :',name(iset))
            else
              call textut (' Even kill :',name(iset))
            end if
            call oddevn (optpar(1),iset,optpar(3),ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
            if (ierr .eq. 0) change (iset) = .true.
            call testrf (iset,maxset,maxhkl,rfree)
          end if
        end do
c
c ... TEMP_FACTOR
c
      else if (optpar(1)(1:2) .eq. 'TE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = '0.0'
          call textin (' Temperature factor ?',optpar(3))
        end if
        call str2r (optpar(3),xdum2,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Temperature factor :',name(iset))
c
            if (know(kreso,iset)) then
              do i=1,numhkl(iset)
                xdum = 0.5 / reso(i,iset)
                xdum = -xdum2*xdum*xdum
                xdum = exp(xdum)
                fobs (i,iset) = fobs (i,iset) * xdum
              end do
            else
              call errcon ('Resolution not calculated')
            end if
c
            change (iset) = .true.
c              
          end if
        end do
c
c ... PROD_PLUS
c
      else if (optpar(1)(1:2) .eq. 'PR') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'FOB'
          call textin (' Which [Fobs|Sigfo|Both] ?',optpar(3))
        end if
        call upcase (optpar(3))
        if (optpar(3)(1:1) .ne. 'F' .and.
     +      optpar(3)(1:1) .ne. 'S' .and.
     +      optpar(3)(1:1) .ne. 'B') then
          call errcon ('Invalid selection')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = '1.0'
          call textin (' Prod ?',optpar(4))
        end if
        call str2r (optpar(4),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 5) then
          optpar(5) = '0.0'
          call textin (' Plus ?',optpar(5))
        end if
        call str2r (optpar(5),xdum2,ierr)
        if (ierr .ne. 0) goto 10
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Prod_plus :',name(iset))
c
            if (optpar(3)(1:1) .eq. 'F' .or.
     +          optpar(3)(1:1) .eq. 'B') then
              do i=1,numhkl(iset)
                fobs(i,iset) = xdum*fobs(i,iset) + xdum2
              end do
            end if
c
            if (optpar(3)(1:1) .eq. 'S' .or.
     +          optpar(3)(1:1) .eq. 'B') then
              do i=1,numhkl(iset)
                sigfob(i,iset) = xdum*sigfob(i,iset) + xdum2
              end do
            end if
c
            change (iset) = .true.
c              
          end if
        end do
c
c ... HISTOGRAM
c
      else if (optpar(1)(1:2) .eq. 'HI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'FOB'
          call textin (' Which [FOB|SIG|F/S|RES] ?',optpar(3))
        end if
        call upcase (optpar(3))
        if (optpar(3)(1:3) .ne. 'FOB' .and.
     +      optpar(3)(1:3) .ne. 'SIG' .and.
     +      optpar(3)(1:3) .ne. 'F/S' .and.
     +      optpar(3)(1:3) .ne. 'RES') then
          call errcon ('Invalid selection')
          goto 10
        end if
c
        numhis = 0
        do i=1,3
          if (nopt .lt. (i+3)) then
            write (optpar(i+3),*) float(i)
            call textin (' Histogram value ?',optpar(i+3))
          end if
          call str2r (optpar(i+3),xdum,ierr)
          if (ierr .ne. 0) goto 10
          numhis = numhis + 1
          rhis (numhis) = xdum
        end do
c
        if (nopt .gt. 6) then
          do i=7,nopt
            call str2r (optpar(i),xdum,ierr)
            if (ierr .ne. 0) goto 10
            numhis = numhis + 1
            rhis (numhis) = xdum
          end do
        end if
c
ccc        call qsortg (rhis,numhis)
        call hsortr (numhis,rhis)
        call rvalut (' Histogram limits :',numhis,rhis)
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Histogram :',name(iset))
c
            if (optpar(3)(1:3) .eq. 'FOB') then
              call histo (numhkl(iset),fobs(1,iset),
     +                    numhis,rhis,nhis)
            else if (optpar(3)(1:3) .eq. 'SIG') then
              call histo (numhkl(iset),sigfob(1,iset),
     +                    numhis,rhis,nhis)
            else if (optpar(3)(1:3) .eq. 'F/S') then
              do i=1,numhkl(iset)
                if (sigfob(i,iset) .ge. small) then
                  buffer (i) = fobs(i,iset)/sigfob(i,iset)
                else
                  buffer (i) = 0.0
                end if
              end do
              call histo (numhkl(iset),buffer,
     +                    numhis,rhis,nhis)
            else if (optpar(3)(1:3) .eq. 'RES') then
              if (know(kreso,iset)) then
                call histo (numhkl(iset),reso(1,iset),
     +                      numhis,rhis,nhis)
              else
                call errcon ('Resolution not calculated')
              end if
            end if
          end if
        end do
c
c ... SYMMOP
c
      else if (optpar(1)(1:2) .eq. 'SY') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'symop.o'
          call textin (' O file with symmetry operators ?',
     +      optpar(3))
        end if
c
        call osymop (iunit,optpar(3),ierr)
        if (ierr .ne. 0) then
          call errcon ( 'While opening O datablock file')
          goto 10
        end if
        close (iunit)
c
        call opoodb (iunit,optpar(3),parnam,partyp,numpar,
     +    parfmt,ierr)
        if (ierr .ne. 0) then
          close (iunit)
          goto 10
        end if
c
        i = numpar/12
        call jvalut (' Nr of symmetry operators :',1,i)
        if (i .lt. 1 .or. i .gt. maxsym) then
          call errcon ('Invalid number of symmetry operators')
          close (iunit)
          goto 10
        end if
c
        nsymop (iset) = i
        read (iunit,parfmt,err=69,end=69)
     +    ((symmop(i,j,iset),i=1,12),j=1,nsymop(iset))
        close (iunit)
c
        call anasgs (nsymop(iset),symmop(1,1,iset),.true.,ierr)
        if (ierr .ne. 0) then
          call errcon ('In spacegroup symmetry operators')
          nsymop (iset) = 0
          goto 10
        end if
c
        know (ksymm,iset) = .true.
c
cc        do i=1,nsymop(iset)
cc          write (*,6200) i,symmop(1,i,iset),symmop(4,i,iset),
cc     +      symmop(7,i,iset),symmop(10,i,iset),
cc     +      symmop(2,i,iset),symmop(5,i,iset),
cc     +      symmop(8,i,iset),symmop(11,i,iset),
cc     +      symmop(3,i,iset),symmop(6,i,iset),
cc     +      symmop(9,i,iset),symmop(12,i,iset)
cc        end do
c
c ... find primitive symm-ops (unique rotation matrix)
c
        nuniq (iset) = 0
        do i=1,nsymop(iset)
          if (nuniq(iset) .le. 0) goto 110
          do j=1,nuniq(iset)
            do k=1,9
              if (symmop(k,i,iset) .ne. symmop(k,j,iset)) goto 111
            end do
c
c ... not unique
c
            write (*,6210) i,j
            goto 112
c
  111       continue
          end do
c
c ... unique; create it's transpose
c
  110     continue
          nuniq (iset) = nuniq (iset) + 1
          uniqso (nuniq(iset),iset) = i
c
          transp (1,nuniq(iset),iset) = nint(symmop (1,i,iset))
          transp (2,nuniq(iset),iset) = nint(symmop (4,i,iset))
          transp (3,nuniq(iset),iset) = nint(symmop (7,i,iset))
          transp (4,nuniq(iset),iset) = nint(symmop (2,i,iset))
          transp (5,nuniq(iset),iset) = nint(symmop (5,i,iset))
          transp (6,nuniq(iset),iset) = nint(symmop (8,i,iset))
          transp (7,nuniq(iset),iset) = nint(symmop (3,i,iset))
          transp (8,nuniq(iset),iset) = nint(symmop (6,i,iset))
          transp (9,nuniq(iset),iset) = nint(symmop (9,i,iset))
c
  112     continue
c
        end do
c
        call ivalut (' Unique symmops :',nuniq(iset),uniqso(1,iset))
c
        goto 10
c
 6200 format (' # ',i2,1x,4f10.3,2(/6x,4f10.3))
 6210 format (' Rot matrix of symmop ',i2,' identical to ',i2)
c
   69   continue
        call errcon ('While reading symmetry operators')
        close (iunit)
        goto 10
c
c ... CELL
c
      else if (optpar(1)(1:2) .eq. 'CE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        do i=1,6
          if (nopt .lt. (i+2)) then
            write (optpar(i+2),'(f15.3)') cell(i,iset)
            call remspa (optpar(i+2))
            write (line,'(a,i1,a)') ' Value for cell ',i,' ?'
            call textin (line,optpar(i+2))
          end if
          call str2r(optpar(i+2),xdum,ierr)
          if (ierr .ne. 0) goto 10
          cell (i,iset) = xdum
        end do
c
        call fvalut (' Cell :',6,cell(1,iset))
        call celvol (cell(1,iset),cvol)
        call rvalut (' Volume (A3) :',1,cvol)
        know (kcell,iset) = .true.
c
c ... ANnotate
c
      else if (optpar(1)(1:2) .eq. 'AN') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = coment(iset)
          call textin (' Label ?',optpar(3))
        end if
c
        coment (iset) = optpar (3)
c
c ... CALC
c
      else if (optpar(1)(1:2) .eq. 'CA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        call selecm (optpar(2),ierr,maxset)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'R'
          call textin (' What [Res|Cent|Orb|F2I|I2F] ?',optpar(3))
        end if
        call upcase (optpar(3))
        if (optpar(3)(1:1) .ne. 'R' .and.
     +      optpar(3)(1:1) .ne. 'C' .and.
     +      optpar(3)(1:1) .ne. 'F' .and.
     +      optpar(3)(1:1) .ne. 'I' .and.
     +      optpar(3)(1:1) .ne. 'O') then
          call errcon ('Invalid selection')
          goto 10
        end if
c
        do iset=1,maxset
          if (select(iset)) then
            write (*,*)
            call textut (' Calc :',name(iset))
c
c ... RESOLUTION
c
            if (optpar(3)(1:1) .eq. 'R') then
              if (.not. know(kcell,iset)) then
                call errcon ('I don''t know the cell constants')
              else
                call calres (iset,ierr,
     +                   maxset,maxhkl,reso,hkl)
                if (ierr .eq. 0) know(kreso,iset) = .true.
              end if
c
c ... F => GO FROM Fs TO Is
c
            else if (optpar(3)(1:1) .eq. 'F') then
              nsneg = 0
              do i=1,numhkl(iset)
                if (fobs(i,iset) .gt. 0.0) then
                  sigfob(i,iset)=2.0*fobs(i,iset)*sigfob(i,iset)
                  fobs(i,iset)=fobs(i,iset)**2
                else
                  nsneg = nsneg + 1
                end if
              end do
              call prompt (' Fnow =  "I"   = Fold ^ 2')
              call prompt (' Snow = "S(I)" = 2 * Fold * Sold')
              if (nsneg .gt. 0) then
                call errcon ('There were zero and/or negative Fs')
                call jvalut (' Unchanged Fs <= 0 :',1,nsneg)
              end if
              change (iset) = .true.
c
c ... I => GO FROM Is TO Fs
c
            else if (optpar(3)(1:1) .eq. 'I') then
              nfneg = 0
              do i=1,numhkl(iset)
                if (fobs(i,iset) .gt. 0.0) then
                  sigfob(i,iset) = sigfob(i,iset) /
     +                         (2.0*sqrt(fobs(i,iset)))
                  fobs(i,iset)=sqrt(fobs(i,iset))
                else
                  nfneg = nfneg + 1
                end if
              end do
              call prompt (' Fnow = SQRT (Fold)')
              call prompt (' Snow = Sold / ( 2 * SQRT(Fold) )')
              if (nfneg .gt. 0) then
                call errcon ('There were zero and/or negative Is')
                call jvalut (' Unchanged Is <= 0 :',1,nfneg)
              end if
              change (iset) = .true.
c
c ... CENTRIC/ACENTRIC
c
            else if (optpar(3)(1:1) .eq. 'C') then
              if (.not. know(ksymm,iset)) then
                call errcon ('I don''t know the symmops')
              else
                call calcen (iset,ierr,
     +                   maxset,maxhkl,hkl,centri)
                if (ierr .eq. 0) know(kcent,iset) = .true.
              end if
c
c ... ORBITAL MULTIPLICITY
c
            else if (optpar(3)(1:1) .eq. 'O') then
              if (.not. know(ksymm,iset)) then
                call errcon ('I don''t know the symmops')
              else
                call calorb (iset,ierr,
     +                   maxset,maxhkl,hkl,morbit)
                if (ierr .eq. 0) know(korbi,iset) = .true.
              end if
            end if
c
          end if
        end do
c
c ... SCATTER_PLOT
c
      else if (optpar(1)(1:2) .eq. 'SC') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'scatter.plt'
          call textin (' O2D plot file ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = '?'
          call textin (' Horizontal ?',optpar(4))
        end if
c
        if (optpar(4)(1:1) .eq. '?') then
          write (*,'(99(a:/))')
     +      ' Select one of : FOB = Fobs, SIG = Sigma,',
     +      '   F/S = Fobs/Sigma, INT = F^2,',
     +      '   I/S = F^2/Sigma^2, RES = resolution,',
     +      '   1/R = 1/resolution, STL = sin(theta)/lambda,',
     +      '   DST = 4(STL^2), LNI = ln(F^2)'
          goto 10
        end if
c
        call upcase (optpar(4))
        if (index ('FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|LNI',
     +             optpar(4)(1:3)) .le. 0) then
          call errcon ('Invalid plot variable')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar(5) = '?'
          call textin (' Vertical ?',optpar(5))
        end if
c
        if (optpar(5)(1:1) .eq. '?') then
          write (*,'(99(a:/))')
     +      ' Select one of : FOB = Fobs, SIG = Sigma,',
     +      '   F/S = Fobs/Sigma, INT = F^2,',
     +      '   I/S = F^2/Sigma^2, RES = resolution,',
     +      '   1/R = 1/resolution, STL = sin(theta)/lambda,',
     +      '   DST = 4(STL^2), LNI = ln(F^2)'
          goto 10
        end if
c
        call upcase (optpar(5))
        if (index ('FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|LNI',
     +             optpar(5)(1:3)) .le. 0) then
          call errcon ('Invalid plot variable')
          goto 10
        end if
c
        if (index ('RES|1/R|STL|DST',optpar(4)(1:3)) .gt. 0 .or.
     +      index ('RES|1/R|STL|DST',optpar(5)(1:3)) .gt. 0) then
          if (.not. know(kcell,iset)) then
            call errcon ('I don''t know the cell constants')
            goto 10
          end if
        end if
c
        call scattr (iset,optpar(3),optpar(4),optpar(5),ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,rfree,buffer)
c
c ... HKl_aniso_plot
c
      else if (optpar(1)(1:2) .eq. 'HK') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'aniso.plt'
          call textin (' O2D plot file ?',optpar(3))
        end if
c
        if (.not. know(kreso,iset)) then
          call errcon ('I don''t know the resolution')
          goto 10
        end if
c
        call anisop (iset,optpar(3),ierr,
     +               maxset,maxhkl,fobs,reso,hkl)
c
c ... BIN_PLOT
c
      else if (optpar(1)(1:2) .eq. 'BI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'bin.plt'
          call textin (' O2D plot file ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = '?'
          call textin (' Horizontal ?',optpar(4))
        end if
c
        if (optpar(4)(1:1) .eq. '?') then
          write (*,'(99(a:/))')
     +      ' Select one of : FOB = Fobs, SIG = Sigma,',
     +      '   F/S = Fobs/Sigma, INT = F^2,',
     +      '   I/S = F^2/Sigma^2, RES = resolution,',
     +      '   1/R = 1/resolution, STL = sin(theta)/lambda,',
     +      '   DST = 4(STL^2), LNI = ln(F^2), H, K, L,',
     +      '   H/A, K/B, L/C'
          goto 10
        end if
c
        call upcase (optpar(4))
        if (index ('FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|LNI',
     +             optpar(4)(1:3)) .le. 0 .and.
     +      index ('H  |K  |L  |H/A|K/B|L/C',
     +             optpar(4)(1:3)) .le. 0) then
          call errcon ('Invalid plot variable')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar(5) = '?'
          call textin (' Vertical ?',optpar(5))
        end if
c
        if (optpar(5)(1:1) .eq. '?') then
          write (*,'(99(a:/))')
     +      ' Select one of : FOB = Fobs, SIG = Sigma,',
     +      '   F/S = Fobs/Sigma, INT = F^2,',
     +      '   I/S = F^2/Sigma^2, RES = resolution,',
     +      '   1/R = 1/resolution, STL = sin(theta)/lambda,',
     +      '   DST = 4(STL^2), NRF = nr of reflections',
     +      '   LNI = ln(F^2)'
          goto 10
        end if
c
        call upcase (optpar(5))
        if (index ('FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|NRF|LNI',
     +             optpar(5)(1:3)) .le. 0) then
          call errcon ('Invalid plot variable')
          goto 10
        end if
c
        if (index ('RES|1/R|STL|DST',optpar(4)(1:3)) .gt. 0 .or.
     +      index ('RES|1/R|STL|DST',optpar(5)(1:3)) .gt. 0 .or.
     +      index ('H/A|K/B|L/C',optpar(4)(1:3)) .gt. 0) then
          if (.not. know(kcell,iset)) then
            call errcon ('I don''t know the cell constants')
            goto 10
          end if
        end if
c
        if (nopt .lt. 6) then
          optpar(6) = '-15'
          call textin (' Bin size (>0) or Nr bins (<0) ?',optpar(6))
        end if
        call str2r (optpar(6),dpbin,ierr)
        if (ierr .ne. 0) goto 10
c
        call binplt (iset,optpar(3),optpar(4),optpar(5),dpbin,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,buffer)
c
c ... EO_PLOT
c
      else if (optpar(1)(1:2) .eq. 'EO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar(3) = 'even_odd.plt'
          call textin (' O2D plot file ?',optpar(3))
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'L'
          call textin (' H, K, or L ?',optpar(4))
        end if
c
        call upcase (optpar(4))
        if (index ('HKL',optpar(4)(1:1)) .le. 0) then
          call errcon ('Invalid selection; assuming L')
          optpar(4) = 'L'
        end if
c
        call eoplot (iset,optpar(3),optpar(4),ierr,
     +               maxset,maxhkl,maxbuf,fobs,hkl)
c
c ... DOUBLE_PLOT (used to be called DUO_PLOT)
c
      else if (optpar(1)(1:2) .eq. 'DO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which set1 ?',optpar(2))
        end if
        iset = whichm (optpar(2),maxset)
c
        if (iset .le. 0 .or. iset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(iset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Which set2 ?',optpar(3))
        end if
        jset = whichm (optpar(3),maxset)
c
        if (jset .le. 0 .or. jset .gt. maxset) then
          call errcon ('Invalid set selection')
          goto 10
        end if
c
        if (.not. incore(jset)) then
          call errcon ('Set not in memory')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'duo.plt'
          call textin (' O2D plot file ?',optpar(4))
        end if
c
        if (nopt .lt. 5) then
          optpar(5) = '?'
          call textin (' Horizontal ?',optpar(5))
        end if
c
        if (optpar(5)(1:1) .eq. '?') then
          write (*,'(99(a:/))')
     +      ' Select one of : FOB = Fobs, SIG = Sigma,',
     +      '   F/S = Fobs/Sigma, INT = F^2,',
     +      '   I/S = F^2/Sigma^2, RES = resolution,',
     +      '   1/R = 1/resolution, STL = sin(theta)/lambda,',
     +      '   DST = 4(STL^2), LNI = ln(F^2)'
          goto 10
        end if
c
        call upcase (optpar(5))
        if (index ('FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|LNI',
     +             optpar(5)(1:3)) .le. 0) then
          call errcon ('Invalid plot variable')
          goto 10
        end if
c
        if (nopt .lt. 6) then
          optpar(6) = '?'
          call textin (' Vertical ?',optpar(6))
        end if
c
        if (optpar(6)(1:1) .eq. '?') then
          write (*,'(99(a:/))')
     +      ' Select one of : FOB = Fobs, SIG = Sigma,',
     +      '   F/S = Fobs/Sigma, INT = F^2,',
     +      '   I/S = F^2/Sigma^2, RES = resolution,',
     +      '   1/R = 1/resolution, STL = sin(theta)/lambda,',
     +      '   DST = 4(STL^2), NRF = nr of reflections',
     +      '   LNI = ln(F^2)'
          goto 10
        end if
c
        call upcase (optpar(6))
        if (index ('FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|NRF|LNI',
     +             optpar(6)(1:3)) .le. 0) then
          call errcon ('Invalid plot variable')
          goto 10
        end if
c
        if (index ('RES|1/R|STL|DST',optpar(5)(1:3)) .gt. 0 .or.
     +      index ('RES|1/R|STL|DST',optpar(6)(1:3)) .gt. 0) then
          if (.not. know(kcell,iset)) then
            call errcon ('I don''t know the cell constants')
            goto 10
          end if
        end if
c
        if (nopt .lt. 7) then
          optpar(7) = '-15'
          call textin (' Bin size (>0) or Nr bins (<0) ?',optpar(7))
        end if
        call str2r (optpar(7),dpbin,ierr)
        if (ierr .ne. 0) goto 10
c
        call duoplt (iset,jset,optpar(4),optpar(5),optpar(6),
     +               dpbin,ierr,maxset,maxhkl,maxbuf,
     +               fobs,sigfob,reso,buffer)
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
      do i=1,maxset
        if (incore(i)) then
          if (change(i)) then
              call textut (' There are unsaved changes to :',
     +          name(i))
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
      return
c
      end
c
c
c
      subroutine allocm (string,imsk,ierr,maxset)
c
c ... allocate a dataset
c
      include 'dataman.incl'
c
      integer maxset
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
      do i=1,maxset
        if (.not. incore(i)) then
          imsk = i
          goto 910
        end if
      end do
c
      call errcon ('No more sets available')
      return
c
  910 continue
      name (imsk) = string
c
      call upcase (name(imsk))
      if (whichm(name(imsk),maxset) .ne. imsk) then
        call errcon ('Invalid set name (empty or not unique)')
        name (imsk) = '*&^$'
        return
      end if
c
      if (incore(imsk)) then
        call errcon ('Set in use; DELETE it first')
        return
      end if
c
      numhkl (imsk) = 0
      nsymop (imsk) = 0
      file   (imsk) = 'not_saved_yet'
      coment (imsk) = 'No comment'
      incore (imsk) = .false.
      change (imsk) = .false.
      do i=1,3
        cell (i,imsk)   = 100.0
        cell (i+3,imsk) = 90.0
      end do
      do i=1,maxitm
        know (i,imsk) = .false.
      end do
c
      ierr = 0
c
      return
      end
c
c
c
      integer function whichm (nam,maxset)
c
c ... which dataset does the name "nam" correspond to ?
c
c ... if "*", then return 0, meaning ALL sets
c     if okay, return index of set
c     otherwise:
c     -1 if length = 0
c     -2 if duplicate name
c     -3 if not found
c
      include 'dataman.incl'
c
      integer maxset
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
      do i=1,maxset
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
      subroutine selecm (nam,ierr,maxset)
c
c ... figure out which set(s) are to be selected
c
      include 'dataman.incl'
c
      integer maxset
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
      do i=1,maxset
        select (i) = .false.
        if (incore(i)) nmask = nmask + 1
      end do
c
      if (nmask .le. 0) then
        call errcon ('No sets in memory')
        return
      end if
c
      imsk = whichm(nam,maxset)
c
      if (imsk .lt. 0 .or. imsk .gt. maxset) then
        call errcon ('Invalid set name')
        return
      end if
c
      if (imsk .eq. 0) then
        do i=1,maxset
          select (i) = incore (i)
        end do
        ierr = 0
      else
        if (incore(imsk)) then
          select (imsk) = .true.
          ierr = 0
        else
          call errcon ('Selected set not in memory')
        end if
      end if
c
      return
      end
c
c
c
      subroutine kilhkl (mode,iset,crit,oper,val,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
c
c ... delete certain HKLs if mode = 'K'
c     show them if mode = 'S'
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer morbit(maxhkl,maxset),rfree(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      real val,xval,xdum
c
      integer iset,i,new,ierr,j,leng1,nold
c
      logical lf,ls,lz,lr,lg,ll,lkill
c
      character crit*(*),oper*(*),line*128,mode*(*)
c
code ...
c
      ierr = -1
c
      lkill = (mode(1:1).eq.'k' .or. mode(1:1).eq.'K')
c
      call upcase (crit)
      if (crit(1:3) .ne. 'FOB' .and.
     +    crit(1:3) .ne. 'SIG' .and.
     +    crit(1:3) .ne. 'RES' .and.
     +    crit(1:3) .ne. 'F/S') then
        call errcon ('Invalid criterion')
        return
      end if
      lf = (crit(1:3) .eq. 'FOB')
      ls = (crit(1:3) .eq. 'SIG')
      lr = (crit(1:3) .eq. 'RES')
      lz = (crit(1:3) .eq. 'F/S')
c
      if (lr .and. (.not. know(kreso,iset))) then
        call errcon ('I don''t know the resolution')
        if (.not. know(kcell,iset)) then
          call errcon ('I don''t know the cell constants')
        end if
        return
      end if
c
      call upcase (oper)
      if (oper(1:1) .ne. '<' .and.
     +    oper(1:1) .ne. '>') then
        call errcon ('Invalid operand')
        return
      end if
      ll = (oper(1:1) .eq. '<')
      lg = (oper(1:1) .eq. '>')
c
      if (lkill) then
        call jvalut (' Nr of reflections before :',1,numhkl(iset))
      else
        call jvalut (' Nr of reflections :',1,numhkl(iset))
      end if
c
      nold = numhkl(iset)
c
      write (line,*) crit(1:leng1(crit)),' ',
     +  oper(1:leng1(oper)),' ',val
      call pretty (line)
      if (lkill) then
        call textut (' Kill reflection if :',line)
      else
        call textut (' Show reflection if :',line)
        write (*,6901) 'Reflxn','   H','   K','   L',
     +    'Resoln','   Fobs   ',' Sigma(Fobs)',' F/Sig',
     +    'Test','A/C','Orb'
      end if
c
 6900 format (1x,i6,1x,3i4,1x,f6.2,1x,1p,2e12.4,1x,0p,f6.2,
     +  1x,i4,3x,a1,1x,i3)
 6901 format (1x,a6,1x,3a4,1x,a6,1x,2a12,1x,a6,1x,a4,1x,a3,1x,a3)
c
      new = 0
      do i=1,numhkl(iset)
c
        if (lf) then
          xval = fobs(i,iset)
        else if (ls) then
          xval = sigfob(i,iset)
        else if (lr) then
          xval = reso(i,iset)
        else
          if (sigfob(i,iset) .ge. small) then
            xval = fobs(i,iset) / sigfob(i,iset)
          else
            goto 100
          end if
        end if
c
        if (ll) then
          if (xval .lt. val) goto 100
        else
          if (xval .gt. val) goto 100
        end if
c
c ... keep this one
c
        if (lkill) then
          new = new + 1
          ibuff (new) = i
        end if
        goto 110
c
  100   continue
        if (.not. lkill) then
          new = new + 1
          xdum = -1.0
          if (sigfob(i,1) .ne. 0.0) then
            xdum = fobs(i,1)/sigfob(i,1)
          end if
          write (*,6900) i,(hkl(j,i,1),j=1,3),reso(i,1),
     +      fobs(i,1),sigfob(i,1),xdum,rfree(i,1),
     +      centri(i,1),morbit(i,1)
ccc          write (*,6000) i,(hkl(j,i,iset),j=1,3),
ccc     +      fobs(i,iset),sigfob(i,iset),rfree(i,iset)
        end if
c
  110   continue
      end do
c
 6000 format (' # ',i8,' HKL ',3i6,' Fo, S(Fo) = ',1p,2e12.4,
     +  ' Test ',i1)
c
      if (.not. lkill) then
        call jvalut (' Nr of reflections listed :',1,new)
        return
      end if
c
      if (new .eq. 0) then
        call errcon ('This would kill all reflections; aborted')
        return
      end if
c
      if (new .eq. numhkl(iset)) then
        call prompt (' NOTE --- No reflections killed')
        return
      end if
c
      do i=1,new
        j=ibuff(i)
        if (j .ne. i) then
          hkl(1,i,iset) = hkl(1,j,iset)
          hkl(2,i,iset) = hkl(2,j,iset)
          hkl(3,i,iset) = hkl(3,j,iset)
          fobs(i,iset)   = fobs(j,iset)
          sigfob(i,iset) = sigfob(j,iset)
          reso(i,iset)   = reso(j,iset)
          centri(i,iset) = centri(j,iset)
          morbit(i,iset) = morbit(j,iset)
          rfree(i,iset)  = rfree(j,iset)
        end if
      end do
c
      numhkl (iset) = new
      call jvalut (' Nr of reflections after :',1,new)
      ierr = 0
c
      i = nold - new
      if (i .gt. 0) then
        call jvalut (' WARNING - Reflections killed :',1,i)
      end if
c
      return
      end
c
c
c
      subroutine calres (iset,ierr,
     +                   maxset,maxhkl,reso,hkl)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
c
      real ca,cb,cg,sa,sb,sg,vol,astar,bstar,cstar,car,cbr,cgr
      real ah,bk,cl,sintol,res,xmin,xmax
c
      integer iset,ierr,i
c
code ...
c
      ierr = -1
c
      ca = cos(cell(4,iset)*degtor)
      sa = sin(cell(4,iset)*degtor)
      cb = cos(cell(5,iset)*degtor)
      sb = sin(cell(5,iset)*degtor)
      cg = cos(cell(6,iset)*degtor)
      sg = sin(cell(6,iset)*degtor)
c
      vol = cell(1,iset)*cell(2,iset)*cell(3,iset) *
     +      sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg)
      call rvalut (' Cell volume :',1,vol)
      if (vol .le. small) return
c
      astar = cell(2,iset)*cell(3,iset)*sa/vol
      bstar = cell(1,iset)*cell(3,iset)*sb/vol
      cstar = cell(1,iset)*cell(2,iset)*sg/vol
c
      car = (cb*cg-ca)/(sb*sg)
      cbr = (ca*cg-cb)/(sa*sg)
      cgr = (ca*cb-cg)/(sa*sb)
c
      xmin = 9999.9
      xmax = -999.9
c
      do i=1,numhkl(iset)
        ah = astar*float(hkl(1,i,iset))
        bk = bstar*float(hkl(2,i,iset))
        cl = cstar*float(hkl(3,i,iset))
        sintol = 0.5 * sqrt (ah*ah + bk*bk + cl*cl +
     +    2.0 * (bk*cl*car + ah*cl*cbr + ah*bk*cgr) )
        res = 0.50/sintol
        reso (i,iset) = res
        xmin = min (xmin,res)
        xmax = max (xmax,res)
      end do
c
      call fvalut (' Lowest  resolution :',1,xmax)
      call fvalut (' Highest resolution :',1,xmin)
      ierr = 0
c
      return
      end
c
c
c
      subroutine calcen (iset,ierr,
     +                   maxset,maxhkl,hkl,centri)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      integer hkl(3,maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      integer iset,ierr,i,nc,na
c
      logical iscent
c
code ...
c
      ierr = -1
c
      call jvalut (' Nr of reflections    :',1,numhkl(iset))
      call jvalut (' Nr of unique symmops :',1,nuniq(iset))
c
      na = 0
      nc = 0
c
      do i=1,numhkl(iset)
        centri (i,iset) = 'A'
        if (iscent(hkl(1,i,iset),nuniq(iset),transp(1,1,iset))) then
          centri (i,iset) = 'C'
          nc = nc + 1
cc          print *,hkl(1,i,iset),hkl(2,i,iset),hkl(3,i,iset)
        else
          na = na + 1
        end if
      end do
c
      call jvalut (' Nr of acentric reflections :',1,na)
      call jvalut (' Nr of  centric reflections :',1,nc)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine calorb (iset,ierr,
     +                   maxset,maxhkl,hkl,morbit)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      integer hkl(3,maxhkl,maxset),morbit(maxhkl,maxset)
c
      integer twosym
      parameter (twosym = 2 * maxsym)
c
      integer iset,ierr,i,done(3,twosym),cnt(twosym),om
c
code ...
c
      ierr = -1
c
      call jvalut (' Nr of reflections    :',1,numhkl(iset))
      call jvalut (' Nr of unique symmops :',1,nuniq(iset))
c
      do i=1,twosym
        cnt (i) = 0
      end do
c
      do i=1,numhkl(iset)
        call mulorb (hkl(1,i,iset),nuniq(iset),transp(1,1,iset),
     +    twosym,done,om)
        morbit (i,iset) = om
        cnt (om) = cnt (om) + 1
      end do
c
      do i=twosym,1,-1
        if (cnt(i) .gt. 0) then
          write (*,6000) cnt(i),i
        end if
      end do
c
 6000 format (' There are ',i8,' reflections with O.M. ',i3)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine gerard (iset,jset,pltf1,pltf2,step,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,reso,morbit,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer morbit(maxhkl,maxset)
c
      integer nbinmx, nbinmn
      parameter (nbinmx = 100, nbinmn = 10)
c
      real xdum,xmin,xmax,step,x1,x2,xx(2),yy(2)
      real sumsq1(nbinmx),sumsq2(nbinmx),scale,btemp
      real stol2(nbinmx),lgi2i1(nbinmx),ymin,ymax
      real rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2
c
      integer iset,jset,ierr,i,j,nbins,jerr,npnt,leng1
      integer nrefl1(nbinmx),nrefl2(nbinmx)
      integer nfull1(nbinmx),nfull2(nbinmx)
c
      logical xinter
c
      character pltf1*(*),pltf2*(*),line*80
c
code ...
c
      ierr = -1
c
      call prompt (' Wilson scaling - G. Bricogne')
      call jvalut (' Max nr of bins :',1,nbinmx)
      call jvalut (' Min nr of bins :',1,nbinmn)
      call rvalut (' Bin width      :',1,step)
c
      call textut (' Unchanged Set 1 =',name(iset))
      call textut ('    Scaled Set 2 =',name(jset))
c
      xmin = 9999.9
      xmax = -999.9
c
c ... calculate 4 * (sin(theta)/lambda)**2
c
      do i=1,numhkl(iset)
        xdum = 0.5/reso(i,iset)
        buffer (i) = 4.0 * xdum * xdum
        xmin = min (xmin,buffer(i))
        xmax = max (xmax,buffer(i))
      end do
c
      do i=1,numhkl(jset)
        xdum = 0.5/reso(i,jset)
        j = maxhkl + i
        buffer (j) = 4.0 * xdum * xdum
        xmin = min (xmin,buffer(j))
        xmax = max (xmax,buffer(j))
      end do
c
      call rvalut (' Min 4 * (sin(theta)/lambda)**2 :',1,xmin)
      call rvalut (' Max 4 * (sin(theta)/lambda)**2 :',1,xmax)
c
      call w_scale_2 (
     +  fobs(1,iset),buffer(1),morbit(1,iset),numhkl(iset),
     +  fobs(1,jset),buffer(maxhkl+1),morbit(1,jset),numhkl(jset),
     +  xmin,xmax,step,nbinmx,nbinmn,nbins,
     +  nrefl1,nfull1,nrefl2,nfull2,sumsq1,sumsq2,
     +  scale,btemp)
c
      call jvalut (' Nr of bins used    :',1,nbins)
      call rvalut (' Step size for bins :',1,step)
c
      if (nbins .lt. 2) then
        call errcon ('Less than two bins; aborting')
        ierr = -1
        return
      end if
c
      call prompt (' Applying scale to set 2')
      do i=1,numhkl(jset)
        xdum = 0.5 / reso(i,jset)
        xdum = -btemp*xdum*xdum
        xdum = exp(xdum)
        fobs (i,jset) = scale * fobs (i,jset) * xdum
      end do
c
      write (*,*)
      write (*,6000) 'Bin','NF1','NF2','SSQ1','SSQ2','<I1>',
     +               '<I2>','LOG<I2>/<I1>','4(SIN(T)/L)^2'
      npnt = 0
      xmin = 9.99E19
      xmax = -xmin
      ymin = xmin
      ymax = xmax
      do i=1,nbins
        if (nrefl1(i) .gt. 0 .and.
     +      nrefl2(i) .gt. 0) then
          npnt = npnt + 1
          x1 = sumsq1(i)
          x2 = sumsq2(i)
          sumsq1 (npnt) = x1/float(nfull1(i))
          sumsq2 (npnt) = x2/float(nfull2(i))
          stol2  (npnt) = (float(i)-0.5)*step
          lgi2i1 (npnt) = alog (sumsq2(npnt)/sumsq1(npnt))
          write (*,6010) i,nfull1(i),
     +      nfull2(i),x1,x2,sumsq1(npnt),sumsq2(npnt),
     +      lgi2i1(npnt),stol2(npnt)
c
          xmin = min (xmin, sumsq2(npnt),sumsq1(npnt))
          xmax = max (xmax, sumsq2(npnt),sumsq1(npnt))
          ymin = min (ymin, lgi2i1(npnt))
          ymax = max (ymax, lgi2i1(npnt))
        end if
      end do
      write (*,*)
c
      if (npnt .lt. 2) then
        call errcon ('Not enough filled bins')
        ierr = -1
        return
      end if
c
      call xystat (sumsq1,sumsq2,npnt,
     +  rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2)
      call prompt (' Comparison of <I1> and <I2> :')
      call fvalut (' Correlation coefficient :',1,corr)
      call rvalut (' Scaled R w.r.t. <I1>    :',1,rf3)
      call rvalut (' Scaled R w.r.t. <I2>    :',1,rf4)
      call rvalut (' RMS difference          :',1,rmsd)
      write (*,*)
c
c ... make first plot
c
      call xopxua (1,pltf1,xinter(),jerr)
      if (jerr .ne. 0) then
        call errcon ('Cannot produce first plot')
        goto 100
      end if
c
      call stamp (line)
      write (1,'(a,a)') 'REMARK ',line(1:leng1(line))
      write (1,6100) 'REMARK Wilson plot file #1:'
      write (1,6100) 'REMARK <I> vs (sin(t)/l)**2'
      write (1,6102) 'REMARK Filename = ',pltf1(1:leng1(pltf1))
      write (1,6102) 'REMARK Dataset 1 = ',
     +  name(iset)(1:leng1(name(iset))),' File = ',
     +  file(iset)(1:leng1(file(iset)))
      write (1,6102) 'REMARK Comment 1 = ',
     +  coment(iset)(1:leng1(coment(iset)))
      write (1,6102) 'REMARK Dataset 2 = ',
     +  name(jset)(1:leng1(name(jset))),' File = ',
     +  file(jset)(1:leng1(file(jset)))
      write (1,6102) 'REMARK Comment 2 = ',
     +  coment(jset)(1:leng1(coment(jset)))
      write (1,6100) '!'
      write (1,6100) 'NPOINT',npnt
      write (1,6100) 'XLABEL 4 * (sin(theta)/lambda)**2'
      write (1,6100) 'YLABEL <I>bin = SUM (F**2)/Nfull'
      write (1,6100) 'MRKTYP',1
      write (1,6100) 'COLOUR',4
      write (1,6105) 'XYVIEW',0.0,1.05*stol2(npnt),xmin,xmax
      write (1,6100) 'XVALUE (1p,6e12.4)'
      write (1,6110) (stol2(i),i=1,npnt)
      write (1,6100) '!'
      write (1,6100) 'YVALUE (1p,6e12.4)'
      write (1,6110) (sumsq1(i),i=1,npnt)
      write (1,6100) '!'
      write (1,6100) 'MORE'
      write (1,6100) '!'
      write (1,6100) 'MRKTYP',2
      write (1,6100) 'COLOUR',1 
      write (1,6100) '!'
      write (1,6100) 'YVALUE (1p,6e12.4)'
      write (1,6110) (sumsq2(i),i=1,npnt)
      write (1,6100) '!'
      write (1,6100) 'END'
c
c ... make second plot
c
  100 continue
c
      close (1)
      call xopxua (1,pltf2,xinter(),jerr)
      if (jerr .ne. 0) then
        call errcon ('Cannot produce second plot')
        goto 200
      end if
c
      x1 = 2.0 * btemp
      x2 = -2.0 * alog(scale)
c
      xx(1) = stol2(1)
      yy(1) = x2 + x1*xx(1)
      xx(2) = stol2(npnt)
      yy(2) = x2 + x1*xx(2)
      ymin = min (ymin,yy(1),yy(2))
      ymax = max (ymax,yy(1),yy(2))
c
      call stamp (line)
      write (1,'(a,a)') 'REMARK ',line(1:leng1(line))
      write (1,6100) 'REMARK Wilson plot file #2:'
      write (1,6100) 'REMARK log<I2>/<I1> vs (sin(t)/l)**2'
      write (1,6102) 'REMARK Filename = ',pltf2(1:leng1(pltf2))
      write (1,6102) 'REMARK Dataset 1 = ',
     +  name(iset)(1:leng1(name(iset))),' File = ',
     +  file(iset)(1:leng1(file(iset)))
      write (1,6102) 'REMARK Comment 1 = ',
     +  coment(iset)(1:leng1(coment(iset)))
      write (1,6102) 'REMARK Dataset 2 = ',
     +  name(jset)(1:leng1(name(jset))),' File = ',
     +  file(jset)(1:leng1(file(jset)))
      write (1,6102) 'REMARK Comment 2 = ',
     +  coment(jset)(1:leng1(coment(jset)))
c
      write (line,*) 'Temperature factor = ',btemp
      write (1,'(a,a)') 'REMARK ',line(1:leng1(line))
      write (line,*) 'Scale = ',scale
      write (1,'(a,a)') 'REMARK ',line(1:leng1(line))
c
      write (1,6100) '!'
      write (1,6100) 'NPOINT',npnt
      write (1,6100) 'XLABEL 4 * (sin(theta)/lambda)**2'
      write (1,6100) 'YLABEL log<I2>/<I1>'
      write (1,6100) 'MRKTYP',1
      write (1,6100) 'COLOUR',4 
      write (1,6105) 'XYVIEW',0.0,1.05*stol2(npnt),ymin,ymax
      write (1,6100) 'XVALUE (1p,6e12.4)'
      write (1,6110) (stol2(i),i=1,npnt)
      write (1,6100) '!'
      write (1,6100) 'YVALUE (1p,6e12.4)'
      write (1,6110) (lgi2i1(i),i=1,npnt)
      write (1,6100) '!'
      write (1,6100) 'MORE'
      write (1,6100) '!'
      write (1,6100) 'COLOUR',1
      write (1,6100) 'NPOINT',2
      write (1,6100) 'VALUES (1p,2e12.4)'
      write (1,6110) xx(1),yy(1)
      write (1,6110) xx(2),yy(2)
c
      write (1,6100) '!'
      write (1,6100) 'END'
c
 6000 format (1x,a3,2(1x,a6),5(1x,a12),1x,a13)
 6010 format (1x,i3,2(1x,i6),1p,6(1x,e12.4))
c
 6100 format (a,10(1x,i8))
 6102 format (10a)
 6105 format (a,1p,4(1x,e12.4))
 6110 format (1p,6e12.4)
c
c ... done
c
  200 continue
c
      close (1)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine deltaf (kset,iset,jset,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,
     +                   morbit,centri,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer rfree(maxhkl,maxset),morbit(maxhkl,maxset)
      character*1 centri(maxhkl,maxset)
c
      integer iset,jset,kset,ierr,i,new,icode,m
c
code ...
c
      ierr = -1
c
      call textut (' Delta-F Set 1 =',name(iset))
      call textut ('     and Set 2 =',name(jset))
c
      call prompt (' Encoding reflections of set 2 ...')
      do i=1,numhkl(jset)
        call packin (hkl(1,i,jset),hkl(2,i,jset),
     +    hkl(3,i,jset),0,icode)
        ibuff (3*maxhkl+i) = icode
        ibuff (4*maxhkl+i) = i
      end do
c
      call prompt (' Sorting reflections of set 2 ...')
      call shell (ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),numhkl(jset))
c
      new = 0
      call prompt (' Locating reflections of set 1 in set 2 ...')
c
      do i=1,numhkl(iset)
c
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,icode)
c
        call bindex (icode,ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),
     +    numhkl(jset),m)
c
        if (m .le. 0) goto 190
c
c ... found its mate !
c
        new = new + 1
        hkl (1,new,kset) = hkl(1,i,iset)
        hkl (2,new,kset) = hkl(2,i,iset)
        hkl (3,new,kset) = hkl(3,i,iset)
        fobs (new,kset)   = fobs (i,iset) - fobs (m,jset)
        sigfob (new,kset) = sqrt ( sigfob(i,iset)**2 +
     +                             sigfob(m,jset)**2 )
        reso (new,kset)   = reso (i,iset)
        morbit (new,kset) = morbit (i,iset)
        centri (new,kset) = centri (i,iset)
        rfree (new,kset)  = rfree (i,iset)
c
  190   continue
c
      end do
c
      call jvalut (' HKLs in native     set 1:',1,numhkl(iset))
      call jvalut (' HKLs in derivative set 2:',1,numhkl(jset))
      call jvalut (' HKLs in new nat-der set :',1,new)
c
      if (new .le. 0) return
c
      numhkl (kset) = new
      ierr = 0
c
      return
      end
c
c
c
      subroutine lauegr (kset,iset,laue,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,
     +                   morbit,centri)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
      integer rfree(maxhkl,maxset),morbit(maxhkl,maxset)
      character*1 centri(maxhkl,maxset)
c
      integer twosym
      parameter (twosym = 2 * maxsym)
c
      integer iset,kset,ierr,i,new,laue,j
      integer done(3,twosym),om
c
      logical oklaue
c
code ...
c
      ierr = -1
c
      new = 0
c
      do i=1,numhkl(iset)
        call mulorb (hkl(1,i,iset),nuniq(iset),transp(1,1,iset),
     +    twosym,done,om)
ccc        call ivalut (' OLD :',3,hkl(1,i,iset))
        if (om .gt. 0) then
          do j=1,om
            if (oklaue(done(1,j),laue)) then
              if (new .lt. maxhkl) then
                new = new + 1
                hkl (1,new,kset) = done (1,j)
                hkl (2,new,kset) = done (2,j)
                hkl (3,new,kset) = done (3,j)
                fobs (new,kset) = fobs (i,iset)
                sigfob (new,kset) = sigfob (i,iset)
                reso (new,kset) = reso (i,iset)
                morbit (new,kset) = morbit (i,iset)
                centri (new,kset) = centri (i,iset)
                rfree (new,kset) = rfree (i,iset)
ccc                call ivalut (' ... :',3,hkl(1,new,kset))
              else
                call errcon ('Too many reflections - aborting !!!')
                call jvalut (' Max nr allowed :',1,maxhkl)
                new = 0
                goto 100
              end if
            end if
          end do
        else
          call errcon ('Orbital multiplicity = zero ???')
          call ivalut (' HKL :',3,hkl(1,i,iset))
        end if
      end do
c
  100 continue
c
      call jvalut (' HKLs in old set :',1,numhkl(iset))
      call jvalut (' HKLs in new set :',1,new)
c
      if (new .le. 0) return
c
      call prompt (' HKLs in the new set are UNSORTED !')
      ierr = 0
      numhkl (kset) = new
c
      return
      end
c
c
c
      subroutine sortem (kset,iset,i123,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,
     +                   morbit,centri,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer rfree(maxhkl,maxset),morbit(maxhkl,maxset)
      character*1 centri(maxhkl,maxset)
c
      integer iset,kset,ierr,i,j,i123(3),i1,i2,i3,ioff
c
code ...
c
      ierr = -1
c
      i1 = i123(1)
      i2 = i123(2)
      i3 = i123(3)
c
      ioff = numhkl(iset)
      call prompt (' Encoding reflections of old set ...')
      do i=1,numhkl(iset)
        call packin (hkl(i1,i,iset),hkl(i2,i,iset),
     +    hkl(i3,i,iset),0,ibuff(ioff+i))
        ibuff (i) = i
      end do
c
      call prompt (' Sorting reflections ...')
      call shell (ibuff(ioff+1),ibuff(1),numhkl(iset))
c
      do i=1,numhkl(iset)
        j = ibuff(i)
        hkl (1,i,kset) = hkl (1,j,iset)
        hkl (2,i,kset) = hkl (2,j,iset)
        hkl (3,i,kset) = hkl (3,j,iset)
        fobs (i,kset) = fobs (j,iset)
        sigfob (i,kset) = sigfob (j,iset)
        reso (i,kset) = reso (j,iset)
        morbit (i,kset) = morbit (j,iset)
        centri (i,kset) = centri (j,iset)
        rfree (i,kset) = rfree (j,iset)
      end do
c
      numhkl (kset) = numhkl (iset)
      ierr = 0
c
      return
      end
c
c
c
      subroutine duplic (kset,iset,maxset,maxhkl,
     +                   fobs,sigfob,reso,hkl,rfree,morbit,centri)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
      integer rfree(maxhkl,maxset),morbit(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      integer iset,kset,i,j
c
code ...
c
      do i=1,numhkl(iset)
        j = i
        hkl (1,i,kset) = hkl (1,j,iset)
        hkl (2,i,kset) = hkl (2,j,iset)
        hkl (3,i,kset) = hkl (3,j,iset)
        fobs (i,kset) = fobs (j,iset)
        sigfob (i,kset) = sigfob (j,iset)
        reso (i,kset) = reso (j,iset)
        rfree (i,kset) = rfree (j,iset)
        morbit (i,kset) = morbit (j,iset)
        centri (i,kset) = centri (j,iset)
      end do
c
      numhkl (kset) = numhkl (iset)
c
      know (kdata,kset) = know(kdata,iset)
      know (kreso,kset) = know(kreso,iset)
      know (kcent,kset) = know(kcent,iset)
      know (korbi,kset) = know(korbi,iset)
      know (kfree,kset) = know(kfree,iset)
c
      return
      end
c
c
c
      subroutine cpknow (iset,kset)
c
      include 'dataman.incl'
c
      integer iset,kset,i,j
c
code ...
c
      know (kcell,kset) = know(kcell,iset)
      if (know(kcell,iset)) then
        do i=1,6
          cell (i,kset) = cell (i,iset)
        end do
      end if
c
      know (ksymm,kset) = know(ksymm,iset)
      if (know(ksymm,iset)) then
        nsymop (kset) = nsymop (iset)
        do i=1,nsymop (iset)
          do j=1,12
            symmop (j,i,kset) = symmop (j,i,iset)
          end do
        end do
c
        nuniq (kset) = nuniq (iset)
        do i=1,nuniq (iset)
          uniqso (i,kset) = uniqso (i,iset)
          do j=1,9
            transp (j,i,kset) = transp (j,i,iset)
          end do
        end do
      end if
c
      return
      end
c
c
c
      subroutine twin (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,centri,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),buffer(maxbuf)
      character centri(maxhkl,maxset)*1
c
      real xdum,ave,sdv,xmin,xmax,afx,afa,afc,aix,aia,aic
      real aisqx,aisqa,aisqc
c
      integer iset,i,ni1,ni2,nn1,nn2
c
code ...
c
      ni1 = maxhkl
      ni2 = 2*maxhkl
c
      nn1 = 0
      nn2 = 0
      do i=1,numhkl(iset)
        xdum = abs (fobs(i,iset))
        buffer(i) = xdum
        if (centri(i,iset).eq.'A') then
          nn1 = nn1 + 1
          buffer (ni1+nn1) = xdum
        else
          nn2 = nn2 + 1
          buffer (ni2+nn2) = xdum
        end if
      end do
c
      call jvalut (' Nr of reflections :',1,numhkl(iset))
      call jvalut (' Acentrics         :',1,nn1)
      call jvalut (' Centrics          :',1,nn2)
c
      write (*,*)
      write (*,6000) 'Item  ',  'Average ', 'StDev  ',
     +  'Min   ',   'Max   '
      write (*,6000) '========','=========','=========',
     +  '=========','========='
c
 6000 format (1x,a10,4a15)
 6010 format (1x,a10,1p,4e15.5)
c
      call xstats (buffer(1),numhkl(iset),ave,sdv,xmin,xmax,xdum)
      write (*,6010) ' |F| all ',ave,sdv,xmin,xmax
      afx = ave
c
      if (nn1 .gt. 0) then
        call xstats (buffer(ni1+1),nn1,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' |F| acn ',ave,sdv,xmin,xmax
        afa = ave
      else
        call prompt (' No acentric reflexions')
        afa = 0.0
      end if
c
      if (nn2 .gt. 0) then
        call xstats (buffer(ni2+1),nn2,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' |F| cen ',ave,sdv,xmin,xmax
        afc = ave
      else
        call prompt (' No centric reflexions')
        afc = 0.0
      end if
c
      nn1 = 0
      nn2 = 0
      do i=1,numhkl(iset)
        xdum = fobs(i,iset) * fobs(i,iset)
        buffer(i) = xdum
        if (centri(i,iset).eq.'A') then
          nn1 = nn1 + 1
          buffer (ni1+nn1) = xdum
        else
          nn2 = nn2 + 1
          buffer (ni2+nn2) = xdum
        end if
      end do
c
      call xstats (buffer(1),numhkl(iset),ave,sdv,xmin,xmax,xdum)
      write (*,6010) ' |I| all ',ave,sdv,xmin,xmax
      aix = ave
c
      if (nn1 .gt. 0) then
        call xstats (buffer(ni1+1),nn1,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' |I| acn ',ave,sdv,xmin,xmax
        aia = ave
      else
        call prompt (' No acentric reflexions')
        aia = 0.0
      end if
c
      if (nn2 .gt. 0) then
        call xstats (buffer(ni2+1),nn2,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' |I| cen ',ave,sdv,xmin,xmax
        aic = ave
      else
        call prompt (' No centric reflexions')
        aic = 0.0
      end if
c
      nn1 = 0
      nn2 = 0
      do i=1,numhkl(iset)
        xdum = fobs(i,iset)
        xdum = xdum*xdum*xdum*xdum
        buffer(i) = xdum
        if (centri(i,iset).eq.'A') then
          nn1 = nn1 + 1
          buffer (ni1+nn1) = xdum
        else
          nn2 = nn2 + 1
          buffer (ni2+nn2) = xdum
        end if
      end do
c
      call xstats (buffer(1),numhkl(iset),ave,sdv,xmin,xmax,xdum)
      write (*,6010) ' I*I all ',ave,sdv,xmin,xmax
      aisqx = ave
c
      if (nn1 .gt. 0) then
        call xstats (buffer(ni1+1),nn1,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' I*I acn ',ave,sdv,xmin,xmax
        aisqa = ave
      else
        call prompt (' No acentric reflexions')
        aisqa = 0.0
      end if
c
      if (nn2 .gt. 0) then
        call xstats (buffer(ni2+1),nn2,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' I*I cen ',ave,sdv,xmin,xmax
        aisqc = ave
      else
        call prompt (' No centric reflexions')
        aisqc = 0.0
      end if
c
      write (*,*)
      if (aic .ne. 0.0)
     +  call fvalut (' <I**2>/<I>**2 for CENTRO       :',
     +    1,(aisqc/(aic**2)))
      if (aia .ne. 0.0)
     +  call fvalut (' <I**2>/<I>**2 for NON-CS       :',
     +  1,(aisqa/(aia**2)))
      if (aix .ne. 0.0)
     +  call fvalut (' <I**2>/<I>**2 for ALL          :',
     +  1,(aisqx/(aix**2)))
c
      write (*,*)
      if (afc .ne. 0.0)
     +  call fvalut (' <F**2>/<F>**2 for CENTRO       :',
     +  1,(aic/(afc**2)))
      if (afa .ne. 0.0)
     +  call fvalut (' <F**2>/<F>**2 for NON-CS       :',
     +  1,(aia/(afa**2)))
      if (afx .ne. 0.0)
     +  call fvalut (' <F**2>/<F>**2 for ALL          :',
     +  1,(aix/(afx**2)))
c
      write (*,*)
      if (aic .ne. 0.0)
     +  call fvalut (' Wilson ratio for CENTRO        :',
     +  1,(afc**2/aic))
      if (aia .ne. 0.0)
     +  call fvalut (' Wilson ratio for NON-CS        :',
     +  1,(afa**2/aia))
      if (aix .ne. 0.0)
     +  call fvalut (' Wilson ratio for ALL           :',
     +  1,(afx**2/aix))
c
      write (*,*)
      call fvalut (' Wilson ratio non-twinned CENTRO :',1,0.637)
      call fvalut (' Wilson ratio non-twinned NON-CS :',1,0.785)
      call fvalut (' Wilson ratio 1:1-twinned NON-CS :',1,0.885)
c
      return
      end
c
c
c
      subroutine gemini (iset,psf1,psf2,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,centri,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),buffer(maxbuf)
      character centri(maxhkl,maxset)*1
c
      integer nnz, nna
      parameter (nnz = 10, nna = 5)
c
      real pint
      integer npint
      parameter (pint = 0.05, npint = 2 + (1/pint))
c
      integer f1
      parameter (f1=11)
c
      real xdum,ave,sdv,xmin,xmax,aix,aia,aic,x1,x2,x3,x4,x5
      real z(nnz),a(nna),nzanon(nnz,nna),nzacen(nnz,nna)
      real x(nnz),y(nnz),dumx,dumc,dumn,dum,ba,best
      real x6,x7,x8,x9
c
      integer iset,i,j,ni1,ni2,nn1,nn2
c
      character psf1*(*),psf2*(*),myline*120
c
      data z /0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      data a /0.0,0.1,0.2,0.3,0.5/
c
      data nzanon/
     +  0.10,0.18,0.26,0.33,0.39,0.45,0.50,0.55,0.59,0.63,
     +  0.04,0.12,0.20,0.28,0.36,0.42,0.48,0.54,0.59,0.63,
     +  0.03,0.08,0.16,0.24,0.31,0.39,0.45,0.52,0.57,0.62,
     +  0.02,0.07,0.14,0.21,0.28,0.36,0.43,0.49,0.55,0.61,
     +  0.02,0.06,0.12,0.19,0.26,0.34,0.41,0.48,0.54,0.59/
c
      data nzacen/
     +  0.25,0.35,0.42,0.47,0.52,0.56,0.60,0.63,0.66,0.68,
     +  0.15,0.26,0.35,0.42,0.48,0.53,0.58,0.62,0.65,0.68,
     +  0.12,0.22,0.30,0.38,0.44,0.50,0.55,0.59,0.63,0.66,
     +  0.10,0.19,0.28,0.35,0.41,0.47,0.52,0.57,0.61,0.65,
     +  0.10,0.18,0.26,0.33,0.39,0.45,0.50,0.55,0.59,0.63/
c
code ...
c
      write (*,'(20(/,1x,a))')
     +  ' REFERENCES:',
     +  ' ',
     +  ' (1) D.C. Rees, "The Influence of Twinning by Merohedry on',
     +  ' Intensity Statistics", Acta Cryst A36, 578-581 (1980)',
     +  ' (2) E. Stanley, "The Identification of Twins from',
     +  ' Intensity Statistics", J Appl Cryst 5, 191-194 (1972)'
c
      call textut (' Gemini Set =',name(iset))
c
      call fvalut (' Z sampled at :',nnz,z)
      call fvalut (' ALPHA sampled at :',nna,a)
      write (*,*)
c
      ni1 = maxhkl
      ni2 = 2*maxhkl
c
      nn1 = 0
      nn2 = 0
      do i=1,numhkl(iset)
        xdum = fobs(i,iset) * fobs(i,iset)
        buffer(i) = xdum
        if (centri(i,iset).eq.'A') then
          nn1 = nn1 + 1
          buffer (ni1+nn1) = xdum
        else
          nn2 = nn2 + 1
          buffer (ni2+nn2) = xdum
        end if
      end do
c
      call jvalut (' Nr of reflections :',1,numhkl(iset))
      call jvalut (' Acentrics         :',1,nn1)
      call jvalut (' Centrics          :',1,nn2)
c
      write (*,*)
      write (*,6000) 'Item  ',  'Average ', 'StDev  ',
     +  'Min   ',   'Max   '
      write (*,6000) '========','=========','=========',
     +  '=========','========='
c
 6000 format (1x,a10,4a15)
 6010 format (1x,a10,1p,4e15.5)
c
      call xstats (buffer(1),numhkl(iset),ave,sdv,xmin,xmax,xdum)
      write (*,6010) ' |I| all ',ave,sdv,xmin,xmax
      aix = ave
c
      if (nn1 .gt. 0) then
        call xstats (buffer(ni1+1),nn1,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' |I| acn ',ave,sdv,xmin,xmax
        aia = ave
      else
        call prompt (' No acentric reflexions')
        aia = 1.0
      end if
c
      if (nn2 .gt. 0) then
        call xstats (buffer(ni2+1),nn2,ave,sdv,xmin,xmax,xdum)
        write (*,6010) ' |I| cen ',ave,sdv,xmin,xmax
        aic = ave
      else
        call prompt (' No centric reflexions')
        aic = 1.0
      end if
c
      dumx = 1.0 / aix
      dumc = 1.0 / aic
      dumn = 1.0 / aia
c
      do i=1,numhkl(iset)
        buffer (i) = dumx * buffer (i)
      end do
c
      do i=1,nn1
        buffer (ni1+i) = dumn * buffer (ni1+i)
      end do
c
      do i=1,nn2
        buffer (ni2+i) = dumc * buffer (ni2+i)
      end do
c
      call xstats (buffer(1),numhkl(iset),ave,sdv,xmin,xmax,xdum)
      write (*,6010) '  z  all ',ave,sdv,xmin,xmax
c
      if (nn1 .gt. 0) then
        call xstats (buffer(ni1+1),nn1,ave,sdv,xmin,xmax,xdum)
        write (*,6010) '  z  acn ',ave,sdv,xmin,xmax
      else
        call prompt (' No acentric reflexions')
      end if
c
      if (nn2 .gt. 0) then
        call xstats (buffer(ni2+1),nn2,ave,sdv,xmin,xmax,xdum)
        write (*,6010) '  z  cen ',ave,sdv,xmin,xmax
      else
        call prompt (' No centric reflexions')
      end if
c
c ... initialise counters
c
      do i=1,nnz
        x (i) = 0.0
        y (i) = 0.0
      end do
c
c ... fill bins
c
      if (nn2 .gt. 0) then
      do i=1,nn2
        do j=1,nnz
          if (buffer(ni2+i) .lt. z(j)) then
            y (j) = y (j) + 1
            goto 1297
          end if
        end do
 1297   continue
      end do
      end if
c
      if (nn1 .gt. 0) then
      do i=1,nn1
        do j=1,nnz
          if (buffer(ni1+i) .lt. z(j)) then
            x (j) = x (j) + 1
            goto 1298
          end if
        end do
 1298   continue
      end do
      end if
c
c ... make distribution cumulative
c
      do i=2,nnz
        if (nn1 .gt. 0) x (i) = x (i) + x (i-1)
        if (nn2 .gt. 0) y (i) = y (i) + y (i-1)
      end do
c
c ... compute fractions
c
      dumc = 1.0
      dumn = 1.0
      if (nn2 .gt. 0)
     +  dumc = 1.0 / float(nn2)
      if (nn1 .gt. 0)
     +  dumn = 1.0 / float(nn1)
      do i=1,nnz
        if (nn1 .gt. 0) x (i) = x (i) * dumn
        if (nn2 .gt. 0) y (i) = y (i) * dumc
      end do
c
      if (nn1 .gt. 0) then
        write (*,*)
        call fvalut (' DIST NON-CS :',nnz,x)
c
        best = 999.9
        ba   = 0.0
        do i=1,nna
          call xystat (x,nzanon(1,i),nnz,
     +      x1,x2,x3,x4,x5,x6,x7,x8,x9)
          write (*,'(a,f4.2,a,f6.3,a,f6.3)')
     +      ' For ALPHA = ',a(i),' RMSD to curve = ',x1,
     +      ' and CORR COEFF = ',x3
          if (x1*(1.0-x3) .lt. best) then
            best = x1*(1.0-x3)
            ba   = a(i)
          end if
        end do
        call fvalut (' Most likely twin fraction :',1,ba)
      else
        write (*,*)
        call prompt (' No acentric reflections !')
      end if
c
      if (nn2 .gt. 0) then
        write (*,*)
        call fvalut (' DIST CENTRO :',nnz,y)
c
        best = 999.9
        ba   = 0.0
        do i=1,nna
          call xystat (y,nzanon(1,i),nnz,
     +      x1,x2,x3,x4,x5,x6,x7,x8,x9)
          write (*,'(a,f4.2,a,f6.3,a,f6.3)')
     +      ' For ALPHA = ',a(i),' RMSD to curve = ',x1,
     +      ' and CORR COEFF = ',x3
          if (x1*(1.0-x3) .lt. best) then
            best = x1*(1.0-x3)
            ba   = a(i)
          end if
        end do
        call fvalut (' Most likely twin fraction :',1,ba)
      else
        write (*,*)
        call prompt (' No centric reflections !')
      end if
c
      write (*,*)
c
c ... make PostScript plot 1
c
      if (nn1 .gt. 0) then
      call xps_init ()
      call xps_open (f1,psf1,'DATAMAN')
c
      call xps_legend ('DATAMAN Gemini plot')
      call xps_legend ('1N(z,alpha) non-centrosymmetric vs. Z')
      myline = 'Dataset '//name(iset)//' File '//file(iset)
      call xps_legend (myline)
      myline = 'Comment '//coment(iset)
      call xps_legend (myline)
      write (myline,*) 'Most likely twin fraction ',ba
      call xps_legend (myline)
c
      call xps_dash ()
      call xps_colour (1)
      call xps_scale (0.0,1.0,0.0,1.0)
      call xps_label ('Z','1N(z,alpha) non-centrosymmetric')
c
      call xps_ps_comment ('Curves')
      do i=1,nna
        call xps_move (0.0,0.0)
        do j=1,nnz
          call xps_draw (z(j),nzanon(j,i))
        end do
      end do
c
      call xps_ps_comment ('Points')
      dum = 0.01
      call xps_symbol (1,-dum,dum,-dum,dum)
      call xps_colour (4)
      do j=1,nnz
        call xps_symbol (1,z(j)-dum,z(j)+dum,
     +    x(j)-dum,x(j)+dum)
      end do
      call xps_move (0.0,0.0)
c      do j=1,npint
c        call xps_draw (qnz(j),qnx(j))
c      end do
c
      call xps_close ()
      end if
c
c ... make PostScript plot 2
c
      if (nn2 .gt. 0) then
      call xps_init ()
      call xps_open (f1,psf2,'DATAMAN')
c
      call xps_legend ('DATAMAN Gemini plot')
      call xps_legend ('-1N(z,alpha) centrosymmetric vs. Z')
      myline = 'Dataset '//name(iset)//' File '//file(iset)
      call xps_legend (myline)
      myline = 'Comment '//coment(iset)
      call xps_legend (myline)
      write (myline,*) 'Most likely twin fraction ',ba
      call xps_legend (myline)
c
      call xps_dash ()
      call xps_colour (1)
      call xps_scale (0.0,1.0,0.0,1.0)
      call xps_label ('Z','-1N(z,alpha) centrosymmetric')
c
      call xps_ps_comment ('Curves')
      do i=1,nna
        call xps_move (0.0,0.0)
        do j=1,nnz
          call xps_draw (z(j),nzacen(j,i))
        end do
      end do
c
      call xps_ps_comment ('Points')
      dum = 0.01
      call xps_symbol (1,-dum,dum,-dum,dum)
      call xps_colour (4)
      do j=1,nnz
        call xps_symbol (1,z(j)-dum,z(j)+dum,
     +    y(j)-dum,y(j)+dum)
      end do
      call xps_move (0.0,0.0)
c      do j=1,npint
c        call xps_draw (qnz(j),qny(j))
c      end do
c
      call xps_close ()
      end if
c
      return
      end
c
c
c
      subroutine hemisp (kset,iset,maxset,maxhkl,resolu,laue,
     +                   fobs,sigfob,reso,hkl,rfree,type,ierr)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
      integer rfree(maxhkl,maxset)
c
      real ca,cb,cg,sa,sb,sg,vol,astar,bstar,cstar,car,cbr,cgr
      real ah,bk,cl,sintol,res,xmin,xmax,resolu
c
      integer myhkl(3)
      integer iset,kset,i,ierr,ih,ik,il,new,hmax,kmax,lmax,laue
c
      logical oklaue
c
      character type*(*)
c
      equivalence (ih,myhkl(1))
      equivalence (ik,myhkl(2))
      equivalence (il,myhkl(3))
c
code ...
c
      ierr = -1
c
      call textut (' Generate   :',type)
      call fvalut (' Resolution :',1,resolu)
      call ivalut (' Laue group :',1,laue)
      if (resolu .lt. 0.1) then
        call errcon ('Silly resolution limit')
        return
      end if
c
      do i=1,6
        cell (i,kset) = cell (i,iset)
      end do
      call fvalut (' Cell :',6,cell(1,kset))
c
      ca = cos(cell(4,kset)*degtor)
      sa = sin(cell(4,kset)*degtor)
      cb = cos(cell(5,kset)*degtor)
      sb = sin(cell(5,kset)*degtor)
      cg = cos(cell(6,kset)*degtor)
      sg = sin(cell(6,kset)*degtor)
c
      vol = cell(1,kset)*cell(2,kset)*cell(3,kset) *
     +      sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg)
      call rvalut (' Cell volume :',1,vol)
      if (vol .le. small) then
        call errcon ('Silly cell volume')
        return
      end if
c
      astar = cell(2,kset)*cell(3,kset)*sa/vol
      bstar = cell(1,kset)*cell(3,kset)*sb/vol
      cstar = cell(1,kset)*cell(2,kset)*sg/vol
c
      car = (cb*cg-ca)/(sb*sg)
      cbr = (ca*cg-cb)/(sa*sg)
      cgr = (ca*cb-cg)/(sa*sb)
c
      hmax = 1 + nint(cell(1,kset)/resolu)
      kmax = 1 + nint(cell(2,kset)/resolu)
      lmax = 1 + nint(cell(3,kset)/resolu)
c
      call jvalut (' Hmax :',1,hmax)
      call jvalut (' Kmax :',1,kmax)
      call jvalut (' Lmax :',1,lmax)
c
      xmin = 9999.9
      xmax = -999.9
c
      new = 0
      do ih=-hmax,hmax
        ah = astar*float(ih)
        do ik=-kmax,kmax
          bk = bstar*float(ik)
          do il=0,lmax
c
           if (laue .eq. 3) then
c
c ... use: HKL: L>=0
c          HK0: H>=0
c          0K0: K>=0
c
              if (il.eq.0) then
                if (ih .lt. 0) goto 90
                if (ih .eq. 0) then
                  if (ik .lt. 0) goto 90
                end if
              end if
            else
              if (.not. oklaue(myhkl,laue)) goto 90
            end if
c
            cl = cstar*float(il)
            sintol = 0.5 * sqrt (ah*ah + bk*bk + cl*cl +
     +               2.0 * (bk*cl*car + ah*cl*cbr + ah*bk*cgr) )
c
c ... skip F(000)
c
            if (sintol .eq. 0.0) goto 90
c
            res = 0.50/sintol
            if (res .ge. resolu) then
              new = new + 1
              if (new .gt. maxhkl) then
                new = new - 1
                call jvalut (' Max nr of reflections :',1,maxhkl)
                call errcon ('Too many reflections; rest skipped')
                goto 1000
              end if
              hkl (1,new,kset) = ih
              hkl (2,new,kset) = ik
              hkl (3,new,kset) = il
              fobs (new,kset) = 1.0
              reso (new,kset) = res
              rfree (new,kset) = 0
              sigfob (new,kset) = 0.0
              if (res .lt. xmin) xmin = res
              if (res .gt. xmax) xmax = res
            else
              goto 100
            end if
   90       continue
          end do
  100     continue
        end do
      end do
c
 1000 continue
      call jvalut (' Nr of reflections generated :',1,new)
      if (new .le. 0) return
c
      numhkl (kset) = new
c
      call fvalut (' Lowest  resolution :',1,xmax)
      call fvalut (' Highest resolution :',1,xmin)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine locint (iset,psf1,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,centri,
     +                   ibuff,jbuff,kbuff,buffer)
c
c ... Produce Padilla-Yeates statistics and plot
c
c ... Gerard J Kleywegt 2004
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),buffer(maxhkl)
      integer hkl(3,maxhkl,maxset)
      integer ibuff(maxhkl),jbuff(maxhkl),kbuff(maxhkl)
      character centri(maxhkl,maxset)*1
c
      real pint,zeroto
      integer npint
      parameter (pint = 0.05, zeroto=1.0E-6)
      parameter (npint = (1/pint))
c
      integer f1
      parameter (f1=11)
c
      real nlabs(-1:npint)
      real sumabs,sumsqr,bigl,aveabs,avesqr,xvalue,yvalue,dum
c
      integer bincnt(0:npint)
      integer iset,i,j,ibin,icode,ndelta,npair,nn1,nn2,k,l,m
c
      character psf1*(*),myline*120
c
      data ndelta /2/
c
code ...
c
      write (*,'(20(/,1x,a))')
     +  'REFERENCE:',
     +  'JE Padilla & TO Yeates, Acta Cryst D59, 1124 (2003).'
c
      call textut (' PY stats for Set =',name(iset))
      write (*,*)
c
      call prompt (' Encoding reflections ...')
      nn1 = 0
      nn2 = 0
      do i=1,numhkl(iset)
        if (centri(i,iset) .eq. 'A') then
          nn1 = nn1 + 1
        else
          nn2 = nn2 + 1
        end if
c
c ... generate a unique integer hash code for every reflection
c
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +               hkl(3,i,iset),0,ibuff(i))
        kbuff (i) = ibuff (i)
        jbuff (i) = i
c
c ... approximate intensity as I = F^2
c
        buffer (i) = fobs (i,iset) * fobs (i,iset)
      end do
c
      call jvalut (' Nr of reflections :',1,numhkl(iset))
      call jvalut (' Acentrics         :',1,nn1)
      call jvalut (' Centrics          :',1,nn2)
      write (*,*)
c
      if (nn1 .lt. 100) then
        call errcon ('Fewer than 100 acentric reflections')
        return
      end if
c
c ... sort the reflections by their hash code (this is
c     necessary so we can use bisectioning during the
c     search for neighbouring reflections later on which
c     in turn is *critical* for the performance of the
c     algorithm !!!)
c
      call prompt (' Sorting reflections ...')
c      call jvalut (' ICODE BEFORE :',10,kbuff)
c      call jvalut (' INDEX BEFORE :',10,jbuff)
      call shell (kbuff,jbuff,numhkl(iset))
c      call jvalut (' ICODE AFTER  :',10,kbuff)
c      call jvalut (' INDEX AFTER :',10,jbuff)
c
      call prompt (' Calculating local intensity statistics ...')
      npair = 0
      do i=0,npint
        bincnt (i) = 0
      end do
      sumabs = 0.0
      sumsqr = 0.0
c
c ... loop over the reflections
c
      do i=1,numhkl(iset)-1
c
c ... skip if centric
c
        if (centri(i,iset) .ne. 'A') goto 199
c
c ... skip if intensity close to zero
c
        if (buffer(i) .le. zeroto) goto 199
c
c ... loop over the "nearby" reflections (in terms of "h,k,l distance")
c
        do j=-ndelta,ndelta,ndelta
          do k=-ndelta,ndelta,ndelta
            do l=-ndelta,ndelta,ndelta
c
c ... skip the reflection itself
c
              if (j.eq.0 .and. k.eq.0 .and. l.eq.0) goto 190
c
c ... generate hash code for this reflection
c
              call packin (hkl(1,i,iset)+j,hkl(2,i,iset)+k,
     +                     hkl(3,i,iset)+l,0,icode)
c
c ... is the reflection in our list of reflections ???
c     we use bisectioning to determine if this is the case
c     the "naive" approach of looping over all reflections,
c     even when using hash codes, will take N/2 steps, whereas
c     bisectioning takes 2LOG(N) - the difference is enormous
c     for large numbers of reflections and is the difference
c     between this routine taking a couple of seconds or many
c     minutes !!! the entire loop is now O(N.LOG(N)) instead
c     of O(N^2)
c
              call bindex (icode,kbuff,jbuff,numhkl(iset),m)
c
c ... if m < 0 the reflection doesn't occur in our list of reflections
c
              if (m .le. 0) goto 190
c
c ... consider each pair only once (unless the other guy is centric,
c     otherwise the pair would never get used)
c
              if (centri(m,iset) .eq. 'A' .and.
     +            m .lt. i) goto 190
c
c ... add new pair
c
              npair = npair + 1
              bigl = (buffer (i) - buffer (m)) /
     +               (buffer (i) + buffer (m))
              if (bigl .lt. 0.0) bigl = -bigl
c
              sumabs = sumabs + bigl
              sumsqr = sumsqr + (bigl * bigl)
              ibin = int (bigl / pint)
              if (ibin .gt. npint - 1) ibin = npint - 1
ccc              print *,i,m,bigl,ibin,buffer(i),buffer(m)
              bincnt (ibin) = bincnt (ibin) + 1
c
c              write (98,'(2i12,6i6)') i,m,
c     +          hkl(1,i,iset),hkl(1,m,iset),
c     +          hkl(2,i,iset),hkl(2,m,iset),
c     +          hkl(3,i,iset),hkl(3,m,iset)
c
  190         continue
            end do
          end do
        end do
  199   continue
      end do
c
c ... note: originally the inner loop contained a loop over all
c           reflections to look for ICODE - this made the operation
c           of order N^2, whereas now it is of order N.log(N).
c           tested that the numerical results are identical
c
c              if (j.eq.0 .and. k.eq.0 .and. l.eq.0) goto 190
c              call packin (hkl(1,i,iset)+j,hkl(2,i,iset)+k,
c     +                     hkl(3,i,iset)+l,0,icode)
c              do m=i+1,numhkl(iset)
c                if (ibuff(m) .eq. icode) then
c                  npair = npair + 1
c                  bigl = (buffer (i) - buffer (m)) /
c     +                   (buffer (i) + buffer (m))
c                  if (bigl .lt. 0.0) bigl = -bigl
c                  sumabs = sumabs + bigl
c                  sumsqr = sumsqr + (bigl * bigl)
c                  ibin = int (bigl / pint)
ccc                  print *,i,m,bigl,ibin,buffer(i),buffer(m)
c                  bincnt (ibin) = bincnt (ibin) + 1
c                  write (97,'(2i12,6i6)') i,m,
c     +              hkl(1,i,iset),hkl(1,m,iset),
c     +              hkl(2,i,iset),hkl(2,m,iset),
c     +              hkl(3,i,iset),hkl(3,m,iset)
c                  goto 190
c                end if
c              end do
c  190         continue
c
c ... print overall statistics
c
      write (*,*)
      aveabs = sumabs / float(npair)
      avesqr = sumsqr / float(npair)
      write (*,6000) aveabs,0.5,0.375
      write (*,6010) avesqr,1.0/3.0,0.2
      write (*,6020) npair,sumabs,sumsqr
      write (*,*)
c
 6000 format (' <|L|> = ',f6.3,' Untwinned = ',f6.3,
     +  ' Perfectly twinned = ',f6.3)
 6010 format (' <L^2> = ',f6.3,' Untwinned = ',f6.3,
     +  ' Perfectly twinned = ',f6.3)
 6020 format (' Npair = ',i12,' SUM |L| = ',f12.3,
     +  ' SUM L^2 = ',f12.3)
c
c ... generate data for plot
c
      nlabs (-1) = 0.0
      do i=0,npint
        nlabs (i) = nlabs (i-1) + float(bincnt(i))
      end do
      do i=0,npint
        nlabs (i) = nlabs (i) / float(npair)
        write (*,6100) float(i)*pint,float(i+1)*pint,
     +    bincnt(i),nlabs(i)
        if (float(i+1)*pint .ge. 1.0) goto 299
      end do
  299 continue
c
 6100 format (' Bin ',f6.3,' - ',f6.3,' Count = ',i12,
     +  ' N(|L|) = ',f6.3)
c
c ... make PostScript plot
c
      call xps_init ()
      call xps_open (f1,psf1,'DATAMAN')
c
      call xps_legend ('DATAMAN Local Intensity Statistics plot')
      call xps_legend ('Cumulative N(|L|) vs. |L| (acentrics)')
      myline = 'Dataset '//name(iset)//' File '//file(iset)
      call xps_legend (myline)
      myline = 'Comment '//coment(iset)
      call xps_legend (myline)
      write (myline,6000) aveabs,0.5,0.375
      call xps_legend (myline)
      write (myline,6010) avesqr,1.0/3.0,0.2
      call xps_legend (myline)
c
      call xps_dash ()
      call xps_colour (1)
      call xps_scale (0.0,1.0,0.0,1.0)
      call xps_label ('|L|','N(|L|) acentrics')
c
      call xps_ps_comment ('Curves')
c
      call xps_move (0.0,0.0)
      call xps_draw (1.0,1.0)
c
      call xps_move (0.0,0.0)
      do i=0,npint
        xvalue = float(i+1) * pint
        yvalue = 0.5*xvalue*(3.0-xvalue*xvalue)
        if (xvalue .le. 1.001) call xps_draw (xvalue,yvalue)
      end do
c
      call xps_ps_comment ('Points')
      dum = 0.01
      call xps_symbol (1,-dum,dum,-dum,dum)
      call xps_colour (4)
c
      do i=0,npint
        xvalue = float(i+1) * pint
        if (xvalue .le. 1.001) then
          call xps_symbol (1,xvalue-dum,xvalue+dum,
     +      nlabs(i)-dum,nlabs(i)+dum)
        end if
      end do
c
      call xps_move (0.0,0.0)
      do i=0,npint
        xvalue = float(i+1) * pint
        if (xvalue .le. 1.001) then
          call xps_draw (xvalue,nlabs(i))
        end if
      end do
c
      call xps_legend (
     +  'See: JE Padilla & TO Yeates, Acta Cryst D59, 1124 (2003).')
      call xps_legend ('Red line = theoretical untwinned')
      call xps_legend ('Red curve = theoretical perfectly twinned')
      call xps_legend ('Blue curve + points = observed')
c
      call xps_close ()
      call prompt (' PostScript plot file written')
c
      return
      end
