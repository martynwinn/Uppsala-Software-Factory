      program lsqman
c
c ... LSQMAN - Least-squares superposition and improvement
c
c ... Gerard Kleywegt @ 931007
c
c ... 0.1 @ 931007 - first version
c
      character*12 prognm,vers
      parameter (prognm = 'LSQMAN', vers = '081126/9.7.9')
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr
      integer fmalloc
#endif
c
      integer maxbuf,nb
code ...
c
      call gkinit (prognm,vers)
c
      maxbuf = 1000000
      call extint ('GKBUFFER',maxbuf)
      maxbuf = max ( maxbuf , 10000 )
      call jvalut (' Allocate buffer arrays of size  :',1,maxbuf)
c
      nb = 4 * maxbuf
      iaptr = fmalloc (nb)
      ibptr = fmalloc (nb)
      icptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call dolsq (maxbuf,%val(iaptr),%val(ibptr),%val(icptr),
     +  prognm,vers)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine dolsq (maxbuf,pl1buf,pl2buf,pl3buf,prognm,vers)
c
      include 'lsqman.incl'
c
      integer maxbuf
      real pl1buf(maxbuf),pl2buf(maxbuf),pl3buf(maxbuf)
      character*12 prognm,vers
c
      integer maxopt,maxaat,maxncs,maxsis,maxnuc
ccc      parameter (maxopt = max (25, maxtyp+2))
      parameter (maxopt = 200, maxaat=20, maxnuc=6)
      parameter (maxncs = 200, maxsis=1000)
c
      real xmat(maxchn,maxchn),cmpmat(maxaat,maxaat),rmsrms(maxchn)
      real dumncs(12,maxncs),duminv(12,maxncs),sirmsd(maxsis)
      real rtdum (12),fra2ca(3,3),ca2fra(3,3),dummy(3),rtbest(12)
      real total,user,sys,xdum,det3,xbest,qqq,cacut,delta,rmslow
      real vrdist,rgbbg1,rgbbg2,rgbbg3,rgbfg1,rgbfg2,rgbfg3
      real rr,gg,bb,judis,juphi,juchi,tarcut,sistart,siend,sistep
      real qqq1,qqq2,qqq3,qqq4,qqq5,qqq6,xmax,ymax,ranamp,q3,qqq9
      real castrt,castop,castep,cuttor,cutcv,cutbf,omdist,gappen
      real nwcut,cutmp,gaploc,radius,distce,qqq7,dmax,xold,qqq8
      real minfra,mamcdf,qqq10,qqq11,qqq12,qqq13,xtr
c
      integer sinali(maxsis),iseq1(maxres),iseq2(maxres)
      integer i,j,k,length,nempty,nopt,ierr,imol,frleng,frstep,nmlow
      integer whichm,jmol,idum,ibest,nat1,nat2,i1low,j1low,minmat
      integer munit,leng1,iii,ichn,jchn,i1,i2,j1,j2,ii,ivrml,ncsdum
      integer kmol,kchn,nmorph,k1,k2,jdum,iptr,jptr,jj,ncolen,inew
      integer ires,nplot,m,l,nmat,jres,kk,nout,maxdp,iter,nold,f2step
      integer iremin,ireshi,ntries,mtries,nimped,ltarg,nt,mtr1,mtr2
c
      logical xinter,linter
      logical ldone,linit,lvrml,lecho,lbsame,lrtf
c
      character typ3lc(maxaat)*3,typ1lc(maxaat)*1
      character nuc1lc(0:maxnuc)*1
      character seq1(maxres)*1,seq2(maxres)*1
      character line*1024,optpar(maxopt)*1024,sdum*120,cavrml*4
      character parnam*40,partyp*1,parfmt*40,pro*12,prev*10
      character zone1*80,zone2*80,zone3*2,zone4*2,inimac*128
      character vrmlbg*25,vrmldc*25,vrfile*128,omcent*4,omconn*20
      character chncol(maxchn)*20,nwmode*10
c
c ... DATA
c
      data typ3lc /'ALA','ARG','ASN','ASP','CYS','GLU','GLN',
     +             'GLY','HIS','ILE','LEU','LYS','MET','PHE',
     +             'PRO','SER','THR','TRP','TYR','VAL'/
c
      data typ1lc /'A',  'R',  'N',  'D',  'C',  'E',  'Q',
     +             'G',  'H',  'I',  'L',  'K',  'M',  'F',
     +             'P',  'S',  'T',  'W',  'Y',  'V'/
c
      data nuc1lc /'?', 'a', 'g', 'c', 't', 'u', '?'/
c
code ...
c
c ... initialise history
c
      call dohist ('*INIT*',ldone)
c
      write (*,*)
      call jvalut (' Max nr of molecules             :',1,maxmol)
      call jvalut (' Max nr of residues per molecule :',1,maxres)
      call jvalut (' Max nr of atoms per molecule    :',1,maxatm)
      call jvalut (' Max nr of atom types            :',1,maxtyp)
      call jvalut (' Max nr of chains/models per mol :',1,maxchn)
      write (*,*)
c
      if (maxopt .lt. (maxtyp+2) ) then
        call errstp ('MAXOPT too small; ask Gerard to recompile !')
        goto 9000
      end if
c
      call inicmp (cmpmat,maxaat)
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
c
      lecho = .false.
c
      ounit = -1
      omol = -1
      punit = -1
c
      nmorph = 10
      cutwat = 1.5
      cuthis = 5.0
      cacut  = 3.5
      delta  = 0.01
      bcutlo = -999.999
      bcuthi = 9999.999
      frleng = 50
      frstep = 25
      f2step = 1
      minmat = 100
      minfra = 0.9
      judis = 1.0
      juphi = 10.0
      juchi = 30.0
      tarcut = 999.0
      sistart = 0.1
      siend = 10.0
      sistep = 0.1
      ranamp = 3.0
      cuttor = 10.0
      cutcv = 0.2
      cutbf = 5.0
      cutmp = 1.0
      ncolen = 2
      gappen = 5.0
      gaploc = 3.5
      nwcut = 3.5
      nwmode = 'squaredist'
      dmax = 3.5
      maxdp = 10
      mamcdf = 0.25
c
      castrt = 3.0
      castop = 3.0
      castep = 1.0
c
      lvrml = .false.
      ivrml = 99
      cavrml = ' CA '
      vrdist = 4.5
      vrmlbg = 'black'
      vrmldc = 'white'
      vrfile = 'lsqman.wrl'
      rgbbg1 = 0.0
      rgbbg2 = 0.0
      rgbbg2 = 0.0
      rgbfg1 = 1.0
      rgbfg2 = 1.0
      rgbfg2 = 1.0
      call xvrml_init ()
c
      chncol (1) = 'yellow'
      chncol (2) = 'green'
      chncol (3) = 'cyan'
      chncol (4) = 'magenta'
      chncol (5) = 'blue'
      chncol (6) = 'red'
      chncol (7) = 'purple'
      chncol (8) = 'orange'
      chncol (9) = 'honeydew'
      chncol (10) = 'peru'
      chncol (11) = 'light_salmon'
      chncol (12) = 'tomato3'
      chncol (13) = 'darkorchid2'
      chncol (14) = 'dark_khaki'
      chncol (15) = 'thistle'
      chncol (16) = 'navy_blue'
      chncol (17) = 'violetred'
      chncol (18) = 'olivedrab1'
      chncol (19) = 'crimson'
      chncol (20) = 'ivory1'
      chncol (21) = 'beige'
      chncol (22) = 'forestgreen'
      chncol (23) = 'peachpuff3'
      chncol (24) = 'sienna2'
      chncol (25) = 'mint_cream'
      chncol (26) = 'bisque'
c
c ... user input unit (5=interactive; other=macro)
c
      munit = 5
c
      linit = .false.
c
      linter = xinter()
      if (linter) then
        pro='$LSQMAN > '
      else
        pro=' LSQMAN > '
      end if
c
      call setdef (0)
c
      do i=1,maxmol
        name   (i) = '$#@%^&@^62'
        file   (i) = ' '
        coment (i) = 'No comment'
        incore (i) = .false.
        mmodel (i) = .false.
        natoms (i) = 0
        nchain (i) = 0
        do j=1,3
          cell (j,i) = 1.0
          cell (j+3,i) = 90.0
        end do
        do j=1,maxmol
          call initrt (i,j)
        end do
      end do
c
      natype = 1
      atypes (1) = ' CA '
      do i=2,maxtyp
        atypes (i) = ' '
      end do
c
      chamod = 'RE'
      lhkeep = .false.
      ltkeep = .true.
      lnmral = .true.
c
      omcent = ' CA '
      omdist = 4.5
      omconn = 'all.dat'
c
c ... formats
c
c ... TO DO: COnvolution mol1 chain1 mol2 chain2 frag_length 2d_plot_file
c ... TO DO: LOcal_nw mol1 chain1 mol2 chain2 frag_length gap_penalty
c ... TO DO: GLobal register-error-detection parameters
c
 6000 format (/
     +  ' LSQMAN options :'//
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
     + /' REad mol pdb_file [chain] [atom]    ',
     +              ' WRite mol pdb_file [chain] [first] [last]'/
     +  ' DElete mol                          ',
     +              ' ANnotate mol comment_string'/
     +  ' LIst [mol]                          ',
     +              ' CHain_mode mode'/
     +  ' TYpe_residues mol                   ',
     +              ' BFactor_range b_lo b_hi'/
     +  ' FRactionalise mol                   ',
     +              ' ORthogonalise mol'/
     +  ' CEll mol a b c al be ga             ',
     +              ' SUbtract_ave_b mol'/
     +  ' HYdrogens keep_or_strip             ',
     +              ' HEtatm keep_or_strip'/
     +  ' NMr_model_mode all_or_first         ',
     +              ' AA_substitution_matrix filename'/
     + /' ALter CHain_id mol chain new_chain  ',
     +              ' ALter SEgid mol segid new_segid'/
     +  ' ALter FOrce mol chain new_segid     ',
     +              ' ALter SAme mol chain'/
     +  ' ALter REnumber mol chain [first]    ',
     +              ' '/
     + /' EXplicit mol1 range1 mol2 range2    ',
     +              ' NWunsch mol1 chain1 mol2 chain2 gap'/
     +  ' BRute_force mol1 chain1 mol2 chain2 frag_length ',
     +              'frag_step min_match [S|D] [slide_step]'/
     +  ' FAst_force  mol1 chain1 mol2 chain2 frag_length ',
     +              'frag_step min_match [S|D] [slide_step]'/
     +  ' XAlignment mol1 chain1 mol2 chain2 p',
     +              'ir_alignment_file'/
     + /' IMprove mol1 range1 mol2 range2     ',
     +              ' '/
     +  ' DP_improve mol1 chain1 mol2 chain2 m',
     +              'ode cut_off max_cycles [verbose]'/
     + /' GLobal_nw mol1 chain1 mol2 chain2 cu',
     +              't_off [log_file]'/
     + /' EDit_operator mol1 mol2 val1 ...    ',
     +              ' SHow_operator mol1 mol2'/
     +  ' SAve_operator mol1 mol2 file [name] ',
     +              ' PErturb_operator mol1 mol2 [amplitude]'/
     +  ' APply_operator mol1 mol2_to_move [ch',
     +              'ain] [first] [last]'/
     +  ' OLd_o_operator mol1 mol2 file       ',
     +              ' RMsd_calc mol1 range1 mol2 range2 [Ltarget]'/
     +  ' WAters mol1 mol2 cut_off plot_file  ',
     +              ' HIsto_dist mol1 mol2 cut_off bin'/
     + /' MOrph mol1 range1 mol2 range2 nsteps',
     +              ' basename type oid range3 cutoff'/
     +  ' SImilarity_plot mol1 chain1 mol2 cha',
     +              'in2 plot_file [start] [end] [step]'/
     +  ' LEsk_plot mol1 chain1 mol2 chain2 pl',
     +              'ot_file'/
     +  ' PHipsi_plot mol1 range1 mol2 range2 ',
     +              'plot_file [cut-off] [hist_bin] [hist_max]'/
     +  ' ETa_theta_plot mol1 range1 mol2 rang',
     +              'e2 plot_file [cut-off] [hist_bin] [hist_max]'/
     +  ' DIstance_plot mol1 range1 mol2 range',
     +              '2 plot_file [cut-off] [hist_bin] [hist_max]'/
     +  ' DDihe_plot mol1 range1 mol2 range2 pl',
     +              'ot_file [cut-off] [hist_bin] [hist_max]'/
     +  ' D1_D2_plot mol1 range1 mol2 range2 pl',
     +              'ot_file [cut-off] [hist_bin] [hist_max]'/
     +  ' QDiff_dist_plot mol1 range1 mol2 rang',
     +              'e2 2d_plot_file'/
     +  ' DChi mol1 range1 mol2 range2 [cut-off',
     +              '] '/
     +  ' SOap_film mol1 chain1 mol2 chain2 odl',
     +              '_file [verbose]'/
     + /' MCentral mol residue_range exp_imp [',
     +              'rt_file]  '/
     +  ' MAlign mol residue_range exp_imp cha',
     +              'in'/
     +  ' MDihedral mol chain plot_file [cut] ',
     +              ' MRamachandran mol chain ps_file [cut] [how]'/
     +  ' MSide_ch mol chain plot_file [cut]  ',
     +              ' MTorsion mol chain ps_file [cut] [how]'/
     +  ' VMain_ch mol chain plot_file [cut]  ',
     +              ' VSide_ch mol chain plot_file [cut]'/
     +  ' MPlot mol chain plot_file ps_file [d',
     +              'max_black] [cut_dist_print]'/
     +  ' MBfactors mol chain plot_file [cut] ',
     +              ' '/
     + /' JUdge target tchn parent pchn model ',
     +              'mchn dist phi chi '/
     +  ' CAsp target tchn model mchn [start] ',
     +              '[end] [step]'/
     + /' GEt XYz mol chain x y z radius symbo',
     +              'l_name [O_macro]'/
     + /' FIx_atom_names mol1 range1 mol2 rang',
     +              'e2 mode how what [min_gain] [cut_off]'/
     +  ' NOmenclature mol                    ',
     +              ' INvert_ncs infile outfile'/
     +  ' NUcleic_acid_pdb_nomenclature mol   ',
     +              ' '/
     + /' ATom_types ?                        ',
     +              ' ATom_types CA'/
     +  ' ATom_types MAin_chain               ',
     +              ' ATom_types SIde_chain'/
     +  ' ATom_types EXtended_main_chain      ',
     +              ' ATom_types ALl'/
     +  ' ATom_types NOn_hydrogen             ',
     +              ' ATom_types DEfine type1 [type2 ...]'/
     +  ' ATom_types PHosphorous              ',
     +              ' ATom_types TRace_and_side_chain'/
     +  ' ATom_types C4*                      ',
     +              ' ATom_types NUcleic_acid_backbone'/
     + /' SEt ?                               ',
     +              ' SEt REset_defaults'/
     +  ' SEt COarse_6A_fit_defaults          ',
     +              ' SEt INtermediate_4A_fit_defaults'/
     +  ' SEt FIne_tune_3A_fit_defaults       ',
     +              ' SEt SImilar_mols_2A_fit_defaults'/
     + /' SEt MAx_nr_improve_cycles value     ',
     +              ' SEt DIst_max value'/
     +  ' SEt MIn_fragment_length value       ',
     +              ' SEt DEcay value'/
     +  ' SEt OPtimisation_criterion value    ',
     +              ' SEt SEquential_hits on_off'/
     +  ' SEt RMs_weight value                ',
     +              ' SEt FRagment_length_decay value'/
     +  ' SEt SHift_correction on_off         ',
     +              ' SEt NUcleic_acid_defaults'/
     + /' OMacro INit mol1 file               ',
     +              ' OMacro APpend mol2'/
     +  ' OMacro WRite o_command_string       ',
     +              ' OMacro CLose_file'/
     +  ' OMacro DEfine central_atom max_dist ',
     +              'connect_file'/
     + /' VRml SEtup central_atom max_dist backgr_col default_col'/
     +  ' VRml INit [vrml_file]               ',
     +              ' VRml COlour_list'/
     +  ' VRml ADd mol [chain] [colour]       ',
     +              ' VRml ALl_chains mol'/
     +  )
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6000)
c
      call jvalut (' Max nr of molecules             :',1,maxmol)
      call jvalut (' Max nr of residues per molecule :',1,maxres)
      call jvalut (' Max nr of atoms per molecule    :',1,maxatm)
      call jvalut (' Max nr of atom types            :',1,maxtyp)
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
        call gknval ('GKLSQMAN',inimac,ierr)
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
c ... READ
c
      else if (optpar(1)(1:2) .eq. 'RE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'm9999'
          call textin (' New mol ?',optpar(2))
        end if
        call allocm (optpar(2),imol,ierr)
        if (ierr .ne. 0) goto 10
        prev = name(imol)
c
        if (nopt .lt. 3) then
          optpar (3) = file(imol)
          call textin (' File name ?',optpar(3))
        end if
        file (imol) = optpar(3)
c
        if (nopt .lt. 4) optpar(4) = '*'
        call upcase (optpar(4))
c
        if (nopt .lt. 5) optpar(5) = '*'
        call upcase (optpar(5))
c
        call molin (imol,optpar(4)(1:1),optpar(5)(1:4),ierr)
        close (iunit)
c
        if (ierr .ne. 0) goto 10
c
        incore (imol) = .true.
        coment (imol) = 'Read from '//file(imol)
c
        goto 10
c
c ... BFACTOR_RANGE
c
      else if (optpar(1)(1:2) .eq. 'BF') then
c
        if (nopt .lt. 2) then
          write (optpar (2),*) bcutlo
          call textin (' Lower B-factor cut-off ?',optpar(2))
        end if
        call str2r (optpar(2),xdum,ierr)
        if (ierr .eq. 0) bcutlo = xdum
c
        if (nopt .lt. 3) then
          write (optpar (3),*) bcuthi
          call textin (' Upper B-factor cut-off ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .eq. 0) bcuthi = xdum
c
        call rlohi (bcutlo,bcuthi)
        call fvalut (' Lower B cut-off :',1,bcutlo)
        call fvalut (' Upper B cut-off :',1,bcuthi)
c
c ... NOMENCLATURE
c
      else if (optpar(1)(1:2) .eq. 'NO') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        call textut (' Enforce proper nomenclature for :',name(imol))
        call nomen (imol)
c
c ... NUCLEIC ACID PDB NOMENCLATURE
c
      else if (optpar(1)(1:2) .eq. 'NU') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        call textut (' PDB NA nomenclature for :',name(imol))
        call nomenu (imol)
c
c ... SUBTRACT_AVE_B
c
      else if (optpar(1)(1:2) .eq. 'SU') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        call subavb (imol)
c
c ... FRACTIONALISE
c
      else if (optpar(1)(1:2) .eq. 'FR') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        call orthog (cell(1,imol),ca2fra,1)
        call fvalut (' Operator :',12,ca2fra)
c
        call fvalut (' Atom #1 before :',3,atmxyz(1,1,imol))
        do i=1,natoms(imol)
          call mulmtx (ca2fra,atmxyz(1,i,imol),dummy,3,3,1)
          atmxyz (1,i,imol) = dummy (1)
          atmxyz (2,i,imol) = dummy (2)
          atmxyz (3,i,imol) = dummy (3)
        end do
        call fvalut (' Atom #1 after  :',3,atmxyz(1,1,imol))
c
        call textut (' Fractionalised :',name(imol))
c
c ... ORTHOGONALISE
c
      else if (optpar(1)(1:2) .eq. 'OR') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        call orthog (cell(1,imol),fra2ca,0)
        call fvalut (' Operator :',12,fra2ca)
c
        call fvalut (' Atom #1 before :',3,atmxyz(1,1,imol))
        do i=1,natoms(imol)
          call mulmtx (fra2ca,atmxyz(1,i,imol),dummy,3,3,1)
          atmxyz (1,i,imol) = dummy (1)
          atmxyz (2,i,imol) = dummy (2)
          atmxyz (3,i,imol) = dummy (3)
        end do
        call fvalut (' Atom #1 after  :',3,atmxyz(1,1,imol))
c
        call textut (' Orthogonalised :',name(imol))
c
c ... CELL
c
      else if (optpar(1)(1:2) .eq. 'CE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          write (optpar(3),*) cell(1,imol)
          call textin (' A axis (A) ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .eq. 0 .and. xdum .ge. 1.0) cell(1,imol) = xdum
c
        if (nopt .lt. 4) then
          write (optpar(4),*) cell(2,imol)
          call textin (' B axis (A) ?',optpar(4))
        end if
        call str2r (optpar(4),xdum,ierr)
        if (ierr .eq. 0 .and. xdum .ge. 1.0) cell(2,imol) = xdum
c
        if (nopt .lt. 5) then
          write (optpar(5),*) cell(3,imol)
          call textin (' C axis (A) ?',optpar(5))
        end if
        call str2r (optpar(5),xdum,ierr)
        if (ierr .eq. 0 .and. xdum .ge. 1.0) cell(3,imol) = xdum
c
        if (nopt .lt. 6) then
          write (optpar(6),*) cell(4,imol)
          call textin (' Alpha (deg) ?',optpar(6))
        end if
        call str2r (optpar(6),xdum,ierr)
        if (ierr .eq. 0 .and. xdum .ge. 1.0 .and. xdum .le. 179.0)
     +    cell(4,imol) = xdum
c
        if (nopt .lt. 7) then
          write (optpar(7),*) cell(5,imol)
          call textin (' Beta (deg) ?',optpar(7))
        end if
        call str2r (optpar(7),xdum,ierr)
        if (ierr .eq. 0 .and. xdum .ge. 1.0 .and. xdum .le. 179.0)
     +    cell(5,imol) = xdum
c
        if (nopt .lt. 8) then
          write (optpar(8),*) cell(6,imol)
          call textin (' Gamma (deg) ?',optpar(8))
        end if
        call str2r (optpar(8),xdum,ierr)
        if (ierr .eq. 0 .and. xdum .ge. 1.0 .and. xdum .le. 179.0)
     +    cell(6,imol) = xdum
c
        call textut (' Molecule :',name(imol))
        call fvalut (' Cell axes (A) :',3,cell(1,imol))
        call fvalut (' Angles (deg)  :',3,cell(4,imol))
c
c ... WRITE
c
      else if (optpar(1)(1:2) .eq. 'WR') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = file(imol)
          call textin (' File name ?',optpar(3))
        end if
        file (imol) = optpar(3)
c
        if (nopt .lt. 4) optpar (4) = '*'
        call upcase (optpar(4))
c
        i1 = -9999
        i2 = -9999
        if (nopt .ge. 5) then
          call str2i (optpar(5),idum,ierr)
          if (ierr .eq. 0) i1 = idum
          if (nopt .ge. 6) then
            call str2i (optpar(6),idum,ierr)
            if (ierr .eq. 0) i2 = idum
          end if
        end if
c
        call textut (' Write mol :',name(imol))
        call textut (' Chain id  :',optpar(4))
        call textut (' PDB file  :',file(imol))
        if (i1 .ne. -9999) call ivalut (' First res :',1,i1)
        if (i2 .ne. -9999) call ivalut (' Last  res :',1,i2)
c
        call molut (imol,ierr,optpar(4),i1,i2)
c
        if (ierr .ne. 0) goto 10
c
        goto 10
c
c ... DELETE
c
      else if (optpar(1)(1:2) .eq. 'DE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        call selecm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do imol=1,maxmol
          if (select(imol)) then
c
            call textut (' Deleted :',name(imol))
c
            incore (imol) = .false.
            name (imol) = '&^%$%'
            file (imol) = ' '
            mmodel (imol) = .false.
            natoms (imol) = 0
            nchain (imol) = 0
            do j=1,maxmol
              call initrt (imol,j)
              call initrt (j,imol)
            end do
c
 1210       continue
          end if
        end do
c
c ... PERTURB
c
      else if (optpar(1)(1:2) .eq. 'PE') then
c
        call prompt (' Perturb operator ...')
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Mol 2 ?',optpar(3))
        end if
        jmol = whichm (optpar(3))
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f8.1)') ranamp
          call remspa (optpar(4))
          call textin (' Amplitude of random perturbations ?',
     +      optpar(4))
        end if
        call str2r (optpar(4),xdum,ierr)
        if (ierr .eq. 0) ranamp = xdum
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (imol .eq. jmol) then
          call prompt (' Warning: mol1 == mol2 !')
cc          call errcon ('Unit operator: mol1 == mol2 !')
cc          goto 10
        end if
c
        call textut (' Molecule 1 :',name(imol))
        call textut (' Molecule 2 :',name(jmol))
        call fvalut (' Amplitude  :',1,ranamp)
c
        call perrot (rtlsq(1,imol,jmol),ranamp,.false.)
c
        call inista (imol,jmol)
c
        write (last(imol,jmol),'(99(a,1x))',err=10)
     +    (optpar(i)(1:leng1(optpar(i))),i=1,4)
c
c ... APPLY
c
      else if (optpar(1)(1:2) .eq. 'AP') then
c
        call prompt (' Bring Mol 2 on top of Mol 1 ...')
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 (operator to use) ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Mol 2 (to be moved) ?',optpar(3))
        end if
        jmol = whichm (optpar(3))
c
        if (nopt .lt. 4) optpar (4) = '*'
        call upcase (optpar(4))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (imol .eq. jmol) then
          call prompt (' Warning: mol1 == mol2 !')
cc          call errcon ('Unit operator: mol1 == mol2 !')
cc          goto 10
        end if
c
        i1 = -9999
        i2 = -9999
        if (nopt .ge. 5) then
          call str2i (optpar(5),idum,ierr)
          if (ierr .eq. 0) i1 = idum
          if (nopt .ge. 6) then
            call str2i (optpar(6),idum,ierr)
            if (ierr .eq. 0) i2 = idum
          end if
        end if
c
c ... apply operator
c
        call textut (' Molecule 1 :',name(imol))
        call textut (' Molecule 2 :',name(jmol))
        call textut (' Apply to mol 2 chain :',optpar(4))
        if (i1 .ne. -9999) call ivalut (' First res :',1,i1)
        if (i2 .ne. -9999) call ivalut (' Last  res :',1,i2)
c
        call prompt (' Applying operator to mol 2 ...')
        call vecrtv (atmxyz(1,1,jmol),buffi,natoms(jmol),
     +               rtlsq(1,imol,jmol),rtlsq(10,imol,jmol))
        call prompt (' Updating selected chain(s)/zone ...')
c
        k = 0
        do i=1,natoms(jmol)
c
          if (optpar(4) .ne. '*') then
            if (optpar(4) .ne. achain(i,jmol)) goto 1938
          end if
c
          if (i1 .gt. -9999) then
            if (iresid(i,jmol) .lt. i1) goto 1938
          end if
c
          if (i2 .gt. -9999) then
            if (iresid(i,jmol) .gt. i2) goto 1938
          end if
c
          atmxyz (1,i,jmol) = buffi (1,i)
          atmxyz (2,i,jmol) = buffi (2,i)
          atmxyz (3,i,jmol) = buffi (3,i)
          k = k + 1
c
 1938     continue
        end do
c
        call jvalut (' Nr of atoms moved :',1,k)
c
c ... reset all operators to other molecules
c
        if (optpar(4) .eq. '*' .and. i1 .eq. -9999 .and.
     +      i2 .eq. -9999) then
          call prompt (' Resetting ALL operators of mol 2 ...')
          do j=1,maxmol
            call initrt (j,jmol)
            call initrt (jmol,j)
          end do
        else
          call prompt (' NOT resetting operators of mol 2 !')
        end if
c
        write (last(imol,jmol),'(99(a,1x))',err=10)
     +    (optpar(i)(1:leng1(optpar(i))),i=1,4)
c
ccc        write (last(jmol,imol),'(99(a,1x))',err=10)
ccc     +    (optpar(i)(1:leng1(optpar(i))),i=1,4)
c
c ... CHAIN_MODE
c
      else if (optpar(1)(1:2) .eq. 'CH') then
c
        if (nopt .lt. 2) then
          optpar (2) = chamod
          write (*,'(1x,a)')
     +      'Select one of the following modes:',
     +      'REname    = chains are renamed A, B, .. Z',
     +      'ORiginal  = chain names are not altered',
     +      'NOn-blank = keep chain names but replace " " by "_"',
     +      'XPlor     = rename; SEGIds delineate chains',
     +      'BReak     = rename; use breaks in residue numbers',
     +      'LOwer     = rename; use drop in residue numbers'
          call textin (' Chain mode ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'RE') then
          chamod = 'RE'
          call prompt (' Chain-mode REname')
        else if (optpar(2)(1:2) .eq. 'OR') then
          chamod = 'OR'
          call prompt (' Chain-mode ORiginal')
        else if (optpar(2)(1:2) .eq. 'NO') then
          chamod = 'NO'
          call prompt (' Chain-mode NOn-blank')
        else if (optpar(2)(1:2) .eq. 'XP') then
          chamod = 'XP'
          call prompt (' Chain-mode XPlor')
        else if (optpar(2)(1:2) .eq. 'BR') then
          chamod = 'BR'
          call prompt (' Chain-mode BReak')
        else if (optpar(2)(1:2) .eq. 'LO') then
          chamod = 'LO'
          call prompt (' Chain-mode LOwer')
        else
          call errcon ('Invalid chain-mode (not RE/OR/NO/XP/BR/LO)')
        end if
c
c ... HYDROGENS
c
      else if (optpar(1)(1:2) .eq. 'HY') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'STrip'
          write (*,'(1x,a)')
     +      'Select one of the following modes:',
     +      'KEep   = retain hydrogens on read/write',
     +      'STrip  = strip  hydrogens on read/write'
          call textin (' Hydrogen mode ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'KE') then
          lhkeep = .true.
          call prompt (' Keep hydrogens')
        else if (optpar(2)(1:2) .eq. 'ST') then
          lhkeep = .false.
          call prompt (' Strip hydrogens')
        else
          call errcon ('Invalid hydrogen mode (not KE/ST)')
        end if
c
c ... NMR
c
      else if (optpar(1)(1:2) .eq. 'NM') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'ALl'
          write (*,'(1x,a)')
     +      'Select one of the following modes:',
     +      'ALl   = keep all NMR models on read',
     +      'FIrst = only keep first NMR model on read'
          call textin (' NMR model mode ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'AL') then
          lnmral = .true.
          call prompt (' Keep all NMR models')
        else if (optpar(2)(1:2) .eq. 'FI') then
          lnmral = .false.
          call prompt (' Only keep first NMR model')
        else
          call errcon ('Invalid NMR model mode (not AL/FI)')
        end if
c
c ... AA_substitution_matrix
c
      else if (optpar(1)(1:2) .eq. 'AA') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'sbin_blosum45.lib'
          call gklibf (optpar(2))
          call textin (' Library file with matrix ?',optpar(2))
        end if
        call remspa (optpar(2))
        call textut (' Library file with matrix :',optpar(2))
        call xopxoa (iunit,optpar(2),linter,ierr)
c
        if (ierr .ne. 0) then
          call errcon ('While opening library file; using default')
          call inicmp (cmpmat,maxaat)
          goto 1841
        end if
c
        call rdlib (iunit,maxaat,typ1lc,cmpmat,ierr)
c
        if (ierr .ne. 0) then
          call errstp ('While reading library file; using default')
          call inicmp (cmpmat,maxaat)
          goto 1841
        end if
c
        call prompt (' Matrix read successfully !')
c
 1841   continue
        close (iunit)
c
c ... HETATMS
c
      else if (optpar(1)(1:2) .eq. 'HE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'STrip'
          write (*,'(1x,a)')
     +      'Select one of the following modes:',
     +      'KEep   = retain HETATMs on read',
     +      'STrip  = strip  HETATMs on read'
          call textin (' HETATM mode ?',optpar(2))
        end if
        call upcase (optpar(2))
        call remspa (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'KE') then
          ltkeep = .true.
          call prompt (' Keep HETATMs')
        else if (optpar(2)(1:2) .eq. 'ST') then
          ltkeep = .false.
          call prompt (' Strip HETATMs')
        else
          call errcon ('Invalid HETATM mode (not KE/ST)')
        end if
c
c ... LIST
c
      else if (optpar(1)(1:2) .eq. 'LI') then
c
        if (nopt .lt. 2) optpar (2) = '*'
        call selecm (optpar(2),ierr)
        if (ierr .ne. 0) goto 10
        prev = optpar (2)
c
        do imol=1,maxmol
          if (select(imol)) then
c
            write (*,*)
            call textut (' List    :',name(imol))
c
            call textut (' File    :',file(imol))
            call textut (' Comment :',coment(imol))
            call fvalut (' Cell    :',6,cell(1,imol))
            call jvalut (' Nr of atoms in mol  :',1,natoms(imol))
            call logiut (' Multiple NMR models ?',1,mmodel(imol))
            call jvalut (' Nr of chains/models :',1,nchain(imol))
c
            do i=1,nchain(imol)
              write (*,6009) i,chname(i,imol),
     +          (chnptr(2,i,imol)-chnptr(1,i,imol)+1)
            end do
c
          end if
        end do
c
 6009 format (' Chain/Model # ',i2,' - Name |',a1,'| Nr of atoms ',
     +  i8)
c
c ... TYPE_RESIDUES
c
      else if (optpar(1)(1:2) .eq. 'TY') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        call textut (' List of residues in :',name(imol))
        call seqnce (imol)
c
c ... ANNOTATE
c
      else if (optpar(1)(1:2) .eq. 'AN') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Which mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = coment(imol)
          call textin (' Label ?',optpar(3))
        end if
c
        coment (imol) = optpar (3)
c
c ... EXPLICIT or IMPROVE or RMSD_CALC or PHIPSI or DIST or DDIHE or D1D2
c     or FIX_ATOM_NAMES or MORPH or DChi or ETa_theta
c
      else if (optpar(1)(1:2) .eq. 'EX' .or.
     +         optpar(1)(1:2) .eq. 'IM' .or.
     +         optpar(1)(1:2) .eq. 'RM' .or.
     +         optpar(1)(1:2) .eq. 'DI' .or.
     +         optpar(1)(1:2) .eq. 'DD' .or.
     +         optpar(1)(1:2) .eq. 'DC' .or.
     +         optpar(1)(1:2) .eq. 'QD' .or.
     +         optpar(1)(1:2) .eq. 'D1' .or.
     +         optpar(1)(1:2) .eq. 'PH' .or.
     +         optpar(1)(1:2) .eq. 'FI' .or.
     +         optpar(1)(1:2) .eq. 'MO' .or.
     +         optpar(1)(1:2) .eq. 'ET') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'A1-999'
          call textin (' Range 1 ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) optpar(3) = 'A1-999'
        call upcase (optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Mol 2 ?',optpar(4))
        end if
        jmol = whichm (optpar(4))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'A1'
          call textin (' Range 2 ?',optpar(5))
        end if
        if (length(optpar(5)) .lt. 1) optpar(5) = 'A1'
        call upcase (optpar(5))
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - mol1 == mol2 !')
        end if
c
        if (optpar(1)(1:2) .eq. 'EX') then
          call lsqexp (imol,optpar(3),jmol,optpar(5),1,.true.,ierr)
c
        else if (optpar(1)(1:2) .eq. 'RM') then
c
          if (nopt .ge. 6) then
            call str2i (optpar(6),ltarg,ierr)
            if (ierr .ne. 0 .or. ltarg .le. 0) ltarg = 1
            ltarg = -ltarg
          else
            ltarg = -1
          end if
c
          call lsqexp (imol,optpar(3),jmol,optpar(5),ltarg,.true.,ierr)
c
        else if (optpar(1)(1:2) .eq. 'IM') then
          call lsqimp (imol,optpar(3),jmol,optpar(5),.true.,ierr)
c
        else if (optpar(1)(1:2) .eq. 'MO') then
c
          if (nopt .lt. 6) then
            write (optpar(6),*) nmorph
            call textin (' Number of steps ?',optpar(6))
          end if
          call str2i (optpar(6),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .lt. 3 .or. idum .gt. 999) then
            call errcon ('Invalid number of steps')
            goto 10
          end if
          nmorph = idum
c
          if (nopt .lt. 7) then
            optpar (7) = 'morphy'
            call textin (' Basename for output PDB files ?',
     +                   optpar(7))
          end if
          call remspa (optpar(7))
c
          if (nopt .lt. 8) then
            optpar (8) = 'internal'
            call textin (' Morph type (Internal/Cartesian) ?',optpar(8))
          end if
          call remspa (optpar(8))
          call upcase (optpar(8))
          if (optpar(8)(1:1) .ne. 'C') optpar (8)='I'
c
          if (nopt .lt. 9) then
            optpar (9) = 'm'
            call textin (' O mol name suffix ?',optpar(9))
          end if
          call remspa (optpar(9))
          call upcase (optpar(9))
          if (ichar(optpar(9)(1:1)) .lt. ichar('A') .or.
     +        ichar(optpar(9)(1:1)) .gt. ichar('Z')) then
            optpar (9) = 'M'
          end if
c
          if (nopt .lt. 10) then
            optpar (10) = optpar(3)
            call textin (' Superpositioning range ?',optpar(10))
          end if
          if (length(optpar(10)) .lt. 1) optpar(10) = optpar(3)
          call upcase (optpar(10))
          if (index(optpar(10),'*') .gt. 0) then
            call errcon ('Range may not contain wildcard (*)')
            optpar(10) = optpar(3)
          end if
c
          if (nopt .lt. 11) then
            write (optpar(11),*) tarcut
            call textin (' Torsion range cut-off ?',optpar(11))
          end if
          call str2r (optpar(11),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .gt. 0.0) tarcut = xdum
c
          call morphm (imol,optpar(3),jmol,optpar(5),nmorph,
     +                 optpar(7),optpar(8),optpar(9),
     +                 optpar(10),tarcut,ierr)
c
        else if (optpar(1)(1:2) .eq. 'FI') then
c
          if (nopt .lt. 6) then
            optpar (6) = 'S'
            call textin (' Mode (Strict/All) ?',optpar(6))
          end if
          call remspa (optpar(6))
          call upcase (optpar(6))
          if (optpar(6)(1:1) .ne. 'A') optpar(6) = 'Strict'
c
          if (nopt .lt. 7) then
            optpar (7) = 'S'
            call textin (' How (Sequential/Nearest) ?',optpar(7))
          end if
          call remspa (optpar(7))
          call upcase (optpar(7))
          if (optpar(7)(1:1) .ne. 'N') optpar(7) = 'Sequential'
c
          if (nopt .lt. 8) then
            optpar (8) = 'T'
            call textin (' Minimise what (Torsions/Rmsd) ?',optpar(8))
          end if
          call remspa (optpar(8))
          call upcase (optpar(8))
          if (optpar(8)(1:1) .ne. 'T') optpar(8) = 'Rmsd'
c
          if (nopt .ge. 9) then
            call str2r (optpar(9),xdum,ierr)
            if (ierr .eq. 0) delta = xdum
          end if
c
          if (nopt .ge. 10) then
            call str2r (optpar(10),xdum,ierr)
            if (ierr .eq. 0) cacut = xdum
          end if
c
          call fixatm (imol,optpar(3),jmol,optpar(5),
     +                 cacut,optpar(6)(1:1),optpar(7)(1:1),
     +                 optpar(8)(1:1),delta,ierr)
c
        else if (optpar(1)(1:2) .eq. 'QD') then
c
          if (nopt .lt. 6) then
            optpar (6) = name(imol)//'_'//name(jmol)//
     +                   '_diff_dist.pl2'
            call remspa (optpar(6))
            call locase (optpar(6))
            call textin (' Plot file ?',optpar(6))
          end if
          call remspa (optpar(6))
          call plotqd (imol,optpar(3),jmol,optpar(5),optpar(6),ierr,
     +                 maxbuf,pl1buf)
c
        else if (optpar(1)(1:2) .eq. 'DC') then
c
          qqq = 10.0
          if (nopt .ge. 6) then
            call str2r (optpar(6),xdum,ierr)
            if (ierr .eq. 0) qqq = xdum
          end if
c
          call delchi (imol,optpar(3),jmol,optpar(5),qqq,ierr)
c
        else
c
c ... various types of plots
c
          if (nopt .lt. 6) then
            optpar (6) = name(imol)//'_'//name(jmol)//'_'//
     +        optpar(1)(1:2)//'.plt'
            call remspa (optpar(6))
            call locase (optpar(6))
            call textin (' Plot file ?',optpar(6))
          end if
          call remspa (optpar(6))
c
          xdum = -1.0
          if (nopt .ge. 7) then
            call str2r (optpar(7),xdum,ierr)
            if (ierr .ne. 0) xdum = -1.0
          end if
c
          qqq = -1.0
          if (nopt .ge. 8) then
            call str2r (optpar(8),qqq,ierr)
            if (ierr .ne. 0) qqq = -1.0
          end if
c
          q3 = -1.0
          if (nopt .ge. 8) then
            call str2r (optpar(9),q3,ierr)
            if (ierr .ne. 0) q3 = -1.0
          end if
c
          call plotem (optpar(1)(1:2),imol,optpar(3),
     +                 jmol,optpar(5),optpar(6),q3,qqq,xdum,ierr)
        end if
c
        if (ierr .ne. 0) goto 10
c
        if (optpar(1)(1:2) .eq. 'EX' .or.
     +      optpar(1)(1:2) .eq. 'IM') then
          write (last(imol,jmol),'(99(a,1x))',err=10)
     +      (optpar(i)(1:leng1(optpar(i))),i=1,5)
        end if
c
c ... MULTIPLE DIHEDRALS/B-FACTORS/RMSD
c
      else if (optpar(1)(1:2) .eq. 'MD' .or.
     +         optpar(1)(1:2) .eq. 'MB' .or.
     +         optpar(1)(1:2) .eq. 'MR' .or.
     +         optpar(1)(1:2) .eq. 'VM' .or.
     +         optpar(1)(1:2) .eq. 'VS' .or.
     +         optpar(1)(1:2) .eq. 'MS' .or.
     +         optpar(1)(1:2) .eq. 'MT' .or.
     +         optpar(1)(1:2) .eq. 'MP') then
c
        if (optpar(1)(1:2) .eq. 'MD') then
          call prompt (' Multiple chain/model dihedral plot')
        else if (optpar(1)(1:2) .eq. 'MB') then
          call prompt (' Multiple chain/model B-factor plot')
        else if (optpar(1)(1:2) .eq. 'MR') then
          call prompt (' Multiple Ramachandran plot')
        else if (optpar(1)(1:2) .eq. 'MS') then
          call prompt (' Multiple side-chain torsion plot')
        else if (optpar(1)(1:2) .eq. 'MT') then
          call prompt (' Multiple side-chain Chi-1/Chi-2 plot')
        else if (optpar(1)(1:2) .eq. 'VM') then
          call prompt (' Main-chain circular variance plot')
        else if (optpar(1)(1:2) .eq. 'VS') then
          call prompt (' Side-chain circular variance plot')
        else if (optpar(1)(1:2) .eq. 'MP') then
          call prompt (' Multiple chain/model RMS distance plot')
        end if
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nchain(imol) .lt. 2) then
          call errcon ('Fewer than 2 chains/models in this mol')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = 'A'
          call textin (' Reference chain ?',optpar(3))
        end if
        call upcase (optpar(3))
        call remspa (optpar(3))
        if (optpar(3) .eq. ' ') optpar(3) = 'A'
c
c ... DIHEDRALS
c
        if (optpar(1)(1:2) .eq. 'MD') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_phi_psi_sigma.plt'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' Plot file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'phi_psi_sigma.plt'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cuttor = xdum
          end if
c
          call muldih (imol,optpar(3),optpar(4),1,'C',cuttor,ierr,
     +      prognm)
c
c ... RAMACHANDRAN
c
        else if (optpar(1)(1:2) .eq. 'MR') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_multi_rama.ps'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' PostScript file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'multi_rama.ps'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cuttor = xdum
          end if
c
          if (nopt .lt. 6) then
            optpar(6) = 'Cartesian'
          end if
          call remspa (optpar(6))
          call upcase (optpar(6))
          if (optpar(6)(1:1).ne.'P') optpar(6)(1:1)='C'
c
          call muldih (imol,optpar(3),optpar(4),2,optpar(6),cuttor,ierr,
     +      prognm)
c
c ... SIDE CHAINS
c
        else if (optpar(1)(1:2) .eq. 'MS') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_chi12_sigma.plt'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' Plot file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'chi12_sigma.plt'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cuttor = xdum
          end if
c
          call mulsid (imol,optpar(3),optpar(4),1,'C',cuttor,ierr)
c
c ... CIRC VAR MAIN
c
        else if (optpar(1)(1:2) .eq. 'VM') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_phi_psi_cv.plt'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' Plot file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'phi_psi_cv.plt'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cutcv = xdum
          end if
c
          call muldih (imol,optpar(3),optpar(4),3,'C',cutcv,ierr,
     +      prognm)
c
c ... CIRC VAR SIDE
c
        else if (optpar(1)(1:2) .eq. 'VS') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_chi12_cv.plt'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' Plot file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'chi12_cv.plt'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cutcv = xdum
          end if
c
          call mulsid (imol,optpar(3),optpar(4),3,'C',cutcv,ierr)
c
c ... TORSIONS
c
        else if (optpar(1)(1:2) .eq. 'MT') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_chi12_dist.ps'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' PostScript file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'chi12_dist.ps'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cuttor = xdum
          end if
c
          if (nopt .lt. 6) then
            optpar(6) = 'Cartesian'
          end if
          call remspa (optpar(6))
          call upcase (optpar(6))
          if (optpar(6)(1:1).ne.'P') optpar(6)(1:1)='C'
c
          call mulsid (imol,optpar(3),optpar(4),2,optpar(6),cuttor,ierr)
c
c ... B-FACTORS
c
        else if (optpar(1)(1:2) .eq. 'MB') then
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_bfac_multi.plt'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' Plot file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'bfac_multi.plt'
c
          if (nopt .ge. 5) then
            call str2r (optpar(5),xdum,ierr)
            if (ierr .eq. 0) cutbf = xdum
          end if
c
          call mulbfs (imol,optpar(3),optpar(4),cutbf,ierr)
c
c ... MULTIPLE RMSD PLOT
c
        else if (optpar(1)(1:2) .eq. 'MP') then
c
          if (nopt .lt. 4) then
            optpar (4) = name(imol)(1:leng1(name(imol)))//
     +                   '_multi_rms_dist.plt'
            call remspa (optpar(4))
            call locase (optpar(4))
            call textin (' Plot file ?',optpar(4))
          end if
          call remspa (optpar(4))
          if (optpar(4) .eq. ' ') optpar(4) = 'multi_rms_dist.plt'
c
          if (nopt .lt. 5) then
            optpar (5) = name(imol)(1:leng1(name(imol)))//
     +                   '_cdplot.ps'
            call remspa (optpar(5))
            call locase (optpar(5))
            call textin (' CD plot PostScript file ?',optpar(5))
          end if
          call remspa (optpar(5))
          if (optpar(5) .eq. ' ') optpar(5) = 'cdplot.ps'
c
          if (nopt .ge. 6) then
            call str2r (optpar(6),xdum,ierr)
            if (ierr .eq. 0) dmax = xdum
          end if
c
          if (nopt .ge. 7) then
            call str2r (optpar(7),xdum,ierr)
            if (ierr .eq. 0) cutmp = xdum
          end if
c
          call mulrms (imol,optpar(3),optpar(4),optpar(5),
     +                 dmax,cutmp,ierr,prognm)
c
        end if
c
c ... MULTIPLE ALIGNMENT
c
      else if (optpar(1)(1:2) .eq. 'MC' .or.
     +         optpar(1)(1:2) .eq. 'MA') then
c
        if (optpar(1)(1:2) .eq. 'MC') then
          call prompt (' Find central chain/model')
        else if (optpar(1)(1:2) .eq. 'MA') then
          call prompt (' Align all to one chain/model')
        end if
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nchain(imol) .lt. 3) then
          call errcon ('Fewer than 3 chains/models in this mol')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = '1-999'
          call textin (' Residue range (NO CHAIN !) ?',optpar(3))
        end if
        if (length(optpar(3)) .lt. 1) optpar(3) = '1-999'
        call upcase (optpar(3))
        call pretty (optpar(3))
        if (index('1234567890',optpar(3)(1:1)) .le. 0) then
          call errcon ('Chain name removed')
          optpar (3) = optpar(3)(2:)
        end if
        idum = index (optpar(3),' ')
        optpar(3)(idum+1:) = ' '
        idum = length(optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = 'Improve'
          call textin (' Explicit or Improve (E/I) ?',optpar(4))
        end if
        call upcase (optpar(4))
        call remspa (optpar(4))
        if (optpar(4)(1:1) .ne. 'E') optpar(4) = 'I'
c
        if (optpar(1)(1:2) .eq. 'MA') then
c
c ... MA
c
          if (nopt .lt. 5) then
            optpar (5) = 'A'
            call textin (' Central chain ?',optpar(5))
          end if
          call upcase (optpar(5))
          call remspa (optpar(5))
          if (optpar(5) .eq. ' ') optpar(5) = 'A'
          goto 1741
c
        else
c
c ... MC
c
          if (nopt .ge. 5) then
            lrtf = .true.
            close (iunit)
            call xopxua (iunit,optpar(5),linter,ierr)
            if (ierr .ne. 0) then
              call errcon ('Cannot open operator file')
              lrtf = .false.
            else
              call textut (' Save operators in file :',optpar(5))
              call stamp (sdum)
              write (iunit,'(a,a)') '! ',sdum(1:leng1(sdum))
            end if
          else
            lrtf = .false.
          end if
        end if
c
c ... MC === do explicit alignment of similar chains/models
c
        do i=1,nchain(imol)-1
          xmat (i,i) = 0.0
          nat1 = chnptr(2,i,imol)-chnptr(1,i,imol)+1
          zone1 = chname(i,imol)//optpar(3)(1:idum)
          zone3 = chname(i,imol)//'*'
c
          do j=i+1,nchain(imol)
c
            xmat (j,i) = 999.99
c
            if (nat1 .le. 3) then
ccc              call errcon ('Too few atoms in chain 1')
              goto 1638
            end if
c
            nat2 = chnptr(2,j,imol)-chnptr(1,j,imol)+1
            if (nat2 .le. 3) then
ccc              call errcon ('Too few atoms in chain 2')
              goto 1638
            end if
c
            xdum = float(nat1-nat2)/float(min(nat1,nat2))
c
ccc            print *,' I,J,NAT1, NAT2, DIF% = ',i,j,nat1,nat2,xdum
c
            if (abs(xdum) .gt. mamcdf) then
ccc              call errcon ('Nr of atoms in chains differs too much')
              goto 1638
            end if
c
            write (*,'(/1x,a,a1,a,a1)') 'Aligning ',chname(j,imol),
     +        ' to ',chname(i,imol)
            zone2 = chname(j,imol)//optpar(3)(1:idum)
            rmsd (imol,imol) = 999.99
            call lsqexp (imol,zone1,imol,zone2,1,.false.,ierr)
            xmat (j,i) = rmsd (imol,imol)
c
c ... improve alignment if requested
c
            if (optpar(4)(1:1) .eq. 'I') then
              zone4 = chname(j,imol)//'*'
              rmsd (imol,imol) = 999.99
              call lsqimp (imol,zone3,imol,zone4,.false.,ierr)
              xmat (j,i) = rmsd (imol,imol)
            end if
c
            if (xmat(j,i) .gt. 900.) then
ccc              call errcon ('Calculated RMSD > 900 A')
              goto 1638
            end if
c
            write (*,6010) nmatch(imol,imol),rmsd(imol,imol),
     +        simind(imol,imol),matchi(imol,imol),
     +        cripp(imol,imol),rrmsd(imol,imol),
     +        normsd(imol,imol),
     +        sas1(imol,imol),sas2(imol,imol),
     +        sas3(imol,imol),sas4(imol,imol),
     +        rmsdna(imol,imol),
     +        rmsb(imol,imol),corb(imol,imol),
     +        rtlsq(1,imol,imol),rtlsq(4,imol,imol),
     +        rtlsq(7,imol,imol),rtlsq(2,imol,imol),
     +        rtlsq(5,imol,imol),rtlsq(8,imol,imol),
     +        rtlsq(3,imol,imol),rtlsq(6,imol,imol),
     +        rtlsq(9,imol,imol),rtlsq(10,imol,imol),
     +        rtlsq(11,imol,imol),rtlsq(12,imol,imol)
c
            if (lrtf) then
              sdum = '.lsq_rt_' //
     +               name(imol)(1:leng1(name(imol))) //
     +               chname(j,imol)// '_to_' //
     +               name(imol)(1:leng1(name(imol))) //
     +               chname(i,imol)
              call locase (sdum)
              call textut (' Datablock name :',sdum)
              write (iunit,'(a,a)',err=3917)
     +          sdum(1:leng1(sdum)),' r 12 (3f15.7)'
              write (iunit,'(3f15.7)',err=3917)
     +          (rtlsq(k,imol,imol),k=1,12)
              goto 3919
c
 3917         continue
              call errcon ('While writing datablock file')
c
 3919         continue
            end if
c
 1638       continue
c
            xmat (i,j) = xmat (j,i)
          end do
        end do
c
        close (iunit)
c
 6047   format (' Chain/model ',a1,' - RMSDs (A) to the others:'/
     +    4(7f10.3,:,/))
c
c ... for each chain, calc RMS(RMSD)
c
        ibest = -1
        xbest = 9999.
        iii = 0
        qqq = 0.0
        do i=1,nchain(imol)
          idum = 0
          rmsrms (i) = 0.0
          write (*,*)
          write (*,6047) chname(i,imol),(xmat(i,j),j=1,nchain(imol))
          do j=1,nchain(imol)
            if (i .ne. j .and. xmat(i,j) .lt. 999.) then
              idum = idum + 1
              rmsrms (i) = rmsrms (i) + (xmat(i,j)*xmat(i,j))
              if (j .gt. i) then
                iii = iii + 1
                qqq = qqq + xmat(i,j)
              end if
            end if
          end do
          if (idum .le. 0) then
            write (*,'(1x,a,a1)') 'NO matches for chain/model ',
     +        chname(i,imol)
          else
            rmsrms (i) = sqrt (rmsrms(i) / float(idum))
            write (*,'(1x,a,a1,a,f8.3)') 'RMS(RMSD) for chain/model ',
     +        chname(i,imol),' = ',rmsrms(i)
            if (rmsrms(i) .lt. xbest) then
              xbest = rmsrms(i)
              ibest = i
            end if
          end if
        end do
c
c ... central chain is the one with the lowest RMS(RMSD)
c
        write (*,*)
        if (ibest .gt. 0) then
          write (*,'(1x,a,a1)') '==> Central chain is ',
     +      chname(ibest,imol)
        else
          call prompt (' ==> NO central chain/model !!??')
        end if
        if (iii .gt. 0) then
          write (*,*)
          write (*,'(1x,a,f8.3,a)')
     +      'Average RMSD between chains = ',qqq/float(iii),' A'
        end if
        write (*,*)
        write (last(imol,imol),'(9(a,1x))',err=10)
     +    (optpar(i)(1:leng1(optpar(i))),i=1,4)
c
        goto 10
c
c ... MA === alignment to a chain
c
 1741   continue
c
c ... do explicit alignment of chains/models to central one
c
        ichn = -1
        do i=1,nchain(imol)
          if (chname(i,imol) .eq. optpar(5)(1:1)) then
            ichn = i
            goto 1743
          end if
        end do
        call errcon ('Requested chain not found')
        goto 10
c
 1743   continue
        i = ichn
        nat1 = chnptr(2,ichn,imol)-chnptr(1,ichn,imol)+1
        zone1 = chname(ichn,imol)//optpar(3)(1:idum)
        zone3 = chname(ichn,imol)//'*'
c
        do j=1,nchain(imol)
c
ccc          print *,' chain ',j
ccc          print *,' ptr i ',chnptr(1,ichn,imol),chnptr(2,ichn,imol)
ccc          print *,' ptr j ',chnptr(1,j,imol),chnptr(2,j,imol)
c
          if (j .eq. ichn) goto 1738
c
          if (nat1 .le. 3)  then
ccc            call errcon ('Too few atoms in chain 1')
            goto 1738
          end if
c
          nat2 = chnptr(2,j,imol)-chnptr(1,j,imol)+1
          if (nat2 .le. 3)  then
ccc            call errcon ('Too few atoms in chain 2')
            goto 1738
          end if
c
          xdum = float(nat1-nat2)/float(min(nat1,nat2))
c
ccc          print *,' I,J,NAT1, NAT2, DIF% = ',i,j,nat1,nat2,xdum
c
          if (abs(xdum) .gt. mamcdf) then
ccc            call errcon ('Nr of atoms in chains differs too much')
            goto 1738
          end if
c
          zone2 = chname(j,imol)//optpar(3)(1:idum)
          call lsqexp (imol,zone1,imol,zone2,1,.false.,ierr)
c
c ... improve alignment if requested
c
          if (optpar(4)(1:1) .eq. 'I') then
            zone4 = chname(j,imol)//'*'
            call lsqimp (imol,zone3,imol,zone4,.false.,ierr)
          end if
c
c ... align the chain
c
          write (*,'(1x,a,a1,a,f8.3)') 'Chain/model ',chname(j,imol),
     +      ' - RMSD (A) to the selected one = ',rmsd(imol,imol)
          call prompt (' Applying operator to chain ...')
c
c ... it's not efficient to do this for all atoms, but it's simple
c     (no index-juggling in the next DO-loop ;-)
c
          call vecrtv (atmxyz(1,1,imol),buffi,natoms(imol),
     +                 rtlsq(1,imol,imol),rtlsq(10,imol,imol))
c
          do k=chnptr(1,j,imol),chnptr(2,j,imol)
            atmxyz (1,k,imol) = buffi (1,k)
            atmxyz (2,k,imol) = buffi (2,k)
            atmxyz (3,k,imol) = buffi (3,k)
          end do
c
 1738     continue
c
        end do
c
c ... reset operator of IMOL to itself
c
        call initrt (imol,imol)
c
        write (last(imol,imol),'(9(a,1x))',err=10)
     +    (optpar(i)(1:leng1(optpar(i))),i=1,5)
c
        goto 10
c
c ... JUDGE
c
      else if (optpar(1)(1:2) .eq. 'JU') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 = TARGET ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'A'
          call textin (' Chain 1 ?',optpar(3))
        end if
        call upcase (optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Mol 2 = PARENT ?',optpar(4))
        end if
        jmol = whichm (optpar(4))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'A'
          call textin (' Chain 2 ?',optpar(5))
        end if
        call upcase (optpar(5))
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - TARGET == PARENT !')
        end if
c
        if (nopt .lt. 6) then
          optpar (6) = prev
          call textin (' Mol 3 = MODEL ?',optpar(6))
        end if
        kmol = whichm (optpar(6))
c
        if (kmol .le. 0 .or. kmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(kmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 7) then
          optpar (7) = 'A'
          call textin (' Chain 3 ?',optpar(7))
        end if
        call upcase (optpar(7))
c
        if (imol .eq. kmol) then
          call prompt (' WARNING - TARGET == MODEL !')
c          call errcon ('TARGET == MODEL !')
c          goto 10
        end if
c
        if (jmol .eq. kmol) then
          call prompt (' WARNING - MODEL == PARENT !')
        end if
c
        if (nopt .lt. 8) then
          write (optpar (8),*) judis
          call textin (' CA-CA distance cut-off ?',optpar(8))
        end if
        call str2r (optpar(8),xdum,ierr)
        if (ierr .ne. 0) goto 10
        judis = xdum
c
        if (nopt .lt. 9) then
          write (optpar (9),*) juphi
          call textin (' Phi-Psi distance cut-off ?',optpar(9))
        end if
        call str2r (optpar(9),xdum,ierr)
        if (ierr .ne. 0) goto 10
        juphi = xdum
c
        if (nopt .lt. 10) then
          write (optpar (10),*) juchi
          call textin (' Chi1-2 distance cut-off ?',optpar(10))
        end if
        call str2r (optpar(10),xdum,ierr)
        if (ierr .ne. 0) goto 10
        juchi = xdum
c
c ... do it
c
        ichn = -1
        do i=1,nchain(imol)
          if (chname(i,imol) .eq. optpar(3)(1:1)) then
            ichn = i
            goto 4343
          end if
        end do
        call errcon ('TARGET chain not found')
        goto 10
 4343   continue
c
        jchn = -1
        do i=1,nchain(jmol)
          if (chname(i,jmol) .eq. optpar(5)(1:1)) then
            jchn = i
            goto 4345
          end if
        end do
        call errcon ('PARENT chain not found')
        goto 10
 4345   continue
c
        kchn = -1
        do i=1,nchain(kmol)
          if (chname(i,kmol) .eq. optpar(7)(1:1)) then
            kchn = i
            goto 4347
          end if
        end do
        call errcon ('MODEL chain not found')
        goto 10
 4347   continue
c
        call casp3 (imol,ichn,jmol,jchn,kmol,kchn,judis,juphi,juchi)
c
        goto 10
c
c ... CASP
c
      else if (optpar(1)(1:2) .eq. 'CA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 = TARGET ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'A'
          call textin (' Chain 1 ?',optpar(3))
        end if
        call upcase (optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Mol 2 = MODEL ?',optpar(4))
        end if
        kmol = whichm (optpar(4))
c
        if (kmol .le. 0 .or. kmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(kmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'A'
          call textin (' Chain 2 ?',optpar(5))
        end if
        call upcase (optpar(5))
c
        if (imol .eq. kmol) then
          call prompt (' WARNING - TARGET == MODEL !')
c          call errcon ('TARGET == MODEL !')
c          goto 10
        end if
c
        if (nopt .ge. 6) then
          call str2r (optpar(6),xdum,ierr)
          if (ierr .ne. 0) goto 10
          castrt = xdum
        end if
c
        if (nopt .ge. 7) then
          call str2r (optpar(7),xdum,ierr)
          if (ierr .ne. 0) goto 10
          castop = xdum
        else
          castop = castrt
        end if
c
        if (nopt .ge. 8) then
          call str2r (optpar(8),xdum,ierr)
          if (ierr .ne. 0) goto 10
          castep = xdum
        end if
c
        if (nopt .lt. 9) optpar (9) = 'E'
        call remspa (optpar(9))
        call upcase (optpar(9))
        if (optpar(9)(1:1) .ne. 'C') optpar (9) = 'E'
c
c ... do it
c
        ichn = -1
        do i=1,nchain(imol)
          if (chname(i,imol) .eq. optpar(3)(1:1)) then
            ichn = i
            goto 4348
          end if
        end do
        call errcon ('TARGET chain not found')
        goto 10
 4348   continue
c
        kchn = -1
        do i=1,nchain(kmol)
          if (chname(i,kmol) .eq. optpar(5)(1:1)) then
            kchn = i
            goto 4349
          end if
        end do
        call errcon ('MODEL chain not found')
        goto 10
 4349   continue
c
        call casp2 (imol,ichn,kmol,kchn,castrt,castop,castep,
     +              optpar(9))
c
        goto 10
c
c ... GET
c
      else if (optpar(1)(1:2) .eq. 'GE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'XYz'
          call textin (' Get option ?',optpar(2))
          nopt = 2
        end if
c
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'XY') then
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Mol ?',optpar(3))
          end if
          imol = whichm (optpar(3))
c
          if (imol .le. 0 .or. imol .gt. maxmol) then
            call errcon ('Invalid mol selection')
            goto 10
          end if
c
          if (.not. incore(imol)) then
            call errcon ('Mol not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) then
            optpar (4) = 'A'
            call textin (' Chain ?',optpar(4))
          end if
          call upcase (optpar(4))
c
          if (nopt .lt. 5) then
            optpar (5) = '0.0'
            call textin (' X ?',optpar(5))
          end if
          call str2r (optpar(5),dummy(1),ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            optpar (6) = '0.0'
            call textin (' Y ?',optpar(6))
          end if
          call str2r (optpar(6),dummy(2),ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 7) then
            optpar (7) = '0.0'
            call textin (' Z ?',optpar(7))
          end if
          call str2r (optpar(7),dummy(3),ierr)
          if (ierr .ne. 0) goto 10
          call fvalut (' X, Y, Z :',3,dummy)
c
          if (nopt .lt. 8) then
            optpar (8) = '5.0'
            call textin (' Radius (A) ?',optpar(8))
          end if
          call str2r (optpar(8),radius,ierr)
          call fvalut (' Radius (A) :',1,radius)
          if (ierr .ne. 0) goto 10
          if (radius .le. 0.001) goto 10
c
          if (nopt .lt. 9) then
            optpar (9) = 'xyzsel'
            call textin (' Symbol name ?',optpar(9))
          end if
          call remspa (optpar(9))
          if (length(optpar(9)) .lt. 1) goto 10
c
          punit = -1
          if (nopt .lt. 10) then
            optpar (10) = ' '
          else
            punit = 50
            call xopxua (punit,optpar(10),linter,ierr)
            if (ierr .ne. 0) then
              punit = -1
            else
              call textut (' Creating O macro :',optpar(10))
              write (punit,'(20(a,1x))')
     +          '! LSQMAN-generated O macro',
     +          optpar(10)(1:length(optpar(10)))
              write (punit,'(20(a,1x))')
     +          '! Command :',
     +          (optpar(i)(1:length(optpar(i))),i=1,10)
              write (punit,'(20(a,1x))') '!'
              write (punit,'(20(a,1x))') '! molec #Molecule ?#'
              sdum = optpar(9)(1:4)//'_'//optpar(4)(1:1)
              call remspa (sdum)
              call upcase (sdum)
              call textut (' Object name in O :',sdum)
              write (punit,'(20(a,1x))') 'object',
     +          sdum(1:length(sdum))
            end if
          end if
          if (punit .le. 0) call prompt (' No O macro generated')
c
c ... do it
c
          ichn = -1
          do i=1,nchain(imol)
            if (chname(i,imol) .eq. optpar(4)(1:1)) then
              ichn = i
              goto 3453
            end if
          end do
          call errcon ('Requested chain not found')
          goto 10
 3453     continue
c
 4713 format (' # ',i5,' @ ',f6.2,' A -> ',1x,a4,1x,a3,1x,a1,
     +        i4,1x,3x,3f8.3,2f6.2,6x,a4)
c
          optpar (1) = '&'
          optpar (2) = optpar (9)
          optpar (5) = optpar (4)
          optpar (3) = ' '
          jres = 0
          nout = 0
          do i=chnptr(1,ichn,imol),chnptr(2,ichn,imol)
            if (iresid(i,imol) .ne. jres) then
              xdum = distce (dummy,atmxyz(1,i,imol))
              if (xdum .le. radius) then
                nout = nout + 1
                jres = iresid(i,imol)
                write (optpar(4),'(1x,a1,i6)') optpar(5)(1:1),jres
                call remspa (optpar(4)(2:))
                kk = i
          write (*,4713) nout,xdum,
     +      atmnam(kk,imol),resnam(kk,imol),achain(kk,imol),
     +      iresid(kk,imol),atmxyz(1,kk,imol),atmxyz(2,kk,imol),
     +      atmxyz(3,kk,imol),qatom(kk,imol),batom(kk,imol),
     +      axplor(kk,imol)
ccc                call textut (' Add :',optpar(4))
                if (length(optpar(3)) .lt. 1) then
                  optpar (3) = optpar (4)
                else
                  call appstr (optpar(3),optpar(4))
                end if
c
                if (punit .gt. 0) then
                  write (punit,'(20(a,1x))') 'zone',
     +              optpar(4)(1:length(optpar(4)))
                end if
              end if
            end if
          end do
c
          call textut (' Selection :',optpar(3))
          call dosymb (3,optpar,ldone)
c
          if (punit .gt. 0) then
            write (punit,'(20(a,1x))') 'end'
            call prompt (' O macro written')
            close (punit)
            punit = -1
          end if
c
          goto 10
c
        else
          call errcon ('Invalid GEt option')
          goto 10
        end if
c
c ... CONVOLUTION/LOCAL_NW/GLOBAL_NW/NWUNSCH/DP_IMPROVE/SOAP_FILM/XALIGN
c
      else if (optpar(1)(1:2) .eq. 'CO' .or.
     +         optpar(1)(1:2) .eq. 'LO' .or.
     +         optpar(1)(1:2) .eq. 'GL' .or.
     +         optpar(1)(1:2) .eq. 'DP' .or.
     +         optpar(1)(1:2) .eq. 'SO' .or.
     +         optpar(1)(1:2) .eq. 'NW' .or.
     +         optpar(1)(1:2) .eq. 'XA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'A'
          call textin (' Chain 1 ?',optpar(3))
        end if
        call upcase (optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Mol 2 ?',optpar(4))
        end if
        jmol = whichm (optpar(4))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'A'
          call textin (' Chain 2 ?',optpar(5))
        end if
        call upcase (optpar(5))
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - mol1 == mol2 !')
        end if
c
c ... OPTPAR 6
c
        if (optpar(1)(1:2) .eq. 'CO' .or.
     +      optpar(1)(1:2) .eq. 'LO') then
          if (nopt .lt. 6) then
            write (optpar (6),*) ncolen
            call textin (' Fragment length (one-sided) ?',optpar(6))
          end if
          call str2i (optpar(6),idum,ierr)
          if (ierr .ne. 0) goto 10
          ncolen = idum
        else if (optpar(1)(1:2) .eq. 'NW') then
          if (nopt .lt. 6) then
            write (optpar (6),*) gappen
            call textin (' Gap penalty ?',optpar(6))
          end if
          call str2r (optpar(6),xdum,ierr)
          if (ierr .ne. 0) goto 10
          gappen = xdum
        else if (optpar(1)(1:2) .eq. 'GL') then
          if (nopt .lt. 6) then
            write (optpar (6),*) nwcut
            call textin (' Cut-off distance (A) ?',optpar(6))
          end if
          call str2r (optpar(6),xdum,ierr)
          if (ierr .ne. 0) goto 10
          nwcut = xdum
        else if (optpar(1)(1:2) .eq. 'DP') then
          if (nopt .lt. 6) then
            optpar (6) = nwmode
            call textin (' Matrix mode (S/D) ?',optpar(6))
          end if
          call remspa (optpar(6))
          call upcase (optpar(6))
          nwmode = optpar(6)
        else if (optpar(1)(1:2) .eq. 'SO') then
          if (nopt .lt. 6) then
            optpar (6) = 'soap_film.odl'
            call textin (' ODL file name ?',optpar(6))
          end if
          call remspa (optpar(6))
        else if (optpar(1)(1:2) .eq. 'XA') then
          if (nopt .lt. 6) then
            optpar (6) = 'alignment.pir'
            call textin (' Alignment file name ?',optpar(6))
          end if
          call remspa (optpar(6))
        end if
c
c ... OPTPAR 7
c
        if (optpar(1)(1:2) .eq. 'CO') then
          if (nopt .lt. 7) then
            optpar (7) = 'convolution.pl2'
            call textin (' Plot file ?',optpar(7))
          end if
        else if (optpar(1)(1:2) .eq. 'LO') then
          if (nopt .lt. 7) then
            write (optpar (7),*) gaploc
            call textin (' Gap penalty ?',optpar(7))
          end if
          call str2r (optpar(7),xdum,ierr)
          if (ierr .ne. 0) goto 10
          gaploc = xdum
        else if (optpar(1)(1:2) .eq. 'GL') then
          if (nopt .lt. 7) optpar (7) = ' '
        else if (optpar(1)(1:2) .eq. 'SO') then
          if (nopt .lt. 7) then
            optpar (7) = 'No'
          end if
          call upcase (optpar(7))
          call remspa (optpar(7))
        else if (optpar(1)(1:2) .eq. 'NW') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'XA') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'DP') then
          if (nopt .lt. 7) then
            write (optpar (7),*) nwcut
            call textin (' Cut-off distance (A) ?',optpar(7))
          end if
          call str2r (optpar(7),xdum,ierr)
          if (ierr .ne. 0) goto 10
          nwcut = xdum
        end if
c
c ... OPTPAR 8
c
        if (optpar(1)(1:2) .eq. 'CO') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'LO') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'GL') then
          if (nopt .lt. 8) optpar (8) = '-1'
          call str2i (optpar(8),idum,ierr)
          if (ierr .ne. 0) goto 10
          iremin = idum
        else if (optpar(1)(1:2) .eq. 'SO') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'NW') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'XA') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'DP') then
          if (nopt .lt. 8) then
            write (optpar (8),*) maxdp
            call textin (' Max nr of cycles ?',optpar(8))
          end if
          call str2i (optpar(8),idum,ierr)
          if (ierr .ne. 0) goto 10
          maxdp = idum
        end if
c
c ... OPTPAR 9
c
        if (optpar(1)(1:2) .eq. 'CO') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'LO') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'GL') then
          if (nopt .lt. 9) optpar (9) = '10'
          call str2i (optpar(9),idum,ierr)
          if (ierr .ne. 0) goto 10
          ireshi = abs(idum)
        else if (optpar(1)(1:2) .eq. 'SO') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'NW') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'XA') then
c ... nothing
        else if (optpar(1)(1:2) .eq. 'DP') then
          if (nopt .lt. 9) then
            optpar (9) = 'No'
          end if
          call upcase (optpar(9))
          call remspa (optpar(9))
        end if
c
c ... do it
c
        ichn = -1
        do i=1,nchain(imol)
          if (chname(i,imol) .eq. optpar(3)(1:1)) then
            ichn = i
            goto 3243
          end if
        end do
        call errcon ('Requested chain in mol 1 not found')
        goto 10
 3243   continue
c
        jchn = -1
        do i=1,nchain(jmol)
          if (chname(i,jmol) .eq. optpar(5)(1:1)) then
            jchn = i
            goto 3245
          end if
        end do
        call errcon ('Requested chain in mol 2 not found')
        goto 10
 3245   continue
c
        if (optpar(1)(1:2) .eq. 'CO') then
          write (*,6022) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),ncolen,optpar(7)(1:length(optpar(7)))
        else if (optpar(1)(1:2) .eq. 'LO') then
          write (*,6024) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),ncolen,gaploc
        else if (optpar(1)(1:2) .eq. 'GL') then
          write (*,6226) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),nwcut,optpar(7)(1:leng1(optpar(7)))
ccc     +      iremin,ireshi
        else if (optpar(1)(1:2) .eq. 'DP') then
          write (*,6230) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),nwcut,nwmode,maxdp,optpar(9)
        else if (optpar(1)(1:2) .eq. 'NW') then
          write (*,6228) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),gappen
        else if (optpar(1)(1:2) .eq. 'XA') then
          write (*,6232) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),optpar(6)(1:leng1(optpar(6)))
        else if (optpar(1)(1:2) .eq. 'SO') then
          write (*,6224) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      atypes(1),optpar(6)(1:leng1(optpar(6))),
     +      optpar(7)
        end if
c
 6022 format (' Convolution of  ',a,1x,a/
     +        ' And             ',a,1x,a/
     +        ' Atom type      |',a4,'|'/
     +        ' Fragment length ',i4/
     +        ' Plot file       ',a)
c
 6024 format (
     + ' Local-RMSD-based Needleman-Wunsch alignment'/
     +        ' Of              ',a,1x,a/
     +        ' And             ',a,1x,a/
     +        ' Atom type      |',a4,'|'/
     +        ' Fragment length ',i4/
     +        ' Gap penalty     ',f8.2)
c
 6224 format (
     + ' Soap-film ODL file generation'/
     +        ' Of               ',a,1x,a/
     +        ' And              ',a,1x,a/
     +        ' Atom type       |',a4,'|'/
     +        ' ODL file name    ',a/
     +        ' Verbose output   ',a10)
c
 6226 format (
     + ' Global-superposition-distance-based',
     +   ' Needleman-Wunsch alignment'/
     +        ' Of                   ',a,1x,a/
     +        ' And                  ',a,1x,a/
     +        ' Atom type           |',a4,'|'/
     +        ' Cut-off distance     ',f8.2/
     +        ' Log file (optional)  ',a/)
ccc     +        ' Register-shift par 1 ',i8/
ccc     +        ' Register-shift par 2 ',i8)
c
 6228 format (
     + ' Sequence-based Needleman-Wunsch alignment'/
     +        ' Of              ',a,1x,a/
     +        ' And             ',a,1x,a/
     +        ' Atom type      |',a4,'|'/
     +        ' Gap penalty     ',f8.2)
c
 6230 format (
     + ' Dynamic-Programming-based',
     +   ' operator improvement (Needleman-Wunsch)'/
     +        ' Of               ',a,1x,a/
     +        ' And              ',a,1x,a/
     +        ' Atom type       |',a4,'|'/
     +        ' Cut-off distance ',f8.2/
     +        ' Matrix mode      ',a10/
     +        ' Max nr of cycles ',i8/
     +        ' Verbose output   ',a10)
c
 6232 format (
     + ' External alignment'/
     +        ' Of               ',a,1x,a/
     +        ' And              ',a,1x,a/
     +        ' Atom type       |',a4,'|'/
     +        ' Alignment file   ',a)
c
        k1 = 0
        do i=chnptr(1,ichn,imol),chnptr(2,ichn,imol)
          if (atmnam(i,imol) .eq. atypes(1)) then
            k1 = k1 + 1
            buffi (1,k1) = atmxyz (1,i,imol)
            buffi (2,k1) = atmxyz (2,i,imol)
            buffi (3,k1) = atmxyz (3,i,imol)
            if (optpar(1)(1:2) .eq. 'LO' .or.
     +          optpar(1)(1:2) .eq. 'GL' .or.
     +          optpar(1)(1:2) .eq. 'DP' .or.
     +          optpar(1)(1:2) .eq. 'NW' .or.
     +          optpar(1)(1:2) .eq. 'XA') then
              l = -1
              do k=1,maxaat
                if (resnam(i,imol) .eq. typ3lc(k)) then
                  l = k
                  goto 8310
                end if
              end do
 8310         continue
              seq1(k1) = '?'
              iseq1(k1) = l
              if (l .gt. 0) then
                seq1(k1) = typ1lc(l)
              else
                call nuctyp (resnam(i,imol),nt)
                nt = max (0, min (nt, maxnuc))
                seq1 (k1) = nuc1lc (nt)
              end if
              ibufi(k1) = iresid(i,imol)
            end if
          end if
        end do
c
        k2 = 0
        do j=chnptr(1,jchn,jmol),chnptr(2,jchn,jmol)
          if (atmnam(j,jmol) .eq. atypes(1)) then
            k2 = k2 + 1
            buffj (1,k2) = atmxyz (1,j,jmol)
            buffj (2,k2) = atmxyz (2,j,jmol)
            buffj (3,k2) = atmxyz (3,j,jmol)
            if (optpar(1)(1:2) .eq. 'LO' .or.
     +          optpar(1)(1:2) .eq. 'GL' .or.
     +          optpar(1)(1:2) .eq. 'DP' .or.
     +          optpar(1)(1:2) .eq. 'NW' .or.
     +          optpar(1)(1:2) .eq. 'XA') then
              l = -1
              do k=1,maxaat
                if (resnam(j,jmol) .eq. typ3lc(k)) then
                  l = k
                  goto 8320
                end if
              end do
 8320         continue
              seq2(k2) = '?'
              iseq2(k2) = l
              if (l .gt. 0) then
                seq2(k2) = typ1lc(l)
              else
                call nuctyp (resnam(j,jmol),nt)
                nt = max (0, min (nt, maxnuc))
                seq2 (k2) = nuc1lc (nt)
              end if
              ibufj(k2) = iresid(j,jmol)
            end if
          end if
        end do
c
        call jvalut (' Central atoms mol 1 :',1,k1)
        call jvalut (' Central atoms mol 2 :',1,k2)
c
        if (k1 .lt. 3 .or. k2. lt. 3) then
          call errcon ('Not enough atoms')
          goto 10
        end if
c
        if (optpar(1)(1:2) .eq. 'CO' .or.
     +      optpar(1)(1:2) .eq. 'LO') then
          idum = 3 + (2*ncolen+1)
          if (k1 .lt. idum .or. k2 .lt. idum) then
            call errcon (' Not enough residues')
            call jvalut (' Minimum required :',1,idum)
            goto 10
          end if
        end if
c
c ... call the appropriate subroutine
c
        if (optpar(1)(1:2) .eq. 'CO') then
          call convol (optpar(7),imol,ichn,jmol,jchn,k1,k2,
     +                 ncolen,buffi,buffj,pl1buf,maxbuf)
        else if (optpar(1)(1:2) .eq. 'LO') then
          idum = maxbuf
          jdum = (k1+1) * (k2+1)
          if (jdum .gt. idum) then
            call errcon ('Not enough memory - sorry !')
            call jvalut (' Required  :',1,jdum)
            call jvalut (' Available :',1,idum)
          else
            call nwlocs (imol,ichn,jmol,jchn,k1,k2,
     +                   ncolen,buffi,buffj,pl1buf,
     +                   gaploc,pl2buf,pl3buf,seq1,seq2,
     +                   rtlsq(1,imol,jmol))
            write (last(imol,jmol),'(99(a,1x))',err=10)
     +        (optpar(i)(1:leng1(optpar(i))),i=1,7)
          end if
        else if (optpar(1)(1:2) .eq. 'GL') then
          idum = maxbuf
          jdum = (k1+1) * (k2+1)
          if (jdum .gt. idum) then
            call errcon ('Not enough memory - sorry !')
            call jvalut (' Required  :',1,jdum)
            call jvalut (' Available :',1,idum)
          else
            call nwglob (imol,ichn,jmol,jchn,k1,k2,
     +                   nwcut,buffi,buffj,pl1buf,
     +                   pl2buf,pl3buf,seq1,seq2,
     +                   rtlsq(1,imol,jmol),optpar(7),
     +                   iremin,ireshi,chname(ichn,imol),
     +                   chname(jchn,jmol),ibufi,ibufj)
          end if
        else if (optpar(1)(1:2) .eq. 'DP') then
          idum = maxbuf
          jdum = (k1+1) * (k2+1)
          if (jdum .gt. idum) then
            call errcon ('Not enough memory - sorry !')
            call jvalut (' Required  :',1,jdum)
            call jvalut (' Available :',1,idum)
          else
            nold = nmatch (imol,jmol)
            xold = rmsd (imol,jmol)
            iter = 1
 2257       continue
            write (*,*)
            call jvalut (' DP_improve iteration :',1,iter)
c
ccc            call fvalut (' RT before :',12,rtlsq(1,imol,jmol))
c
            call nwimpr (imol,ichn,jmol,jchn,k1,k2,
     +                   nwcut,buffi,buffj,buffk,buffl,
     +                   pl1buf,pl2buf,pl3buf,seq1,seq2,
     +                   rtlsq(1,imol,jmol),nwmode,optpar(9))
c
ccc            call fvalut (' RT after  :',12,rtlsq(1,imol,jmol))
ccc            print *,nold,nmatch(imol,jmol),xold,rmsd(imol,jmol)
c
            if (iter .lt. maxdp .and.
     +          (nold .ne. nmatch(imol,jmol) .or.
     +           abs(xold-rmsd(imol,jmol)) .ge. 0.0005)) then
              iter = iter + 1
              nold = nmatch (imol,jmol)
              xold = rmsd (imol,jmol)
              goto 2257
            end if
            write (last(imol,jmol),'(99(a,1x))',err=10)
     +        (optpar(i)(1:leng1(optpar(i))),i=1,9)
          end if
        else if (optpar(1)(1:2) .eq. 'SO') then
          idum = maxbuf
          jdum = (k1+1) * (k2+1)
          if (jdum .gt. idum) then
            call errcon ('Not enough memory - sorry !')
            call jvalut (' Required  :',1,jdum)
            call jvalut (' Available :',1,idum)
          else
            call nwsoap (k1,k2,
     +                   optpar(6),buffi,buffj,pl1buf,pl2buf,
     +                   pl3buf,rtlsq(1,imol,jmol),optpar(7))
          end if
        else if (optpar(1)(1:2) .eq. 'NW') then
          call nwseq (imol,ichn,jmol,jchn,k1,k2,
     +                buffi,buffj,cmpmat,maxaat,
     +                gappen,pl2buf,pl3buf,seq1,seq2,
     +                iseq1,iseq2,rtlsq(1,imol,jmol))
          write (last(imol,jmol),'(99(a,1x))',err=10)
     +      (optpar(i)(1:leng1(optpar(i))),i=1,6)
        else if (optpar(1)(1:2) .eq. 'XA') then
ccc          write (*,'(a,1x,10a1)') 'SEQ1',(seq1(i),i=1,10)
ccc          write (*,'(a,1x,10a1)') 'SEQ2',(seq2(i),i=1,10)
          call extali (imol,ichn,jmol,jchn,k1,k2,
     +                 buffi,buffj,optpar(6),seq1,seq2,
     +                 rtlsq(1,imol,jmol))
          write (last(imol,jmol),'(99(a,1x))',err=10)
     +      (optpar(i)(1:leng1(optpar(i))),i=1,6)
        end if
c
        goto 10
c
c ... BRUTE_FORCE/FAST_FORCE
c
      else if (optpar(1)(1:2) .eq. 'BR' .or.
     +         optpar(1)(1:2) .eq. 'FA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'A'
          call textin (' Chain 1 ?',optpar(3))
        end if
        call upcase (optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Mol 2 ?',optpar(4))
        end if
        jmol = whichm (optpar(4))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'A'
          call textin (' Chain 2 ?',optpar(5))
        end if
        call upcase (optpar(5))
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - mol1 == mol2 !')
        end if
c
        if (nopt .lt. 6) then
          write (optpar (6),*) frleng
          call textin (' Fragment length ?',optpar(6))
        end if
        call str2i (optpar(6),idum,ierr)
        if (ierr .ne. 0) goto 10
        frleng = idum
c
        if (nopt .lt. 7) then
          write (optpar (7),*) frstep
          call textin (' Fragment step ?',optpar(7))
        end if
        call str2i (optpar(7),idum,ierr)
        if (ierr .ne. 0) goto 10
        frstep = idum
c
        if (nopt .lt. 8) then
          write (optpar (8),*) minmat
          call textin (' Min nr or fraction of residues to match ?',
     +      optpar(8))
        end if
        call str2r (optpar(8),xdum,ierr)
        if (ierr .ne. 0) goto 10
        if (xdum .le. 0.0) then
          call errcon ('Non-positive value')
          goto 10
        end if
        if (xdum .le. 1.0) then
          minfra = xdum
          minmat = -1
        else
          minmat = nint(xdum)
        end if
c
        if (nopt .lt. 9) optpar (9) = 'D'
        call upcase (optpar(9))
        call remspa (optpar(9))
        lbsame = (optpar(9)(1:1) .eq. 'S')
c
        if (nopt .lt. 10) optpar (10) = '1'
        call str2i (optpar(10),idum,ierr)
        if (ierr .ne. 0) goto 10
        idum = max (1, min (idum, frstep))
        f2step = idum
c
c ... do it
c
        ichn = -1
        do i=1,nchain(imol)
          if (chname(i,imol) .eq. optpar(3)(1:1)) then
            ichn = i
            goto 2943
          end if
        end do
        call errcon ('Requested chain not found')
        goto 10
 2943   continue
c
        jchn = -1
        do i=1,nchain(jmol)
          if (chname(i,jmol) .eq. optpar(5)(1:1)) then
            jchn = i
            goto 2945
          end if
        end do
        call errcon ('Requested chain not found')
        goto 10
 2945   continue
c
c ... count number of trials and refined trials
c
        ntries = 0
        nimped = 0
c
        if (optpar(1)(1:2) .eq. 'BR') goto 2940
c
c ... FAST_FORCE
c
        write (*,6026) name(imol)(1:leng1(name(imol))),
     +    chname(ichn,imol),
     +    name(jmol)(1:leng1(name(jmol))),
     +    chname(jchn,jmol),
     +    atypes(1)
        write (*,6027) frleng,frstep,f2step
c
        k1 = 0
        do i=chnptr(1,ichn,imol),chnptr(2,ichn,imol)
          if (atmnam(i,imol) .eq. atypes(1)) then
            k1 = k1 + 1
            buffz (1,k1) = atmxyz (1,i,imol)
            buffz (2,k1) = atmxyz (2,i,imol)
            buffz (3,k1) = atmxyz (3,i,imol)
          end if
        end do
c
        k2 = 0
        do j=chnptr(1,jchn,jmol),chnptr(2,jchn,jmol)
          if (atmnam(j,jmol) .eq. atypes(1)) then
            k2 = k2 + 1
            buffz (1,k1+k2) = atmxyz (1,j,jmol)
            buffz (2,k1+k2) = atmxyz (2,j,jmol)
            buffz (3,k1+k2) = atmxyz (3,j,jmol)
          end if
        end do
c
        call jvalut (' Central atoms mol 1 :',1,k1)
        call jvalut (' Central atoms mol 2 :',1,k2)
        if (k1 .lt. frleng .or. k2 .lt. frleng) then
          call errcon (' Not enough atoms')
          goto 10
        end if
c
        if (minmat .lt. 0) then
          write (*,6031) minfra
          minmat = nint (minfra*float(min(k1,k2)))
        end if
        minmat = max (3,minmat)
        write (*,6025) minmat
c
        nmlow  = -1
        rmslow = 999.99
c
        do i=1,12
          rtbest (i) = 0.0
        end do
        rtbest (1) = 1.0
        rtbest (5) = 1.0
        rtbest (9) = 1.0
c
        zone3 = chname(ichn,imol)//'*'
        zone4 = chname(jchn,jmol)//'*'
c
        mtries = (((k1-frleng+1)/frstep)+1) *
     +           (((k2-frleng+1)/f2step)+1)
        call jvalut (' Max number of trials :',1,mtries)
        if (mtries .gt. 10000) then
          call prompt (' WARNING - Many trials to do !')
          call prompt (
     +      ' ... Be patient or re-do with bigger step sizes !')
        end if
        mtr1 = mtries/10
        mtr2 = mtr1
c
        do i=1,k1-frleng+1,frstep
c
          j1 = 1
          j2 = k2-frleng+1
          if (lbsame) then
            j1 = i
            j2 = i
          end if
c
          do j=j1,j2,f2step
c
            call lsqgjk (buffz(1,i),buffz(1,k1+j),frleng,
     +                   xdum,rtdum,ierr)
            ntries = ntries + 1
c
            if (ntries .eq. mtr2) then
              if (mtries .gt. 10000) then
                xtr = 100.0 * float(ntries) / float(mtries)
                call fvalut (' Trials done  (%) :',1,xtr)
                mtr2 = mtr2 + mtr1
              end if
            end if
c
            if (ierr .eq. 0) then
              if (xdum .lt. 10.0) then
                nimped = nimped + 1
                call copyrt (rtlsq(1,imol,jmol),rtdum)
                rmsd (imol,jmol) = 999.99
                call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
                if (nmatch(imol,jmol) .gt. nmlow .and.
     +            ierr .eq. 0) then
                  nmlow = nmatch(imol,jmol)
                  rmslow = rmsd(imol,jmol)
                  call copyrt (rtbest,rtdum)
c
                  call jvalut (' Max match so far :',1,nmlow)
                  call fvalut (' RMSD (A)         :',1,rmslow)
                  if (nmlow .ge. minmat) goto 2930
c
                end if
              end if
            end if
          end do
        end do
c
 2930   continue
c
        write (*,*)
        call jvalut (' Number of trials :',1,ntries)
        call jvalut (' Number IMproved  :',1,nimped)
c
        write (*,*)
        call jvalut (' Max match :',1,nmlow)
        call fvalut (' RMSD (A)  :',1,rmslow)
c
        if (nmlow .le. 0) then
          call errcon ('Could not align; operator etc. reset !')
          call initrt (imol,jmol)
          goto 10
        end if
c
        call prompt (' Regenerating best alignment ...')
c
        call copyrt (rtlsq(1,imol,jmol),rtbest)
        call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
c
        write (*,6010) nmatch(imol,jmol),rmsd(imol,jmol),
     +    simind(imol,jmol),matchi(imol,jmol),
     +    cripp(imol,jmol),rrmsd(imol,jmol),
     +    normsd(imol,jmol),
     +    sas1(imol,jmol),sas2(imol,jmol),
     +    sas3(imol,jmol),sas4(imol,jmol),
     +    rmsdna(imol,jmol),
     +    rmsb(imol,jmol),corb(imol,jmol),
     +    rtlsq(1,imol,jmol),rtlsq(4,imol,jmol),
     +    rtlsq(7,imol,jmol),rtlsq(2,imol,jmol),
     +    rtlsq(5,imol,jmol),rtlsq(8,imol,jmol),
     +    rtlsq(3,imol,jmol),rtlsq(6,imol,jmol),
     +    rtlsq(9,imol,jmol),rtlsq(10,imol,jmol),
     +    rtlsq(11,imol,jmol),rtlsq(12,imol,jmol)
c
        write (last(imol,jmol),'(99(a,1x))',err=10)
     +    (optpar(i)(1:leng1(optpar(i))),i=1,10)
c
        goto 10
c
c ... BRUTE_FORCE
c
 2940   continue
        write (*,6028) name(imol)(1:leng1(name(imol))),
     +    chname(ichn,imol),
     +    name(jmol)(1:leng1(name(jmol))),
     +    chname(jchn,jmol),
     +    (atypes(i),i=1,natype)
        write (*,6029) bcutlo,bcuthi
        write (*,6027) frleng,frstep,f2step
c
 6026 format (' Fast-force fit of  ',a,1x,a/
     +        ' And                ',a,1x,a/
     +        ' Atom type      |',15(a4,'|'))
 6028 format (' Brute-force fit of ',a,1x,a/
     +        ' And                ',a,1x,a/
     +        ' Atom types     |',15(a4,'|'))
 6029 format (' B-factor range used  ',f8.2,' - ',f8.2,' A2')
 6027 format (' Fragment length      ',i8/
     +        ' Fragment step size   ',i8/
     +        ' Sliding step size    ',i8)
 6025 format (' Min matched residues ',i8)
 6031 format (' Min fraction matched ',f8.2)
c
        nmlow  = -1
        rmslow = 999.99
        i1low  = -1
        j1low  = -1
c
c ... find start and end residue numbers, but ignore residues
c     with negative or zero residue number
c
        i1 = iresid ( chnptr(1,ichn,imol) , imol )
        i1 = max (1, i1)
        i2 = iresid ( chnptr(2,ichn,imol) , imol )
        i2 = max (1, i2)
        j1 = iresid ( chnptr(1,jchn,jmol) , jmol )
        j1 = max (1, j1)
        j2 = iresid ( chnptr(2,jchn,jmol) , jmol )
        j2 = max (1, j2)
c
c ... with HETATMs, lower numbered residues may come after higher
c     numbered ones; therefore, loop over all atoms in the chains
c     to find the highest residue number
c
        do i=chnptr(1,ichn,imol),chnptr(2,ichn,imol)
          i2 = max (i2,iresid(i,imol))
        end do
        do j=chnptr(1,jchn,jmol),chnptr(2,jchn,jmol)
          j2 = max (j2,iresid(j,jmol))
        end do
c
        write (zone1,'(a1,i6,a1,i6)') chname(ichn,imol),i1,'-',i2
        call remspa (zone1)
        write (zone2,'(a1,i6,a1,i6)') chname(jchn,jmol),j1,'-',j2
        call remspa (zone2)
        call textut (' Mol 1 zone to try :',zone1)
        call textut (' Mol 2 zone to try :',zone2)
        zone3 = chname(ichn,imol)//'*'
        zone4 = chname(jchn,jmol)//'*'
c
        if (minmat .lt. 0) then
          write (*,6031) minfra
          k1 = i2-i1+1
          k2 = j2-j1+1
          minmat = nint (minfra*float(min(k1,k2)))
        end if
        minmat = max (3,minmat)
        write (*,6025) minmat
c
        write (*,*)
c
        do i=i1,i2-frleng+1,frstep
          ii = i + frleng - 1
c
          write (zone1,'(a1,i6,a1,i6)') chname(ichn,imol),i,'-',ii
          call remspa (zone1)
          call textut (' Try zone :',zone1)
c
c ... check that all residues exist and have a central atom !
c
          do j=i,ii
c
            call getptr (imol,chname(ichn,imol),j,atypes(1),
     +        -1,chname(jchn,jmol),jdum,iptr,jptr,ierr)
c
c      print *,' => ',imol,chname(ichn,imol),j,atypes(1),
c     +        -1,chname(jchn,jmol),jdum,iptr,jptr,ierr
c
            if (ierr .ne. 0) then
              call prompt (' Skip - missing residue(s) in zone')
              goto 2948
            end if
          end do
c
          k1 = j1
          k2 = j2-frleng+1
          if (lbsame) then
            k1 = i1
            k2 = i1
          end if
c
          do j=k1,k2,f2step
c
c ... check if last residue has a central atom (otherwise skip)
c
            jj = j + frleng - 1
            call getptr (-1,chname(ichn,imol),i,atypes(1),
     +        jmol,chname(jchn,jmol),jj,iptr,jptr,ierr)
            if (ierr .ne. 0) goto 2946
c
            write (zone2,'(a1,i6)') chname(jchn,jmol),j
            call remspa (zone2)
            rmsd (imol,jmol) = 999.99
ccc            print *,' Mol 1 ',zone1(1:10),' - Mol 2 ',zone2(1:10)
            call lsqexp (imol,zone1,jmol,zone2,1,.false.,ierr)
            ntries = ntries + 1
c
            if (rmsd(imol,jmol) .lt. 10.0 .and.
     +          nmatch(imol,jmol) .gt. 3 .and.
     +          ierr .eq. 0) then
              rmsd (imol,jmol) = 999.99
              call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
              nimped = nimped + 1
              if (nmatch(imol,jmol) .gt. nmlow .and.
     +            ierr .eq. 0) then
                nmlow = nmatch(imol,jmol)
                rmslow = rmsd(imol,jmol)
                i1low = i
                j1low = j
c
                call jvalut (' Max match so far :',1,nmlow)
                call fvalut (' RMSD (A)         :',1,rmslow)
                if (nmlow .ge. minmat) goto 2950
c
              end if
            end if
c
 2946       continue
c
          end do
c
 2948     continue
c
        end do
c
 2950   continue
c
        write (*,*)
        call jvalut (' Number of trials :',1,ntries)
        call jvalut (' Number IMproved  :',1,nimped)
c
        write (*,*)
        call jvalut (' Max match :',1,nmlow)
        call fvalut (' RMSD (A)  :',1,rmslow)
        call jvalut (' Mol 1 res :',1,i1low)
        call jvalut (' Mol 2 res :',1,j1low)
c
        if (nmlow .le. 0) then
          call errcon ('Could not align; operator etc. reset !')
          call initrt (imol,jmol)
          goto 10
        end if
c
        call prompt (' Regenerating best alignment ...')
c
        i = i1low
        j = j1low
        ii = i + frleng - 1
        write (zone1,'(a1,i6,a1,i6)') chname(ichn,imol),i,'-',ii
        call remspa (zone1)
        write (zone2,'(a1,i6)') chname(jchn,jmol),j
        call remspa (zone2)
        call lsqexp (imol,zone1,jmol,zone2,1,.false.,ierr)
        call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
c
        write (*,6010) nmatch(imol,jmol),rmsd(imol,jmol),
     +    simind(imol,jmol),matchi(imol,jmol),
     +    cripp(imol,jmol),rrmsd(imol,jmol),
     +    normsd(imol,jmol),
     +    sas1(imol,jmol),sas2(imol,jmol),
     +    sas3(imol,jmol),sas4(imol,jmol),
     +    rmsdna(imol,jmol),
     +    rmsb(imol,jmol),corb(imol,jmol),
     +    rtlsq(1,imol,jmol),rtlsq(4,imol,jmol),
     +    rtlsq(7,imol,jmol),rtlsq(2,imol,jmol),
     +    rtlsq(5,imol,jmol),rtlsq(8,imol,jmol),
     +    rtlsq(3,imol,jmol),rtlsq(6,imol,jmol),
     +    rtlsq(9,imol,jmol),rtlsq(10,imol,jmol),
     +    rtlsq(11,imol,jmol),rtlsq(12,imol,jmol)
c
        write (last(imol,jmol),'(99(a,1x))',err=10)
     +    (optpar(i)(1:leng1(optpar(i))),i=1,10)
c
        goto 10
c
c ... SIMILARITY_PLOT/LESK_PLOT
c
      else if (optpar(1)(1:2) .eq. 'SI' .or.
     +         optpar(1)(1:2) .eq. 'LE') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = 'A'
          call textin (' Chain 1 ?',optpar(3))
        end if
        call upcase (optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = prev
          call textin (' Mol 2 ?',optpar(4))
        end if
        jmol = whichm (optpar(4))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'A'
          call textin (' Chain 2 ?',optpar(5))
        end if
        call upcase (optpar(5))
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - mol1 == mol2 !')
        end if
c
        if (nopt .lt. 6) then
          if (optpar(1)(1:2) .eq. 'SI') then
            write (optpar (6),'(9a)')
     +        name(imol)(1:leng1(name(imol))),'_',
     +        name(jmol)(1:leng1(name(jmol))),'_simil.plt'
          else
            write (optpar (6),'(9a)')
     +        name(imol)(1:leng1(name(imol))),'_',
     +        name(jmol)(1:leng1(name(jmol))),'_lesk.plt'
          end if
          call locase (optpar(6))
          call textin (' Plot file name ?',optpar(6))
        end if
        call remspa (optpar(6))
c
        if (optpar(1)(1:2) .eq. 'SI') then
c
          if (nopt .lt. 7) then
            write (optpar (7),*) sistart
          end if
          call str2r (optpar(7),xdum,ierr)
          if (ierr .ne. 0) goto 10
          sistart = xdum
c
          if (nopt .lt. 8) then
            write (optpar (8),*) siend
          end if
          call str2r (optpar(8),xdum,ierr)
          if (ierr .ne. 0) goto 10
          siend = xdum
c
          if (nopt .lt. 9) then
            write (optpar (9),*) sistep
          end if
          call str2r (optpar(9),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (abs(xdum) .le. 0.01) then
            call errcon (' Step size too small')
            goto 10
          end if
          sistep = xdum
c
          sistart = abs (sistart)
          siend = abs(siend)
          call rlohi (sistart,siend)
          sistep = abs (sistep)
c
        end if
c
c ... do it
c
        ichn = -1
        do i=1,nchain(imol)
          if (chname(i,imol) .eq. optpar(3)(1:1)) then
            ichn = i
            goto 2843
          end if
        end do
        call errcon ('Requested chain not found')
        goto 10
 2843   continue
c
        jchn = -1
        do i=1,nchain(jmol)
          if (chname(i,jmol) .eq. optpar(5)(1:1)) then
            jchn = i
            goto 2845
          end if
        end do
        call errcon ('Requested chain not found')
        goto 10
 2845   continue
c
        if (optpar(1)(1:2) .eq. 'SI') then
c
          write (*,6128) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      (atypes(i),i=1,natype)
          write (*,6029) bcutlo,bcuthi
          idum = 1 + nint ((siend-sistart)/sistep)
          write (*,6127) sistart,siend,sistep,idum
c
          if (idum .lt. 3) then
            call errcon ('Too few steps !')
            goto 10
          end if
c
          if (idum .gt. maxsis) then
            call errcon ('Too many steps !')
            call jvalut (' Max nr of steps :',1,maxsis)
            goto 10
          end if
c
        else
          write (*,6129) name(imol)(1:leng1(name(imol))),
     +      chname(ichn,imol),
     +      name(jmol)(1:leng1(name(jmol))),
     +      chname(jchn,jmol),
     +      (atypes(i),i=1,natype)
        end if
c
 6129 format (' Lesk plot of ',a,1x,a/
     +        ' And          ',a,1x,a/
     +        ' Atom types   |',15(a4,'|'))
 6128 format (' Similarity plot of ',a,1x,a/
     +        ' And                ',a,1x,a/
     +        ' Atom types         |',15(a4,'|'))
 6127 format (' Cut-off start ',f6.2/
     +        ' Cut-off end   ',f6.2/
     +        ' Cut-off step  ',f6.2/
     +        ' Nr of steps   ',i6)
 6126 format (' # ',i4,' Cut-off ',f8.3,' Nmatch ',i8,
     +        ' RMSD ',f8.3,' A')
c
        zone3 = chname(ichn,imol)//'*'
        zone4 = chname(jchn,jmol)//'*'
        write (*,*)
c
        xdum = dismax
        do i=1,12
          rtdum (i) = rtlsq (i,imol,jmol)
        end do
        qqq1 = simind(imol,jmol)
        qqq2 = rmsd(imol,jmol)
        qqq3 = matchi(imol,jmol)
        qqq4 = rmsb(imol,jmol)
        qqq5 = cripp(imol,jmol)
        qqq6 = corb(imol,jmol)
        qqq7 = rrmsd(imol,jmol)
        qqq8 = normsd(imol,jmol)
        qqq9 = rmsdna(imol,jmol)
        qqq10 = sas1(imol,jmol)
        qqq11 = sas2(imol,jmol)
        qqq12 = sas3(imol,jmol)
        qqq13 = sas4(imol,jmol)
        ii = nmatch (imol,jmol)
c
        k = 0
        xmax = 0
        ymax = 0
c
        if (optpar(1)(1:2) .eq. 'SI') then
c
          do i=1,idum
            dismax = sistart + float(i-1)*sistep
            rmsd (imol,jmol) = 999.99
            nmatch (imol,jmol) = -1
            do j=1,12
              rtlsq(j,imol,jmol) = rtdum(j)
            end do
            call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
            write (*,6126) i,dismax,nmatch(imol,jmol),rmsd(imol,jmol)
            if (nmatch(imol,jmol) .gt. 0 .and. ierr .eq. 0) then
              if (k .ge. 1) then
                if (nmatch(imol,jmol) .eq. sinali(k) .and.
     +              rmsd(imol,jmol) .eq. sirmsd (k)) goto 2828
              end if
              if (k .ge. maxsis) goto 2828
c
              k = k + 1
              sinali (k) = nmatch(imol,jmol)
              sirmsd (k) = rmsd(imol,jmol)
              xmax = max (xmax,float(nmatch(imol,jmol)))
              ymax = max (ymax,rmsd(imol,jmol))
c
 2828         continue
            end if
          end do
c
        else
c
          nplot = 0
c
          dismax = 3.0
          call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
          write (*,6126) 0,dismax,nmatch(imol,jmol),rmsd(imol,jmol)
          dismax = 4.0
          call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
          write (*,6126) 0,dismax,nmatch(imol,jmol),rmsd(imol,jmol)
          dismax = 5.0
          call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
          write (*,6126) 0,dismax,nmatch(imol,jmol),rmsd(imol,jmol)
          dismax = 6.0
          call lsqimp (imol,zone3,jmol,zone4,.false.,ierr)
          write (*,6126) 0,dismax,nmatch(imol,jmol),rmsd(imol,jmol)
c
          if (nmatch(imol,jmol) .gt. 0 .and. ierr .eq. 0) then
            nplot = 1
            sinali (nplot) = nmatch(imol,jmol)
            sirmsd (nplot) = rmsd(imol,jmol)
            xmax = max (xmax,float(nmatch(imol,jmol)))
            ymax = max (ymax,rmsd(imol,jmol))
          end if
c
          m = 0
          do i=1,nuse1
            if (sptrij(i) .gt. 0) then
              m = m + 1
              j = ptri (i)
              k = ptrj (sptrij (i))
              l = sptrij (i)
              buffl (1,m) = buffi(1,i)
              buffm (1,m) = buffj(1,l)
              buffl (2,m) = buffi(2,i)
              buffm (2,m) = buffj(2,l)
              buffl (3,m) = buffi(3,i)
              buffm (3,m) = buffj(3,l)
            end if
          end do
          nmat = m
c
 2824     continue
c
c ... apply current operator
c
          call vecrtv (buffm,buffk,nmat,rtlsq(1,imol,jmol),
     +                 rtlsq(10,imol,jmol))
c
c ... find worst-fitting pair
c
          j = 1
          qqq = -999.99
          do i=1,nmat
            xbest = (buffl(1,i)-buffk(1,i)) ** 2 +
     +              (buffl(2,i)-buffk(2,i)) ** 2 +
     +              (buffl(3,i)-buffk(3,i)) ** 2
            if (xbest .gt. qqq) then
              j = i
              qqq = xbest
            end if
          end do
c
          do k=j,nmat-1
            buffl(1,k) = buffl(1,k+1)
            buffl(2,k) = buffl(2,k+1)
            buffl(3,k) = buffl(3,k+1)
            buffm(1,k) = buffm(1,k+1)
            buffm(2,k) = buffm(2,k+1)
            buffm(3,k) = buffm(3,k+1)
          end do
c
          nmat = nmat - 1
          call lsqgjk (buffl,buffm,nmat,qqq,rtlsq(1,imol,jmol),ierr)
c
          if (nmat .gt. 0 .and. ierr .eq. 0) then
            nplot = nplot + 1
            sinali (nplot) = nmat
            sirmsd (nplot) = qqq
            xmax = max (xmax,float(nmat))
            ymax = max (ymax,qqq)
          end if
c
          if (nmat  .gt. 9 .and.
     +        qqq   .gt. 0.2 .and.
     +        nplot .lt. maxsis) goto 2824
c
          k = nplot
c
        end if
c
        dismax = xdum
        do i=1,12
          rtlsq(i,imol,jmol) = rtdum(i)
        end do
        simind(imol,jmol) = qqq1
        rmsd(imol,jmol) = qqq2
        matchi(imol,jmol) = qqq3
        rmsb(imol,jmol) = qqq4
        cripp(imol,jmol) = qqq5
        corb(imol,jmol) = qqq6
        rrmsd(imol,jmol) = qqq7
        normsd(imol,jmol) = qqq8
        rmsdna(imol,jmol) = qqq9
        sas1(imol,jmol) = qqq10
        sas2(imol,jmol) = qqq11
        sas3(imol,jmol) = qqq12
        sas4(imol,jmol) = qqq13
        nmatch (imol,jmol) = ii
c
        call ivalut (' Nr of unique valid plot points :',1,k)
        if (k .lt. 2) then
          call errcon ('Not enough valid points for plot')
          goto 10
        end if
c
        close (iunit)
        call xopxua (iunit,optpar(6),linter,ierr)
        if (ierr .ne. 0) goto 10
c
        call stamp (line)
        write (iunit,'(a,a)',err=2899) '! ',line(1:leng1(line))
        write (iunit,'(a)',err=2899)
     +    'XLABEL Number of superimposed residues'
        write (iunit,'(a)',err=2899)
     +    'YLABEL RMSD of superimposed residues'
c
        if (optpar(1)(1:2) .eq. 'SI') then
          write (iunit,'(a)',err=2899)
     +      'REMARK Similarity plot (Sanchez & Sali, 1997)'
        else
          write (iunit,'(a)',err=2899)
     +      'REMARK Lesk plot (Irving, Whisstock & Lesk, 2001)'
        end if
c
        write (iunit,'(9a)',err=2899) 'REMARK',
     +    ' Compared mol ',name(imol)
     +    (1:leng1(name(imol))),' = file ',
     +    file(imol)(1:leng1(file(imol)))
        write (iunit,'(9a)',err=2899) 'REMARK',
     +    '      and mol ',name(jmol)
     +    (1:leng1(name(jmol))),' = file ',
     +    file(jmol)(1:leng1(file(jmol)))
c
        write (iunit,'(a,1x,i6)',err=2899) 'NPOINT',k
        write (iunit,'(a)',err=2899) 'COLOUR 4'
        write (iunit,'(a,4f8.2)',err=2899) 'XYVIEW',0.0,
     +    1.1*xmax,0.0,1.1*ymax
        write (iunit,'(a)',err=2899) 'XVALUE *'
        write (iunit,'(9i8)',err=2899) (sinali(i),i=1,k)
        write (iunit,'(a)',err=2899) 'YVALUE *'
        write (iunit,'(9f8.3)',err=2899) (sirmsd(i),i=1,k)
        write (iunit,'(a)',err=2899) 'END'
        call prompt (' Plot file written')
        close (iunit)
c
        goto 10
c
 2899   continue
        call errcon ('While writing plot file')
        close (iunit)
c
        goto 10
c
c ... WATERS
c
      else if (optpar(1)(1:2) .eq. 'WA') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Mol 2 ?',optpar(3))
        end if
        jmol = whichm (optpar(3))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - mol1 == mol2 !')
        end if
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f6.2)') cutwat
          call remspa (optpar(4))
          call textin (' Distance cut-off (A) ?',optpar(4))
        end if
        call str2r (optpar(4),xdum,ierr)
        if (ierr .ne. 0) goto 10
        if (xdum .le. 0.0) then
          call errcon ('Cut-off must be positive')
          goto 10
        end if
        cutwat = xdum
c
        if (nopt .lt. 5) then
          optpar (5) = name(imol)//'_'//name(jmol)//
     +      '_dist_db.plt'
          call remspa (optpar(5))
          call locase (optpar(5))
          call textin (' Plot file ?',optpar(5))
        end if
        call remspa (optpar(5))
c
        call waters (imol,jmol,optpar(5))
c
c ... HISTO_DISTO
c
      else if (optpar(1)(1:2) .eq. 'HI') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Mol 2 ?',optpar(3))
        end if
        jmol = whichm (optpar(3))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (imol .eq. jmol) then
          call prompt (' WARNING - mol1 == mol2 !')
        end if
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f6.2)') cuthis
          call remspa (optpar(4))
          call textin (' Distance cut-off (A) ?',optpar(4))
        end if
        call str2r (optpar(4),xdum,ierr)
        if (ierr .ne. 0) goto 10
        if (xdum .le. 0.0) then
          call errcon ('Cut-off must be positive')
          goto 10
        end if
        cuthis = xdum
c
        if (nopt .lt. 5) then
          write (optpar(5),'(f6.2)') 0.10
          call remspa (optpar(5))
          call textin (' Bin size (A) ?',optpar(5))
        end if
        call str2r (optpar(5),xdum,ierr)
        if (ierr .ne. 0) goto 10
        if (xdum .le. 0.0) then
          call errcon ('Bin size must be positive')
          goto 10
        end if
c
        call hisdis (imol,jmol,xdum)
c
c ... ALTER
c
      else if (optpar(1)(1:2) .eq. 'AL') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'CHain'
          call textin (' Alter option ?',optpar(2))
          nopt = 2
        end if
c
        call upcase (optpar(2))
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Mol ?',optpar(3))
        end if
        imol = whichm (optpar(3))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (3)
c
        if (optpar(2)(1:2) .eq. 'CH') then
c
          if (nopt .lt. 4) then
            optpar (4) = 'X'
            call textin (' Chain ID to alter ?',optpar(4))
          end if
          call remspa (optpar(4))
          call upcase (optpar(4))
c
          if (nopt .lt. 5) then
            optpar (5) = 'X'
            call textin (' New chain ID ?',optpar(5))
          end if
          call remspa (optpar(5))
          call upcase (optpar(5))
c
          call textut (' Chain ID to alter :',optpar(4)(1:1))
          call textut (' New chain ID      :',optpar(5)(1:1))
c
          k = 0
          do j=1,nchain(imol)
            if (chname(j,imol) .eq. optpar(4)(1:1)) then
              do i=chnptr(1,j,imol),chnptr(2,j,imol)
                achain(i,imol) = optpar(5)(1:1)
                k = k + 1
              end do
              chname(j,imol) = optpar(5)(1:1)
            end if
          end do
          call jvalut (' Nr of atoms changed :',1,k)
c
        else if (optpar(2)(1:2) .eq. 'SE') then
c
          if (nopt .lt. 4) then
            optpar (4) = 'XXXX'
            call textin (' Segment ID to alter ?',optpar(4))
          end if
          call remspa (optpar(4))
          call upcase (optpar(4))
c
          if (nopt .lt. 5) then
            optpar (5) = 'XXXX'
            call textin (' New segment ID ?',optpar(5))
          end if
          call remspa (optpar(5))
          call upcase (optpar(5))
c
          call textut (' Segment ID to alter :',optpar(4)(1:4))
          call textut (' New segment ID      :',optpar(5)(1:4))
          k = 0
          do i=1,natoms(imol)
            if (optpar(4)(1:4) .eq. axplor(i,imol)) then
              axplor(i,imol) = optpar(5)(1:4)
              k = k + 1
            end if
          end do
          call jvalut (' Nr of atoms changed :',1,k)
c
        else if (optpar(2)(1:2) .eq. 'FO') then
c
          if (nopt .lt. 4) then
            optpar (4) = 'X'
            call textin (
     +        ' Chain to set segment ID for ?',optpar(4))
          end if
          call remspa (optpar(4))
          call upcase (optpar(4))
c
          if (nopt .lt. 5) then
            optpar (5) = 'XXXX'
            call textin (' New segment ID ?',optpar(5))
          end if
          call remspa (optpar(5))
          call upcase (optpar(5))
c
          call textut (' Chain to alter :',optpar(4)(1:1))
          call textut (' New segment ID :',optpar(5)(1:4))
c
          k = 0
          do j=1,nchain(imol)
            if (chname(j,imol) .eq. optpar(4)(1:1)) then
              do i=chnptr(1,j,imol),chnptr(2,j,imol)
                axplor(i,imol) = optpar(5)(1:4)
                k = k + 1
              end do
            end if
          end do
          call jvalut (' Nr of atoms changed :',1,k)
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          if (nopt .lt. 4) then
            optpar (4) = 'X'
            call textin (
     +        ' Chain to renumber ?',optpar(4))
          end if
          call remspa (optpar(4))
          call upcase (optpar(4))
c
          inew = 1
          if (nopt .ge. 5) then
            call str2i (optpar(5),idum,ierr)
            if (ierr .eq. 0) inew = idum
          end if
c
          call textut (' Chain to renumber :',optpar(4)(1:1))
          call ivalut (' First new residue :',1,idum)
c
          k = 0
          inew = inew - 1
          ires = -99999
          do j=1,nchain(imol)
            if (chname(j,imol) .eq. optpar(4)(1:1)) then
              do i=chnptr(1,j,imol),chnptr(2,j,imol)
                if (iresid(i,imol) .ne. ires) then
                  inew = inew + 1
                  ires = iresid(i,imol)
                end if
                iresid(i,imol) = inew
                k = k + 1
              end do
            end if
          end do
          call jvalut (' Nr of atoms changed :',1,k)
          call jvalut (' Last new residue nr :',1,inew)
c
        else if (optpar(2)(1:2) .eq. 'SA') then
c
          if (nopt .lt. 4) then
            optpar (4) = 'X'
            call textin (
     +        ' Chain to set segment ID for ?',optpar(4))
          end if
          call remspa (optpar(4))
          call upcase (optpar(4))
c
          optpar(5) = optpar(4)(1:1)//'   '
          call textut (' Chain to alter :',optpar(4)(1:1))
          call textut (' New segment ID :',optpar(5)(1:4))
c
          k = 0
          do j=1,nchain(imol)
            if (chname(j,imol) .eq. optpar(4)(1:1)) then
              do i=chnptr(1,j,imol),chnptr(2,j,imol)
                axplor(i,imol) = optpar(5)(1:4)
                k = k + 1
              end do
            end if
          end do
          call jvalut (' Nr of atoms changed :',1,k)
c
        else
          call errcon ('Invalid ALter option')
          goto 10
        end if
c
c ... VRML
c
      else if (optpar(1)(1:2) .eq. 'VR') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'SEtup'
          call textin (' Sub-command ?',optpar(2))
          nopt = 2
        end if
c
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'IN') then
c
c ... VRML INIT
c
          if (nopt .ge. 3) vrfile = optpar(3)
c
          if (lvrml) then
            call xvrml_close ()
            lvrml = .false.
          end if
c
          call textut (' Open VRML file :',vrfile)
          call xopxua (ivrml,vrfile,linter,ierr)
          if (ierr .ne. 0) then
            call errcon (' Could not open VRML file')
            goto 10
          end if
c
          call xvrml_open (ivrml,rgbbg1,rgbbg2,rgbbg3)
          lvrml = .true.
          call xvrml_colour (rgbfg1,rgbfg2,rgbfg3)
c
        else if (optpar(2)(1:2) .eq. 'AD') then
c
c ... VRML ADD
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Mol ?',optpar(3))
          end if
          imol = whichm (optpar(3))
c
          if (imol .le. 0 .or. imol .gt. maxmol) then
            call errcon ('Invalid mol selection')
            goto 10
          end if
c
          if (.not. incore(imol)) then
            call errcon ('Mol not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (nopt .lt. 4) optpar(4) = '*'
          call upcase (optpar(4))
c
          if (nopt .ge. 5) then
            call xvrml_rgb_name (optpar(5),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          write (*,'(99a)') ' VRML - Add mol ',name(imol),
     +      ' chain ',optpar(4)(1:1)
c
          call vrmlca (imol,optpar(4),vrdist,cavrml)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'AL') then
c
c ... VRML ALL_CHAINS
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Mol ?',optpar(3))
          end if
          imol = whichm (optpar(3))
c
          if (imol .le. 0 .or. imol .gt. maxmol) then
            call errcon ('Invalid mol selection')
            goto 10
          end if
c
          if (.not. incore(imol)) then
            call errcon ('Mol not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          call textut (' VRML traces for mol :',name(imol))
c
          do i=1,nchain(imol)
c
            call xvrml_rgb_name (chncol(i),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
c
            write (*,'(99a)') ' VRML - Add chain ',chname(i,imol),
     +        ' colour ',chncol(i)
c
            call vrmlca (imol,chname(i,imol),vrdist,cavrml)
c
            call flusho (ivrml)
c
          end do
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
c ... VRML COLOURS
c
          call xvrml_col_list ()
c
        else if (optpar(2)(1:2) .eq. 'SE') then
c
c ... VRML SETUP
c
          if (nopt .lt. 3) then
            optpar(3) = cavrml
            call textin (' Central atom type ?',optpar(3))
          end if
          cavrml = optpar(3)
          call upcase (cavrml)
          call textut (' Central atom type :',cavrml)
c
          if (nopt .lt. 4) then
            write (optpar(4),*) vrdist
            call textin (' Max central atom distance ?',optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .eq. 0) vrdist = xdum
          call fvalut (' Max central atom distance :',1,vrdist)
c
          if (nopt .lt. 5) then
            write (optpar(5),*) rgbbg1,rgbbg2,rgbbg3
            call pretty (optpar(5))
            call textin (' Background colour ?',optpar(5))
          end if
          call xvrml_rgb_name (optpar(5),rgbbg1,rgbbg2,rgbbg3)
          write (optpar(5),*) rgbbg1,rgbbg2,rgbbg3
          call pretty (optpar(5))
          call textut (' Background colour :',optpar(5))
c
          if (nopt .lt. 6) then
            write (optpar(6),*) rgbfg1,rgbfg2,rgbfg3
            call pretty (optpar(6))
            call textin (' Default colour ?',optpar(6))
          end if
          call xvrml_rgb_name (optpar(6),rgbfg1,rgbfg2,rgbfg3)
          write (optpar(6),*) rgbfg1,rgbfg2,rgbfg3
          call pretty (optpar(6))
          call textut (' Default colour :',optpar(6))
c
        end if
c
c ... OMACRO
c
      else if (optpar(1)(1:2) .eq. 'OM') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'APpend'
          call textin (' Sub-command ?',optpar(2))
          nopt = 2
        end if
c
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'CL') then
c
c ... CLOSE
c
          if (ounit .gt. 0) then
            close (ounit)
            ounit = -1
            omol = -1
            call prompt (' O macro file closed')
          else
            call errcon ('No O macro file open')
          end if
c
        else if (optpar(2)(1:2) .eq. 'DE') then
c
c ... DEFINE
c
          if (nopt .lt. 3) then
            optpar (3) = omcent
            call textin (' Central atom type ?',optpar(3))
          end if
          omcent = optpar(3)(1:4)
c
          if (nopt .lt. 4) then
            write (optpar(4),'(f8.2)') omdist
            call textin (' Max central atom dist (A) ?',optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .eq. 0) omdist = xdum
c
          if (nopt .lt. 5) then
            optpar (5) = omconn
            call textin (' Connectivity file ?',optpar(5))
          end if
          omconn = optpar(5)
c
        else if (optpar(2)(1:2) .eq. 'IN') then
c
c ... INIT
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Mol ?',optpar(3))
          end if
          imol = whichm (optpar(3))
c
          if (imol .le. 0 .or. imol .gt. maxmol) then
            call errcon ('Invalid mol selection')
            goto 10
          end if
c
          if (.not. incore(imol)) then
            call errcon ('Mol not in memory')
            goto 10
          end if
          prev = optpar (3)
c
          if (ounit .gt. 0) then
            call prompt (' Closing previous O macro file')
            close (ounit)
            ounit = -1
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = 'lsq_' // 
     +                   name(imol)(1:leng1(name(imol))) //
     +                   '.omac'
            call locase (optpar(4))
            call textin (' File name ?',optpar(4))
          end if
c
          ounit = iunit + 1
c
          call xopxua (ounit,optpar(4),linter,ierr)
          if (ierr .ne. 0) then
            call errcon ('While opening O macro file')
            ounit = -1
            omol = -1
            goto 10
          end if
c
          dummy (1) = 0.0
          dummy (2) = 0.0
          dummy (3) = 0.0
          do i=1,natoms(imol)
            dummy (1) = dummy(1) + atmxyz(1,i,imol)
            dummy (2) = dummy(2) + atmxyz(2,i,imol)
            dummy (3) = dummy(3) + atmxyz(3,i,imol)
          end do
          dummy (1) = dummy(1) / float(natoms(imol))
          dummy (2) = dummy(2) / float(natoms(imol))
          dummy (3) = dummy(3) / float(natoms(imol))
c
          call stamp (sdum)
c
          write (ounit,6020,err=3933)
     +      '! O macro',optpar(4)(1:leng1(optpar(4)))
          write (ounit,6020,err=3933)
     +      '!',sdum(1:leng1(sdum))
          write (ounit,6020,err=3933)
     +      '!'
          write (ounit,6020,err=3933)
     +      '! o_setup off off on'
          write (ounit,6020,err=3933)
     +      '!'
          write (ounit,6020,err=3933)
     +      'conn',omconn(1:leng1(omconn))
          write (ounit,'(a,3f10.2)',err=3933)
     +      'centre_xyz ',dummy(1),dummy(2),dummy(3)
          write (ounit,6020,err=3933)
     +      '!'
          write (ounit,6020,err=3933)
     +      'print ... Analysing',
     +      name(imol)(1:leng1(name(imol)))
          write (ounit,6020,err=3933)
     +      'print ... From file',
     +      file(imol)(1:leng1(file(imol)))
          write (ounit,6020,err=3933)
     +      '!'
          write (ounit,6020,err=3933)
     +      'read stereo_chem.odb'
c          write (ounit,6020,err=3933)
c     +      'sam_at_in',
c     +      file(imol)(1:leng1(file(imol))),
c     +      name(imol)(1:leng1(name(imol))),' pdb'
          write (ounit,6020,err=3933)
     +      'pdb_read',
     +      file(imol)(1:leng1(file(imol))),
     +      name(imol)(1:leng1(name(imol))),' y y'
          write (ounit,6020,err=3933)
     +      'mol',
     +      name(imol)(1:leng1(name(imol)))
          write (ounit,6020,err=3933)
     +      'pa_case atom_z 5 6 7 8 15 16',
     +      'green cyan magenta orange yellow'
c
 6746 format (4a)
 6748 format (3a,f8.2)
c
          write (ounit,6746,err=3933)
     +      'db_set_dat ',name(imol)(1:leng1(name(imol))),
     +      '_MOLECULE_CA 1 1 ',omcent(1:leng1(omcent))
          write (ounit,6748,err=3933)
     +      'db_set_dat ',name(imol)(1:leng1(name(imol))),
     +      '_MOLECULE_CA_MXDST 1 1 ',omdist
c
          write (ounit,6020,err=3933)
     +      ( 'obj c'//name(imol)(1:leng1(name(imol))) ),
     +      'ca ; end_obj'
c
          write (ounit,6020,err=3933)
     +      'sketch_setup stick smooth 0.3 8'
          write (ounit,6020,err=3933)
     +      'sketch_stick',
     +      'c'//name(imol)(1:leng1(name(imol)))
          write (ounit,6020,err=3933)
     +      'sketch_setup stick smooth 0.1 8'
c
          write (ounit,6020,err=3933)
     +      '!'
          write (ounit,6020,err=3933)
     +      'paint_colour red'
          write (ounit,6020,err=3933)
     +      '!'
c
          omol = imol
          call flusho (ounit)
          call prompt (' O macro initialised')
c
          goto 10
c
 3933     continue
          call errcon ('While writing file; closing it')
          close (ounit)
          ounit = -1
          omol = -1
c
        else if (optpar(2)(1:2) .eq. 'AP') then
c
c ... APPEND
c
          if (ounit .le. 0) then
            call errcon ('No O macro file open')
            goto 10
          end if
c
          if (omol .le. 0 .or. omol .gt. maxmol) then
            call errcon ('Comparison molecule not found')
            call ivalut (' Value of "OMOL" :',1,omol)
            goto 10
          end if
c
          if (.not. incore(omol)) then
            call errcon ('Comparison molecule not in memory')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = prev
            call textin (' Mol ?',optpar(3))
          end if
          imol = whichm (optpar(3))
c
          if (imol .le. 0 .or. imol .gt. maxmol) then
            call errcon ('Invalid mol selection')
            goto 10
          end if
c
          if (.not. incore(imol)) then
            call errcon ('Mol not in memory')
            goto 10
          end if
          prev = optpar (3)
c
c ... create operator name
c
          sdum = '.lsq_rt_' //
     +           name(imol)(1:leng1(name(imol))) //
     +           '_to_' //
     +           name(omol)(1:leng1(name(omol)))
c
          write (ounit,6020,err=3943)
     +      '!'
          write (ounit,6020,err=3943)
     +      'print =========================================='
          write (ounit,6020,err=3943)
     +      '!'
          write (ounit,6020,err=3943)
     +      'print ... Comparing',
     +      name(imol)(1:leng1(name(imol)))
          write (ounit,6020,err=3943)
     +      'print ... From file',
     +      file(imol)(1:leng1(file(imol)))
          write (ounit,6020,err=3943)
     +      '!'
          write (ounit,'(a,1x,i12)',err=3943)
     +      'print ... Nr of matched residues',
     +      nmatch(omol,imol)
c
          if (nmatch(omol,imol) .lt. 3) goto 6321
c
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... RMS distance (A)          ',
     +      rmsd(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... RMS delta B (A2)          ',
     +      rmsb(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... Similarity index          ',
     +      simind(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... Match index               ',
     +      matchi(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... Crippen RHO               ',
     +      cripp(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... Relative RMSD             ',
     +      rrmsd(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... Normalised RMSD (100) (A) ',
     +      normsd(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... SAS(1) (A)                ',
     +      sas1(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... SAS(2) (A)                ',
     +      sas2(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... SAS(3) (A)                ',
     +      sas3(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... SAS(4) (A)                ',
     +      sas4(omol,imol)
          write (ounit,'(a,1x,f12.5)',err=3943)
     +      'print ... RMSD / Nalign (A)         ',
     +      rmsdna(omol,imol)
          write (ounit,6020,err=3943)
     +      '!'
c          write (ounit,6020,err=3943)
c     +      'sam_at_in',
c     +      file(imol)(1:leng1(file(imol))),
c     +      name(imol)(1:leng1(name(imol))),' pdb'
          write (ounit,6020,err=3943)
     +      'pdb_read',
     +      file(imol)(1:leng1(file(imol))),
     +      name(imol)(1:leng1(name(imol))),' y y'
          write (ounit,6020,err=3943)
     +      'mol',
     +      name(imol)(1:leng1(name(imol)))
          write (ounit,6020,err=3943)
     +      'pa_case atom_z 5 6 7 8 15 16',
     +      'green cyan magenta orange yellow'
c
          write (ounit,6746,err=3943)
     +      'db_set_dat ',name(imol)(1:leng1(name(imol))),
     +      '_MOLECULE_CA 1 1 ',omcent(1:leng1(omcent))
          write (ounit,6748,err=3943)
     +      'db_set_dat ',name(imol)(1:leng1(name(imol))),
     +      '_MOLECULE_CA_MXDST 1 1 ',omdist
c
          write (ounit,6020,err=3943)
     +      '!'
          write (ounit,6020,err=3943)
     +      'print ... Operator in datablock ',
     +      sdum(1:leng1(sdum))
c
          write (ounit,6020,err=3943)
     +      'db_create',
     +      sdum(1:leng1(sdum)),
     +      '12 R'
          write (ounit,6020,err=3943)
     +      '!'
c
          do i=1,12
            write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +        'db_set_data',
     +        sdum(1:leng1(sdum)),
     +        i,i,rtlsq(i,omol,imol)
          end do
c
          sdum = name(imol)(1:leng1(name(imol))) //
     +           '_to_' //
     +           name(omol)(1:leng1(name(omol)))
          write (ounit,6020,err=3943)
     +      '!'
          write (ounit,6020,err=3943)
     +      'lsq_mol',
     +      sdum(1:leng1(sdum)),
     +      name(imol)(1:leng1(name(imol))),
     +      ';'
c
          write (ounit,6020,err=3943)
     +      ( 'obj c'//name(imol)(1:leng1(name(imol))) ),
     +      'ca ; end_obj'
          write (ounit,6020,err=3943)
     +      'paint_object',
     +      'c'//name(imol)(1:leng1(name(imol)))
          write (ounit,6020,err=3943)
     +      'sketch_stick',
     +      'c'//name(imol)(1:leng1(name(imol)))
c
          sdum = '.lsq_stats_' //
     +           name(imol)(1:leng1(name(imol))) //
     +           '_to_' //
     +           name(omol)(1:leng1(name(omol)))
c
          write (ounit,6020,err=3943)
     +      'print ... Statistics in datablock ',
     +      sdum(1:leng1(sdum))
c
          write (ounit,6020,err=3943)
     +      'db_create',
     +      sdum(1:leng1(sdum)),
     +      '12 R'
          write (ounit,6020,err=3943)
     +      '!'
c
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      1,1,float(nmatch(omol,imol))
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      2,2,rmsd(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      3,3,simind(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      4,4,matchi(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      5,5,cripp(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      6,6,rrmsd(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      7,7,normsd(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      8,8,rmsdna(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      9,9,sas1(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      10,10,sas2(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      11,11,sas3(omol,imol)
          write (ounit,'(a,1x,a,1x,i2,1x,i2,f15.8)',err=3943)
     +      'db_set_data',sdum(1:leng1(sdum)),
     +      12,12,sas4(omol,imol)
c
          write (ounit,6020,err=3943)
     +      '!'
          write (ounit,6020,err=3943)
     +      '! del_obj',
     +      'c'//name(imol)(1:leng1(name(imol))),
     +      ';'
          write (ounit,6020,err=3943)
     +      '! db_kill',
     +      '*'//name(imol)(1:leng1(name(imol)))//'*'
c
 6321     continue
          call flusho (ounit)
          call prompt (' O macro extended')
c
          goto 10
c
 3943     continue
          call errcon ('While writing file; closing it')
          close (ounit)
          ounit = -1
          omol = -1
c
        else if (optpar(2)(1:2) .eq. 'WR') then
c
c ... WRITE
c
          if (ounit .le. 0) then
            call errcon ('No O macro file open')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = 'bell message NEXT'
            call textin (' O command(s) ?',optpar(3))
          end if
c
          write (ounit,6020,err=3953)
     +      optpar(3)(1:leng1(optpar(3)))
c
          call flusho (ounit)
          call textut (' Written to O macro :',optpar(3))
c
          goto 10
c
 3953     continue
          call errcon ('While writing file; closing it')
          close (ounit)
          ounit = -1
          omol = -1
c
        else
          call errcon ('Invalid OMacro option')
          call textut (' Option :',optpar(2))
        end if
c
 6020 format (20(a,1x))
c
c ... SET
c
      else if (optpar(1)(1:2) .eq. 'SE') then
c
        if (nopt .lt. 2) optpar (2) = '?'
        call upcase (optpar(2))
c
        if (optpar(2)(1:1) .eq. '?') then
c
          call prompt (' Current parameters:')
          call fvalut (' (DI) Max matching distance (A) :',1,dismax)
          call fvalut (' (DE) Decay factor              :',1,decay)
          call ivalut (' (MI) Min fragment length (res) :',1,minlen)
          call ivalut (' (FR) Fragment length decay     :',1,lendec)
          call ivalut (' (MA) Max nr of improve cycles  :',1,maxcyc)
          call textut (' (OP) Criterion                 :',optcri)
          call fvalut (' (RM) RMS weight (MI only)      :',1,rmswgt)
          call textut (' (SE) Sequential hits only      :',seqcri)
          call textut (' (SH) Frameshift correction     :',ashift)
c
        else if (optpar(2)(1:2) .eq. 'SH') then
c
          if (nopt .lt. 3) then
            optpar(3) = ashift
            call textin (' Frameshift correction (ON/OFf) ?',optpar(3))
          end if
c
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3)(1:2) .eq. 'ON' .or.
     +        optpar(3)(1:2) .eq. 'OF') then
            ashift = optpar(3)(1:2)
            lshift = (ashift .eq. 'ON')
            call textut (' Frameshift correction :',ashift)
          else
            call errcon (
     +        'Invalid selection (must be ON or OFf)')
            call textut (' Value entered :',optpar(3))
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'SE') then
c
          if (nopt .lt. 3) then
            optpar(3) = seqcri
            call textin (' Sequential hits (ON/OFf) ?',optpar(3))
          end if
c
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3)(1:2) .eq. 'ON' .or.
     +        optpar(3)(1:2) .eq. 'OF') then
            seqcri = optpar(3)(1:2)
            call textut (' Sequential hits :',seqcri)
          else
            call errcon (
     +        'Invalid selection (must be ON or OFf)')
            call textut (' Value entered :',optpar(3))
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'OP') then
c
          if (nopt .lt. 3) then
            optpar(3) = optcri
            call prompt (
     +        ' Allowed values: SI/MI/RMs/NMatch/CRippen/RRmsd/')
            call prompt (
     +        '                 NRmsd/S1/S2/S3/S4')
            call textin (
     +        ' Criterion ?',
     +        optpar(3))
          end if
c
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3)(1:2) .eq. 'RM' .or.
     +        optpar(3)(1:2) .eq. 'NM' .or.
     +        optpar(3)(1:2) .eq. 'MI' .or.
     +        optpar(3)(1:2) .eq. 'SI' .or.
     +        optpar(3)(1:2) .eq. 'S1' .or.
     +        optpar(3)(1:2) .eq. 'S2' .or.
     +        optpar(3)(1:2) .eq. 'S3' .or.
     +        optpar(3)(1:2) .eq. 'S4' .or.
     +        optpar(3)(1:2) .eq. 'CR' .or.
     +        optpar(3)(1:2) .eq. 'RR' .or.
     +        optpar(3)(1:2) .eq. 'NR') then
            optcri = optpar(3)(1:2)
            call textut (' Criterion :',optcri)
          else
            call errcon ('Invalid criterion entered')
            call textut (' Value entered :',optpar(3))
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          call prompt (' Resetting program defaults')
          call setdef (0)
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          call prompt (' Setting coarse 6 A fit defaults')
          call setdef (1)
c
        else if (optpar(2)(1:2) .eq. 'NU') then
c
          call prompt (' Setting nucleic acid defaults')
          call setdef (5)
          seqcri = 'OF'
c
        else if (optpar(2)(1:2) .eq. 'IN') then
c
          call prompt (' Setting intermediate 4 A fit defaults')
          call setdef (2)
c
        else if (optpar(2)(1:2) .eq. 'FI') then
c
          call prompt (' Setting fine-tune 3 A fit defaults')
          call setdef (3)
c
        else if (optpar(2)(1:2) .eq. 'SI') then
c
          call prompt (
     +      ' Setting defaults for 2 A fit of similar molecules')
          call setdef (4)
c
        else if (optpar(2)(1:2) .eq. 'RM') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) rmswgt
            call textin (' RMS weight ?',optpar(3))
          end if
c
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .gt. 0.0 .and. xdum .lt. 1000.0) then
            rmswgt = xdum
            call fvalut (' RMS weight :',1,rmswgt)
          else
            call errcon ('Value outside range <0,1000>')
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'DI') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) dismax
            call textin (' Max matching distance (A) ?',optpar(3))
          end if
c
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .gt. 0.0 .and. xdum .lt. 1000.0) then
            dismax = xdum
            call fvalut (' Max matching distance (A) :',1,dismax)
          else
            call errcon ('Value outside range <0,1000>')
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'DE') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) decay
            call textin (' Decay factor ?',optpar(3))
          end if
c
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .gt. 0.0 .and. xdum .lt. 10.0) then
            decay = xdum
            call fvalut (' Decay factor :',1,decay)
          else
            call errcon ('Value outside range <0,10>')
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'FR') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) lendec
            call textin (' Fragment length decay (res) ?',optpar(3))
          end if
c
          call str2i (optpar(3),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .ge. -5 .and. idum .le. 5) then
            lendec = idum
            call ivalut (' Fragment length decay (res) :',1,lendec)
          else
            call errcon ('Value outside range [-5,+5]')
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'MI') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) minlen
            call textin (' Min fragment length (res) ?',optpar(3))
          end if
c
          call str2i (optpar(3),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .gt. 0 .and. idum .lt. 1000) then
            minlen = idum
            call ivalut (' Min fragment length (res) :',1,minlen)
          else
            call errcon ('Value outside range <0,1000>')
            goto 10
          end if
c
        else if (optpar(2)(1:2) .eq. 'MA') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) maxcyc
            call textin (' Max nr of improve cycles ?',optpar(3))
          end if
c
          call str2i (optpar(3),idum,ierr)
          if (ierr .ne. 0) goto 10
          if (idum .gt. 0 .and. idum .lt. 1000) then
            maxcyc = idum
            call ivalut (' Max nr of improve cycles :',1,maxcyc)
          else
            call errcon ('Value outside range <0,1000>')
            goto 10
          end if
c
        else
          call errcon ('Invalid SEt option')
          call textut (' Option :',optpar(2))
        end if
c
c ... ATOM_TYPES
c
      else if (optpar(1)(1:2) .eq. 'AT') then
c
        if (nopt .lt. 2) optpar (2) = '?'
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'CA') then
c
          natype = 1
          atypes (1) = ' CA '
c
        else if (optpar(2)(1:2) .eq. 'PH') then
c
          natype = 1
          atypes (1) = ' P  '
c
        else if (optpar(2)(1:2) .eq. 'C4') then
c
          natype = 1
          atypes (1) = ' C4*'
c
        else if (optpar(2)(1:2) .eq. 'AL') then
c
          natype = 1
          atypes (1) = 'ALL '
c
        else if (optpar(2)(1:2) .eq. 'NO') then
c
          natype = 1
          atypes (1) = 'NONH'
c
        else if (optpar(2)(1:2) .eq. 'SI') then
c
          natype = 1
          atypes (1) = 'SIDE'
c
        else if (optpar(2)(1:2) .eq. 'TR') then
c
          natype = 1
          atypes (1) = 'TRAC'
c
        else if (optpar(2)(1:2) .eq. 'MA') then
c
          natype = 3
          atypes (1) = ' CA '
          atypes (2) = ' N  '
          atypes (3) = ' C  '
c
        else if (optpar(2)(1:2) .eq. 'EX') then
c
          natype = 5
          atypes (1) = ' CA '
          atypes (2) = ' N  '
          atypes (3) = ' C  '
          atypes (4) = ' O  '
          atypes (5) = ' CB '
c
        else if (optpar(2)(1:2) .eq. 'NU') then
c
          natype = 8
          atypes (1) = ' C4*'
          atypes (2) = ' P  '
          atypes (3) = ' C1*'
          atypes (4) = ' C2*'
          atypes (5) = ' C3*'
          atypes (6) = ' O2*'
          atypes (7) = ' O3*'
          atypes (8) = ' O4*'
c
        else if (optpar(2)(1:2) .eq. 'DE') then
c
          if (nopt .lt. 3) then
            call errcon ('No atom types entered')
            goto 10
          end if
c
          natype = 0
c
          do i=3,nopt
            if (length(optpar(i)) .gt. 0) then
              call upcase (optpar(i))
              natype = natype + 1
              atypes (natype) = optpar(i)
              if (natype .ge. maxtyp) then
                call prompt (' Max nr of atom types reached')
                goto 7382
              end if
            end if
          end do
c
          if (natype .lt. 1) then
            call errcon ('No atom types entered; reset to CA')
            natype = 1
            atypes (1) = ' CA '
            goto 10
          end if
c
 7382     continue
c
        else if (optpar(2)(1:1) .ne. '?') then
          call errcon ('Invalid ATom_types option')
          call textut (' Option :',optpar(2))
        end if
c
        call ivalut (' Nr of atom types :',1,natype)
        call asciut (' Types :',natype,atypes)
c
c ... EDIT_OPERATOR or SHOW_OPERATOR or SAVE_OPERATOR
c     or OLD_O_OPERATOR
c
      else if (optpar(1)(1:2) .eq. 'ED' .or.
     +         optpar(1)(1:2) .eq. 'SH' .or.
     +         optpar(1)(1:2) .eq. 'SA' .or.
     +         optpar(1)(1:2) .eq. 'OL') then
c
        if (nopt .lt. 2) then
          optpar (2) = prev
          call textin (' Mol 1 ?',optpar(2))
        end if
        imol = whichm (optpar(2))
c
        if (imol .le. 0 .or. imol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(imol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
        prev = optpar (2)
c
        if (nopt .lt. 3) then
          optpar (3) = prev
          call textin (' Mol 2 ?',optpar(3))
        end if
        jmol = whichm (optpar(3))
c
        if (jmol .le. 0 .or. jmol .gt. maxmol) then
          call errcon ('Invalid mol selection')
          goto 10
        end if
c
        if (.not. incore(jmol)) then
          call errcon ('Mol not in memory')
          goto 10
        end if
c
        if (optpar(1)(1:2) .eq. 'SH') then
c
c ... SHOW
c
          call textut (' Operator bringing :',name(jmol))
          call textut (' on top of         :',name(imol))
c
          call textut (' Last command was  :',last(imol,jmol))
c
          write (*,6010) nmatch(imol,jmol),rmsd(imol,jmol),
     +      simind(imol,jmol),matchi(imol,jmol),
     +      cripp(imol,jmol),rrmsd(imol,jmol),
     +      normsd(imol,jmol),
     +      sas1(imol,jmol),sas2(imol,jmol),
     +      sas3(imol,jmol),sas4(imol,jmol),
     +      rmsdna(imol,jmol),
     +      rmsb(imol,jmol),corb(imol,jmol),
     +      rtlsq(1,imol,jmol),rtlsq(4,imol,jmol),
     +      rtlsq(7,imol,jmol),rtlsq(2,imol,jmol),
     +      rtlsq(5,imol,jmol),rtlsq(8,imol,jmol),
     +      rtlsq(3,imol,jmol),rtlsq(6,imol,jmol),
     +      rtlsq(9,imol,jmol),rtlsq(10,imol,jmol),
     +      rtlsq(11,imol,jmol),rtlsq(12,imol,jmol)
c
ccc          call anamat (rtlsq(1,imol,jmol),xdum)
c
          call anancs (1,rtlsq(1,imol,jmol),.true.,ierr)
c
          if (imol .eq. jmol) then
            if (rmsd(imol,jmol) .le. 0.01) then
              call prompt ('0 *GOOD* - NCS constrained')
            else if (rmsd(imol,jmol) .le. 0.2) then
              call prompt ('0 *GOOD* - NCS tightly restrained')
            else if (rmsd(imol,jmol) .le. 0.3) then
              call prompt ('0 *OKAY* - NCS restrained')
            else if (rmsd(imol,jmol) .le. 0.5) then
              call prompt ('0 *FAIR* - NCS weakly restrained ?')
            else
              call prompt ('0 *POOR* - NCS not restrained ?')
            end if
c
            if (rmsb(imol,jmol) .le. 0.01) then
              call prompt ('  *GOOD* - NCS Bs constrained')
            else if (rmsb(imol,jmol) .le. 1.0) then
              call prompt ('  *GOOD* - NCS Bs tightly restrained')
            else if (rmsb(imol,jmol) .le. 3.0) then
              call prompt ('  *OKAY* - NCS Bs restrained')
            else if (rmsb(imol,jmol) .le. 6.0) then
              call prompt ('  *FAIR* - NCS Bs weakly restrained ?')
            else
              call prompt ('  *POOR* - NCS Bs not restrained ?')
            end if
          end if
c
        else if (optpar(1)(1:2) .eq. 'ED') then
c
c ... EDIT
c
          call textut (' Operator bringing :',name(jmol))
          call textut (' on top of         :',name(imol))
c
          do i=1,12
            rtdum (i) = rtlsq(i,imol,jmol)
            if (nopt .lt. (i+3)) then
              write (optpar(i+3),*) rtdum(i)
              write (sdum,'(a,i2,a)')
     +          ' Operator element ',i,' ?'
              call textin (sdum,optpar(i+3))
            end if
c
            call str2r (optpar(i+3),xdum,ierr)
            if (ierr .ne. 0) goto 10
            if (i .le. 9) then
              if (xdum .lt. -1.0 .or. xdum .gt. 1.0) then
                call errcon (
     +            ' Matrix element outside range [-1,+1]')
                goto 10
              end if
            end if
c
            rtdum (i) = xdum
          end do
c
ccc          call anamat (rtdum,xdum)
c
          call anancs (1,rtdum,.true.,ierr)
          xdum = det3 (rtdum)
c
          if (xdum .lt. 0.9995 .or. xdum .gt. 1.0005) then
            call errcon (
     +        'Determinant outside range [0.9995,1.0005]')
            goto 10
          end if
c
          do i=1,12
            rtlsq (i,imol,jmol) = rtdum (i)
          end do
c
          nmatch (imol,jmol) = 0
          rmsd   (imol,jmol) = 999.99999
          rmsb   (imol,jmol) = 999.99999
          corb   (imol,jmol) = 999.99999
          simind (imol,jmol) = 999.99999
          matchi (imol,jmol) = 999.99999
          cripp  (imol,jmol) = 999.99999
          rrmsd  (imol,jmol) = 999.99999
          normsd (imol,jmol) = 999.99999
          rmsdna (imol,jmol) = 999.99999
          sas1   (imol,jmol) = 999.99999
          sas2   (imol,jmol) = 999.99999
          sas3   (imol,jmol) = 999.99999
          sas4   (imol,jmol) = 999.99999
c
          write (last(imol,jmol),'(99(a,1x))',err=10)
     +      (optpar(i)(1:leng1(optpar(i))),i=1,15)
c
        else if (optpar(1)(1:2) .eq. 'SA') then
c
c ... SAVE
c
          call textut (' Operator bringing :',name(jmol))
          call textut (' on top of         :',name(imol))
c
          if (nopt .lt. 4) then
            optpar (4) = 'rt_' // 
     +                   name(jmol)(1:leng1(name(jmol))) //
     +                   '_to_' //
     +                   name(imol)(1:leng1(name(imol))) //
     +                   '.odb'
            call locase (optpar(4))
            call textin (' File name ?',optpar(4))
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = '.lsq_rt_' // 
     +                   name(jmol)(1:leng1(name(jmol))) //
     +                   '_to_' //
     +                   name(imol)(1:leng1(name(imol)))
            call locase (optpar(5))
          end if
c
          call textut (' Save in file   :',optpar(4))
          call textut (' Datablock name :',optpar(5))
c
          close (iunit)
          call xopxua (iunit,optpar(4),linter,ierr)
          if (ierr .ne. 0) goto 10
c
          call stamp (sdum)
          write (iunit,'(a,a)',err=3913) '! ',sdum(1:leng1(sdum))
          write (iunit,'(a,a)',err=3913)
     +      optpar(5)(1:leng1(optpar(5))),' r 12 (3f15.7)'
          write (iunit,'(3f15.7)',err=3913)
     +      (rtlsq(i,imol,jmol),i=1,12)
c
          close (iunit)
          goto 10
c
 3913     continue
          call errcon ('While writing file')
          close (iunit)
          goto 10
c
        else if (optpar(1)(1:2) .eq. 'OL') then
c
c ... OLD
c
          call textut (' Operator bringing :',name(jmol))
          call textut (' on top of         :',name(imol))
c
          if (nopt .lt. 4) then
            optpar (4) = 'rt_' // 
     +                   name(jmol)(1:leng1(name(jmol))) //
     +                   '_to_' //
     +                   name(imol)(1:leng1(name(imol))) //
     +                   '.odb'
            call locase (optpar(4))
            call textin (' File name ?',optpar(4))
          end if
c
          close (iunit)
          call opoodb (iunit,optpar(4),parnam,partyp,
     +      idum,parfmt,ierr)
          if (ierr .ne. 0) goto 10
c
          call upcase (partyp)
          if (partyp .ne. 'R') then
            call errcon ('Not a REAL datablock')
            goto 10
          end if
c
          if (idum .ne. 12) then
            call errcon ('Datablock MUST have TWELVE values')
            goto 10
          end if
c
          read (iunit,fmt=parfmt,err=3923) (rtdum(i),i=1,12)
c
          close (iunit)
          call rvalut (' Operator :',12,rtdum)
c
ccc          call anamat (rtdum,xdum)
c
          call anancs (1,rtdum,.true.,ierr)
          xdum = det3 (rtdum)
c
          if (xdum .lt. 0.9995 .or. xdum .gt. 1.0005) then
            call errcon (
     +        'Determinant outside range [0.9995,1.0005]')
            goto 10
          end if
c
          do i=1,9
            xdum = rtdum(i)
            if (xdum .lt. -1.0 .or. xdum .gt. 1.0) then
              call errcon (
     +          ' Matrix element outside range [-1,+1]')
              goto 10
            end if
          end do
c
          call prompt (' Operator looks okay !')
          do i=1,12
            rtlsq (i,imol,jmol) = rtdum (i)
          end do
c
          nmatch (imol,jmol) = 0
          rmsd   (imol,jmol) = 999.99999
          rmsb   (imol,jmol) = 999.99999
          corb   (imol,jmol) = 999.99999
          simind (imol,jmol) = 999.99999
          matchi (imol,jmol) = 999.99999
          cripp  (imol,jmol) = 999.99999
          rrmsd  (imol,jmol) = 999.99999
          normsd (imol,jmol) = 999.99999
          rmsdna (imol,jmol) = 999.99999
          sas1   (imol,jmol) = 999.99999
          sas2   (imol,jmol) = 999.99999
          sas3   (imol,jmol) = 999.99999
          sas4   (imol,jmol) = 999.99999
c
          write (last(imol,jmol),'(99(a,1x))',err=10)
     +      (optpar(i)(1:leng1(optpar(i))),i=1,4)
c
          goto 10
c
 3923     continue
          call errcon ('While reading file')
          close (iunit)
          goto 10
c
        end if
c
c ... INVERT_NCS
c
      else if (optpar(1)(1:2) .eq. 'IN') then
c
        if (nopt .lt. 2) then
          optpar (2) = ' '
          call textin (' Input RT file ?',optpar(2))
        end if
c
        if (length(optpar(2)) .lt. 1) then
          call errcon ('No file name provided')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = ' '
          call textin (' Output RT file ?',optpar(3))
        end if
c
        if (length(optpar(3)) .lt. 1) then
          call errcon ('No file name provided')
          goto 10
        end if
c
        call rdoncs (iunit,optpar(2),ncsdum,maxncs,dumncs,ierr)
c
        if (ierr .ne. 0) then
          call errcon ('While reading RT operator file')
          goto 10
        end if
c
        call jvalut (' Nr of RT operators to invert :',1,ncsdum)
c
        call invncs (ncsdum,dumncs,duminv)
c
        call wroncs (iunit,optpar(3),ncsdum,maxncs,duminv,ierr)
c
        if (ierr .ne. 0) then
          call errcon ('While writing RT operator file')
        end if
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
 6010 format (
     +  ' The ',i6,' atoms have an RMS distance of ',f8.3,' A'/
     +  ' SI = RMS * Nmin / Nmatch             = ',f12.5/
     +  ' MI = (1+Nmatch)/{(1+W*RMS)*(1+Nmin)} = ',f12.5/
     +  ' CR = Maiorov-Crippen RHO (0-2)       = ',f12.5/
     +  ' RR = relative RMSD                   = ',f12.5/
     +  ' NR = normalised RMSD (100)           = ',f10.3,' A'/
     +  ' SAS(1) = cRMS * (100/Nmatch)         = ',f10.3,' A'/
     +  ' SAS(2) = cRMS * (100/Nmatch)^2       = ',f10.3,' A'/
     +  ' SAS(3) = cRMS * (100/Nmatch)^3       = ',f10.3,' A'/
     +  ' SAS(4) = cRMS * (100/Nmatch)^4       = ',f10.3,' A'/
     +  ' RMSD / Nalign                        = ',f10.5,' A'/
     +  ' RMS delta B for matched atoms        = ',f9.3,' A2'/
     +  ' Corr. coefficient matched atom Bs    = ',f12.3/
     +  ' Rotation     : ',3f12.8/16x,3f12.8/16x,3f12.8/
     +  ' Translation  : ',3f12.4)
c
c --- END OF MAIN EVENT LOOP
c
c ... check if all changes saved
c
 9000 continue
c
      call gkquit
c
      end
