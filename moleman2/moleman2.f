      program mole2
c
c ... MOLEMAN2 - GJ Kleywegt @ 951107
c
      include 'moleman2.incl'
c
      integer maxopt
      parameter (maxopt=50)
c
      real total,user,sys,dummy,volume,border,vrdist,rr,gg,bb
      real rgbbg1,rgbbg2,rgbbg3,rgbfg1,rgbfg2,rgbfg3,xdum,ydum
      real tramul,veclen,rcocut
c
      integer iunit,ierr,nopt,nempty,i,j,k,length,junit,kunit,idum
      integer munit,ivrml,jdum,ncolen
c
      logical xinter,ldone,linit,lvrml,lmono,lecho,ltrueb,ldist
c
      character line*256,optpar(maxopt)*80,pro*12,inimac*128
      character cavrml*4,vrmlbg*25,vrmldc*25,vrfile*128,acotyp*4
c
code ...
c
      call gkinit (prognm,vers)
c
c ... initialise history
c
      call dohist ('*INIT*',ldone)
c
c ... initialise some variables
c
      do i=1,maxopt
        optpar (i) = ' '
      end do
      nempty = 0
      nopt   = 0
c
      lecho = .false.
c
      iunit  = 10
      junit  = 11
      kunit  = 12
c
c ... user input unit (5=interactive; other=macro)
c
      munit = 5
c
      linit = .false.
c
      linter = xinter()
c
      if (linter) then
        pro='$MOLEMAN2 > '
      else
        pro=' MOLEMAN2 > '
      end if
c
      natoms = 0
      nres   = 0
      nrem   = 0
      nother = 0
      pdbfil = 'm1.pdb'
c
      do i=1,3
        cell (i) = 1.0
        cell (i+3) = 90.0
      end do
      zmol = 1
      spgrp = 'P 1'
      border = 1.0
      veclen = 3.0
      tramul = 1.0
      ltrueb = .true.
c
      ncolen = 5
      acotyp = ' CA '
c
      rcocut = 6.0
c
      lvrml  = .false.
      ivrml  = 99
      lmono  = .true.
      cavrml = ' CA '
      vrdist = 4.5
      vrmlbg = 'black'
      vrmldc = 'white'
      vrfile = 'moleman2.wrl'
      rgbbg1 = 0.0
      rgbbg2 = 0.0
      rgbbg2 = 0.0
      rgbfg1 = 1.0
      rgbfg2 = 1.0
      rgbfg2 = 1.0
      call xvrml_init ()
c
      do i=1,maxpln
        lplane (i) = .false.
        write (nplane(i),'(a,i3)') 'Plane ',i
        vplane (1,i) = 1.0
        vplane (2,i) = 1.0
        vplane (3,i) = 1.0
      end do
c
      ssenam (-1) = 'Non-protein'
      ssenam (0)  = 'Loop or turn'
      ssenam (1)  = 'Alpha helix'
      ssenam (2)  = 'Beta strand'
      ssenam (3)  = 'Left-handed helix'
c
c ... print dimensioning
c
      call prompt ('0Array dimensioning:')
c
      call prompt ('01) Library:')
      call jvalut ('    Max nr of residue types       :',
     +  1,mxrtyp)
      call jvalut ('    Max nr of atom types          :',
     +  1,mxrtat)
      call jvalut ('    Max nr of residue aliases     :',
     +  1,mxrtal)
      call jvalut ('    Nr of defined residue classes :',
     +  1,mxrtal)
c
      call prompt ('02) Molecule:')
      call jvalut ('    Max nr of atoms               :',
     +  1,maxatm)
      call jvalut ('    Max nr of residues            :',
     +  1,maxres)
      call jvalut ('    Max nr of REMARK records      :',
     +  1,maxcom)
      call jvalut ('    Max nr of other records       :',
     +  1,maxcom)
c
      call prompt ('03) Program:')
      call jvalut ('    Max buffer size               :',
     +  1,maxbuf)
      call jvalut ('    Max nr of atoms per residue   :',
     +  1,maxapr)
      call jvalut ('    Max nr of residue torsions    :',
     +  1,maxtor)
c
c ... check buffer size
c
      if (maxbuf .lt. max (512*1024, 6*maxatm, 5*maxres)) then
        call errstp ('Buffer too small; ask Gerard to recompile !')
        goto 9000
      end if
c
c ... set default values
c
      write (*,*)
      call set_defaults ()
c
c ... define some symbols
c
      write (*,*)
      nopt = 3
      optpar(1) = '&'
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
      optpar(1) = ' '
      optpar(2) = ' '
      optpar(3) = ' '
c
c ... get name of library file
c
      line = 'moleman2.lib'
      call gklibf (line)
c
      write (*,*)
      call textin (' Name of library file ?',line)
      call xopxoa (iunit,line,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening MOLEMAN2 library file')
        goto 9000
      end if
c
      call prompt ('0Reading library ...')
      call readdb (iunit,ierr)
      close (iunit)
      if (ierr .ne. 0) then
        call errcon ('While reading MOLEMAN2 library file')
        goto 9000
      end if
c
      call gkdcpu (total,user,sys)
      if (total .ge. 1.0)
     +  write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
c ... formats
c
 6000 format (/
     +  ' MOLEMAN2 commands :'//
     +  ' ? [command] (list (sub-)commands)   ',
     +              ' ! (comment)'/
     +  ' QUit                                ',
     +              ' $ shell_command'/
     +  ' & symbol value                      ',
     +              ' & ? (list symbols)'/
     +  ' @ macro_file                        ',
     +              ' ECho on_off'/
     +  ' # parameter(s) (command history)    ',
     +              ' '/
     +  ' COnstants REset                     ',
     +              ' COnstants SEt name value'/
     +  ' COnstants LIst                      ',
     +              ' STatistics'/
     +  ' GEometry_selected                   ',
     +              ' LIst_selected [which]'/
     +  ' MUlti_geometry which                ',
     +              ' LS_plane'/
     +  ' RContact_order [cut_off]            ',
     +              ' '/
     + /' REad file [format] [hydro]          ',
     +              ' WRite file [format] [which]'/,
     +  ' APpend file [format] [hydro]        ',
     +              ' SPlit file_prefix'/
     +  ' DELETE_molecule                     ',
     +              ' BOok_keeping')
c
 6001 format (/
     +  ' Commands with sub-commands:'/
     +  ' SElect  BFactor  OCcupancy  CHain    PDb '/
     +  ' XYz     ONo      DIstance   SQuence  AUto'/
     +  ' PRotein NUcleic  VRml'/
     +  ' To see sub-commands, use for instance: ? xy')
c
 6010 format (/
     +  ' SElect All                          ',
     +              ' SElect NOne'/
     +  ' SElect HYdrogen                     ',
     +              ' SElect EXhydrogen'/
     +  ' SElect OR what which                ',
     +              ' SElect ANd what which'/
     +  ' SElect NEgate                       ',
     +              ' SElect NUmeric and_or_butnot what lo hi'/
     +  ' SElect DIst_to_sel lo hi            ',
     +              ' SElect POint and_or_butnot x y z lo hi'/
     +  ' SElect BY_residue                   ',
     +              ' SElect BUtnot what which'/
     +  ' SElect ?                            ',
     +              ' ')
c
 6020 format (/
     +  ' BFactor STats [how]                 ',
     +              ' BFactor LImit lo hi'/
     +  ' BFactor PLot file [what]            ',
     +              ' BFactor SMooth'/
     +  ' BFactor BOnded                      ',
     +              ' BFactor GRoup how'/
     +  ' BFactor PRod_plus [prod] [plus]     ',
     +              ' BFactor NO_anisou'/
     +  ' BFactor SAve                        ',
     +              ' BFactor REstore'/
     +  ' BFactor PSeudo what [parameters]    ',
     +              ' BFactor SCale [min] [max]'/
     +  ' BFactor ODb [filename] [molname]    ',
     +              ' ')
c
 6030 format (/
     +  ' OCcupancy STats [how]               ',
     +              ' OCcupancy LImit lo hi'/
     +  ' OCcupancy PRod_plus [prod] [plus]   ',
     +              ' OCcupancy PLot file [what]')
c
 6040 format (/
     +  ' CHain AUto                          ',
     +              ' CHain ASk_auto'/
     +  ' CHain FRom_segid [how]              ',
     +              ' CHain TO_segid [how]'/
     +  ' CHain REname old new                ',
     +              ' CHain SEgid_rename old new'/
     +  ' CHain NAme_selection chain segid    ',
     +              ' CHain OT2_suggest')
c
 6050 format (/
     +  ' PDb HEtero option                   ',
     +              ' PDb CRystal cell1..6 zmol spacegroup'/
     +  ' PDb SSbond what [how]               ',
     +              ' PDb REmark text'/
     +  ' PDb LIst_remark                     ',
     +              ' PDb DElete_remark which'/
     +  ' PDb NAme which old new              ',
     +              ' PDb NUmber first_new'/
     +  ' PDb CHemical+charge                 ',
     +              ' PDb NO_atomic_numbers'/
     +  ' PDb INdonesia                       ',
     +              ' PDb SAnity_check'/
     +  ' PDb SEqres                          ',
     +              ' ')
c
 6060 format (/
     +  ' PRotein MC_analysis [file] [what] [c',
     +              'hain]'/
     +  ' PRotein CA_analysis [file] [what] [c',
     +              'hain]'/
     +  ' PRotein SC_analysis [file]          ',
     +              ' ')
c
 6070 format (/
     +  ' XYz FRactionalise                   ',
     +              ' XYz ORthogonalise'/
     +  ' XYz ROtate how ang1 ang2 ang3       ',
     +              ' XYz TRanslate tx ty tz'/
     +  ' XYz AXis_rotate axis_xyz angle      ',
     +              ' XYz RAndom_rotation'/
     +  ' XYz MAtrix [r11 .. r33 tx ty tz]    ',
     +              ' XYz RT file'/
     +  ' XYz CEntre_origin [use_masses]      ',
     +              ' XYz PErturb [dx] [dy] [dz] [db] [dq]'/
     +  ' XYz MIrror [xyz] [value]            ',
     +              ' XYz INvert [icx] [icy] [icz]'/
     +  ' XYZ DIstort [r11 .. r33 tx ty tz]   ',
     +              ' XYz ALign_inertia_axes')
c
 6080 format (/
     +  ' ONo RSr res_type [file] [res_nr]    ',
     +              ' ONo FIt res_type [file] [res_nr]'/
     +  ' ONo TOrs res_type [file] [res_nr]   ',
     +              ' ONo COnnect res_type [file] [res_nr]'/
     +  ' ONo CEll [file] [mode] [colour]     ',
     +              ' ONo WAter_fit_macro [file]'/
     +  ' ONo OOps_macro o_mol [file]         ',
     +              ' ONo XPlor_hydrogens'/
     +  ' ONo LS_plane_odl [file] [obj] [col] ',
     +              '[border] [nr] [name]'/
     +  ' ONo INertia_axes_odl [file] [obj] [c',
     +              'ol] [length]'/
     +  ' ONo MOlray [file] [mode] [mult] [nea',
     +              'r_far] [for_rev] [resi] [atom]'/
     +  ' ONo ANgle_ls_planes [nr1] [nr2]     ',
     +              ' ONo RIngs [file] [obj] [col]'/
     +  ' ONo DIsulfide_odl o_mol [file] [how]',
     +              ' [object]')
c
 6090 format (/
     +  ' DIstance PLot 2d_plot_file          ',
     +              ' DIstance DIstribution [bin_size]'/
     +  ' DIstance SHort [cut_off]            ',
     +              ' DIstance SElect [cut-off] [mode]'/
     +  ' DIstance LIst [lower] [upper]       ',
     +              ' DIstance CHains chain1 chain2 [cut-off]')
c
c ... to do: DIstance COnvolution 2d_plot_file [fragment_length] [atom_type]
c
 6100 format (/
     +  ' SQuence LIst how                    ',
     +              ' SQuence PIr file [seq_name] [title]'/
     +  ' SQuence GLyco_sites                 ',
     +              ' SQuence MOtif motif'/
     +  ' SQuence COunt_residue_types         ',
     +              ' SQuence EXtinction_280')
c
 6110 format (/
     +  ' AUto SPink type nr_residues         ',
     +              ' AUto BOnes type x1 y1 z1 x2 y2 z2'/
     +  ' AUto SSe filename                   ',
     +              ' ')
c
 6120 format (/
     +  ' VRml SEtup central_atom max_dist backgr_col default_col'/
     +  ' VRml INit [filename]                ',
     +              ' VRml NAmed_colours'/
     +  ' VRml COlour_selection colour        ',
     +              ' VRml CRamp_selection how'/
     +  ' VRml REset_colours_radii            ',
     +              ' VRml RAdii mult add'/
     +  ' VRml CEll [colour]                  ',
     +              ' VRml CLose_file'/
     +  ' VRml TRace [colour]                 ',
     +              ' VRml FAt_trace fatness [colour]'/
     +  ' VRml CPk radius [colour]            ',
     +              ' VRml SPhere [colour]'/
     +  ' VRml CYlinder radius [colour]       ',
     +              ' VRml STick [colour]'/
     +  ' VRml LIquorice radius [colour]      ',
     +              ' VRml BAll_cylinder radius [colour]')
c
 6130 format (/
     +  ' NUcleic DUarte_pyle [file] [what] [c',
     +              'hain]')
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6000)
      write (*,6010)
      write (*,6020)
      write (*,6030)
      write (*,6040)
      write (*,6050)
      write (*,6060)
      write (*,6070)
      write (*,6080)
      write (*,6090)
      write (*,6100)
      write (*,6110)
      write (*,6120)
      write (*,6130)
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
        call gknval ('GKMOLEMAN2',inimac,ierr)
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
        call textut (' ERROR - No valid command :',line)
        goto 100
      end if
c
c ... handle symbols
c
      call dosymb (nopt,optpar,ldone)
      if (ldone) goto 10
c
      if (optpar(1)(1:1) .eq. '?') then
        if (nopt .lt. 2) goto 100
        call upcase (optpar(2))
        if (optpar(2)(1:2) .eq. '? ') then
          write (*,6000)
          write (*,6001)
        else if (optpar(2)(1:2) .eq. 'SE') then
          write (*,6010)
        else if (optpar(2)(1:2) .eq. 'BF') then
          write (*,6020)
        else if (optpar(2)(1:2) .eq. 'OC') then
          write (*,6030)
        else if (optpar(2)(1:2) .eq. 'CH') then
          write (*,6040)
        else if (optpar(2)(1:2) .eq. 'PD') then
          write (*,6050)
        else if (optpar(2)(1:2) .eq. 'PR') then
          write (*,6060)
        else if (optpar(2)(1:2) .eq. 'XY') then
          write (*,6070)
        else if (optpar(2)(1:2) .eq. 'ON') then
          write (*,6080)
        else if (optpar(2)(1:2) .eq. 'DI') then
          write (*,6090)
        else if (optpar(2)(1:2) .eq. 'SQ') then
          write (*,6100)
        else if (optpar(2)(1:2) .eq. 'AU') then
          write (*,6110)
        else if (optpar(2)(1:2) .eq. 'VR') then
          write (*,6120)
        else if (optpar(2)(1:2) .eq. 'NU') then
          write (*,6130)
        else
          goto 100
        end if
        write (*,*)
        goto 10
      end if
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
c ... STatistics
c
      else if (optpar(1)(1:2) .eq. 'ST') then
c
        call mole_stats ()
        goto 10
c
c ... REad file [format] [hydro]
c
      else if (optpar(1)(1:2) .eq. 'RE') then
c
        if (nopt .lt. 2) then
          optpar(2) = pdbfil
          call textin (' PDB file ?',optpar(2))
        end if
c
        if (optpar(2) .eq. '?') then
          call prompt (' file   = PDB file name')
          call prompt (' format = Pdb | Alwyn')
          call prompt (' hydro  = No  | Yes')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = 'PDB'
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'NO'
        end if
c
        call read_pdb (iunit,0,optpar(2),optpar(3),optpar(4))
        write (*,*)
        call jvalut (' Nr of atoms now :',1,natoms)
        call jvalut (' Nr of residues  :',1,nres)
        call do_select ('ALL',' ',' ')
        close (iunit)
        pdbfil = optpar(2)
        ltrueb = .true.
c
        goto 10
c
c ... AUto commands
c
      else if (optpar(1)(1:2) .eq. 'AU') then
c
        if (nopt .lt. 2) then
          optpar(2) = 'LIst'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'SP') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'beta'
            call textin (' Type ?',optpar(3))
          end if
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' type    = Alpha | Beta')
            call prompt (' nr_res  = number of residues to generate')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = '10'
            call textin (' Nr of residues ?',optpar(4))
          end if
c
          call auto_spink (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'BO') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'beta'
            call textin (' Type ?',optpar(3))
          end if
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' type    = Alpha | Beta')
            call prompt (
     +        ' x1..z2  = coordinates of first and last CA atom')
            goto 10
          end if
c
          do i=1,6
            if (nopt .lt. (i+3) ) then
              optpar(i+3) = '0.0'
              call textin (' Coordinate ?',optpar(i+3))
            end if
          end do
c
          call auto_bones (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'SS') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'test.sse'
            call textin (' DEJAVU SSE file name ?',optpar(3))
          end if
c
          call auto_sse (iunit,optpar(3))
c
        else
          call errcon ('Invalid AUto option')
        end if
c
        goto 10
c
c ... COnstants commands
c
      else if (optpar(1)(1:2) .eq. 'CO') then
c
        if (nopt .lt. 2) then
          optpar(2) = 'LIst'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'RE') then
          call prompt (' Reset program constants to defaults')
          call set_defaults ()
        else if (optpar(2)(1:2) .eq. 'LI') then
          call user_defaults (0,optpar(3),optpar(4))
        else if (optpar(2)(1:2) .eq. 'SE') then
          if (nopt .lt. 3) then
            optpar(3) = 'CANCEL'
            call textin (' Name ?',optpar(3))
          end if
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' name   = parameter name')
            call prompt (' value  = new value for this parameter')
            goto 10
          end if
          if (nopt .lt. 4) then
            optpar(4) = '0.0'
            call textin (' Value ?',optpar(4))
          end if
          call user_defaults (1,optpar(3),optpar(4))
        else
          call errcon ('Invalid COnstants option')
        end if
c
        goto 10
c
      end if
c
c ... the following commands require that atoms have been read ...........
c
      if (natoms .le. 0) then
        call errcon ('No molecule in memory')
        goto 10
      end if
c
c ... WRite file [format] [hydro]
c     format = Pdb | Xplor | Ccp4
c     which  = ALl | NO_hydro | SElected | PAla | PGly | PSer | CAlpha
c
      if (optpar(1)(1:2) .eq. 'WR') then
c
        if (nopt .lt. 2) then
          optpar(2) = pdbfil
          call textin (' PDB file ?',optpar(2))
        end if
c
        if (optpar(2) .eq. '?') then
          call prompt (' file   = PDB file name')
          call prompt (' format = Pdb | Xplor | Ccp4')
          call prompt (' which  = ALl | NO_hydro | SElected | PAla |')
          call prompt ('          PGly | PSer | CAlpha')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = 'PDB'
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'ALl'
        end if
c
        call write_pdb (iunit,optpar(2),optpar(3),optpar(4))
        write (*,*)
        close (iunit)
        pdbfil = optpar(2)
c
        goto 10
c
c ... APpend file [format] [hydro]
c
      else if (optpar(1)(1:2) .eq. 'AP') then
c
        if (nopt .lt. 2) then
          optpar(2) = pdbfil
          call textin (' PDB file ?',optpar(2))
        end if
c
        if (optpar(2) .eq. '?') then
          call prompt (' file   = PDB file name')
          call prompt (' format = Pdb | Alwyn')
          call prompt (' hydro  = No  | Yes')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = 'PDB'
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'NO'
        end if
c
        call read_pdb (iunit,1,optpar(2),optpar(3),optpar(4))
        write (*,*)
        call jvalut (' Nr of atoms now :',1,natoms)
        call jvalut (' Nr of residues  :',1,nres)
        call do_select ('ALL',' ',' ')
        close (iunit)
c
        goto 10
c
c ... BOok_keeping
c
      else if (optpar(1)(1:2) .eq. 'BO') then
c
        call book_keep (.true.)
        goto 10
c
c ... DELETE
c
      else if (optpar(1)(1:6) .eq. 'DELETE') then
c
        call prompt (' ALL ATOMS AND RESIDUES DELETED !!!')
        natoms  = 0
        nres    = 0
        nrem    = 0
        nother  = 0
        nselect = 0
        do i=1,3
          cell (i) = 1.0
          cell (i+3) = 90.0
        end do
        zmol = 1
        spgrp = 'P 1'
c
        goto 10
c
c ... LS_plane
c
      else if (optpar(1)(1:2) .eq. 'LS') then
c
        call xyz_ls_plane (iunit,' ',optpar(4),optpar(5),border,
     +       -1,' ','L')
        goto 10
c
c ... GEometry_selected
c
      else if (optpar(1)(1:2) .eq. 'GE') then
c
        call geom_select ()
        goto 10
c
c ... MUlti_geometry
c
      else if (optpar(1)(1:2) .eq. 'MU') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'GLY'
          call textin (' Residue type ?',optpar(2))
        end if
        call upcase (optpar(2))
c
        if (length(optpar(2)) .lt. 1) then
          call errcon ('No residue type')
          goto 10
        end if
c
        call geom_multi (optpar(2))
        goto 10
c
c ... RContact_order
c
      else if (optpar(1)(1:2) .eq. 'RC') then
c
        if (nopt .ge. 2) then
          call str2r (optpar(2),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .lt. 0.0 .or. xdum .gt. 100000.0) then
            call errcon ('Invalid value for cut-off')
            goto 10
          end if
          rcocut = xdum
        end if
c
        call rel_cont_order (rcocut)
        goto 10
c
c ... LIst_selected
c
      else if (optpar(1)(1:2) .eq. 'LI') then
c
        if (nopt .lt. 2) optpar(2) = 'R'
        call upcase (optpar(2))
c
        if (optpar(2) .eq. '?') then
          call prompt (' which  = Residues | Atoms')
          goto 10
        end if
c
        call list_select (optpar(2))
        goto 10
c
c ... BFactor commands
c
      else if (optpar(1)(1:2) .eq. 'BF') then
c
        if (nselect .le. 0) then
          call errcon ('No atoms selected')
          goto 10
        end if
c
        if (nopt .lt. 2) then
          optpar(2) = 'STats'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'ST') then
c
c ... BFactor STats [how]
c     how    = Chain | Type
c
          if (nopt. lt. 3) then
            optpar(3) = 'Chain'
          else if (optpar(3) .eq. '?') then
            call prompt (' how    = Chain | Type')
            goto 10
          end if
c
          call bq_stats (optpar(3),'B')
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'NO') then
c
          do i=1,natoms
            laniso (i) = .false.
          end do
          naniso = 0
          call prompt (' All ANISOUs deleted')
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
          if (nopt .lt. 3) then
            write (optpar(3),'(f8.2)') blimlo
            call textin (' Lower B limit ?',optpar(3))
          end if
          call str2r (optpar(3),dummy,ierr)
          if (ierr .ne. 0) goto 10
          blimlo = dummy
c
          if (nopt .lt. 4) then
            write (optpar(4),'(f8.2)') blimhi
            call textin (' Upper B limit ?',optpar(4))
          end if
          call str2r (optpar(4),dummy,ierr)
          if (ierr .ne. 0) goto 10
          blimhi = dummy
c
          call bq_limit ('B')
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'PL') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'b.plt'
            call textin (' Plot file ?',optpar(3))
          end if
c
          if (optpar(3) .eq. '?') then
            call prompt (' file   = O2D plot file name')
            call prompt (' what   = Average | Radial')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = 'Average'
            call textin (' What ?',optpar(4))
          end if
c
          call bplot (iunit,optpar(3),optpar(4))
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'GR') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'Mc_sc'
            call textin (' How ?',optpar(3))
          end if
c
          if (optpar(3) .eq. '?') then
            call prompt (' how    = Mc_sc | Residue | Overall')
            goto 10
          end if
c
          call bgroup (optpar(3))
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'SM') then
c
          call bsmooth ()
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'BO') then
c
          call bbond ()
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'SA') then
c
          if (ltrueb) then
            do i=1,natoms
              bsave (i) = batom (i)
            end do
            call prompt (' Saved all atomic B-factors')
          else
            call errcon (' Cannot save Pseudo-B-factors !')
          end if
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          do i=1,natoms
            batom (i) = bsave (i)
          end do
          call prompt (' Restored all atomic B-factors')
          ltrueb = .true.
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'OD') then
c
          if (nopt .lt. 3) optpar(3) = 'bfactors.odb'
          if (nopt .lt. 4) optpar(4) = 'm1'
c
          call bodb (iunit,optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'SC') then
c
          if (nopt .lt. 3) optpar(3) = '1.0'
          if (nopt .lt. 4) optpar(4) = '100.0'
c
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          call str2r (optpar(4),ydum,ierr)
          if (ierr .ne. 0) goto 10
c
          call bscale (xdum,ydum)
c
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'PS') then
c
          call prompt (' Replace B-factors by Pseudo-B-factors ...')
c
          if (nopt .lt. 3) then
            optpar (3) = 'G'
            call textin (' What (X|Y|Z|Q|A|I|R|G|D|N|C|H) ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
          if (index('XYZQAIRGDNCH',optpar(3)(1:1)) .le. 0) then
            call errcon ('Invalid property type')
            goto 10
          end if
c
          if (optpar(3)(1:1) .eq. 'D') then
            if (nopt .lt. 4) optpar(4) = '0.0'
            if (nopt .lt. 5) optpar(5) = '0.0'
            if (nopt .lt. 6) optpar(6) = '0.0'
          else if (optpar(3)(1:1) .eq. 'N') then
            if (nopt .lt. 4) optpar(4) = '10.0'
          else if (optpar(3)(1:1) .eq. 'H') then
            if (nopt .lt. 4) optpar(4) = '7.35'
            if (.not. ltrueb) then
              call errcon ('Pseudo-B-values used !!!')
            end if
          else if (optpar(3)(1:1) .eq. 'C') then
            if (nopt .lt. 4) optpar(4) = '10.0'
            if (nopt .lt. 5) optpar(5) = '20.0'
          end if
c
          call bpseudo (optpar(3)(1:1),optpar(4),optpar(5),
     +                  optpar(6),ierr)
c
          if (ierr .eq. 0) ltrueb = .false.
c
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'PR') then
c
          if (nopt .lt. 3) optpar (3) = '1.0'
          if (nopt .lt. 4) optpar (4) = '0.0'
c
          call bq_prod_plus ('B',optpar(3),optpar(4))
          goto 10
c
        else
          call errcon ('Invalid BFactor command')
        end if
c
        goto 10
c
c ... OCcupancy commands
c
      else if (optpar(1)(1:2) .eq. 'OC') then
c
        if (nselect .le. 0) then
          call errcon ('No atoms selected')
          goto 10
        end if
c
        if (nopt .lt. 2) then
          optpar(2) = 'STats'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'ST') then
c
c ... OCcupancy STats [how]
c     how    = Chain | Type
c
          if (nopt. lt. 3) then
            optpar(3) = 'Chain'
          else if (optpar(3) .eq. '?') then
            call prompt (' how    = Chain | Type')
            goto 10
          end if
c
          call bq_stats (optpar(3),'Q')
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
          if (nopt .lt. 3) then
            write (optpar(3),'(f8.2)') qlimlo
            call textin (' Lower Q limit ?',optpar(3))
          end if
          call str2r (optpar(3),dummy,ierr)
          if (ierr .ne. 0) goto 10
          qlimlo = dummy
c
          if (nopt .lt. 4) then
            write (optpar(4),'(f8.2)') qlimhi
            call textin (' Upper Q limit ?',optpar(4))
          end if
          call str2r (optpar(4),dummy,ierr)
          if (ierr .ne. 0) goto 10
          qlimhi = dummy
c
          call bq_limit ('Q')
c
        else if (optpar(2)(1:2) .eq. 'PR') then
c
          if (nopt .lt. 3) optpar (3) = '1.0'
          if (nopt .lt. 4) optpar (4) = '0.0'
c
          call bq_prod_plus ('Q',optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'PL') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'q.plt'
            call textin (' Plot file ?',optpar(3))
          end if
c
          if (optpar(3) .eq. '?') then
            call prompt (' file   = O2D plot file name')
            call prompt (' what   = Average | Radial')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = 'Average'
            call textin (' What ?',optpar(4))
          end if
c
          call qplot (iunit,optpar(3),optpar(4))
          goto 10
c
        else
          call errcon ('Invalid OCcupancy command')
        end if
c
        goto 10
c
c ... SPlit file_prefix
c
      else if (optpar(1)(1:2) .eq. 'SP') then
c
        if (nopt .lt. 2) then
          optpar (2) = './m1'
          call textin (' File prefix ?',optpar(2))
        end if
c
        call split_pdb (iunit,junit,optpar(2))
c
        goto 10
c
c ... NUcleic commands
c
      else if (optpar(1)(1:2) .eq. 'NU') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'DU') then
c
          if (nopt .gt. 2) then
            if (optpar(3) .eq. '?') then
              call prompt (' file   = name of output file or blank')
              call prompt (
     +          ' what   = Plot | ACGTU_plot | Text_file')
              call prompt (
     +          ' chain  = * | _ | chain_name')
              goto 10
            end if
          else
            optpar (3) = ' '
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = 'R'
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = '*'
          end if
c
          if (optpar(5) .ne. '_') then
            call na_duarte (iunit,optpar(3),optpar(4),optpar(5))
          else
            line = achain(atmptr(1,1)) // '_'
            k = 1
            do i=1,nres
              j = index (line(1:k),achain(atmptr(1,i)))
              if (j .le. 0) then
                line = achain(atmptr(1,i)) // line(1:k)
                k = k + 1
              end if
            end do
            call upcase (line)
            write (*,'(/a,i2,a,a,a)') ' List of ',k,
     +        ' chain IDs found: |',line(1:k),'|'
            do i=k,1,-1
              optpar(6) = optpar(3)(1:length(optpar(3))) //
     +                    '_chain_' // line(i:i) // '.ps'
              call remspa (optpar(6))
              optpar (7) = line(i:i)
              write (*,*)
              call textut (' ==> Process chain   :',line(i:i))
              call textut (' ==> PostScript file :',optpar(6))
              call na_duarte (iunit,optpar(6),optpar(4),optpar(7))
            end do
          end if
c
          goto 10
c
        else
          call errcon ('Invalid NUcleic command')
        end if
c
c ... PRotein commands
c
      else if (optpar(1)(1:2) .eq. 'PR') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'MC') then
c
          if (nopt .gt. 2) then
            if (optpar(3) .eq. '?') then
              call prompt (' file   = name of output file or blank')
              call prompt (
     +          ' what   = Rama_plot | Labelled_rama | Text_file')
              call prompt (
     +          ' chain  = * | _ | chain_name')
              goto 10
            end if
          else
            optpar (3) = ' '
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = 'R'
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = '*'
          end if
c
          if (optpar(5) .ne. '_') then
            call mc_check (iunit,optpar(3),optpar(4),optpar(5))
          else
            line = achain(atmptr(1,1)) // '_'
            k = 1
            do i=1,nres
              j = index (line(1:k),achain(atmptr(1,i)))
              if (j .le. 0) then
                line = achain(atmptr(1,i)) // line(1:k)
                k = k + 1
              end if
            end do
            call upcase (line)
            write (*,'(/a,i2,a,a,a)') ' List of ',k,
     +        ' chain IDs found: |',line(1:k),'|'
            do i=k,1,-1
              optpar(6) = optpar(3)(1:length(optpar(3))) //
     +                    '_chain_' // line(i:i) // '.ps'
              call remspa (optpar(6))
              optpar (7) = line(i:i)
              write (*,*)
              call textut (' ==> Process chain   :',line(i:i))
              call textut (' ==> PostScript file :',optpar(6))
              call mc_check (iunit,optpar(6),optpar(4),optpar(7))
            end do
          end if
c
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'SC') then
c
          if (nopt .gt. 2) then
            if (optpar(3) .eq. '?') then
              call prompt (' file   = name of output file or blank')
              goto 10
            end if
          else
            optpar (3) = ' '
          end if
c
          call sc_check (iunit,optpar(3))
          close (iunit)
c
          goto 10
c
        else if (optpar(2)(1:2) .eq. 'CA') then
c
          if (nopt .gt. 2) then
            if (optpar(3) .eq. '?') then
              call prompt (' file   = name of output file or blank')
              call prompt (
     +          ' what   = Rama_plot | Labelled_rama | Text_file')
              call prompt (
     +          ' chain  = * | _ | chain_name')
              goto 10
            end if
          else
            optpar (3) = ' '
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = 'R'
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = '*'
          end if
c
          if (optpar(5) .ne. '_') then
            call ca_check (iunit,optpar(3),optpar(4),optpar(5))
          else
            line = achain(atmptr(1,1)) // '_'
            k = 1
            do i=1,nres
              j = index (line(1:k),achain(atmptr(1,i)))
              if (j .le. 0) then
                line = achain(atmptr(1,i)) // line(1:k)
                k = k + 1
              end if
            end do
            call upcase (line)
            write (*,'(/a,i2,a,a,a)') ' List of ',k,
     +        ' chain IDs found: |',line(1:k),'|'
            do i=k,1,-1
              optpar(6) = optpar(3)(1:length(optpar(3))) //
     +                    '_chain_' // line(i:i) // '.ps'
              call remspa (optpar(6))
              optpar (7) = line(i:i)
              write (*,*)
              call textut (' ==> Process chain   :',line(i:i))
              call textut (' ==> PostScript file :',optpar(6))
              call ca_check (iunit,optpar(6),optpar(4),optpar(7))
            end do
          end if
c
          goto 10
c
        else
          call errcon ('Invalid PRotein command')
        end if
c
        goto 10
c
c ... CHain commands
c
      else if (optpar(1)(1:2) .eq. 'CH') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'AU') then
c
          call chain_auto (0)
c
        else if (optpar(2)(1:2) .eq. 'AS') then
c
          call chain_auto (1)
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          if (nopt .lt. 3) then
            optpar (3) = ' '
            call textin (' Old chain name ?',optpar(3))
          end if
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' old    = old chain name (* means ANY)')
            call prompt (' new    = new chain name')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = ' '
            call textin (' New chain name ?',optpar(4))
          end if
c
          call chain_rename (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'SE') then
c
          if (nopt .lt. 3) then
            optpar (3) = ' '
            call textin (' Old segment id ?',optpar(3))
          end if
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' old    = old segment id (* means ANY)')
            call prompt (' new    = new segment id')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = ' '
            call textin (' New segment id ?',optpar(4))
          end if
c
          call chain_reseg (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'NA') then
c
          if (nopt .lt. 3) then
            optpar (3) = ' '
            call textin (' New chain name ?',optpar(3))
          end if
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' chain  = new chain name (= means no change)')
            call prompt (' segid  = new segment id (= means no change)')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = ' '
            call textin (' New segment id ?',optpar(4))
          end if
c
          call chain_select (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'FR') then
c
          if (nopt .lt. 3) optpar (3) = 'AUto'
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' how    = AUto | ASk')
            goto 10
          end if
c
          call chain_from_seg (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'TO') then
c
          if (nopt .lt. 3) optpar (3) = 'AUto'
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' how    = AUto | ASk')
            goto 10
          end if
c
          call chain_to_seg (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'OT') then
c
          call chain_ot2 ()
c
        else
          call errcon ('Invalid CHain command')
        end if
c
        goto 10
c
c ... PDb commands
c
      else if (optpar(1)(1:2) .eq. 'PD') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'HE') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'Deduce'
            call textin (' Option ?',optpar(3))
          end if
c
          if (optpar(3) .eq. '?') then
            call prompt (' option = Atom_all | Deduce')
            goto 10
          end if
c
          call pdb_hetero (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'NO') then
c
          do i=1,natoms
            inote (i) (3:4) = '  '
          end do
          call prompt (' O-style atomic numbers removed')
c
        else if (optpar(2)(1:2) .eq. 'SS') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'List'
            call textin (' What ?',optpar(3))
          end if
c
          if (optpar(3) .eq. '?') then
            call prompt (' what   = List | Delete | Generate')
            goto 10
          end if
c
          call pdb_ssbond (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          if (nopt .lt. 3) then
            optpar(3) = '???'
            call textin (' Text ?',optpar(3))
          end if
c
          optpar(4) = 'ADd'
          call pdb_remark (optpar(4),optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'NU') then
c
          if (nopt .lt. 3) then
            optpar (3) = '1000'
            call textin (' First new residue number ?',optpar(3))
          end if
c
          call pdb_number (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'NA') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'Atom'
            call textin (' Which ?',optpar(3))
          end if
          if (optpar(3) .eq. '?') then
            call prompt (' which  = Atom | Residue')
            call prompt (' old    = old name to replace')
            call prompt (' new    = new name')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar (4) = ' O1 '
            call textin (' Old ?',optpar(4))
          end if
c
          if (nopt .lt. 5) then
            optpar (5) = ' O  '
            call textin (' New ?',optpar(5))
          end if
c
          call pdb_name (optpar(3),optpar(4),optpar(5))
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
          optpar(3) = 'LIst'
          optpar(4) = '?'
          call pdb_remark (optpar(3),optpar(4))
c
c ... PDb FArout === PDb INdonesia
c
        else if (optpar(2)(1:2) .eq. 'FA') then
c
          call pdb_farout ()
c
c ... PDb SEqres
c
        else if (optpar(2)(1:2) .eq. 'SE') then
c
          call pdb_seqres ()
c
        else if (optpar(2)(1:2) .eq. 'IN') then
c
c ... PDb INdonesia === PDb FArout
c
          call pdb_farout ()
c
        else if (optpar(2)(1:2) .eq. 'CH') then
c
          call pdb_chem_charge (.true.,.false.)
c
        else if (optpar(2)(1:2) .eq. 'SA') then
c
          call pdb_sanity_check ()
c
        else if (optpar(2)(1:2) .eq. 'DE') then
c
          if (nopt .lt. 3) then
            optpar(3) = '-1'
            call textin (' Which ?',optpar(3))
          end if
c
          if (optpar(3) .eq. '?') then
            call prompt (' which  = number of the remark (* means ALL)')
            goto 10
          end if
c
          optpar(4) = 'DElete'
          call pdb_remark (optpar(4),optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'CR') then
c
          if (nopt .lt. 3) then
            write (optpar(3),'(f8.2)') cell(1)
            call textin (' A axis (A) ?',optpar(3))
          end if
          call str2r (optpar(3),dummy,ierr)
          if (ierr .ne. 0) goto 10
          cell (1) = max (1.0, dummy)
c
          if (nopt .lt. 4) then
            write (optpar(4),'(f8.2)') cell(2)
            call textin (' B axis (A) ?',optpar(4))
          end if
          call str2r (optpar(4),dummy,ierr)
          if (ierr .ne. 0) goto 10
          cell (2) = max (1.0, dummy)
c
          if (nopt .lt. 5) then
            write (optpar(5),'(f8.2)') cell(3)
            call textin (' C axis (A) ?',optpar(5))
          end if
          call str2r (optpar(5),dummy,ierr)
          if (ierr .ne. 0) goto 10
          cell (3) = max (1.0, dummy)
c
          if (nopt .lt. 6) then
            write (optpar(6),'(f8.2)') cell(4)
            call textin (' Alpha angle (deg) ?',optpar(6))
          end if
          call str2r (optpar(6),dummy,ierr)
          if (ierr .ne. 0) goto 10
          cell (4) = max (1.0, min (179.0, dummy))
c
          if (nopt .lt. 7) then
            write (optpar(7),'(f8.2)') cell(5)
            call textin (' Beta  angle (deg) ?',optpar(7))
          end if
          call str2r (optpar(7),dummy,ierr)
          if (ierr .ne. 0) goto 10
          cell (5) = max (1.0, min (179.0, dummy))
c
          if (nopt .lt. 8) then
            write (optpar(8),'(f8.2)') cell(6)
            call textin (' Gamma angle (deg) ?',optpar(8))
          end if
          call str2r (optpar(8),dummy,ierr)
          if (ierr .ne. 0) goto 10
          cell (6) = max (1.0, min (179.0, dummy))
c
          if (nopt .lt. 9) then
            write (optpar(9),'(i8)') zmol
            call textin (' Nr of molecules in cell ?',optpar(9))
          end if
          call str2i (optpar(9),idum,ierr)
          if (ierr .ne. 0) goto 10
          zmol = max (1, idum)
c
          if (nopt .lt. 10) then
            optpar(10) = spgrp
            call textin (' Spacegroup symbol ?',optpar(10))
          end if
          spgrp = optpar(10)
c
          call fvalut (' Unit-cell axes (A)      :',3,cell)
          call fvalut (' Unit-cell angles (deg)  :',3,cell(4))
          call celvol (cell,volume)
          call rvalut (' Unit-cell volume (A3)   :',1,volume)
          call ivalut (' Nr of molecules in cell :',1,zmol)
          call textut (' Spacegroup symbol       :',spgrp)
c
        else
          call errcon ('Invalid PDb command')
        end if
c
        goto 10
c
c ... SElection commands
c
      else if (optpar(1)(1:2) .eq. 'SE') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:1) .eq. '?') then
c
          call do_select ('? ',' ',' ')
c
        else if (optpar(2)(1:2) .eq. 'AL') then
c
          call do_select ('AL',' ',' ')
c
        else if (optpar(2)(1:2) .eq. 'NO') then
c
          call do_select ('NO',' ',' ')
c
        else if (optpar(2)(1:2) .eq. 'HY') then
c
          call do_select ('HY',' ',' ')
c
        else if (optpar(2)(1:2) .eq. 'EX') then
c
          call do_select ('EX',' ',' ')
c
        else if (optpar(2)(1:2) .eq. 'BY') then
c
          call select_by_residue ()
c
        else if (optpar(2)(1:2) .eq. 'NE') then
c
          call do_select ('NE',' ',' ')
c
        else if (optpar(2)(1:2) .eq. 'AN' .or.
     +           optpar(2)(1:2) .eq. 'BU' .or.
     +           optpar(2)(1:2) .eq. 'OR') then
c
          if (nopt .lt. 3) then
            optpar(3) = '?'
            call textin (' What ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3)(1:1) .eq. '?') then
        call prompt (' what   = CHain | SEgid | TYpe | CLass |')
        call prompt ('          REsidue | ATom | ALtloc | ANisou |')
        call prompt ('          STruc_sec')
        call prompt (' which  = which chain, segid, etc. to use')
        call prompt ('        - type can prot, nucl, wate, etc.')
        call prompt ('        - class can be main or side chain')
        call prompt ('        - residue can be Ala, HOH, BGL etc.')
        call prompt ('        - atom can be " CA ", " O1 " etc.')
        call prompt ('        - altloc can be " ", "A", "X" etc.')
        call prompt ('        - anisou can be T(rue) or F(alse)')
        call prompt ('        - struc_sec can be LOop, TUrn, ALpha,')
        call prompt ('          BEta, LEft-handed, NOn-protein')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = '?'
            call textin (' Which ?',optpar(4))
          end if
          call upcase (optpar(4))
c
          call do_select (optpar(2),optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'NU') then
c
          if (nopt .lt. 3) then
            optpar(3) = '?'
            call textin (' And/Or/Butnot ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' and_or_butnot = And | Or | Butnot')
            call prompt (
     +        ' what   = Residue_nr | B-factor | Occupancy |')
            call prompt ('          X-coord | Y-coord | Z-coord |')
            call prompt ('          Mass | Element | Coval_radius |')
            call prompt ('          Sec_struc')
            call prompt (' lo     = minimum value to select')
            call prompt (' hi     = maximum value to select')
            call prompt (' Sec_str: 0 = loop/turn, 1 = alpha,')
            call prompt ('          2 = beta, 3 = left-handed,')
            call prompt ('          -1 = non-protein residues')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = 'Residue_nr'
            call textin (' What ?',optpar(4))
          end if
          call upcase (optpar(4))
c
          if (nopt .lt. 5) then
            optpar(5) = '99998'
            call textin (' Lo ?',optpar(5))
          end if
          call upcase (optpar(5))
c
          if (nopt .lt. 6) then
            optpar(6) = '99999'
            call textin (' Hi ?',optpar(6))
          end if
          call upcase (optpar(6))
c
          call do_numeric_select (optpar(3),optpar(4),
     +                            optpar(5),optpar(6))
c
        else if (optpar(2)(1:2) .eq. 'PO') then
c
          if (nopt .lt. 3) then
            optpar(3) = '?'
            call textin (' And/Or/Butnot ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' and_or_butnot = And | Or | Butnot')
            call prompt (' x y z  = coordinates of point')
            call prompt (' lo     = minimum distance to select')
            call prompt (' hi     = maximum distance to select')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = '0.0'
            call textin (' X-coordinate of point ?',optpar(4))
          end if
          call upcase (optpar(4))
c
          if (nopt .lt. 5) then
            optpar(5) = '0.0'
            call textin (' Y-coordinate of point ?',optpar(5))
          end if
          call upcase (optpar(5))
c
          if (nopt .lt. 6) then
            optpar(6) = '0.0'
            call textin (' Z-coordinate of point ?',optpar(6))
          end if
          call upcase (optpar(6))
c
          if (nopt .lt. 7) then
            optpar(7) = '0.0'
            call textin (' Minimum distance to point ?',optpar(7))
          end if
          call upcase (optpar(7))
c
          if (nopt .lt. 8) then
            optpar(8) = '2.4'
            call textin (' Maximum distance to point ?',optpar(8))
          end if
          call upcase (optpar(8))
c
          call do_point_select (optpar(3),optpar(4),optpar(5),
     +         optpar(6),optpar(7),optpar(8))
c
        else if (optpar(2)(1:2) .eq. 'DI') then
c
          if (nopt .lt. 3) then
            optpar(3) = '0.0'
            call textin (' Lo ?',optpar(3))
          end if
          call upcase (optpar(3))
c
          if (nopt .lt. 4) then
            optpar(4) = '2.4'
            call textin (' Hi ?',optpar(4))
          end if
          call upcase (optpar(4))
c
          call select_near (optpar(3),optpar(4))
c
        else
          call errcon ('Invalid SElect command')
        end if
c
        goto 10
c
c ... ONo commands
c
      else if (optpar(1)(1:2) .eq. 'ON') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'RS' .or. optpar(2)(1:2) .eq. 'FI' .or.
     +      optpar(2)(1:2) .eq. 'TO' .or. optpar(2)(1:2) .eq. 'CO') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'XYZ'
            call textin (' Residue type ?',optpar(3))
          end if
          if (optpar(3) .eq. '?') then
            call prompt (
     +        ' res_typ = type of residue (3 characters; e.g. GLC)')
            call prompt (
     +        ' file    = name of output file (default: automatic)')
            call prompt (
     +        ' res_nr  = residue number (default: first)')
            goto 10
          end if
c
          if (nopt .lt. 4) optpar(4) = ' '
          if (nopt .lt. 5) optpar(5) = ' '
c
          call ono_dicts (optpar(2),optpar(3),optpar(4),optpar(5),
     +                    iunit)
c
        else if (optpar(2)(1:2) .eq. 'OO') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'M1'
            call textin (' Molecule name in O ?',optpar(3))
          end if
c
          if (nopt .lt. 4) optpar(4) = 'pre_oops.omac'
c
          call ono_oops (optpar(3),optpar(4),iunit)
c
        else if (optpar(2)(1:2) .eq. 'DI') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'M1'
            call textin (' Molecule name in O ?',optpar(3))
          end if
c
          if (nopt .lt. 4) optpar(4) = 'disulfide.odl'
c
c ... how ? Sticks or Lines ?
c
          if (nopt .lt. 5) optpar(5) = 'sticks'
c
          if (nopt .lt. 6) optpar(6) = '_ssbonds'
c
          call ono_disulfide (optpar(3),optpar(4),optpar(5),
     +                        optpar(6),iunit)
c
        else if (optpar(2)(1:2) .eq. 'WA') then
c
          if (nopt .lt. 3) optpar(3) = 'water_fit.omac'
c
          call ono_water_fit (optpar(3),iunit)
c
        else if (optpar(2)(1:2) .eq. 'LS') then
c
          if (nopt .lt. 3) optpar(3) = 'plane.odl'
          if (nopt .lt. 4) optpar(4) = '_plane'
          if (nopt .lt. 5) optpar(5) = 'cyan'
          if (nopt .lt. 6) write (optpar(6),*) border
          call str2r (optpar(6),dummy,ierr)
          if (ierr .eq. 0) border = dummy
          if (nopt .lt. 7) optpar(7) = '-1'
          call str2i (optpar(7),idum,ierr)
          if (ierr .ne. 0) idum = -1
          if (nopt .lt. 8) optpar(8) = ' '
c
          call xyz_ls_plane (iunit,optpar(3),optpar(4),optpar(5),
     +                       border,idum,optpar(8),'O')
c
        else if (optpar(2)(1:2) .eq. 'IN') then
c
          if (nopt .lt. 3) optpar(3) = 'inertia.odl'
          if (nopt .lt. 4) optpar(4) = '_inert'
          if (nopt .lt. 5) optpar(5) = 'cyan'
          if (nopt .lt. 6) write (optpar(6),*) veclen
          call str2r (optpar(6),dummy,ierr)
          if (ierr .eq. 0) veclen = dummy
c
          call xyz_ls_plane (iunit,optpar(3),optpar(4),optpar(5),
     +                       veclen,-1,' ','I')
c
        else if (optpar(2)(1:2) .eq. 'RI') then
c
          if (nopt .lt. 3) optpar(3) = 'rings.odl'
          if (nopt .lt. 4) optpar(4) = '_rings'
          if (nopt .lt. 5) optpar(5) = 'cyan'
c
          call textut (' ODL file name :',optpar(3))
          call textut (' Object name   :',optpar(4))
          call textut (' Object colour :',optpar(5))
c
          call ono_rings (iunit,optpar(3),optpar(4),optpar(5))
c
        else if (optpar(2)(1:2) .eq. 'CE') then
c
          if (nopt .lt. 3) optpar(3) = 'cell.odl'
          if (nopt .lt. 4) optpar(4) = 'line'
          if (nopt .lt. 5) optpar(5) = 'red'
c
          call textut (' ODL file :',optpar(3))
          call textut (' Mode     :',optpar(4))
          call textut (' Colour   :',optpar(5))
c
          call upcase (optpar(4))
c
          call ono_cell (iunit,optpar(3),optpar(4),optpar(5))
c
        else if (optpar(2)(1:2) .eq. 'MO') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'molray_trace.pdb'
            call textin (' Trace PDB file name ?',optpar(3))
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = 'peptide'
            call textin (' Mode (PEPT) ?',optpar(4))
          end if
c
          if (nopt .lt. 5) then
            write (optpar(5),'(f8.2)') tramul
            call textin (' Multiplier ?',optpar(5))
          end if
          call str2r (optpar(5),tramul,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            optpar(6) = 'nearest'
            call textin (' Near or Far nbrs (NEA/FAR) ?',optpar(6))
          end if
c
          if (nopt .lt. 7) then
            optpar(7) = 'forward'
            call textin (' Forward or reverse (FOR/REV) ?',optpar(7))
          end if
c
          if (nopt .lt. 8) optpar(8) = 'HOH'
          if (nopt .lt. 9) optpar(9) = ' O  '
c
          close (iunit)
c
          call ono_molray_trace (iunit,optpar(3),optpar(4),tramul,
     +      optpar(6),optpar(7),optpar(8),optpar(9),rbuf)
c
          close (iunit)
c
        else if (optpar(2)(1:2) .eq. 'XP') then
c
          call ono_fix_hydro ()
c
        else if (optpar(2)(1:2) .eq. 'AN') then
c
          if (nopt .lt. 3) optpar(3) = '-1'
          call str2i (optpar(3),idum,ierr)
          if (ierr .ne. 0) idum = -1
c
          if (nopt .lt. 4) optpar(4) = '-1'
          call str2i (optpar(4),jdum,ierr)
          if (ierr .ne. 0) jdum = -1
c
          call ono_ls_angles (idum,jdum)
c
        else
          call errcon ('Invalid ONo command')
        end if
c
        goto 10
c
c ... SQuence commands
c
      else if (optpar(1)(1:2) .eq. 'SQ') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'LI') then
c
          if (nopt .lt. 3) then
            optpar (3) = '3'
            call textin (' How ?',optpar(3))
          end if
          if (optpar(3) .eq. '?') then
            call prompt (
     +        ' how    = 3-letter code | 1-letter code | Full')
            goto 10
          end if
c
          call sequence_list (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'MO') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'N?S'
            call textin (' Motif ?',optpar(3))
          end if
          if (optpar(3) .eq. '?') then
            call prompt (
     +        ' motif  = sequence motif in 1-letter code')
            goto 10
          end if
c
          call sequence_motif (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'GL') then
c
          call sequence_glyco ()
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          call sequence_count ()
c
        else if (optpar(2)(1:2) .eq. 'EX') then
c
          call sequence_ex_280 ()
c
        else if (optpar(2)(1:2) .eq. 'PI') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'sequence.pir'
            call textin (' File ?',optpar(3))
          end if
          if (optpar(3) .eq. '?') then
            call prompt (' file    = PIR output file name')
            call prompt (' seq_nam = PIR sequence name')
            call prompt (' title   = any text')
            goto 10
          end if
c
          if (nopt .lt. 4) optpar(4) = 'PIRSEQ'
          if (nopt .lt. 5) optpar(5) = 'No title'
c
          call sequence_pir (iunit,optpar(3),optpar(4),optpar(5))
c
        else
          call errcon ('Invalid SQuence command')
        end if
c
        goto 10
c
c ... DIstance commands
c
      else if (optpar(1)(1:2) .eq. 'DI') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'PL') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'dist.pl2'
            call textin (' File ?',optpar(3))
          end if
c
          call dist_plot (optpar(3),iunit)
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'conv.pl2'
            call textin (' File ?',optpar(3))
          end if
c
          if (nopt .ge. 4) then
            call str2i (optpar(4),idum,ierr)
            if (ierr .ne. 0) goto 10
            ncolen = idum
          end if
c
          if (nopt .ge. 5) then
            acotyp = optpar(5)
          end if
c
          call conv_plot (optpar(3),ncolen,acotyp,iunit)
c
        else if (optpar(2)(1:2) .eq. 'DI') then
c
          if (nopt .lt. 3) optpar(3) = '2.0'
c
          call dist_dist (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'SH') then
c
          if (nopt .lt. 3) optpar(3) = '2.4'
c
          call dist_short (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
          if (nopt .lt. 3) optpar(3) = '0.0'
          if (nopt .lt. 4) optpar(4) = '1.0'
c
          call dist_list (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'CH') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'A'
            call textin (' Chain 1 ?',optpar(3))
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = 'B'
            call textin (' Chain 2 ?',optpar(4))
          end if
c
          if (nopt .lt. 5) optpar(5) = '3.5'
c
          call dist_chain (optpar(3),optpar(4),optpar(5))
c
        else if (optpar(2)(1:2) .eq. 'SE') then
c
          if (nopt .lt. 3) optpar(3) = '2.4'
c
          if (nopt .lt. 4) optpar(4) = 'ALL'
c
          call dist_select (optpar(3),optpar(4))
c
        else
          call errcon ('Invalid DIstance command')
        end if
c
        goto 10
c
c ... XYz commands
c
      else if (optpar(1)(1:2) .eq. 'XY') then
c
        if (nopt .lt. 2) then
          optpar(2) = '?'
          call textin (' Option ?',optpar(2))
        end if
        call remspa (optpar(2))
        call upcase (optpar(2))
c
        if (optpar(2)(1:2) .eq. 'FR') then
c
          call frac_cart (0)
c
        else if (optpar(2)(1:2) .eq. 'OR') then
c
          call frac_cart (1)
c
        else if (optpar(2)(1:2) .eq. 'RO') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'CEuler'
            call textin (' How ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3) .eq. '?') then
            call prompt (
     +        ' how    = CEuler | CPolar | MEuler | MPolar |')
            call prompt (
     +        '          XPolar | XLattmann')
            call prompt (
     +        '          (C=CCP4, M=Merlot, X=X-PLOR)')
            call prompt (
     +        ' ang1-3 = three rotation angles (default 0,0,0)')
            goto 10
          end if
c
          do i=4,6
            if (nopt .lt. i) then
              optpar(i) = '0.0'
              call textin (' Angle ?',optpar(i))
            end if
          end do
c
          call xyz_rotate (optpar(3),optpar(4),optpar(5),optpar(6))
c
        else if (optpar(2)(1:2) .eq. 'AX') then
c
          if (nopt .lt. 3) then
            optpar(3) = 'X'
            call textin (' Axis [X|Y|Z] ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
c
          if (optpar(3) .eq. '?') then
            call prompt (
     +        ' axis  = X | Y | Z')
            call prompt (
     +        ' angle = rotation angle around the selected axis')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = '0.0'
            call textin (' Angle ?',optpar(4))
          end if
c
          call xyz_axis_rotate (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'TR') then
c
          do i=3,5
            if (nopt .lt. i) then
              optpar(i) = '0.0'
              call textin (' Translation component ?',optpar(i))
            end if
          end do
c
          call xyz_trans (optpar(3),optpar(4),optpar(5))
c
        else if (optpar(2)(1:2) .eq. 'RA') then
c
          call xyz_ranrot ()
c
        else if (optpar(2)(1:2) .eq. 'PE') then
c
          do i=1,5
            if (nopt .lt. (i+2)) optpar (i+2) = '0.0'
          end do
c
          call xyz_perturb (optpar(3),optpar(4),optpar(5),
     +      optpar(6),optpar(7))
c
        else if (optpar(2)(1:2) .eq. 'MI') then
c
          if (nopt .lt. 3) optpar(3) = 'X'
          if (nopt .lt. 4) optpar(4) = '0.0'
c
          call xyz_mirror (optpar(3),optpar(4))
c
        else if (optpar(2)(1:2) .eq. 'IN') then
c
          do i=1,3
            if (nopt .lt. (i+2)) optpar (i+2) = '0.0'
          end do
c
          call xyz_invert (optpar(3),optpar(4),optpar(5))
c
        else if (optpar(2)(1:2) .eq. 'MA' .or.
     +           optpar(2)(1:2) .eq. 'DI') then
c
          do i=1,12,4
            if (nopt .lt. (2+i)) optpar(2+i) = '1.0'
            do j=i+1,i+3
              if (nopt .lt. (2+j)) optpar(2+j) = '0.0'
            end do
          end do
c
          ldist = (optpar(2)(1:2) .eq. 'DI')
          call xyz_user_rt (optpar(3),ldist)
c
        else if (optpar(2)(1:2) .eq. 'RT') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'rt.o'
            call textin (' O LSQ operator file ?',optpar(3))
          end if
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' file   = LSQ operator in O format')
            goto 10
          end if
c
          call xyz_o_rt (iunit,optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'CE') then
c
          if (nopt .lt. 3) optpar(3) = 'N'
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (' use_masses   = Yes or No (default)')
            goto 10
          end if
          call upcase (optpar(3))
          call remspa (optpar(3))
          if (optpar(3)(1:1) .eq. 'N') then
            call xyz_origin ()
          else
            call xyz_mass_origin ()
          end if
c
        else if (optpar(2)(1:2) .eq. 'AL') then
c
          call xyz_inertia ()
c
        else
          call errcon ('Invalid XYz command')
        end if
c
        goto 10
c
c ... VRML
c
      else if (optpar(1)(1:2) .eq. 'VR') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'SEtup'
          call textin (' Option ?',optpar(2))
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
          call xopxua (ivrml,vrfile,xinter(),ierr)
          if (ierr .ne. 0) then
            call errcon (' Could not open VRML file')
            goto 10
          end if
c
          call xvrml_open (ivrml,rgbbg1,rgbbg2,rgbbg3)
          lvrml = .true.
          call xvrml_colour (rgbfg1,rgbfg2,rgbfg3)
c
        else if (optpar(2)(1:2) .eq. 'CL') then
c
c ... VRML CLOSE_FILE
c
          if (lvrml) then
            call xvrml_close ()
            lvrml = .false.
ccc            call prompt (' VRML file closed')
          else
            call prompt (' No open VRML file at present')
          end if
c
        else if (optpar(2)(1:2) .eq. 'TR') then
c
c ... VRML TRACE
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .ge. 3) then
            call xvrml_rgb_name (optpar(3),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          call prompt (' VRML Trace ...')
c
          call vtrace (vrdist,cavrml)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'FA') then
c
c ... VRML FAT_TRACE
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = '1.0'
            call textin (' Fatness (A) ?',optpar(3))
          end if
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          lmono = .false.
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML Fat trace ...')
c
          call vturd (lmono,vrdist,cavrml,xdum)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'CP') then
c
c ... VRML CPK
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = '1.4'
            call textin (' Sphere radius (A) ?',optpar(3))
          end if
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          lmono = .false.
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML CPK ...')
c
          call vcpk (lmono,xdum)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'ST') then
c
c ... VRML STICK
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          lmono = .false.
          if (nopt .ge. 3) then
            call xvrml_rgb_name (optpar(3),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML Stick ...')
c
          call vstick (.false.,.false.,.true.,lmono,0.0)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'BA') then
c
c ... VRML BALL_CYLINDER
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = '0.2'
            call textin (' Cylinder radius (A) ?',optpar(3))
          end if
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          lmono = .false.
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML Ball-and-cylinder ...')
c
          call vstick (.true.,.true.,.false.,lmono,xdum)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'CY') then
c
c ... VRML CYLINDER
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = '0.2'
            call textin (' Cylinder radius (A) ?',optpar(3))
          end if
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          lmono = .false.
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML Cylinder ...')
c
          call vstick (.false.,.true.,.false.,lmono,xdum)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
c ... VRML LIQUORICE
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = '0.5'
            call textin (' Radius (A) ?',optpar(3))
          end if
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          lmono = .false.
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML Liquorice ...')
c
          call vstick (.false.,.true.,.false.,lmono,xdum)
          call vcpk (lmono,xdum)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'SP') then
c
c ... VRML SPHERE
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          lmono = .false.
          if (nopt .ge. 3) then
            call xvrml_rgb_name (optpar(3),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
            lmono = .true.
          end if
c
          call prompt (' VRML Sphere ...')
c
          call vstick (.true.,.false.,.false.,lmono,xdum)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'NA') then
c
c ... VRML NAMED_COLOURS
c
          call xvrml_col_list ()
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
c ... VRML RESET
c
          call prompt (' Resetting atom colours and radii ...')
c
          call pdb_chem_charge (.false.,.true.)
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
c ... VRML COLOUR_SELECTION
c
          if (nopt .lt. 3) then
            optpar (3) = 'white'
            call textin (' Colour for selection ?',optpar(3))
          end if
          call xvrml_rgb_name (optpar(3),rr,gg,bb)
          call xvrml_encode_rgb (rr,gg,bb,idum)
c
          j = 0
          do i=1,natoms
            if (select(i)) then
              colour (i) = idum
              j = j + 1
            end if
          end do
c
          call jvalut (' Nr of atoms coloured :',1,j)
c
        else if (optpar(2)(1:2) .eq. 'CR') then
c
c ... VRML CRAMP_SELECTION
c
          if (nopt .lt. 3) then
            optpar (3) = '?'
            call textin (' Colour ramp criterion ?',optpar(3))
          end if
          call remspa (optpar(3))
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. '?') then
            call prompt (
     +        ' how    = REsi_nr | BFactor | OCcupancy | MAss |')
            call prompt (
     +        '          ELement | RAdius | ATom_nr | X | Y | Z')
            goto 10
          end if
c
          call vcramp (optpar(3))
c
        else if (optpar(2)(1:2) .eq. 'RA') then
c
c ... VRML RADII
c
          if (nopt .lt. 3) then
            optpar (3) = '1.0'
            call textin (' Multiply radii by ?',optpar(3))
          end if
          call str2r (optpar(3),xdum,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) then
            optpar (4) = '0.0'
            call textin (' Increase radii by ?',optpar(4))
          end if
          call str2r (optpar(4),ydum,ierr)
          if (ierr .ne. 0) goto 10
c
          j = 0
          do i=1,natoms
            if (select(i)) then
              cvbrad (i) = xdum*cvbrad(i) + ydum
              j = j + 1
            end if
          end do
c
          call jvalut (' Nr of radii changed :',1,j)
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
            write (optpar(4),'(f6.2)') vrdist
            call textin (' Max central atom distance ?',optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .eq. 0) vrdist = xdum
          call fvalut (' Max central atom distance :',1,vrdist)
c
          if (nopt .lt. 5) then
            write (optpar(5),'(3f6.3)') rgbbg1,rgbbg2,rgbbg3
            call pretty (optpar(5))
            call textin (' Background colour ?',optpar(5))
          end if
          call xvrml_rgb_name (optpar(5),rgbbg1,rgbbg2,rgbbg3)
          write (optpar(5),'(3f6.3)') rgbbg1,rgbbg2,rgbbg3
          call pretty (optpar(5))
          call textut (' Background colour :',optpar(5))
c
          if (nopt .lt. 6) then
            write (optpar(6),'(3f6.3)') rgbfg1,rgbfg2,rgbfg3
            call pretty (optpar(6))
            call textin (' Default colour ?',optpar(6))
          end if
          call xvrml_rgb_name (optpar(6),rgbfg1,rgbfg2,rgbfg3)
          write (optpar(6),'(3f6.3)') rgbfg1,rgbfg2,rgbfg3
          call pretty (optpar(6))
          call textut (' Default colour :',optpar(6))
c
        else if (optpar(2)(1:2) .eq. 'CE') then
c
c ... VRML CELL
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .ge. 3) then
            call xvrml_rgb_name (optpar(3),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          call fvalut (' VRML Cell :',6,cell)
c
          call xvrml_cell (cell,0,0,0,.true.)
c
          call flusho (ivrml)
c
        end if
c
c ... INVALID OPTION
c
      else
c
        call errcon ('Invalid command')
        call textut (' ==>',line)
c
      end if
c
      goto 10
c
c --- END OF MAIN EVENT LOOP
c


c
c ... end of program
c
 9000 continue
c
      call gkquit ()
c
      end
c
c
c
      subroutine read_pdb (iunit,mode,file,form,hydro)
c
      include 'moleman2.incl'
c
      integer i,j,iunit,ierr,mode,nline,nhydro,kk
c
      logical lalwyn,lhydro,lhkeep
c
      character*(*) file,form,hydro
      character key*6,line*120,dummy*6,oldatm*21
c
code ...
c
      call upcase (form)
      lalwyn = (form(1:1) .eq. 'A')
c
      call upcase (hydro)
      lhkeep = (hydro(1:1) .eq. 'Y')
c
      call textut (' Reading from file :',file)
c
      if (lalwyn) then
        call prompt (' in Alwyn''s PDB format')
      else
        call prompt (' in normal PDB format')
      end if
c
      if (lhkeep) then
        call prompt (' including hydrogen atoms')
      else
        call prompt (' ignoring hydrogen atoms')
      end if
c
      call xopxoa (iunit,file,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file')
        return
      end if
c
c ... reset counts for REad, don't for APpend
c
      if (mode .eq. 0) then
        natoms = 0
        nres   = 0
        nother = 0
        nrem   = 0
        nter   = 0
        naniso = 0
        do i=1,3
          cell (i) = 1.0
          cell (i+3) = 90.0
        end do
        spgrp = 'P 1'
        zmol = 1
      end if
c
      nline  = 0
      nhydro = 0
      oldatm = '???'
c
 6000 format (a6,a)
c
   10 continue
      read (iunit,6000,err=8000,end=9000) key,line
      call upcase (key)
      nline = nline + 1
c
      if (key .eq. 'ATOM  ' .or. key .eq. 'HETATM') goto 20
c
      if (key .eq. 'ANISOU') goto 25
c
      call pdbinfo (key,line)
c
      if (key .eq. 'END   ') then
        call prompt ('0>>>>> END card encountered <<<<<')
        goto 9000
      end if
c
      if (key .eq. 'TER   ') then
        if (nter .ge. maxter) then
          call errcon ('Too many TER records - skipping')
          goto 10
        end if
        nter = nter + 1
        terrec (nter) = line
        terptr (nter) = natoms
        goto 10
      end if
c
      if (key(1:5) .eq. 'ORIGX') then
        call prompt (' ORIGXn card ignored')
        goto 10
      end if
c
      if (key(1:5) .eq. 'SCALE') then
        call prompt (' SCALEn card ignored')
        goto 10
      end if
c
      if (key(1:5) .eq. 'CRYST') then
        read (line,*,err=10) (cell(i),i=1,6)
        spgrp = line(50:60)
        zmol = 1
        read (line(61:64),*,err=10,end=10) zmol
        goto 10
      end if
c
      if (key .eq. 'REMARK') then
        if (nrem .lt. maxcom) then
          nrem = nrem + 1
          remark (nrem) = key//line
          if (nrem .eq. maxcom) call prompt (
     +      ' WARNING - max nr of REMARK records read')
        end if
        goto 10
      end if
c
      if (nother .lt. maxcom) then
        nother = nother + 1
        other (nother) = key//line
        if (nother .eq. maxcom) call prompt (
     +    ' WARNING - max nr of OTHER records read')
      end if
c
      goto 10
c
c ... new atom/hetatm
c
   20 continue
c
      call upcase (line)
c
c ... strip hydrogens ?
c
      if (.not. lhkeep) then
        if (lhydro(line(7:10))) then
          nhydro = nhydro + 1
          goto 10
        end if
      else
        if (lhydro(line(7:10))) then
          nhydro = nhydro + 1
        end if
      end if
c
      natoms = natoms + 1
      if (natoms .gt. maxatm) then
        call errcon ('Too many atoms; rest skipped')
        call jvalut (' Max nr of atoms :',1,maxatm)
        call textut (' Offending line  :',line)
        natoms = natoms - 1
        goto 9000
      end if
c
 4713 format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a40)
 4715 format (i5,1x,a4,a1,a3,1x,   a6   ,3x,3f8.3,2f6.2,a40)
c
      kk = natoms
      if (lalwyn) then
        read (line,4715,err=8000) atomnr(kk),atmnam(kk),altloc(kk),
     +    resnam(kk),dummy,
     +    (xyz(j,kk),j=1,3),qatom(kk),batom(kk),inote(kk)
        call detaja (dummy,achain(kk),iresid(kk),insert(kk))
      else
        read (line,4713,err=8000) atomnr(kk),atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    (xyz(j,kk),j=1,3),qatom(kk),batom(kk),inote(kk)
      end if
c
      oldatm = line(1:21)
      if (atomnr(kk) .le. 0) then
        atomnr(kk) = natoms
      end if
      resptr (kk) = -1
      lmain (kk) = .false.
      lhet (kk) = (key .eq. 'HETATM')
      weight (kk) = 1.0
      laniso (kk) = .false.
      bsave (kk) = batom (kk)
c
      goto 10
c
cATOM     74  NZ XLYS     8      17.216  31.254 -22.466  0.55 28.30           N
cANISOU   74  NZ XLYS     8     2886   2496   5105    713    930     17       N
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8
c
   25 continue
      if (line(1:21) .eq. oldatm) then
        read (line(23:64),'(6i7)',err=8000) (anisou(j,natoms),j=1,6)
        laniso (natoms) = .true.
        naniso = naniso + 1
      end if
c
      goto 10
c
cATOM    163  CB  LYS    20      62.798   1.269   0.832  1.00 20.00   6
c1234567890123456789012345678901234567890123456789012345678901234567890
c      1234567890123456789012345678901234567890123456789012345678901234567890
c
c ... read error
c
 8000 continue
      call errcon ('While reading PDB record')
      call jvalut (' Near line :',1,nline)
      line = key//line
      call textut (' Line :',line)
      return
c
c ... end of file
c
 9000 continue
      write (*,*)
      call jvalut (' Nr of lines read :',1,nline)
      if (.not. lhkeep) then
        call jvalut (' Nr of hydrogens skipped :',1,nhydro)
      else
        call jvalut (' Nr of hydrogens kept :',1,nhydro)
      end if
c
      terptr (nter+1) = natoms + 1
c
      call pdb_chem_charge (.false.,.true.)
c
      call book_keep (.false.)
c
      return
      end
c
c
c
      subroutine book_keep (lprint)
c
      include 'moleman2.incl'
c
      real cotol,dx,dy,dz
c
      integer i,j,iold,nerr,ires,ityp,ii,jj,nsp,nhb,nlb,nhq,nlq,nqo
c
      logical check1(maxapr),check2(maxapr)
      logical lhydro,lprint
c
      character pres*3,pchn*1,pseg*4,pins*1
c
code ...
c
      if (natoms .lt. 1) then
        nres = 0
        return
      end if
c
      iold = -999
      pres = '???'
      pchn = '?'
      pins = '?'
      pseg = '????'
      nres = 0
c
      do i=1,nrlcat+1
        icnt(i) = 0
      end do
c
c ... delineate residues
c
      do i=1,natoms
c
        if (iresid(i) .ne. iold .or.
     +      resnam(i) .ne. pres .or.
     +      achain(i) .ne. pchn .or.
     +      insert(i) .ne. pins .or.
     +      inote(i)(7:10) .ne. pseg) then
c
          nres = nres + 1
          if (nres .gt. maxres) then
            call errcon ('Too many residues')
            call jvalut (' Max allowed :',1,maxres)
            nres = maxres
            atmptr (2,nres) = i - 1
            goto 2234
          end if
c
          resres (nres) = resnam (i)
          call tellib (resres(nres),typptr(nres),.false.)
          if (typptr(nres) .gt. 0) then
            restyp (nres) = lrtype (typptr(nres))
            onelc (nres) = libolc(typptr(nres))
          else
            restyp (nres) = nrlcat + 1
            onelc (nres) = '?'
          end if
          icnt (restyp(nres)) = icnt (restyp(nres)) + 1
c
          atmptr (1,nres) = i
c
          if (nres .gt. 1) then
            atmptr (2,nres-1) = i - 1
          end if
c
          iold = iresid (i)
          pres = resnam (i)
          pchn = achain (i)
          pins = insert (i)
          pseg = inote(i)(7:10)
c
        end if
c
        resptr (i) = nres
c
      end do
c
      if (nres .gt. 0) then
        atmptr (2,nres) = natoms
      end if
c
 2234 continue
      write (*,*)
      call jvalut (' Total nr of residues      :',1,nres)
c
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
      if (nres .le. 1) return
c
c ... find missing and extra atoms
c
      call prompt (
     +  '0Checking for missing/extra atoms ...')
c
      do ires=1,nres
        ityp = typptr(ires)
        if (ityp .gt. 0) then
          do i=1,maxapr
            check1(i) = .false.
            check2(i) = .false.
          end do
          do i=atmptr(1,ires),atmptr(2,ires)
            ii = i - atmptr(1,ires) + 1
            do j=nmrptr(1,ityp),nmrptr(2,ityp)
              jj = j - nmrptr(1,ityp) + 1
              if (atmnam(i) .eq. lratom(j) .or. (
     + ((atmnam(i)(1:1).eq.lratom(j)(1:1)).or.lratom(j)(1:1).eq.'?')
     +  .and.
     + ((atmnam(i)(2:2).eq.lratom(j)(2:2)).or.lratom(j)(2:2).eq.'?')
     +  .and.
     + ((atmnam(i)(3:3).eq.lratom(j)(3:3)).or.lratom(j)(3:3).eq.'?')
     +  .and.
     + ((atmnam(i)(4:4).eq.lratom(j)(4:4)).or.lratom(j)(4:4).eq.'?')
     +  ) ) then
                if (check1(ii)) then
                  if (lprint) then
                    call errcon ('Duplicate atom in library')
                    call textut (' Library residue type :',
     +                lrname(ityp))
                    call textut (' Atom name :',lratom(j))
                  end if
                else if (check2(jj)) then
                  if (lprint) then
                    call errcon ('Duplicate atom in structure')
                    call print_atom (i)
                  end if
                else
                  check1(ii) = .true.
                  check2(jj) = .true.
                  lmain (i) = ismain (j)
                end if
              end if
            end do
          end do
c
          if (lprint) then
            do i=atmptr(1,ires),atmptr(2,ires)
              ii = i - atmptr(1,ires) + 1
              if (.not. check1(ii)) then
                if (.not. lhydro(atmnam(i))) then
                  call errcon ('Unknown atom in structure')
                  call print_atom(i)
                end if
              end if
            end do
            do j=nmrptr(1,ityp),nmrptr(2,ityp)
              jj = j - nmrptr(1,ityp) + 1
              if (.not. check2(jj)) then
                if (.not. lhydro(atmnam(i))) then
                  call errcon ('Missing atom in structure')
                  call print_res (ires,0)
                  call textut (' Atom name :',lratom(j))
                end if
              end if
            end do
          end if
        else
          call errcon ('Residue not in library')
          call print_res (ires,0)
        end if
      end do
c
 5000 continue
c
c ... find CAs of amino acids
c
      do i=1,nres
        captr (i) = -1
        if (restyp(i) .eq. iprot) then
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' CA ') then
              captr (i) = j
              goto 1234
            end if
          end do
        end if
 1234   continue
      end do
c
c ... check for non-unique residue names
c
      nerr = 0
c
 6000 format (' Non-unique residue names : # ',i4,' = ',a3,1x,a1,
     +  i4,a1,1x,a4,' <-> # ',i4,' = ',a3,1x,a1,i4,a1,1x,a4)
c
      if (.not. lprint) goto 1000
      call prompt ('0Checking uniqueness of residue names ...')
      do i=1,nres-1
        do j=i+1,nres
          if ( iresid(atmptr(1,i)) .eq.
     +         iresid(atmptr(1,j)) ) then
            if ( achain(atmptr(1,i)) .eq.
     +           achain(atmptr(1,j)) ) then
              if ( insert(atmptr(1,i)) .eq.
     +             insert(atmptr(1,j)) ) then
              if ( inote(atmptr(1,i))(7:10) .eq.
     +             inote(atmptr(1,j))(7:10) ) then
          write (*,6000) i,resnam(atmptr(1,i)),achain(atmptr(1,i)),
     +      iresid(atmptr(1,i)),insert(atmptr(1,i)),
     +      inote(atmptr(1,i))(7:10),j,resnam(atmptr(1,j)),
     +      achain(atmptr(1,j)),iresid(atmptr(1,j)),
     +      insert(atmptr(1,j)),inote(atmptr(1,j))(7:10)
                nerr = nerr + 1
                if (nerr .gt. 10) then
                  call errcon ('More than 10 non-unique names')
                  goto 1000
                end if
              end if
              end if
            end if
          end if
        end do
      end do
c
      if (nerr .eq. 0) then
        call prompt (' All residue names are unique')
      end if
c
 1000 continue
c
c ... check "special" positions (e.g., 0,0,0) and unusual Bs and Qs
c
      cotol = 0.05
      nsp = 0
      nhb = 0
      nlb = 0
      nhq = 0
      nlq = 0
      nqo = 0
      write (*,*)
      call prompt (
     +  ' Checking "special" positions (X~Y~Z), Bs and Qs ...')
c
      naniso = 0
c
      do i=1,natoms
        if (laniso(i)) naniso = naniso + 1
        if (batom(i) .lt. 2.0) nlb = nlb + 1
        if (qatom(i) .lt. 0.01) nlq = nlq + 1
        if (batom(i) .gt. 100.0) nhb = nhb + 1
        if (qatom(i) .gt. 1.0) nhq = nhq + 1
        if (abs(qatom(i)-1.0) .gt. 0.0001) then
          nqo = nqo + 1
        end if
        dx = abs(xyz(1,i)-xyz(2,i))
        if (dx .le. cotol) then
          dy = abs (xyz(1,i)-xyz(3,i))
          dz = abs (xyz(2,i)-xyz(3,i))
          if (dy .le. cotol .and. dz .le. cotol) then
            if (nsp .le. 10) call print_atom (i)
ccc            call fvalut (' X ~ Y ~ Z :',3,xyz(1,i))
            nsp = nsp + 1
          end if
        end if
      end do
c
      if (nsp .eq. 0) then
        call prompt (' No suspicious coordinates encountered')
      else
        call jvalut (' WARNING - Nr of atoms with X~Y~Z :',1,nsp)
        if (nsp .gt. 10) call prompt (
     +    ' (Only the first 10 were listed)')
      end if
c
      if (nhb .gt. 0) then
        call jvalut (
     +  ' WARNING - Nr of atoms with B > 100 A2 :',1,nhb)
      else
        call prompt (' All atoms have B <= 100 A2')
      end if
c
      if (nlb .gt. 0) then
        call jvalut (
     +  ' WARNING - Nr of atoms with B < 2.0 A2 :',1,nlb)
      else
        call prompt (' All atoms have B >= 2.0 A2')
      end if
c
      if (nhq .gt. 0) then
        call jvalut (
     +  ' WARNING - Nr of atoms with Q > 1.0 :',1,nhq)
      else
        call prompt (' All atoms have Q <= 1.0')
      end if
c
      if (nlq .gt. 0) then
        call jvalut (
     +  ' WARNING - Nr of atoms with Q < 0.01 :',1,nlq)
      else
        call prompt (' All atoms have Q >= 0.01')
      end if
c
      if (nqo .gt. 0) then
        call jvalut (
     +  ' NOTE - Nr of atoms with Q not equal to 1.0 :',1,nqo)
      else
        call prompt (' All atoms have Q = 1.0')
      end if
c
      if (naniso .gt. 0) then
        call jvalut (
     +  ' Nr of atoms with ANISOU cards :',1,naniso)
      else
        call prompt (' No atoms with ANISOU cards')
      end if
c
      if (nter .gt. 0) then
        call jvalut (
     +  ' Nr of TER cards :',1,nter)
      else
        call prompt (' No TER cards found')
      end if
c
c ... get secondary structure assignment
c
      call yasspa (rbuf)
c
      return
      end
c
c
c
      subroutine write_pdb (iunit,file,form,which)
c
      include 'moleman2.incl'
c
      real ca2fa(12)
c
      integer i,j,iunit,ierr,nlines,nok,leng1,iter,niter
c
      logical lhydro,mainch
c
      character*(*) file,form,which
      character key*6,line*256,myatm*4,myres*3
c
code ...
c
      call upcase (form)
      call upcase (which)
      call remspa (form)
      call remspa (which)
c
      if (form(1:1) .ne. 'X' .and. form(1:1) .ne. 'C') form = 'Pdb'
      if (which(1:2) .ne. 'SE' .and. which(1:2) .ne. 'NO' .and.
     +    which(1:2) .ne. 'PA' .and. which(1:2) .ne. 'PG' .and.
     +    which(1:2) .ne. 'PS' .and. which(1:2) .ne. 'CA')
     +  which = 'ALl'
c
      call textut (' Output PDB file :',file)
      call textut (' Format :',form)
      call textut (' Atoms  :',which)
c
      call xopxna (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      nok = 0
      if (which(1:2) .eq. 'AL') then
        do i=1,natoms
          lbuf(i) = .true.
        end do
        nok = natoms
      else if (which(1:2) .eq. 'NO') then
        do i=1,natoms
          lbuf(i) = (.not. lhydro(atmnam(i)))
          if (lbuf(i)) nok = nok + 1
        end do
      else if (which(1:2) .eq. 'SE') then
        do i=1,natoms
          lbuf(i) = (select(i))
          if (lbuf(i)) nok = nok + 1
        end do
      else if (which(1:2) .eq. 'CA') then
        do i=1,natoms
          lbuf(i) = (atmnam(i) .eq. ' CA ')
          if (lbuf(i)) nok = nok + 1
        end do
      else if (which(1:2) .eq. 'PG') then
        do i=1,natoms
          if (restyp(resptr(i)) .eq. iprot) then
            lbuf(i) = (mainch(atmnam(i)) .or. lmain(i))
            if (lbuf(i)) nok = nok + 1
          end if
        end do
      else if (which(1:2) .eq. 'PA') then
        do i=1,natoms
          if (restyp(resptr(i)) .eq. iprot) then
            lbuf(i) = (mainch(atmnam(i)) .or. lmain(i) .or.
     +                 atmnam(i) .eq. ' CB ')
            if (lbuf(i)) nok = nok + 1
          end if
        end do
      else if (which(1:2) .eq. 'PS') then
        do i=1,natoms
          if (restyp(resptr(i)) .eq. iprot) then
            lbuf(i) = (mainch(atmnam(i)) .or. lmain(i) .or.
     +        atmnam(i) .eq. ' CB ' .or. atmnam(i) .eq. ' CG ' .or.
     +        atmnam(i) .eq. ' OG ' .or. atmnam(i) .eq. ' SG ' .or.
     +        atmnam(i) .eq. ' OG1' .or. atmnam(i) .eq. ' CG1')
            if (lbuf(i)) nok = nok + 1
          end if
        end do
      end if
c
      call jvalut (' Number of atoms to write :',1,nok)
      if (nok .le. 0) then
        call errcon ('No atoms to write')
        return
      end if
c
 6000 format (a6,a)
 6010 format (a)
 4713 format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a40)
c
      nlines = 0
      nok = 0
c
      call stamp (line)
      write (iunit,'(a6,1x,a)',err=9000) 'REMARK',line(1:leng1(line))
      nlines = 1
c
c ... if PDB, write comments etc.
c
      if (form(1:1) .eq. 'P') then
        if (nrem .gt. 0) then
          do i=1,nrem
            write (iunit,6010,err=9000) remark(i)(1:leng1(remark(i)))
          end do
        end if
        if (nother .gt. 0) then
          do i=1,nother
            write (iunit,6010,err=9000) other(i)(1:leng1(other(i)))
          end do
        end if
        nlines = nlines + nrem + nother
      end if
c
c ... PDB/CCP4 - write CRYST1 etc. cards
c
      if (form(1:1) .eq. 'P' .or. form(1:1) .eq. 'C') then
c
        write (iunit,'(a6,3f9.3,3f7.2,1x,a11,i4)',err=9000) 
     +      'CRYST1',(cell(i),i=1,6),spgrp,zmol
c
        write (iunit,'(a)',err=9000)
     +      'ORIGX1      1.000000  0.000000  0.000000        0.00000'
        write (iunit,'(a)',err=9000)
     +      'ORIGX2      0.000000  1.000000  0.000000        0.00000'
        write (iunit,'(a)',err=9000)
     +      'ORIGX3      0.000000  0.000000  1.000000        0.00000'
c
c
        do i=10,12
          ca2fa (i) = 0.0
        end do
        call orthog (cell,ca2fa,1)
        write (iunit,'(a6,4x,3f10.6,f15.5)',err=9000)
     +      'SCALE1',(ca2fa(i),i=1,10,3)
        write (iunit,'(a6,4x,3f10.6,f15.5)',err=9000)
     +      'SCALE2',(ca2fa(i),i=2,11,3)
        write (iunit,'(a6,4x,3f10.6,f15.5)',err=9000)
     +      'SCALE3',(ca2fa(i),i=3,12,3)
c
        nlines = nlines + 7
      end if
c
      iter  = 1
      niter = terptr (iter)
c
c ... write atoms
c
      do i=1,natoms
c
c ... check if TER record needs to be written
c
        if (i .gt. niter) then
          write (iunit,6000,err=9000) 'TER   ',
     +      terrec(iter)(1:leng1(terrec(iter)))
          nlines = nlines + 1
          iter   = iter + 1
          niter  = terptr (iter)
          if (iter .gt. nter) niter = natoms + 1
        end if
c
        if (lbuf(i)) then
          myres = resnam(i)
          myatm = atmnam(i)
          if (which(1:2) .eq. 'PG') then
            myres = 'GLY'
          else if (which(1:2) .eq. 'PA') then
            if (myres .ne. 'GLY') myres = 'ALA'
          else if (which(1:2) .eq. 'PS') then
            if (myres .ne. 'GLY' .and. myres .ne. 'ALA') then
              myres = 'SER'
              if (atmnam(i) .eq. ' CG ' .or.
     +            atmnam(i) .eq. ' OG ' .or.
     +            atmnam(i) .eq. ' SG ' .or.
     +            atmnam(i) .eq. ' CG1' .or.
     +            atmnam(i) .eq. ' OG1') then
                myatm = ' OG '
              end if
            end if
          end if
c
          key = 'ATOM  '
          if (form(1:1) .eq. 'P') then
            if (lhet(i)) key = 'HETATM'
          end if
c
          nok = nok + 1
          write (line,4713) atomnr(i),myatm,altloc(i),
     +      myres,achain(i),iresid(i),insert(i),(xyz(j,i),j=1,3),
     +      qatom(i),batom(i),inote(i)
          write (iunit,6000,err=9000) key,line(1:leng1(line))
          nlines = nlines + 1
c
          if (laniso(i)) then
            write (line,4715) atomnr(i),myatm,altloc(i),
     +        myres,achain(i),iresid(i),insert(i),
     +        (anisou(j,i),j=1,6),inote(i)(5:)
            write (iunit,6000,err=9000) 'ANISOU',line(1:leng1(line))
            nlines = nlines + 1
          end if
 4715 format (i5,1x,a4,a1,a3,1x,a1,i4,a1,1x,6i7,a)
c
cATOM     74  NZ XLYS     8      17.216  31.254 -22.466  0.55 28.30           N
cANISOU   74  NZ XLYS     8     2886   2496   5105    713    930     17       N
cATOM     74  NZ XLYS     8      17.216  31.254 -22.466  0.55 28.30           N
cANISOU   74  NZ XLYS     8     2886   2496   5105    713    930     17           N
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8
c
        end if
      end do
c
 8000 continue
c
c ... END card
c
      write (iunit,6000,err=9000) 'END   '
      nlines = nlines + 1
c
      call jvalut (' Nr of atoms written :',1,nok)
      call jvalut (' Nr of lines written :',1,nlines)
c
      return
c
 9000 continue
      call errcon ('While writing file')
      goto 8000
c
      end
c
c
c
      subroutine split_pdb (iunit,junit,prefix)
c
      include 'moleman2.incl'
c
      real myxyz(3)
      real dist,dd
c
      integer i,j,iunit,junit,ierr,nlines,nok,length,nfiles,inow
      integer ncys,nss,ii,jj,leng1,ip,k
c
      character*(*) prefix
      character key*6,line*256,file*256,pseg*4
c
code ...
c
      call remspa (prefix)
      if (length(prefix) .lt. 1) prefix = './'
c
      call prompt (' Split PDB files for X-PLOR')
      call textut (' File prefix :',prefix)
c
      file = prefix(1:leng1(prefix))//'_generate.inp'
      call textut (' X-PLOR generate input file :',file)
      call xopxua (junit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      write (junit,6030,err=9000)
     +  ' remarks File',file(1:leng1(file)),
     +  '- generate pdb/psf file'
      call stamp (line)
      write (junit,6030,err=9000)
     +  ' remarks',line(1:leng1(line))
c
      write (junit,6030,err=9000)
      write (junit,6030,err=9000) ' topology'
      write (junit,6030,err=9000) '   @tophcsdx.pro'
      write (junit,6030,err=9000) '   @toph19.sol'
      write (junit,6030,err=9000) '   { add more topology files here }'
      write (junit,6030,err=9000) ' end'
c
      write (junit,6030,err=9000)
      write (junit,6030,err=9000) ' parameter'
      write (junit,6030,err=9000) '   @parhcsdx.pro'
      write (junit,6030,err=9000) '   @param19.sol'
      write (junit,6030,err=9000) '   { add more parameter files here }'
      write (junit,6030,err=9000)
      write (junit,6030,err=9000) '   nbonds'
      write (junit,6030,err=9000)
     +  '     atom cdie shift eps=8.0  e14fac=0.4'
      write (junit,6030,err=9000)
     +  '     cutnb=7.5 ctonnb=6.0 ctofnb=6.5'
      write (junit,6030,err=9000) '     nbxmod=5 vswitch wmin=0.5'
      write (junit,6030,err=9000) '   end'
      write (junit,6030,err=9000)
     +  '   { dielectric constant set to 8.0 (EPS) }'
      write (junit,6030,err=9000)
     +  '   { close contacts printed only if dist < 0.5 A (WMIN) }'
      write (junit,6030,err=9000) ' end'
c
      do i=1,natoms
        lbuf(i) = .false.
      end do
c
      inow = 1
      nfiles = 0
c
 6000 format (a6,a)
 6010 format (a)
 6020 format (/' File nr ',i6,' Segment |',a4,'| File name : ',a)
 6030 format (20(a,1x))
 4713 format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a40)
c
   10 continue
      pseg = inote(inow)(7:10)
c
      file = '_'//pseg//'.pdb'
      call remspa (file)
      call locase (file)
      file = prefix(1:leng1(prefix))//file(1:leng1(file))
      nfiles = nfiles + 1
      write (*,6020) nfiles,pseg,file(1:leng1(file))
      call print_atom (inow)
c
      call xopxua (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      line = 'REMARK '//line
      call pretty (line)
      write (iunit,6010,err=9000) line(1:leng1(line))
c
      line = 'REMARK File name : '//file
      call pretty (line)
      write (iunit,6010,err=9000) line(1:leng1(line))
      nlines = 2
c
c ... update X-PLOR file
c
      j = resptr(inow)
      if (restyp(j) .eq. iprot) then
        write (junit,6010,err=9000)
        write (junit,6010,err=9000) ' { protein }'
        write (junit,6010,err=9000) ' segment name="'//pseg//'"'
        write (junit,6010,err=9000) '   chain @toph19.pep'
        write (junit,6010,err=9000) '     coordinates @'//
     +    file(1:leng1(file))
        write (junit,6010,err=9000) '   end'
        write (junit,6010,err=9000) ' end'
        write (junit,6010,err=9000)
     +    ' vector do (name="CD1") ( name CD and resname ile )'
        write (junit,6010,err=9000) ' coordinates @'//
     +    file(1:leng1(file))
      else
        write (junit,6010,err=9000)
        write (junit,6010,err=9000) ' { '//libcat(restyp(j))//' }'
        write (junit,6010,err=9000) ' segment name="'//pseg//'"'
        write (junit,6010,err=9000) '   chain'
        write (junit,6010,err=9000) '     coordinates @'//
     +    file(1:leng1(file))
        write (junit,6010,err=9000) '   end'
        write (junit,6010,err=9000) ' end'
        write (junit,6010,err=9000) ' coordinates @'//
     +    file(1:leng1(file))
      end if
c
      nok = 0
c
      do i=1,natoms
        if (inote(i)(7:10) .eq. pseg) then
          key = 'ATOM  '
c
c ... check OT1/OT2
c
          if (atmnam(i) .eq. ' O  ' .or.
     +        atmnam(i) .eq. ' OT1') then
            j = resptr(i)
            if (restyp(j) .eq. iprot) then
ccc       print *,'i,j,natoms',i,j,natoms
c
c ... gjk@960606 - fix next line to avoid INOTE(0)
c
              if (inote(atmptr(2,j)+1)(7:10) .ne. pseg) then
c
                write (line,4713) atomnr(i),' OT1',altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            (xyz(k,i),k=1,3),qatom(i),batom(i),inote(i)
                write (iunit,6000,err=9000) key,line(1:leng1(line))
                nok = nok + 1
                nlines = nlines + 1
                lbuf (i) = .true.
                call prompt (' Found & written OT1')
c
                k = atmptr(1,j)
ccc                print *,i,j,k,atmptr(1,j),atmptr(2,j)
                call suggot (atmptr(2,j)-atmptr(1,j)+1,
     +            atmnam(k),xyz(1,k),ip,myxyz,.false.)
c
                if (ip .le. 0) then
                  call fvalut (' ... Adding OT2 at :',3,myxyz)
                  write (line,4713) atomnr(i)+1,' OT2',altloc(i),
     +              resnam(i),achain(i),iresid(i),insert(i),
     +              (myxyz(k),k=1,3),qatom(i),batom(i),inote(i)
                  write (iunit,6000,err=9000) key,line(1:leng1(line))
                  nok = nok + 1
                  nlines = nlines + 1
                end if
c
                goto 1234
              end if
            end if
          else if (atmnam(i) .eq. ' OT2' .or.
     +             atmnam(i) .eq. ' OTX' .or.
     +             atmnam(i) .eq. ' OXT') then
            write (line,4713) atomnr(i),' OT2',altloc(i),
     +        resnam(i),achain(i),iresid(i),insert(i),
     +        (xyz(k,i),k=1,3),qatom(i),batom(i),inote(i)
            write (iunit,6000,err=9000) key,line(1:leng1(line))
            nok = nok + 1
            nlines = nlines + 1
            lbuf (i) = .true.
            call prompt (' Found & written OT2')
            goto 1234
          end if
c
          write (line,4713) atomnr(i),atmnam(i),altloc(i),
     +      resnam(i),achain(i),iresid(i),insert(i),(xyz(j,i),j=1,3),
     +      qatom(i),batom(i),inote(i)
          write (iunit,6000,err=9000) key,line(1:leng1(line))
          nok = nok + 1
          nlines = nlines + 1
          lbuf (i) = .true.
c
 1234     continue
        end if
      end do
c
c ... X-PLOR doesn't like REMARK cards in between ATOM cards
c
c      write (line,*) 'REMARK Nr of atoms in file : ',nok
c      call pretty (line)
c      write (iunit,6010,err=9000) line(1:leng1(line))
c      nlines = nlines + 1
c
      write (iunit,6000,err=9000) 'END   '
      nlines = nlines + 1
      close (iunit)
c
      call jvalut (' Nr of lines written :',1,nlines)
      call jvalut (' Nr of atoms written :',1,nok)
      if (nok .lt. 1) then
        call errcon ('Empty file ???!!!')
      end if
c
c ... find next segment id (if any)
c
      do i=1,natoms
        if (.not. lbuf(i)) then
          inow = i
          goto 10
        end if
      end do
c
c ... all done
c
      write (*,*)
      call jvalut (' Nr of PDB files generated :',1,nfiles)
c
c ... find disulfides (if any)
c
      call prompt ('0Looking for disulfides ...')
      call find_type (natoms,resnam,atmnam,'CYS',' SG ',
     +  ncys,ibuf,.true.)
c
 6040 format (' Disulfide # ',i4,1x,i5,1x,a4,' <-> ',i5,1x,a4,
     +  ' @ ',f6.2,' A')
c
      call jvalut (' Nr of CYS SG atoms :',1,ncys)
c
      if (ncys .gt. 1) then
        call fvalut (' Max SG-SG distance for link :',1,mxcyss)
        nss = 0
        do i=1,ncys-1
          ii = ibuf(i)
          do j=i+1,ncys
            jj = ibuf(j)
            dd = dist (ii,jj,xyz)
            if (dd .le. mxcyss) then
              nss = nss + 1
              write (*,6040) nss,iresid(ii),inote(ii)(7:10),
     +          iresid(jj),inote(jj)(7:10),dd
              if (nss .eq. 1) then
                write (junit,6010,err=9000)
                write (junit,6010,err=9000) ' { the disulfides }'
              else
                write (junit,6010,err=9000)
              end if
              write (junit,6010,err=9000) ' patch DISU'
              write (line,*) '   refer=1=(segid="'//
     +          inote(ii)(7:10)//'" and resid ',iresid(ii),')'
              write (junit,6010,err=9000) line(1:leng1(line))
              write (line,*) '   refer=2=(segid="'//
     +          inote(jj)(7:10)//'" and resid ',iresid(jj),')'
              write (junit,6010,err=9000) line(1:leng1(line))
              write (junit,6010,err=9000) ' end'
            end if
          end do
        end do
        call jvalut (' Nr of disulfides :',1,nss)
      end if
c
c ... handle hydrogens
c
      write (junit,6030,err=9000)
      write (junit,6030,err=9000)
     +  '{ do not use hydrogens and remove unknown atoms }'
       write (junit,6030,err=9000)
     +  ' delete selection=( hydrogen  ) end'
       write (junit,6030,err=9000)
     +  ' delete selection=( not known ) end'
c
c      write (junit,6030,err=9000)
c      write (junit,6030,err=9000) '! flags exclude vdw end'
c      write (junit,6030,err=9000) '! hbuild'
c      write (junit,6030,err=9000)
c     +  '!   selection=(hydrogen and not known)'
c      write (junit,6030,err=9000) '!   phistep=45'
c      write (junit,6030,err=9000) '! end'
c
c      write (junit,6030,err=9000)
c      write (junit,6030,err=9000)
c     +  ' { optimise hydrogens to get rid of clashes }'
c      write (junit,6030,err=9000)
c     +  '! constraints fix=( not hydrogen ) end'
c      write (junit,6030,err=9000)
c     +  '! flags include vdw end'
c      write (junit,6030,err=9000)
c     +  '! minimize powell nstep=50 end'
c      write (junit,6030,err=9000)
c     +  '! constraints fix=( not all ) end'
c
      file = prefix(1:leng1(prefix))//'_gen.pdb'
      write (junit,6030,err=9000)
      write (junit,6030,err=9000) ' write coordinates output='//
     +  file(1:leng1(file)),' end'
c
      file = prefix(1:leng1(prefix))//'.psf'
      write (junit,6030,err=9000)
      write (junit,6030,err=9000) ' write structure output='//
     +  file(1:leng1(file)),' end'
      write (junit,6030,err=9000)
      write (junit,6030,err=9000) ' stop'
c
      write (*,*)
      call prompt (' X-PLOR generate input file written')
      close (junit)
c
      return
c
 9000 continue
      call errcon ('While writing file')
      close (iunit)
      close (junit)
c
      return
c
      end
c
c
c
      subroutine set_defaults ()
c
      include 'moleman2.incl'
c
      real dx
c
code ...
c
      blimlo = 2.0
      blimhi = 50.0
      qlimlo = 0.0
      qlimhi = 1.0
      mxcaca = 4.5
      mxpp   = 8.0
      mxbond = 2.0
      mxnonb = 3.6
      mxcyss = 2.2
      tortol = 6.0
      largeb = 0.05
      largea = 5.0
      iseed  = 12345
      yalcut = 0.5
      ybecut = 0.8
      ylhcut = 0.3
      call gkrand (dx,0.0,1.0,iseed)
c
      return
      end
c
c
c
      subroutine user_defaults (imode,par,val)
c
      include 'moleman2.incl'
c
      real dummy,dx
c
      integer mode,imode,ierr
c
      character par*(*),val*(*),line*256
c
code ...
c
      mode = imode
c
      if (mode .ne. 1) mode = 0
c
      if (mode .eq. 1) then
        call remspa (par)
        call upcase (par)
        call str2r (val,dummy,ierr)
        if (ierr .ne. 0) return
        write (line,'(3a,f10.3)') ' Set ',par(1:6),' to ',dummy
        call pretty (line(2:))
        call prompt (line)
      end if
c
      if (mode .eq. 0) then
        call fvalut (
     +    ' BLIMLO = B-factor default minimum     :',1,blimlo)
        call fvalut (
     +    ' BLIMHI = B-factor default maximum     :',1,blimhi)
        call fvalut (
     +    ' QLIMLO = Occupancy default minimum    :',1,qlimlo)
        call fvalut (
     +    ' QLIMHI = Occupancy default maximum    :',1,qlimhi)
        call fvalut (
     +    ' MXCACA = Max connected CA-CA distance :',1,mxcaca)
        call fvalut (
     +    ' MXPP   = Max connected P-P distance   :',1,mxpp)
        call fvalut (
     +    ' TORTOL = Torsion/improper tolerance   :',1,tortol)
        call fvalut (
     +    ' MXBOND = Max bonded-atom distance     :',1,mxbond)
        call fvalut (
     +    ' MXNONB = Max non-bonded distance      :',1,mxnonb)
        call fvalut (
     +    ' MXCYSS = Max disulfide S-S distance   :',1,mxcyss)
        call fvalut (
     +    ' LARGEB = Large bond distane range     :',1,largeb)
        call fvalut (
     +    ' LARGEA = Large bond angle range       :',1,largea)
        call ivalut (
     +    ' ISEED  = Random-number seed           :',1,iseed)
        call fvalut (
     +    ' YALCUT = YASSPA alpha-helix cut-off   :',1,yalcut)
        call fvalut (
     +    ' YBECUT = YASSPA beta-strand cut-off   :',1,ybecut)
        call fvalut (
     +    ' YLHCUT = YASSPA left-handed cut-off   :',1,ylhcut)
      else
        if (par(1:6) .eq. 'BLIMLO') then
          blimlo = dummy
        else if (par(1:6) .eq. 'BLIMHI') then
          blimhi = dummy
        else if (par(1:6) .eq. 'QLIMLO') then
          qlimlo = dummy
        else if (par(1:6) .eq. 'QLIMHI') then
          qlimhi = dummy
        else if (par(1:6) .eq. 'MXCACA') then
          mxcaca = dummy
        else if (par(1:6) .eq. 'MXPP  ') then
          mxpp = dummy
        else if (par(1:6) .eq. 'TORTOL') then
          tortol = dummy
        else if (par(1:6) .eq. 'MXBOND') then
          mxbond = dummy
        else if (par(1:6) .eq. 'MXNONB') then
          mxnonb = dummy
        else if (par(1:6) .eq. 'MXCYSS') then
          mxcyss = dummy
        else if (par(1:6) .eq. 'LARGEB') then
          largeb = dummy
        else if (par(1:6) .eq. 'LARGEA') then
          largea = dummy
        else if (par(1:6) .eq. 'ISEED ') then
          iseed = nint(dummy)
          call gkrand (dx,0.0,1.0,iseed)
        else if (par(1:6) .eq. 'YALCUT') then
          yalcut = dummy
        else if (par(1:6) .eq. 'YBECUT') then
          ybecut = dummy
        else if (par(1:6) .eq. 'YLHCUT') then
          ylhcut = dummy
        else
          call errcon ('Name not recognised')
        end if
      end if
c
      return
      end
c
c
c
      subroutine yasspa (resxyz)
c
      include 'moleman2.incl'
c
      real resxyz(3,maxres),calpha(3,5),cbeta(3,5),clefth(3,4)
      real rt(12),rmsd1,rmsd2,rmsd3,x
c
      integer cnt(-1:3),i,j,ierr
c
      character line*128
c
c ...
c
      call salpha (calpha)
      call sbeta (cbeta)
      call slefth (clefth)
c
      do i=1,nres
        if (captr (i) .gt. 0) then
          sstype (i) = 0
          j = captr (i)
          resxyz (1,i) = xyz (1,j)
          resxyz (2,i) = xyz (2,j)
          resxyz (3,i) = xyz (3,j)
        else
          sstype (i) = -1
          resxyz (1,i) = 999.99
          resxyz (2,i) = 999.99
          resxyz (3,i) = 999.99
        end if
      end do
c
      do i=1,nres
c
        if (captr(i) .le. 0) goto 10
c
        if (i .lt. 3 .or. i .gt. (nres-2)) goto 10
c
        if (captr(i-2) .le. 0) goto 10
        if (captr(i-1) .le. 0) goto 10
        if (captr(i+1) .le. 0) goto 10
        if (captr(i+2) .le. 0) goto 10
c
        call lsqgjk (calpha,resxyz(1,i-2),5,rmsd1,rt,ierr)
        call lsqgjk (cbeta, resxyz(1,i-2),5,rmsd2,rt,ierr)
        call lsqgjk (clefth,resxyz(1,i-2),4,rmsd3,rt,ierr)
        if (rmsd1 .le. yalcut) then
          sstype (i-1) = 1
          sstype (i)   = 1
          sstype (i+1) = 1
        else if (rmsd2 .le. ybecut) then
          sstype (i-1) = 2
          sstype (i)   = 2
          sstype (i+1) = 2
        else if (rmsd3 .le. ylhcut) then
          sstype (i-1) = 3
          sstype (i)   = 3
        end if
c
   10   continue
c
      end do
c
      do i=-1,3
        cnt (i) = 0
      end do
c
      do i=1,nres
        j = sstype (i)
        cnt (j) = cnt (j) + 1
      end do
c
      call prompt ('0YASSPA results:')
      do i=-1,3
        line = ' Nr residues of type '//ssenam(i)//' : '
        call jvalut (line,1,cnt(i))
      end do
c
      j = cnt(0)+cnt(1)+cnt(2)+cnt(3)
      if (j .gt. 0) then
        x = 100.0 * float(cnt(1)) / float(j)
        call fvalut (' %-age alpha :',1,x)
        x = 100.0 * float(cnt(2)) / float(j)
        call fvalut (' %-age beta  :',1,x)
        x = 100.0 * float(cnt(3)) / float(j)
        call fvalut (' %-age lefth :',1,x)
      end if
c
      return
      end
