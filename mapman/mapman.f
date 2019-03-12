      program mapman
c
c ... MAP MANipulation
c
c ... Gerard Kleywegt @ 930402
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'MAPMAN', vers = '080625/7.8.5')
c
      include 'maxdim.incl'
c
      integer maxsiz, maxmap
      parameter (maxsiz = maxgk1)
      parameter (maxmap = maxgk4)
c
c      pointer (iaptr,mapa)
c      pointer (ibptr,buffer)
c
c      real mapa(1),buffer(1)
c      integer malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr
      integer fmalloc
#endif
c
      integer nb, mapsize, nummaps, i1, i2
c
      logical lretry,ldone
c
code ...
c
      call gainit (prognm,vers)
c
c ... initialise history
c
      call dohist ('*INIT*',ldone)
c
      mapsize = maxsiz
      call extint ('MAPSIZE',mapsize)
      mapsize = max ( mapsize , minsiz )
      call jvalut (' Allocate maps of size :',1,mapsize)
c
      nummaps = 2
      call extint ('NUMMAPS',nummaps)
      nummaps = max ( min ( maxmap, nummaps ), 1 )
      call jvalut (' Max number of maps    :',1,nummaps)
      write (*,*)
c
c ... WRDBYT accounts for 4 or 8 bytes per word
c
   10 continue
      nb = wrdbyt*mapsize*nummaps
      iaptr = fmalloc (nb)
      nb = wrdbyt*mapsize
      ibptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
c ... BUFFER and SHADOW were equivalenced in the old version of MAPMAN
c     without dynamic memory allocation (they should never both be used
c     at the same time); passing BUFFER for both has the same net effect
c     (i.e., using the same memory addresses for two different arrays)
c
      lretry = .false.
      i1 = mapsize
      i2 = nummaps
c
      call domapm (%val(iaptr),%val(ibptr),%val(ibptr),
     +             mapsize, nummaps, lretry, i1, i2)
c
      call ffree (iaptr)
      call ffree (ibptr)
c
      if (lretry) then
        mapsize = i1
        nummaps = i2
        mapsize = max ( mapsize , minsiz )
        nummaps = max ( min ( maxmap, nummaps ), 1 )
        write (*,*)
        call jvalut (' Allocate maps of size :',1,mapsize)
        call jvalut (' Max number of maps    :',1,nummaps)
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
      subroutine domapm (map,buffer,shadow,maxpnt,maxmap,lretry,i1,i2)
c
c ... MAP MANipulation
c
c ... Gerard Kleywegt @ 930402
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap),buffer(maxpnt)
      integer shadow(maxpnt)
c
      integer maxopt,ninfmt,nutfmt,maxhis
      parameter (maxopt = 15, ninfmt=23, nutfmt=15)
      parameter (maxhis = maxopt - 2)
c
      real total,user,sys,rdum,cellv,voxv,gridv,perc,zdum
      real xdum,ydum,pdum,rhis(maxhis),xpr,xpl,alpha,beta,gamma
      real a(3,3),b(3,3),flo(3),fhi(3),xlo(3),xhi(3),frlim(3)
      real frcar(3),fract(3),newcel(6)
      real vrdist,rgbbg1,rgbbg2,rgbbg3,rgbfg1,rgbfg2,rgbfg3
      real rr,gg,bb,vrlo,vrhi,vrlev,set1st,setdx,setdy,setdz
      real rsfrad,mfrad
c
      integer length,whichm,iindex
      integer frind(3),gklim(2,3),newgrd(3),newori(3),newext(3)
      integer i,nopt,ierr,imap,j,jmap,ndum,iunit,nempty,ifmt,kmap
      integer nhis(maxhis+1),numhis,itype,iplane,idum,jdum,ivrml
      integer jplane,ncube,kuse,leng1,munit,i1,i2,npos,nneg,ibox
      integer iox,ioy,ioz,junit,itra
c
      logical xinter,linter,simap,unsave,samgrd
      logical ldone,linit,lretry,lecho,lvrml,lline
c
      character line*256,pro*10,optpar(maxopt)*80,reply*1
      character inpfmt(ninfmt)*10,outfmt(nutfmt)*10
      character gkuplo(2)*5,gkxyz(3)*1,frtype*1,inimac*128
      character cavrml*4,vrmlbg*25,vrmldc*25,vrfile*128
c
      data inpfmt /'PROTEIN   ','FFT-Y     ','TENEYCK2  ', 
     +             'CCP4      ','X-PLOR    ','OLDEZD    ',
     +             'MASK      ','NEWEZD    ','BINXPLOR  ',
     +             'BRICK     ','DSN6      ','3DMATRIX  ',
     +             'TNT       ','PHASES    ','FSMASK    ',
     +             'BRIX      ','XPLOR     ','CNS       ',
     +             'EZD       ','EM08      ','OMAP      ',
     +             'MPI       ','AMBER     '/
c
      data outfmt /'CCP4      ','OLDEZD    ','MASK      ',
     +             'NEWEZD    ','ENVELOPE  ','X-PLOR    ',
     +             'DSN6      ','BRIX      ','XPLOR     ',
     +             'CNS       ','EZD       ','AMORE     ',
     +             'OMAP      ','TURBO     ','MPI       '/
c
      data gkuplo /'Lower','Upper'/
      data gkxyz  /'X','Y','Z'/
c
      data newori /3*0/, newgrd /3*100/, newext /3*100/
      data newcel /3*100.0, 3*90.0/
c
code ...
c
      lretry = .false.
      lecho = .false.
c
c --- INITIALISATION
c
ccc      call gainit (prognm,vers)
c
      write (*,6000) maxmap,maxpnt,maxpk,maxrho,
     +  maxbop,maxbob,maxboc,krecl
c
 6000 format (
     +  ' Max nr of maps in memory     : ',i10/
     +  ' Max nr of points in map      : ',i10/
     +  ' Max nr of picked peaks       : ',i10/
     +  ' Max size of density buffer   : ',i10/
     +  ' Max nr of BONES points       : ',i10/
     +  ' Max nr of BONES branches     : ',i10/
     +  ' Max nr of BONES connections  : ',i10/
     +  ' Record length for DSN6/BRIX  : ',i10)
c
      write (*,*)
      call asciut (' Input  formats :',ninfmt,inpfmt)
      write (*,*)
      call asciut (' Output formats :',nutfmt,outfmt)
c
      nmap  = 0
c
      iunit = 9
      junit = 8
      nempty = 0
c
c ... user input unit (5=interactive; other=macro)
c
      munit = 5
c
      do i=1,maxmap
        incore (i) = .false.
        change (i) = .false.
c
        commnt (i) = 'No comment'
c
c        write (line,'(a,i2,a)') 'm',i,'.map'
c        call remspa (line)
        file (i) = 'not_defined'
c
        write (line,'(a,i2)') 'M#@$!',i
        call remspa (line)
        name (i) = line
c
        spaceg (i) = 1
        uvw (1,i) = 2
        uvw (2,i) = 1
        uvw (3,i) = 3
        maprod (i) = 1.0
        maplus (i) = 0
c
      end do
c
      do i=1,3
c
        do j=1,maxmap
          cell (i,j)   = 100.0
          cell (i+3,j) =  90.0
          origin (i,j) = 0
          extent (i,j) = 100
          grid (i,j)   = 100
        end do
c
      end do
c
c ... PICK stuff
c
      do i=1,3
        do j=1,2
          pklim(j,i) = 0
        end do
      end do
      npks = 0
      pklev = 0.0
      pkres  = 'HOH'
      pkatom = ' O  '
      pfirst = 1001
      pkfile = 'peaks.pdb'
c
c ... BONES stuff
c
      bonmap = -1
      bobase = 10.0
      bostep = 5.0
      bonex = 100
      bonlen = 5
      bobfac = 5.0
      bonfil = 'bones.odb'
      bonpdb = 'bones.pdb'
      bonid = 'skel'
      bonpcn = 0
      bonccn = 0
c
c ... FILTER stuff
c
      alpha = 0.0
      beta = 0.0
      gamma = 0.0
      ncube = 1
      kuse = 5
c
c ... ROUGHNESS stuff
c
      npos = 0
      nneg = 0
      ibox = 2
c
      do i=1,3
        gklim(1,i) = 0
        gklim(2,i) = 10
      end do
c
      frtype = 'I'
      do i=1,3
        frlim(i) = 0
      end do
c
c ... VRML stuff
c
      lvrml = .false.
      ivrml = 99
      cavrml = ' CA '
      vrdist = 4.5
      vrlo = 99.0
      vrhi = 100.0
      vrlev = 100.0
      vrmlbg = 'black'
      vrmldc = 'white'
      vrfile = 'mapman.wrl'
      rgbbg1 = 0.0
      rgbbg2 = 0.0
      rgbbg2 = 0.0
      rgbfg1 = 1.0
      rgbfg2 = 1.0
      rgbfg2 = 1.0
      call xvrml_init ()
c
      set1st = 0.0
      setdx = 0.0
      setdy = 0.0
      setdz = 0.0
c
      rsfrad = 1.5
      mfrad = 1.5
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
      do i=1,maxopt
        optpar (i) = ' '
      end do
c
      linter = xinter()
      if (linter) then
        pro='$MAPMAN > '
      else
        pro=' MAPMAN > '
      end if
c
c ... to be documented:
c
c     +  ' RS_fit obs_map calc_map pdbfile listfile [rad]'/
c     +  ' LC_map new_map map_1 map_2 [cc|r] [sph|par] [rad]'/
c
 6010 format (/
     +  ' MAPMAN''s options :'//
     +  ' ? (list options)            ! (comment)'/
     +  ' QUit                        $ shell_command'/
     +  ' & symbol value              & ? (list symbols)'/
     +  ' @ macro_file                ZP_restart mapsize nummaps'/
     +  ' ECho on_off                 # parameter(s) (command history)'/
     + /
     +  ' REad new_map file fmt       WRite map file fmt [cutoff] [op]'/
     +  ' DElete map                  LIst [map]'/
     +  ' HIsto map x1 x2 x3 [...]    EXtract map file type plane'/
     +  ' MOments map                 DH map plot_file nr_bins'/
     +  ' 1D_proj map file type       2D_proj map file type first last'/
     + /
     +  ' MAppage map file            BRix map file'/
     +  ' SWapbytes brick_file        '/
     + /
     +  ' NEw ? (list defaults)       NEw ORigin o1 o2 o3'/
     +  ' NEw CEll a b c al be ga     NEw EXtent e1 e2 e3'/
     +  ' NEw GRid g1 g2 g3           NEw SAme map'/
     +  ' NEw ENcompass map           NEw MAke new_map'/
     + /
     +  ' UVw map u v w               SPacegroup map ispcgr'/
     +  ' PRod_plus map prod plus     CEll map a b c al be ga'/
     +  ' RAnge map low high          DRange map low high'/
     +  ' MUltiply map factor         DIvide map factor'/
     +  ' PLus map value              ODl map odl_file type mode'/
     +  ' SCale map low high          ZEro map low high'/
     +  ' NOrmalise map               ERase map e1 .. e6'/
     +  ' INtegrate map e1 .. e6      FRame map type e1 e2 e3'/
     +  ' VAlue_at map e1 .. e6       CHange_value map lo hi new'/
     +  ' WHere_value map value toler SEt_value map start [dx dy dz]'/
     +  ' TRanslate map tx ty tz      GTranslate map gx gy gz'/
     +  ' LAbel map text              '/
     + /
     +  ' SImilarity map1 map2        DUplicate new_map map'/
     +  ' ADd map map_to_add          MIn_max map1 map2'/
     +  ' OPerate map oper map2       COrrelate map1 map2'/
     +  ' MStats map mask_map cutoff  QInvert new_map map'/
     +  ' PAste map1 e1 .. e6 map2 e1 e2 e3'/
     +  ' ROughness new_map map npos nneg box plot_file nr_bins'/
     + /
     +  ' BOnes SKel map base step nr_spares'/
     +  ' BOnes COnnect file id mc    BOnes ?   '/
     +  ' BOnes PRune file mc bfac    '/
     + /
     +  ' PIck PDb res atom first     PIck LImits e1 .. e6'/
     +  ' PIck LEvel threshold        PIck ?'/
     +  ' PIck PEaks map file format  PIck INtegrate map file format'/
     +  ' PIck HIgh map file format   '/
     + /
     +  ' PEek VAlue map pdb_in pdb_out mode'/
     +  ' PEek CUbe map pdb_in pdb_out mode nr_points'/
     +  ' PEek SPhere map pdb_in pdb_out mode radius'/
     + /
     +  ' FIlter ?                    FIlter EDge map a b'/
     +  ' FIlter MEdian map size      FIlter AVerage map size'/
     +  ' FIlter MInimum map size     FIlter MAximum map size'/
     +  ' FIlter SMooth map           FIlter SIgnal/noise map size'/
     +  ' FIlter GRadient map         FIlter LAplace map'/
     +  ' FIlter KHighest map size k  FIlter KLowest map size k'/
     +  ' FIlter VAr_thr map size a b FIlter STat_diff map size a b c'/
     + /
     +  ' VRml SEtup central_atom max_dist backgr_col default_col'/
     +  ' VRml INit [filename]                ',
     +              ' VRml COlour_list'/
     +  ' VRml TRace pdb_file [colour]        ',
     +              ' VRml BOx map [colour] [line_solid]'/
     +  ' VRml CEll map [colour] [line_solid] ',
     +              ' [x_offset] [y_offset] [z_offset]'/
     +  ' VRml DOts map low high [colour]     ',
     +              ' VRml DRaw map level [colour]'/
     +  ' VRml CLose_file                     ',
     +              ' '/
     + )
c
 6040 format (' WARNING - unsaved changes to map ',a,' !!!')
c
      linit = .false.
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6010)
      write (*,6000) maxmap,maxpnt,maxpk,maxrho,
     +  maxbop,maxbob,maxboc,krecl
      write (*,*)
c
c ... 950118 - check if CCP4_OPEN is defined; print warning if not
c
      if (.not. linit) call gkccp4 (.false.)
c
c --- MAIN EVENT LOOP
c
   10 continue
c
      call gkdcpu (total,user,sys)
      if (total .ge. 1.0)
     +  write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
c ... execute initialisation macro ?
c
      if (.not. linit) then
        linit = .true.
        call gknval ('GKMAPMAN',inimac,ierr)
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
          write (optpar(2),*) i1
          call textin (' New MAPSIZE ?',optpar(2))
        end if
        call str2i (optpar(2),idum,ierr)
        if (ierr .ne. 0) goto 10
        i1 = idum
        if (i1 .lt. 10) then
          call errcon ('Silly MAPSIZE')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          write (optpar(3),*) i2
          call textin (' New NUMMAPS ?',optpar(3))
        end if
        call str2i (optpar(3),idum,ierr)
        if (ierr .ne. 0) goto 10
        i2 = idum
        if (i2 .lt. 1) then
          call errcon ('Silly NUMMAPS')
          goto 10
        end if
c
        lretry = .true.
        return
c
c ... SWAPBYTES
c
      else if (optpar(1)(1:2) .eq. 'SW') then
c
        if (nopt .lt. 2) then
          call textin (' Brick file ?',optpar(2))
        end if
c
        call swappy (iunit,krecl,optpar(2),ierr)
        close (iunit)
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
        else if (optpar(2)(1:2) .eq. 'DO') then
c
c ... VRML DOTS
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
ccc            optpar (3) = prev
            call textin (' Map ?',optpar(3))
          end if
          imap = whichm (optpar(3),maxmap,maxpnt)
c
          if (imap .le. 0 .or. imap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(imap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            write (optpar(4),'(1pe12.4)') vrlo
            call textin (' Lowest level ?',optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .ne. 0) xdum = vrlo
          vrlo = xdum
c
          if (nopt .lt. 5) then
            write (optpar(5),'(1pe12.4)') vrhi
            call textin (' Highest level ?',optpar(5))
          end if
          call str2r (optpar(5),xdum,ierr)
          if (ierr .ne. 0) xdum = vrhi
          vrhi = xdum
c
          if (nopt .ge. 6) then
            call xvrml_rgb_name (optpar(6),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          call rlohi (vrlo,vrhi)
c
          call vrdots (map(1,imap),extent(1,imap),
     +             extent(2,imap),extent(3,imap),grid(1,imap),
     +             origin(1,imap),cell(1,imap),vrlo,vrhi,
     +             ivrml,ierr)
c
        else if (optpar(2)(1:2) .eq. 'DR') then
c
c ... VRML DRAW
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
ccc            optpar (3) = prev
            call textin (' Map ?',optpar(3))
          end if
          imap = whichm (optpar(3),maxmap,maxpnt)
c
          if (imap .le. 0 .or. imap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(imap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            write (optpar(4),'(1pe12.4)') vrlev
            call textin (' Level ?',optpar(4))
          end if
          call str2r (optpar(4),xdum,ierr)
          if (ierr .ne. 0) xdum = vrlev
          vrlev = xdum
c
          if (nopt .ge. 5) then
            call xvrml_rgb_name (optpar(5),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          call vrdraw (map(1,imap),extent(1,imap),
     +             extent(2,imap),extent(3,imap),grid(1,imap),
     +             origin(1,imap),cell(1,imap),vrlev,
     +             ivrml,buffer,ierr)
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
c
          if (nopt .lt. 3) then
            optpar (3) = 'm1.pdb'
            call textin (' PDB file name ?',optpar(3))
          end if
c
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          write (*,'(99a)') ' VRML - Trace of ',cavrml,
     +      ' atoms for PDB file ',optpar(3)(1:leng1(optpar(3)))
c
          call vrmlca (iunit,optpar(3),vrdist,ivrml,cavrml,ierr)
          close (iunit)
c
          if (ierr .ne. 0) then
            lvrml = .false.
            call xvrml_close ()
            call errcon ('Closed VRML file')
            goto 10
          end if
c
          call flusho (ivrml)
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
          if (nopt .lt. 3) then
ccc            optpar (3) = prev
            call textin (' Map ?',optpar(3))
          end if
          imap = whichm (optpar(3),maxmap,maxpnt)
c
          if (imap .le. 0 .or. imap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(imap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          lline = .true.
          if (nopt .ge. 5) then
            call upcase (optpar(5))
            lline = (optpar(5)(1:1).ne.'S')
          end if
c
          iox = 0
          ioy = 0
          ioz = 0
          if (nopt .ge. 6) then
            call str2i (optpar(6),idum,ierr)
            if (ierr .eq. 0) iox = idum
          end if
          if (nopt .ge. 7) then
            call str2i (optpar(7),idum,ierr)
            if (ierr .eq. 0) ioy = idum
          end if
          if (nopt .ge. 8) then
            call str2i (optpar(8),idum,ierr)
            if (ierr .eq. 0) ioz = idum
          end if
c
          write (*,'(99a)') ' VRML - Cell for map ',
     +      name(imap)
c
          call xvrml_cell (cell(1,imap),iox,ioy,ioz,lline)
c
          call flusho (ivrml)
c
        else if (optpar(2)(1:2) .eq. 'BO') then
c
c ... VRML BOX
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
ccc            optpar (3) = prev
            call textin (' Map ?',optpar(3))
          end if
          imap = whichm (optpar(3),maxmap,maxpnt)
c
          if (imap .le. 0 .or. imap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(imap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          lline = .true.
          if (nopt .ge. 5) then
            call upcase (optpar(5))
            lline = (optpar(5)(1:1).ne.'S')
          end if
c
          write (*,'(99a)') ' VRML - Box for map ',
     +      name(imap)
c
          call vbox (imap,lline)
c
          call flusho (ivrml)
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
c ... NEW
c
      else if (optpar(1)(1:2) .eq. 'NE') then
c
        call upcase (optpar(2))
c
 6050 format (/' Current defaults for the next NEW map:'/
     +  ' Grid    = ',3i10/
     +  ' Origin  = ',3i10/
     +  ' Extent  = ',3i10/
     +  ' Cell    = ',6f10.3/
     +  ' Nr of points = ',i10,' Max = ',i10)
c
c ... NEW ?
c
        if (nopt .lt. 2 .or. optpar(2)(1:1) .eq. '?') then
c
          ndum = newext(1)*newext(2)*newext(3)
          write (*,6050) (newgrd(i),i=1,3),(newori(i),i=1,3),
     +      (newext(i),i=1,3),(newcel(i),i=1,6),ndum,maxpnt
          if (ndum .gt. maxpnt)
     +      call errcon ('NEW map is too big !')
c
c ... NEW ORIGIN
c
        else if (optpar(2)(1:2) .eq. 'OR') then
c
          do i=1,3
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(i10)') newori(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for origin ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2i(optpar(i+2),ndum,ierr)
            if (ierr .ne. 0) goto 10
            newori (i) = ndum
          end do
          call ivalut (' NEW origin :',3,newori)
c
c ... NEW GRID
c
        else if (optpar(2)(1:2) .eq. 'GR') then
c
          do i=1,3
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(i10)') newgrd(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for grid ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2i(optpar(i+2),ndum,ierr)
            if (ierr .ne. 0) goto 10
            newgrd (i) = max (1, ndum)
          end do
          call ivalut (' NEW grid :',3,newgrd)
c
c ... NEW EXTENT
c
        else if (optpar(2)(1:2) .eq. 'EX') then
c
          do i=1,3
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(i10)') newext(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for extent ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2i(optpar(i+2),ndum,ierr)
            if (ierr .ne. 0) goto 10
            newext (i) = ndum
          end do
          call ivalut (' NEW extent :',3,newext)
c
c ... NEW CELL
c
        else if (optpar(2)(1:2) .eq. 'CE') then
c
          do i=1,6
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(f15.3)') newcel(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for cell ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2r(optpar(i+2),rdum,ierr)
            if (ierr .ne. 0) goto 10
            newcel (i) = rdum
          end do
          call fvalut (' NEW cell :',6,newcel)
c
c ... NEW SAME
c
        else if (optpar(2)(1:2) .eq. 'SA') then
c
          if (nopt .lt. 3) call textin (' Old map ?',optpar(3))
          jmap = whichm(optpar(3),maxmap,maxpnt)
c
          if (jmap .le. 0 .or. jmap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(jmap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          do i=1,3
            newcel (i)   = cell (i,jmap)
            newcel (i+3) = cell (i+3,jmap)
            newori (i)   = origin (i,jmap)
            newgrd (i)   = grid (i,jmap)
            newext (i)   = extent (i,jmap)
          end do
c
c ... NEW ENCOMPASS
c
        else if (optpar(2)(1:2) .eq. 'EN') then
c
          if (nopt .lt. 3) call textin (' Old map ?',optpar(3))
          jmap = whichm(optpar(3),maxmap,maxpnt)
c
          if (jmap .le. 0 .or. jmap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(jmap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          do i=1,3
            if (abs(newcel(i)-cell(i,jmap)).gt.0.01) then
              call errcon ('Map cell differs from NEw CEll')
              goto 10
            end if
            if (abs(newcel(i+3)-cell(i+3,jmap)).gt.0.01) then
              call errcon ('Map cell differs from NEw CEll')
              goto 10
            end if
            if (newgrd(i) .ne. grid(i,jmap)) then
              call errcon ('Map grid differs from NEw GRid')
              goto 10
            end if
          end do
c
          call jvalut (' Current NEw ORigin :',3,newori)
          call jvalut (' Current NEw EXtent :',3,newext)
          call jvalut (' Current map origin :',3,origin(1,jmap))
          call jvalut (' Current map extent :',3,extent(1,jmap))
c
          do i=1,3
            newcel (i)   = cell (i,jmap)
            newcel (i+3) = cell (i+3,jmap)
            newgrd (i)   = grid (i,jmap)
c
            j = max (newori(i)+newext(i)-1,
     +               origin(i,jmap)+extent(i,jmap)-1)
            newori (i)   = min (newori(i),origin(i,jmap))
            newext (i)   = j - newori(i) + 1
          end do
c
          call jvalut (' New origin :',3,newori)
          call jvalut (' New extent :',3,newext)
c
c ... NEW MAKE
c
        else if (optpar(2)(1:2) .eq. 'MA') then
c
          if (nopt .lt. 3) call textin (' New map ?',optpar(3))
          call allocm (optpar(3),imap,ierr,maxmap,maxpnt)
          if (ierr .ne. 0) goto 10
c
          ndum = newext(1)*newext(2)*newext(3)
          if (ndum .gt. maxpnt) then
            call errcon ('NEW map is too big !')
            call jvalut (' Requested map size :',1,ndum)
            call jvalut (' Max available size :',1,maxpnt)
            goto 10
          end if
c
          incore (imap) = .true.
          select (imap) = .false.
          change (imap) = .true.
          file (imap)   = 'not_defined'
          commnt (imap) = 'Created from scratch with NEw MAke'
c
          do i=1,3
            cell (i,imap)   = newcel (i)
            cell (i+3,imap) = newcel (i+3)
            grid (i,imap)   = newgrd (i)
            origin (i,imap) = newori (i)
            extent (i,imap) = newext (i)
            uvw (i,imap)    = i
          end do
          spaceg (imap) = 1
          npnt (imap) = extent(1,imap)*extent(2,imap)*extent(3,imap)
c
          call inimap (map(1,imap),extent(1,imap),extent(2,imap),
     +      extent(3,imap),0.0)
c
          call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +      edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
          edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c
          call proplu (edstat(1,imap),edstat(2,imap),maprod(imap),
     +      maplus(imap))
c
          call prompt (' New (empty) map created (P1) !')
c
        else
c
          call errcon ('Unknown NEw sub-command')
c
        end if
c
        goto 10
c
c ... FILTER
c
      else if (optpar(1)(1:2) .eq. 'FI') then
c
        call upcase (optpar(2))
c
c ... FILTER ?
c
        if (nopt .lt. 2 .or. optpar(2)(1:1) .eq. '?') then
c
          write (*,'(1x,a)')
     +  ' ','Filter options:','---------------',
     +  'MEdian  - replace point by median  value in local cube',
     +  'MInimum - replace point by minimum value in local cube',
     +  'MAximum - replace point by maximum value in local cube',
     +  'AVerage - replace point by average value in local cube',
     +  '          i.e., convolution with     1 1 1   1 1 1   1 1 1',
     +  '                            (1/27) * 1 1 1   1 1 1   1 1 1',
     +  '                                     1 1 1   1 1 1   1 1 1',
     +  'SMooth  - convolution with           1 2 1   2 3 2   1 2 1',
     +  '                            (1/55) * 2 3 2   3 5 3   2 3 2',
     +  '                                     1 2 1   2 3 2   1 2 1',
     +  'EDge    - convolution with       A A A   A   A   A   A A A',
     +  '                         (1/B) * A A A   A B-26A A   A A A',
     +  '          (A<0, B>0)             A A A   A   A   A   A A A',
     +  'SIgnal/noise - replace point by average/sigma in local cube',
     +  'GRadient     - convolution with (-1 0 +1) in 3 directions',
     +  'LAplace      - convolution with (-1 2 -1) in 3 directions',
     +  'KLowest      - replace point by average of K lowest  nbrs',
     +  'KHighest     - replace point by average of K highest nbrs',
     +  'VAr_thr      - variable threshold (A, B in <0,1]):',
     +  '                  NEW = old - (A*sigma + B*ave)',
     +  'STat_diff    - statistical differencing (A, B, C in <0,1]):',
     +  '                  NEW = (1-A)*old + (old-ave)*B*C/(C+B*sigma)',
     +  ' ','Parameters:',
     +  'Size  = one-sided size of small, local cube (usually, 1)',
     +  'A,B,C = parameters that determine the convolution matrix',
     +  '        (the program suggests sensible defaults)',
     +  'K     = number of neigbours to use in averaging',
     +  ' ','Further reading :',
     +  'Wayne Niblack, "An Introduction to Digital Image Processing",',
     +  'Prentice-Hall International, London, 1986.',' '
c
          goto 10
        end if
c
        if (nopt .lt. 3) call textin
     +    (' Map to filter ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
c ... FILTER MEDIAN MINIMUM MAXIMUM AVERAGE SIGNAL/NOISE
c
        if (optpar(2)(1:2) .eq. 'ME' .or.
     +      optpar(2)(1:2) .eq. 'MA' .or.
     +      optpar(2)(1:2) .eq. 'MI' .or.
     +      optpar(2)(1:2) .eq. 'SI' .or.
     +      optpar(2)(1:2) .eq. 'AV') then
c
          if (nopt .lt. 4) then
            optpar(4) = '1'
            call textin (' Cube size ?',optpar(4))
          end if
          call str2i (optpar(4),ncube,ierr)
          if (ierr .ne. 0) goto 10
c
          call filter (optpar(2)(1:2),map(1,jmap),buffer,ncube,
     +      extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +      alpha,beta,gamma,kuse,ierr)
c
          if (ierr. ne. 0) goto 10
c
c ... FILTER KLOWEST KHIGHEST
c
        else if (optpar(2)(1:2) .eq. 'KL' .or.
     +           optpar(2)(1:2) .eq. 'KH') then
c
          if (nopt .lt. 4) then
            optpar(4) = '1'
            call textin (' Cube size ?',optpar(4))
          end if
          call str2i (optpar(4),ncube,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            optpar(5) = '5'
            call textin (' Value for K ?',optpar(5))
          end if
          call str2i (optpar(5),kuse,ierr)
          if (ierr .ne. 0) goto 10
c
          call filter (optpar(2)(1:2),map(1,jmap),buffer,ncube,
     +      extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +      alpha,beta,gamma,kuse,ierr)
c
          if (ierr. ne. 0) goto 10
c
c ... FILTER SMOOTH GRADIENT LAPLACE
c
        else if (optpar(2)(1:2) .eq. 'SM' .or.
     +           optpar(2)(1:2) .eq. 'GR' .or.
     +           optpar(2)(1:2) .eq. 'LA') then
c
          ncube = 1
          call filter (optpar(2)(1:2),map(1,jmap),buffer,ncube,
     +      extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +      alpha,beta,gamma,kuse,ierr)
c
          if (ierr. ne. 0) goto 10
c
c ... FILTER STAT_DIFF
c
        else if (optpar(2)(1:2) .eq. 'ST') then
c
          if (nopt .lt. 4) then
            optpar(4) = '1'
            call textin (' Cube size ?',optpar(4))
          end if
          call str2i (optpar(4),ncube,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            write (optpar(5),*) 0.7
            call textin (' Value for A (<0,1]) ?',optpar(5))
          end if
          call str2r (optpar(5),alpha,ierr)
          if (ierr .ne. 0) goto 10
          if (alpha .le. 0.0 .or. alpha .gt. 1.0) then
            call errcon ('Invalid value for A')
            goto 10
          end if
c
          if (nopt .lt. 6) then
            write (optpar(6),*) 0.5
            call textin (' Value for B (<0,1]) ?',optpar(6))
          end if
          call str2r (optpar(6),beta,ierr)
          if (ierr .ne. 0) goto 10
          if (beta .le. 0.0 .or. beta .gt. 1.0) then
            call errcon ('Invalid value for B')
            goto 10
          end if
c
          if (nopt .lt. 7) then
            write (optpar(7),*) 1.0
            call textin (' Value for C (<0,1]) ?',optpar(7))
          end if
          call str2r (optpar(7),gamma,ierr)
          if (ierr .ne. 0) goto 10
          if (gamma .le. 0.0 .or. gamma .gt. 1.0) then
            call errcon ('Invalid value for C')
            goto 10
          end if
c
          call filter (optpar(2)(1:2),map(1,jmap),buffer,ncube,
     +      extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +      alpha,beta,gamma,kuse,ierr)
c
          if (ierr. ne. 0) goto 10
c
c ... FILTER VAR_THRES
c
        else if (optpar(2)(1:2) .eq. 'VA') then
c
          if (nopt .lt. 4) then
            optpar(4) = '1'
            call textin (' Cube size ?',optpar(4))
          end if
          call str2i (optpar(4),ncube,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            write (optpar(5),*) 0.18
            call textin (' Value for A (<0,1]) ?',optpar(5))
          end if
          call str2r (optpar(5),alpha,ierr)
          if (ierr .ne. 0) goto 10
          if (alpha .le. 0.0 .or. alpha .gt. 1.0) then
            call errcon ('Invalid value for A')
            goto 10
          end if
c
          if (nopt .lt. 6) then
            write (optpar(6),*) 1.0
            call textin (' Value for B (<0,1]) ?',optpar(6))
          end if
          call str2r (optpar(6),beta,ierr)
          if (ierr .ne. 0) goto 10
          if (beta .le. 0.0 .or. beta .gt. 1.0) then
            call errcon ('Invalid value for B')
            goto 10
          end if
c
          call filter (optpar(2)(1:2),map(1,jmap),buffer,ncube,
     +      extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +      alpha,beta,gamma,kuse,ierr)
c
          if (ierr. ne. 0) goto 10
c
c ... FILTER EDGE
c
        else if (optpar(2)(1:2) .eq. 'ED') then
c
          if (nopt .lt. 4) then
            write (optpar(4),*) -1.0
            call textin (' Value for A (<0) ?',optpar(4))
          end if
          call str2r (optpar(4),alpha,ierr)
          if (ierr .ne. 0) goto 10
          if (alpha .ge. 0.0) then
            call errcon ('Invalid value for A')
            goto 10
          end if
c
          if (nopt .lt. 5) then
            write (optpar(5),*) -alpha
            call textin (' Value for B (>0) ?',optpar(5))
          end if
          call str2r (optpar(5),beta,ierr)
          if (ierr .ne. 0) goto 10
          if (beta .le. 0.0) then
            call errcon ('Invalid value for B')
            goto 10
          end if
c
          ncube = 1
          call filter (optpar(2)(1:2),map(1,jmap),buffer,ncube,
     +      extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +      alpha,beta,gamma,kuse,ierr)
c
          if (ierr. ne. 0) goto 10
c
c ... more options
c
        else
          call errcon ('Invalid FILTER option')
          call textut (' Option :',optpar(2))
          goto 10
        end if
c
c ... update map statistics
c
        call xstats (map(1,jmap),npnt(jmap),edstat(3,jmap),
     +        edstat(4,jmap),edstat(1,jmap),edstat(2,jmap),pdum)
        edstat(5,jmap) = edstat(4,jmap)*edstat(4,jmap)
        change (jmap) = .true.
c
c ... changed
c
c        call dynran (maprod(jmap),maplus(jmap),pdum,rdum)
c        pdum = max (pdum, edstat(1,jmap))
c        rdum = min (rdum, edstat(2,jmap))
c        call proplu (pdum,rdum,maprod(jmap),maplus(jmap))
c
        call proplu (edstat(1,jmap),edstat(2,jmap),
     +    maprod(jmap),maplus(jmap))
c
        goto 10
c
c ... BONES
c
      else if (optpar(1)(1:2) .eq. 'BO') then
c
        call upcase (optpar(2))
c
c ... BONES ?
c
        if (nopt .lt. 2 .or. optpar(2)(1:1) .eq. '?') then
c
          if (bonmap .le. 0) then
            call prompt (' No map has been skeletonised yet !')
            write (*,6820) '***NONE***','***NONE***',
     +        bobase,bostep,
     +        bonex,bonpcn,bonccn,bonfil(1:leng1(bonfil)),
     +        bonid(1:leng1(bonid)),bonlen,bobfac
          else
            write (*,6820) name(bonmap)(1:leng1(name(bonmap))),
     +        file(bonmap)(1:leng1(file(bonmap))),bobase,bostep,
     +        bonex,bonpcn,bonccn,bonfil(1:leng1(bonfil)),
     +        bonid(1:leng1(bonid)),bonlen,bobfac
          end if
c
 6820 format (/
     +  ' Current BONES parameters and defaults :'/
     +  ' Current BONES map : ',a/
     +  ' From file         : ',a/
     +  ' Skeletonised with base level ',f8.3,' and step ',f8.3/
     +  ' And including ',i8,' spare points'/
     +  ' Yielding ',i8,' points and ',i8,' connections'/
     +  ' BONES file        : ',a/
     +  ' BONES ID for O    : ',a/
     +  ' Minimum length for main-chain fragments = ',i3/
     +  ' Default temperature factor = ',f8.3/
     +  )
c
c ... BONES SKELETONISE
c
        else if (optpar(2)(1:2) .eq. 'SK') then
c
          if (nopt .lt. 3) call textin
     +      (' Map to skeletonise ?',optpar(3))
          jmap = whichm(optpar(3),maxmap,maxpnt)
c
          if (jmap .le. 0 .or. jmap .gt. maxmap) then
            call errcon ('Invalid map selection')
            goto 10
          end if
c
          if (.not. incore(jmap)) then
            call errcon ('Map not in memory')
            goto 10
          end if
c
          if (nopt .lt. 4) then
            write (optpar(4),*) bobase
            call textin (' Base level ?',optpar(4))
          end if
c
          call str2r (optpar(4),xdum,ierr)
          if (ierr .ne. 0) goto 10
          bobase = xdum
c
          if (nopt .lt. 5) then
            write (optpar(5),*) bostep
            call textin (' Step size ?',optpar(5))
          end if
c
          call str2r (optpar(5),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .le. 0.001) then
            call errcon (' Step size must be > 0.001')
            goto 10
          end if
          bostep = xdum
c
          if (nopt .lt. 6) then
            write (optpar(6),*) bonex
            call textin (' Nr of spare points ?',optpar(6))
          end if
c
          call str2i (optpar(6),ndum,ierr)
          if (ierr .ne. 0) goto 10
          bonex = max (0, ndum)
c
          call xbones (jmap,ierr,
     +    maxmap,maxpnt,map,shadow)
          if (ierr .eq. 0) then
            bonmap = jmap
          end if
c
c ... BONES CONNECTIVITY
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          if (nopt .lt. 3) then
            optpar(3) = bonfil
            call textin (' Output BONES file ?',optpar(3))
          end if
c
          if (nopt .lt. 4) then
            optpar(4) = bonid
            call textin (' BONES ID for O ?',optpar(4))
          end if
c
          if (nopt .lt. 5) then
            write(optpar(5),*) bonlen
            call textin (' Min length main-chain fragments ?',optpar(5))
          end if
c
          if (length(optpar(3)) .lt. 1) then
            call errcon ('No BONES filename provided')
            goto 10
          end if
          bonfil = optpar(3)
c
          if (length(optpar(4)) .lt. 1) then
            call errcon ('No BONES ID provided')
            goto 10
          end if
          bonid = optpar(4)
c
          call str2i (optpar(5),ndum,ierr)
          if (ierr .ne. 0) goto 10
          if (ndum .le. 1) then
            call errcon ('Main-chain length must be > 1')
            goto 10
          end if
          bonlen = ndum
c
          call xopxua (iunit,bonfil,xinter(),ierr)
          if (ierr .ne. 0) goto 10
c
          call ybones (iunit,ierr)
c
          close (iunit)
c
c ... BONES PRUNE
c
        else if (optpar(2)(1:2) .eq. 'PR') then
c
          if (nopt .lt. 3) then
            optpar(3) = bonpdb
            call textin (' BONES PDB file ?',optpar(3))
          end if
c
          if (nopt .lt. 4) then
            write(optpar(4),*) bonlen
            call textin (' Min length main-chain fragments ?',optpar(4))
          end if
c
          if (length(optpar(3)) .lt. 1) then
            call errcon ('No BONES filename provided')
            goto 10
          end if
          bonpdb = optpar(3)
c
          call str2i (optpar(4),ndum,ierr)
          if (ierr .ne. 0) goto 10
          if (ndum .le. 1) then
            call errcon ('Main-chain length must be > 1')
            goto 10
          end if
          bonlen = ndum
c
          if (nopt .lt. 5) then
            write (optpar(5),*) bobfac
            call textin (' Temperature factor ?',optpar(5))
          end if
c
          call str2r (optpar(5),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (xdum .le. 0.001) then
            call errcon (' Temperature factor must be > 0.001')
            goto 10
          end if
          bobfac = xdum
c
          call xopxua (iunit,bonpdb,xinter(),ierr)
          if (ierr .ne. 0) goto 10
c
          call zbones (iunit,ierr)
c
          close (iunit)
c
c ... invalid sub-command
c
        else
c
          call errcon ('Invalid BONES option')
          call textut (' Selected option :',optpar(2))
c
        end if         
c
c ... PEEK
c
      else if (optpar(1)(1:2) .eq. 'PE') then
c
        if (nopt .lt. 2) then
          optpar (2) = 'VA'
          call textin (' Peek option (VAlue|CUbe|SPhere) ?',
     +      optpar(2))
        end if
c
        call upcase (optpar(2))
        if (optpar(2)(1:2) .ne. 'VA' .and.
     +      optpar(2)(1:2) .ne. 'CU' .and.
     +      optpar(2)(1:2) .ne. 'SP') then
          call errcon ('Invalid PEEK option')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which map ?',optpar(3))
        imap = whichm(optpar(3),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = 'in.pdb'
          call textin (' Input PDB file ?',optpar(4))
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'out.pdb'
          call textin (' Output PDB file ?',optpar(5))
        end if
c
        if (optpar(2)(1:2) .eq. 'VA') then
          if (nopt .lt. 6) then
            optpar(6) = 'Spline'
            call textin (' Mode (Nearest|Interpolate|Spline) ?',
     +        optpar(6))
          end if
          optpar (7) = ' '
        else if (optpar(2)(1:2) .eq. 'CU') then
          if (nopt .lt. 6) then
            optpar(6) = 'Rms'
            call textin (' Mode (Int.|Abs.int.|Rms|Mean) ?',
     +        optpar(6))
          end if
          if (nopt .lt. 7) then
            optpar(7) = '2'
            call textin (' Grid points ?',
     +        optpar(7))
          end if
        else if (optpar(2)(1:2) .eq. 'SP') then
          if (nopt .lt. 6) then
            optpar(6) = 'Rms'
            call textin (' Mode (Int.|Abs.int.|Rms|Mean) ?',
     +        optpar(6))
          end if
          if (nopt .lt. 7) then
            optpar(7) = '3.0'
            call textin (' Radius (A) ?',
     +        optpar(7))
          end if
        end if
c
        call upcase (optpar(6))
        call remspa (optpar(6))
c
        call dopeek (imap,optpar(2),optpar(4),optpar(5),
     +               optpar(6),optpar(7),maxpnt,maxmap,map)
c
c ... PICK
c
      else if (optpar(1)(1:2) .eq. 'PI') then
c
        call upcase (optpar(2))
c
c ... PICK ?
c
        if (nopt .lt. 2 .or. optpar(2)(1:1) .eq. '?') then
c
c
          write (*,6810) pkres,pkatom,pfirst,pklev,pklim
c
 6810 format (/
     +  ' Peak residue & atom names : |',a3,'|',a4,'|'/
     +  ' First peak residue number : ',i4/
     +  ' Peak threshold level      : ',1pe12.4/
     +  ' Pick limits               : ',6i6/)
c
c ... PICK LEVEL
c
        else if (optpar(2)(1:2) .eq. 'LE') then
c
          if (nopt .lt. 3) then
            write (optpar(3),*) pklev
            call remspa (optpar(3))
            call textin (' Lowest peak pick level ?',optpar(3))
          end if
          call str2r (optpar(3),rdum,ierr)
          if (ierr .ne. 0) goto 10
          pklev = rdum
c
c ... PICK PDB
c
        else if (optpar(2)(1:2) .eq. 'PD') then
c
          if (nopt .lt. 3) then
            optpar(3) = pkres
            call textin (' Residue name ?',optpar(3))
          end if
          call upcase (optpar(3))
          if (length(optpar(3)) .gt. 0) pkres = optpar(3)(1:3)
c
          if (nopt .lt. 4) then
            optpar(4) = pkatom
            call textin (' Atom name ?',optpar(4))
          end if
          call upcase (optpar(4))
          if (length(optpar(4)) .gt. 0) pkatom = optpar(4)(1:4)
c
          if (nopt .lt. 5) then
            write (optpar(5),*) pfirst
            call remspa (optpar(5))
            call textin (
     +        ' Residue number for first peak ?',optpar(5))
          end if
          call str2i (optpar(5),ndum,ierr)
          if (ierr .ne. 0) goto 10
          pfirst = ndum
c
c ... PICK LIMITS
c
        else if (optpar(2)(1:2) .eq. 'LI') then
c
          do i=1,3
            do j=1,2
              idum = 2*(i-1)+j+2
              if (nopt .lt. idum) then
                line = ' '//gkuplo(j)//' bound for '//
     +                 gkxyz(i)//' ?'
                write (optpar(idum),*) pklim(j,i)
                call remspa (optpar(idum))
                call textin (line,optpar(idum))
              end if
              call str2i (optpar(idum),jdum,ierr)
              if (ierr .eq. 0) then
                pklim (j,i) = jdum
              end if
            end do
          end do
c
c ... PICK PEAKS/INTEGRATE/HIGH
c
        else if (optpar(2)(1:2) .eq. 'PE' .or.
     +           optpar(2)(1:2) .eq. 'IN' .or.
     +           optpar(2)(1:2) .eq. 'HI') then
c
        if (nopt .lt. 3) call textin (' Which map ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = pkfile
          call textin (' Output file ?',optpar(4))
        end if
        if (length(optpar(4)) .gt. 0) pkfile = optpar(4)
c
        if (nopt .lt. 5) then
          optpar (5) = 'PDB'
          call textin (' Format (Pdb|Amore) ?',optpar(5))
        end if
        call upcase (optpar(5))
        if (optpar(5)(1:1) .ne. 'A') optpar(5)='PDB'
c
        if (optpar(2)(1:2) .eq. 'PE') then
          call pickem (jmap,0,optpar(5),
     +      maxmap,maxpnt,map,buffer)
        else if (optpar(2)(1:2) .eq. 'IN') then
          call pickem (jmap,1,optpar(5),
     +      maxmap,maxpnt,map,buffer)
        else if (optpar(2)(1:2) .eq. 'HI') then
          call pickem (jmap,2,optpar(5),
     +      maxmap,maxpnt,map,buffer)
        end if
c
c ... invalid sub-command
c
        else
c
          call errcon ('Invalid PICK option')
          call textut (' Selected option :',optpar(2))
c
        end if         
c
c ... DUPLICATE
c
      else if (optpar(1)(1:2) .eq. 'DU') then
c
        if (nopt .lt. 2) call textin (' New map ?',optpar(2))
        call allocm (optpar(2),imap,ierr,maxmap,maxpnt)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) call textin (' Old map ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        call copymp (imap,jmap,maxmap,maxpnt,map,.true.)
c
c ... ROUGHNESS
c
      else if (optpar(1)(1:2) .eq. 'RO') then
c
        if (nopt .lt. 2) call textin (' New map ?',optpar(2))
        call allocm (optpar(2),imap,ierr,maxmap,maxpnt)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) call textin (' Old map ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          write (optpar(4),*) npos
          call textin (' Nr of positive peaks to mask out ?',
     +      optpar(4))
        end if
        call str2i (optpar(4),idum,ierr)
        if (ierr .ne. 0) goto 10
        npos = max (0, min (idum, 1000))
c
        if (nopt .lt. 5) then
          write (optpar(5),*) nneg
          call textin (' Nr of negative peaks to mask out ?',
     +      optpar(5))
        end if
        call str2i (optpar(5),idum,ierr)
        if (ierr .ne. 0) goto 10
        nneg = max (0, min (idum, 1000))
c
        if (nopt .lt. 6) then
          write (optpar(6),*) ibox
          call textin (' Half box size ?',optpar(6))
        end if
        call str2i (optpar(6),idum,ierr)
        if (ierr .ne. 0) goto 10
        ibox = max (1, min (idum, 9))
c
        if (nopt .lt. 7) then
          optpar(7) = 'dens_hist.plt'
          call textin (' O2D plot file ?',optpar(7))
          if (length(optpar(7)) .lt. 1) goto 10
        end if
c
        if (nopt .lt. 8) then
          optpar(8) = '100'
          call textin (' Number of bins ?',optpar(8))
        end if
        call upcase (optpar(8))
        call str2i (optpar(8),idum,ierr)
        if (ierr .ne. 0) goto 10
        idum = max (10, idum)
c
        call copymp (imap,jmap,maxmap,maxpnt,map,.false.)
c
        call rough (imap,jmap,shadow,maxmap,maxpnt,map,
     +    npos,nneg,ibox,optpar(7),idum)
c
c ... QINVERT
c
      else if (optpar(1)(1:2) .eq. 'QI') then
c
        if (nopt .lt. 2) call textin (' New map ?',optpar(2))
        call allocm (optpar(2),imap,ierr,maxmap,maxpnt)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) call textin (' Old map ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        call invmap (imap,jmap,maxmap,maxpnt,map)
c
c ... ERASE/INTEGRATE/VALUES_AT
c
      else if (optpar(1)(1:2) .eq. 'ER' .or.
     +         optpar(1)(1:2) .eq. 'IN' .or.
     +         optpar(1)(1:2) .eq. 'VA') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        do i=1,3
          do j=1,2
            idum = 2*(i-1)+j+2
            if (nopt .lt. idum) then
              line = ' '//gkuplo(j)//' bound for '//
     +               gkxyz(i)//' ?'
              write (optpar(idum),*) gklim(j,i)
              call remspa (optpar(idum))
              call textin (line,optpar(idum))
            end if
            call str2i (optpar(idum),jdum,ierr)
            if (ierr .eq. 0) then
              gklim (j,i) = jdum
            end if
          end do
        end do
c
        do i=1,3
          gklim (1,i) = max (origin(i,imap),
     +         min(origin(i,imap)+extent(i,imap)-1,gklim(1,i)))
          gklim (2,i) = max (origin(i,imap),
     +         min(origin(i,imap)+extent(i,imap)-1,gklim(2,i)))
          call ilohi (gklim(1,i),gklim(2,i))
        end do
        call ivalut (' Limits :',6,gklim(1,1))
c
        if (optpar(1)(1:2) .eq. 'ER') then
c
c ... erase a block
c
          call textut (' Erase :',name(imap))
          call eramap (map(1,imap),origin(1,imap),extent(1,imap),
     +      extent(2,imap),extent(3,imap),gklim)
c
          call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +        edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
          edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
          change (imap) = .true.
c
c          call dynran (maprod(imap),maplus(imap),pdum,rdum)
c          call proplu (pdum,rdum,maprod(imap),maplus(imap))
c
           call proplu (edstat(1,imap),edstat(2,imap),
     +       maprod(imap),maplus(imap))
c
        else if (optpar(1)(1:2) .eq. 'IN') then
c
c ... integrate a block
c
          call textut (' Integrate :',name(imap))
          call intmap (map(1,imap),origin(1,imap),extent(1,imap),
     +      extent(2,imap),extent(3,imap),gklim)
c
        else if (optpar(1)(1:2) .eq. 'VA') then
c
c ... list values in a block
c
          call textut (' List values :',name(imap))
          call values (map(1,imap),origin(1,imap),extent(1,imap),
     +      extent(2,imap),extent(3,imap),gklim)
c
        end if
c
c ... UVW
c
      else if (optpar(1)(1:2) .eq. 'UV') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        do i=1,3
          if (nopt .lt. (i+2)) then
            write (optpar(i+2),'(i1)') uvw(i,imap)
            call remspa (optpar(i+2))
            write (line,'(a,i1,a)') ' Value for uvw ',i,' ?'
            call textin (line,optpar(i+2))
          end if
          call str2i(optpar(i+2),ndum,ierr)
          if (ierr .ne. 0) goto 10
          uvw(i,imap) = max (1, ndum)
        end do
c
        call ivalut (' UVW :',3,uvw(1,imap))
c
c ... check
c
        do i=1,3
          ndum = iindex (i,0,3,uvw(1,imap))
          if (ndum .le. 0) then
            call errcon (
     +        'Invalid UVW; not a permutation of 1 2 3')
            uvw (1,imap) = 2
            uvw (2,imap) = 1
            uvw (3,imap) = 3
            call ivalut (' Reset uvw :',3,uvw(1,imap))
            goto 10
          end if
        end do
c
c ... SPACEGROUP
c
      else if (optpar(1)(1:2) .eq. 'SP') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          write (optpar(3),'(i3)') spaceg(imap)
          call textin (' Spacegroup ?',optpar(3))
        end if
        call str2i (optpar(3),ndum,ierr)
        if (ierr .ne. 0) goto 10
c
        ndum = max (1,ndum)
        spaceg (imap) = ndum
        call ivalut (' Spacegroup :',1,spaceg(imap))
        call prompt (' NOTE: THIS DOES NOT CHANGE OR EVEN ADD')
        call prompt ('       SYMMETRY OPERATORS TO CCP4 MAPS !')
c
c ... LABEL
c
      else if (optpar(1)(1:2) .eq. 'LA') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        if (nopt .lt. 3) call textin (' Label ?',optpar(3))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        commnt (imap) = optpar(3)
c
c ... TRANSLATE
c
      else if (optpar(1)(1:2) .eq. 'TR') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        call jvalut (' Old origin :',3,origin(1,imap))
        call jvalut (' Grid       :',3,grid(1,imap))
c
        if (nopt .lt. 3) then
          optpar(3) = '0'
          call textin (
     +      ' Nr of UNIT CELLs translation in X ?',optpar(3))
        end if
        call str2i (optpar(3),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (1,imap) = origin(1,imap) + itra*grid(1,imap)
c
        if (nopt .lt. 4) then
          optpar(4) = '0'
          call textin (
     +      ' Nr of UNIT CELLs translation in Y ?',optpar(4))
        end if
        call str2i (optpar(4),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (2,imap) = origin(2,imap) + itra*grid(2,imap)
c
        if (nopt .lt. 5) then
          optpar(5) = '0'
          call textin (
     +      ' Nr of UNIT CELLs translation in Z ?',optpar(5))
        end if
        call str2i (optpar(5),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (3,imap) = origin(3,imap) + itra*grid(3,imap)
c
        call jvalut (' New origin :',3,origin(1,imap))
c
c ... GTRANSLATE
c
      else if (optpar(1)(1:2) .eq. 'GT') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        write (*,'(/10(1x,a/))')
     +    'WARNING !!! This option will change the origin of your',
     +    'map by any number of grid points that you specify.',
     +    'This is bound to upset the compatibility between this',
     +    'map and your model/data. Only use this option if you',
     +    'understand this and know what you are doing ...'
c
        call jvalut (' Old origin :',3,origin(1,imap))
        call jvalut (' Grid       :',3,grid(1,imap))
c
        if (nopt .lt. 3) then
          optpar(3) = '0'
          call textin (
     +      ' Nr of GRID POINTs translation in X ?',optpar(3))
        end if
        call str2i (optpar(3),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (1,imap) = origin(1,imap) + itra
c
        if (nopt .lt. 4) then
          optpar(4) = '0'
          call textin (
     +      ' Nr of GRID POINTs translation in Y ?',optpar(4))
        end if
        call str2i (optpar(4),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (2,imap) = origin(2,imap) + itra
c
        if (nopt .lt. 5) then
          optpar(5) = '0'
          call textin (
     +      ' Nr of GRID POINTs translation in Z ?',optpar(5))
        end if
        call str2i (optpar(5),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (3,imap) = origin(3,imap) + itra
c
        call jvalut (' New origin :',3,origin(1,imap))
c
c ... CELL
c
      else if (optpar(1)(1:2) .eq. 'CE') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        call fvalut (' CURRENT Cell :',6,cell(1,imap))
        do i=1,6
          if (nopt .lt. (i+2)) then
            write (optpar(i+2),'(f15.3)') cell(i,imap)
            call remspa (optpar(i+2))
            write (line,'(a,i1,a)') ' Value for cell ',i,' ?'
            call textin (line,optpar(i+2))
          end if
          call str2r(optpar(i+2),rdum,ierr)
          if (ierr .ne. 0) goto 10
          cell (i,imap) = rdum
        end do
        call fvalut (' CHANGED Cell :',6,cell(1,imap))
c
c ... READ
c
      else if (optpar(1)(1:2) .eq. 'RE') then
c
        if (nopt .lt. 2) call textin (' New map ?',optpar(2))
        call allocm (optpar(2),imap,ierr,maxmap,maxpnt)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          optpar (3) = file(imap)
          call textin (' File name ?',optpar(3))
        end if
        file (imap) = optpar(3)
        if (length(file(imap)) .lt. 1) then
          call errcon ('No file name provided')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = 'CCP4'
          call textin (' Format ?',optpar(4))
        end if
        call upcase (optpar(4))
        if (optpar(4)(1:1) .eq. '?') then
          call asciut (' Supported formats :',ninfmt,inpfmt)
          goto 10
        end if
c
        do i=1,ninfmt
          if (optpar(4)(1:10) .eq. inpfmt(i)) then
            ifmt=i
            goto 6501
          end if
        end do
        call errcon ('Invalid map format: '//optpar(4))
        goto 10
c
 6501   continue
c
        incore (imap) = .false.
        spaceg (imap) = 1
        uvw (1,imap) = 2
        uvw (2,imap) = 1
        uvw (3,imap) = 3
c
c ... XPLOR = CNS = X-PLOR
c
        if (ifmt .eq. 17 .or. ifmt .eq. 18) ifmt = 5
c
c ... EZD = NEWEZD
c
        if (ifmt .eq. 19) ifmt = 8
c
c ... MAPPAGE/DSN6
c
        if (ifmt .eq. 10 .or. ifmt .eq. 11 .or. ifmt .eq. 21) then
          call undsn6 (imap,ierr,maxmap,maxpnt,map,.false.)
c
c ... BRIX
c
        else if (ifmt .eq. 16) then
          call undsn6 (imap,ierr,maxmap,maxpnt,map,.true.)
        else if (ifmt .eq. 12) then
          call xpl3dm (imap,ierr,maxmap,maxpnt,map)
        else if (ifmt .eq. 20) then
          call rdem08 (imap,ierr,maxmap,maxpnt,map)
        else if (ifmt .eq. 22) then
          call rdempi (imap,ierr,maxmap,maxpnt,map)
        else if (ifmt .eq. 14 .or. ifmt .eq. 15) then
          call rfurey (imap,ifmt,ierr,maxmap,maxpnt,map)
        else
          call maprd (imap,ifmt,ierr,
     +                maxmap,maxpnt,map,buffer)
        end if
c
        if (ierr .ne. 0) goto 10
c
        call prompt (' Map read into memory - calculating statistics')
        npnt (imap) = extent(1,imap)*extent(2,imap)*extent(3,imap)
        call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +    edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
        edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
        call rvalut (' Sum of density in map :',1,pdum)
c
        call proplu (edstat(1,imap),edstat(2,imap),maprod(imap),
     +    maplus(imap))
c
        change (imap) = .false.
        incore (imap) = .true.
        commnt (imap) = 'Read from '//file(imap)
c
        goto 10
c
c ... FRAME
c
      else if (optpar(1)(1:2) .eq. 'FR') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = frtype
          call textin (
     +      ' Input Indices/Fractional/Cartesian (I/F/C) ?',optpar(3))
        end if
c
        frtype = optpar(3)(1:1)
        call upcase (frtype)
        if (index ('I|F|C',frtype) .le. 0) then
          call errcon (' Invalid coordinate type')
          goto 10
        end if
c
        call prompt (' Note: indices assumed to start at 0 !')
c
        do i=1,3
          idum = i+3
          if (nopt .lt. idum) then
            line = ' Value for '//gkxyz(i)//' ?'
            write (optpar(idum),*) frlim(i)
            call remspa (optpar(idum))
            call textin (line,optpar(idum))
          end if
          call str2r (optpar(idum),xdum,ierr)
          if (ierr .eq. 0) then
            frlim (i) = xdum
          end if
        end do
c
        call orthog (cell(1,imap), a, 0)
        call orthog (cell(1,imap), b, 1)
c
        if (frtype .eq. 'I') then
c
c ... Indices
c
          do i=1,3
            frind (i) = nint(frlim (i))
            fract (i) = float(frind(i))/float(grid(i,imap))
          end do
          call mulmtx (a,fract,frcar, 3, 3, 1)
          call ivalut (' Indices    :',3,frind)
          call fvalut (' Fractional :',3,fract)
          call fvalut (' Cartesian  :',3,frcar)
c
        else if (frtype .eq. 'F') then
c
c ... Fractional
c
          do i=1,3
            fract (i) = frlim (i)
            frind (i) = nint (grid(i,imap)*fract(i))
          end do
          call mulmtx (a,fract,frcar, 3, 3, 1)
          call fvalut (' Fractional :',3,fract)
          call ivalut (' Indices    :',3,frind)
          call fvalut (' Cartesian  :',3,frcar)
c
        else if (frtype .eq. 'C') then
c
c ... Cartesian
c
          do i=1,3
            frcar (i) = frlim (i)
          end do
          call mulmtx (b,frcar,fract, 3, 3, 1)
          do i=1,3
            frind (i) = nint (grid(i,imap)*fract(i))
          end do
          call fvalut (' Cartesian  :',3,frcar)
          call fvalut (' Fractional :',3,fract)
          call ivalut (' Indices    :',3,frind)
c
        end if
c
        if (frind(1) .lt. origin(1,imap) .or.
     +      frind(1) .gt. (origin(1,imap)+extent(1,imap)-1) .or.
     +      frind(2) .lt. origin(2,imap) .or.
     +      frind(2) .gt. (origin(2,imap)+extent(2,imap)-1) .or.
     +      frind(3) .lt. origin(3,imap) .or.
     +      frind(3) .gt. (origin(3,imap)+extent(3,imap)-1)) then
          call prompt (' Point lies outside the map boundaries')
        else
          call prompt (' Point lies inside the map')
        end if
c
c ... LIST
c
      else if (optpar(1)(1:2) .eq. 'LI') then
c
        nmap = 0
        do i=1,maxmap
          if (incore(i)) nmap = nmap + 1
        end do
        call ivalut (' Nr of maps in memory :',1,nmap)
c
        if (nopt .lt. 2)  optpar(2) = '*'
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nmap .gt. 0) then
c
          do i=1,maxmap
            if (select(i)) then
              if (incore(i)) then
c
                call voxvol (cell(1,i),grid(1,i),cellv,voxv)
                gridv = voxv * float(npnt(i))
                perc = 100.0*gridv/cellv
c
                call dynran (maprod(i),maplus(i),xdum,rdum)
c
                call orthog (cell(1,i), a, 0)
c
c ... calc origin and top in fractionals
c
                do j=1,3
                  flo (j) = float(origin(j,i)) * 
     +             (cell(j,i)/float(grid(j,i))) / cell(j,i)
                  fhi (j) = float(extent(j,i)+origin(j,i)-1) *
     +             (cell(j,i)/float(grid(j,i))) / cell(j,i)
                end do
c
                call mulmtx (a, flo, xlo, 3, 3, 1)
                call mulmtx (a, fhi, xhi, 3, 3, 1)
c
                write (*,6020) i,name(i)(1:leng1(name(i))),
     +            file(i)(1:leng1(file(i))),(grid(j,i),j=1,3),
     +            (origin(j,i),j=1,3),(extent(j,i),j=1,3),
     +            (cell(j,i),j=1,6),(uvw(j,i),j=1,3),
     +            spaceg(i),npnt(i),(edstat(j,i),j=1,5),
     +            maprod(i),maplus(i),xdum,rdum
                write (*,6023) cellv,voxv,gridv,perc
                write (*,6027) flo,xlo,fhi,xhi
                write (*,6025)
     +            (cell(j,i)/float(grid(j,i)),j=1,3),
     +            (origin(j,i)+extent(j,i)-1,j=1,3),
     +            change(i),commnt(i)(1:leng1(commnt(i)))
              else
                write (*,6030) i
              end if
            end if
          end do
c
        end if
c
        goto 10
c
 6020 format (/' Map ',i3,'     = ',a/
     +  ' File        = ',a/
     +  ' Grid        = ',3i10/
     +  ' Origin      = ',3i10/
     +  ' Extent      = ',3i10/
     +  ' Cell        = ',6f10.3/
     +  ' UVW         = ',3i10/
     +  ' Spcgrp      = ',i10,  5x,' Nr of points = ',i10/1p,
     +  ' Density min = ',e10.3,5x,' Max          = ',e10.3/
     +  ' Average     = ',e10.3,5x,' Sigma        = ',e10.3/
     +  ' Variance    = ',e10.3/
     +  ' Mappage Prod= ',e10.3,5x,' Plus         = ',i10/
     +  ' Dyn. range  = ',e10.3,5x,' To           = ',e10.3)
c
 6023 format (1p,
     +  ' Cell vol.   = ',e10.3,5x,' Voxel vol.   = ',e10.3/
     +  ' Grid vol.   = ',e10.3,5x,' %Cell vol.   = ',0p,f10.2)
c
 6027 format (
     +  ' Origin frac = ',3f10.5/
     +  ' Origin Cart = ',3f10.3/
     +  ' Top frac    = ',3f10.5/
     +  ' Top Cart    = ',3f10.3)
c
 6025 format (
     +  ' Spacing     = ',3f10.3/
     +  ' Top         = ',3i10/
     +  ' Changes     = ',l1/
     +  ' Label       = ',a)
c
 6030 format (/' Map ',i2,' available')
c
c ... MOMENTS
c
      else if (optpar(1)(1:2) .eq. 'MO') then
c
        if (nopt .lt. 2)  call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Moments :',name(i))
            call mpstat (npnt(i),map(1,i))
          end if
        end do
c
c ... HISTOGRAM
c
      else if (optpar(1)(1:2) .eq. 'HI') then
c
        if (nopt .lt. 2)  call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
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
ccc        call qsortg (rhis,numhis)
        call hsortr (numhis,rhis)
        call rvalut (' Histogram limits :',numhis,rhis)
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Histogram :',name(i))
            call histo (npnt(i),map(1,i),numhis,rhis,nhis)
          end if
        end do
c
c ... DELETE
c
      else if (optpar(1)(1:2) .eq. 'DE') then
c
        if (nopt .lt. 2)  call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            if (change(i)) then
              write (*,6040) name(i)(1:leng1(name(i)))
              if (linter) then
                reply = 'N'
                call textin (
     +            ' Really DELETE this map (Y/N) ?',reply)
                call upcase (reply)
                if (reply .ne. 'Y') goto 1210
              end if
            end if
c
            incore (i) = .false.
            change (i) = .false.
            commnt (i) = 'No comment'
            call textut (' Deleted :',name(i))
            name (i) = '!$%#?@'
c
 1210       continue
          end if
        end do
c
c ... ODL
c
      else if (optpar(1)(1:2) .eq. 'OD') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = 'map.odl'
          call textin (' ODL file ?',optpar(3))
        end if
c
        if (length(optpar(3)) .lt. 1) then
          call errcon ('Invalid ODL file name')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = 'box'
          call textin (' Draw Cell or Box (C/B) ?',optpar(4))
        end if
        call upcase (optpar(4))
c
        if (nopt .lt. 5) then
          optpar (5) = 'solid'
          call textin (' Draw in Line or Solid mode (L/S) ?',
     +      optpar(5))
        end if
        call upcase (optpar(5))
c
        call odlmap (imap,optpar(3),optpar(4),optpar(5))
c
c ... PROD_PLUS
c
      else if (optpar(1)(1:2) .eq. 'PR') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f3.1)') 1.0
          call textin (' PROD ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(i1)') 1
          call textin (' PLUS ?',optpar(4))
        end if
        call str2i (optpar(4),ndum,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Prod_plus :',name(i))
            maprod (i) = xdum
            maplus (i) = ndum
            call dynran (maprod(i),maplus(i),pdum,rdum)
            call proplu (pdum,rdum,maprod(i),maplus(i))
          end if
        end do
c
c ... RANGE
c
      else if (optpar(1)(1:2) .eq. 'RA') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f6.1)') -100.0
          call textin (' Lower bound of dynamic range ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f6.1)') 100.0
          call textin (' Upper bound of dynamic range ?',optpar(4))
        end if
        call str2r (optpar(4),ydum,ierr)
        if (ierr .ne. 0) goto 10
c
        call rlohi (xdum,ydum)
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Range :',name(i))
            call proplu (xdum,ydum,maprod(i),maplus(i))
            call dynran (maprod(i),maplus(i),pdum,rdum)
          end if
        end do
c
c ... DRANGE
c
      else if (optpar(1)(1:2) .eq. 'DR') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f6.1)') -10.0
          call textin (' Lower bound of dynamic range ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f6.1)') 10.0
          call textin (' Upper bound of dynamic range ?',optpar(4))
        end if
        call str2r (optpar(4),ydum,ierr)
        if (ierr .ne. 0) goto 10
c
        call rlohi (xdum,ydum)
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Cap dynamic range for map :',name(i))
            pdum = edstat (1,i)
            rdum = edstat (2,i)
            call rvalut (' ... Lower before :',1,pdum)
            call rvalut (' ... Upper before :',1,rdum)
            if (xdum .gt. pdum) pdum = xdum
            if (ydum .lt. rdum) rdum = ydum
            call proplu (pdum,rdum,maprod(i),maplus(i))
            call dynran (maprod(i),maplus(i),pdum,rdum)
            call rvalut (' ... Lower after  :',1,pdum)
            call rvalut (' ... Upper after  :',1,rdum)
          end if
        end do
c
c ... SCALE
c
      else if (optpar(1)(1:2) .eq. 'SC') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f6.1)') -100.0
          call textin (' Lower bound ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f6.1)') 100.0
          call textin (' Upper bound ?',optpar(4))
        end if
        call str2r (optpar(4),ydum,ierr)
        if (ierr .ne. 0) goto 10
c
        call rlohi (xdum,ydum)
        if ( (ydum-xdum).lt.0.001 ) then
          call errcon ('Range too narrow')
          goto 10
        end if
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Scale :',name(i))
            xpr = (ydum-xdum)/(edstat(2,i)-edstat(1,i))
            xpl = xdum - edstat(1,i)*xpr
            call rvalut (' Multiply map by :',1,xpr)
            call rvalut (' ... and add     :',1,xpl)
            call mapmat (map(1,i),npnt(i),xpr,xpl)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),pdum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            call proplu ( (xpr*pdum+xpl), (xpr*rdum+xpl),
c     +        maprod(i),maplus(i))
c
             call proplu (edstat(1,i),edstat(2,i),
     +         maprod(i),maplus(i))
c
          end if
        end do
c
c ... ZERO
c
      else if (optpar(1)(1:2) .eq. 'ZE') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f6.1)') -100.0
          call textin (' Lower bound ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f6.1)') 100.0
          call textin (' Upper bound ?',optpar(4))
        end if
        call str2r (optpar(4),ydum,ierr)
        if (ierr .ne. 0) goto 10
c
        call rlohi (xdum,ydum)
        if ( (ydum-xdum).lt.0.001 ) then
          call errcon ('Range too narrow')
          goto 10
        end if
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Zero :',name(i))
            call zeromp (map(1,i),npnt(i),xdum,ydum)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),pdum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            pdum = max (pdum, edstat(1,i))
c            rdum = min (rdum, edstat(2,i))
c            call proplu (pdum,rdum,maprod(i),maplus(i))
c
             call proplu (edstat(1,i),edstat(2,i),
     +         maprod(i),maplus(i))
c
          end if
        end do
c
c ... CHANGE_VALUE
c
      else if (optpar(1)(1:2) .eq. 'CH') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(1pe12.4)') -10000.0
          call textin (' Lower bound ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(1pe12.4)') xdum+1.0
          call textin (' Upper bound ?',optpar(4))
        end if
        call str2r (optpar(4),ydum,ierr)
        if (ierr .ne. 0) goto 10
c
        call rlohi (xdum,ydum)
c
        if (nopt .lt. 5) then
          write (optpar(5),'(f6.1)') 0.0
          call textin (' New value ?',optpar(5))
        end if
        call str2r (optpar(5),zdum,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Change :',name(i))
            call chaval (map(1,i),npnt(i),xdum,ydum,zdum)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),pdum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            pdum = max (pdum, edstat(1,i))
c            rdum = min (rdum, edstat(2,i))
c            call proplu (pdum,rdum,maprod(i),maplus(i))
c
             call proplu (edstat(1,i),edstat(2,i),
     +         maprod(i),maplus(i))
c
          end if
        end do
c
c ... SET_VALUE
c
      else if (optpar(1)(1:2) .eq. 'SE') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(1pe12.4)') set1st
          call textin (' Start value ?',optpar(3))
        end if
        call str2r (optpar(3),set1st,ierr)
        if (ierr .ne. 0) goto 10
c
        setdx = 0.0
        if (nopt .ge. 4) then
          call str2r (optpar(4),setdx,ierr)
          if (ierr .ne. 0) goto 10
        end if
c
        setdy = 0.0
        if (nopt .ge. 5) then
          call str2r (optpar(5),setdy,ierr)
          if (ierr .ne. 0) goto 10
        end if
c
        setdz = 0.0
        if (nopt .ge. 6) then
          call str2r (optpar(6),setdz,ierr)
          if (ierr .ne. 0) goto 10
        end if
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Set values :',name(i))
            call rvalut (' Start value  :',1,set1st)
            call rvalut (' DX increment :',1,setdx)
            call rvalut (' DY increment :',1,setdy)
            call rvalut (' DZ increment :',1,setdz)
            call setval (map(1,i),extent(1,i),extent(2,i),
     +                   extent(3,i),set1st,setdx,setdy,setdz)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),pdum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            pdum = max (pdum, edstat(1,i))
c            rdum = min (rdum, edstat(2,i))
c            call proplu (pdum,rdum,maprod(i),maplus(i))
c
             call proplu (edstat(1,i),edstat(2,i),
     +         maprod(i),maplus(i))
c
          end if
        end do
c
c ... WHERE_VALUE
c
      else if (optpar(1)(1:2) .eq. 'WH') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f10.3)') 1.0
          call textin (' Density value ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 4) then
          write (optpar(4),'(f10.6)') 0.00001
          call textin (' Tolerance ?',optpar(4))
        end if
        call str2r (optpar(4),ydum,ierr)
        if (ierr .ne. 0) goto 10
c
        call rvalut (' Density value :',1,xdum)
        call rvalut (' Tolerance     :',1,ydum)
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Checking map :',name(i))
            call wherev (map(1,i),origin(1,i),extent(1,i),
     +                   extent(2,i),extent(3,i),xdum,ydum)
          end if
        end do
c
c ... MULTIPLY
c
      else if (optpar(1)(1:2) .eq. 'MU') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f3.1)') 1.0
          call textin (' Multiply by ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Multiply :',name(i))
            call mapmat (map(1,i),npnt(i),xdum,0.0)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),ydum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            call proplu (xdum*pdum,xdum*rdum,maprod(i),maplus(i))
c
            call proplu (edstat(1,i),edstat(2,i),
     +        maprod(i),maplus(i))
c
          end if
        end do
c
c ... NORMALISE
c
      else if (optpar(1)(1:2) .eq. 'NO') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Normalise :',name(i))
            xdum = edstat(3,i)
            call rvalut (' Subtract average :',1,xdum)
            call mapmat (map(1,i),npnt(i),1.0,-1.0*xdum)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),ydum)
            xdum = edstat(4,i)
            call rvalut (' Divide by sigma  :',1,xdum)
            xdum = 1.0 / xdum
            call mapmat (map(1,i),npnt(i),xdum,0.0)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),ydum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            call proplu (xdum*pdum,xdum*rdum,maprod(i),maplus(i))
c
             call proplu (edstat(1,i),edstat(2,i),
     +         maprod(i),maplus(i))
c
          end if
        end do
c
c ... DIVIDE
c
      else if (optpar(1)(1:2) .eq. 'DI') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f3.1)') 1.0
          call textin (' Divide by ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
        if (xdum .eq. 0.0) then
          call errcon ('Attempted division by zero')
          goto 10
        end if
        xdum = 1.0 / xdum
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Divide :',name(i))
            call mapmat (map(1,i),npnt(i),xdum,0.0)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),ydum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            call proplu (xdum*pdum,xdum*rdum,maprod(i),maplus(i))
c
            call proplu (edstat(1,i),edstat(2,i),
     +        maprod(i),maplus(i))
c
          end if
        end do
c
c ... PLUS
c
      else if (optpar(1)(1:2) .eq. 'PL') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        call selecm (optpar(2),maxmap,maxpnt,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          write (optpar(3),'(f3.1)') 0.0
          call textin (' Add ?',optpar(3))
        end if
        call str2r (optpar(3),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmap
          if (select(i)) then
            call textut (' Plus :',name(i))
            call mapmat (map(1,i),npnt(i),1.0,xdum)
            call xstats (map(1,i),npnt(i),edstat(3,i),
     +        edstat(4,i),edstat(1,i),edstat(2,i),ydum)
            edstat(5,i) = edstat(4,i)*edstat(4,i)
            change (i) = .true.
c
c            call dynran (maprod(i),maplus(i),pdum,rdum)
c            call proplu (xdum+pdum,xdum+rdum,maprod(i),maplus(i))
c
            call proplu (edstat(1,i),edstat(2,i),
     +        maprod(i),maplus(i))
c
          end if
        end do
c
c ... WRITE
c
      else if (optpar(1)(1:2) .eq. 'WR') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        if (nopt .lt. 3) call textin (' Which file ?',optpar(3))
c
        if (nopt .lt. 4) then
          optpar (4) = 'CCP4'
          call textin (' Format ?',optpar(4))
        end if
        call upcase (optpar(4))
        if (optpar(4)(1:1) .eq. '?') then
          call asciut (' Supported formats :',nutfmt,outfmt)
          goto 10
        end if
c
        do i=1,nutfmt
          if (optpar(4)(1:10) .eq. outfmt(i)) then
            ifmt=i
            goto 6401
          end if
        end do
        call errcon ('Invalid map format: '//optpar(4))
        goto 10
c
 6401   continue
c
c ... XPLOR = CNS = X-PLOR
c
        if (ifmt .eq. 9 .or. ifmt .eq. 10) ifmt = 6
c
c        if (ifmt .eq. 7) then
c          call prompt (
c     +      ' Use the MAppage command to write DSN6 maps !')
c          goto 10
c        end if
c        if (ifmt .eq. 8) then
c          call prompt (
c     +      ' Use the BRix command to write BRIX maps !')
c          goto 10
c        end if
c
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        file (imap) = optpar(3)
        call upcase (optpar(4))
c
c ... if output is a mask, ask for cut-off
c
        if (ifmt .eq. 3 .or. ifmt .eq. 5) then
          if (nopt .lt. 5) then
            write (optpar(5),'(f6.4)') 0.001
            call textin (' Cut-off between non-mask and mask ?',
     +        optpar(5))
          end if
          call str2r (optpar(5),xdum,ierr)
          if (ierr .ne. 0) goto 10
          if (nopt .lt. 6) optpar(6) = '>'
c            call textin (
c     +        ' Set mask point >, >=, <, or <= cut-off ?',
c     +        optpar(6))
          call remspa (optpar(6))
          call upcase (optpar(6))
        end if
c
        if (ifmt .eq. 7 .or. ifmt .eq. 13) then
          call prompt (' Doing MAPPAGE for you')
          call mappag (imap,file(imap),'MAPPAGE',
     +                 maxmap,maxpnt,map,shadow)
          ierr = 0
        else if (ifmt .eq. 8) then
          call prompt (' Doing BRIX for you')
          call mappag (imap,file(imap),'BRIX',
     +                 maxmap,maxpnt,map,shadow)
          ierr = 0
        else if (ifmt .eq. 14) then
          call prompt (' Doing TURBO-MAPPAGE for you')
          call mappag (imap,file(imap),'TURBO',
     +                 maxmap,maxpnt,map,shadow)
          ierr = 0
        else
          call writem (imap,ifmt,xdum,optpar(6),ierr,
     +                 maxmap,maxpnt,map,buffer,shadow)
        end if
c
        if (ierr .eq. 0) change(imap) = .false.
c
        goto 10
c
c ... ADD
c
      else if (optpar(1)(1:2) .eq. 'AD') then
c
        if (nopt .lt. 2) call textin (' Which map1 ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which map2 ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (.not. simap(imap,jmap)) goto 10
c
        call textut (' Add to Map :',name(imap))
        call textut (' Map to Add :',name(jmap))
c
        call addmap (imap,jmap,ierr,maxmap,maxpnt,map)
        if (ierr .ne. 0) goto 10
c
        change (imap) = .true.
        call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +        edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
        edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c
c        call dynran (maprod(imap),maplus(imap),pdum,rdum)
c        pdum = max (pdum, edstat(1,imap))
c        rdum = min (rdum, edstat(2,imap))
c        call proplu (pdum,rdum,maprod(imap),maplus(imap))
c
        call proplu (edstat(1,imap),edstat(2,imap),
     +    maprod(imap),maplus(imap))
c
c ... OPERATE
c
      else if (optpar(1)(1:2) .eq. 'OP') then
c
        if (nopt .lt. 2) call textin (' Which map1 ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = '+'
          call textin (' Operator (+|-|*) ?',optpar(3))
        end if
c
        call upcase (optpar(3))
        call remspa (optpar(3))
        if (index ('+-*',optpar(3)(1:1)) .le. 0) then
          call errcon ('Invalid operator')
          goto 10
        end if
c
        if (nopt .lt. 4) call textin (' Which map2 ?',optpar(4))
        jmap = whichm(optpar(4),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call prompt (' WARNING - Map1 and Map2 are identical !')
        end if
c
        if (.not. simap(imap,jmap)) goto 10
c
        call textut (' Store Map :',name(imap))
        call textut (' Other Map :',name(jmap))
        call textut (' Operator  :',optpar(3))
c
        call opemap (imap,jmap,optpar(3),ierr,maxmap,maxpnt,map)
        if (ierr .ne. 0) goto 10
c
        change (imap) = .true.
        call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +        edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
        edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c
c        call dynran (maprod(imap),maplus(imap),pdum,rdum)
c        pdum = max (pdum, edstat(1,imap))
c        rdum = min (rdum, edstat(2,imap))
c        call proplu (pdum,rdum,maprod(imap),maplus(imap))
c
        call proplu (edstat(1,imap),edstat(2,imap),
     +    maprod(imap),maplus(imap))
c
c ... MIN_MAX
c
      else if (optpar(1)(1:2) .eq. 'MI') then
c
        if (nopt .lt. 2) call textin (' Which map1 ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which map2 ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (.not. simap(imap,jmap)) goto 10
c
        call textut (' Minimum Map :',name(imap))
        call textut (' Maximum Map :',name(jmap))
c
        call minmax (imap,jmap,ierr,maxmap,maxpnt,map)
        if (ierr .ne. 0) goto 10
c
        change (imap) = .true.
        call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +        edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
        edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c
c        call dynran (maprod(imap),maplus(imap),pdum,rdum)
c        pdum = max (pdum, edstat(1,imap))
c        rdum = min (rdum, edstat(2,imap))
c        call proplu (pdum,rdum,maprod(imap),maplus(imap))
c
        call proplu (edstat(1,imap),edstat(2,imap),
     +    maprod(imap),maplus(imap))
c
        change (jmap) = .true.
        call xstats (map(1,jmap),npnt(jmap),edstat(3,jmap),
     +        edstat(4,jmap),edstat(1,jmap),edstat(2,jmap),pdum)
        edstat(5,jmap) = edstat(4,jmap)*edstat(4,jmap)
c
c        call dynran (maprod(jmap),maplus(jmap),pdum,rdum)
c        pdum = max (pdum, edstat(1,jmap))
c        rdum = min (rdum, edstat(2,jmap))
c        call proplu (pdum,rdum,maprod(jmap),maplus(jmap))
c
        call proplu (edstat(1,jmap),edstat(2,jmap),
     +    maprod(jmap),maplus(jmap))
c
c ... CORRELATE
c
      else if (optpar(1)(1:2) .eq. 'CO') then
c
        if (nopt .lt. 2) call textin (' Which map1 ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which map2 ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (.not. simap(imap,jmap)) goto 10
c
        call textut (' Correlate Map1 :',name(imap))
        call textut ('           Map2 :',name(jmap))
c
        call cormap (imap,jmap,ierr,maxmap,maxpnt,map)
        if (ierr .ne. 0) goto 10
c
        change (imap) = .true.
        call xstats (map(1,imap),npnt(imap),edstat(3,imap),
     +        edstat(4,imap),edstat(1,imap),edstat(2,imap),pdum)
        edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c
c        call dynran (maprod(imap),maplus(imap),pdum,rdum)
c        pdum = max (pdum, edstat(1,imap))
c        rdum = min (rdum, edstat(2,imap))
c        call proplu (pdum,rdum,maprod(imap),maplus(imap))
c
        call proplu (edstat(1,imap),edstat(2,imap),
     +    maprod(imap),maplus(imap))
c
        change (jmap) = .true.
        call xstats (map(1,jmap),npnt(jmap),edstat(3,jmap),
     +        edstat(4,jmap),edstat(1,jmap),edstat(2,jmap),pdum)
        edstat(5,jmap) = edstat(4,jmap)*edstat(4,jmap)
c
c        call dynran (maprod(jmap),maplus(jmap),pdum,rdum)
c        pdum = max (pdum, edstat(1,jmap))
c        rdum = min (rdum, edstat(2,jmap))
c        call proplu (pdum,rdum,maprod(jmap),maplus(jmap))
c
        call proplu (edstat(1,jmap),edstat(2,jmap),
     +    maprod(jmap),maplus(jmap))
c
c ... PASTE
c
      else if (optpar(1)(1:2) .eq. 'PA') then
c
c ... get map FROM which to copy
        if (nopt .lt. 2) call textin (
     +    ' Copy FROM which map1 ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
c ... get limits to copy
c
        do i=1,3
          do j=1,2
            idum = 2*(i-1)+j+2
            if (nopt .lt. idum) then
              line = ' '//gkuplo(j)//' bound for copy '//
     +               gkxyz(i)//' ?'
              write (optpar(idum),*) gklim(j,i)
              call remspa (optpar(idum))
              call textin (line,optpar(idum))
            end if
            call str2i (optpar(idum),jdum,ierr)
            if (ierr .eq. 0) then
              gklim (j,i) = jdum
            end if
          end do
        end do
c
        do i=1,3
          gklim (1,i) = max (origin(i,imap),
     +         min(origin(i,imap)+extent(i,imap)-1,gklim(1,i)))
          gklim (2,i) = max (origin(i,imap),
     +         min(origin(i,imap)+extent(i,imap)-1,gklim(2,i)))
          call ilohi (gklim(1,i),gklim(2,i))
        end do
        call ivalut (' Limits :',6,gklim(1,1))
c
        if (nopt .lt. 9) call textin (
     +    ' Copy *TO* which map2 ?',optpar(9))
        jmap = whichm(optpar(9),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (.not. simap(imap,jmap)) goto 10
c
        do i=1,3
          frind (i) = nint(frlim (i))
          idum = i+9
          if (nopt .lt. idum) then
            line = ' Lower bound for paste '//gkxyz(i)//' ?'
            write (optpar(idum),*) frind(i)
            call remspa (optpar(idum))
            call textin (line,optpar(idum))
          end if
          call str2i (optpar(idum),jdum,ierr)
          if (ierr .eq. 0) then
            frind (i) = jdum
          end if
        end do
        call ivalut (' Paste start :',3,frind)
c
        call pastem (
     +    map(1,imap),extent(1,imap),extent(2,imap),
     +    extent(3,imap),origin(1,imap),gklim(1,1),
     +    map(1,jmap),extent(1,jmap),extent(2,jmap),
     +    extent(3,jmap),origin(1,jmap),frind(1),ierr)
c
        if (ierr .eq. 0) then
          call xstats (map(1,jmap),npnt(jmap),edstat(3,jmap),
     +        edstat(4,jmap),edstat(1,jmap),edstat(2,jmap),pdum)
          edstat(5,jmap) = edstat(4,jmap)*edstat(4,jmap)
          change (jmap) = .true.
c
c          call dynran (maprod(jmap),maplus(jmap),pdum,rdum)
c          call proplu (pdum,rdum,maprod(jmap),maplus(jmap))
c
          call proplu (edstat(1,jmap),edstat(2,jmap),
     +      maprod(jmap),maplus(jmap))
c
        end if
c
c ... SIMILARITY
c
      else if (optpar(1)(1:2) .eq. 'SI') then
c
        if (nopt .lt. 2) call textin (' Which map1 ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which map2 ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (.not. simap(imap,jmap)) goto 10
c
        call textut (' Similarity map :',name(imap))
        call textut ('        and map :',name(jmap))
c
        call simmap (imap,jmap,maxmap,maxpnt,map)
c
c ... LC_map
c
      else if (optpar(1)(1:2) .eq. 'LC') then
c
        if (nopt .lt. 2) call textin (' New map ?',optpar(2))
        call allocm (optpar(2),kmap,ierr,maxmap,maxpnt)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) call textin (' Which map 1 ?',optpar(3))
        imap = whichm(optpar(3),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. kmap) then
          call errcon ('Map 1 and new map are identical')
          goto 10
        end if
c
        call copymp (kmap,imap,maxmap,maxpnt,map,.false.)
c
        if (nopt .lt. 4) call textin (' Which map 2 ?',optpar(4))
        jmap = whichm(optpar(4),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map 1 and map 2 are identical')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map 2 and new map are identical')
          goto 10
        end if
c
        if (.not. samgrd(imap,jmap)) goto 10
c
        if (nopt .lt. 5) then
          optpar (5) = 'cc'
ccc          call textin (' Map type: CC or R-factor (C|R) ?',optpar(5))
        end if
        call upcase (optpar(5))
        call remspa (optpar(5))
        if (optpar(5)(1:1) .ne. 'R') optpar(5) = 'C'
c
        if (nopt .lt. 6) then
          optpar (6) = 'sph'
ccc          call textin (' Sphere or parallellopiped (S|P) ?',optpar(6))
        end if
        call upcase (optpar(6))
        call remspa (optpar(6))
        if (optpar(6)(1:1) .ne. 'P') optpar(6) = 'S'
c
        if (nopt .lt. 7) then
          write (optpar(7),'(f6.2)') mfrad
ccc          call textin (' Radius (A) ?',optpar(7))
        end if
        call str2r (optpar(7),xdum,ierr)
        if (ierr .eq. 0) then
          if (xdum .gt. 0.1 .and. xdum .lt. 9.9) then
            mfrad = xdum
          else
            call errcon ('Unreasonable radius ignored')
          end if
        end if
c
        call mfmap (imap,jmap,kmap,maxmap,maxpnt,map,shadow,
     +              optpar(5),optpar(6),mfrad)
c
        call xstats (map(1,kmap),npnt(kmap),edstat(3,kmap),
     +    edstat(4,kmap),edstat(1,kmap),edstat(2,kmap),pdum)
        edstat(5,kmap) = edstat(4,kmap)*edstat(4,kmap)
        call rvalut (' Sum of density in map :',1,pdum)
c
        call proplu (edstat(1,kmap),edstat(2,kmap),maprod(kmap),
     +    maplus(kmap))
c
        change (kmap) = .true.
        incore (kmap) = .true.
        commnt (kmap) = 'LC_map'
c
c ... RS_FIT
c
      else if (optpar(1)(1:2) .eq. 'RS') then
c
        if (nopt .lt. 2) call textin (' Which obs_map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which calc_map ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Obs_map and calc_map are identical')
          goto 10
        end if
c
        if (.not. samgrd(imap,jmap)) goto 10
c
        if (nopt .lt. 4) then
          optpar (4) = 'm1.pdb'
          call textin (' Model PDB file ?',optpar(4))
        end if
c
        if (nopt .lt. 5) then
          optpar (5) = 'rs_fit.list'
          call textin (' RS-fit list output file ?',optpar(5))
        end if
c
        if (nopt .ge. 6) then
          call str2r (optpar(6),xdum,ierr)
          if (ierr .eq. 0) then
            if (xdum .gt. 0.1 .and. xdum .lt. 9.9) then
              rsfrad = xdum
              call fvalut (' Entered radius (A) :',1,rsfrad)
            else if (xdum .lt. -0.1 .and. xdum .gt. -9.9) then
              xdum = -1.0 * xdum
              call fvalut (' Assuming resolution (A) :',1,xdum)
              if (xdum .le. 0.6) then
                rsfrad = 0.7
              else if (xdum .ge. 3.0) then
                rsfrad = 0.5 * xdum
              else
                rsfrad = 0.7 + ( (xdum-0.6) / 3.0 )
              end if
              call fvalut (' Calculated radius (A) :',1,rsfrad)
            else
              call errcon ('Unreasonable radius ignored')
              call fvalut (' Keep old radius (A) :',1,rsfrad)
            end if
          end if
        end if
c
        if (nopt .lt. 7) optpar (7) = 'M'
c
        call rsfit (imap,jmap,maxmap,maxpnt,map,shadow,
     +              optpar(4),iunit,optpar(5),junit,rsfrad,
     +              optpar(7))
c
c ... MSTATS
c
      else if (optpar(1)(1:2) .eq. 'MS') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which mask ?',optpar(3))
        jmap = whichm(optpar(3),maxmap,maxpnt)
c
        if (jmap .le. 0 .or. jmap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(jmap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (imap .eq. jmap) then
          call errcon ('Map1 and Map2 are identical')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = '0.5'
          call textin (' Cut-off for mask ?',optpar(4))
        end if
c
        call str2r (optpar(4),xdum,ierr)
        if (ierr .ne. 0) goto 10
c
        if (.not. simap(imap,jmap)) goto 10
c
        call textut (' Masked stats Map  :',name(imap))
        call textut ('              Mask :',name(jmap))
c
        call mstats (imap,jmap,xdum,buffer,maxpnt,maxmap,map)
c
c ... MAPPAGE
c
      else if (optpar(1)(1:2) .eq. 'MA') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Brick file ?',optpar(3))
c
        call mappag (imap,optpar(3),optpar(1),
     +    maxmap,maxpnt,map,shadow)
c
        change(imap) = .false.
c
c ... BRIX
c
      else if (optpar(1)(1:2) .eq. 'BR') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Brix file ?',optpar(3))
c
        call mappag (imap,optpar(3),optpar(1),
     +    maxmap,maxpnt,map,shadow)
c
        change(imap) = .false.
c
c ... EXTRACT
c
      else if (optpar(1)(1:2) .eq. 'EX') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = 'extract.pl2'
          call textin (' O2D plot file ?',optpar(3))
          if (length(optpar(3)) .lt. 1) goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'Z'
          call textin (' Type of plane (X/Y/Z) ?',optpar(4))
        end if
        call upcase (optpar(4))
        if (length(optpar(4)) .lt. 1) goto 10
        if (index('XYZ',optpar(4)(1:1)) .le. 0) then
          call errcon ('Type must be X, Y or Z')
          goto 10
        end if
        itype = index('XYZ',optpar(4)(1:1))
c
        if (nopt .lt. 5) then
          write (optpar(5),*) origin(itype,imap)
          call remspa (optpar(5))
          call textin (' Which plane ?',optpar(5))
        end if
        call str2i (optpar(5),iplane,ierr)
        if (ierr .ne. 0) goto 10
c
        call o2dctr (imap,optpar(3),itype,iplane,
     +    maxmap,maxpnt,map)
c
c ... 1D_PROJECT
c
      else if (optpar(1)(1:2) .eq. '1D') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = ' '
          call textin (' O2D plot file ?',optpar(3))
          if (length(optpar(3)) .lt. 1) goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'Z'
          call textin (' ONTO which axis (X/Y/Z) ?',optpar(4))
        end if
        call upcase (optpar(4))
        if (length(optpar(4)) .lt. 1) goto 10
        if (index('XYZ',optpar(4)(1:1)) .le. 0) then
          call errcon ('Type must be X, Y or Z')
          goto 10
        end if
        itype = index('XYZ',optpar(4)(1:1))
c
        call proj1d (imap,optpar(3),itype,
     +    maxmap,maxpnt,map,buffer)
c
c ... 2D_PROJECT
c
      else if (optpar(1)(1:2) .eq. '2D') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = 'project.pl2'
          call textin (' O2D plot file ?',optpar(3))
          if (length(optpar(3)) .lt. 1) goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = 'Z'
          call textin (' ALONG which axis (X/Y/Z) ?',optpar(4))
        end if
        call upcase (optpar(4))
        if (length(optpar(4)) .lt. 1) goto 10
        if (index('XYZ',optpar(4)(1:1)) .le. 0) then
          call errcon ('Type must be X, Y or Z')
          goto 10
        end if
        itype = index('XYZ',optpar(4)(1:1))
c
        if (nopt .lt. 5) then
          write (optpar(5),*) origin(itype,imap)
          call remspa (optpar(5))
          call textin (' Start plane ?',optpar(5))
        end if
        call str2i (optpar(5),iplane,ierr)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 6) then
          write (optpar(6),*) origin(itype,imap)+
     +      extent(itype,imap)-1
          call remspa (optpar(6))
          call textin (' End plane ?',optpar(6))
        end if
        call str2i (optpar(6),jplane,ierr)
        if (ierr .ne. 0) goto 10
c
        call proj2d (imap,optpar(3),itype,iplane,jplane,
     +    maxmap,maxpnt,map,buffer)
c
c ... DH (density histogram)
c
      else if (optpar(1)(1:2) .eq. 'DH') then
c
        if (nopt .lt. 2) call textin (' Which map ?',optpar(2))
        imap = whichm(optpar(2),maxmap,maxpnt)
c
        if (imap .le. 0 .or. imap .gt. maxmap) then
          call errcon ('Invalid map selection')
          goto 10
        end if
c
        if (.not. incore(imap)) then
          call errcon ('Map not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar(3) = 'dens_hist.plt'
          call textin (' O2D plot file ?',optpar(3))
          if (length(optpar(3)) .lt. 1) goto 10
        end if
c
        if (nopt .lt. 4) then
          optpar(4) = '100'
          call textin (' Number of bins ?',optpar(4))
        end if
        call upcase (optpar(4))
        call str2i (optpar(4),idum,ierr)
        if (ierr .ne. 0) goto 10
c
        call distri (imap,optpar(3),idum,maxmap,maxpnt,map)
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
      do i=1,maxmap
        if (incore(i)) then
          if (change(i)) then
            write (*,6040) name(i)(1:leng1(name(i)))
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
      integer function whichm (nam,maxmap,maxpnt)
c
c ... which map does the name "nam" correspond to ?
c
c ... if "*", then return 0, meaning ALL maps
c     if okay, return index of map
c     otherwise:
c     -1 if length = 0
c     -2 if duplicate name
c     -3 if not found
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
c
      integer i,ll,imap,length
c
      character nam*(*)
c
code ...
c
      whichm = -1
      ll = length(nam)
      if (ll .le. 0) return
c
cc      print *,'WHICHM ',maxmap,maxpnt
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
      imap = 0
      do i=1,maxmap
        if (nam(1:ll) .eq. name(i)(1:ll)) then
          if (imap .ne. 0) return
          imap = i
        end if
      end do
c
      whichm = -3
      if (imap .eq. 0) return
c
      whichm = imap
c
      return
      end
c
c
c
      subroutine selecm (nam,maxmap,maxpnt,ierr)
c
c ... figure out which map(s) are to be selected
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
c
      integer ierr,i,imap,whichm
c
      character nam*(*)
c
code ...
c
      ierr = -1
c
cc      print *,'SELECM ',maxmap,maxpnt
c
      nmap = 0
      do i=1,maxmap
        select (i) = .false.
        if (incore(i)) nmap = nmap + 1
      end do
c
      if (nmap .le. 0) then
        call errcon ('No maps in memory')
        ierr = -2
        return
      end if
c
      imap = whichm(nam,maxmap,maxpnt)
c
      if (imap .lt. 0 .or. imap .gt. maxmap) then
        call errcon ('Invalid map name')
        ierr = -3
        return
      end if
c
      if (imap .eq. 0) then
        do i=1,maxmap
          select (i) = incore (i)
        end do
        ierr = 0
      else
        if (incore(imap)) then
          select (imap) = .true.
          ierr = 0
        else
          ierr = -4
          call errcon ('Selected map not in memory')
        end if
      end if
c
      return
      end
c
c
c
      subroutine writem (imap,ifmt,cutoff,how,ierr,
     +    maxmap,maxpnt,map,buffer,shadow)
c
c ... write map IMAP to file
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap),buffer(maxpnt)
      integer shadow(maxpnt)
c
      real xscale,x1,x2,cutoff
c
      integer ierr,imap,iunit,ifmt,i,icnt
c
      logical xinter
c
      character how*(*)
c
code ...
c
      ierr = -1
      iunit = 1
c
      if (ifmt .eq. 1) then
c
        call prompt (' Writing CCP4 map ...')
c
        call edout (file(imap),iunit,map(1,imap),
     +    origin(1,imap),extent(1,imap),grid(1,imap),
     +    uvw(1,imap),cell(1,imap),spaceg(imap))
c
        close (iunit)
c
      else if (ifmt .eq. 2) then
c
        call prompt (' Writing OLDEZD map ...')
        call xopxua (iunit,file(imap),xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... get scaling constant to fit format (F6.1)
c
        if (edstat(1,imap) .lt. 0.0 .and.
     +      edstat(2,imap) .gt. 0.0) then
          x1 = -999.0/(edstat(1,imap))
          x2 = 9999.0/(edstat(2,imap))
          xscale = min(x1,x2)
        else if (edstat(2,imap) .gt. 0.0) then
          xscale = 9999.0/(edstat(2,imap))
        else
          xscale = 9999.0/(edstat(1,imap))
        end if
c
        do i=1,npnt(imap)
          buffer (i) = xscale*map(i,imap)
        end do
c
        call wroezd (buffer,extent(1,imap),extent(2,imap),
     +    extent(3,imap),origin(1,imap),grid(1,imap),
     +    cell(1,imap),iunit,xscale,'(13f6.1)',ierr)
c
        close (iunit)
        if (ierr .ne. 0) return
c
      else if (ifmt .eq. 3) then
c
        call prompt (' Writing MASK ...')
        call xopxua (iunit,file(imap),xinter(),ierr)
        if (ierr .ne. 0) return
c
ccc        call fvalut (' Border between 0 and 1 :',1,cutoff)
        write (*,'(1x,a,a2,1pe12.4)')
     +   'Set mask points if map value ',how(1:2),cutoff
c
        icnt = 0
        if (how(1:2) .eq. '>=' .or. how(1:2) .eq. 'GE') then
          do i=1,npnt(imap)
            shadow (i) = 0
            if (map(i,imap) .ge. cutoff) then
              shadow(i) = 1
              icnt = icnt + 1
            end if
          end do
        else if (how(1:2) .eq. '<=' .or. how(1:2) .eq. 'LE') then
          do i=1,npnt(imap)
            shadow (i) = 0
            if (map(i,imap) .le. cutoff) then
              shadow(i) = 1
              icnt = icnt + 1
            end if
          end do
        else if (how(1:2) .eq. '< ' .or. how(1:2) .eq. 'LT') then
          do i=1,npnt(imap)
            shadow (i) = 0
            if (map(i,imap) .lt. cutoff) then
              shadow(i) = 1
              icnt = icnt + 1
            end if
          end do
        else
          do i=1,npnt(imap)
            shadow (i) = 0
            if (map(i,imap) .gt. cutoff) then
              shadow(i) = 1
              icnt = icnt + 1
            end if
          end do
        end if
        call jvalut (' Nr of mask points set :',1,icnt)
c
        call wrmask (iunit,shadow,extent(1,imap),
     +    extent(2,imap),extent(3,imap),origin(1,imap),
     +    grid(1,imap),cell(1,imap),'OMASK',ierr)
c
        close (iunit)
        if (ierr .ne. 0) return
c
      else if (ifmt .eq. 4 .or. ifmt .eq. 11) then
c
        call prompt (' Writing (NEW)EZD map ...')
        call xopxua (iunit,file(imap),xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... get scaling constant to fit format (F6.0)
c
        if (edstat(1,imap) .lt. 0.0 .and.
     +      edstat(2,imap) .gt. 0.0) then
          x1 = -999.0/(edstat(1,imap))
          x2 = 9999.0/(edstat(2,imap))
          xscale = min(x1,x2)
        else if (edstat(2,imap) .gt. 0.0) then
          xscale = 9999.0/(edstat(2,imap))
        else
          xscale = 9999.0/(edstat(1,imap))
        end if
c
        call rvalut (' Scale :',1,xscale)
        do i=1,npnt(imap)
          buffer (i) = xscale*map(i,imap)
        end do
c
        call ezdput (iunit,
     +    cell(1,imap),cell(2,imap),cell(3,imap),
     +    cell(4,imap),cell(5,imap),cell(6,imap),
     +    origin(1,imap),origin(2,imap),origin(3,imap),
     +    extent(1,imap),extent(2,imap),extent(3,imap),
     +    grid(1,imap),grid(2,imap),grid(3,imap),
     +    npnt(imap),buffer,xscale,10,1)
c
        close (iunit)
        if (ierr .ne. 0) return
c
      else if (ifmt .eq. 5) then
c
        call prompt (' Writing ENVELOPE (CCP4 MASK) ...')
c
        call fvalut (' Border between 0 and 1 :',1,cutoff)
        icnt = 0
        do i=1,npnt(imap)
          buffer (i) = 0.0
          if (map(i,imap) .ge. cutoff) then
            buffer(i) = 1.0
            icnt = icnt + 1
          end if
        end do
        call jvalut (' Nr of mask points set :',1,icnt)
c
        call edmout (file(imap),iunit,buffer,
     +    origin(1,imap),extent(1,imap),grid(1,imap),
     +    uvw(1,imap),cell(1,imap),spaceg(imap))
c
        close (iunit)
c
      else if (ifmt .eq. 6) then
c
c ... oops ! (20061124) - maps are not scaled, so the numbers
c     in the last line should be 0 and 1 !
c
        x1 = 0.0
        x2 = 1.0
c
        call prompt (' Writing ASCII X-PLOR map ...')
        call wxplmp (file(imap),iunit,map(1,imap),
     +    extent(1,imap),extent(2,imap),extent(3,imap),
     +    origin(1,imap),grid(1,imap),cell(1,imap),
     +    x1,x2,ierr)
c
cccc     +    edstat(3,imap),edstat(4,imap),ierr)
c
        close (iunit)
        if (ierr .ne. 0) return
c
      else if (ifmt .eq. 12) then
c
        call prompt (' Writing binary AMORE map ...')
        call wamore (file(imap),iunit,map(1,imap),
     +    extent(1,imap),extent(2,imap),extent(3,imap),
     +    origin(1,imap),grid(1,imap),cell(1,imap),ierr)
c
        close (iunit)
        if (ierr .ne. 0) return
c
      else if (ifmt .eq. 15) then
c
        call prompt (' Writing binary MPI (EM) map ...')
        call wrimpi (file(imap),iunit,map(1,imap),
     +    extent(1,imap),extent(2,imap),extent(3,imap),ierr)
c
        close (iunit)
        if (ierr .ne. 0) return
c
      else
        call errcon ('Invalid output format')
        ierr = -2
        return
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      logical function simap (imap,jmap)
c
c ... check if IMAP and JMAP are on same grid etc.
c
      include 'mapman.incl'
c
      integer imap,jmap,i
c
code ...
c
      simap = .false.
c
      do i=1,6
        if (abs(cell(i,imap)-cell(i,jmap)) .ge. 0.01) then
          call errcon ('Maps have different cell constants')
          return
        end if
      end do
c
      do i=1,3
        if (grid(i,imap) .ne. grid(i,jmap)) then
          call errcon ('Maps are on different grids')
          return
        end if
      end do
c
      simap = .true.
c
      return
      end
c
c
c
      logical function samgrd (imap,jmap)
c
c ... check if IMAP and JMAP are on exactly same grid etc.
c
      include 'mapman.incl'
c
      integer imap,jmap,i
c
code ...
c
      samgrd = .false.
c
      do i=1,6
        if (abs(cell(i,imap)-cell(i,jmap)) .ge. 0.01) then
          call errcon ('Maps have different cell constants')
          return
        end if
      end do
c
      do i=1,3
        if (grid(i,imap) .ne. grid(i,jmap)) then
          call errcon ('Maps are on different grids')
          return
        end if
        if (origin(i,imap) .ne. origin(i,jmap)) then
          call errcon ('Maps have different origins')
          return
        end if
        if (extent(i,imap) .ne. extent(i,jmap)) then
          call errcon ('Maps have different extents')
          return
        end if
      end do
c
      samgrd = .true.
c
      return
      end
c
c
c
      subroutine rsfit (imap,jmap,maxmap,maxpnt,map,shadow,
     +                  pdbfil,iunit,lstfil,junit,rsfrad,howscl)
c
c ... calc real-space fit of model
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
      integer shadow(maxpnt)
c
      real rsfrad
c
      integer ierr,imap,jmap,iunit,junit,leng1,idum,jdum,kdum
      integer length,ii1,ii2,jj1,jj2
c
      logical xinter,proscl
c
      character pdbfil*(*),lstfil*(*),line*80,howscl*(*)
c
code ...
c
      write (*,6000) 
     +  name(imap)(1:leng1(name(imap))),
     +  name(jmap)(1:leng1(name(jmap))),
     +  pdbfil(1:leng1(pdbfil)),
     +  lstfil(1:leng1(lstfil)),
     +  rsfrad,howscl
c
 6000 format (' Real-space fit calculation:'/
     +  ' Obs_map (e.g., 2mFo-DFc,PHIc)   : ',a/
     +  ' Calc_map (Fc,PHIc)              : ',a/
     +  ' Model PDB file                  : ',a/
     +  ' Result list output file         : ',a/
     +  ' Masking radius around atoms (A) : ',f8.3/
     +  ' Scaling method (Mask|Quick)     : ',a)
c
      call xopxoa (iunit,pdbfil,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening PDB file')
        return
      end if
c
      call xopxua (junit,lstfil,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening list file')
        close (iunit)
        return
      end if
c
      proscl = .true.
      call upcase (howscl)
      if (howscl(1:1) .eq. 'Q') proscl = .false.
c
      call stamp (line)
      write (junit,6010,err=100)
     +  '!'
      write (junit,6010,err=100)
     +  '! ',line(1:leng1(line))
      write (junit,6010,err=100)
     +  '!'
      write (junit,6010,err=100)
     +  '! RS-fit list file ',lstfil(1:leng1(lstfil))
      write (junit,6010,err=100)
     +  '!'
      write (junit,6010,err=100)
     +  '! Obs_map  ',name(imap)(1:leng1(name(imap))),' = ',
     +  file(imap)(1:leng1(file(imap)))
      write (junit,6010,err=100)
     +  '! Calc_map ',name(jmap)(1:leng1(name(jmap))),' = ',
     +  file(jmap)(1:leng1(file(jmap)))
      write (junit,6010,err=100)
     +  '! Model = ',pdbfil(1:leng1(pdbfil))
c
      write (line,'(a,f8.3,a)') 'Masking radius around atoms = ',
     +  rsfrad,' A'
      write (junit,6010,err=100)
     +  '! ',line(1:leng1(line))
c
      if (proscl) then
        write (junit,6010,err=100)
     +    '! Using masked scaling of obs_map and cal_map'
      else
        write (junit,6010,err=100)
     +    '! Using quick-n-dirty scaling of obs_map and cal_map'
      end if
c
 6010 format (10a)
c
  100 continue
      idum = 2*maxpk - 3
      jdum = maxpk/2
      kdum = extent(1,imap)*extent(2,imap)*extent(3,imap)
      ii1 = maxpnt / 3
      ii2 = 2 * ii1
      jj1 = maxrho / 3
      jj2 = 2 * jj1
c
      call grsfit (
     +  map(1,imap),kdum,
     +  extent(1,imap),extent(2,imap),extent(3,imap),
     +  origin(1,imap),grid(1,imap),cell(1,imap),
     +  map(1,jmap),shadow,
     +  shadow(1),shadow(ii1+1),shadow(ii2+1),ii1,
     +  peaks(1,1),peaks(1,jdum+1),idum,
     +  rho(1),rho(jj1+1),rho(jj2+1),jj1,
     +  rsfrad,iunit,junit,proscl,ierr)
c
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine mfmap (imap,jmap,kmap,maxmap,maxpnt,map,shadow,
     +                  what,how,mfrad)
c
c ... calc real-space fit map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
      integer shadow(maxpnt)
c
      real mfrad
c
      integer ierr,imap,jmap,kmap,leng1,idum,jdum,ii1
c
      character what*(*),how*(*)
c
code ...
c
      write (*,6000) 
     +  name(kmap)(1:leng1(name(kmap))),
     +  name(imap)(1:leng1(name(imap))),
     +  name(jmap)(1:leng1(name(jmap))),
     +  what(1:leng1(what)),
     +  how(1:leng1(how)),
     +  mfrad
c
 6000 format (' Local correlation map calculation:'/
     +  ' Correlation map                 : ',a/
     +  ' Map 1                           : ',a/
     +  ' Map 2                           : ',a/
     +  ' Corr Coeff of R-value     (C|R) : ',a/
     +  ' Sphere or parallellopiped (S|P) : ',a/
     +  ' Radius                      (A) : ',f8.3)
c
      idum = 2*maxpk - 3
      jdum = maxpk/2
      ii1 = maxpnt / 100
c
      call gmfmap (
     +  map(1,kmap),
     +  extent(1,imap),extent(2,imap),extent(3,imap),
     +  origin(1,imap),grid(1,imap),cell(1,imap),
     +  map(1,imap),map(1,jmap),shadow(1),ii1,
     +  peaks(1,1),peaks(1,jdum+1),jdum,
     +  what,how,mfrad,ierr)
c
      return
      end
c
c
c
      subroutine simmap (imap,jmap,maxmap,maxpnt,map)
c
c ... similarity of maps IMSK and JMSK
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer ierr,imap,jmap,ilim(2,3),jlim(2,3)
c
code ...
c
      call comap (imap,jmap,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('SIM',
     +  map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +  map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +  ilim,jlim)
c
      return
      end
c
c
c
      subroutine comap (imap,jmap,ilim,jlim,ierr)
c
c ... find common ground between two maps
c
      include 'mapman.incl'
c
      integer ierr,i,imap,jmap,ilim(2,3),jlim(2,3),ncom
      integer imin(3),jmin(3),imax(3),jmax(3),lmin(3),lmax(3)
c
code ...
c
      ierr = -1
c
c ... find absolute margins for both maps
c
      do i=1,3
        imin (i) = origin(i,imap)
        jmin (i) = origin(i,jmap)
        imax (i) = origin(i,imap) + extent(i,imap) - 1
        jmax (i) = origin(i,jmap) + extent(i,jmap) - 1
      end do
c
c ... determine common box
c
      do i=1,3
        lmin (i) = max (imin(i),jmin(i))
        lmax (i) = min (imax(i),jmax(i))
      end do
c
c ... check if they actually overlap
c     2004-06-02 - changed ".ge." to ".gt." to allow overlaps
c                  of only 1 (esp. for 2D maps)
c
      do i=1,3
        if (lmin(i) .gt. lmax(i)) then
          call errcon ('Maps do not overlap')
          ierr = -2
          return
        end if
      end do
c
      ierr = 0
c
c ... find common limits
c
      do i=1,3
        ilim (1,i) = lmin(i) - origin(i,imap) + 1
        ilim (2,i) = lmax(i) - origin(i,imap) + 1
        jlim (1,i) = lmin(i) - origin(i,jmap) + 1
        jlim (2,i) = lmax(i) - origin(i,jmap) + 1
      end do
c
      call ivalut (' Lower limits common volume :',3,lmin)
      call ivalut (' Upper limits common volume :',3,lmax)
      call ivalut (' Limits first  map  :',6,ilim)
      call ivalut (' Limits second map  :',6,jlim)
c
      ncom = (lmax(1)-lmin(1)+1) * (lmax(2)-lmin(2)+1) *
     +       (lmax(3)-lmin(3)+1)
      call jvalut (' Number of common map points :',1,ncom)
c
      return
      end
c
c
c
      subroutine copymp (imap,jmap,maxmap,maxpnt,map,lcopy)
c
c ... copy JMAP into IMAP; if LCOPY = .TRUE. copy density,
c     otherwise, just copy map parameters
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer i,imap,jmap
c
      logical lcopy
c
code ...
c
      do i=1,3
        cell (i,imap)   = cell (i,jmap)
        cell (i+3,imap) = cell (i+3,jmap)
        origin (i,imap) = origin (i,jmap)
        extent (i,imap) = extent (i,jmap)
        grid (i,imap)   = grid (i,jmap)
        uvw (i,imap)    = uvw (i,jmap)
      end do
c
      select (imap) = .false.
      incore (imap) = .true.
      change (imap) = .true.
      file (imap)   = 'not_defined'
      npnt (imap)   = npnt (jmap)
      spaceg (imap) = spaceg (jmap)
      commnt (imap) = 'Copied from '//name(jmap)
      maprod (imap) = maprod (jmap)
      maplus (imap) = maplus (jmap)
c
      do i=1,5
        edstat (i,imap) = edstat (i,jmap)
      end do
c
      if (lcopy) then
        call copmap (map(1,imap),map(1,jmap),extent(1,imap),
     +    extent(2,imap),extent(3,imap))
      end if
c
      return
      end
c
c
c
      subroutine rough (imap,jmap,shadow,maxmap,maxpnt,map,
     +    npos,nneg,ibox,filnam,nbin)
c
c ... calculate local roughness map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
      integer shadow(maxpnt)
c
      integer imap,jmap,npos,nneg,ibox,nbin
c
      character filnam*(*)
c
code ...
c
      call irough (map(1,imap),map(1,jmap),shadow,
     +  extent(1,jmap),extent(2,jmap),
     +  extent(3,jmap),npos,nneg,ibox,edstat(1,imap),edstat(2,imap),
     +  edstat(3,imap),edstat(4,imap))
c
      edstat(5,imap) = edstat(4,imap) * edstat(4,imap)
c
      change (imap) = .true.
c
      call distri (imap,filnam,nbin,maxmap,maxpnt,map)
c
      return
      end
c
c
c
      subroutine invmap (imap,jmap,maxmap,maxpnt,map)
c
c ... invert JMAP into IMAP
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer i,imap,jmap
c
code ...
c
      do i=1,3
        cell (i,imap)   = cell (i,jmap)
        cell (i+3,imap) = cell (i+3,jmap)
        origin (i,imap) = origin (i,jmap)
        extent (i,imap) = extent (i,jmap)
        grid (i,imap)   = grid (i,jmap)
        uvw (i,imap)    = uvw (i,jmap)
      end do
c
      select (imap) = .false.
      incore (imap) = .true.
      change (imap) = .true.
      file (imap)   = 'not_defined'
      npnt (imap)   = npnt (jmap)
      spaceg (imap) = spaceg (jmap)
      commnt (imap) = 'QInverted from '//name(jmap)
      maprod (imap) = maprod (jmap)
      maplus (imap) = maplus (jmap)
c
      do i=1,5
        edstat (i,imap) = edstat (i,jmap)
      end do
c
      call dominv (map(1,imap),map(1,jmap),extent(1,imap),
     +  extent(2,imap),extent(3,imap))
c
      return
      end
c
c
c
      subroutine allocm (string,imap,ierr,maxmap,maxpnt)
c
c ... allocate a map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
c
      integer length,whichm
      integer i,imap,ierr
c
      character string*(*)
c
code ...
c
      ierr = -1
c
      if (length(string) .lt. 1) return
c
      do i=1,maxmap
        if (.not. incore(i)) then
          imap = i
          goto 910
        end if
      end do
c
      call errcon ('No more maps available')
      ierr = -2
      return
c
  910 continue
      name (imap) = string
c
      call upcase (name(imap))
      if (whichm(name(imap),maxmap,maxpnt) .ne. imap) then
        call errcon ('Invalid map name (empty or not unique)')
        name (imap) = '$%#@'
        ierr = -3
        return
      end if
c
      if (incore(imap)) then
        call errcon ('Map in use; DELETE it first')
        ierr = -4
        return
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine odlmap (imsk,odlfil,type,mode)
c
c ... write ODL file to draw box around map IMSK
c
      include 'mapman.incl'
c
      real a(3,3),x(3),flo(3),fhi(3)
c
      integer i,imsk,iunit,ierr
c
      logical xinter,lcell,lsolid
c
      character odlfil*(*),type*(*),mode*(*)
c
code ...
c
      iunit = 12
c
      call xopxua (iunit,odlfil,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening ODL file')
        return
      end if
c
      lcell  = (type(1:1) .ne. 'B')
      lsolid = (mode(1:1) .ne. 'L')
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell(1,imsk), a, 0)
c
c ... calc origin and top in fractionals
c
      if (lcell) then
        do i=1,3
          flo(i) = 0.0
          fhi (i) = 1.0
        end do
      else
        do i=1,3
          flo (i) = float(origin(i,imsk)) * 
     +             (cell(i,imsk)/float(grid(i,imsk))) / cell(i,imsk)
          fhi (i) = float(extent(i,imsk)+origin(i,imsk)-1) *
     +             (cell(i,imsk)/float(grid(i,imsk))) / cell(i,imsk)
        end do
      end if
c
      call fvalut (' Origin (fract.) :',3,flo)
      call fvalut (' Top    (fract.) :',3,fhi)
c
      call mulmtx (a, flo, x, 3, 3, 1)
      call fvalut (' Origin (cart.)  :',3,x)
      call mulmtx (a, fhi, x, 3, 3, 1)
      call fvalut (' Top    (cart.)  :',3,x)
c
c ... now write the ODL object
c
      if (lcell) then
        write (iunit,6000) 'begin mapcell'
        if (lsolid) write (iunit,6000) 'mode solid'
        write (iunit,6000) 'colour magenta'
      else
        write (iunit,6000) 'begin mapbox'
        if (lsolid) write (iunit,6000) 'mode solid'
        write (iunit,6000) 'colour cyan'
      end if
c
      if (lsolid) then
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'mode line'
c
      else
c
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
c
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
c
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
c
      end if
c
      write (iunit,6000) 'end_object'
c
 6000 format (a,3f10.3)
c
      close (iunit)
      write (*,*) 'ODL file written'
c
      return
      end
c
c
c
      subroutine maprd (imap,ifmt,ierr,
     +    maxmap,maxpnt,map,buffer)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap),buffer(maxpnt)
c
      real xdum,ave,sdv
c
      integer imap,ifmt,ierr,i,jfmt,level,iunit,imax,jmax,kmax
c
code ...
c
      ierr = -1
      jfmt = -ifmt
      level = 0
      iunit = 1
c
      call prompt (' Read header')
      incore (imap) = .false.
      call maphdr (file(imap), iunit, jfmt,
     +  origin(1,imap), extent(1,imap), grid(1,imap),
     +  uvw(1,imap), cell(1,imap), spaceg(imap),
     +  rho, maxrho, map(1,imap), maxpnt, incore(imap))
c
      if (spaceg(imap).le.0) then
        call errcon ('While opening map file')
        ierr = -2
        return
      end if
c
      call prompt (' Header done')
c
      if (.not. incore(imap)) then
        call prompt (' Not in core')
        imax = extent( uvw(3,imap), imap)
        jmax = extent( uvw(1,imap), imap)
        kmax = extent( uvw(2,imap), imap)
        call prompt (' Reading levels ...')
        do i=1,imax
          if (i .eq. i/10*10) 
     +      call ivalut (' Level number :', 1, i)
          call mapin (file(imap), iunit, jfmt,
     +      origin(1,imap), extent(1,imap), grid(1,imap),
     +      uvw(1,imap), cell(1,imap),
     +      rho, maxrho, buffer, maxrho, incore(imap), level)
c
          if (level.lt.0) goto 9020
c
          call pckrho (map(1,imap), extent(1,imap),
     +      extent(2,imap),extent(3,imap), i,
     +      rho, jmax, kmax, uvw(1,imap))
        end do
      end if
c
      npnt (imap) = extent(1,imap)*extent(2,imap)*extent(3,imap)
c
c ... if X-PLOR check if it's X-PLOR version 4
c
      if (ifmt .eq. 5) then
c
        ave = 0.0
        sdv = -1.0
c
        read (iunit,'(i8)',err=200,end=200) i
        if (i .ne. -9999) goto 200
        read (iunit,'(1x,2e12.5)',err=200,end=200) ave,sdv
        if (ave .eq. 0.0 .or. sdv .le. 0.0) goto 200
c
        call prompt (' This is a CNS (or X-PLOR version 4) map !')
        call rvalut (' Average density :',1,ave)
        call rvalut (' Sigma level     :',1,sdv)
        call prompt (' Restoring original density values (e-/A3)')
        do i=1,npnt(imap)
          map(i,imap) = (sdv * map(i,imap)) + ave
        end do
c
  200   continue
c
      end if
c
      ierr = 0
c
c      print *,' first  = ',map(1,imap)
c      print *,' second = ',map(2,imap)
c
c      call xstats (map(1,imap),npnt(imap),edstat(3,imap),
c     +  edstat(4,imap),edstat(1,imap),edstat(2,imap),xdum)
c      edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c      call rvalut (' Sum of density in map :',1,xdum)
c
 9999 continue
c
c ... close map properly
c
      call mapclo (iunit)
c
      return
c
c ... error traps
c
 9020 call errcon ('While reading map')
      ierr = -3
      goto 9999
c
      end
c
c
c
      subroutine mappag (imap,brick,how,
     +    maxmap,maxpnt,map,shadow)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
      integer shadow(maxpnt)
c
      integer imap,iunit,length,nxyz(3),i1,i2,i,indxra,ierr
      integer*2 first(256)
c
      logical litend,little,lbrix,lturbo
c
      character brick*(*),how*(*),str*512
c
code ...
c
      lturbo = (how(1:2) .eq. 'TU')
      lbrix  = (how(1:2) .eq. 'BR')
c
      if (lbrix) then
        call textut (' Opening BRIX file :',brick)
      else if (lturbo) then
        call textut (' Opening TURBO-DSN6 file :',brick)
      else
        call textut (' Opening DSN6 file :',brick)
      end if
c
      if (length(brick) .le. 0) then
        call errcon ('No filename supplied')
        return
      end if
c
      if (extent(1,imap)*extent(2,imap) .gt. maxrho) then
        call errcon ('Slices of this map are too big')
        call jvalut (' Size :',1,(extent(1,imap)*extent(2,imap)))
        call jvalut (' Max  :',1,maxrho)
        return
      end if
c
      little = litend()
      if (little) then
        call prompt (' Little-endian machine; will swap bytes')
      else
        call prompt (' Big-endian machine; will NOT swap bytes')
      end if
c
      iunit = 12
      close (iunit)
c
C+PRE
C  Note that there is a dispute between different computers about how
C  to count record-lengths in unformatted direct access files, either
C  in bytes or words (=4 bytes).  The record length here needs to be
C  512 bytes == 128 words. KRECL is set in the MAPMAN include file
c
c      open (iunit, file=brick, status='unknown',
c     +      form='unformatted', access='direct',
c     +      recl=krecl, maxrec=mmxrec, err=13)
c
      open (iunit, file=brick, status='unknown',
     +      form='unformatted', access='direct',
     +      recl=krecl, err=13)
c
C-PRE
c
      if (lbrix) then
c
        write (str, 10) (origin(i,imap),i=1,3), 
     +    (extent(i,imap),i=1,3), (grid(i,imap),i=1,3),
     +    (cell(i,imap),i=1,6), maprod(imap), maplus(imap),
     +    edstat(4,imap)
c
 10     format (":-) Origin", 3i5," Extent", 3i5, " Grid", 3i5,
     $        " Cell ", 6f10.3, " Prod", f12.5, " Plus",i8, 
     $        " Sigma ", f12.5)
c
        write (iunit, rec=1, iostat=ierr) str
        indxra = 1
c
      else
c
c ... fill up first record
c
        i1 = 80
        i2 = 100
c
c ... check if I1*MAX(CELL) fits in INTEGER*2 !!!
c
        do i=1,6
          if ( (i1*cell(i,imap)) .gt. 32760) then
            i1 = min (i1, int
     +            ( 32760.0 / float(int(0.999 + cell(i,imap)))))
            i1 = max (1, i1)
          end if
        end do
        call ivalut (' Scale constant for cell :',1,i1)
c
        do i=1,3
          first (i  ) = origin (i,imap)
          first (i+3) = extent (i,imap)
          first (i+6) = grid (i,imap)
          if (lturbo) then
            first (i+9) = i2*cell(i,imap)
            first (i+12) = i2*cell(i+3,imap)
          else
            first (i+9) = i1*cell(i,imap)
            first (i+12) = i1*cell(i+3,imap)
          end if
        end do
        first (16) = i2*maprod(imap)
        first (17) = maplus(imap)
        first (18) = i1
        first (19) = i2
c
c ... 930514 - store RHOMIN, RHOMAX, RHOSIGMA
c
        first (20) = edstat(1,imap)*maprod(imap) + maplus(imap)
        first (21) = edstat(2,imap)*maprod(imap) + maplus(imap)
        first (22) = edstat(4,imap)*maprod(imap) + maplus(imap)
c
        do i=23,256
          first (i) = 0
        end do
c
        indxra = 1
        if (little) call bytswp (first)
        write (iunit,rec=indxra,err=14) first
c
      end if
c
      call rvalut (' Prod :',1,maprod(imap))
      call jvalut (' Plus :',1,maplus(imap))
c
      do i=1,3
        nxyz(i) = extent(i,imap)/8
        if ( mod(extent(i,imap),8) .ge. 1) nxyz(i) = nxyz(i) + 1
      end do
      call ivalut (' Pages along X, Y and Z  :',3,nxyz)
      i = nxyz(1)*nxyz(2)*nxyz(3)
      call jvalut (' Total nr of pages :',1,i)
c
c      if (i .gt. mmxrec) then
c        call errcon ('Too many pages required !')
c        call jvalut (' Max nr of records :',1,mmxrec)
c        goto 15
c      end if
c
      if (lbrix) then
        do i=1,extent(3,imap)
          call storem (map(1,imap),extent(1,imap),
     +      extent(2,imap),extent(3,imap),i,rho)
          call pageb (iunit,rho,shadow,extent(1,imap),extent(2,imap),
     +      nxyz(1),nxyz(2),i,indxra,maprod(imap),maplus(imap),ierr)
          if (ierr .ne. 0) goto 15
c
c          if ( (indxra+1) .ge. mmxrec ) then
c            call errcon ('Too many records required !')
c            call jvalut (' Max nr of pages :',1,mmxrec)
c            goto 15
c          end if
c 
       end do
        i = 0
        call pageb (iunit,rho,shadow,extent(1,imap),extent(2,imap),
     +    nxyz(1),nxyz(2),i,indxra,maprod(imap),maplus(imap),ierr)
        if (ierr .ne. 0) goto 15
      else
        do i=1,extent(3,imap)
          call storem (map(1,imap),extent(1,imap),
     +      extent(2,imap),extent(3,imap),i,rho)
          call paged (rho,extent(1,imap),extent(2,imap),i,
     +      maprod(imap),maplus(imap),shadow,maxrho,iunit,
     +      indxra,nxyz(1),nxyz(2),nxyz(3),ierr)
          if (ierr .ne. 0) goto 15
c
c          if ( (indxra+1) .ge. mmxrec ) then
c            call errcon ('Too many records required !')
c            call jvalut (' Max nr of pages :',1,mmxrec)
c            goto 15
c          end if
c
        end do
        call rest (extent(1,imap),extent(2,imap),extent(3,imap),
     +    maprod(imap),maplus(imap),shadow,maxrho,iunit,
     +    indxra,nxyz(1),nxyz(2),nxyz(3),ierr)
        if (ierr .ne. 0) goto 15
      end if
c
        close (iunit)
c
ccc      call jvalut (' Nr of records written :',1,indxra)
c
      return
c
c ... error traps
c
   13 continue
      call errcon ('While opening file')
      close (iunit)
      return
c
   14 continue
      call errcon ('While writing header')
      close (iunit)
      return
c
   15 continue
      call errcon ('While writing data')
      close (iunit)
      return
c
      end
c
c
c
      subroutine o2dctr (imap,filenm,itype,iplane,
     +    maxmap,maxpnt,map)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer nlevel
      parameter (nlevel = 8)
c
      integer imap,itype,iplane,iunt,ierr,j,i1,i2,i,leng1
c
      logical xinter
c
      character filenm*(*),myline*256,plname(3)*1
c
      data plname /'X','Y','Z'/
c
code ...
c
      if (itype .lt. 1 .or. itype .gt. 3) then
        call errcon ('Invalid plane type')
        return
      end if
c
      i1 = origin(itype,imap)
      i2 = (origin(itype,imap)+extent(itype,imap)-1)
      if (iplane .lt. i1 .or.
     +    iplane .gt. i2) then
        call errcon ('Plane index out of range')
        call ivalut (' Requested   :',1,iplane)
        call ivalut (' Lower limit :',1,i1)
        call ivalut (' Upper limit :',1,i2)
        return
      end if
c
c ... account for offset !!!
c
      iplane = iplane - origin(itype,imap) + 1
c
      iunt = 87
      close (iunt)
      call xopxua (iunt,filenm,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call stamp (myline)
      myline = 'REMARK '//myline
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a,i3)',err=9000) 'NLEVEL ',nlevel
      write (iunt,'(a)',err=9000) 'LEVELS'
      write (myline,'(1p,20(e12.4,1x))',err=9000)
     +  (edstat(3,imap)+float(j+1)*0.5*edstat(4,imap),
     +  j=1,nlevel)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a)',err=9000) 'COLOUR'
      write (myline,'(20i5)',err=9000) ( (j+1)/2 ,j=1,nlevel)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      if (itype .eq. 1) then
        i1 = 2
        i2 = 3
      else if (itype .eq. 2) then
        i1 = 1
        i2 = 3
      else if (itype .eq. 3) then
        i1 = 1
        i2 = 2
      end if
c
      write (myline,*,err=9000)
     +  'REMARK Contour plot of ',plname(i1),'-',
     +  plname(i2),' plane nr ',(iplane + origin(itype,imap) - 1)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Map ',name(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK File ',file(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Comment ',commnt(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,6f10.2)',err=9000)
     +  'REMARK Cell ',(cell(i,imap),i=1,6)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,3i8)',err=9000)
     +  'REMARK Grid ',(grid(i,imap),i=1,3)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000)
     +  'XLABEL ',plname(itype),'-plane',
     +  iplane+origin(itype,imap)-1,
     +  '; horizontal = ',plname(i1),'-axis'
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000)
     +  'YLABEL vertical = ',plname(i2),'-axis'
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a,i3)',err=9000) 'XPOINT ',extent(i1,imap)
      write (iunt,'(a,i3)',err=9000) 'YPOINT ',extent(i2,imap)
c
      write (myline,*,err=9000) 'XLIMIT ',
     +  float(origin(i1,imap))/float(grid(i1,imap)),
     +  float(origin(i1,imap)+extent(i1,imap)-1)/float(grid(i1,imap))
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000) 'YLIMIT ',
     +  float(origin(i2,imap))/float(grid(i2,imap)),
     +  float(origin(i2,imap)+extent(i2,imap)-1)/float(grid(i2,imap))
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      call gjk001 (map(1,imap),i1,i2,itype,iplane,extent(1,imap),
     +  extent(2,imap),extent(3,imap),iunt,ierr)
      if (ierr .ne. 0) goto 9000
c
      write (iunt,'(a)',err=9000) 'END'
c
      call prompt (' O2D file written')
      close (iunt)
      return
c
 9000 continue
      call errcon ('While writing O2D file')
      close (iunt)
      return
c
      end
c
c
c
      subroutine addmap (imap,jmap,ierr,maxmap,maxpnt,map)
c
c ... add two maps
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer imap,jmap,ierr,ilim(2,3),jlim(2,3)
c
code ...
c
      ierr = -1
c
      call comap (imap,jmap,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('ADD',
     +  map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +  map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +  ilim,jlim)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine opemap (imap,jmap,oper,ierr,maxmap,maxpnt,map)
c
c ... combine two maps
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer imap,jmap,ierr,ilim(2,3),jlim(2,3)
c
      character oper*(*)
c
code ...
c
      ierr = -1
c
      call comap (imap,jmap,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      ierr = 0
c
      if (oper(1:1) .eq. '+') then
        call combim ('GK+',
     +    map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +    map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +    ilim,jlim)
      else if (oper(1:1) .eq. '-') then
        call combim ('GK-',
     +    map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +    map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +    ilim,jlim)
      else if (oper(1:1) .eq. '*') then
        call combim ('GK*',
     +    map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +    map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +    ilim,jlim)
      else
        call errcon ('Invalid operator in OPEMAP !')
        ierr = -2
      end if
c
      return
      end
c
c
c
      subroutine minmax (imap,jmap,ierr,maxmap,maxpnt,map)
c
c ... generate MIN(map1,map2) and MAX(map1,map2)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer imap,jmap,ierr,ilim(2,3),jlim(2,3)
c
code ...
c
      ierr = -1
c
      call comap (imap,jmap,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('MIN',
     +  map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +  map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +  ilim,jlim)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine cormap (imap,jmap,ierr,maxmap,maxpnt,map)
c
c ... generate correlated maps
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer imap,jmap,ierr,ilim(2,3),jlim(2,3)
c
code ...
c
      ierr = -1
c
      call comap (imap,jmap,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('COR',
     +  map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +  map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +  ilim,jlim)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine pickem (imap,jhow,form,
     +  maxmap,maxpnt,map,buffer)
c
c ... pick peaks from map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap),buffer(maxpnt)
c
      real f2c(3,3),xdum(3)
c
      integer imap,i,j,k,icnt,iunit,ierr,mylim(2,3)
      integer jhow,ihow,leng1
c
      logical xinter,lpdb
c
      character line*128,form*(*)
c
code ...
c
      ihow = jhow
c
      iunit = 32
      close (iunit)
      call xopxua (iunit,pkfile,xinter(),ierr)
      if (ierr .ne. 0) return
c
      lpdb = (form(1:1) .eq. 'P')
c
      call rvalut (' Pick level  :',1,pklev)
c
cc      print *,' IHOW = ',ihow
c
c ... clean up borders
c
      do i=1,3
c
        if (pklim(1,i) .ge. pklim(2,i)) then
          pklim (1,i) = origin(i,imap)
          pklim (2,i) = origin(i,imap)+extent(i,imap)-1
        end if
c
        pklim (1,i) = max (pklim(1,i),
     +                     origin(i,imap)+1)
        pklim (2,i) = min (pklim(2,i),
     +                     origin(i,imap)+extent(i,imap)-2)
c
        call ivalut (' Pick limits :',2,pklim(1,i))
c
        mylim (1,i) = pklim(1,i)-origin(i,imap)+1
        mylim (2,i) = pklim(2,i)-origin(i,imap)+1
c
      end do
c
c ... initialise scratch map to -1.0
c
      call inimap (buffer,extent(1,imap),extent(2,imap),
     +             extent(3,imap),-1.0)
c
      call setpkm (buffer,map(1,imap),extent(1,imap),
     +             extent(2,imap),extent(3,imap),
     +             mylim,pklev,1.0,icnt)
c
      if (icnt .le. 0) then
        call errcon ('No potential peaks; no picking')
        return
      end if
c
      if (ihow .eq. 1) then
        call prompt (' Picking ... Integrate intensity')
      else if (ihow .eq. 2) then
        call prompt (' Picking ... High-density values')
      else
        ihow = 0
        call prompt (' Picking ... Interpolate intensity')
      end if
c
      call pick3d (buffer,map(1,imap),extent(1,imap),
     +             extent(2,imap),extent(3,imap),
     +             mylim,pklev,npks,maxpk,peaks,ihow)
c
      if (npks .le. 0) return
c
 6000 format ('REMARK ',a)
c
      call stamp (line)
      write (iunit,6000,err=996) line(1:leng1(line))
c
      write (line,*) npks,' peaks picked above level ',pklev
      call pretty (line)
      write (iunit,6000,err=996) line(1:leng1(line))
c
      line = 'Map file '//file(imap)
      call pretty (line)
      write (iunit,6000,err=996) line(1:leng1(line))
c
      line = 'Map name '//name(imap)
      call pretty (line)
      write (iunit,6000,err=996) line(1:leng1(line))
c
      line = 'Map text '//commnt(imap)
      call pretty (line)
      write (iunit,6000,err=996) line(1:leng1(line))
c
      if (lpdb) then
        write (iunit,6000,err=996)
        write (iunit,'(a,3f9.3,3f7.2)',err=996)
     +    'CRYST1',(cell(i,imap),i=1,6)
      end if
c
      call orthog (cell(1,imap),f2c,0)
c
      j = pfirst - 1
c
      do i=1,npks
c
        do k=1,3
          peaks(k,i) = peaks(k,i) + origin(k,imap) -1
        end do
c
        if (lpdb) then
          write (line,'(a,i8,a,3f10.3,a,1pe15.6)')
     +      'Peak # ',i,' @ ',peaks(1,i),peaks(2,i),
     +      peaks(3,i),' $ ',peaks(4,i)
          write (iunit,6000,err=996) line(1:leng1(line))
        end if
c
c ... convert to Angstrom & write ATOM record
c
        j = j + 1
c
        do k=1,3
          xdum(k) = peaks(k,i) / float(grid(k,imap))
        end do
c
cc        print *,' GRID = ',(peaks(k,i),k=1,3)
cc        print *,' FRAC = ',(xdum(k),k=1,3)
c
        call mulmtx (f2c,xdum,peaks(1,i),3,3,1)
c
cc        print *,' CART = ',(peaks(k,i),k=1,3)
c
        if (lpdb) then
          write (iunit,6010) i,pkatom,pkres,j,peaks(1,i),
     +      peaks(2,i),peaks(3,i),1.0,20.0
        else
          write (iunit,6020) peaks(1,i),peaks(2,i),peaks(3,i),
     +      peaks(4,i)
        end if
c
      end do
c
 6010 format ('ATOM  ',i5,1x,a4,1x,a3,2x,i4,4x,3f8.3,2f6.2)
 6020 format (' SOLUTIONRC     1',3f9.2,'  0.00000  0.00000  0.00000',
     +  f5.1,'  0.0')
c
c SOLUTIONRC     1   116.50   162.44   144.00  0.00000  0.00000  0.00000 22.9  0.0
c123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      if (lpdb) write (iunit,'(a)',err=996) 'END'
c
      close (iunit)
      return
c
c ... fatal write error
c
  996 continue
      call errcon ('While writing peak file')
      close (iunit)
c
      return
      end
c
c
c
      subroutine undsn6 (imap,ierr,maxmap,maxpnt,map,lbrix)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      real xdum
c
      integer imap,ierr,iunit
c
      logical lbrix
c
code ...
c
      ierr = -1
      iunit = 1
c
      if (lbrix) then
        call prompt (' Reconstructing BRIX file')
      else
        call prompt (' Reconstructing DSN6 ("BRICK") file')
      end if
c
      incore (imap) = .false.
c
      call omapin (file(imap), iunit, map(1,imap),
     +  origin(1,imap), extent(1,imap), grid(1,imap),
     +  cell(1,imap), maxpnt, krecl, lbrix, ierr)
c
      if (ierr .ne. 0) return
c
      call prompt (' Map read into memory')
      call telmap (grid(1,imap),origin(1,imap),extent(1,imap),
     +             cell(1,imap))
      call prompt (' Spacegroup set to P1')
c
c      incore (imap) = .true.
c      spaceg (imap) = 1
c      uvw (1,imap) = 2
c      uvw (2,imap) = 1
c      uvw (3,imap) = 3
c
c      npnt (imap) = extent(1,imap)*extent(2,imap)*extent(3,imap)
c      call xstats (map(1,imap),npnt(imap),edstat(3,imap),
c     +  edstat(4,imap),edstat(1,imap),edstat(2,imap),xdum)
c      edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c      call rvalut (' Sum of density in map :',1,xdum)
c
      ierr = 0
c
 9999 continue
      return
c
      end
c
c
c
      subroutine xbones (imap,ierr,
     +  maxmap,maxpnt,map,shadow)
c
c ... skeletonise a map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
      integer shadow(maxpnt)
c
      integer imap,ierr,maxval,ncount
c
code ...
c
      ierr = -1
c
c ... construct integer envelope from the electron density
c
      write (*,6010) bobase,bostep
 6010 format (' Constructing envelope with base level ',f8.3/
     +        ' Level increment ..................... ',f8.3)
c
      call envset (shadow,map(1,imap),extent(1,imap),
     +  extent(2,imap),extent(3,imap),bostep,bobase,
     +  maxval,ierr)
      if (ierr .ne. 0) return
c
      if (maxval .gt. 10) maxval = 10
      write (*,6020) maxval
 6020 format (' Number of levels (max. 10) .......... ',i8)
c
      call edchck (shadow,extent(1,imap),
     +  extent(2,imap),extent(3,imap),ncount,ierr)
      if (ierr .ne. 0) return
c
      write (*,6030) ncount
 6030 format (' Nr of points (level > 0) in input map ',i8)
c
      call sklton (shadow,extent(1,imap),
     +  extent(2,imap),extent(3,imap),ierr)
      if (ierr .ne. 0) return
c
      call edchck (shadow,extent(1,imap),
     +  extent(2,imap),extent(3,imap),ncount,ierr)
      if (ierr .ne. 0) return
c
      write (*,6040) ncount
 6040 format (' Ditto, in skeletonised map .......... ',i8)
c
      call chncon (shadow,extent(1,imap),
     +  extent(2,imap),extent(3,imap),bonex,cell(1,imap),
     +  grid(1,imap),origin(1,imap),
     +  maxbop,maxbob,maxboc,pntlis,bonesa,levlis,brnlis,
     +  conlis,ierr)
      if (ierr .ne. 0) return
c
c ... save LEVLIS
c
cc      do i=1,bonpcn
cc        conlis (1,i) = levlis(i)
cc      end do
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine envsav (i1, i2, i3, i4, i5)
c
c --- Pass skeleton data structure from original bones->mainlv
c
      include 'mapman.incl'
c
      integer i5(3,*)
      integer i1, i2, i3, i4
c
ccc      common /cons/ npcn, nccn, conlis(20000)
ccc      integer npcn, nccn, conlis
c
      integer i
c
code ...
c
      bonpcn = i1
      bonccn = i2
      do i=1,bonccn
        boncon(i) = i5(i4, i)
cc        print *,i,boncon(i)
      end do
c
cc      print *,'PCN, CCN = ',bonpcn,bonccn
c
      return
      end
c
c
c
      subroutine ybones (iunit,ierr)
c
c ... write BONES file
c
      include 'mapman.incl'
c
      integer iunit,ierr
c
code ...
c
      ierr = -1
c
c ... restore LEVLIS
c
cc      do i=1,bonpcn
cc        levlis (i) = conlis (1,i)
cc      end do
c
      call mainlv (maxbop,maxbob,maxboc,bonpcn,bonccn,
     +  bonlen,levlis,boncon,ierr)
      if (ierr .ne. 0) return
c
      call obones (maxbop,maxbob,maxboc,bonpcn,bonccn,
     +  pntlis,levlis,boncon,iunit,bonid,brnlis,ierr)
      if (ierr .ne. 0) return
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine zbones (iunit,ierr)
c
c ... write pruned BONES to PDB file
c
      include 'mapman.incl'
c
      integer iunit,ierr
c
code ...
c
      ierr = -1
c
      call mainlv (maxbop,maxbob,maxboc,bonpcn,bonccn,
     +  bonlen,levlis,boncon,ierr)
      if (ierr .ne. 0) return
c
      call pbones (maxbop,maxbob,maxboc,bonpcn,bonccn,
     +  pntlis,levlis,boncon,iunit,brnlis,conlis,bonlen,
     +  bobfac,cell(1,bonmap),ierr)
      if (ierr .ne. 0) return
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine xpl3dm (imap,ierr,maxmap,maxpnt,map)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      real xdum,xmin(3),xmax(3),ymin(3),ymax(3),space(3),a(3,3)
c
      logical xinter
c
      integer imap,ierr,iunit,length,i,inum,nnum,j,k
c
      character line*128,reply*1
c
code ...
c
      ierr = -1
      iunit = 1
c
 6000 format (A)
c
      call prompt (' Reconstructing XPLOR 3DMATRIX file')
      incore (imap) = .false.
c
c ... open the file and start reading
c
      call xopxoa (iunit,file(imap),xinter(),ierr)
      if (ierr. ne. 0) then
        call errcon ('While opening file')
        ierr = -2
        return
      end if
c
      read (iunit,6000,err=9000,end=9000) line
      if (line(1:3) .ne. ' (*') then
        call errcon ('Not an XPLOR 3DMATRIX file')
        call prompt (' First line should start with " (*"')
        goto 9100
      end if
      i = length(line)
      call textut (' Title :',line(5:i-2))
c
   10 read (iunit,6000,err=9000,end=9000) line
      if (line(1:11) .ne. '   rnumber=') goto 10
      i = length(line)
      read (line(12:i-1),*,err=9200) inum
      call jvalut (' Nr of points  :',1,inum)
c
      if (inum .lt. 9) then
        call errcon ('File too small')
        goto 9100
      end if
c
      if (inum .gt. maxpnt) then
        call errcon ('Too many points in file')
        call jvalut (' Maximum :',1,maxpnt)
        goto 9100
      end if
c
   20 read (iunit,6000,err=9000,end=9000) line
      if (line(1:8) .ne. '   rave=') goto 20
      i = length(line)
      read (line(9:i-1),*,err=9200) xdum
      call rvalut (' Average value :',1,xdum)
c
   30 read (iunit,6000,err=9000,end=9000) line
      if (line(1:10) .ne. '   rsigma=') goto 30
      i = length(line)
      read (line(11:i-1),*,err=9200) xdum
      call rvalut (' Sigma level   :',1,xdum)
c
   50 read (iunit,6000,err=9000,end=9000) line
      if (line(1:8) .ne. '   rmax=') goto 50
      i = length(line)
      read (line(11:i-1),*,err=9200) xdum
      call rvalut (' Maximum value :',1,xdum)
c
   40 read (iunit,6000,err=9000,end=9000) line
      if (line(1:8) .ne. '   rmin=') goto 40
      i = length(line)
      read (line(11:i-1),*,err=9200) xdum
      call rvalut (' Minimum value :',1,xdum)
c
   60 read (iunit,6000,err=9000,end=9000) line
      if (line(1:6) .ne. ' var={') goto 60
      i = length(line)
      call textut (' Axis labels   :',line(8:i-2))
c
   70 read (iunit,6000,err=9000,end=9000) line
      if (line(1:6) .ne. ' min={') goto 70
      i = length(line)
      read (line(7:i-2),*,err=9200) (xmin(j),j=1,3)
      call rvalut (' Origin is at  :',3,xmin)
c
   80 read (iunit,6000,err=9000,end=9000) line
      if (line(1:6) .ne. ' max={') goto 80
      i = length(line)
      read (line(7:i-2),*,err=9200) (xmax(j),j=1,3)
      call rvalut (' Top is at     :',3,xmax)
c
      k = 3
      do i=1,3
        if (xmin(i).eq.xmax(i)) k = k - 1
      end do
      call ivalut (' Dimension of the matrix (1,2,3) :',1,k)
c
      if (k .lt. 2) then
        call errcon ('Not a 2D or 3D matrix')
        goto 9100
      end if
c
c ... guess nr of grid points by assuming an equal nr of points
c     along each axis
c
c ... better would be to assume identical spacings, D, along
c     each axis; then solve:
c     D^3 * (1-N) + D^2 * (dx+dy+dz) + D * (dx*dy + dx*dz + dy*dz) +
c     dx*dy*dz = 0
c     where dx = xmax - xmin etc.
c     this comes from: (1+dx/D)*(1+dy/D)*(1+dz/D) = N
c     if dx or dy or dz = 0 the equation becomes (assume dz = 0):
c     D^3(1-N) + D^2(dx+dy) + D(dx*dy) = 0
c     ==> D=0 (nonsense) OR D^2 * (1-N) + D * (dx+dy) + dx*dy = 0
c     solve for D, then NX=1+dx/D etc. and N must be NX*NY*NZ
c     I'll implement this some day when I have nothing better to do
c
      if (k .eq. 3) then
        j = nint ( float(inum)**(1.0/3.0) )
        nnum = 1
        do i=1,3
          extent (i,imap) = j
          nnum = nnum * j
        end do
      else
        j = nint ( sqrt(float(inum)) )
        nnum = 1
        do i=1,3
          if (xmin(i).eq.xmax(i)) then
            extent (i,imap) = 1
          else
            extent (i,imap) = j
            nnum = nnum * j
          end if
        end do
      end if
c
  100 continue
      call ivalut (' Try grid extent :',3,extent(1,imap))
      call jvalut (' Nr of points in such a grid   :',1,nnum)
      call jvalut (' Nr of points actually in file :',1,inum) 
      if (inum .ne. nnum) then
        call prompt (' Enter 0 0 0 to abort or')
        call ivalin (' Enter nr of grid points :',3,extent(1,imap))
        nnum = 1
        do i=1,3
          if (extent(i,imap) .lt. 1) then
            call prompt (' ... Giving up on this grid ...')
            goto 9100
          end if
          nnum = nnum * extent(i,imap)
        end do
        goto 100
      end if
c
c ... try to read with the current grid
c
      call x3dmin (iunit,map(1,imap),extent(1,imap),
     +  extent(2,imap),extent(3,imap),ierr)
      if (ierr .ne. 0) goto 9100
      call prompt (' Matrix read okay')
      write (*,*)
c
c ... figure out dummy cell etc.
c
      do i=1,3
        cell (i,imap) = 360.0
        cell (i+3,imap) = 90.0
      end do
      call prompt (
     +  ' Use real cell ONLY for Fractional Translation Functions !')
      call fvalin (' Cell   ?',6,cell(1,imap))
c
      call rvalut (' Origin :',3,xmin)
      call rvalut (' Top    :',3,xmax)
      reply='N'
      if (cell(1,imap).ne.360.0) reply = 'Y'
      call textin (' Are origin and top in FRACTIONAL coordinates ?',
     +  reply)
      call upcase (reply)
      if (reply .eq. 'N') then
        call orthog (cell(1,imap), a, 1)
        call mulmtx (a, xmin, ymin, 3, 3, 1)
        call mulmtx (a, xmax, ymax, 3, 3, 1)
        do i=1,3
          xmin (i) = ymin (i)
          xmax (i) = ymax (i)
        end do
      end if
c
      do i=1,3
        if (extent(i,imap) .gt. 1) then
          space(i) = (xmax(i)-xmin(i)) / float(extent(i,imap)-1)
        else
          space(i) = 1.0
        end if
        grid (i,imap) = nint (1.0/space(i))
        origin (i,imap) = nint (grid(i,imap)*xmin(i))
      end do
c
      write (*,*)
      call prompt (' Set up cell etc.')
      call rvalut (' Grid spacing :',3,space)
c
      call telmap (grid(1,imap),origin(1,imap),extent(1,imap),
     +             cell(1,imap))
      call prompt (' Spacegroup set to P1')
c
c      incore (imap) = .true.
c      spaceg (imap) = 1
c      uvw (1,imap) = 2
c      uvw (2,imap) = 1
c      uvw (3,imap) = 3
c
c      npnt (imap) = extent(1,imap)*extent(2,imap)*extent(3,imap)
c      call xstats (map(1,imap),npnt(imap),edstat(3,imap),
c     +  edstat(4,imap),edstat(1,imap),edstat(2,imap),xdum)
c      edstat(5,imap) = edstat(4,imap)*edstat(4,imap)
c      call rvalut (' Sum of density in map :',1,xdum)
c
      ierr = 0
      close (iunit)
      return
c
 9000 continue
      call errcon ('While reading file')
c
 9100 continue
      ierr = -3
      close (iunit)
      return
c
 9200 continue
      call errcon ('While reading a number')
      call textut (' Line >',line)
      goto 9100
c
      end
c
c
c
      subroutine distri (imap,filenm,nbins,maxmap,maxpnt,map)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer maxbin
      parameter (maxbin = 999)
c
      real dmin,dmax,step
c
      integer imap,i,nbins,iunt,cnt(maxbin),ierr,nxyz,nmax,leng1
c
      logical xinter
c
      character filenm*(*),myline*256
c
code ...
c
      if (nbins .lt. 3 .or. nbins .gt. maxbin) then
        call errcon ('Invalid number of bins')
        call jvalut (' Must be >= 3 and <= :',1,maxbin)
        return
      end if
c
      iunt = 87
      close (iunt)
      call xopxua (iunt,filenm,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call stamp (myline)
      myline = 'REMARK '//myline
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Density histogram of ',name(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK File ',file(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Comment ',commnt(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a,i3)',err=9000) 'NPOINT ',nbins
      write (iunt,'(a,i3)',err=9000) 'COLOUR ',4
c
      write (iunt,'(a)',err=9000)
     +  'XLABEL Electron density intervals'
      write (iunt,'(a)',err=9000)
     +  'YLABEL Number of map points'
c
      dmin = edstat(1,imap)
      dmax = edstat(2,imap)
      step = (dmax-dmin)/float(nbins-1)
c
      write (myline,'(a,1p,2e15.6)',err=9000)
     +  'REMARK Density minimum, maximum ',dmin,dmax
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,i15,1p,2e15.6)',err=9000)
     +  'REMARK Nr of bins, bin size ',nbins,step
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a,1p,2e15.6)',err=9000) 'XLIMIT ',dmin,step
c
c ... get histogram
c
      nxyz = extent(1,imap) * extent(2,imap) *
     +       extent(3,imap)
      call prompt (' Calculating histogram ...')
c
      call gethis (map(1,imap),nxyz,dmin,step,cnt,nbins,ierr)
      if (ierr .ne .0) then
        close (iunt)
        return
      end if
c
      nmax = 0
      do i=1,nbins
        nmax = max (nmax,cnt(i))
      end do
c
      write (iunt,'(a,1p,4e15.6)',err=9000) 'XYVIEW ',
     +  (dmin-0.01*step),(dmax+0.01*step),0.0,float(nmax+1)
c
      write (iunt,'(a)',err=9000) 'YVALUE *'
c
      write (iunt,'(5i15)',err=9000) (cnt(i),i=1,nbins)
c
      write (iunt,'(a)',err=9000) 'END'
c
      call prompt (' O2D file written')
      close (iunt)
      return
c
 9000 continue
      call errcon ('While writing O2D file')
      close (iunt)
      return
c
      end
c
c
c
      subroutine rfurey (imap,ifmt,ierr,maxmap,maxpnt,map)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      integer maxlen
      parameter (maxlen=4096)
c
      real xdum(maxlen)
c
      integer imap,ifmt,ierr,i,iunit,ix,iy,iz
c
      logical xinter
c
      byte mybyte(maxlen)
c
      character line*120
c
code ...
c
      ierr = -1
      iunit = 1
c
      if (ifmt .eq. 14) then
        call textut (' Opening PHASES map :',file(imap))
      else
        call textut (' Opening PHASES mask :',file(imap))
      end if
c
      call xopxob (iunit,file(imap),xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening file')
        ierr = -2
        return
      end if
c
      read (iunit,end=9000,err=9000) line
c
      read (iunit,end=9000,err=9000)
     +  (cell(i,imap),i=1,6),(grid(i,imap),i=1,3),
     +  (origin(1,imap),i=1,3),ix,iy,iz
      extent (1,imap) = ix - origin(1,imap) + 1
      extent (2,imap) = iy - origin(2,imap) + 1
      extent (3,imap) = iz - origin(3,imap) + 1
c
      call telmap (grid(1,imap),origin(1,imap),extent(1,imap),
     +  cell(1,imap))
c
      if (extent(1,imap) .gt. maxlen) then
        ierr = -3
        call errcon ('Records too long')
        call jvalut (' Max record length :',1,maxlen)
        return
      end if
c
      do i=1,3
        if (cell(i,imap) .lt. 1.0 .or. cell(i,imap) .gt. 1000000. .or.
     +      cell(i+3,imap) .lt. 10. .or. cell(i+3,imap) .gt. 170.) then
          ierr = -4
          call errcon ('Strange cell')
          return
        end if
      end do
c
      if (ifmt .eq. 15) then
        call rfur2 (map(1,imap),iunit,extent(1,imap),
     +              extent(2,imap),extent(3,imap),mybyte,maxlen,ierr)
      else
        call rfur1 (map(1,imap),iunit,extent(1,imap),
     +              extent(2,imap),extent(3,imap),xdum,maxlen,ierr)
      end if
c
      if (ierr .ne. 0) goto 9000
c
      call prompt (' File read OK')
      call prompt (' Spacegroup set to P1')
      close (iunit)
      ierr = 0
      return
c
 9000 continue
      close (iunit)
      ierr = -5
      call errcon ('While reading file header')
      return
c
      end
c
c
c
      subroutine proj1d (imap,filenm,itype,
     +  maxmap,maxpnt,map,buffer)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap),buffer(maxpnt)
c
      real ave,sdv,xmin,xmax,xtot,dx
c
      integer imap,itype,iunt,ierr,j1,j2,i,leng1
c
      logical xinter
c
      character filenm*(*),myline*256,plname(3)*1
c
      data plname /'X','Y','Z'/
c
code ...
c
      if (itype .lt. 1 .or. itype .gt. 3) then
        call errcon ('Invalid projection axis')
        return
      end if
c
      j1 = origin(itype,imap)
      j2 = (origin(itype,imap)+extent(itype,imap)-1)
c
      call dopro1 (map(1,imap),extent(1,imap),extent(2,imap),
     +  extent(3,imap),buffer,itype) 
c
      call xstats (buffer,extent(itype,imap),
     +  ave,sdv,xmin,xmax,xtot)
      dx = 0.01*(xmax-xmin)
c
      write (*,6000) extent(itype,imap),ave,sdv,xmin,xmax
 6000 format (' Nr of points in projection : ',i12,1p/
     +        ' Average density            : ',e12.4/
     +        ' Standard deviation         : ',e12.4/
     +        ' Minimum density            : ',e12.4/
     +        ' Maximum density            : ',e12.4)
c
      iunt = 87
      close (iunt)
      call xopxua (iunt,filenm,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call stamp (myline)
      myline = 'REMARK '//myline
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a)',err=9000) 'COLOUR 4'
c
      write (myline,*,err=9000)
     +  'REMARK 1D projection onto ',plname(itype),'-axis'
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Map ',name(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK File ',file(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Comment ',commnt(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,6f10.2)',err=9000)
     +  'REMARK Cell ',(cell(i,imap),i=1,6)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,3i8)',err=9000)
     +  'REMARK Grid ',(grid(i,imap),i=1,3)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(3a)',err=9000)
     +  'XLABEL ',plname(itype),'-axis'
c
      write (iunt,'(a)',err=9000)
     +  'YLABEL Integrated (projected) density'
c
      write (iunt,'(a,i8)',err=9000) 'NPOINT',extent(itype,imap)
c
      write (iunt,'(a,1p,4e15.4)',err=9000) 'XYVIEW',
     +  float(j1-1),float(j2+1),xmin-dx,xmax+dx
c
      write (myline,*,err=9000) 'XLIMIT ',j1,1
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a)',err=9000) 'YVALUE * '
c
      write (iunt,'(1p,6e12.4)',err=9000)
     +  (buffer(i),i=1,extent(itype,imap))
c
      write (iunt,'(a)',err=9000) 'END'
c
      call prompt (' O2D file written')
      close (iunt)
      return
c
 9000 continue
      call errcon ('While writing O2D file')
      close (iunt)
      return
c
      end
c
c
c
      subroutine proj2d (imap,filenm,itype,iplane,jplane,
     +  maxmap,maxpnt,map,buffer)
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap),buffer(maxpnt)
c
      integer nlevel
      parameter (nlevel = 8)
c
      real ave,sdv,xmin,xmax,dx,xtot
c
      integer imap,itype,iplane,iunt,ierr,j,i1,i2,i,jplane,leng1
c
      logical xinter
c
      character filenm*(*),myline*256,plname(3)*1
c
      data plname /'X','Y','Z'/
c
code ...
c
      if (itype .lt. 1 .or. itype .gt. 3) then
        call errcon ('Invalid plane type')
        return
      end if
c
      call ilohi (iplane,jplane)
      i1 = origin(itype,imap)
      i2 = (origin(itype,imap)+extent(itype,imap)-1)
      if (iplane .lt. i1 .or.
     +    iplane .gt. i2) then
        call errcon ('Plane index out of range')
        call ivalut (' Requested   :',1,iplane)
        call ivalut (' Lower limit :',1,i1)
        call ivalut (' Upper limit :',1,i2)
        return
      end if
      if (jplane .lt. i1 .or.
     +    jplane .gt. i2) then
        call errcon ('Plane index out of range')
        call ivalut (' Requested   :',1,jplane)
        call ivalut (' Lower limit :',1,i1)
        call ivalut (' Upper limit :',1,i2)
        return
      end if
c
c ... account for offset !!!
c
      iplane = iplane - origin(itype,imap) + 1
      jplane = jplane - origin(itype,imap) + 1
c
      call dopro2 (map(1,imap),extent(1,imap),extent(2,imap),
     +  extent(3,imap),buffer,itype,iplane,jplane) 
c
      if (itype .eq. 1) then
        i1 = 2
        i2 = 3
      else if (itype .eq. 2) then
        i1 = 1
        i2 = 3
      else if (itype .eq. 3) then
        i1 = 1
        i2 = 2
      end if
c
      call xstats (buffer,extent(i1,imap)*extent(i1,imap),
     +  ave,sdv,xmin,xmax,xtot)
      dx = 0.01*(xmax-xmin)
c
      write (*,6000) extent(i1,imap)*extent(i1,imap),
     +  ave,sdv,xmin,xmax
 6000 format (' Nr of points in projection : ',i12,1p/
     +        ' Average density            : ',e12.4/
     +        ' Standard deviation         : ',e12.4/
     +        ' Minimum density            : ',e12.4/
     +        ' Maximum density            : ',e12.4)
c
      iunt = 87
      close (iunt)
      call xopxua (iunt,filenm,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call stamp (myline)
      myline = 'REMARK '//myline
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a,i3)',err=9000) 'NLEVEL ',nlevel
      write (iunt,'(a)',err=9000) 'LEVELS'
      write (myline,'(1p,20(e12.4,1x))',err=9000)
     +  (ave+float(j+1)*0.5*sdv,j=1,nlevel)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a)',err=9000) 'COLOUR'
      write (myline,'(20i5)',err=9000) ( (j+1)/2 ,j=1,nlevel)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000)
     +  'REMARK Contour plot of ',plname(i1),'-',
     +  plname(i2),' projected from plane nr ',
     +  (iplane + origin(itype,imap) - 1),' to ',
     +  (jplane + origin(itype,imap) - 1)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Map ',name(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK File ',file(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(9a)',err=9000)
     +  'REMARK Comment ',commnt(imap)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,6f10.2)',err=9000)
     +  'REMARK Cell ',(cell(i,imap),i=1,6)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,'(a,3i8)',err=9000)
     +  'REMARK Grid ',(grid(i,imap),i=1,3)
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000)
     +  'XLABEL ',plname(i1),'-axis'
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000)
     +  'YLABEL ',plname(i2),'-axis'
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (iunt,'(a,i3)',err=9000) 'XPOINT ',extent(i1,imap)
      write (iunt,'(a,i3)',err=9000) 'YPOINT ',extent(i2,imap)
c
      write (myline,*,err=9000) 'XLIMIT ',
     +  float(origin(i1,imap))/float(grid(i1,imap)),
     +  float(origin(i1,imap)+extent(i1,imap)-1)/float(grid(i1,imap))
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      write (myline,*,err=9000) 'YLIMIT ',
     +  float(origin(i2,imap))/float(grid(i2,imap)),
     +  float(origin(i2,imap)+extent(i2,imap)-1)/float(grid(i2,imap))
      call pretty (myline)
      write (iunt,'(a)',err=9000) myline(1:leng1(myline))
c
      call gjk001 (buffer,1,2,3,1,extent(i1,imap),
     +  extent(i2,imap),1,iunt,ierr)
      if (ierr .ne. 0) goto 9000
c
      write (iunt,'(a)',err=9000) 'END'
c
      call prompt (' O2D file written')
      close (iunt)
      return
c
 9000 continue
      call errcon ('While writing O2D file')
      close (iunt)
      return
c
      end
c
c
c
      subroutine mstats (imap,jmap,cutoff,buffer,maxpnt,maxmap,map)
c
c ... calculate statistics inside/outside mask
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
      real buffer(maxpnt)
c
      real cutoff
c
      integer imap,jmap,ierr,ilim(2,3),jlim(2,3)
c
code ...
c
      ierr = -1
c
      call comap (imap,jmap,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call getmst (
     +  map(1,imap),extent(1,imap),extent(2,imap),extent(3,imap),
     +  map(1,jmap),extent(1,jmap),extent(2,jmap),extent(3,jmap),
     +  ilim,jlim,cutoff,buffer,maxpnt)
c
      return
      end
c
c
c
      subroutine dopeek (imap,option,infile,outfil,mode,par,
     +                   maxpnt,maxmap,map)
c
c ... peek density for atoms in a PDB file
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      real rad
c
      integer imap,iunit,junit,ierr,npt
c
      logical xinter
c
      character*(*) option,infile,outfil,mode,par
      character mymode*1
c
code ...
c
      iunit = 30
      junit = 31
c
      call xopxoa (iunit,infile,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening input PDB file')
        return
      end if
c
      call xopxua (junit,outfil,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening output PDB file')
        close (iunit)
        return
      end if
c
      mymode = mode(1:1)
      if (option(1:2) .eq. 'VA') then
        if (mymode .ne. 'N' .and. mymode .ne. 'I' .and.
     +      mymode .ne. 'S') mymode = 'S'
        npt = 0
        rad = 0.0
      else if (option(1:2) .eq. 'CU') then
        if (mymode .ne. 'I' .and. mymode .ne. 'A' .and.
     +      mymode .ne. 'R' .and. mymode .ne. 'M') mymode = 'R'
        call str2i (par,npt,ierr)
        if (ierr .ne. 0) npt = 2
        if (npt .lt. 1) then
          call errcon ('Nr of points reset to 1')
          npt = 1
        end if
        rad = 0.0
      else if (option(1:2) .eq. 'SP') then
        if (mymode .ne. 'I' .and. mymode .ne. 'A' .and.
     +      mymode .ne. 'R' .and. mymode .ne. 'M') mymode = 'R'
        npt = 0
        call str2r (par,rad,ierr)
        if (ierr .ne. 0) rad = 3.0
        if (rad .lt. 0.1) then
          call errcon ('Radius reset to 0.1 A')
          rad = 0.1
        end if
      end if
c
      call peekit (map(1,imap),grid(1,imap),origin(1,imap),
     +  extent(1,imap),extent(2,imap),extent(3,imap),
     +  cell(1,imap),
     +  option,iunit,junit,mymode,npt,rad)
c
      close (iunit)
      close (junit)
c
      return
      end
c
c
c
      subroutine rdem08 (imap,ierr,maxmap,maxpnt,map)
c
c ... read KI-Stockholm-style EM map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      real pixel
c
      integer imap,ierr,i,j,iunit
c
      logical xinter
c
code ...
c
      ierr = -1
      iunit = 1
c
      call textut (' Opening EM08 map :',file(imap))
c
      call xopxob (iunit,file(imap),xinter(),ierr)
      if (ierr .ne. 0) then
        ierr = -2
        call errcon ('While opening file')
        return
      end if
c
      do i=1,3
        extent(i,imap) = 100
      end do
c
      call jvalin (' Points along X, Y, Z ?',3,extent(1,imap))
c
      do i=1,3
        if (extent(i,imap) .lt. 1) then
          call errcon ('Invalid grid')
          ierr = -3
          return
        end if
      end do
c
      pixel = 7.58
      call fvalin (' Pixel size (A) ?',1,pixel)
      if (pixel .lt. 0.001) then
        call errcon ('Invalid pixel size')
        ierr = -4
        return
      end if
c
      j = 1
      do i=1,3
        j = j * extent(i,imap)
        cell (i,imap) = 100.0*pixel
        cell (i+3,imap) = 90.0
        grid (i,imap) = 100
        origin (i,imap) = 0
      end do
c
      call telmap (grid(1,imap),origin(1,imap),extent(1,imap),
     +  cell(1,imap))
c
      call jvalut (' Points in map :',1,j)
      if (j .gt. maxpnt) then
        call errcon ('Too many points in map')
        ierr = -5
        return
      end if
c
      call rem08 (map(1,imap),iunit,extent(1,imap),
     +            extent(2,imap),extent(3,imap),ierr)
c
      if (ierr .ne. 0) goto 9000
c
      call prompt (' File read OK')
      call prompt (' Spacegroup set to P1')
      close (iunit)
      ierr = 0
c
      return
c
 9000 continue
      close (iunit)
      ierr = -6
      call errcon ('While reading file')
c
      return
      end
c
c
c
      subroutine vbox (imsk,lline)
c
c ... write VRML code to draw box around mask IMSK
c
      include 'mapman.incl'
c
      real a(3,3),xyz(3,8),dum(3)
c
      integer i,imsk,ix,iy,iz
c
      logical lline
c
code ...
c
      call orthog (cell(1,imsk), a, 0)
c
      i = 0
      do ix=0,1
        dum(1) = float (origin(1,imsk)+ix*extent(1,imsk))/
     +           float(grid(1,imsk))
        do iy=0,1
          dum(2) = float (origin(2,imsk)+iy*extent(2,imsk))/
     +             float(grid(2,imsk))
          do iz=0,1
            dum(3) = float (origin(3,imsk)+iz*extent(3,imsk))/
     +               float(grid(3,imsk))
            i = i + 1
            call mulmtx (a,dum,xyz(1,i),3,3,1)
          end do
        end do
      end do
c
      call xvrml_box8 (xyz,lline)
c
      return
      end
c
c
c
      subroutine rdempi (imap,ierr,maxmap,maxpnt,map)
c
c ... read MPI-style EM map
c
      include 'mapman.incl'
c
      integer maxpnt,maxmap
      real map(maxpnt,maxmap)
c
      real*4 emdata(40)
      real pixel(3)
c
      integer imap,ierr,i,j,k,iunit
      integer*4 nb,nx,ny,nz
c
      logical xinter
c
      character coment*80
c
      data pixel /3*7.58/
c
code ...
c
      ierr = -1
      iunit = 1
c
      call textut (' Opening MPI map :',file(imap))
c
      call xopxob (iunit,file(imap),xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening file')
        return
      end if
c
      call prompt (' Reading header ...')
      read (iunit,end=9000,err=9000) nb,nx,ny,nz,
     +  (coment(j:j),j=1,80),
     +  (emdata(k),k=1,40)
c
      extent(1,imap) = nx
      extent(2,imap) = ny
      extent(3,imap) = nz
c
      commnt (imap) = coment
c
      call jvalut (' # X/Y/Z :',3,extent(1,imap))
      call textut (' Comment :',coment)
      call rvalut (' EM data :',40,emdata)
      call jvalut (' Bytes per pixel :',1,nb)
c
      if (nb .ne. 4) then
        call errcon ('Must be 4 bytes per pixel - sorry !')
        ierr = -2
        return
      end if
c
      do i=1,3
        if (extent(i,imap) .lt. 1) then
          call errcon ('Invalid grid')
          ierr = -3
          return
        end if
      end do
c
      call fvalin (' X/Y/Z Pixel size (A) ?',3,pixel)
      do i=1,3
        if (pixel(i) .lt. 0.001) then
          call errcon ('Invalid pixel size')
          ierr = -4
          return
        end if
      end do
c
      j = 1
      do i=1,3
        j = j * extent(i,imap)
        cell (i,imap) = 100.0*pixel(i)
        cell (i+3,imap) = 90.0
        grid (i,imap) = 100
        origin (i,imap) = 0
      end do
c
      call telmap (grid(1,imap),origin(1,imap),extent(1,imap),
     +  cell(1,imap))
c
      call jvalut (' Points in map :',1,j)
      if (j .gt. maxpnt) then
        call errcon ('Too many points in map')
        ierr = -5
        return
      end if
c
      call rempi (map(1,imap),iunit,extent(1,imap),
     +            extent(2,imap),extent(3,imap),ierr)
c
ccc      call rvalut (' DATA :',10,map(1,imap))
c
      if (ierr .ne. 0) goto 9000
c
      call prompt (' File read OK')
      call prompt (' Spacegroup set to P1')
      close (iunit)
      ierr = 0
c
      return
c
 9000 continue
      close (iunit)
      ierr = -6
      call errcon ('While reading file')
c
      return
      end
