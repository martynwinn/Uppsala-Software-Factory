      program mama
c
c ... MAnipulate MAsks
c
c ... Gerard Kleywegt @ 930224
c
      implicit none
c
      character*12 prognm,vers
      parameter (prognm = 'MAMA', vers = '080625/6.1.9')
c
      include 'mama.incl'
c
      integer maxsiz, maxmask
      parameter (maxsiz = maxgk1)
      parameter (maxmask = maxgk3)
c
c      pointer (iaptr,maska)
c      pointer (ibptr,shadow)
c      pointer (icptr,ibuff)
c
c      integer maska(1),shadow(1),ibuff(1)
c      integer malloc
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr
      integer fmalloc
#endif
c
      integer nb, masksize, nummasks, i1, i2
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
      masksize = maxsiz
      call extint ('MASKSIZE',masksize)
      masksize = max ( masksize , minsiz )
      call jvalut (' Allocate masks of size :',1,masksize)
c
      nummasks = 4
      call extint ('NUMMASKS',nummasks)
      nummasks = max ( min ( maxmask, nummasks ), 1 )
      call jvalut (' Max number of masks    :',1,nummasks)
      write (*,*)
c
c ... WRDBYT accounts for 4 or 8 bytes per word
c
   10 continue
      nb = wrdbyt*masksize*nummasks
      iaptr = fmalloc (nb)
      nb = wrdbyt*masksize
      ibptr = fmalloc (nb)
      nb = wrdbyt*masksize
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
c ... BUFFER and IBUFF were equivalenced in the old version of MAMA
c     without dynamic memory allocation (they should never both be used
c     at the same time); passing BUFFER for both has the same net effect
c     (i.e., using the same memory addresses for two different arrays)
c ... ditto for RSHAD and SHADOW
c
      lretry = .false.
      i1 = masksize
      i2 = nummasks
c
      call domama (%val(iaptr),%val(ibptr),%val(ibptr),
     +             %val(icptr),%val(icptr),
     +             masksize, nummasks, lretry, i1, i2)
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
c
      if (lretry) then
        masksize = i1
        nummasks = i2
        masksize = max ( masksize , minsiz )
        nummasks = max ( min ( maxmask, nummasks ), 1 )
        write (*,*)
        call jvalut (' Allocate masks of size :',1,masksize)
        call jvalut (' Max number of masks    :',1,nummasks)
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
      subroutine domama (mask, shadow, rshad, ibuff, buffer,
     +                   maxpnt, maxmsk, lretry, i1, i2)
c
c ... MAnipulate MAsks
c
c ... Gerard Kleywegt @ 930224
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
      integer shadow(maxpnt),ibuff(maxpnt)
      real    buffer(maxpnt),rshad(maxpnt)
c
      real dumncs(12,maxncs),duminv(12,maxncs)
c
      integer maxopt
      parameter (maxopt = 10)
c
      real cog(3)
      real perc,total,user,sys,rdum,cellv,voxv,gridv,maskv
      real xball,yball,zball,rball,vrdist,rr,gg,bb,xdum
      real rgbbg1,rgbbg2,rgbbg3,rgbfg1,rgbfg2,rgbfg3,sphrad
c
      integer ifudge(3)
      integer length,whichm,ivrml,nunmin,nextra,i1,i2,numcyc
      integer i,nopt,ierr,imsk,j,smofac,jmsk,ndum,iunit,kmsk
      integer nrd,munit,itra,idum,leng1,iox,ioy,ioz,ncsdum,jdum
c
      logical xinter,linter,simask,unsave
      logical ldone,linit,lvrml,lline,lretry,lecho
c
      character line*256,pro*8,optpar(maxopt)*80,reply*1
      character opar(4)*40,inimac*128,cavrml*4,vrmlbg*25
      character vrmldc*25,vrfile*128
c
      data ifudge /0,0,0/
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
      write (*,6000) maxmsk,maxpnt,maxsym,maxncs
c
      nmask  = 0
      nsymop = 1
      nrtncs = 0
      numcyc = 1
c
      munit = 5
      iunit = 9
c
      linit = .false.
c
      do i=1,maxmsk
        incore (i) = .false.
        change (i) = .false.
        unt (i) = 19 + i
c
        commnt (i) = 'No comment'
c
c        write (line,'(a,i2,a)') 'm',i,'.mask'
c        call remspa (line)
        file (i) = 'not_defined'
c
        write (line,'(a,i2)') 'M#@$!',i
        call remspa (line)
        name (i) = line
      end do
c
      do i=1,3
c
        do j=1,maxmsk
          cell (i,j)   = 100.0
          cell (i+3,j) =  90.0
          origin (i,j) = 0
          extent (i,j) = 100
          grid (i,j)   = 100
        end do
c
        newcel (i) = 100.0
        newcel (i+3) = 90.0
        newori (i) = 0
        newext (i) = 100
        newgrd (i) = 100
        newpad (i) = 10
c
      end do
c
      newrad = 2.0
      xball = 0.0
      yball = 0.0
      zball = 0.0
      rball = 25.0
c
      nunmin = 5
      nextra = 2
c
      lvrml = .false.
      ivrml = 99
      cavrml = ' CA '
      vrdist = 4.5
      vrmlbg = 'black'
      vrmldc = 'white'
      vrfile = 'mama.wrl'
      rgbbg1 = 0.0
      rgbbg2 = 0.0
      rgbbg2 = 0.0
      rgbfg1 = 1.0
      rgbfg2 = 1.0
      rgbfg2 = 1.0
      call xvrml_init ()
c
      do i=1,12
        newrt (i) = 0.0
        do j=1,maxsym
          symmop (i,j) = 0.0
        end do
        do j=1,maxncs
          rtncs (i,j) = 0.0
        end do
      end do
      newrt (1) = 1.0
      newrt (5) = 1.0
      newrt (9) = 1.0
      do j=1,maxsym
        symmop (1,j) = 1.0
        symmop (5,j) = 1.0
        symmop (9,j) = 1.0
      end do
      do j=1,maxncs
        rtncs (1,j) = 1.0
        rtncs (5,j) = 1.0
        rtncs (9,j) = 1.0
      end do
      do i=1,3
        asuori (i) = 0
        asuext (i) = 100
      end do
      ovmode = 'COUNT'
c
c ... define some symbols
c
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
        pro='$MAMA > '
      else
        pro=' MAMA > '
      end if
c
c --- FORMATS
c
 6000 format (
     +  ' Max nr of masks in memory    : ',i10/
     +  ' Max nr of points in mask     : ',i10/
     +  ' Max nr of symmetry operators : ',i10/
     +  ' Max nr of NCS RT operators   : ',i10/)
c
 6010 format (/
     +  ' MAMA''s options :'//
     +  ' ? (list options)           ! (comment)'/
     +  ' @ macrofile                QUit'/
     +  ' & symbol value             & ? (list symbols)'/
     +  ' $ shell_command            ZP_restart masksize nummasks'/
     +  ' ECho on_off                # parameter(s) (command history)'/
     + /
     +  ' REad maskname filename     WRite maskname filename [how]'/
     +  ' DElete maskname            LIst [maskname]'/
     +  ' DUplicate newmask oldmask  NBr_count maskname'/
     +  ' UNite newmask mask1 mask2  SImilarity mask1 mask2'/
     +  ' ODl maskname filename      BOrder_check maskname'/
     +  ' DOt_odl maskname filename [radius]'/
     +  ' EZd maskname ezdfile       LAbel mask text'/
     +  ' TRanslate mask tx ty tz    GTranslate gx gy gz'/
     +  ' ATom_fit maskname pdbfile [pdb_inside] [pdb_outside]'/
     +  ' INvert_ncs infile outfile'/
     + /
     +  ' NEw ? (list defaults)      NEw ORigin o1 o2 o3'/
     +  ' NEw CEll a b c al be ga    NEw EXtent e1 e2 e3'/
     +  ' NEw GRid g1 g2 g3          NEw SAme maskname'/
     +  ' NEw RAdius value           NEw RT_operator filename'/
     +  ' NEw FActor value           NEw PAd p1 p2 p3'/
     +  ' NEw SPacing value          NEw REset_rt_operator'/
     +  ' NEw ENcompass oldmask      '/
     + /
     +  ' NEw PDb maskname pdbfile   NEw BOnes maskname bonesfile'/
     +  ' NEw COpy mask oldmask      NEw MAke maskname'/
     +  ' NEw UNit_cell mask oldmask NEw OLdbones maskname bonesfile'/
     +  ' NEw BAll mask x y z rad    NEw CUbe mask x y z extent'/
     + /
     +  ' FIll_voids maskname        ISland_erase maskname'/
     +  ' EXpand maskname [cycles]   SMooth maskname factor [cycles]'/
     +  ' COntract maskname [cycles] CUt maskname factor [cycles]'/
     +  ' NOt maskname               BLob_erase maskname min [max]'/
     + /
     +  ' ANd mask1 mask2            OR mask1 mask2'/
     +  ' BUtnot mask1 mask2         XOr mask1 mask2'/
     + /
     +  ' OVerlap ? (list status)    OVerlap SYmm_op filename'/
     +  ' OVerlap ORigin o1 o2 o3    OVerlap EXtent e1 e2 e3'/
     +  ' OVerlap NCs_rt filename    OVerlap REset_ncs'/
     +  ' OVerlap MOde mode_type     '/
     +  ' OVerlap UNit_cell mask factor trim nextra [ezdf]'/
     +  ' OVerlap TRim mask factor [ezdf] '/
     +  ' OVerlap ERase mask factor [ezdf]'/
     +  ' OVerlap EZd mask ezdfile   '/
     +  ' OVerlap ASu_mask newmask oldmask [f1 f2 f3]'/
     + /
     +  ' VRml SEtup central_atom max_dist backgr_col default_col'/
     +  ' VRml INit [filename]                ',
     +              ' VRml COlour_list'/
     +  ' VRml DOt_surface mask [colour]      ',
     +              ' VRml TRace pdb_file [colour]'/
     +  ' VRml CEll mask [colour] [line_solid]',
     +              ' [x_offset] [y_offset] [z_offset]'/
     +  ' VRml BOx mask [colour] [line_solid] ',
     +              ' VRml CLose_file'/
     +  )
c
 6020 format (/' Mask ',i2,' = ',a/
     +  ' File    = ',a/
     +  ' Grid    = ',3i10/
     +  ' Origin  = ',3i10/
     +  ' Extent  = ',3i10/
     +  ' Cell    = ',6f10.3/
     +  ' Nr of points = ',i10,5x,' Set   = ',i10,' (',f6.2,' %)')
c
 6023 format (
     +  ' Cell volume  = ',1pe10.3,5x,' Voxel = ',e10.3/
     +  ' Grid volume  = ',1pe10.3,5x,' Mask  = ',e10.3)
c
 6025 format (
     +  ' Centre-of-gravity (A) = ',3f10.3/
     +  ' Spacing = ',3f10.3/
     +  ' Top     = ',3i10/
     +  ' Changes = ',l1/
     +  ' Label   = ',a)
c
 6030 format (/' Mask ',i2,' available')
c
 6040 format (' WARNING - unsaved changes to mask ',a,' !!!')
c
 6050 format (/' Current defaults for the next NEW mask:'/
     +  ' Grid    = ',3i10/
     +  ' Origin  = ',3i10/
     +  ' Extent  = ',3i10/
     +  ' Padding = ',3i10/
     +  ' Cell    = ',6f10.3/
     +  ' Radius  = ',f10.3/
     +  ' RT-oper = ',3f12.6/3(11x,3f12.6/),
     +  ' Nr of points = ',i10,' Max = ',i10)
c
 6060 format (/' Current settings for overlap calculations:'/
     +  ' Nr of symmetry operators = ',i10/
     +  ' Nr of NCS RT operators   = ',i10/
     +  ' Asymmetric unit origin   = ',3i10/
     +  ' Asymmetric unit extent   = ',3i10/
     +  ' Nr of points in as. unit = ',i10/
     +  ' Maximum allowed          = ',i10/
     +  ' Overlap mode selected    = ',a)
c
 6070 format (
     +  ' Symmetry operator ',i3,' = ',3f12.6,3(/25x,3f12.6))
 6072 format (
     +  ' NCS RT operator   ',i3,' = ',3f12.6,3(/25x,3f12.6))
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6010)
      write (*,6000) maxmsk,maxpnt,maxsym,maxncs
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
        call gknval ('GKMAMA',inimac,ierr)
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
          call textin (' New MASKSIZE ?',optpar(2))
        end if
        call str2i (optpar(2),idum,ierr)
        if (ierr .ne. 0) goto 10
        i1 = idum
        if (i1 .lt. 10) then
          call errcon ('Silly MASKSIZE')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          write (optpar(3),*) i2
          call textin (' New NUMMASKS ?',optpar(3))
        end if
        call str2i (optpar(3),idum,ierr)
        if (ierr .ne. 0) goto 10
        i2 = idum
        if (i2 .lt. 1) then
          call errcon ('Silly NUMMASKS')
          goto 10
        end if
c
        lretry = .true.
        return
c
c ... OVERLAP
c
      else if (optpar(1)(1:2) .eq. 'OV') then
c
        if (nopt .lt. 2) then
          optpar (2) = '?'
          call textin (' Option ?',optpar(2))
        end if
c
        call upcase (optpar(2))
c
c ... OVERLAP ?
c
        if (optpar(2)(1:1) .eq. '?') then
c
          ndum = asuext(1)*asuext(2)*asuext(3)
c
          write (*,6060) nsymop,nrtncs,(asuori(i),i=1,3),
     +      (asuext(i),i=1,3),ndum,maxpnt,ovmode
c
          if (nsymop .gt. 0) then
            do j=1,nsymop
              write (*,6070) j,(symmop(i,j),i=1,12)
            end do
          end if
c
          if (nrtncs .gt. 0) then
            do j=1,nrtncs
              write (*,6072) j,(rtncs(i,j),i=1,12)
            end do
          end if
c
          if (ndum .gt. maxpnt) then
            call errcon ('Asymmetric unit grid too large')
            call jvalut (' Nr of points :',1,ndum)
            call jvalut (' Maximum      :',1,maxpnt)
          end if          
c
c ... OVERLAP ORIGIN
c
        else if (optpar(2)(1:2) .eq. 'OR') then
c
          do i=1,3
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(i10)') asuori(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for origin ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2i(optpar(i+2),ndum,ierr)
            if (ierr .ne. 0) goto 10
            asuori (i) = ndum
          end do
          call ivalut (' Asymmetric unit origin :',3,asuori)
c
c ... OVERLAP EXTENT
c
        else if (optpar(2)(1:2) .eq. 'EX') then
c
          do i=1,3
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(i10)') asuext(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for origin ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2i(optpar(i+2),ndum,ierr)
            if (ierr .ne. 0) goto 10
            asuext (i) = ndum
          end do
          call ivalut (' Asymmetric unit extent :',3,asuext)
c
          ndum = asuext(1)*asuext(2)*asuext(3)
          if (ndum .gt. maxpnt) then
            call errcon ('Asymmetric unit grid too large')
            call jvalut (' Nr of points :',1,ndum)
            call jvalut (' Maximum      :',1,maxpnt)
          end if
c
c ... OVERLAP SYMMETRY_OPERATORS
c
        else if (optpar(2)(1:2) .eq. 'SY') then
c
          if (nopt .lt. 3) then
            optpar (3) = 'symop.o'
            call textin (' File name ?',optpar(3))
          end if
c
          if (length(optpar(3)) .lt. 1) then
            call errcon ('No file name provided')
            goto 10
          end if
c
          call osymop (iunit,optpar(3),ierr)
          if (ierr .ne. 0) then
            call errcon ( 'While opening O datablock file')
            goto 10
          end if
          close (iunit)
c
          call opoodb (iunit,optpar(3),opar(1),opar(2),
     +      nrd,opar(4),ierr)
          if (ierr .ne. 0) then
            call errcon ('While opening RT file')
            if (linter) goto 10
            call errcon ('Unable to open RT file')
            return
          end if
c
          nsymop = nrd / 12
          call jvalut (' Nr of symmetry operators :',1,nsymop)
          if (nsymop .lt. 1 .or. nsymop .gt. maxsym) then
            call errcon ('Invalid nr of symmetry operators')
            goto 10
          end if
c
          if ((12*nsymop) .ne. nrd) then
            call errcon ('Invalid nr of elements; must be N*12')
            goto 10
          end if
c
          do j=1,nsymop
            read (iunit,opar(4)) (symmop(i,j),i=1,12)
            call fratra (symmop(10,j))
          end do
c
          close (iunit)
c
          call anasgs (nsymop,symmop,.true.,ierr)
          if (ierr .ne. 0) then
            call errcon ('In spacegroup symmetry operators')
            nsymop = 0
          end if
c
c ... OVERLAP RESET_NCS
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          nrtncs = 0
          call jvalut (' Nr of NCS RT operators :',1,nrtncs)
c
c ... OVERLAP MODE
c
        else if (optpar(2)(1:2) .eq. 'MO') then
c
          if (nopt .lt. 3) then
            optpar (3) = ovmode
            call textin (' Mode (Count/Label) ?',optpar(3))
          end if
c
          if (length(optpar(3)) .lt. 1) then
            call errcon ('No overlap mode provided')
            goto 10
          end if
c
          call upcase (optpar(3))
          if (optpar(3)(1:1) .eq. 'C') then
            ovmode = 'COUNT'
          else if (optpar(3)(1:1) .eq. 'L') then
            ovmode = 'LABEL'
          else
            call errcon ('Invalid mode :'//optpar(3))
          end if
c
c ... OVERLAP NCS_RT
c
        else if (optpar(2)(1:2) .eq. 'NC') then
c
          if (nrtncs .ge. maxncs) then
            call errcon ('Too many NCS operators')
            call jvalut (' Maximum :',1,maxncs)
            goto 10
          end if
c
          if (nopt .lt. 3) then
            optpar (3) = ' '
            call textin (' File name ?',optpar(3))
          end if
c
          if (length(optpar(3)) .lt. 1) then
            call errcon ('No file name provided')
            goto 10
          end if
c
          call rdoncs (iunit,optpar(3),nrtncs,maxncs,rtncs,ierr)
c
          if (ierr .ne. 0) then
            call errcon ('While reading NCS operator file')
            if (linter) goto 10
            call errcon ('While reading NCS operator file')
            return
          end if
c
c ... OVERLAP TRIM
c
        else if (optpar(2)(1:2) .eq. 'TR') then
c
          if (nopt .lt. 3) call textin (' Which mask ?',optpar(3))
c
          if (nopt .lt. 4) then
            write (optpar(4),*) float(nrtncs+1)
            call remspa (optpar(4))
            call textin (' Factor ?',optpar(4))
          end if
c
          if (nopt .lt. 5) optpar(5) = ' '
c
          jmsk = whichm(optpar(3),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          call str2r (optpar(4),rdum,ierr)
          if (ierr .ne. 0) goto 10
          if (rdum .le. float(nrtncs)) then
            call errcon ('Factor too low')
            call ivalut (' Should be greater than :',1,nrtncs)
            goto 10
          end if
c
          call textut (' Trim Mask :',name(jmsk))
          call countm (jmsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set before :',1,nset(jmsk))
c
          call overlp ('TRim',jmsk,ierr,1,rdum,optpar(5),
     +      maxpnt,maxmsk,mask,shadow,ibuff,rshad,buffer,ifudge)
c
          if (ierr .ne. 0) then
            call errcon ('During overlap trimming; mask not trimmed')
            goto 10
          end if
c
          call countm (jmsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set after :',1,nset(jmsk))
          change (jmsk) = .true.
c
c ... OVERLAP UNIT_CELL
c
        else if (optpar(2)(1:2) .eq. 'UN') then
c
          if (nopt .lt. 3) call textin (' Which mask ?',optpar(3))
c
          if (nopt .lt. 4) then
            write (optpar(4),*) float(nrtncs+1)
            call remspa (optpar(4))
            call textin (' Overlap factor ?',optpar(4))
          end if
c
          if (nopt .lt. 5) then
            write (optpar(5),*) nunmin
            call remspa (optpar(5))
            call textin (' Trim cut-off ?',optpar(5))
          end if
c
          if (nopt .lt. 6) then
            write (optpar(6),*) nextra
            call remspa (optpar(6))
            call textin (' Extra points ?',optpar(6))
          end if
c
          if (nopt .lt. 7) optpar(7) = ' '
c
          jmsk = whichm(optpar(3),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          call str2r (optpar(4),rdum,ierr)
          if (ierr .ne. 0) goto 10
          if (rdum .le. float(nrtncs)) then
            call errcon ('Factor too low')
            call ivalut (' Should be greater than :',1,nrtncs)
            goto 10
          end if
c
          call str2i (optpar(5),idum,ierr)
          if (ierr .ne. 0) goto 10
          nunmin = max (0, min (26, idum))
c
          call str2i (optpar(6),idum,ierr)
          if (ierr .ne. 0) goto 10
          nextra = max (0, min (3, nextra))
c
          call textut (
     +      ' Unit-cell-based overlap removal for mask :',name(jmsk))
          call countm (jmsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set before :',1,nset(jmsk))
c
          call overun (jmsk,ierr,nunmin,rdum,nextra,optpar(7),
     +      maxpnt,maxmsk,mask,shadow,ibuff,rshad,buffer)
c
          if (ierr .ne. 0) then
            call errcon ('During overlap removal; mask not changed')
            goto 10
          end if
c
          call countm (jmsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set after :',1,nset(jmsk))
          change (jmsk) = .true.
c
c ... OVERLAP ERASE
c
        else if (optpar(2)(1:2) .eq. 'ER') then
c
          if (nopt .lt. 3) call textin (' Which mask ?',optpar(3))
c
          if (nopt .lt. 4) then
            write (optpar(4),*) float(nrtncs+1)
            call remspa (optpar(4))
            call textin (' Factor ?',optpar(4))
          end if
          if (nopt .lt. 5) optpar(5) = ' '
c
          jmsk = whichm(optpar(3),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          call str2r (optpar(4),rdum,ierr)
          if (ierr .ne. 0) goto 10
          if (rdum .le. float(nrtncs)) then
            call errcon ('Factor too low')
            call ivalut (' Should be greater than :',1,nrtncs)
            goto 10
          end if
c
          call textut (' Erase Mask :',name(jmsk))
          call countm (jmsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set before :',1,nset(jmsk))
c
          call overlp ('ERase',jmsk,ierr,0,rdum,optpar(5),
     +      maxpnt,maxmsk,mask,shadow,ibuff,rshad,buffer,ifudge)
c
          if (ierr .ne. 0) then
            call errcon ('During overlap trimming; mask not trimmed')
            goto 10
          end if
c
          call countm (jmsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set after :',1,nset(jmsk))
          change (jmsk) = .true.
c
c ... OVERLAP EZD
c
        else if (optpar(2)(1:2) .eq. 'EZ') then
c
          if (nopt .lt. 3) call textin (' Which mask ?',optpar(3))
          if (nopt .lt. 4) call textin (' EZD file ?',optpar(4))
          jmsk = whichm(optpar(3),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          if (length(optpar(4)) .lt. 1) then
            call errcon ('Invalid EZD file name')
            goto 10
          end if
c
          call textut (' EZD Mask :',name(jmsk))
c
          call overlp ('EZd',jmsk,ierr,-1,float(nrtncs+1),optpar(4),
     +      maxpnt,maxmsk,mask,shadow,ibuff,rshad,buffer,ifudge)
c
c ... OVERLAP ASU_MASK
c
        else if (optpar(2)(1:2) .eq. 'AS') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) call textin (' Old mask ?',optpar(4))
          jmsk = whichm(optpar(4),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          if (imsk .eq. jmsk) then
            call errcon ('Mask1 and Mask2 are identical')
            goto 10
          end if
c
          idum = 0
          if (nopt .ge. 5) then
            call str2i (optpar(5),idum,ierr)
            if (ierr .ne. 0) goto 10
            if (idum .ne. 1) idum = 0
          end if
          ifudge (1) = idum
c
          idum = 0
          if (nopt .ge. 6) then
            call str2i (optpar(6),idum,ierr)
            if (ierr .ne. 0) goto 10
            if (idum .ne. 1) idum = 0
          end if
          ifudge (2) = idum
c
          idum = 0
          if (nopt .ge. 7) then
            call str2i (optpar(7),idum,ierr)
            if (ierr .ne. 0) goto 10
            if (idum .ne. 1) idum = 0
          end if
          ifudge (3) = idum
c
c ... set mask parameters (EXTENT will be set by OVERLP routine)
c
          do i=1,3
            cell (i,imsk)   = cell(i,jmsk)
            cell (i+3,imsk) = cell(i+3,jmsk)
            grid (i,imsk)   = grid(i,jmsk)
            origin (i,imsk) = asuori (i)
          end do
c
          call textut (' ASU Mask :',name(imsk))
          call ivalut (' Fudge margins :',3,ifudge)
c
          call overlp ('ASu',jmsk,ierr,imsk,-1.0,optpar(4),
     +      maxpnt,maxmsk,mask,shadow,ibuff,rshad,buffer,ifudge)
c
          if (ierr .ne. 0) then
            call errcon ('In overlap calculation; mask not generated')
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... invalid OVERLAP sub-option
c
        else
          call errcon ('Invalid option for OVERLAP command')
          call textut (' Value :',optpar(2))
        end if
c
c ... NEW
c
      else if (optpar(1)(1:2) .eq. 'NE') then
c
        if (nopt .lt. 2) then
          optpar (2) = '?'
          call textin (' Option ?',optpar(2))
        end if
c
        call upcase (optpar(2))
c
c ... NEW ?
c
        if (optpar(2)(1:1) .eq. '?') then
c
          ndum = newext(1)*newext(2)*newext(3)
          write (*,6050) (newgrd(i),i=1,3),(newori(i),i=1,3),
     +      (newext(i),i=1,3),(newpad(i),i=1,3),
     +      (newcel(i),i=1,6),newrad,
     +      (newrt(i),i=1,12),ndum,maxpnt
          if (ndum .gt. maxpnt)
     +      call errcon ('NEW mask is too big !')
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
c ... NEW PAD
c
        else if (optpar(2)(1:2) .eq. 'PA') then
c
          do i=1,3
            if (nopt .lt. (i+2)) then
              write (optpar(i+2),'(i10)') newpad(i)
              call remspa (optpar(i+2))
              write (line,'(a,i1,a)') ' Value for padding ',i,' ?'
              call textin (line,optpar(i+2))
            end if
            call str2i(optpar(i+2),ndum,ierr)
            if (ierr .ne. 0) goto 10
            newpad (i) = max (1, ndum)
          end do
          call ivalut (' NEW padding :',3,newpad)
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
c ... NEW SPACING
c
        else if (optpar(2)(1:2) .eq. 'SP') then
c
          if (nopt .lt. 3) then
            write (optpar(3),'(f10.2)') 1.0
            call remspa (optpar(3))
            call textin (' Value for spacing ?',optpar(3))
          end if
          call str2r(optpar(3),rdum,ierr)
          if (ierr .ne. 0) goto 10
          rdum = max (0.1, rdum)
          do i=1,3
            newgrd (i) = newcel(i)/rdum
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
c ... NEW RADIUS
c
        else if (optpar(2)(1:2) .eq. 'RA') then
c
          if (nopt .lt. (3)) then
            write (optpar(3),'(f15.3)') newrad
            call remspa (optpar(3))
            call textin (' Value for radius ?',optpar(3))
          end if
          call str2r(optpar(3),rdum,ierr)
          if (ierr .ne. 0) goto 10
          newrad = rdum
          call fvalut (' NEW radius :',1,newrad)
c
c ... NEW RT_OPERATOR
c
        else if (optpar(2)(1:2) .eq. 'RT') then
c
          if (nopt .lt. 3) then
            optpar (3) = ' '
            call textin (' File name ?',optpar(3))
          end if
c
          if (length(optpar(3)) .lt. 1) then
            call errcon ('No file name provided')
            goto 10
          end if
c
          call opoodb (iunit,optpar(3),opar(1),opar(2),
     +      nrd,opar(4),ierr)
          if (ierr .ne. 0) then
            call errcon ('While opening RT file')
            if (linter) goto 10
            call errcon ('Unable to open RT file')
            return
          end if
c
          if (nrd .ne. 12 .or. ierr .ne. 0) then
            call errcon ('Invalid number of data items')
            goto 10
          end if
c
          read (iunit,opar(4)) (newrt(i),i=1,12)
c
          close (iunit)
c
          call fvalut (' NEW RT-operator :',12,newrt)
c
c ... NEW RESET_RT_OPERATOR
c
        else if (optpar(2)(1:2) .eq. 'RE') then
c
          do i=1,12
            newrt (i) = 0.0
          end do
          newrt (1) = 1.0
          newrt (5) = 1.0
          newrt (9) = 1.0
c
          call fvalut (' NEW RT-operator :',12,newrt)
c
c ... NEW MAKE
c
        else if (optpar(2)(1:2) .eq. 'MA') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          ndum = newext(1)*newext(2)*newext(3)
          if (ndum .gt. maxpnt) then
            call errcon ('NEW mask is too big !')
            call jvalut (' Requested mask size :',1,ndum)
            call jvalut (' Max available size  :',1,maxpnt)
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
            origin (i,imsk) = newori (i)
            extent (i,imsk) = newext (i)
          end do
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call inimsk (mask(1,imsk),extent(1,imsk),
     +      extent(2,imsk),extent(3,imsk),0)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW BALL
c
        else if (optpar(2)(1:2) .eq. 'BA') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) then
            write (optpar(4),*) xball
            call textin (' X-coordinate ?',optpar(4))
          end if
          call str2r (optpar(4),xball,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            write (optpar(5),*) yball
            call textin (' Y-coordinate ?',optpar(5))
          end if
          call str2r (optpar(5),yball,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            write (optpar(6),*) zball
            call textin (' Z-coordinate ?',optpar(6))
          end if
          call str2r (optpar(6),zball,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 7) then
            write (optpar(7),*) rball
            call textin (' Radius ?',optpar(7))
          end if
          call str2r (optpar(7),rball,ierr)
          if (ierr .ne. 0) goto 10
          if (rball .le. 0.0) then
            call errcon ('Radius must be positive')
            goto 10
          end if
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
          end do
c
          call balmsk (cell(1,imsk),grid(1,imsk),
     +      newpad,newrt,maxpnt,xball,yball,zball,rball,
     +      mask(1,imsk),origin(1,imsk),extent(1,imsk),ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW CUBE
c
        else if (optpar(2)(1:2) .eq. 'CU') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) then
            write (optpar(4),*) xball
            call textin (' X-coordinate centre ?',optpar(4))
          end if
          call str2r (optpar(4),xball,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 5) then
            write (optpar(5),*) yball
            call textin (' Y-coordinate centre ?',optpar(5))
          end if
          call str2r (optpar(5),yball,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 6) then
            write (optpar(6),*) zball
            call textin (' Z-coordinate centre ?',optpar(6))
          end if
          call str2r (optpar(6),zball,ierr)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 7) then
            write (optpar(7),*) rball
            call textin (' Extent ?',optpar(7))
          end if
          call str2r (optpar(7),rball,ierr)
          if (ierr .ne. 0) goto 10
          if (rball .le. 0.0) then
            call errcon ('Extent must be positive')
            goto 10
          end if
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
          end do
c
          call cubmsk (cell(1,imsk),grid(1,imsk),
     +      newpad,newrt,maxpnt,xball,yball,zball,rball,
     +      mask(1,imsk),origin(1,imsk),extent(1,imsk),ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW PDB
c
        else if (optpar(2)(1:2) .eq. 'PD') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) then
            optpar (4) = ' '
            call textin (' PDB file ?',optpar(4))
          end if
c
          if (length(optpar(4)) .lt. 1) then
            call errcon ('No PDB file provided')
            goto 10
          end if
c
          call xopxoa (iunit,optpar(4),linter,ierr)
          if (ierr .ne. 0) then
            if (linter) goto 10
            call errcon ('Unable to open PDB file')
            return
          end if
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
          end do
c
          call pdbmsk (iunit,cell(1,imsk),grid(1,imsk),
     +      newpad,newrad,newrt,maxpnt,
     +      mask(1,imsk),origin(1,imsk),extent(1,imsk),ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW BONES
c
        else if (optpar(2)(1:2) .eq. 'BO') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) then
            optpar (4) = ' '
            call textin (' BONES file ?',optpar(4))
          end if
c
          if (length(optpar(4)) .lt. 1) then
            call errcon ('No BONES file provided')
            goto 10
          end if
c
          call xopxoa (iunit,optpar(4),linter,ierr)
          if (ierr .ne. 0) then
            if (linter) goto 10
            call errcon ('Unable to open BONES file')
            return
          end if
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
          end do
c
          call bonesm (iunit,cell(1,imsk),grid(1,imsk),
     +      newpad,newrad,newrt,maxpnt,
     +      mask(1,imsk),origin(1,imsk),extent(1,imsk),ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW OLDBONES
c
        else if (optpar(2)(1:2) .eq. 'OL') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) then
            optpar (4) = ' '
            call textin (' OLD BONES file ?',optpar(4))
          end if
c
          if (length(optpar(4)) .lt. 1) then
            call errcon ('No OLD BONES file provided')
            goto 10
          end if
c
          call xopxoa (iunit,optpar(4),linter,ierr)
          if (ierr .ne. 0) then
            if (linter) goto 10
            call errcon ('Unable to open OLD BONES file')
            return
          end if
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
          end do
c
          call boneso (iunit,cell(1,imsk),grid(1,imsk),
     +      newpad,newrad,newrt,maxpnt,
     +      mask(1,imsk),origin(1,imsk),extent(1,imsk),ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW COPY
c
        else if (optpar(2)(1:2) .eq. 'CO') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) call textin (' Old mask ?',optpar(4))
          jmsk = whichm(optpar(4),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          if (imsk .eq. jmsk) then
            call errcon ('Mask1 and Mask2 are identical')
            goto 10
          end if
c
          ndum = newext(1)*newext(2)*newext(3)
          if (ndum .gt. maxpnt) then
            call errcon ('NEW mask is too big !')
            goto 10
          end if
c
          do i=1,3
            cell (i,imsk)   = cell(i,jmsk)
            cell (i+3,imsk) = cell(i+3,jmsk)
            grid (i,imsk)   = newgrd (i)
            origin (i,imsk) = newori (i)
            extent (i,imsk) = newext (i)
          end do
c
          call sub06m (
     +      mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +      origin(1,jmsk),grid(1,jmsk),
     +      mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +      origin(1,imsk),grid(1,imsk),cell(1,imsk),ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW UNIT_CELL
c
        else if (optpar(2)(1:2) .eq. 'UN') then
c
          if (nopt .lt. 3) call textin (' New mask ?',optpar(3))
          call allocm (optpar(3),imsk,ierr,maxmsk)
          if (ierr .ne. 0) goto 10
c
          if (nopt .lt. 4) call textin (' Old mask ?',optpar(4))
          jmsk = whichm(optpar(4),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          if (imsk .eq. jmsk) then
            call errcon ('Mask1 and Mask2 are identical')
            goto 10
          end if
c
          do i=1,3
            cell (i,imsk)   = newcel (i)
            cell (i+3,imsk) = newcel (i+3)
            grid (i,imsk)   = newgrd (i)
            origin (i,imsk) = newori (i)
            extent (i,imsk) = newext (i)
          end do
c
          call newunt (
     +      mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +      origin(1,jmsk),grid(1,jmsk),cell(1,jmsk),
     +      mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +      origin(1,imsk),grid(1,imsk),cell(1,imsk),
     +      ibuff,shadow,newrt,newpad,maxpnt,ierr)
c
          if (ierr .ne. 0) then
            incore (imsk) = .false.
            name (imsk) = '!@#$%^&*()'
            goto 10
          end if
c
          incore (imsk) = .true.
          select (imsk) = .false.
          change (imsk) = .true.
          file (imsk)   = 'not_defined'
          npnt (imsk) = extent(1,imsk)*extent(2,imsk)*extent(3,imsk)
c
          call countm (imsk,maxpnt,maxmsk,mask)
          call jvalut (' Nr of points set :',1,nset(imsk))
c
c ... NEW SAME
c
        else if (optpar(2)(1:2) .eq. 'SA') then
c
          if (nopt .lt. 3) call textin (' Old mask ?',optpar(3))
          jmsk = whichm(optpar(3),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          do i=1,3
            newcel (i)   = cell (i,jmsk)
            newcel (i+3) = cell (i+3,jmsk)
            newori (i)   = origin (i,jmsk)
            newgrd (i)   = grid (i,jmsk)
            newext (i)   = extent (i,jmsk)
          end do
c
c ... NEW ENCOMPASS
c
        else if (optpar(2)(1:2) .eq. 'EN') then
c
          if (nopt .lt. 3) call textin (' Old mask ?',optpar(3))
          jmsk = whichm(optpar(3),maxmsk)
c
          if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(jmsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
c
          do i=1,3
            if (abs(newcel(i)-cell(i,jmsk)).gt.0.01) then
              call errcon ('Mask cell differs from NEw CEll')
              goto 10
            end if
            if (abs(newcel(i+3)-cell(i+3,jmsk)).gt.0.01) then
              call errcon ('Mask cell differs from NEw CEll')
              goto 10
            end if
            if (newgrd(i) .ne. grid(i,jmsk)) then
              call errcon ('Mask grid differs from NEw GRid')
              goto 10
            end if
          end do
c
          call jvalut (' Current NEw ORigin  :',3,newori)
          call jvalut (' Current NEw EXtent  :',3,newext)
          call jvalut (' Current mask origin :',3,origin(1,jmsk))
          call jvalut (' Current mask extent :',3,extent(1,jmsk))
c
          do i=1,3
            newcel (i)   = cell (i,jmsk)
            newcel (i+3) = cell (i+3,jmsk)
            newgrd (i)   = grid (i,jmsk)
c
            j = max (newori(i)+newext(i)-1,
     +               origin(i,jmsk)+extent(i,jmsk)-1)
            newori (i)   = min (newori(i),origin(i,jmsk))
            newext (i)   = j - newori(i) + 1
          end do
c
          call jvalut (' New origin :',3,newori)
          call jvalut (' New extent :',3,newext)
c
c ... NEW FACTOR
c
        else if (optpar(2)(1:2) .eq. 'FA') then
c
          if (nopt .lt. 3) then
            optpar(3) = '1.0'
            call textin (' Factor ?',optpar(3))
          end if
          call str2r (optpar(3),rdum,ierr)
          if (ierr .ne. 0) goto 10
c
          if (rdum .lt. 0.01 .or. rdum .gt. 999.9) then
            call errcon ('Factor should be in range 0.01 - 999.9')
            goto 10
          end if
c
          do i=1,3
            newori (i)   = nint (rdum*newori(i))
            newext (i)   = nint (rdum*newext(i))
            newgrd (i)   = nint (rdum*newgrd(i))
          end do
c
c ... invalid NEW sub-option
c
        else
          call errcon ('Invalid option for NEW command')
          call textut (' Value :',optpar(2))
        end if
c
c ... DUPLICATE
c
      else if (optpar(1)(1:2) .eq. 'DU') then
c
        if (nopt .lt. 2) call textin (' New mask ?',optpar(2))
        call allocm (optpar(2),imsk,ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) call textin (' Old mask ?',optpar(3))
        jmsk = whichm(optpar(3),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        call copym (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... GTRANSLATE
c
      else if (optpar(1)(1:2) .eq. 'GT') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        write (*,'(/10(1x,a/))')
     +    'WARNING !!! This option will change the origin of your',
     +    'map by any number of grid points that you specify.',
     +    'This is bound to upset the compatibility between this',
     +    'mask and your model/data. Only use this option if you',
     +    'understand this and know what you are doing ...'
c
        call jvalut (' Old origin :',3,origin(1,imsk))
        call jvalut (' Grid       :',3,grid(1,imsk))
c
        if (nopt .lt. 3) then
          optpar(3) = '0'
          call textin (
     +      ' Nr of GRID POINTs translation in X ?',optpar(3))
        end if
        call str2i (optpar(3),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (1,imsk) = origin(1,imsk) + itra
c
        if (nopt .lt. 4) then
          optpar(4) = '0'
          call textin (
     +      ' Nr of GRID POINTs translation in Y ?',optpar(4))
        end if
        call str2i (optpar(4),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (2,imsk) = origin(2,imsk) + itra
c
        if (nopt .lt. 5) then
          optpar(5) = '0'
          call textin (
     +      ' Nr of GRID POINTs translation in Z ?',optpar(5))
        end if
        call str2i (optpar(5),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (3,imsk) = origin(3,imsk) + itra
c
        call jvalut (' New origin :',3,origin(1,imsk))
c
c ... TRANSLATE
c
      else if (optpar(1)(1:2) .eq. 'TR') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        call jvalut (' Old origin :',3,origin(1,imsk))
        call jvalut (' Grid       :',3,grid(1,imsk))
c
        if (nopt .lt. 3) then
          optpar(3) = '0'
          call textin (
     +      ' Nr of unit cells translation in X ?',optpar(3))
        end if
        call str2i (optpar(3),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (1,imsk) = origin(1,imsk) + itra*grid(1,imsk)
c
        if (nopt .lt. 4) then
          optpar(4) = '0'
          call textin (
     +      ' Nr of unit cells translation in Y ?',optpar(4))
        end if
        call str2i (optpar(4),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (2,imsk) = origin(2,imsk) + itra*grid(2,imsk)
c
        if (nopt .lt. 5) then
          optpar(5) = '0'
          call textin (
     +      ' Nr of unit cells translation in Z ?',optpar(5))
        end if
        call str2i (optpar(5),itra,ierr)
        if (ierr .ne. 0) goto 10
        origin (3,imsk) = origin(3,imsk) + itra*grid(3,imsk)
c
        call jvalut (' New origin :',3,origin(1,imsk))
c
c ... LABEL
c
      else if (optpar(1)(1:2) .eq. 'LA') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        if (nopt .lt. 3) call textin (' Label ?',optpar(3))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        commnt (imsk) = optpar(3)
c
c ... READ
c
      else if (optpar(1)(1:2) .eq. 'RE') then
c
        if (nopt .lt. 2) call textin (' New mask ?',optpar(2))
        call allocm (optpar(2),imsk,ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          optpar (3) = file(imsk)
          call textin (' File name ?',optpar(3))
        end if
        file (imsk) = optpar(3)
        if (length(file(imsk)) .lt. 1) then
          call errcon ('No file name provided')
          goto 10
        end if
c
        call xopxoa (unt(imsk),file(imsk),linter,ierr)
        if (ierr .ne. 0) then
          if (linter) goto 10
          call errcon ('Unable to open input file')
            return
        end if
c
        call maskin (unt(imsk),mask(1,imsk),origin(1,imsk),
     +               extent(1,imsk),grid(1,imsk),cell(1,imsk),
     +               maxpnt,ierr)
c
        if (ierr .ne. 0) then
          incore (imsk) = .false.
          name (imsk) = '!@#$%^&*()'
          goto 10
        end if
c
        close (unt(imsk))
        change (imsk) = .false.
        incore (imsk) = .true.
        commnt (imsk) = 'Read from '//file(imsk)
        npnt (imsk) = extent(1,imsk) * extent(2,imsk) *
     +                extent(3,imsk)
c
        goto 10
c
c ... LIST
c
      else if (optpar(1)(1:2) .eq. 'LI') then
c
        nmask = 0
        do i=1,maxmsk
          if (incore(i)) nmask = nmask + 1
        end do
        call ivalut (' Nr of masks in memory :',1,nmask)
c
        if (nopt .lt. 2)  optpar(2) = '*'
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nmask .gt. 0) then
c
          do i=1,maxmsk
            if (select(i)) then
              if (incore(i)) then
                call centrm (i,maxpnt,maxmsk,mask,cog)
ccc                call countm (i,maxpnt,maxmsk,mask)
                perc = 100.0 * float(nset(i)) / float(npnt(i))
                call voxvol (cell(1,i),grid(1,i),cellv,voxv)
                gridv = voxv * float(npnt(i))
                maskv = voxv * float(nset(i))
                write (*,6020) i,name(i)(1:leng1(name(i))),
     +            file(i)(1:leng1(file(i))),(grid(j,i),j=1,3),
     +            (origin(j,i),j=1,3),(extent(j,i),j=1,3),
     +            (cell(j,i),j=1,6),npnt(i),nset(i),perc
                write (*,6023) cellv,voxv,gridv,maskv
                write (*,6025)
     +            (cog(j),j=1,3),
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
c ... DELETE
c
      else if (optpar(1)(1:2) .eq. 'DE') then
c
        if (nopt .lt. 2)  call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmsk
          if (select(i)) then
            if (change(i)) then
              write (*,6040) name(i)(1:leng1(name(i)))
              if (linter) then
                reply = 'N'
                call textin (
     +            ' Really DELETE this mask (Y/N) ?',reply)
                call upcase (reply)
                if (reply .ne. 'Y') goto 1210
              end if
            end if
c
            incore (i) = .false.
            change (i) = .false.
            commnt (i) = 'No comment'
            call textut (' Deleted :',name(i))
            name (i) = '!@#$%^&*?'
c
 1210       continue
          end if
        end do
c
c ... ODL
c
      else if (optpar(1)(1:2) .eq. 'OD') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        if (nopt .lt. 3) then
          optpar (3) = 'mask.odl'
          call textin (' ODL file ?',optpar(3))
        end if
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (length(optpar(3)) .lt. 1) then
          call errcon ('Invalid ODL file name')
          goto 10
        end if
c
        call odlmsk (imsk,optpar(3))
c
c ... DOT
c
      else if (optpar(1)(1:2) .eq. 'DO') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        if (nopt .lt. 3) then
          optpar (3) = 'dot.odl'
          call textin (' ODL file ?',optpar(3))
        end if
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (length(optpar(3)) .lt. 1) then
          call errcon ('Invalid ODL file name')
          goto 10
        end if
c
        if (nopt .lt. 4) then
          sphrad = -1.0
        else
          call str2r (optpar(4),xdum,ierr)
          xdum = max(0.001, xdum)
          sphrad = xdum
          call fvalut (' Sphere radius (A):',1,sphrad)
        end if
c
        call odldot (imsk,optpar(3),sphrad,maxpnt,maxmsk,mask,shadow)
c
c ... ATOM_FIT
c
      else if (optpar(1)(1:2) .eq. 'AT') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
c
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) then
          optpar (3) = ' '
          call textin (' PDB file ?',optpar(3))
        end if
c
        if (length(optpar(3)) .lt. 1) then
          call errcon ('Invalid PDB file name')
          goto 10
        end if
c
        if (nopt .lt. 4) optpar(4) = ' '
        if (nopt .lt. 5) optpar(5) = ' '
c
        call chkatm (imsk,optpar(3),optpar(4),optpar(5),
     +               maxpnt,maxmsk,mask)
c
c ... WRITE
c
      else if (optpar(1)(1:2) .eq. 'WR') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        if (nopt .lt. 3) call textin (' Which file ?',optpar(3))
        if (nopt .lt. 4) optpar(4) = 'OMASK'
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        file (imsk) = optpar(3)
        call upcase (optpar(4))
        call writem (imsk,optpar(4),ierr,maxpnt,maxmsk,mask)
c
        if (ierr .eq. 0) change(imsk) = .false.
c
        goto 10
c
c ... EZD
c
      else if (optpar(1)(1:2) .eq. 'EZ') then
c
        if (nopt .lt. 2) call textin (' Which mask ?',optpar(2))
        if (nopt .lt. 3) call textin (' EZD file ?',optpar(3))
        jmsk = whichm(optpar(2),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (length(optpar(3)) .lt. 1) then
          call errcon ('Invalid EZD file name')
          goto 10
        end if
c
        call xopxua (iunit,optpar(3),xinter(),ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening EZD file')
          goto 10
        end if
c
        call textut (' EZD Mask :',name(jmsk))
c
        call coprim (rshad,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +    extent(3,imsk),-1)
        call wroezd (rshad,extent(1,imsk),extent(2,imsk),
     +      extent(3,imsk),origin(1,imsk),grid(1,imsk),
     +      cell(1,imsk),iunit,1.0,'(1X,13F6.2)',ierr)
c
c ... ISLAND_ERASE
c
      else if (optpar(1)(1:2) .eq. 'IS') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Island Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call island (mask(1,i),extent(1,i),extent(2,i),
     +        extent(3,i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
c
c ... BLOB_ERASE
c
      else if (optpar(1)(1:2) .eq. 'BL') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          optpar (3) = '1'
          call textin (' Minimum size of blobs to KEEP ?',optpar(3))
        end if
        call str2i (optpar(3),idum,ierr)
        if (ierr .ne. 0) goto 10
        if (idum .le. 0) goto 10
c
        if (nopt .lt. 4) then
          optpar (4) = '-1'
        end if
        call str2i (optpar(4),jdum,ierr)
        if (ierr .ne. 0) goto 10
        if (jdum .gt. 0 .and. jdum .lt. idum) then
          call errcon ('Max blob size less than min size')
          goto 10
        end if
c
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Blob Mask :',name(i))
            call jvalut (' Min size of blobs to keep :',1,idum)
            if (jdum .ge. idum) then
              call jvalut (' Max size of blobs to keep :',1,jdum)
            else
              jdum = -1
              call prompt (' No max size limit on blobs')
            end if
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call blob (mask(1,i),extent(1,i),extent(2,i),
     +        extent(3,i),idum,jdum)
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
c
c ... FILL
c
      else if (optpar(1)(1:2) .eq. 'FI') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Fill Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call fillm (mask(1,i),extent(1,i),extent(2,i),
     +        extent(3,i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
c
c ... BORDER_CHECK
c
      else if (optpar(1)(1:2) .eq. 'BO') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Border Mask :',name(i))
            call border (mask(1,i),extent(1,i),extent(2,i),
     +        extent(3,i))
          end if
        end do
c
c ... NBR_COUNT
c
      else if (optpar(1)(1:2) .eq. 'NB') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Nbr Mask :',name(i))
            call nbrcnt (mask(1,i),shadow,extent(1,i),extent(2,i),
     +        extent(3,i))
          end if
        end do
c
c ... SMOOTH
c
      else if (optpar(1)(1:2) .eq. 'SM') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          optpar (3) = '15'
          call textin (' Factor ?',optpar(3))
        end if
        call str2i (optpar(3),smofac,ierr)
        if (ierr .ne. 0) goto 10
        if (smofac .lt. 1 .or. smofac .gt. 26) then
          call errcon ('Factor outside range 1 - 26')
          goto 10
        end if
c
        numcyc = 1
        if (nopt .ge. 4) then
          call str2i (optpar(4),idum,ierr)
          if (ierr .eq. 0 .and. idum .gt. 0 .and. idum .le. 99) then
            numcyc = idum
          end if
        end if
        call ivalut (' Nr of cycles :',1,numcyc)
c
        do j=1,numcyc
        call ivalut (' Cycle :',1,j)
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Smooth Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call smooth (mask(1,i),shadow,extent(1,i),extent(2,i),
     +        extent(3,i),smofac)
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
        end do
c
c ... EXPAND (= smooth with a factor of 1)
c
      else if (optpar(1)(1:2) .eq. 'EX') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        numcyc = 1
        if (nopt .ge. 3) then
          call str2i (optpar(3),idum,ierr)
          if (ierr .eq. 0 .and. idum .gt. 0 .and. idum .le. 99) then
            numcyc = idum
          end if
        end if
        call ivalut (' Nr of cycles :',1,numcyc)
c
        do j=1,numcyc
        call ivalut (' Cycle :',1,j)
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Expand Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call smooth (mask(1,i),shadow,extent(1,i),extent(2,i),
     +        extent(3,i),1)
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
        end do
c
c ... CUT
c
      else if (optpar(1)(1:2) .eq. 'CU') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) then
          optpar (3) = '15'
          call textin (' Factor ?',optpar(3))
        end if
        call str2i (optpar(3),smofac,ierr)
        if (ierr .ne. 0) goto 10
        if (smofac .lt. 1 .or. smofac .gt. 26) then
          call errcon ('Factor outside range 1 - 26')
          goto 10
        end if
c
        numcyc = 1
        if (nopt .ge. 4) then
          call str2i (optpar(4),idum,ierr)
          if (ierr .eq. 0 .and. idum .gt. 0 .and. idum .le. 99) then
            numcyc = idum
          end if
        end if
        call ivalut (' Nr of cycles :',1,numcyc)
c
        do j=1,numcyc
        call ivalut (' Cycle :',1,j)
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Cut Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call cutmsk (mask(1,i),shadow,extent(1,i),extent(2,i),
     +        extent(3,i),smofac)
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
        end do
c
c ... CONTRACT (= cut with a factor of 1)
c
      else if (optpar(1)(1:2) .eq. 'CO') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        numcyc = 1
        if (nopt .ge. 3) then
          call str2i (optpar(3),idum,ierr)
          if (ierr .eq. 0 .and. idum .gt. 0 .and. idum .le. 99) then
            numcyc = idum
          end if
        end if
        call ivalut (' Nr of cycles :',1,numcyc)
c
        do j=1,numcyc
        call ivalut (' Cycle :',1,j)
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Contract Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call cutmsk (mask(1,i),shadow,extent(1,i),extent(2,i),
     +        extent(3,i),1)
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
        end do
c
c ... NOT
c
      else if (optpar(1)(1:2) .eq. 'NO') then
c
        if (nopt .lt. 2)
     +    call textin (' Which mask ?',optpar(2))
        call selecm (optpar(2),ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        do i=1,maxmsk
          if (select(i)) then
            call textut (' Not Mask :',name(i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set before :',1,nset(i))
            call notmsk (mask(1,i),extent(1,i),extent(2,i),
     +        extent(3,i))
            call countm (i,maxpnt,maxmsk,mask)
            call jvalut (' Nr of points set after  :',1,nset(i))
            change (i) = .true.
          end if
        end do
c
c ... SIMILARITY
c
      else if (optpar(1)(1:2) .eq. 'SI') then
c
        if (nopt .lt. 2) call textin (' Which mask1 ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which mask2 ?',optpar(3))
        jmsk = whichm(optpar(3),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        if (.not. simask(imsk,jmsk)) goto 10
c
        call textut (' Similarity Mask :',name(imsk))
        call textut (' Similarity Mask :',name(jmsk))
        call simmsk (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... AND
c
      else if (optpar(1)(1:2) .eq. 'AN') then
c
        if (nopt .lt. 2) call textin (' Which mask1 ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which mask2 ?',optpar(3))
        jmsk = whichm(optpar(3),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        if (.not. simask(imsk,jmsk)) goto 10
c
        call textut (' Mask     :',name(imsk))
        call textut (' And Mask :',name(jmsk))
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set before :',1,nset(imsk))
        call andmsk (imsk,jmsk,maxpnt,maxmsk,mask)
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set after  :',1,nset(imsk))
        change (imsk) = .true.
c
c ... OR
c
      else if (optpar(1)(1:2) .eq. 'OR') then
c
        if (nopt .lt. 2) call textin (' Which mask1 ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which mask2 ?',optpar(3))
        jmsk = whichm(optpar(3),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        if (.not. simask(imsk,jmsk)) goto 10
c
        call textut (' Mask    :',name(imsk))
        call textut (' Or Mask :',name(jmsk))
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set before :',1,nset(imsk))
        call ormask (imsk,jmsk,maxpnt,maxmsk,mask)
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set after  :',1,nset(imsk))
        change (imsk) = .true.
c
c ... XOR
c
      else if (optpar(1)(1:2) .eq. 'XO') then
c
        if (nopt .lt. 2) call textin (' Which mask1 ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which mask2 ?',optpar(3))
        jmsk = whichm(optpar(3),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        if (.not. simask(imsk,jmsk)) goto 10
c
        call textut (' Mask     :',name(imsk))
        call textut (' Xor Mask :',name(jmsk))
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set before :',1,nset(imsk))
        call xormsk (imsk,jmsk,maxpnt,maxmsk,mask)
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set after  :',1,nset(imsk))
        change (imsk) = .true.
c
c ... mask1 BUTNOT mask2
c
      else if (optpar(1)(1:2) .eq. 'BU') then
c
        if (nopt .lt. 2) call textin (' Which mask1 ?',optpar(2))
        imsk = whichm(optpar(2),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 3) call textin (' Which mask2 ?',optpar(3))
        jmsk = whichm(optpar(3),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        if (.not. simask(imsk,jmsk)) goto 10
c
        call textut (' Mask        :',name(imsk))
        call textut (' Butnot Mask :',name(jmsk))
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set before :',1,nset(imsk))
        call butnot (imsk,jmsk,maxpnt,maxmsk,mask)
        call countm (imsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set after  :',1,nset(imsk))
        change (imsk) = .true.
c
c ... UNITE
c
      else if (optpar(1)(1:2) .eq. 'UN') then
c
        if (nopt .lt. 2) call textin (' New mask ?',optpar(2))
        call allocm (optpar(2),kmsk,ierr,maxmsk)
        if (ierr .ne. 0) goto 10
c
        if (nopt .lt. 3) call textin (' Which mask1 ?',optpar(3))
        imsk = whichm(optpar(3),maxmsk)
c
        if (imsk .le. 0 .or. imsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(imsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (nopt .lt. 4) call textin (' Which mask2 ?',optpar(4))
        jmsk = whichm(optpar(4),maxmsk)
c
        if (jmsk .le. 0 .or. jmsk .gt. maxmsk) then
          call errcon ('Invalid mask selection')
          goto 10
        end if
c
        if (.not. incore(jmsk)) then
          call errcon ('Mask not in memory')
          goto 10
        end if
c
        if (imsk .eq. jmsk) then
          call errcon ('Mask1 and Mask2 are identical')
          goto 10
        end if
c
        if (.not. simask(imsk,jmsk)) goto 10
c
        npnt (kmsk) = 1
        do i=1,3
          cell (i,kmsk)   = cell (i,imsk)
          cell (i+3,kmsk) = cell (i+3,imsk)
          grid (i,kmsk)   = grid (i,imsk)
          origin (i,kmsk) = min (origin(i,imsk),origin(i,jmsk))
          extent (i,kmsk) = max ( (origin(i,imsk)+extent(i,imsk)),
     +      (origin(i,jmsk)+extent(i,jmsk)) ) - origin (i,kmsk)
          npnt (kmsk) = npnt (kmsk) * extent (i,kmsk)
        end do
c
        if (npnt(kmsk) .gt. maxpnt) then
          call errcon ('United mask is too big')
          call jvalut (' Maximum  nr of points :',1,maxpnt)
          call jvalut (' Required nr of points :',1,npnt(kmsk))
          incore (kmsk) = .false.
          name (kmsk) = '!@#$%^&*()'
          goto 10
        end if
c
        incore (kmsk) = .true.
        select (kmsk) = .false.
        change (kmsk) = .true.
        file (kmsk)   = 'not_defined'
c
        call inimsk (mask(1,kmsk),extent(1,kmsk),
     +    extent(2,kmsk),extent(3,kmsk),0)
        call jvalut (' Nr of points in united mask  :',1,npnt(kmsk))
c
        call ormask (kmsk,imsk,maxpnt,maxmsk,mask)
        call countm (kmsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set after mask1 :',1,nset(kmsk))
c
        call ormask (kmsk,jmsk,maxpnt,maxmsk,mask)
        call countm (kmsk,maxpnt,maxmsk,mask)
        call jvalut (' Nr of points set after mask2 :',1,nset(kmsk))
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
            call prompt (' VRML file closed')
          else
            call prompt (' No open VRML file at present')
          end if
c
        else if (optpar(2)(1:2) .eq. 'DO') then
c
c ... VRML DOT_SURFACE
c
          if (.not. lvrml) then
            call errcon ('No VRML file opened yet !')
            goto 10
          end if
c
          if (nopt .lt. 3) then
ccc            optpar (3) = prev
            call textin (' Mask ?',optpar(3))
          end if
          imsk = whichm (optpar(3),maxmsk)
c
          if (imsk .le. 0 .or. imsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(imsk)) then
            call errcon ('Mask not in memory')
            goto 10
          end if
ccc          prev = optpar (3)
c
          if (nopt .ge. 4) then
            call xvrml_rgb_name (optpar(4),rr,gg,bb)
            call xvrml_colour (rr,gg,bb)
          end if
c
          write (*,'(99a)') ' VRML - Dot surface for mask ',
     +      name(imsk)
c
          call vrmask (imsk,ivrml,maxpnt,maxmsk,mask,shadow,ierr)
c
          if (ierr .ne. 0) then
            lvrml = .false.
            call xvrml_close ()
ccc            call errcon ('Closed VRML file')
            goto 10
          end if
c
          call flusho (ivrml)
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
          if (ierr .ne. 0) then
            lvrml = .false.
            call xvrml_close ()
ccc            call errcon ('Closed VRML file')
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
            call textin (' Mask ?',optpar(3))
          end if
          imsk = whichm (optpar(3),maxmsk)
c
          if (imsk .le. 0 .or. imsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(imsk)) then
            call errcon ('Mask not in memory')
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
          write (*,'(99a)') ' VRML - Cell for mask ',
     +      name(imsk)
c
          call xvrml_cell (cell(1,imsk),iox,ioy,ioz,lline)
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
            call textin (' Mask ?',optpar(3))
          end if
          imsk = whichm (optpar(3),maxmsk)
c
          if (imsk .le. 0 .or. imsk .gt. maxmsk) then
            call errcon ('Invalid mask selection')
            goto 10
          end if
c
          if (.not. incore(imsk)) then
            call errcon ('Mask not in memory')
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
          write (*,'(99a)') ' VRML - Box for mask ',
     +      name(imsk)
c
          call vbox (imsk,lline)
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
          if (linter) goto 10
          return
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
c --- END OF MAIN EVENT LOOP
c
c ... check if all changes saved
c
 9000 continue
c
      write (*,*)
      unsave = .false.
      do i=1,maxmsk
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
      integer function whichm (nam,maxmsk)
c
c ... which mask does the name "nam" correspond to ?
c
c ... if "*", then return 0, meaning ALL masks
c     if okay, return index of mask
c     otherwise:
c     -1 if length = 0
c     -2 if duplicate name
c     -3 if not found
c
      include 'mama.incl'
c
      integer maxmsk
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
      do i=1,maxmsk
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
      subroutine countm (imsk,maxpnt,maxmsk,mask)
c
c ... count nr of points set in mask nr IMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ns,imsk
c
code ...
c
      if (imsk .le. 0 .or. imsk .gt. maxmsk) return
      if (.not. incore(imsk)) return
c
      call cntmsk (mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +  extent(3,imsk),ns)
c
      nset (imsk) = ns
c
      return
      end
c
c
c
      subroutine centrm (imsk,maxpnt,maxmsk,mask,cog)
c
c ... count nr of points set in mask nr IMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      real amat(3,3),cog(3),dum(3),fra(3)
c
      integer ns,imsk,i
c
code ...
c
      if (imsk .le. 0 .or. imsk .gt. maxmsk) return
      if (.not. incore(imsk)) return
c
      call cogmsk (mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +  extent(3,imsk),ns,dum(1),dum(2),dum(3))
c
      nset (imsk) = ns
c
      call orthog (cell(1,imsk),amat,0)
c
      do i=1,3
        fra(i) = (dum(i)+float(origin(i,imsk)-1))/float(grid(i,imsk))
      end do
c
      call mulmtx (amat, fra, cog, 3, 3, 1)
c
ccc      write (*,*) 'COG - GPT = ',dum
ccc      write (*,*) 'COG - FRA = ',fra
ccc      write (*,*) 'COG - XYZ = ',cog
c
      return
      end
c
c
c
      subroutine selecm (nam,ierr,maxmsk)
c
c ... figure out which mask(s) are to be selected
c
      include 'mama.incl'
c
      integer maxmsk
c
      integer ierr,i,imsk,whichm
c
      character nam*(*)
c
code ...
c
      ierr = -1
c
      nmask = 0
      do i=1,maxmsk
        select (i) = .false.
        if (incore(i)) nmask = nmask + 1
      end do
c
      if (nmask .le. 0) then
        call errcon ('No masks in memory')
        return
      end if
c
      imsk = whichm(nam,maxmsk)
c
      if (imsk .lt. 0 .or. imsk .gt. maxmsk) then
        call errcon ('Invalid mask name')
        return
      end if
c
      if (imsk .eq. 0) then
        do i=1,maxmsk
          select (i) = incore (i)
        end do
        ierr = 0
      else
        if (incore(imsk)) then
          select (imsk) = .true.
          ierr = 0
        else
          call errcon ('Selected mask not in memory')
        end if
      end if
c
      return
      end
c
c
c
      subroutine writem (imsk,how,ierr,maxpnt,maxmsk,mask)
c
c ... write mask IMSK to file
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ierr,imsk,iunit
c
      logical xinter
c
      character how*(*)
c
code ...
c
      call xopxua (unt(imsk),file(imsk),xinter(),ierr)
      if (ierr .ne. 0) return
c
      iunit = unt(imsk)
c
      call wrmask (iunit,mask(1,imsk),extent(1,imsk),
     +  extent(2,imsk),extent(3,imsk),origin(1,imsk),
     +  grid(1,imsk),cell(1,imsk),how,ierr)
c
      close (iunit)
c
      return
c
 7000 format (8i5)
 7010 format (8f10.3)
 7020 format (40(1x,a1))
c
      end
c
c
c
      logical function simask (imsk,jmsk)
c
c ... check if IMSK and JMSK on same grid etc.
c
      include 'mama.incl'
c
      integer imsk,jmsk,i
c
code ...
c
      simask = .false.
c
      do i=1,6
        if (abs(cell(i,imsk)-cell(i,jmsk)) .ge. 0.01) then
          call errcon ('Masks have different cell constants')
          return
        end if
      end do
c
      do i=1,3
        if (grid(i,imsk) .ne. grid(i,jmsk)) then
          call errcon ('Masks are on different grids')
          return
        end if
      end do
c
      simask = .true.
c
      return
      end
c
c
c
      subroutine simmsk (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... similarity of masks IMSK and JMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ierr,imsk,jmsk,ilim(2,3),jlim(2,3)
c
code ...
c
      call comask (imsk,jmsk,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('SIM',
     +  mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +  mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +  ilim,jlim)
c
      return
      end
c
c
c
      subroutine andmsk (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... logical AND of masks IMSK and JMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ierr,imsk,jmsk,ilim(2,3),jlim(2,3)
c
code ...
c
      call comask (imsk,jmsk,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('AND',
     +  mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +  mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +  ilim,jlim)
c
      return
      end
c
c
c
      subroutine ormask (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... logical OR of masks IMSK and JMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ierr,imsk,jmsk,ilim(2,3),jlim(2,3)
c
code ...
c
      call comask (imsk,jmsk,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('OR',
     +  mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +  mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +  ilim,jlim)
c
      return
      end
c
c
c
      subroutine xormsk (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... logical XOR of masks IMSK and JMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ierr,imsk,jmsk,ilim(2,3),jlim(2,3)
c
code ...
c
      call comask (imsk,jmsk,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('XOR',
     +  mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +  mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +  ilim,jlim)
c
      return
      end
c
c
c
      subroutine butnot (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... IMSK but NOT JMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer ierr,imsk,jmsk,ilim(2,3),jlim(2,3)
c
code ...
c
      call comask (imsk,jmsk,ilim,jlim,ierr)
      if (ierr .ne. 0) return
c
      call combim ('BUTNOT',
     +  mask(1,imsk),extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +  mask(1,jmsk),extent(1,jmsk),extent(2,jmsk),extent(3,jmsk),
     +  ilim,jlim)
c
      return
      end
c
c
c
      subroutine comask (imsk,jmsk,ilim,jlim,ierr)
c
c ... find common ground between two masks
c
      include 'mama.incl'
c
      integer ierr,i,imsk,jmsk,ilim(2,3),jlim(2,3),ncom
      integer imin(3),jmin(3),imax(3),jmax(3),lmin(3),lmax(3)
c
code ...
c
      ierr = -1
c
c ... find absolute margins for both masks
c
      do i=1,3
        imin (i) = origin(i,imsk)
        jmin (i) = origin(i,jmsk)
        imax (i) = origin(i,imsk) + extent(i,imsk) - 1
        jmax (i) = origin(i,jmsk) + extent(i,jmsk) - 1
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
c
      do i=1,3
        if (lmin(i) .ge. lmax(i)) then
          call errcon ('Masks do not overlap')
          return
        end if
      end do
c
      ierr = 0
c
c ... find common limits
c
      do i=1,3
        ilim (1,i) = lmin(i) - origin(i,imsk) + 1
        ilim (2,i) = lmax(i) - origin(i,imsk) + 1
        jlim (1,i) = lmin(i) - origin(i,jmsk) + 1
        jlim (2,i) = lmax(i) - origin(i,jmsk) + 1
      end do
c
      call ivalut (' Lower limits common volume :',3,lmin)
      call ivalut (' Upper limits common volume :',3,lmax)
      call ivalut (' Limits first  mask  :',6,ilim)
      call ivalut (' Limits second mask  :',6,jlim)
c
      ncom = (lmax(1)-lmin(1)+1) * (lmax(2)-lmin(2)+1) *
     +       (lmax(3)-lmin(3)+1)
      call jvalut (' Number of common mask points :',1,ncom)
c
      return
      end
c
c
c
      subroutine copym (imsk,jmsk,maxpnt,maxmsk,mask)
c
c ... copy JMSK into IMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer i,imsk,jmsk
c
code ...
c
      do i=1,3
        cell (i,imsk)   = cell (i,jmsk)
        cell (i+3,imsk) = cell (i+3,jmsk)
        origin (i,imsk) = origin (i,jmsk)
        extent (i,imsk) = extent (i,jmsk)
        grid (i,imsk)   = grid (i,jmsk)
      end do
c
      select (imsk) = .false.
      incore (imsk) = .true.
      change (imsk) = .false.
      file (imsk)   = 'not_defined'
      npnt (imsk)   = npnt (jmsk)
c
      call copmij (mask(1,imsk),mask(1,jmsk),extent(1,imsk),
     +  extent(2,imsk),extent(3,imsk))
c
      call countm (imsk,maxpnt,maxmsk,mask)
      call jvalut (' Nr of points set :',1,nset(imsk))
c
      return
      end
c
c
c
      subroutine allocm (string,imsk,ierr,maxmsk)
c
c ... allocate a mask
c
      include 'mama.incl'
c
      integer maxmsk
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
      do i=1,maxmsk
        if (.not. incore(i)) then
          imsk = i
          goto 910
        end if
      end do
c
      call errcon ('No more masks available')
      return
c
  910 continue
      name (imsk) = string
c
      call upcase (name(imsk))
      if (whichm(name(imsk),maxmsk) .ne. imsk) then
        call errcon ('Invalid mask name (empty or not unique)')
        return
      end if
c
      if (incore(imsk)) then
        call errcon ('Mask in use; DELETE it first')
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
      subroutine overlp (task,imsk,ierr,nmin,xmin,ezdfil,
     +                   maxpnt,maxmsk,mask,shadow,ibuff,
     +                   rshad,buffer,ifudge)
c
c ... calculate overlap
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
      integer shadow(maxpnt),ibuff(maxpnt)
      real    buffer(maxpnt),rshad(maxpnt)
c
      integer ifudge(3)
      integer imsk,ierr,i,nbad,nmin,nunit,nasun,nunas,iunit,length
      integer margin
c
      real space(3),xmin
c
      logical xinter
c
      character task*(*),ezdfil*(*)
c
code ...
c
      ierr = -1
c
      if (nrtncs .lt. 1 .or. nsymop .lt. 1) then
        call errcon ('Not enough (NCS) symmetry operators')
        return
      end if
c
      iunit = 12
      if (length(ezdfil) .gt. 0) then
        call xopxua (iunit,ezdfil,xinter(),ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening EZD file')
          return
        end if
      end if
c
      do i=1,3
        space (i) = cell(i,imsk) / float(grid(i,imsk))
      end do
c
c ... if #NCS = 1 or 2, then simple
c
c      if (nrtncs .le. 2) then
c
c ... count incidences
c
c        write (*,*) 'Creating overlap map ...'
c        call sub05x (ibuff,extent(1,imsk),extent(2,imsk),
c     +    extent(3,imsk),origin(1,imsk),
c     +    mask(1,imsk),extent(1,imsk),extent(2,imsk),
c     +    extent(3,imsk),origin(1,imsk),cell(1,imsk),
c     +    space,rtncs,nrtncs,
c     +    symmop,nsymop,nbad,ierr)
c
c        if (ierr .ne. 0) return
c
c        call jvalut (' Nr of suspicious points :',1,nbad)
c        if (nbad .le. 0) return
c
c        if (nmin .lt. 0) then
c
c ... output EZD file
c
c          write (*,*) 'Creating EZD file ...'
c          call coprim (rshad,ibuff,extent(1,imsk),extent(2,imsk),
c     +      extent(3,imsk),-1)
c          call wroezd (rshad,extent(1,imsk),extent(2,imsk),
c     +      extent(3,imsk),origin(1,imsk),grid(1,imsk),
c     +      cell(1,imsk),iunit,1.0,'(1X,13F6.2)',ierr)
c
c        else
c
c ... count neighbours for each mask point; store in SHADOW
c
c          write (*,*) 'Counting neighbours ...'
c          call cntnbr (mask(1,imsk),shadow,extent(1,imsk),
c     +      extent(2,imsk),extent(3,imsk),1,0)
c
c          write (*,*) 'Removing overlapping mask points ...'
c          call trimov (mask(1,imsk),ibuff,shadow,
c     +      extent(1,imsk),extent(2,imsk),extent(3,imsk),2,nmin,ierr)
c
c        end if
c
c ... if > 2 NCS, then more complicated
c
c      else
c
        if (asuext(1)*asuext(2)*asuext(3) .gt. maxpnt) then
          call errcon ('Too many points in asymmetric unit')
          call jvalut (' Required :',1,
     +                 asuext(1)*asuext(2)*asuext(3))
          call jvalut (' Maximum  :',1,maxpnt)
          return
        end if
c
c ... check if really an asymmetric unit plus one or two points
c
        nunit = grid(1,imsk)*grid(2,imsk)*grid(3,imsk)
c
        do i=1,3
c
          call ivalut (' Assume margin (extra points) :',1,i)
          nasun = (asuext(1)-i)*(asuext(2)-i)*(asuext(3)-i)
          nunas = nsymop * nasun
c
          if (nunas .eq. nunit) then
c
            call jvalut (' Nr of points in unit cell  :',1,nunit)
            call jvalut (' Nr in your asymmetric unit :',1,nasun)
            call jvalut (' Nr of symmetry operators   :',1,nsymop)
            call jvalut (' Asymm. unit * symm. opers. :',1,nunas)
c
            margin = i
            call jvalut (' Extra points at all sides  :',1,margin)
            goto 6205
          end if
        end do
c
        nasun = asuext(1)*asuext(2)*asuext(3)
        nunas = nsymop * nasun
        call jvalut (' Nr of points in unit cell  :',1,nunit)
        call jvalut (' Nr in your asymmetric unit :',1,nasun)
        call jvalut (' Nr of symmetry operators   :',1,nsymop)
        call jvalut (' Asymm. unit * symm. opers. :',1,nunas)
        if (nunas .eq. nasun) then
          call errcon (
     +    'Add 1, 2, or 3 points to all sides of the asymm. unit !')
        else
          call errcon (
     +    'Not an asymmetric unit (plus 1, 2, or 3 points) !')
        end if
        return
c
c ... count incidences; store "overlap map" in integer IBUFF
c
 6205   continue
c
        if (ovmode(1:5) .eq. 'LABEL') then
          write (*,*) 'Creating overlap map (mode LABEL) ...'
          call sub05x (ibuff,asuext(1),asuext(2),asuext(3),
     +      asuori,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +      extent(3,imsk),origin(1,imsk),cell(1,imsk),
     +      space,rtncs,nrtncs,symmop,nsymop,nbad,ierr)
        else if (ovmode(1:5) .eq. 'COUNT') then
          write (*,*) 'Creating overlap map (mode COUNT) ...'
          call sub15 (rshad,asuext(1),asuext(2),asuext(3),
     +      asuori,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +      extent(3,imsk),origin(1,imsk),cell(1,imsk),
     +      space,rtncs,nrtncs,symmop,nsymop,nbad,ierr)
        else
          call errcon ('Invalid overlap mode :'//ovmode)
          return
        end if
c
        if (ierr .ne. 0) return
c
c ... OVERLAP ASU
c
        if (task(1:2) .eq. 'AS') then
c
c ... define EXTENT of asu-mask
c
          do i=1,3
            extent(i,nmin) = asuext(i) - margin + ifudge(i)
          end do
          call ivalut (' ASU extent :',3,extent(1,nmin))
c
          if (ovmode(1:5) .eq. 'LABEL') then
            i = 1
          else if (ovmode(1:5) .eq. 'COUNT') then
            i = 2
          end if
c
          call asucpy (ibuff,rshad,i,asuext(1),asuext(2),asuext(3),
     +                 mask(1,nmin),extent(1,nmin),extent(2,nmin),
     +                 extent(3,nmin))
c
          ierr = 0
          return
c
        end if
c
cc        call jvalut (' Nr of suspicious points :',1,nbad)
cc        if (nbad .le. 0) then
cc          ierr = -1
cc          return
cc        end if
c
c ... now "average": project all points onto the original mask
c                    and ADD them UP (i.e., DO NOT AVERAGE)
c
c ... copy overlap map to real RSHAD
c
        if (ovmode(1:5) .eq. 'LABEL') then
          write (*,*) 'Copying overlap map ...'
          call coprim (rshad,ibuff,asuext(1),asuext(2),asuext(3),-1)
        end if
c
c ... store "averaged" overlap map in real BUFFER
c
        write (*,*) '"Averaging" overlap map ...'
        call sub01x (rshad,asuext(1),asuext(2),asuext(3),
     +    asuori,buffer,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +    extent(3,imsk),origin(1,imsk),cell(1,imsk),space,
     +    rtncs,nrtncs,symmop,nsymop,ierr)
c
        if (ierr .ne. 0) return
c
        if (length(ezdfil) .gt. 0) then
c
c ... output EZD file
c
          write (*,*) 'Creating EZD file ...'
          call wroezd (buffer,extent(1,imsk),extent(2,imsk),
     +      extent(3,imsk),origin(1,imsk),grid(1,imsk),
     +      cell(1,imsk),iunit,1.0,'(1X,13F6.2)',ierr)
c
          close (iunit)
          if (task(1:2) .eq. 'EZ') return
c
        end if
c
c ... count neighbours for each mask point; store in SHADOW
c
        write (*,*) 'Counting neighbours ...'
        call cntnbr (mask(1,imsk),shadow,extent(1,imsk),
     +    extent(2,imsk),extent(3,imsk),1,0)
c
        write (*,*) 'Removing overlapping mask points ...'
        call fvalut (' Min count for removal       :',1,xmin)
        call ivalut (' Min nr of nbrs outside mask :',1,nmin)
        call trrmov (mask(1,imsk),buffer,shadow,
     +    extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +    xmin,nmin,ierr)
c
c      end if
c
      return
      end
c
c
c
      subroutine overun (imsk,ierr,nmin,xmin,nextra,ezdfil,
     +                   maxpnt,maxmsk,mask,shadow,ibuff,rshad,buffer)
c
c ... calculate overlap
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
      integer shadow(maxpnt),ibuff(maxpnt)
      real    buffer(maxpnt),rshad(maxpnt)
c
      integer imsk,ierr,i,nbad,nmin,nunit,nasun,iunit,length
      integer iori(3),iext(3),nextra
c
      real space(3),xmin
c
      logical xinter
c
      character ezdfil*(*)
c
code ...
c
      ierr = -1
c
      if (nrtncs .lt. 1 .or. nsymop .lt. 1) then
        call errcon ('Not enough (NCS) symmetry operators')
        return
      end if
c
      iunit = 12
      if (length(ezdfil) .gt. 0) then
        call xopxua (iunit,ezdfil,xinter(),ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening EZD file')
          return
        end if
      end if
c
      nunit = 1
      nasun = 1
      do i=1,3
        space (i) = cell(i,imsk) / float(grid(i,imsk))
        iori (i) = 0
        iext (i) = grid(i,imsk) + nextra
        nunit = nunit * grid(i,imsk)
        nasun = nasun * iext(i)
      end do
c
      call jvalut (' Unit cell origin :',3,iori)
      call jvalut (' Unit cell grid   :',3,grid(1,imsk))
      call jvalut (' Unit cell extent :',3,iext)
      call fvalut (' Unit cell spacing (A) :',3,space)
c
      call jvalut (' Nr of points in unit cell   :',1,nunit)
      call jvalut (' Nr of extra points for grid :',1,nextra)
      call jvalut (' Nr of unit cell grid points :',1,nasun)
c
      if (nasun .gt. maxpnt) then
        call errcon ('Too many points in unit cell grid')
        call jvalut (' Required :',1,nasun)
        call jvalut (' Maximum  :',1,maxpnt)
        return
      end if
c
c ... count incidences; store "overlap map" in integer IBUFF
c
      if (ovmode(1:5) .eq. 'LABEL') then
        write (*,*) 'Creating overlap map (mode LABEL) ...'
        call sub05x (ibuff,iext(1),iext(2),iext(3),
     +    iori,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +    extent(3,imsk),origin(1,imsk),cell(1,imsk),
     +    space,rtncs,nrtncs,symmop,nsymop,nbad,ierr)
      else if (ovmode(1:5) .eq. 'COUNT') then
        write (*,*) 'Creating overlap map (mode COUNT) ...'
        call sub15 (rshad,iext(1),iext(2),iext(3),
     +    iori,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +    extent(3,imsk),origin(1,imsk),cell(1,imsk),
     +    space,rtncs,nrtncs,symmop,nsymop,nbad,ierr)
      else
        call errcon ('Invalid overlap mode :'//ovmode)
        return
      end if
c
      if (ierr .ne. 0) return
c
c ... now "average": project all points onto the original mask
c                    and ADD them UP (i.e., DO NOT AVERAGE)
c
c ... copy overlap map to real RSHAD
c
      if (ovmode(1:5) .eq. 'LABEL') then
        write (*,*) 'Copying overlap map ...'
        call coprim (rshad,ibuff,iext(1),iext(2),iext(3),-1)
      end if
c
c ... store "averaged" overlap map in real BUFFER
c
      write (*,*) '"Averaging" overlap map ...'
      call sub01x (rshad,iext(1),iext(2),iext(3),
     +  iori,buffer,mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +  extent(3,imsk),origin(1,imsk),cell(1,imsk),space,
     +  rtncs,nrtncs,symmop,nsymop,ierr)
c
      if (ierr .ne. 0) return
c
      if (length(ezdfil) .gt. 0) then
c
c ... output EZD file
c
        write (*,*) 'Creating EZD file ...'
        call wroezd (buffer,extent(1,imsk),extent(2,imsk),
     +    extent(3,imsk),origin(1,imsk),grid(1,imsk),
     +    cell(1,imsk),iunit,1.0,'(1X,13F6.2)',ierr)
c
        close (iunit)
c
      end if
c
c ... count neighbours for each mask point; store in SHADOW
c
      write (*,*) 'Counting neighbours ...'
      call cntnbr (mask(1,imsk),shadow,extent(1,imsk),
     +  extent(2,imsk),extent(3,imsk),1,0)
c
      write (*,*) 'Removing overlapping mask points ...'
      call fvalut (' Min count for removal       :',1,xmin)
      call ivalut (' Min nr of nbrs outside mask :',1,nmin)
      call trrmov (mask(1,imsk),buffer,shadow,
     +  extent(1,imsk),extent(2,imsk),extent(3,imsk),
     +  xmin,nmin,ierr)
c
      return
      end
c
c
c
      subroutine odlmsk (imsk,odlfil)
c
c ... write ODL file to draw box around mask IMSK
c
      include 'mama.incl'
c
      real a(3,3),x(3),flo(3),fhi(3)
c
      integer i,imsk,iunit,ierr
c
      logical xinter
c
      character odlfil*(*)
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
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell(1,imsk), a, 0)
c
c ... calc origin and top in fractionals
c
      do i=1,3
        flo (i) = float(origin(i,imsk)) * 
     +           (cell(i,imsk)/float(grid(i,imsk))) / cell(i,imsk)
        fhi (i) = float(extent(i,imsk)+origin(i,imsk)-1) *
     +           (cell(i,imsk)/float(grid(i,imsk))) / cell(i,imsk)
      end do
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
      write (iunit,6000) 'begin mskgr'
      write (iunit,6000) 'colour 16799999'
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
      subroutine chkatm (imsk,pdbfil,fil1,fil2,maxpnt,maxmsk,mask)
c
c ... check if atoms from PDB file fit inside mask/grid
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
c
      integer imsk,iunit,ierr,junit,kunit,length
c
      logical xinter
c
      character pdbfil*(*),fil1*(*),fil2*(*)
c
code ...
c
      iunit = 41
      junit = 42
      kunit = 43
c
      call xopxoa (iunit,pdbfil,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening input PDB file')
        return
      end if
c
      if (length(fil1) .gt. 0) then
        call xopxua (junit,fil1,xinter(),ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening "inside" PDB file')
          return
        end if
      else
        junit = -1
      end if
c
      if (length(fil2) .gt. 0) then
        call xopxua (kunit,fil2,xinter(),ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening "outside" PDB file')
          return
        end if
      else
        kunit = -1
      end if
c
      call lmokay (mask(1,imsk),extent(1,imsk),extent(2,imsk),
     +  extent(3,imsk),origin(1,imsk),grid(1,imsk),cell(1,imsk),
     +  newrad,iunit,junit,kunit)
c
      close (iunit)
      if (junit .gt. 0) close (junit)
      if (kunit .gt. 0) close (kunit)
c
      return
      end
c
c
c
      subroutine odldot (imsk,odlfil,sphrad,maxpnt,maxmsk,mask,shadow)
c
c ... write ODL file to draw dot surface of mask IMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
      integer shadow(maxpnt)
c
      real a(3,3),sphrad
c
      integer imsk,iunit,ierr
c
      logical xinter
c
      character odlfil*(*)
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
      call cntnbr (mask(1,imsk),shadow,extent(1,imsk),
     +             extent(2,imsk),extent(3,imsk),1,0)
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell(1,imsk), a, 0)
c
      call dottum (mask(1,imsk),shadow,extent(1,imsk),
     +             extent(2,imsk),extent(3,imsk),grid(1,imsk),
     +             origin(1,imsk),cell(1,imsk),a,sphrad,iunit,ierr)
c
      return
      end
c
c
c
      subroutine vrmask (imsk,iunit,maxpnt,maxmsk,mask,shadow,ierr)
c
c ... write VRML file to draw dot surface of mask IMSK
c
      include 'mama.incl'
c
      integer maxpnt,maxmsk
      integer mask(maxpnt,maxmsk)
      integer shadow(maxpnt)
c
      real a(3,3)
c
      integer imsk,iunit,ierr
c
code ...
c
      call cntnbr (mask(1,imsk),shadow,extent(1,imsk),
     +             extent(2,imsk),extent(3,imsk),1,0)
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell(1,imsk), a, 0)
c
      call dovrml (mask(1,imsk),shadow,extent(1,imsk),
     +             extent(2,imsk),extent(3,imsk),grid(1,imsk),
     +             origin(1,imsk),cell(1,imsk),a,iunit,ierr)
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
      include 'mama.incl'
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
