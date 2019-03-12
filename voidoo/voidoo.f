      program voidoo
c
c === VOIDOO === calculate cavities in or volume of a (protein) structure
c
c Gerard Kleywegt @ 920807
c
c ===========================================================================
c
c Version history:
c ----------------
c
c Version 0.1 @ 920807 - basic algorithm & documentation
c Version 0.2 @ 920810 - altered definition of Vanderwaals radii
c                      - sped up protein finder algorithm by factor 400
c                      - atom fattening to detect "true" nr of cavities
c                      - output proper mask file
c                      - implemented volume refinement
c                      - implemented simple ODL files
c Version 0.3 @ 920811 - also read residue definitions from library file
c                      - removed references to unit cell data (not used)
c                      - add/subtract grid spacing to/from lower/upper limits
c                        to account for the fact that we don't use the grid
c                        points on the borders
c                      - added probe radius
c                      - implemented iterative zapping instead of recursive
c                        (it's a bit faster and it doesn't crash with huge
c                        grids)
c                      - use character*1 parameters for filling the grid
c                      - added volume option
c                      - made MASK and ODL files optional
c                      - implemented trace option to reduce output for norma
c                        non-debugging runs
c Version 0.4 @ 920812 - improved test file vol.pdb for volume testing
c                      - took some Vanderwaals radii from the AMBER paper:
c                        SJ Weiner et al., JACS 106, 765-784 (1984)
c                      - implemented refinement for volume calculation option
c                      - exclude borders in iterative zapping routines
c                      - debugged grid-definition code (... aaarrgghhh ...)
c                      - made it possible to actually perform all N detection
c                        cycles (rather than only the "first best" option)
c                      - made refinement optional
c                      - implemented percentage convergence criterion
c                      - removed several smaller bugs which only showed up in
c                        "unfortunate" cases
c Version 0.5 @ 920813 - implemented 1-sweep contour odl-output
c                      - implemented 3-sweep contour odl-output
c                      - updated documentation
c Version 0.6 @ 920814 - implemented residue/atom combinations in vdw-lib
c Version 0.7 @ 920820 - implemented dots odl-output
c                      - added option to exclude probe from cavity volume
c                      - print average/sigma of volumes calculated
c Version 0.8 @ 920821 - use probe radius in volume calculation (if 1.4 A =>
c                        gives solvent-excluded volume)
c                      - implemented option to look for a specific cavity by
c                        providing XYZ-coordinates of a point you think or know
c                        is inside it
c                      - implemented output of ODL file which contains dots for
c                        the protein surface points
c Version 0.9 @ 920825 - implemented EZD plot output for cavities
c Version 1.0 @ 920826 - implemented EZD plot output for protein itself
c Version 1.1 @ 920922 - don't bail out if primary grid is too small and the
c                        program is run interactively
c                      - write cavity box to log file (for future option to
c                        cut out a certain box to do the calculations on in
c                        case you are looking for a specific cavity and want
c                        to use a very fine primary grid)
c                      - only allow error recovery in file open routines
c                        (XOPxxx) if program is run interactively
c Version 1.2 @ 920923 - implement limits on X,Y,Z to see if it works
c Version 1.3 @ 921020 - problems with limit option -> removed it again
c                        (couldn't always find "outside world")
c Version 1.4 @ 921029 - altered find_spec_cavity so that the program
c                        checks a 11*11*11 cube around the coordinates
c                        that the user entered (used to be 5*5*5);
c                        also changed so that closest point is selected;
c                        doubled buffer size to 2 megawords
c Version 1.5 @ 930202 - print radius of corresponding sphere whenever a
c                        (cavity or molecular) volume is calculated;
c                        doubled buffer size again, now to 4 megawords;
c                        implemented new EZD format
c Version 1.6 @ 930222 - minor changes
c Version 2.0 @ 930302 - attempted to implement Connolly's "solvent-
c                        occupied volume" ... it works !!!
c Version 2.1 @ 930607 - implement molecule pre-rotation
c Version 2.2 @ 930618 - print list of non-included atoms inside cavities
c                        also print list of protein atoms lining it
c                        (maybe this should be optional ?)
c Version 2.3 @ 930630 - removed bugs from these two new options
c Version 2.4 @ 930809 - changed log-file format; new feature: IF plot
c                        files are requested for the cavities, THEN also
c                        O-macros to draw just the residues that are INSIDE
c                        and those that are LINING the cavities are made
c Version 2.4.1 930823 - include 'centre_xyz cavity_center_of_gravity'
c			       in output O macros
c Version 2.4.2 930825 - made NEWEZD format compatible with O 5.9.1
c Version 3.0 @ 931116 - recognise HETATM cards; echo REMARK cards;
c                        add comments to rotated PDB files; updated manual
c
c ===========================================================================
c
      include 'cavity.incl'
c
code ...
c
      call gkinit (prognm,version)
c
c ... initialise program
c
      ierror = 0
c
      call cav_init
      if (ierror .ne. nix) goto 9999
c
cc      call gkdcpu (total,user,sys)
cc      write (tty,'(1x,a,3f10.1)') ' 1 CPU total/user/sys :',
cc     +  total,user,sys
c
c ... get all relevant data
c
      call get_data
      if (ierror .ne. nix) goto 9999
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') ' 2 CPU total/user/sys :',
     +  total,user,sys
c
      if (ctype .eq. 'C') then
c
        if (lspec) then
c
          call find_spec_cav
c
          if (ierror .ne. nix) goto 9999
          if (numcav .le.   0) goto 9999
          if (.not.     lrefi) goto 9999
c
c ... refine cavity volumes
c
          call ref_cav_vol
c
        else
c
c ... get nr of cavities
c
          call get_nr_cav
c
          if (ierror .ne. nix) goto 9999
          if (numcav .le.   0) goto 9999
          if (.not.     lrefi) goto 9999
c
c ... refine cavity volumes
c
          call ref_cav_vol
c
        end if
c
      else if (ctype .eq. 'V') then
c
        call calc_volume
c
      else if (ctype .eq. 'R') then
c
        call rotate_mol
c
      else if (ctype .eq. 'Q') then
c
        goto 9998
c
      end if
c
cc      call gkdcpu (total,user,sys)
cc      write (tty,'(1x,a,3f10.1)') ' 3 CPU total/user/sys :',
cc     +  total,user,sys
c
c ... THE END (almost)
c
 9999 continue
c
c ... create surface-dot file if requested
c
      if (ctype .ne. 'R') then
c
        if (lpdot) call prot_surf
c
        close (iunit)
        close (junit)
        close (flog)
c
      end if
c
 9998 continue
c
      call gkquit
c
      stop
      end
c
c ===========================================================================
c
      subroutine cav_init
c
      include 'cavity.incl'
c
code ...
c
      write (tty,1200)
      write (tty,1000) 'Current version',version
      write (tty,1100) 'Max nr of atoms',maxatm
      write (tty,1100) 'Max nr of elements',maxelm
      write (tty,1100) 'Max nr of atom types',maxaty
      write (tty,1100) 'Max nr of residue/atom types',maxrac
      write (tty,1100) 'Max nr of residue types',maxres
      write (tty,1100) 'Max nr of cavities',maxhol
      write (tty,1100) 'Max nr of detection cycles',maxcyc
      write (tty,1100) 'Max nr of refinement cycles',maxref
      write (tty,1100) 'Max nr of grid points per axis',maxgrd
      write (tty,1100) 'Max nr of points for grid 1',maxbuf
      write (tty,1100) 'Max nr of points for grid 2',maxbuf
      write (tty,1100) 'Max nr of points for grid 3',maxgrd*maxgrd
      write (tty,1100) 'Memory use (Bytes) for major arrays',memuse
      write (tty,1200)
c
c ... Vanderwaals radii
c
      do i=1,maxelm
        namelm (i) = '??'
        vdwrad (i) = 2.00
      end do
c
      do i=1,maxaty
        namaty (i) = '????'
        vdwrat (i) = 2.00
      end do
c
      do i=1,maxrac
        namrac (i) = '??? ????'
        vdwrac (i) = 2.00
      end do
c
c ... Residue types
c
      do i=1,maxres
        resnam (i) = '???'
      end do
c
c ... Other stuff
c
      ctype  = 'C'
      gsize  = 1.0
      numcyc = 10
      fvdw   = 1.1
      minsiz = 1
      mingrd = 10
      numref = 10
      gfact  = 0.9
      vtoler = 0.1
      ptoler = 0.1
      pgrid  = 0.7
      dgrid  = 0.7
      probe  = 1.4
      dtype  = 'F'
      ptype  = 'N'
      cavxyz (1) = 0.0
      cavxyz (2) = 0.0
      cavxyz (3) = 0.0
c
      ierror = 0
      return
c
 1000 format (1x,a40,' : ',a)
 1100 format (1x,a40,' : ',i10)
 1200 format (/' ***** ',5('VOIDOO ***** ')/)
 1300 format (1x,a40,' : ',f10.5)
c
      end
c
c ===========================================================================
c
      subroutine get_data
c
      include 'cavity.incl'
c
      real xat,ave,sdv,xtot,ytot,ztot,savem(2)
      real save1,save2,save3,save4,save5,save6,totvol
c
      integer ndum,n2,length,leng1
c
      character residu*3,elemnt*2,pline*80,mline*80,atnam*4,racnm*8
c
code ...
c
      write (tty,1200) ' ','Data input','----------',' '
c
c ... type of calculation
c
      write (tty,'(10(1x,a/))')
     +  'Select one of the following types of calculation:',
     +  'C = cavity calculations',
     +  'V = volume calculations',
     +  'R = rotate a molecule',
     +  'Q = Quit program'
c
      call textin (
     +    ' Type of calculation (C/V/R/Q)         ?',ctype)
      call upcase (ctype)
c
      if (ctype .ne. 'V' .and. ctype .ne. 'R' .and.
     +    ctype .ne. 'C' .and. ctype .ne. 'Q') then
        call errcon ('Invalid option; aborting')
        call gkquit
      end if
c
      if (ctype .eq. 'R') return
c
      if (ctype .eq. 'Q') return
c
      answer = 'N'
      call textin (
     +    ' Do you want extensive output          ?',answer)
      call upcase (answer)
      ltrace = (answer .eq. 'Y')
c
c ... Vanderwaals radii
c
      write (tty,1200) ' ','(1) Vanderwaals radii and residue types',' '
c
      line = 'cavity.lib'
c
c ... 950118 - check if environment variable GKLIB is defined
c
      call gklibf (line)
c
      call textin (' Library file ?',line)
      close (iunit)
      call xopxoa (iunit,line,xinter(),ierror)
      if (ierror .ne. nix) then
        call errcon ('Could not open library file')
        ierror = -8
        return
      end if
c
      numelm = 0
      numaty = 0
      numrac = 0
      numres = 0
      ndum   = 0
      write (tty,1200) 'Reading your library file ...'
c
  110 continue
      read (iunit,'(a)',end=899) line
      call upcase (line)
      ndum = ndum + 1
      if (line(1:4) .eq. 'REMA') goto 110
      if (line(1:4) .eq. 'END ') goto 899
c
      if (line(1:4) .eq. 'ELEM') then
        if (numelm .lt. maxelm) then
          numelm = numelm + 1
          read (line(5:),*) namelm(numelm),vdwrad(numelm)
        else
          call errcon ('Too many elements')
          call textut (' SKIPPED :',line(5:))
        end if
      else if (line(1:4) .eq. 'ATOM') then
        if (numaty .lt. maxaty) then
          numaty = numaty + 1
          read (line(5:),*) namaty(numaty),vdwrat(numaty)
        else
          call errcon ('Too many atom types')
          call textut (' SKIPPED :',line(5:))
        end if
      else if (line(1:4) .eq. 'SPAT') then
        if (numrac .lt. maxrac) then
          numrac = numrac + 1
          read (line(5:),*) namrac(numrac),vdwrac(numrac)
          namrac(numrac) (4:4) = ' '
        else
          call errcon ('Too many residue/atom types')
          call textut (' SKIPPED :',line(5:))
        end if
      else if (line(1:4) .eq. 'RESI') then
        if (numres .lt. maxres) then
          numres = numres + 1
          read (line(5:),*) resnam(numres)
        else
          call errcon ('Too many residue types')
          call textut (' SKIPPED :',line(5:))
        end if
      else
        call errcon ('Unknown keyword')
        call textut (' Line :',line)
      end if
      goto 110
c
c ... done reading
c
  899 continue
      call jvalut (' Nr of lines in library file :',1,ndum)
      call jvalut (' Nr of elements defined      :',1,numelm)
      call jvalut (' Nr of atom types defined    :',1,numaty)
      call jvalut (' Nr of residue/atom types    :',1,numrac)
      call jvalut (' Nr of residue types defined :',1,numres)
c
      if (numelm .gt. 0) write (tty,1210)
     +  (namelm(i),vdwrad(i),i=1,numelm)
      if (numaty .gt. 0) write (tty,1212)
     +  (namaty(i),vdwrat(i),i=1,numaty)
      if (numrac .gt. 0) write (tty,1214)
     +  (namrac(i),vdwrac(i),i=1,numrac)
      if (numres .gt. 0) write (tty,1220) (resnam(i),i=1,numres)
c
      if (numelm + numaty + numrac .le. 0) then
        call errcon (' No Vanderwaals radii defined')
        ierror = -6
        return
      end if
c
      if (numres .le. 0) then
        call errcon (' No residue types defined')
        ierror = -7
        return
      end if
c
c ... PDB file
c
      write (tty,1200) ' ','(2) PDB file',' '
      pline = 'in.pdb'
      call textin (' PDB file name ?',pline)
      close (iunit)
      call xopxoa (iunit,pline,xinter(),ierror)
      if (ierror .ne. nix) then
        call errcon ('Could not open PDB file')
        ierror = -1
        return
      end if
c
      ndum   = 0
      natoms = 0
      nrejec = 0
      n2 = numres
      vdwmax = 0.0
      totvol = 0.0
      write (tty,1200) 'Reading your PDB file ...'
c
   10 continue
      if (natoms .ge. maxatm .or. nrejec .ge. maxatm) then
        call errcon ('Too many atoms - rest skipped !!!')
        goto 999
      end if
c
      read (iunit,'(a)',end=999) line
      call upcase (line)
c
ccc      write (tty,'(1x,a)') line(1:leng1(line))
c
      if (line(1:5) .eq. 'ATOM ') goto 11
      if (line(1:6) .eq. 'HETATM') goto 11
      if (line(1:6) .eq. 'REMARK') 
     +  write (tty,'(1x,a)') line(1:leng1(line))
c
      goto 10
c
   11 continue
      ndum = ndum + 1
c
c ... residue okay ?
c
      residu = line (18:20)
      do i=1,numres
        if (residu .eq. resnam(i)) goto 20
      end do
c
c ... store rejected atoms
c
      nrejec = nrejec + 1
ccc      read (line(29:),*) xr(nrejec),yr(nrejec),zr(nrejec)
      read (line(31:),'(3f8.3)')
     +  xr(nrejec),yr(nrejec),zr(nrejec)
      rejnam (nrejec) = line (13:27)
c
c ... keep track of rejected residue types
c
      if (n2 .gt. numres) then
        do i=numres+1,n2
          if (residu .eq. resnam(i)) goto 10
        end do
      end if
      n2 = n2 + 1
      resnam (n2) = residu
      goto 10
c
c ... radius defined for this atom/element ?
c
   20 continue
c
      if (numrac .gt. 0) then
        racnm = residu//' '//line(13:16)
        do i=1,numrac
          if (racnm .eq. namrac(i)) then
            natoms = natoms + 1
            dvdw (natoms) = vdwrac (i)
            vdwmax = max (vdwmax, dvdw(natoms))
            totvol = totvol + (dvdw (natoms))**3
ccc            read (line(29:),*) xd(natoms),yd(natoms),zd(natoms)
            read (line(31:),'(3f8.3)')
     +        xd(natoms),yd(natoms),zd(natoms)
            atmnam (natoms) = line (13:27)
            goto 10
          end if
        end do
      end if
c
      if (numaty .gt. 0) then
        atnam = line(13:16)
        do i=1,numaty
          if (atnam .eq. namaty(i)) then
            natoms = natoms + 1
            dvdw (natoms) = vdwrat (i)
            vdwmax = max (vdwmax, dvdw(natoms))
            totvol = totvol + (dvdw (natoms))**3
ccc            read (line(29:),*) xd(natoms),yd(natoms),zd(natoms)
            read (line(31:),'(3f8.3)')
     +        xd(natoms),yd(natoms),zd(natoms)
            atmnam (natoms) = line (13:27)
            goto 10
          end if
        end do
      end if
c
      if (numelm .gt. 0) then
        elemnt = line (13:14)
        do i=1,numelm
          if (namelm(i) .eq. elemnt) then
            natoms = natoms + 1
            dvdw (natoms) = vdwrad (i)
            vdwmax = max (vdwmax, dvdw(natoms))
            totvol = totvol + (dvdw (natoms))**3
ccc            read (line(29:),*) xd(natoms),yd(natoms),zd(natoms)
            read (line(31:),'(3f8.3)')
     +        xd(natoms),yd(natoms),zd(natoms)
            atmnam (natoms) = line (13:27)
            goto 10
          end if
        end do
      end if
c
c ... unknown atom/element
c
      call errcon ('Element without Vanderwaals radius')
      call textut (' Element :',elemnt)
      call textut (' Atom    :',line(1:28))
      ierror = -1
      return
c
c ... end of PDB file
c
  999 continue
c
      call jvalut (' Number of atoms read       :',1,ndum)
      call jvalut (' Number of atoms kept       :',1,natoms)
      call jvalut (' Number of atoms rejected   :',1,nrejec)
      call fvalut (' Max Vanderwaals radius (A) :',1,vdwmax)
c
      if (natoms .lt. 1) then
        call errcon ('No atoms kept !')
        ierror = -1
        return
      end if
c
      do i=1,natoms
        dvdw2(i)=dvdw(i)*dvdw(i)
      end do
c
      totvol = 4.0 * pi * totvol / 3.0
      call rvalut (' Sum of atomic volumes (A3) :',1,totvol)
c
      if (n2 .gt. numres) then
        call asciut (' Rejected residue types :',n2-numres,
     +    resnam(numres+1))
      else
        write (tty,1200) 'No residue types rejected'
      end if
c
c ... enough atoms ?
c
      if (ctype .eq. 'C' .and. natoms .lt. 10) then
        call errcon ('Not enough atoms (min = 10)')
        ierror = -2
        return
      end if
c
c ... Box and grid
c
      write (tty,1200) ' ','(3) Primary grid',' '
c
      call xstats (xd,natoms,ave,sdv,xmin,xmax,xtot)
      call xstats (yd,natoms,ave,sdv,ymin,ymax,ytot)
      call xstats (zd,natoms,ave,sdv,zmin,zmax,ztot)
      xat = 1.0/float(natoms)
c
      write (tty,1230) 'X',xmin,xmax,xtot*xat
      write (tty,1230) 'Y',ymin,ymax,ytot*xat
      write (tty,1230) 'Z',zmin,zmax,ztot*xat
c
cc      savem (1) = xmin
cc      savem (2) = xmax
cc      call rvalin (' Limits on X ?',2,savem)
cc      call rlohi (savem(1),savem(2))
cc      xmin = max (xmin, savem(1))
cc      xmax = min (xmax, savem(2))
c
cc      savem (1) = ymin
cc      savem (2) = ymax
cc      call rvalin (' Limits on Y ?',2,savem)
cc      call rlohi (savem(1),savem(2))
cc      ymin = max (ymin, savem(1))
cc      ymax = min (ymax, savem(2))
c
cc      savem (1) = zmin
cc      savem (2) = zmax
cc      call rvalin (' Limits on Z ?',2,savem)
cc      call rlohi (savem(1),savem(2))
cc      zmin = max (zmin, savem(1))
cc      zmax = min (zmax, savem(2))
c
cc      write (tty,1230) 'X',xmin,xmax
cc      write (tty,1230) 'Y',ymin,ymax
cc      write (tty,1230) 'Z',zmin,zmax
c
      save1 = xmin
      save2 = xmax
      save3 = ymin
      save4 = ymax
      save5 = zmin
      save6 = zmax
c
 6996 continue
      call fvalin (' Primary grid spacing (A) ?',1,gsize)
      gsize = max (0.01, gsize)
c
      probe = 1.4
      if (ctype .eq. 'V') probe = 0.0
      call fvalin (' Probe radius (1.4 A for water) ?',1,probe)
      probe = max (0.0, probe)
c
c ... NOTE: the following box sizes will make that if the
c           Vanderwaals growth radius becomes too large,
c           you will get spurious cavities, which are merely
c           parts of the outside world no longer connected
c           to the outside world as detected by the program !
c
      xmin = save1 - vdwmax - probe
      ymin = save3 - vdwmax - probe
      zmin = save5 - vdwmax - probe
      xmax = save2 + vdwmax + probe
      ymax = save4 + vdwmax + probe
      zmax = save6 + vdwmax + probe
c
      call def_grd (gsize,xmin,xmax,ymin,ymax,zmin,zmax,
     +              ngrid(1),ngrid(2),ngrid(3))
c
      write (tty,1230) 'X',xmin,xmax
      write (tty,1230) 'Y',ymin,ymax
      write (tty,1230) 'Z',zmin,zmax
c
      call jvalut (' Number of grid points :',3,ngrid)
      nxy  = ngrid(1)*ngrid(2)
      nxyz = ngrid(1)*ngrid(2)*ngrid(3)
c
      if (nxyz     .gt. maxbuf .or.
     +    ngrid(1) .gt. maxgrd .or.
     +    ngrid(2) .gt. maxgrd .or.
     +    ngrid(3) .gt. maxgrd ) then
        call errcon ('Grid too big for buffer; reduce spacing')
        if (xinter()) goto 6996
        ierror = -3
        return
      end if
c
      do i=1,ngrid(1)
        xg (i) = xmin + float(i-1)*gsize
      end do
      do i=1,ngrid(2)
        yg (i) = ymin + float(i-1)*gsize
      end do
      do i=1,ngrid(3)
        zg (i) = zmin + float(i-1)*gsize
      end do
c
      volppt = (xmax-xmin)*(ymax-ymin)*(zmax-zmin) /
     +  float( (ngrid(1)-1)*(ngrid(2)-1)*(ngrid(3)-1) )
      call rvalut (' Volume per voxel (A3) :',1,volppt)
c
c ... Parameters
c
      write (tty,1200) ' ','(4) Various parameters',' '
c
      if (ctype .eq. 'C') then
c
        call jvalin (
     +    ' Nr of detection cycles                ?',1,numcyc)
        numcyc = max (1, min (numcyc,maxcyc))
c
        call fvalin (
     +    ' Growth factor for Vanderwaals radii   ?',1,fvdw)
        fvdw = max (1.01, min (10.0, fvdw))
c
      end if
c
      if (ctype .eq. 'C') then
c
        call jvalin (
     +    ' Min size of "real" cavities (voxels)  ?',1,minsiz)
        minsiz = max (1, minsiz)
c
        answer = 'N'
        call textin (
     +    ' Are you looking for a specific cavity ?',answer)
        call upcase (answer)
        lspec = (answer .eq. 'Y')
c
        if (lspec) then
          call fvalin (
     +      ' Enter coordinates of a cavity point   :',3,cavxyz)
          if (cavxyz(1) .le. xmin .or. cavxyz(1) .gt. xmax .or.
     +        cavxyz(2) .le. ymin .or. cavxyz(2) .gt. ymax .or.
     +        cavxyz(3) .le. zmin .or. cavxyz(3) .gt. zmax) then
            call errcon ('Cavity seed point outside grid')
            ierror = -63
            return
          end if
        else
          call textin (
     +      ' All cycles or First decline (A/F)     ?',dtype)
          call upcase (dtype)
          if (numcyc .eq. 1 .or. dtype .ne. 'A') dtype = 'F'
        end if
c
        answer = 'N'
        call textin (
     +    ' Do you want a MASK file               ?',answer)
        call upcase (answer)
        lmask = (answer .eq. 'Y')
c
        if (lmask) then
c
          mline = 'out.mask'
          call textin (
     +      ' MASK file name                        ?',mline)
          close (iunit)
          call xopxua (iunit,mline,xinter(),ierror)
          if (ierror .ne. nix) then
            call errcon ('Could not open MASK file')
            ierror = -1
            return
          end if
c
        end if
c
        line = 'cavity.log'
        call textin (
     +    ' LOG file name                         ?',line)
        call xopxua (flog,line,xinter(),ierror)
        if (ierror .ne. nix) then
          call errcon ('Could not open LOG file')
          ierror = -1
          return
        end if
c
        answer = 'N'
        if (dtype .ne. 'A') answer = 'Y'
        call textin (
     +    ' Do you want to refine the cavities    ?',answer)
        call upcase (answer)
        lrefi = (answer .eq. 'Y')
c
      end if
c
      if (lrefi .or. ctype .eq. 'V') then
c
        call jvalin (
     +    ' Nr of volume-refinement cycles        ?',1,numref)
        numref = max (1, min (maxref, numref))
c
        if (ctype .eq. 'C') then
          vtype = 'A'
          write (tty,1200)
     +      'You may calculate one of these type of cavity volumes:',
     +      '(V) Vanderwaals: not occupied by the molecule',
     +      '(A) Probe-accessible: available to the probe centre',
     +      '(O) Probe-occupied: occupied by rolling probe'
          call textin (
     +      ' Type of cavity volume (V/A/O)         ?',vtype)
          call upcase (vtype)
          if (vtype .eq. 'V') then
            lprob = .false.
          else if (vtype .eq. 'A') then
            lprob = .true.
          else
            vtype = 'O'
            lprob = .true.
          end if
        end if
c
        call fvalin (
     +    ' Grid-shrink factor                    ?',1,gfact)
        gfact = max (0.1, min (0.999, gfact))
c
        if (ctype .eq. 'C') then
c
          call jvalin (
     +      ' Min size of secondary grid            ?',1,mingrd)
          mingrd = max (3, min (maxgrd, mingrd))
c
        end if
c
        call fvalin (
     +    ' Convergence criterion (A3)            ?',1,vtoler)
        vtoler = max (0.001, vtoler)
c
        call fvalin (
     +    ' Convergence criterion (%)             ?',1,ptoler)
        ptoler = max (0.001, min (99.999, ptoler))
c
      end if
c
      answer = 'N'
      call textin (
     +  ' Create protein-surface plot file      ?',answer)
      call upcase (answer)
      lpdot = (answer .eq. 'Y')
c
      if (lpdot) then
c
        apdot = 'N'
        call textin (
     +      ' Dots, Old-EZD or New-EZD (D/O/N)      ?',apdot)
        call upcase (apdot)
        if (apdot .ne. 'D' .and. apdot .ne. 'O') apdot = 'N'
c
        dfile = 'prot_surf.odl'
        if (apdot .eq. 'O') dfile = 'protein.oldezd'
        if (apdot .eq. 'N') dfile = 'protein.ezd'
        call textin (
     +      ' Name of this file                     ?',dfile)
c
        call fvalin (
     +      ' Grid spacing to use                   ?',1,dgrid)
        dgrid = max (0.01, dgrid)
c
      end if
c
c ... if only a volume calculation, return
c
      if (ctype .eq. 'V') then
c
        ierror = 0
        return
      end if
c
      answer = 'Y'
      call textin (
     +    ' Do you want plot files                ?',answer)
      call upcase (answer)
      lodl = (answer .eq. 'Y')
c
      if (lodl) then
c
        ofile = 'cavity'
        call textin (
     +    ' First part of plot file names         ?',ofile)
c
        call fvalin (
     +    ' Grid for plot files                   ?',1,pgrid)
        pgrid = max (0.01, pgrid)
c
        write (tty,1200)
     +    'You may choose from the following types of',
     +    'graphical representations for your cavities:',
     +    ' O * generate Old-EZD files',
     +    ' N * generate New-EZD files',
     +    ' D * draw dots for all cavity points',
     +    ' 3 * 3-sweep contour (fairly quick)',
     +    ' 1 * 1-sweep contour (quick and dirty)',
     +    ' C * connect all surface points (fast/big objects)',
     +    ' T * tiles (ie, poly-triangles; not implemented yet)'
        call textin (
     +    ' Graphical representation (C/1/3/D/T/O/N)?',ptype)
        call upcase (ptype)
        if (index ('C13TDON',ptype) .le. 0) ptype = 'N'
c
      end if
c
      write (flog,1200) ' ',' *** VOIDOO ***',' '
      write (flog,1200) ' ','- Vanderwaals radii'
      if (numelm .gt. 0) write (flog,1210)
     +  (namelm(i),vdwrad(i),i=1,numelm)
      if (numaty .gt. 0) write (flog,1212)
     +  (namaty(i),vdwrat(i),i=1,numaty)
      if (numrac .gt. 0) write (flog,1214)
     +  (namrac(i),vdwrac(i),i=1,numrac)
      write (flog,1200) ' ','- Residue types'
      write (flog,1220) (resnam(i),i=1,numres)
      write (flog,1200) ' ','- PDB file',pline(1:leng1(pline))
      if (lmask) write (flog,1200) ' ','- MASK file',
     +  mline(1:leng1(mline))
      write (flog,1200) ' ','- LOG file',line(1:leng1(line))
      if (lodl) write (flog,1200) ' ','- Plot files',
     +  (ofile(1:leng1(ofile))//'_#.o')
      write (flog,1200) ' ','- Grid & cavity parameters'
      write (flog,1202) 'Max nr of cavity-detection cycles   : ',numcyc
      if (lspec) then
        write (flog,1204) 'Specific cavity to be found near    : ',
     +    (cavxyz(i),i=1,3)
      else
        if (dtype .eq. 'A') then
          write (flog,1200) 'Will try ALL detection cycles to find best'
        else
          write (flog,1200) 'Will do detection cycles until first best'
        end if
      end if
      write (flog,1204) 'Vanderwaals growth factor           : ',fvdw
      write (flog,1204) 'Probe radius                        : ',probe
      write (flog,1202) 'Min nr of voxels in "real" cavities : ',minsiz
      write (flog,1202) 'Max nr of volume-refinement cycles  : ',numref
      if (vtype .eq. 'V') then
        write (flog,1200)
     +    'Will calculate Vanderwaals volumes'
      else if (vtype .eq. 'A') then
        write (flog,1200)
     +    'Will calculate probe-accessible volumes'
      else
        write (flog,1200)
     +    'Will calculate probe-occupied volumes'
      end if
      if (lprob) then
        write (flog,1200)
     +    'Will use probe in cavity-volume calculations'
      else
        write (flog,1200)
     +    'Will NOT use probe in cavity-volume calculations'
      end if
      write (flog,1204) 'Grid shrink factor                  : ',gfact
      write (flog,1202) 'Min size of secondary grid          : ',mingrd
      if (lodl) write (flog,1204)
     +  'Grid for plot files                 : ',pgrid
      write (flog,1204) 'Volume-convergence tolerance (A3)   : ',vtoler
      write (flog,1204) 'Volume-convergence tolerance (%)    : ',ptoler
      write (flog,1204) 'Primary grid spacing                : ',gsize
      write (flog,1202) 'Grid sizes                          : ',
     +  ngrid(1),ngrid(2),ngrid(3)
      write (flog,1204) 'Lower X/Y/Z limits                  : ',
     +  xmin,ymin,zmin
      write (flog,1204) 'Upper X/Y/Z limits                  : ',
     +  xmax,ymax,zmax
      if (lpdot) then
        if (apdot .eq. 'D') then
          write (flog,1200)
     +      'Creating ODL file of protein surface as dots'
        else if (apdot .eq. 'N') then
          write (flog,1200)
     +      'Creating New-EZD file of protein surface'
        else
          write (flog,1200)
     +      'Creating Old-EZD file of protein surface'
        end if
        write (flog,1204) 'Plot grid spacing                   : ',dgrid
        write (flog,1200) 'Plot file name                      : '//
     +    (dfile(1:leng1(dfile)))
      else
        write (flog,1200)
     +    'No plot file for protein surface will be made'
      end if
      if (lodl) then
        write (flog,1204) 'Plot grid spacing                   : ',pgrid
        if (ptype .eq. 'C') then
          write (flog,1200)
     +      'Representation: all surface points connected'
        else if (ptype .eq. '1') then
          write (flog,1200) 'Representation: 1-sweep contour'
        else if (ptype .eq. '3') then
          write (flog,1200) 'Representation: 3-sweep contour'
        else if (ptype .eq. 'T') then
          write (flog,1200) 'Representation: tiles'
        else if (ptype .eq. 'D') then
          write (flog,1200) 'Representation: dots'
        else if (ptype .eq. 'O') then
          write (flog,1200) 'Representation: Old-EZD'
        else if (ptype .eq. 'N') then
          write (flog,1200) 'Representation: New-EZD'
        end if
      end if
c
      call flusho (flog)
c
      write (tty,1200) ' ','Data input done','---------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
 1202 format (99(1x,a,3i10:/))
 1204 format (99(1x,a,3f10.3:/))
 1210 format (/99(' ELEM ',4(5x,a2,3x,f5.2,1x)/))
 1212 format (/99(' ATOM ',4(5x,a4,1x,f5.2,1x)/))
 1214 format (/99(' SPAT ',4(1x,a8,1x,f5.2,1x)/))
 1220 format (/99(' RESI ',4(5x,a3,5x,a3)/))
 1230 format (' Min, max, cog for ',a1,' : ',3f10.3)
c
      end
c
c ===========================================================================
c
      subroutine get_nr_cav
c
      include 'cavity.incl'
c
      real fnow
c
      integer ibest,nbest,nscale
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Find number of cavities',
     +  '-----------------------',' '
c
      write (flog,1200) ' ','Determining number of cavities'
c
c ... set up detection cycles
c
      fnow = 1.00
c
      if (numcyc .eq. 1) fnow = fvdw
      nbest = 0
c
      do i=1,numcyc
c
        if (ltrace) write (tty,1200) ' ','NEW CYCLE'
        fstore (i) = fnow
c
        write (tty,*)
        call jvalut (' DETECTION CYCLE :',1,i)
        call fvalut (' Vanderwaals factor :',1,fnow)
c
        call set_up_grid (1,tty,fnow,probe,natoms,ngrid,gsize,
     +    xmin,xmax,ymin,ymax,zmin,zmax,xg,yg,zg,nxyz,nxy,
     +    xd,yd,zd,dvdw,ibuff,ierror,prot,notp,cavi,ltrace)
c
        if (ierror .eq. 0) then
          call count_cav
        else
          numcav = 0
          ierror = 0
        end if
c
        write (flog,1299) i,fnow,numcav
        call flusho (flog)
 1299   format (' ... Cycle ',i3,' | Vdw-factor ',f10.5,
     +    ' | Cavities ',i6)
c
        call jvalut (' Nr of cavities :',1,numcav)
        numhol (i) = numcav
        if (dtype .eq. 'F' .and. i.gt.1) then
          if (numhol(i-1) .gt. 0 .and.
     +        numhol(i)   .le. numhol(i-1)) then
            fnow = fstore (i-1)
            numcav = numhol (i-1)
            goto 100
          end if
        end if
c
        fnow = fnow * fvdw
c
      end do
c
      if (dtype .eq. 'A') then
        ibest = 1
        nbest = numhol (1)
        do i=2,numcyc
          if (numhol(i) .gt. nbest) then
            ibest = i
            nbest = numhol (i)
          end if
        end do
c
        fnow   = fstore (ibest)
        numcav = nbest
        goto 100
      end if
c
      if (numhol(numcyc) .eq. 0) then
        call errcon ('Sorry - no cavities found')
        ierror = -1
        return
      else if (numcyc .gt. 1) then
        call errcon ('Nr of cavities still increasing !!!')
        call errcon ('Increase nr of cycles or growth factor')
        ierror = -2
        return
      else
        fnow = fnow / fvdw
        goto 100
      end if
c
c ... OKAY
c
  100 continue
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') ' 4 CPU total/user/sys :',
     +  total,user,sys
c
c ... restore best grid
c
      if (lrefi .or. lmask) then
        call set_up_grid (1,tty,fnow,probe,natoms,ngrid,gsize,
     +    xmin,xmax,ymin,ymax,zmin,zmax,xg,yg,zg,nxyz,nxy,
     +    xd,yd,zd,dvdw,ibuff,ierror,prot,notp,cavi,ltrace)
c
      end if
c
      fused = fnow
c
      write (flog,1290) fnow,numcav
      call flusho (flog)
 1290 format (' Best Vdw-factor ',f10.5,' gives ',i6,' cavities')
c
c ... Write mask
c
      if (lmask) then
c
        write (tty,1200) 'Writing mask'
c
        nscale = 100
 6511   continue
        if (max(ngrid(1),ngrid(2),ngrid(3)).gt.(nscale-20)) then
          nscale = nscale * 2
          goto 6511
        end if
c
        write (iunit,7000,err=998)
     +    nint(xmin/gsize),nint(ymin/gsize),nint(zmin/gsize)
        write (iunit,7000,err=998)
     +    ngrid(1),ngrid(2),ngrid(3)
        write (iunit,7000,err=998) nscale,nscale,nscale
        write (iunit,7010,err=998) float(nscale)*gsize,
     +    float(nscale)*gsize,float(nscale)*gsize,90.0,90.0,90.0
c
        write (iunit,7020,err=998) (ibuff(i),i=1,nxyz)
        close (iunit)
c
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') ' 5 CPU total/user/sys :',
     +    total,user,sys
c
      end if
c
      if (ltrace) write (tty,1200) ' ',
     +  'Find number of cavities done',
     +  '----------------------------',' '
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') ' 6 CPU total/user/sys :',
     +  total,user,sys
c
      ierror = 0
      return
c
c ... write error
c
  998 continue
      call errcon ('While writing MASK file')
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') ' 7 CPU total/user/sys :',
     +  total,user,sys
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
 7000 format (8i5)
 7010 format (8f10.3)
 7020 format (40(1x,a1))
c
      end
c
c ===========================================================================
c
      subroutine count_cav
c
      include 'cavity.incl'
c
      real voltot
c
      integer j0,k0,istart,jstart,kstart,nzap,ntot
c
code ...
c
      if (ltrace) write (tty,1200) ' ','Count cavities',
     +  '--------------',' '
c
c ... do it !
c
      numcav = 0
c
c ... find a new starting point with a '1' in it
c
   10 continue
c
      do i=2,ngrid(1)-1
        do j=2,ngrid(2)-1
          j0 = (j-1)*ngrid(1) + i
          do k=2,ngrid(3)-1
            k0 = (k-1)*nxy + j0
            if (ibuff(k0) .eq. notp) then
              istart = i
              jstart = j
              kstart = k
              goto 100
            end if
          end do
        end do
      end do
c
      goto 999
c
  100 continue
      if (ltrace) then
        write (tty,1200) 'Potential new cavity'
        write (tty,'(1x,a,3i10)')   'First point = ',
     +    istart,jstart,kstart
        write (tty,'(1x,a,3f10.3)') 'Coordinates = ',
     +    xg(istart),yg(jstart),zg(kstart)
      end if
c
      nzap = 0
      call reczap (istart,jstart,kstart,nzap,notp,cavi)
c
      if (ltrace) call jvalut (' Nr of points "zapped" :',1,nzap)
      if (nzap .lt. minsiz) then
        if (ltrace) write (tty,1200) 'NOT a "real" cavity !'
        goto 200
      end if
c
c ... it is a cavity; process it
c
      numcav = numcav + 1
      volcav (numcav) = float(nzap)*volppt
      nptcav (numcav) = nzap
      call jvalut (' Cavity found ! Nr of points :',1,nzap)
      if (ltrace) call rvalut (' Approximate volume (A3) :',1,
     +  volcav(numcav))
c
c ... clean up and find next
c
  200 continue
      do i=1,nxyz
        if (ibuff(i) .eq. cavi) ibuff(i) = prot
      end do
c
      goto 10
c
c ... finished !
c
  999 continue
c
      call jvalut (' Nr of cavities found :',1,numcav)
c
      if (numcav .gt. 0) then
        voltot = 0.0
        ntot = 0
        do i=1,numcav
          voltot = voltot + volcav (i)
          ntot   = ntot   + nptcav (i)
        end do
c
        call jvalut (' Total nr of points in cavities : ',1,ntot)
        call rvalut (' Total cavity volume (A3)       : ',1,voltot)
      end if
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') ' 8 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ','Count cavities done',
     +  '-------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine ref_cav_vol
c
      include 'cavity.incl'
c
      real voltot,xh(maxgrd),yh(maxgrd),zh(maxgrd),oldvol,newvol
      real xlo,ylo,zlo,xhi,yhi,zhi,xxlo,xxhi,yylo,yyhi,zzlo,zzhi
      real voxvol,hsize,xtot,ytot,ztot,xzap,xdum
      real ave,sdv,vmin,vmax,vtot,xprob,yprob
c
      integer j0,k0,istart,jstart,kstart,nzap,ntot,iref,nex
      integer ilo,jlo,klo,ihi,jhi,khi,mgrid(3),mxy,mxyz,leng1
      integer length,sav1,sav2,sav3,icnt,i1,i2,i3,j1,j2,j3
c
      logical conver
c
      character prev*11
c
code ...
c
      if (ltrace) write (tty,1200) ' ','Refine cavity volumes',
     +  '---------------------',' '
c
      write (flog,1200) ' ','Refining cavity volumes'
c
c ... do it !
c
      numcav = 0
cc      margin = (fused - 1.0) * vdwmax
      margin = fused * vdwmax
      if (vtype .eq. 'O') margin = 2.0 * margin
c
      call fvalut (' Vanderwaals factor employed :',1,fused)
      call fvalut (' Max Vanderwaals radius used :',1,vdwmax)
      call fvalut (' Margin for adding to grid   :',1,margin)
c
c ... find a new starting point with a '1' in it
c
   10 continue
c
      do i=2,ngrid(1)-1
        do j=2,ngrid(2)-1
          j0 = (j-1)*ngrid(1) + i
          do k=2,ngrid(3)-1
            k0 = (k-1)*nxy + j0
            if (ibuff(k0) .eq. notp) then
              istart = i
              jstart = j
              kstart = k
              goto 100
            end if
          end do
        end do
      end do
c
      goto 999
c
  100 continue
      if (ltrace) then
        write (tty,1200) 'Potential new cavity'
        write (tty,'(1x,a,3i10)')   'First point = ',
     +    istart,jstart,kstart
        write (tty,'(1x,a,3f10.3)') 'Coordinates = ',
     +    xg(istart),yg(jstart),zg(kstart)
      end if
c
      nzap = 0
      ilo = istart
      jlo = jstart
      klo = kstart
      ihi = istart
      jhi = jstart
      khi = kstart
      call re3zap (istart,jstart,kstart,nzap,notp,cavi,
     +  ilo,jlo,klo,ihi,jhi,khi)
c
      if (ltrace) call jvalut (' Nr of points "zapped" :',1,nzap)
      if ( (.not.lspec) .and. nzap .lt. minsiz) then
        if (ltrace) write (tty,1200) 'NOT a "real" cavity !'
        goto 200
      end if
c
c ... it is a cavity; process it
c
      numcav = numcav + 1
      volcav (numcav) = float(nzap)*volppt
      nptcav (numcav) = nzap
      write (tty,'(//1x,70a//)') ('-',i=1,70)
      call jvalut (' NEW CAVITY ! Nr of points :',1,nzap)
      call rvalut (' Zero-order volume (A3) :',1,volcav(numcav))
c
      write (flog,1200) ' ','Cavity'
      write (flog,1210) numcav,xg(istart),yg(jstart),zg(kstart)
      call flusho (flog)
 1210 format (' # ',i6,' Starts at ',3f10.3)
c
      if (ltrace) then
        write (tty,'(1x,a,3i6)') 'Lowest  i/j/k ',ilo,jlo,klo
        write (tty,'(1x,a,3i6)') 'Highest i/j/k ',ihi,jhi,khi
      end if
c
      xlo = xg(ilo) - margin - probe
      ylo = yg(jlo) - margin - probe
      zlo = zg(klo) - margin - probe
      xhi = xg(ihi) + margin + probe
      yhi = yg(jhi) + margin + probe
      zhi = zg(khi) + margin + probe
c
cxyz      xlo = xg(ilo)
cxyz      ylo = yg(jlo)
cxyz      zlo = zg(klo)
cxyz      xhi = xg(ihi)
cxyz      yhi = yg(jhi)
cxyz      zhi = zg(khi)
c
      write (flog,'(a,6f10.3)')
     +  ' Cavity box ',xlo,xhi,ylo,yhi,zlo,zhi
c
      if (ltrace) then
        write (tty,'(1x,a,3f10.3)')
     +    'Lowest  x/y/z ',xg(ilo),yg(jlo),zg(klo)
        write (tty,'(1x,a,3f10.3)')
     +    'Highest x/y/z ',xg(ihi),yg(jhi),zg(khi)
        write (tty,'(a,6f10.3)')
     +    ' BOX ',xlo,xhi,ylo,yhi,zlo,zhi
      end if
c
      hsize = gsize
      oldvol = 0.0
      newvol = 0.0
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') ' 9 CPU total/user/sys :',
     +  total,user,sys
c
      write (flog,1218)
      call flusho (flog)
c
      icnt = 0
      do iref=1,numref
c
        if (ltrace) write (tty,1200) ' ','NEW CYCLE'
        write (tty,*)
        call jvalut (' REFINEMENT CYCLE :',1,iref)
c
 1999   continue
        call fvalut (' Grid spacing :',1,hsize)
c
c ... define secondary grid
c
        xxlo = xlo
        yylo = ylo
        zzlo = zlo
        xxhi = xhi
        yyhi = yhi
        zzhi = zhi
c
        call def_grd (hsize,xxlo,xxhi,yylo,yyhi,zzlo,zzhi,
     +                mgrid(1),mgrid(2),mgrid(3))
c
        call jvalut (' Nr of secondary grid points :',3,mgrid)
        if (ltrace) write (tty,'(a,6f10.3)')
     +    ' BOX ',xxlo,xxhi,yylo,yyhi,zzlo,zzhi
c
cc        write (flog,'(a,6f10.3)') ' BOX ',xxlo,xxhi,yylo,yyhi,zzlo,zzhi
c
        mxy  = mgrid(1)*mgrid(2)
        mxyz = mgrid(1)*mgrid(2)*mgrid(3)
c
        if (mgrid(1) .lt. mingrd .or.
     +      mgrid(2) .lt. mingrd .or.
     +      mgrid(3) .lt. mingrd) then
          call errcon ('Grid too course')
          hsize = hsize * gfact
          goto 1999
        end if
c
        if (mxyz     .gt. maxbuf .or.
     +      mgrid(1) .gt. maxgrd .or.
     +      mgrid(2) .gt. maxgrd .or.
     +      mgrid(3) .gt. maxgrd ) then
c
          call errcon ('Grid too big for buffer')
          if (iref.eq.1) goto 200
c
          mgrid (1) = sav1
          mgrid (2) = sav2
          mgrid (3) = sav3
          mxy       = sav1*sav2
          mxyz      = mxy *sav3
c
          goto 800
        end if
c
        sav1 = mgrid(1)
        sav2 = mgrid(2)
        sav3 = mgrid(3)
c
        do i=1,mgrid(1)
          xh (i) = xxlo + float(i-1)*hsize
        end do
        do i=1,mgrid(2)
          yh (i) = yylo + float(i-1)*hsize
        end do
        do i=1,mgrid(3)
          zh (i) = zzlo + float(i-1)*hsize
        end do
c
        voxvol = (xxhi-xxlo)*(yyhi-yylo)*(zzhi-zzlo) /
     +    float( (mgrid(1)-1)*(mgrid(2)-1)*(mgrid(3)-1) )
        call rvalut (' Volume per voxel (A3) :',1,voxvol)
c
c ... set up grid values
c
        if (lprob) then
          call set_up_grid (0,tty,1.0,probe,natoms,mgrid,hsize,
     +      xxlo,xxhi,yylo,yyhi,zzlo,zzhi,xh,yh,zh,mxyz,mxy,
     +      xd,yd,zd,dvdw,ibuff2,ierror,prot,notp,cavi,ltrace)
        else
          call set_up_grid (0,tty,1.0,0.0,natoms,mgrid,hsize,
     +      xxlo,xxhi,yylo,yyhi,zzlo,zzhi,xh,yh,zh,mxyz,mxy,
     +      xd,yd,zd,dvdw,ibuff2,ierror,prot,notp,cavi,ltrace)
        end if
c
c ... find a point that is part of this cavity
c
        do i=mgrid(1)/2,mgrid(1)-1
          do j=mgrid(2)/2,mgrid(2)-1
            j0 = (j-1)*mgrid(1) + i
            do k=mgrid(3)/2,mgrid(3)-1
c
              k0 = (k-1)*mxy + j0
              if (ibuff2(k0) .eq. notp) then
                istart = i
                jstart = j
                kstart = k
                goto 400
              end if
c
              if ((mgrid(3)-k) .ge. 2) then
                k0 = (mgrid(3)-k)*mxy + j0
                if (ibuff2(k0) .eq. notp) then
                  istart = i
                  jstart = j
                  kstart = mgrid(3)-k+1
                  goto 400
                end if
              end if
c
            end do
          end do
        end do
c
        do i=mgrid(1)/2,2,-1
          do j=mgrid(2)/2,2,-1
            j0 = (j-1)*mgrid(1) + i
            do k=mgrid(3)/2,mgrid(3)-1
c
              k0 = (k-1)*mxy + j0
              if (ibuff2(k0) .eq. notp) then
                istart = i
                jstart = j
                kstart = k
                goto 400
              end if
c
              if ((mgrid(3)-k) .ge. 2) then
                k0 = (mgrid(3)-k)*mxy + j0
                if (ibuff2(k0) .eq. notp) then
                  istart = i
                  jstart = j
                  kstart = mgrid(3)-k+1
                  goto 400
                end if
              end if
c
            end do
          end do
        end do
c
        call errcon ('??? No cavity points found in secondary grid ???')
        goto 200
c
c ... okay, found a point of the cavity
c
  400   continue
c
        if (ltrace) then
          write (tty,'(1x,a,3i10,1x,a1)')
     +      'First = ',istart,jstart,kstart,
     +      ibuff2(istart+(jstart-1)*mgrid(1)+(kstart-1)*mxy)
c
          call gkdcpu (total,user,sys)
          write (tty,'(1x,a,3f10.1)') '10 CPU total/user/sys :',
     +      total,user,sys
        end if
c
        nzap = 0
        call re5zap (istart,jstart,kstart,nzap,notp,cavi,
     +    mgrid(1),mgrid(2),mgrid(3),mxy)
        call jvalut (' Nr of points "zapped" :',1,nzap)
c
c ... calculate solvent-occupied cavity volume by adding all points which
c     are within the probe radius from the accessible cavity to the cavity
c
        if (vtype .eq. 'O') then
c
          call gkdcpu (total,user,sys)
          write (tty,'(1x,a,3f10.1)') '35 CPU total/user/sys :',
     +      total,user,sys
c
          write (tty,*) 'Rolling probe over cavity surface ...'
c
          call add_probe (ibuff2,mgrid,probe,hsize,cavi,temp,nex)
c
          nzap = nzap + nex
          call jvalut (' Nr of points "zapped" now :',1,nzap)
c
        end if
c
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '11 CPU total/user/sys :',
     +    total,user,sys
c
        volcav (numcav) = float(nzap)*voxvol
        call rvalut (' Volume (A3) :',1,volcav(numcav))
        call sphere_volume (volcav(numcav))
c
        icnt = icnt + 1
        volums (icnt) = volcav (numcav)
c
        write (flog,1220) iref,hsize,voxvol,nzap,volcav(numcav)
        call flusho (flog)
 1218   format (' Cycle','       Grid','V(vox) (A3)','       #Vox',
     +          ' Cavity vol (A3)')
 1220   format (1x,i5,1x,f10.5,1x,f10.5,1x,i10,1x,1pe15.5)
c
        if (iref.gt.1) then
          if (conver(oldvol,volcav(numcav),
     +               vtoler,ptoler,vdif(1),vdif(2))) then
            write (tty,*)
            write (tty,*) '>>> CONVERGENCE <<<'
            write (tty,*)
            call rvalut (' Last change (A3/%) :',2,vdif)
            goto 800
          end if
        end if
c
        oldvol = volcav(numcav)
        hsize  = gfact * hsize
c
      end do
c
c ... determine centre-of-gravity
c
  800 continue
      if (oldvol .le. 1.0e-6) goto 200
      if (nzap .le. 0) goto 200
c
      xtot = 0.0
      ytot = 0.0
      ztot = 0.0
      do i=1,mgrid(1)
        do j=1,mgrid(2)
          j0 = (j-1)*mgrid(1) + i
          do k=1,mgrid(3)
            k0 = (k-1)*mxy + j0
            if (ibuff2(k0) .eq. cavi) then
              xtot = xtot + xh(i)
              ytot = ytot + yh(j)
              ztot = ztot + zh(k)
            end if
          end do
        end do
      end do
      xzap = 1.0 / float (nzap)
      xtot = xtot * xzap
      ytot = ytot * xzap
      ztot = ztot * xzap
      write (tty, '(a,3f10.3)')
     +  ' Centre of cavity gravity ',xtot,ytot,ztot
      write (flog,'(a,3f10.3)')
     +  ' Centre of cavity gravity ',xtot,ytot,ztot
c
c ... compute average and sigma
c
      if (icnt .gt. 1) then
        call xstats (volums,icnt,ave,sdv,vmin,vmax,vtot)
        call jvalut (' Nr of volume calculations :',1,icnt)
        call rvalut (' Average volume       (A3) :',1,ave)
        call sphere_volume (ave)
        call rvalut (' Standard deviation   (A3) :',1,sdv)
        write (flog,'(a,i3,2f10.3)')
     +  ' Nr of calcns/average/sigma volume (A3) ',icnt,ave,sdv
      end if
c
      if (ltrace) then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '12 CPU total/user/sys :',
     +    total,user,sys
      end if
c
c ... open O macro file
c
      if (lodl) then
c
        write (line,'(a,a,i8,a)')
     +    ofile(1:leng1(ofile)),'_',numcav,'.omac'
        call remspa (line)
        close (junit)
        call xopxua (junit,line,xinter(),ierror)
c
        call stamp (line)
        write (junit,'(9a)') '! ',line(1:leng1(line))
        write (junit,'(a)') '!'
        write (junit,'(a)') 'mol #Cavity-molecule name ?#'
        write (junit,'(a)') '!'
c
        write (line,*) 'centre_xyz ',xtot,ytot,ztot
        call pretty (line)
        write (junit,'(a)') line(1:leng1(line))
        write (junit,'(a)') '!'
c
        write (line,'(a,a,i8)') 'ins','_',numcav
        call remspa (line)
        write (junit,'(9a)') 'obj ',line(1:leng1(line))
c
      end if
c
c ... find out which rejected atoms are inside the cavity
c
      if (nrejec .gt. 0) then
c
        write (tty,*)
        write (tty,*) 'Looking for atoms inside this cavity ...'
        j0 = 0
        prev = ' '
c
        do i=1,nrejec
          i1 = int (xr(i)-xxlo)/hsize
          i2 = int (yr(i)-yylo)/hsize
          i3 = int (zr(i)-zzlo)/hsize
c
c ... check 3*3*3 grid points near atom's position
c
          do j1=i1-1,i1+1
            do j2=i2-1,i2+1
              do j3=i3-1,i3+1
                if (j1.ge.1 .and. j1.le.mgrid(1) .and.
     +              j2.ge.1 .and. j2.le.mgrid(2) .and.
     +              j3.ge.1 .and. j3.le.mgrid(3)) then
                  k0 = j1 + mgrid(1)*(j2-1) + (j3-1)*mxy
                  if (ibuff2(k0) .eq. cavi) then
                    write (tty,6983)  rejnam(i),xr(i),yr(i),zr(i)
                    write (flog,6983) rejnam(i),xr(i),yr(i),zr(i)
                    j0 = j0 + 1
c
c ... if a new residue, add instruction to O macro
c
                    if (.not. lodl) goto 6980
                    if (prev .eq. rejnam(i)(5:15)) goto 6980
c
                    prev = rejnam(i)(5:15)
                    line = rejnam(i)(9:15)
                    call remspa(line)
                    write (junit,'(a,a,a,a)') 'zone ',
     +                line(1:leng1(line)),' ',line(1:leng1(line))
c
                    goto 6980
                  end if
                end if
              end do
            end do
          end do
 6980     continue
        end do
c
        if (lodl) then
          write (junit,'(a)') 'end_object'
          write (junit,'(a)') '!'
        end if
c
        if (j0 .eq. 0) then
          write (tty,*) 'There were none ...'
        else
          call jvalut (' Nr found :',1,j0)
        end if
c
        if (ltrace) then
          call gkdcpu (total,user,sys)
          write (tty,'(1x,a,3f10.1)') '35 CPU total/user/sys :',
     +      total,user,sys
        end if
c
      end if
c
 6983 format (' Inside > ',a,' @ ',3f10.3)
c
c ... find out which atoms line the surface of the cavity
c
      write (tty,*)
      write (tty,*) 'Looking for atoms lining this cavity ...'
c
      if (lodl) then
c
        write (line,'(a,a,i8)') 'lin','_',numcav
        call remspa (line)
        write (junit,'(9a)') 'obj ',line(1:leng1(line))
c
      end if
c
c ... only check atoms inside or near the box around this cavity
c
      j0 = 0
      do i=1,natoms
        linecs (i) = .false.
        test (i) = .false.
        if (xd(i) .lt. xh(1)-dvdw(i)) goto 6975
        if (yd(i) .lt. yh(1)-dvdw(i)) goto 6975
        if (zd(i) .lt. zh(1)-dvdw(i)) goto 6975
        if (xd(i) .gt. xh(mgrid(1))+dvdw(i)) goto 6975
        if (yd(i) .gt. yh(mgrid(2))+dvdw(i)) goto 6975
        if (zd(i) .gt. zh(mgrid(3))+dvdw(i)) goto 6975
        test (i) = .true.
        j0 = j0 + 1
 6975   continue
      end do
      call jvalut (' Nr of candidates :',1,j0)
      if (j0 .le. 0) goto 6970
c
      xprob = 0.0
      yprob = 0.0
      if (lprob) then
        xprob = probe*probe
        yprob = probe
      end if
c
      prev = ' '
c
      do i1=2,mgrid(1)-1
        do i2=2,mgrid(2)-1
          do i3=2,mgrid(3)-1
            k0 = i1 + mgrid(1)*(i2-1) + (i3-1)*mxy
            if (ibuff2(k0) .eq. cavi) then
ccc         print *,i1,i2,i3
              do j1=i1-1,i1+1
                do j2=i2-1,i2+1
                  do j3=i3-1,i3+1
                    j0 = j1 + mgrid(1)*(j2-1) + (j3-1)*mxy
                    if (ibuff2(j0) .ne. cavi) then
c
                      xtot = xh(j1)
                      ytot = yh(j2)
                      ztot = zh(j3)
ccc                   print *,xtot,ytot,ztot
c
                      do i=1,natoms
                        if (test(i) .and. (.not.linecs(i))) then
c
c ... XZAP = dist(atom,grid_point)
c
                          xzap = (xd(i)-xtot)**2 +
     +                           (yd(i)-ytot)**2 +
     +                           (zd(i)-ztot)**2
c
c ... XDUM = (Rad + Prob)**2 = Rad**2 + 2*Rad*prob + Prob**2
c
                          xdum = dvdw2(i) +
     +                           2.0*dvdw(i)*yprob +
     +                           xprob
c
ccc   if (xzap.lt.10.0) print *,i,xd(i),yd(i),zd(i),xzap,dvdw2(i)
c
                          if (xzap .le. xdum) then
                            linecs (i) = .true.
                            test (i) = .false.
                          end if
c
                        end if
                      end do
c
                    end if
c
                  end do
                end do
              end do
c
            end if
c
          end do
        end do
      end do
c
      j0 = 0
      do i=1,natoms
        if (linecs(i)) then
          write (tty,6993)  atmnam(i),xd(i),yd(i),zd(i)
          write (flog,6993) atmnam(i),xd(i),yd(i),zd(i)
          j0 = j0 + 1
c
c ... if a new residue, add instruction to O macro
c
          if (lodl) then
            if (prev .ne. atmnam(i)(5:15)) then
              prev = atmnam(i)(5:15)
              line = atmnam(i)(9:15)
              call remspa(line)
              write (junit,'(a,a,a,a)') 'zone ',
     +          line(1:leng1(line)),' ',line(1:leng1(line))
            end if
          end if
c
        end if
      end do
c
 6970 continue
c
      if (lodl) then
        write (junit,'(a)') 'end_object'
        write (junit,'(a)') '!'
        write (junit,'(a)') 'bell'
        write (junit,'(a)') 'message Done'
        write (junit,'(a)') '!'
        close (junit)
      end if
c
      if (j0 .eq. 0) then
        write (tty,*) 'There were none ... ??? !!!'
      else
        call jvalut (' Nr found :',1,j0)
      end if
c
      if (ltrace) then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '36 CPU total/user/sys :',
     +    total,user,sys
      end if
c
 6993 format (' Lining > ',a,' @ ',3f10.3)
c
c ... draw the bastard (clumsily; ODL file)
c
      if (.not. lodl) goto 200
c
      write (tty,*)
c
      if (ptype .eq. 'O') then
        write (line,'(a,a,i8,a)')
     +    ofile(1:leng1(ofile)),'_',numcav,'.oldezd'
      else if (ptype .eq. 'N') then
        write (line,'(a,a,i8,a)')
     +    ofile(1:leng1(ofile)),'_',numcav,'.ezd'
      else
        write (line,'(a,a,i8,a)')
     +    ofile(1:leng1(ofile)),'_',numcav,'.o'
      end if
      call remspa (line)
      close (junit)
      call xopxua (junit,line,xinter(),ierror)
c
      if (ptype .ne. 'O' .and. ptype .ne. 'N') then
        write (line,'(a,a,i8)') 'cav','_',numcav
        call remspa (line)
        write (junit,'(9a)') 'begin ',line(1:leng1(line))
        write (junit,'(a)') 'colour 10000000'
        write (junit,'(a,3f10.3,1x,a)') 't ',xtot,ytot,ztot,
     +    line(1:leng1(line))
        write (junit,'(a)') 'colour 16799999'
      end if
c
c ... first find the cavity on the requested grid
c
      xxlo = xlo
      yylo = ylo
      zzlo = zlo
      xxhi = xhi
      yyhi = yhi
      zzhi = zhi
c
      call def_grd (pgrid,xxlo,xxhi,yylo,yyhi,zzlo,zzhi,
     +              mgrid(1),mgrid(2),mgrid(3))
c
      call jvalut (' Nr of plot grid points :',3,mgrid)
      mxy  = mgrid(1)*mgrid(2)
      mxyz = mgrid(1)*mgrid(2)*mgrid(3)
c
      if (mxyz     .gt. maxbuf .or.
     +    mgrid(1) .gt. maxgrd .or.
     +    mgrid(2) .gt. maxgrd .or.
     +    mgrid(3) .gt. maxgrd ) then
        call errcon ('Plot grid too big for buffer')
        goto 200
      end if
c
      do i=1,mgrid(1)
        xh (i) = xxlo + float(i-1)*pgrid
      end do
      do i=1,mgrid(2)
        yh (i) = yylo + float(i-1)*pgrid
      end do
      do i=1,mgrid(3)
        zh (i) = zzlo + float(i-1)*pgrid
      end do
c
cc      call gkdcpu (total,user,sys)
cc      write (tty,'(1x,a,3f10.1)') '13 CPU total/user/sys :',
cc     +  total,user,sys
c
c ... set up grid values
c
      if (lprob) then
        call set_up_grid (0,tty,1.0,probe,natoms,mgrid,pgrid,
     +    xxlo,xxhi,yylo,yyhi,zzlo,zzhi,xh,yh,zh,mxyz,mxy,
     +    xd,yd,zd,dvdw,ibuff2,ierror,prot,notp,cavi,ltrace)
      else
        call set_up_grid (0,tty,1.0,0.0,natoms,mgrid,pgrid,
     +    xxlo,xxhi,yylo,yyhi,zzlo,zzhi,xh,yh,zh,mxyz,mxy,
     +    xd,yd,zd,dvdw,ibuff2,ierror,prot,notp,cavi,ltrace)
      end if
c
c ... find a point that is part of this cavity
c
      do i=mgrid(1)/2,mgrid(1)
        do j=mgrid(2)/2,mgrid(2)
          j0 = (j-1)*mgrid(1) + i
          do k=mgrid(3)/2,mgrid(3)
c
            k0 = (k-1)*mxy + j0
            if (ibuff2(k0) .eq. notp) then
              istart = i
              jstart = j
              kstart = k
              goto 1400
            end if
c
            k0 = (mgrid(3)-k)*mxy + j0
            if (ibuff2(k0) .eq. notp) then
              istart = i
              jstart = j
              kstart = mgrid(3)-k+1
              goto 1400
            end if
c
          end do
        end do
      end do
c
      call errcon ('No cavity point found on plot grid ???')
      goto 839
c
 1400 continue
c
cc      call gkdcpu (total,user,sys)
cc      write (tty,'(1x,a,3f10.1)') '14 CPU total/user/sys :',
cc     +  total,user,sys
c
      nzap = 0
      call re5zap (istart,jstart,kstart,nzap,notp,cavi,
     +  mgrid(1),mgrid(2),mgrid(3),mxy)
      call jvalut (' Nr of points "zapped" :',1,nzap)
c
      if (vtype .eq. 'O') then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '37 CPU total/user/sys :',
     +    total,user,sys
        write (tty,*) 'Rolling probe over cavity surface ...'
        call add_probe (ibuff2,mgrid,probe,pgrid,cavi,temp,nex)
        nzap = nzap + nex
        call jvalut (' Nr of points "zapped" now :',1,nzap)
      end if
c
      voxvol = (xxhi-xxlo)*(yyhi-yylo)*(zzhi-zzlo) /
     +    float( (mgrid(1)-1)*(mgrid(2)-1)*(mgrid(3)-1) )
c
      call jvalut (' Nr of points "zapped" for plot   :',1,nzap)
      call rvalut (' Volume per voxel (A3)            :',1,voxvol)
      call rvalut (' Cavity volume on plot grid (A3)  :',1,
     +  (voxvol*float(nzap)))
      call sphere_volume (voxvol*float(nzap))
      write (flog,'(a,i10,f10.7,f10.3)')
     +  ' Nr voxels/Voxel volume/Cavity volume on plot grid (A3) ',
     +  nzap,voxvol,(voxvol*float(nzap))
c
      if (ptype .eq. 'C') then
c
        call odl_con_surf (mgrid,xh,yh,zh)
c
      else if (ptype .eq. '1') then
c
        call odl_1sw_cont (mgrid,xh,yh,zh)
c
      else if (ptype .eq. '3') then
c
        call odl_3sw_cont (mgrid,xh,yh,zh)
c
      else if (ptype .eq. 'D') then
c
        call odl_dots (mgrid,xh,yh,zh)
c
      else if (ptype .eq. 'O') then
c
        call ezd_write (mgrid,xh,yh,zh)
c
      else if (ptype .eq. 'N') then
c
        call new_ezd_write (mgrid,xh,yh,zh)
c
      else if (ptype .eq. 'T') then
c
        call errcon ('Tiles not implemented yet')
c
      end if
c
  839 continue
c
      if (ptype .ne. 'O' .and. ptype .ne. 'N') then
        write (junit,'(a)') 'end_object'
      end if
c
c ... clean up primary buffer and find next cavity
c
  200 continue
      close (junit)
      do i=1,nxyz
        if (ibuff(i) .eq. cavi) ibuff(i) = prot
      end do
c
      if (ltrace) then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '15 CPU total/user/sys :',
     +    total,user,sys
      end if
c
      goto 10
c
c ... finished !
c
  999 continue
c
      write (tty,1200) ' ','Summary :'
      call jvalut (' Nr of cavities found :',1,numcav)
c
      if (numcav .gt. 0) then
        voltot = 0.0
        ntot = 0
        do i=1,numcav
          voltot = voltot + volcav (i)
          ntot   = ntot   + nptcav (i)
        end do
c
        call jvalut (
     +    ' Nr of original grid points in cavities : ',1,ntot)
        call rvalut (
     +    ' Total cavity volume (A3)               : ',1,voltot)
      end if
c
cc      call gkdcpu (total,user,sys)
cc      write (tty,'(1x,a,3f10.1)') '16 CPU total/user/sys :',
cc     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Refine cavity volumes done',
     +  '--------------------------',' '
c
      write (flog,1200) ' ',' *** VOIDOO ***',' '
      call flusho (flog)
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine calc_volume
c
      include 'cavity.incl'
c
      real hsize,newvol,oldvol,xlo,xhi,ylo,yhi,zlo,zhi
      real ave,sdv,vmin,vmax,vtot
c
      integer iref,icnt
c
      logical conver
c
code ...
c
      if (ltrace) write (tty,1200) ' ','Calculate volume',
     +  '----------------',' '
c
c ... set up refinement cycles
c
      hsize = gsize
c
      xlo = xmin
      xhi = xmax
      ylo = ymin
      yhi = ymax
      zlo = zmin
      zhi = zmax
c
      icnt = 0
      do iref=1,numref
c
        if (ltrace) write (tty,1200) ' ','NEW CYCLE'
c
        write (tty,*)
        call jvalut (' CYCLE :',1,iref)
        call fvalut (' Grid spacing :',1,hsize)
c
        call set_up_grid (0,tty,1.0,probe,natoms,ngrid,hsize,
     +    xlo,xhi,ylo,yhi,zlo,zhi,xg,yg,zg,nxyz,nxy,
     +    xd,yd,zd,dvdw,ibuff,ierror,prot,notp,cavi,ltrace)
c
        j = 0
        do i=1,nxyz
          if (ibuff(i) .eq. prot) j = j + 1
        end do
c
        call jvalut (' Nr of voxels in protein :',1,j)
        call rvalut (' Volume per voxel (A3)   :',1,volppt)
        newvol = volppt*float(j)
        call rvalut (' Protein volume (A3)     :',1,newvol)
        call sphere_volume (newvol)
c
        icnt = icnt + 1
        volums (icnt) = newvol
c
        if (conver(oldvol,newvol,vtoler,ptoler,
     +            vdif(1),vdif(2))) then
          write (tty,*)
          write (tty,*) '>>> CONVERGENCE <<<'
          write (tty,*)
          call rvalut (' Last change (A3/%) :',2,vdif)
          goto 999
        end if
c
        if (numref .eq.      1) goto 999
        if (iref   .eq. numref) goto 999
c
c ... set up for next cycle
c
        oldvol = newvol
        hsize  = hsize * gfact
c
        xlo = xmin
        xhi = xmax
        ylo = ymin
        yhi = ymax
        zlo = zmin
        zhi = zmax
c
        call def_grd (hsize,xlo,xhi,ylo,yhi,zlo,zhi,
     +                ngrid(1),ngrid(2),ngrid(3))
c
        call jvalut (' Nr of new grid points :',3,ngrid)
        nxy  = ngrid(1)*ngrid(2)
        nxyz = ngrid(1)*ngrid(2)*ngrid(3)
c
        if (nxyz     .gt. maxbuf .or.
     +      ngrid(1) .gt. maxgrd .or.
     +      ngrid(2) .gt. maxgrd .or.
     +      ngrid(3) .gt. maxgrd ) then
c
          call errcon ('New grid too big for buffer')
          goto 999
        end if
c
        do i=1,ngrid(1)
          xg (i) = xlo + float(i-1)*hsize
        end do
        do i=1,ngrid(2)
          yg (i) = ylo + float(i-1)*hsize
        end do
        do i=1,ngrid(3)
          zg (i) = zlo + float(i-1)*hsize
        end do
c
        volppt = (xhi-xlo)*(yhi-ylo)*(zhi-zlo) /
     +    float( (ngrid(1)-1)*(ngrid(2)-1)*(ngrid(3)-1) )
c
      end do
c
  999 continue
c
c ... compute average and sigma
c
      if (icnt .gt. 1) then
        call xstats (volums,icnt,ave,sdv,vmin,vmax,vtot)
        call jvalut (' Nr of volume calculations :',1,icnt)
        call rvalut (' Average volume       (A3) :',1,ave)
        call sphere_volume (ave)
        call rvalut (' Standard deviation   (A3) :',1,sdv)
      end if
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '17 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ','Calculate volume done',
     +  '---------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine find_spec_cav
c
      include 'cavity.incl'
c
      real fnow
c
      integer ibest,nbest
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Find specific cavity',
     +  '--------------------',' '
c
      write (flog,1200) ' ','Determining number of cavities'
c
c ... set up detection cycles
c
      fnow = 1.00
c
      if (numcyc .eq. 1) fnow = fvdw
      nbest = 0
c
      do i=1,numcyc
c
        if (ltrace) write (tty,1200) ' ','NEW CYCLE'
        fstore (i) = fnow
c
        write (tty,*)
        call jvalut (' DETECTION CYCLE :',1,i)
        call fvalut (' Vanderwaals factor :',1,fnow)
c
        call set_up_grid (1,tty,fnow,probe,natoms,ngrid,gsize,
     +    xmin,xmax,ymin,ymax,zmin,zmax,xg,yg,zg,nxyz,nxy,
     +    xd,yd,zd,dvdw,ibuff,ierror,prot,notp,cavi,ltrace)
c
        if (ierror .eq. 0) then
          call check_cav
        else
          numcav = 0
          ierror = 0
        end if
c
        write (flog,1299) i,fnow,numcav
        call flusho (flog)
 1299   format (' ... Cycle ',i3,' | Vdw-factor ',f10.5,
     +    ' | Cavities ',i6)
c
        call jvalut (' Nr of cavities :',1,numcav)
        if (numcav .eq. 1) goto 100
c
        fnow = fnow * fvdw
c
      end do
c
      call errcon ('Sorry - your cavity was not found')
      ierror = -1
      return
c
c ... OKAY
c
  100 continue
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '27 CPU total/user/sys :',
     +  total,user,sys
c
c ... restore best grid
c
      do i=1,nxyz
        if (ibuff(i) .eq. notp) ibuff (i) = prot
      end do
c
      do i=1,nxyz
        if (ibuff(i) .eq. cavi) ibuff (i) = notp
      end do
c
      fused = fnow
c
      write (flog,1290) fnow,numcav
      call flusho (flog)
 1290 format (' Best Vdw-factor ',f10.5,' gives ',i6,' cavities')
c
c ... Write mask
c
      if (lmask) then
c
        write (tty,1200) 'Writing mask'
c
        write (iunit,7000,err=998)
     +    nint(xmin/gsize),nint(ymin/gsize),nint(zmin/gsize)
        write (iunit,7000,err=998)
     +    ngrid(1),ngrid(2),ngrid(3)
        write (iunit,7000,err=998)
     +    100,100,100
        write (iunit,7010,err=998)
     +    100.0*gsize,100.0*gsize,100.0*gsize,90.0,90.0,90.0
c
        write (iunit,7020,err=998) (ibuff(i),i=1,nxyz)
        close (iunit)
c
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '28 CPU total/user/sys :',
     +    total,user,sys
c
      end if
c
      if (ltrace) write (tty,1200) ' ',
     +  'Find specific cavity done',
     +  '-------------------------',' '
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '29 CPU total/user/sys :',
     +  total,user,sys
c
      ierror = 0
      return
c
c ... write error
c
  998 continue
      call errcon ('While writing MASK file')
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '30 CPU total/user/sys :',
     +  total,user,sys
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
 7000 format (8i5)
 7010 format (8f10.3)
 7020 format (40(1x,a1))
c
      end
c
c ===========================================================================
c
      subroutine check_cav
c
      include 'cavity.incl'
c
      real disx,disxy,disxyz,dclose
c
      integer j0,k0,istart,jstart,kstart,nzap,nlook
      integer i1,i2,j1,j2,k1,k2,io,jo,ko
c
code ...
c
      if (ltrace) write (tty,1200) ' ','Check cavities',
     +  '--------------',' '
c
c ... do it !
c
      numcav = 0
c
c ... NLOOK determines how far around the user-specified
c     point the program checks for cavity points;
c     in any event, the one that is closest is selected
c
      nlook  = 5
c
      io = int ( (cavxyz(1) - xmin) / gsize ) + 1
      i1 = max (1, min (ngrid(1), io - nlook))
      i2 = max (1, min (ngrid(1), io + nlook))
      jo = int ( (cavxyz(2) - ymin) / gsize ) + 1
      j1 = max (1, min (ngrid(2), jo - nlook))
      j2 = max (1, min (ngrid(2), jo + nlook))
      ko = int ( (cavxyz(3) - zmin) / gsize ) + 1
      k1 = max (1, min (ngrid(3), ko - nlook))
      k2 = max (1, min (ngrid(3), ko + nlook))
c
      dclose = 9.999E19
      istart = -99
      jstart = -99
      kstart = -99
c
cc      write (*,'(3f10.3)') cavxyz(1),cavxyz(2),cavxyz(3)
cc      write (*,*) dclose
cc      write (*,'(7i10)') io,i1,i2
cc      write (*,'(7i10)') jo,j1,j2
cc      write (*,'(7i10)') ko,k1,k2
c
      do i=i1,i2
        disx = (i-io)**2
        do j=j1,j2
          j0 = (j-1)*ngrid(1) + i
          disxy = disx + (j-jo)**2
          do k=k1,k2
            k0 = (k-1)*nxy + j0
            if (ibuff(k0) .eq. notp) then
              disxyz = disxy + (k-ko)**2
cc        print *,i,j,k,disxyz,dclose
              if (disxyz .lt. dclose) then
                dclose = disxyz
                istart = i
                jstart = j
                kstart = k
              end if
            end if
          end do
        end do
      end do
c
      if (istart .gt. 0 .and. jstart .gt. 0 .and.
     +    kstart .gt. 0) goto 100
      call errcon ('Could not locate your cavity (yet ?)')
      goto 999
c
  100 continue
      if (ltrace) then
        write (tty,1200) 'Found your cavity'
        write (tty,'(1x,a,3i10)')   'First point = ',
     +    istart,jstart,kstart
        write (tty,'(1x,a,3f10.3)') 'Coordinates = ',
     +    xg(istart),yg(jstart),zg(kstart)
      end if
c
      nzap = 0
      call reczap (istart,jstart,kstart,nzap,notp,cavi)
c
      if (ltrace) call jvalut (' Nr of points "zapped" :',1,nzap)
c
c ... it is a cavity; process it
c
      numcav = 1
      volcav (numcav) = float(nzap)*volppt
      nptcav (numcav) = nzap
      call jvalut (' Cavity found ! Nr of points :',1,nzap)
      if (ltrace) call rvalut (' Approximate volume (A3) :',1,
     +  volcav(numcav))
c
c ... finished !
c
  999 continue
c
      call jvalut (' Nr of cavities found :',1,numcav)
c
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '31 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ','Check cavities done',
     +  '-------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine rotate_mol
c
      include 'cavity.incl'
c
      real dx,x,y,z,xn,yn,zn,rotmat(3,3)
c
      integer iseed,length,ndum,leng1
c
      character pline*80
c
code ...
c
      call gkrand (dx,0.0,0.0,-1)
c
      iseed = 0
      call jvalin (' Seed for random number generator     ?',1,iseed)
      call gkrand (dx,0.0,0.0,iseed)
c
      call ranrot (rotmat)
c
      write (tty,*)
      pline = 'in.pdb'
      call textin (' Input  PDB file name ?',pline)
      close (iunit)
      call xopxoa (iunit,pline,xinter(),ierror)
      if (ierror .ne. nix) then
        call errcon ('Could not open input PDB file')
        ierror = -1
        return
      end if
      line = pline
c
      write (tty,*)
      pline = 'out.pdb'
      call textin (' Output PDB file name ?',pline)
      close (junit)
      call xopxua (junit,pline,xinter(),ierror)
      if (ierror .ne. nix) then
        call errcon ('Could not open output PDB file')
        ierror = -1
        return
      end if
c
      ndum   = 0
      natoms = 0
      write (tty,*)
      call prompt (' Reading, rotating and writing ...')
c
      call stamp (pline)
      write (junit,'(9a)',err=998) 'REMARK ',
     +  pline(1:leng1(pline))
c
      write (junit,'(9a)',err=998) 'REMARK ',
     +  'Created by rotating atoms in file ',
     +  line(1:leng1(line))
c
      write (line,'(a,f10.7,a,f10.7,a,f10.7,a)')
     +  'REMARK Xnew = ',rotmat(1,1),' * Xold + ',
     +  rotmat(1,2),' * Yold + ',rotmat(1,3),' * Zold'
      call pretty (line)
      write (junit,'(a)',err=998)
     +  line(1:leng1(line))
c
      write (line,'(a,f10.7,a,f10.7,a,f10.7,a)')
     +  'REMARK Ynew = ',rotmat(2,1),' * Xold + ',
     +  rotmat(2,2),' * Yold + ',rotmat(2,3),' * Zold'
      call pretty (line)
      write (junit,'(a)',err=998)
     +  line(1:leng1(line))
c
      write (line,'(a,f10.7,a,f10.7,a,f10.7,a)')
     +  'REMARK Znew = ',rotmat(3,1),' * Xold + ',
     +  rotmat(3,2),' * Yold + ',rotmat(3,3),' * Zold'
      call pretty (line)
      write (junit,'(a)',err=998)
     +  line(1:leng1(line))
c
      write (line,'(a,i15,a)')
     +  'REMARK Seed ',iseed,
     +  ' used for random number generator' 
      call pretty (line)
      write (junit,'(a)',err=998)
     +  line(1:leng1(line))
c
   10 continue
      read (iunit,'(a)',err=997,end=999) line
      ndum = ndum + 1
      call upcase (line)
      if (line(1:5) .eq. 'ATOM ') goto 11
      if (line(1:6) .eq. 'HETATM') goto 11
c
      write (junit,'(a)',err=998) line(1:leng1(line))
      goto 10
c
   11 continue
      natoms = natoms + 1
      read (line(31:),'(3f8.3)') x,y,z
      xn = rotmat(1,1)*x + rotmat(1,2)*y + rotmat(1,3)*z
      yn = rotmat(2,1)*x + rotmat(2,2)*y + rotmat(2,3)*z
      zn = rotmat(3,1)*x + rotmat(3,2)*y + rotmat(3,3)*z
      write (line(31:54),'(3f8.3)') xn,yn,zn
      write (junit,'(a)',err=998) line(1:leng1(line))
c
      goto 10
c
c ... end of PDB file
c
  999 continue
c
      if (line(1:3) .ne. 'END')
     +  write (junit,'(a)',err=998) 'END'
c
      write (tty,*)
      call jvalut (' Number of lines read   :',1,ndum)
      call jvalut (' Number of atoms read   :',1,natoms)
c
      close (iunit)
      close (junit)
c
      return
c
  998 continue
      call errcon ('While writing new PDB file')
      goto 999
c
  997 continue
      call errcon ('While reading PDB file')
      goto 999
c
      end
