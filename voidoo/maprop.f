      program maprop
c
c === MAPROP === calculate grid properties
c
c Gerard Kleywegt @ 940317
c
c ===========================================================================
c
c Version history:
c ----------------
c
c Version 0.1 @ 940317 - initial version & documentation; 3 libraries
c Version 0.2 @ 940318 - removed a bug; minor changes
c Version 0.3 @ 940320 - added RESET constant; made 'radii' library
c
c ===========================================================================
c
      include 'maprop.incl'
c
code ...
c
      call gkinit (prognm,version)
c
c ... initialise program
c
      ierror = 0
c
      call prog_init
      if (ierror .ne. 0) goto 9999
c
c ... get all relevant data
c
      call get_data
      if (ierror .ne. 0) goto 9999
c
      call gkdcpu (total,user,sys)
      write (*,'(1x,a,3f10.1)') ' 1 CPU total/user/sys :',
     +  total,user,sys
c
c ... initialise map
c
      call init_map
      if (ierror .ne. 0) goto 9999
c
c ... calculate map
c
      call calc_map
      if (ierror .ne. 0) goto 9999
c
c ... create map file
c
      call put_surf
c
 9999 continue
      call gkquit
c
      stop
      end
c
c ===========================================================================
c
      subroutine prog_init
c
      include 'maprop.incl'
c
      integer i
c
code ...
c
      write (*,1200)
      write (*,1000) 'Current version',version
      write (*,1100) 'Max nr of atoms',maxatm
      write (*,1100) 'Max nr of elements',maxelm
      write (*,1100) 'Max nr of atom types',maxaty
      write (*,1100) 'Max nr of residue/atom types',maxrac
      write (*,1100) 'Max nr of residue types',maxres
      write (*,1100) 'Max nr of points for map',maxbuf
      write (*,1100) 'Ditto, for distance buffer',maxbuf
      write (*,1100) 'Memory use (Bytes) for major arrays',memuse
      write (*,1200)
c
c ... Property values
c
      do i=1,maxelm
        namelm (i) = '??'
        prorad (i) = 0.00
      end do
c
      do i=1,maxaty
        namaty (i) = '????'
        prorat (i) = 0.00
      end do
c
      do i=1,maxrac
        namrac (i) = '??? ????'
        prorac (i) = 0.00
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
      gsize  = 1.0
      cuton  = 1.0
      cutoff = 7.0
      const  = 1.0
      power  = 1
      corest = 0.0
      combi  = 'SUM '
c
      ierror = 0
      return
c
 1000 format (1x,a40,' : ',a)
 1100 format (1x,a40,' : ',i10)
 1200 format (/' ***** ',5('MAPROP ***** ')/)
 1300 format (1x,a40,' : ',f10.5)
c
      end
c
c ===========================================================================
c
      subroutine get_data
c
      include 'maprop.incl'
c
      real xat,ave,sdv,xtot,ytot,ztot,perc,margin
      real save1,save2,save3,save4,save5,save6
c
      integer ndum,n2,length,iunit,i,leng1
c
      logical xinter
c
      character residu*3,elemnt*2,pline*80,atnam*4,racnm*8
      character line*256
c
code ...
c
      iunit = 11
c
      write (*,1200) 'Data input','----------',' '
c
c ... Property values
c
      write (*,1200) '(1) Property values and residue types',' '
c
      line = 'maprop.radii'
c
c ... 950118 - check if environment variable GKLIB is defined
c
      call gklibf (line)
c
      call textin (' Library file ?',line)
      close (iunit)
      call xopxoa (iunit,line,xinter(),ierror)
      if (ierror .ne. 0) then
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
      write (*,1200) ' ','Reading your library file ...'
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
          read (line(5:),*) namelm(numelm),prorad(numelm)
        else
          call errcon ('Too many elements')
          call textut (' SKIPPED :',line(5:))
        end if
      else if (line(1:4) .eq. 'ATOM') then
        if (numaty .lt. maxaty) then
          numaty = numaty + 1
          read (line(5:),*) namaty(numaty),prorat(numaty)
        else
          call errcon ('Too many atom types')
          call textut (' SKIPPED :',line(5:))
        end if
      else if (line(1:4) .eq. 'SPAT') then
        if (numrac .lt. maxrac) then
          numrac = numrac + 1
          read (line(5:),*) namrac(numrac),prorac(numrac)
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
      if (numelm .gt. 0) write (*,1210)
     +  (namelm(i),prorad(i),i=1,numelm)
      if (numaty .gt. 0) write (*,1212)
     +  (namaty(i),prorat(i),i=1,numaty)
      if (numrac .gt. 0) write (*,1214)
     +  (namrac(i),prorac(i),i=1,numrac)
      if (numres .gt. 0) write (*,1220) (resnam(i),i=1,numres)
c
      if (numelm + numaty + numrac .le. 0) then
        call errcon (' No property values defined')
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
      write (*,1200) ' ','(2) PDB file',' '
      pline = 'in.pdb'
      call textin (' PDB file name ?',pline)
      close (iunit)
      call xopxoa (iunit,pline,xinter(),ierror)
      if (ierror .ne. 0) then
        call errcon ('Could not open PDB file')
        ierror = -1
        return
      end if
c
      ndum   = 0
      natoms = 0
      nrejec = 0
      n2 = numres
      write (*,1200) ' ','Reading your PDB file ...'
c
   10 continue
      read (iunit,'(a)',end=999) line
      call upcase (line)
      if (line(1:6) .eq. 'ATOM  ') goto 11
      if (line(1:6) .eq. 'HETATM') goto 11
      if (line(1:6) .eq. 'REMARK') 
     +  write (*,'(1x,a)') line(1:leng1(line))
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
c ... count rejected atoms
c
      nrejec = nrejec + 1
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
c ... property defined for this atom/element ?
c
   20 continue
c
      if (numrac .gt. 0) then
        racnm = residu//' '//line(13:16)
        do i=1,numrac
          if (racnm .eq. namrac(i)) then
            natoms = natoms + 1
            dpro (natoms) = prorac (i)
            read (line(29:),*) xd(natoms),yd(natoms),zd(natoms)
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
            dpro (natoms) = prorat (i)
            read (line(29:),*) xd(natoms),yd(natoms),zd(natoms)
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
            dpro (natoms) = prorad (i)
            read (line(29:),*) xd(natoms),yd(natoms),zd(natoms)
            atmnam (natoms) = line (13:27)
            goto 10
          end if
        end do
      end if
c
c ... unknown atom/element
c
      call errcon ('Element without property value')
      call textut (' Element :',elemnt)
      call textut (' Atom    :',line(1:28))
      ierror = -1
      return
c
c ... end of PDB file
c
  999 continue
c
      close (iunit)
      call jvalut (' Number of atoms read     :',1,ndum)
      call jvalut (' Number of atoms kept     :',1,natoms)
      call jvalut (' Number of atoms rejected :',1,nrejec)
c
      if (n2 .gt. numres) then
        call asciut (' Rejected residue types :',n2-numres,
     +    resnam(numres+1))
      else
        write (*,1200) 'No residue types rejected'
      end if
c
c ... enough atoms ?
c
      if (natoms .lt. 10) then
        call errcon ('Not enough atoms (min = 10)')
        ierror = -2
        return
      end if
c
c ... Parameters
c
      write (*,1200) ' ','(3) Various parameters'
c
      write (*,'(20(/1x,a))')
     +  'The following map P(X) will be calculated:',
     +  ' ',
     +  'For each atom A and each grid point X,',
     +  'R = distance(A,X):',
     +  ' ',
     +  'IF R > Cutoff => no contribution',
     +  'IF R < Cuton  => P(X) = P(A) * C',
     +  'ELSE          => P(X) = P(A) * C / (R^N)',
     +  ' ',
     +  'All points which are NOT set, will be given',
     +  'a value Q',
     +  ' ',
     +  'Values for different atoms A may be combined:',
     +  'SUM  -> sum all values from contributing atoms',
     +  'PROD -> multiply them',
     +  'MIN  -> take the minimum',
     +  'MAX  -> take the maximum',
     +  'NEAR -> only keep value of the nearest atom',
     +  'NEST -> ditto, but require at least two nearby atoms'
c
      call fvalin (
     +  ' Minimal radius CUTON  (>0)     ?',1,cuton)
      cuton = max (0.001, cuton)
      call rvalut (
     +  ' Minimal radius CUTON           :',1,cuton)
c
      call fvalin (
     +  ' Maximal radius CUTOFF (>CUTON) ?',1,cutoff)
      cutoff = max (cuton+0.001, cutoff)
      call rvalut (
     +  ' Maximal radius CUTOFF          :',1,cutoff)
c
      call fvalin (
     +  ' Constant C (<> 0)              ?',1,const)
      const = max (0.001, const)
      call rvalut (
     +  ' Constant C                     :',1,const)
c
      call ivalin (
     +  ' Radius power N                 ?',1,power)
      call ivalut (
     +  ' Radius power N                 :',1,power)
c
      call fvalin (
     +  ' Uniform value for unset points ?',1,corest)
      call rvalut (
     +  ' Uniform value for unset points :',1,corest)
c
      call textin (
     +  ' Combine SUM/PROD/MIN/MAX/NEAR/NEST ?',combi)
      call upcase (combi)
      call remspa (combi)
      if (index ('SUM |PROD|MIN |MAX |NEAR|NEST',combi(1:4))
     +    .le. 0) then
        call errcon (' Invalid combination selected')
        combi = 'SUM '
      end if
      call textut (
     +  ' Combination method             :',combi)
c
      dfile = 'maprop.nezd'
      call textin (
     +  ' Name of map file (NEWEZD)      ?',dfile)
      call remspa (dfile)
      call textut (
     +  ' Name of map file (NEWEZD)      :',dfile)
c
c ... Box and grid
c
      write (*,1200) ' ','(4) Map grid',' '
c
      call xstats (xd,natoms,ave,sdv,xmin,xmax,xtot)
      call xstats (yd,natoms,ave,sdv,ymin,ymax,ytot)
      call xstats (zd,natoms,ave,sdv,zmin,zmax,ztot)
      xat = 1.0/float(natoms)
c
      write (*,1230) 'X',xmin,xmax,xtot*xat
      write (*,1230) 'Y',ymin,ymax,ytot*xat
      write (*,1230) 'Z',zmin,zmax,ztot*xat
c
      save1 = xmin
      save2 = xmax
      save3 = ymin
      save4 = ymax
      save5 = zmin
      save6 = zmax
c
      margin = max (3.0, min (10.0, cutoff+1.0))
c
 6996 continue
      call fvalin (' Grid spacing (A)        ?',1,gsize)
      gsize = max (0.01, gsize)
      call fvalut (' Grid spacing (A)        :',1,gsize)
      call fvalut (' Margin on all sides (A) :',1,margin)
c
      xmin = save1 - margin
      ymin = save3 - margin
      zmin = save5 - margin
      xmax = save2 + margin
      ymax = save4 + margin
      zmax = save6 + margin
c
      call def_grd (gsize,xmin,xmax,ymin,ymax,zmin,zmax,
     +              ngrid(1),ngrid(2),ngrid(3))
c
      write (*,1230) 'X',xmin,xmax
      write (*,1230) 'Y',ymin,ymax
      write (*,1230) 'Z',zmin,zmax
c
      nxy  = ngrid(1)*ngrid(2)
      nxyz = ngrid(1)*ngrid(2)*ngrid(3)
      perc = float(nxyz) * 100.0 / float(maxbuf)
      call jvalut (' Number of grid points :',3,ngrid)
      call jvalut (' Required buffer       :',1,nxyz)
      call jvalut (' Buffer size           :',1,maxbuf)
      call fvalut (' % of buffer needed    :',1,perc)
c
      if (nxyz .gt. maxbuf) then
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
      write (*,1200) ' ','Data input done','---------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
 1202 format (99(1x,a,3i10:/))
 1204 format (99(1x,a,3f10.3:/))
 1210 format (/99(' ELEM ',4(5x,a2,3x,f5.2,1x:)/))
 1212 format (/99(' ATOM ',4(5x,a4,1x,f5.2,1x:)/))
 1214 format (/99(' SPAT ',4(1x,a8,1x,f5.2,1x:)/))
 1220 format (/99(' RESI ',4(5x,a3,5x,a3:)/))
 1230 format (' Min, max, cog for ',a1,' : ',3f10.3)
c
      end
c
c ===========================================================================
c
      subroutine def_grd (gsize,xmin,xmax,ymin,ymax,zmin,zmax,
     +                    ngrid1,ngrid2,ngrid3)
c
c ... def_grd (...) - define a proper grid
c
      implicit none
c
      real gsize,xmin,xmax,ymin,ymax,zmin,zmax
c
      integer ngrid1,ngrid2,ngrid3
c
code ...
c
c ... NOTE: make the grid 2 points wider than strictly necessary !
c
      xmin  = gsize * float (int (xmin/gsize) - 1)
      xmax  = gsize * float (int (xmax/gsize) + 1)
      ymin  = gsize * float (int (ymin/gsize) - 1)
      ymax  = gsize * float (int (ymax/gsize) + 1)
      zmin  = gsize * float (int (zmin/gsize) - 1)
      zmax  = gsize * float (int (zmax/gsize) + 1)
c
      ngrid1 = 1 + nint ((xmax - xmin)/gsize)
      ngrid2 = 1 + nint ((ymax - ymin)/gsize)
      ngrid3 = 1 + nint ((zmax - zmin)/gsize)
c
      xmax = xmin + gsize * float (ngrid1-1)
      ymax = ymin + gsize * float (ngrid2-1)
      zmax = zmin + gsize * float (ngrid3-1)
c
      return
      end
c
c ===========================================================================
c
      subroutine calc_map
c
      include 'maprop.incl'
c
c ... fill the map
c
      real dx,dy,dz,radius,sqrrad,radmin,r,t,tmin,tmax
c
      integer atom,j0,k0,i1,i2,j1,j2,k1,k2,i,j,k,ip,np,ires
c
code ...
c
      write (*,*)
      call prompt (' Calculating map ...')
c
      radius = cutoff
      sqrrad = radius**2
      radmin = cuton**2
c
      nxyz = ngrid(1)*ngrid(2)*ngrid(3)
      np = 500
      ip = 0
c
      do atom=1,natoms
c
c ... find max limits of grid to check for this atom
c
        i1 = int ( (xd(atom) - radius - xmin) / gsize ) + 1
        i1 = max (1, min (ngrid(1), i1))
        i2 = int ( (xd(atom) + radius - xmin) / gsize ) + 2
        i2 = max (1, min (ngrid(1), i2))
        j1 = int ( (yd(atom) - radius - ymin) / gsize ) + 1
        j1 = max (1, min (ngrid(2), j1))
        j2 = int ( (yd(atom) + radius - ymin) / gsize ) + 2
        j2 = max (1, min (ngrid(2), j2))
        k1 = int ( (zd(atom) - radius - zmin) / gsize ) + 1
        k1 = max (1, min (ngrid(3), k1))
        k2 = int ( (zd(atom) + radius - zmin) / gsize ) + 2
        k2 = max (1, min (ngrid(3), k2))
c
cATOM      1  CB  PYR A   1      36.358  61.252  40.207  1.00 20.16   6
c123456789012345678901234567
c
        if (combi .eq. 'NEST') then
          read (atmnam(atom)(11:14),'(i4)') ires
        end if
c
cc        write (*,'(7i10)') atom,i1,i2,j1,j2,k1,k2
c
        if (combi .eq. 'SUM ') then
          do i=i1,i2
            dx = (xg(i)-xd(atom))**2
            do j=j1,j2
              dy = dx + (yg(j)-yd(atom))**2
              j0 = (j-1)*ngrid(1) + i
              do k=k1,k2
                k0 = (k-1)*nxy + j0
                  dz = dy + (zg(k)-zd(atom))**2
                  if (dz .le. radmin) then
                    rbuff (k0) = rbuff(k0) + const * dpro(atom)
                    dbuff (k0) = 0.0
                  else if (dz .le. sqrrad) then
                    t = const * dpro(atom)
                    if (power .ne. 0) then
                      r = sqrt(dz)
                      if (power .eq. 1) then
                        t = t / r
                      else if (power .eq. -1) then
                        t = t * r
                      else if (power .eq. 2) then
                        t = t / (r*r)
                      else if (power .eq. -2) then
                        t = t * r * r
                      else
                        t = t / (r**power)
                      end if
                    end if
                    rbuff(k0) = rbuff(k0) + t
                    dbuff (k0) = 0.0
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
                end if
              end do
            end do
          end do
c
        else if (combi .eq. 'PROD') then
          do i=i1,i2
            dx = (xg(i)-xd(atom))**2
            do j=j1,j2
              dy = dx + (yg(j)-yd(atom))**2
              j0 = (j-1)*ngrid(1) + i
              do k=k1,k2
                k0 = (k-1)*nxy + j0
                  dz = dy + (zg(k)-zd(atom))**2
                  if (dz .le. radmin) then
                    rbuff (k0) = rbuff(k0) * const * dpro(atom)
                    dbuff (k0) = 0.0
                  else if (dz .le. sqrrad) then
                    t = const * dpro(atom)
                    if (power .ne. 0) then
                      r = sqrt(dz)
                      if (power .eq. 1) then
                        t = t / r
                      else if (power .eq. -1) then
                        t = t * r
                      else if (power .eq. 2) then
                        t = t / (r*r)
                      else if (power .eq. -2) then
                        t = t * r * r
                      else
                        t = t / (r**power)
                      end if
                    end if
                    rbuff(k0) = rbuff(k0) * t
                    dbuff (k0) = 0.0
                  end if
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
              end do
            end do
          end do
c
        else if (combi .eq. 'MIN ') then
          do i=i1,i2
            dx = (xg(i)-xd(atom))**2
            do j=j1,j2
              dy = dx + (yg(j)-yd(atom))**2
              j0 = (j-1)*ngrid(1) + i
              do k=k1,k2
                k0 = (k-1)*nxy + j0
                  dz = dy + (zg(k)-zd(atom))**2
                  if (dz .le. radmin) then
                    t = const * dpro(atom)
                    if (t.lt.rbuff(k0)) then
                      rbuff (k0) = t
                      dbuff (k0) = 0.0
                    end if
                  else if (dz .le. sqrrad) then
                    t = const * dpro(atom)
                    if (power .ne. 0) then
                      if (power .eq. 1) then
                        t = t / r
                      else if (power .eq. -1) then
                        t = t * r
                      else if (power .eq. 2) then
                        t = t / (r*r)
                      else if (power .eq. -2) then
                        t = t * r * r
                      else
                        t = t / (r**power)
                      end if
                      t = t / (r**power)
                    end if
                    if (t.lt.rbuff(k0)) then
                      rbuff (k0) = t
                      dbuff (k0) = 0.0
                    end if
                  end if
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
              end do
            end do
          end do
c
        else if (combi .eq. 'MAX ') then
          do i=i1,i2
            dx = (xg(i)-xd(atom))**2
            do j=j1,j2
              dy = dx + (yg(j)-yd(atom))**2
              j0 = (j-1)*ngrid(1) + i
              do k=k1,k2
                k0 = (k-1)*nxy + j0
                  dz = dy + (zg(k)-zd(atom))**2
                  if (dz .le. radmin) then
                    t = const * dpro(atom)
                    if (t.gt.rbuff(k0)) then
                      rbuff (k0) = t
                      dbuff (k0) = 0.0
                    end if
                  else if (dz .le. sqrrad) then
                    t = const * dpro(atom)
                    if (power .ne. 0) then
                      r = sqrt(dz)
                      if (power .eq. 1) then
                        t = t / r
                      else if (power .eq. -1) then
                        t = t * r
                      else if (power .eq. 2) then
                        t = t / (r*r)
                      else if (power .eq. -2) then
                        t = t * r * r
                      else
                        t = t / (r**power)
                      end if
                    end if
                    if (t.gt.rbuff(k0)) then
                      rbuff (k0) = t
                      dbuff (k0) = 0.0
                    end if
                  end if
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
              end do
            end do
          end do
c
        else if (combi.eq.'NEST') then
c
c xyz
c this is not correct yet
c it requires two loops over the atoms,
c one to find the nearest atom and hence the grid value
c and one to see if there is an atom from another
c residue nearby
c
          do i=i1,i2
            dx = (xg(i)-xd(atom))**2
            do j=j1,j2
              dy = dx + (yg(j)-yd(atom))**2
              j0 = (j-1)*ngrid(1) + i
              do k=k1,k2
                k0 = (k-1)*nxy + j0
                  dz = dy + (zg(k)-zd(atom))**2
                  r = sqrt(dz)
                  if (dz .le. radmin) then
                      rbuff(k0) = const * dpro(atom)
                      dbuff(k0) = float(ires)
                      nbrcnt(k0) = -99999
                  else if (dz .le. sqrrad) then
c
                     if (nint(dbuff(k0)) .eq. ires) goto 6423
c
                      t = const * dpro(atom)
                      if (power .ne. 0) then
                        if (power .eq. 1) then
                          t = t / r
                        else if (power .eq. -1) then
                          t = t * r
                        else if (power .eq. 2) then
                          t = t / (r*r)
                        else if (power .eq. -2) then
                          t = t * r * r
                        else
                          t = t / (r**power)
                        end if
                      end if
                      rbuff(k0) = t
                      dbuff(k0) = float(ires)
                      nbrcnt(k0) = nbrcnt(k0)+1
                  end if
c
 6423             continue
c
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
              end do
            end do
          end do
c
        else if (combi.eq.'NEAR') then
          do i=i1,i2
            dx = (xg(i)-xd(atom))**2
            do j=j1,j2
              dy = dx + (yg(j)-yd(atom))**2
              j0 = (j-1)*ngrid(1) + i
              do k=k1,k2
                k0 = (k-1)*nxy + j0
                  dz = dy + (zg(k)-zd(atom))**2
                  r = sqrt(dz)
                  if (dz .le. radmin) then
                    if (r.lt.dbuff(k0)) then
                      rbuff(k0) = const * dpro(atom)
                      dbuff(k0) = r
                    end if
                  else if (dz .le. sqrrad) then
                    if (r.lt.dbuff(k0)) then
                      t = const * dpro(atom)
                      if (power .ne. 0) then
                        if (power .eq. 1) then
                          t = t / r
                        else if (power .eq. -1) then
                          t = t * r
                        else if (power .eq. 2) then
                          t = t / (r*r)
                        else if (power .eq. -2) then
                          t = t * r * r
                        else
                          t = t / (r**power)
                        end if
                      end if
                      rbuff(k0) = t
                      dbuff(k0) = r
                    end if
                  end if
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
              end do
            end do
          end do
c
        end if
c
        ip = ip + 1
        if (ip .eq. np) then
          call jvalut (' Nr of atoms done :',1,atom)
          ip = 0
        end if
c
      end do
c
      write (*,*)
      call gkdcpu (total,user,sys)
      write (*,'(1x,a,3f10.1)')
     +  ' 4 CPU total/user/sys :',total,user,sys
c
c ... for NEST, only keep if point has > 1 neighbour atom
c
      if (combi .eq. 'NEST') then
        write (*,*)
        call prompt (' Resetting points with < 2 neighbours')
        do i=1,nxyz
          if (nbrcnt(i) .lt. 2) dbuff(i)=999999.999
        end do
      end if
c
c ... reset points which haven't received a value
c
      write (*,*)
      call rvalut (' Resetting unset points to :',1,corest)
      i1 = 0
      do i=1,nxyz
        if (dbuff(i).gt.999000.0) then
          i1 = i1 + 1
          rbuff (i) = corest
        end if
      end do
      call jvalut (' Nr of points reset :',1,i1)
      call jvalut (' Out of a total of  :',1,nxyz)
c
      write (*,*)
      call xstats (rbuff,nxyz,t,r,tmin,tmax,dz)
      call jvalut (' Nr of points in grid :',1,nxyz)
      call rvalut (' Sum of values in map :',1,dz)
      call rvalut (' Average value in map :',1,t)
      call rvalut (' St. Deviation        :',1,r)
      call rvalut (' Minimum value in map :',1,tmin)
      call rvalut (' Maximum value in map :',1,tmax)
c
      write (*,*)
      call gkdcpu (total,user,sys)
      write (*,'(1x,a,3f10.1)') ' 5 CPU total/user/sys :',
     +  total,user,sys
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
      subroutine init_map
c
      include 'maprop.incl'
c
      real start
c
      integer i,j
c
code ...
c
      j = ngrid(1)*ngrid(2)*ngrid(3)
c
c ... clear NBRCNT buffer for NEST option
c
      if (combi .eq. 'NEST') then
        write (*,*)
        call prompt (' Clearing neighbour count buffer ...')
        do i=1,j
          nbrcnt (i) = 0
        end do
      end if
c
      write (*,*)
      call prompt (' Clearing nearest atom buffer ...')
      start = 999999.999
      do i=1,j
        dbuff (i) = start
      end do
c
      write (*,*)
      call gkdcpu (total,user,sys)
      write (*,'(1x,a,3f10.1)') ' 2 CPU total/user/sys :',
     +  total,user,sys
c
      if (combi .eq. 'NEAR') then
        start = 0.0
      else if (combi .eq. 'NEST') then
        start = 0.0
      else if (combi .eq. 'SUM ') then
        start = 0.0
      else if (combi .eq. 'PROD') then
        start = 1.0
      else if (combi .eq. 'MIN ') then
        start = 999999.999
      else if (combi .eq. 'MAX ') then
        start = -999999.999
      end if
c
      write (*,*)
      call rvalut (' Initialising map with value :',1,start)
      do i=1,j
        rbuff (i) = start
      end do
c
      write (*,*)
      call gkdcpu (total,user,sys)
      write (*,'(1x,a,3f10.1)') ' 3 CPU total/user/sys :',
     +  total,user,sys
c
      return
      end
c
c ===========================================================================
c
      subroutine put_surf
c
      include 'maprop.incl'
c
      integer i,j,k,junit,nscale,leng1
c
      logical xinter
c
      character myline*256
c
code ...
c
      junit = 12
c
      write (*,*)
      call prompt (' Opening New-EZD map file ...')
c
      call xopxua (junit,dfile,xinter(),ierror)
      if (ierror .ne. 0) return
c
      nxyz = ngrid(1) * ngrid(2) * ngrid(3)
c
      nscale = 100
 6511 continue
      if (max(ngrid(1),ngrid(2),ngrid(3)) .gt. (nscale-20)) then
        nscale = nscale * 2
        goto 6511
      end if
c
      write (*,*)
      call prompt (' Writing map ...')
c
      write (junit,'(a)') 'EZD_MAP'
c
      call stamp (myline)
      myline = '! '//myline
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,*) '! CUTON = ',cuton,' A'
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,*) '! CUTOFF = ',cutoff,' A'
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,*) '! CONSTANT = ',const
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,*) '! RADIUS POWER = ',power
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,*) '! COMBINATION = ',combi
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,*) '! UNSET POINTS = ',corest
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,8010) 'CELL',float(nscale)*gsize,
     +  float(nscale)*gsize,float(nscale)*gsize,90.0,90.0,90.0
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
      write (*,*) myline(1:leng1(myline))
c
      write (myline,8000) 'ORIGIN',nint(xg(1)/gsize),
     +  nint(yg(1)/gsize),nint(zg(1)/gsize)
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
      write (*,*) myline(1:leng1(myline))
c
      write (myline,8000) 'EXTENT',ngrid(1),ngrid(2),ngrid(3)
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
      write (*,*) myline(1:leng1(myline))
c
      write (myline,8000) 'GRID',nscale,nscale,nscale
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
      write (*,*) myline(1:leng1(myline))
c
      write (myline,8010) 'SCALE',1.0
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
      write (*,*) myline(1:leng1(myline))
c
      write (junit,'(a)') 'MAP'
      do i=1,nxyz,7
        j = min(nxyz,i+6)
        write (myline,8030) (rbuff(k),k=i,j)
c
c ... replace "0 " with "  " to save disk space
c     (we only do this calculation once so it's
c     wothwhile)
c
  123   k = index (myline,'0 ')
        if (k .le. 0) goto 125
        myline (k:k) = ' '
        goto 123
c
  125   continue
        call pretty (myline)
        write (junit,'(a)') myline(1:leng1(myline))
      end do
c
 8000 format (a,1x,8i5)
 8010 format (a,1x,8f10.3)
 8030 format (7(1x,f10.3))
c
c ... end
c
  999 continue
      write (*,*)
      call gkdcpu (total,user,sys)
      write (*,'(1x,a,3f10.1)') ' 6 CPU total/user/sys :',
     +  total,user,sys
c
      return
      end
