      program xplo2d
c
c ... XPLO2D - X-PLOR utilities
c
c ... compile and link with:
c
c     f77 -Olimit 3000 -v -w0 -u -check_bounds -c xplo2d.f
c     f77 -o XPLO2D xplo2d.o ~gerard/progs/gklib/kleylib
c     strip XPLO2D
c
c ... TO DO:
c     - use CONECT records (if they exist) to figure out bonds
c
      implicit none
c
      character*20 prognm, vers
      parameter (prognm = 'XPLO2D', vers = '050802/3.3.2')
c
      integer maxang, inter
      parameter (maxang = 37, inter = 10)
c
      real xvoid
      parameter (xvoid = -999.99)
c
      real pcval(maxang,maxang,maxang),pcmin,pcmax,pctot,pc
      real pcave,pcsdv,dx,levels(15)
c
      integer nfile,length,ierr,ifile,nval,ntot,i1,i2,i3,leng1
      integer ityp,ival,ilim(5,3),i,nl,ijkmax(3)
c
      logical lforce,lnoh
c
      character filenm*80,fmt*20,which*6,myline*128
c
code ...
c
      call gkinit (prognm,vers)
      which = 'auto'
c
c ... check for "-force" and "-noh" flags (used in AUTODI)
c
      lforce = .false.
      lnoh = .false.
      call gknarg (nl)
      if (nl .gt. 0) then
        do i=1,nl
          call gkgarg (i,myline,ierr)
          if (ierr .eq. 0) then
            call locase (myline)
            if (myline .eq. '-force') then
              lforce = .true.
              call prompt (' -force flag found !')
            else if (myline .eq. '-noh') then
              lnoh = .true.
              call prompt (' -noh flag found !')
            end if
          end if
        end do
      end if
 6666 continue
c
      write (*,*)
      write (*,*) 'XSEARCH - convert X-PLOR search.dat files'
      write (*,*) '          from a direct rotation search'
      write (*,*) '          into O2D 2D contour plot files'
      write (*,*)
      write (*,*) 'X1DTRAZ - convert X-PLOR xxx.3dmatrix files'
      write (*,*) '          from a 1D translation function'
      write (*,*) '          along the Z-axis (C-axis) to an'
      write (*,*) '          O2D 1D plot file'
      write (*,*)
      write (*,*) 'X1DYTRA - convert X-PLOR xxx.3dmatrix files'
      write (*,*) '          from a 1D translation function'
      write (*,*) '          along the Y-axis (B-axis) to an'
      write (*,*) '          O2D 1D plot file'
      write (*,*)
      write (*,*) 'XFILTER - analyse filter.list files from an'
      write (*,*) '          X-PLOR PC-refinement OR search.dat'
      write (*,*) '          from a direct rotation search'
      write (*,*)
      write (*,*) 'LUZZATI - convert a luzzati.list file into'
      write (*,*) '          an O2D plot file (or two, if you'
      write (*,*) '          are using Rfree)'
      write (*,*)
      write (*,*) 'XBADCON - list residues with bad contacts'
      write (*,*) '          during a refinement'
      write (*,*)
      write (*,*) 'AUTODIC - auto-generate X-PLOR dictionary'
      write (*,*) '          from a PDB file'
      write (*,*)
      write (*,*) 'DIHE    - generate DIHEdral statements to'
      write (*,*) '          restrain a model to another one'
c
c ... "secret" phipsi command
c
c      write (*,*)
c      write (*,*) 'PHIPSI  - generate Phi/Psi restraint file'
c
    1 continue
      write (*,*)
      write (*,*) 'QUIT XSEA X1DT X1DY XFIL LUZZ XBAD AUTO'
      call textin (' What would you like to do ?',which)
      call upcase (which)
      write (*,*)
c
      if (which (1:4) .eq. 'QUIT') then
        goto 900
      else if (which (1:4) .eq. 'XFIL') then
        write (*,*) 'XFILTER'
        write (*,*)
        call xfil
        goto 1
      else if (which (1:4) .eq. 'X1DT') then
        write (*,*) 'X1DTRAZ'
        write (*,*)
        call x1dz
        goto 1
      else if (which (1:4) .eq. 'X1DY') then
        write (*,*) 'X1DYTRA'
        write (*,*)
        call x1dy
        goto 1
      else if (which (1:4) .eq. 'XSEA') then
        write (*,*) 'XSEARCH'
        write (*,*)
        goto 2
      else if (which (1:4) .eq. 'LUZZ') then
        write (*,*) 'LUZZATI'
        write (*,*)
        call luzat
        goto 1
      else if (which (1:4) .eq. 'XBAD') then
        write (*,*) 'XBADCON'
        write (*,*)
        call badcon
        goto 1
      else if (which (1:4) .eq. 'AUTO') then
        write (*,*) 'AUTODIC'
        write (*,*)
        call autodi (lforce,lnoh)
        goto 1
      else if (which (1:4) .eq. 'PHIP') then
        write (*,*) 'PHIPSI'
        write (*,*)
        call phipsi
        goto 1
      else if (which (1:4) .eq. 'DIHE') then
        write (*,*) 'DIHE'
        write (*,*)
        call dihere
        goto 1
      else
        call errcon ('Invalid option : '//which)
        goto 1
      end if
c
    2 continue
      call jvalut (' Assumed angle increment  :',1,inter)
      call jvalut (' Max nr of angles sampled :',1,maxang)
      write (*,*)
c
      ifile = 10
      nfile = 0
      ntot  = 0
      ityp  = 3
      ival  = 0
      fmt = '(8f10.3)'
      do i=1,3
        ilim (1,i) = 0
        ilim (2,i) = 360
      end do
c
      do i1=1,maxang
        do i2=1,maxang
          do i3=1,maxang
            pcval (i1,i2,i3) = xvoid
          end do
        end do
      end do
c
   10 continue
      filenm = ' '
      write (*,*)
      call textin (' Enter name of search.dat file :',filenm)
      if (length(filenm) .lt. 1) goto 100
c
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) goto 100
c
      nval = 0
   20 continue
      read (ifile,*,end=30) i1,i2,i3,pc
      pcval (1+i1/inter,1+i2/inter,1+i3/inter) = pc
      nval = nval + 1
      goto 20
c
   30 continue
      nfile = nfile + 1
      close (ifile)
      call jvalut (' Nr of values read :',1,nval)
      ntot = ntot + nval
      call jvalut (' Total             :',1,ntot)
      goto 10
c
c ... done reading
c
  100 continue
      write (*,*)
      call jvalut (' Nr of files read :',1,nfile)
      call jvalut (' Nr of PC values  :',1,ntot)
c
      if (ntot .lt. 10) call errstp ('Not enough values')
c
      pcmin =  999.9
      pcmax = -999.9
      pctot = 0.0
      pcsdv = 0.0
      ntot  = 0
      do i=1,3
        ijkmax(i) = -1
      end do
c
      do i1=1,maxang
        do i2=1,maxang
          do i3=1,maxang
            if (pcval (i1,i2,i3) .gt. xvoid) then
              ntot = ntot + 1
              pctot = pctot + pcval (i1,i2,i3)
              pcsdv = pcsdv + pcval (i1,i2,i3)**2
              pcmin = min (pcmin,pcval (i1,i2,i3))
              if (pcval (i1,i2,i3) .gt. pcmax) then
                pcmax = pcval (i1,i2,i3)
                ijkmax (1) = i1
                ijkmax (2) = i2
                ijkmax (3) = i3
              end if
            end if
          end do
        end do
      end do
c
      write (*,*)
      pcave = pctot/float(ntot)
      pcsdv = pcsdv/float(ntot) - pcave*pcave
      pcsdv = sqrt (pcsdv)
      call jvalut (' Actual nr of unique data points  :',1,ntot)
      call rvalut (' Average PC value   :',1,pcave)
      call rvalut (' Standard deviation :',1,pcsdv)
      call rvalut (' Minimum PC value   :',1,pcmin)
      call rvalut (' Maximum PC value   :',1,pcmax)
      call ivalut (' Max at point :',3,ijkmax)
c
      nl = 10
      dx = (pcmax-pcave)/float(nl)
      do i=1,nl
        levels (i) = pcave + dx*float(i-1)
      end do
c
  110 continue
c
      call jvalin (' Type of THETA plane (1/2/3; 0 to quit) ?',
     +  1,ityp)
      if (ityp .lt. 1 .or. ityp .gt. 3) goto 1
c
      filenm = 'pcplot.pl2'
      write (*,*)
      call textin (' Name for O2D plot file ?',filenm)
c      if (length(filenm) .lt. 1) goto 1
c
      close (ifile)
      call xopxua (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) goto 100
c
      call jvalin (' Type of THETA plane (1/2/3; 0 to quit) ?',1,ityp)
      if (ityp .lt. 1 .or. ityp .gt. 3) goto 110
c
      call stamp (myline)
      write (ifile,'(2a)') 'REMARK ','XPLO2D - XSEARCH option'
      write (ifile,'(2a)') 'REMARK ',myline(1:leng1(myline))
      write (ifile,'(2a)') 'REMARK ',filenm(1:leng1(filenm))
      write (ifile,'(2a,i10)') 'REMARK ',
     +  'Nr of datapoints ',ntot
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Average PC value ',pcave
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Minimum PC value ',pcmin
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Maximum PC value ',pcmax
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Sigma ',pcsdv
c
      write (ifile,'(a,i3)') 'NLEVEL ',nl
      write (ifile,'(a)') 'LEVELS'
      write (ifile,'(15(f6.3,1x))') (levels(i),i=1,nl)
      write (ifile,'(a)') 'COLOUR'
      write (ifile,'(a)') '1 1 5 5 2 2 6 6 4 4 '
c
      if (ityp .eq. 1) then
        call jvalin (' Which THETA1 value (0-360) ?',1,ival)
        ilim (1,1) = ival
        ilim (2,1) = ival
        call jvalin (' Limits for THETA2 ?',2,ilim(1,2))
        call jvalin (' Limits for THETA3 ?',2,ilim(1,3))
        write (ifile,'(a,i3)') 'XLABEL Theta 2 ; Theta 1 = ',ival
        write (ifile,'(a)') 'YLABEL Theta 3'
      else if (ityp .eq. 2) then
        call jvalin (' Which THETA2 value (0-360) ?',1,ival)
        ilim (1,2) = ival
        ilim (2,2) = ival
        call jvalin (' Limits for THETA1 ?',2,ilim(1,1))
        call jvalin (' Limits for THETA3 ?',2,ilim(1,3))
        write (ifile,'(a,i3)') 'XLABEL Theta 1 ; Theta 2 = ',ival
        write (ifile,'(a)') 'YLABEL Theta 3'
      else if (ityp .eq. 3) then
        call jvalin (' Which THETA3 value (0-360) ?',1,ival)
        ilim (1,3) = ival
        ilim (2,3) = ival
        call jvalin (' Limits for THETA1 ?',2,ilim(1,1))
        call jvalin (' Limits for THETA2 ?',2,ilim(1,2))
        write (ifile,'(a,i3)') 'XLABEL Theta 1 ; Theta 3 = ',ival
        write (ifile,'(a)') 'YLABEL Theta 2'
      end if
c
      do i=1,3
        ilim (4,i) = 1+ilim(1,i)/inter
        ilim (5,i) = 1+ilim(2,i)/inter
        ilim (3,i) = ilim(5,i) - ilim(4,i) + 1
      end do
c
      if (ityp .eq. 1) then
        write (ifile,'(a,i3)') 'XPOINT ',ilim(3,2)
        write (ifile,'(a,i3)') 'YPOINT ',ilim(3,3)
        write (ifile,'(a,i3,1x,i3)') 'XLIMIT ',ilim(1,2),ilim(2,2)
        write (ifile,'(a,i3,1x,i3)') 'YLIMIT ',ilim(1,3),ilim(2,3)
        write (ifile,'(9a)') 'ZVALUE ',fmt
        i1 = ilim(4,1)
        write (ifile,fmt) ((pcval(i1,i2,i3),
     +    i2=ilim(4,2),ilim(5,2)),i3=ilim(4,3),ilim(5,3))
      else if (ityp .eq. 2) then
        write (ifile,'(a,i3)') 'XPOINT ',ilim(3,1)
        write (ifile,'(a,i3)') 'YPOINT ',ilim(3,3)
        write (ifile,'(a,i3,1x,i3)') 'XLIMIT ',ilim(1,1),ilim(2,1)
        write (ifile,'(a,i3,1x,i3)') 'YLIMIT ',ilim(1,3),ilim(2,3)
        write (ifile,'(9a)') 'ZVALUE ',fmt
        i2 = ilim(4,2)
        write (ifile,fmt) ((pcval(i1,i2,i3),
     +    i1=ilim(4,1),ilim(5,1)),i3=ilim(4,3),ilim(5,3))
      else if (ityp .eq. 3) then
        write (ifile,'(a,i3)') 'XPOINT ',ilim(3,1)
        write (ifile,'(a,i3)') 'YPOINT ',ilim(3,2)
        write (ifile,'(a,i3,1x,i3)') 'XLIMIT ',ilim(1,1),ilim(2,1)
        write (ifile,'(a,i3,1x,i3)') 'YLIMIT ',ilim(1,2),ilim(2,2)
        write (ifile,'(9a)') 'ZVALUE ',fmt
        i3 = ilim(4,3)
        write (ifile,fmt) ((pcval(i1,i2,i3),
     +    i1=ilim(4,1),ilim(5,1)),i2=ilim(4,2),ilim(5,2))
      end if
c
      write (ifile,'(a)') '!'
      write (ifile,'(a)') 'END'
c
      close (ifile)
      goto 110
c
c ... done
c
  900 continue
      call gkquit ()
c
      end
c
c
c
      subroutine x1dz
c
      implicit none
c
      integer maxdat
      parameter (maxdat = 16*1024)
c
      real xdata (maxdat),ave,sig,xmi,xma,xlo,xhi,dx,xi,xa
c
      integer ifile,ierr,ndata,length,i,leng1
c
      character filenm*80,filepl*80,line*250,fmt*20,myline*128
c
code ...
c
      write (*,*)
      call jvalut (' Max nr of data points :',1,maxdat)
      write (*,*)
      fmt = '(11f7.4)'
      ifile = 11
c
   10 continue
      filenm = ' '
      write (*,*)
      call textin (' Enter name of xxx.3dmatrix file :',filenm)
      if (length(filenm) .lt. 1) return
c
      close (ifile)
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) return
c
 6000 format (a)
c
      read (ifile,6000) line
      if (line(5:15) .ne. 'Translation') then
        call errcon ('Not a translation function 3dmatrix file')
        return
      end if
c
      read (ifile,6000) line
      read (line(12:18),*) ndata
      call jvalut (' Nr of data points :',1,ndata)
      if (ndata .gt. maxdat) then
        call errcon ('Too many data points in this file')
        return
      end if
c
      read (ifile,6000) line
      read (line(9:15),*) ave
      call fvalut (' Average           :',1,ave)
c
      read (ifile,6000) line
      read (line(11:17),*) sig
      call fvalut (' Sigma             :',1,sig)
c
      read (ifile,6000) line
      read (line(9:15),*) xma
      call fvalut (' Maximum           :',1,xma)
c
      read (ifile,6000) line
      read (line(9:15),*) xmi
      call fvalut (' Minimum           :',1,xmi)
c
      read (ifile,6000) line
      read (ifile,6000) line
      read (line(33:43),*) xlo
      call fvalut (' Search start      :',1,xlo)
c
      read (ifile,6000) line
      read (line(33:43),*) xhi
      call fvalut (' Search end        :',1,xhi)
      if (xlo .ge. xhi) then
        call errcon ('Not a search along the C-axis')
        return
      end if
c
      dx = (xhi - xlo) / (ndata-1)
      call fvalut (' Search step       :',1,dx)
c
      read (ifile,6000) line
      read (ifile,6000) line
      read (ifile,6000) line
c
 7000 format (6(f11.4,1x))
      read (ifile,7000) (xdata(i),i=1,ndata)
c
      close (ifile)
      filepl = ' '
      write (*,*)
      call textin (' Enter name of O2D plot file :',filepl)
      if (length(filepl) .lt. 1) return
c
      call xopxua (ifile,filepl,.true.,ierr)
      if (ierr .ne. 0) return
c
      xi = xmi
      xa = xma
      call rvalin (' Min TF value to display ?',1,xi)
      call rvalin (' Max TF value to display ?',1,xa)
c
      call stamp (myline)
      write (ifile,'(2a)') 'REMARK ','XPLO2D - X1DTRAZ option'
      write (ifile,'(2a)') 'REMARK ',myline(1:leng1(myline))
      write (ifile,'(2a)') 'REMARK ',filepl(1:leng1(filepl))
      write (ifile,'(2a)') 'REMARK ',filenm(1:leng1(filenm))
      write (ifile,'(2a,i10)') 'REMARK ',
     +  'Nr of datapoints ',ndata
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Average PC value ',ave
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Minimum PC value ',xmi
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Maximum PC value ',xma
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Sigma ',sig
      write (ifile,'(a,i10)') 'COLOUR ',4
      write (ifile,'(a,i10)') 'NPOINT ',ndata
      write (ifile,'(a,4f10.5)') 'XYVIEW ',xlo,xhi,xi,xa
      write (ifile,'(a,2f10.5)') 'XLIMIT ',xlo,dx
c
      write (line,'(a,f8.4,a,f8.4,a,f8.4,a,a)')
     +  'XLABEL Transl C-axis ',
     +  xlo,' - ',xhi,', ',dx,'; file ',
     +  filenm(1:leng1(filenm))
      call pretty (line)
      write (ifile,6000) (line(1:leng1(line)))
c
      write (line,'(a,f8.4,a,f8.4,a,f8.4,a,f8.4)')
     +  'YLABEL TF values ',xmi,' - ',xma,
     +  '; ave ',ave,'; sigma ',sig
      call pretty (line)
      write (ifile,6000) (line(1:leng1(line)))
c
      write (ifile,'(9a)') 'YVALUE ',fmt
      write (ifile,fmt) (xdata(i),i=1,ndata)
c
      close (ifile)
c
      goto 10
c
      end
c
c
c
      subroutine xfil
c
      implicit none
c
      integer maxdat
      parameter (maxdat = 256*1024)
c
      real ra(3,maxdat),rf,pc1,pc2,pc(maxdat),period(3)
      real ox,oy,oz,rx,ry,rz,pc3,toler,pcmax,xn
      real xave,xsdv,xmin,xmax,xtot,pcut,rms
c
      integer ifile,ierr,ndata,length,i,ind,code(maxdat)
      integer nmax,j,nn,imax,iper(3),ifmt,mxpk,jfile,kfile,leng1
c
      character filenm*80,line*80
c
code ...
c
      write (*,*)
      call jvalut (' Max nr of data points :',1,maxdat)
      write (*,*)
      ifmt = 1
      mxpk = 50
      ifile = 11
      jfile = 12
      kfile = 13
      pcut = 0.095
      toler = 5.0
      rms = 0.0
      iper(1) = 1
      iper(2) = 1
      iper(3) = 1
c
      write (*,*)
     +  'Periodicity = 360.0 / (searched range of angles)'
      write (*,*)
     +  'Ex: if you searched 0 - 90 => periodicity = 4'
      call ivalin (' Periodicity of rotation angles ?',3,iper)
      do i=1,3
        iper (i) = max (1,iper(i))
        period (i) = 360.0/float(iper(i))
      end do
      call fvalut (' Intervals :',3,period)
c
      filenm = ' '
      call textin (' Enter name of input file :',filenm)
      if (length(filenm) .lt. 1) return
c
      close (ifile)
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) return
c
      write (*,*) 'Select one of the following format types:'
      write (*,*) '(0) = direct rotation search output'
      write (*,*) '(1) = PC-refine output with 1 PC value'
      write (*,*) '(2) = PC-refine output with 2 PC values'
      write (*,*) '(3) = PC-refine output with 3 PC values'
      call ivalin (' Format type ?',1,ifmt)
      if (ifmt .lt. 0 .or. ifmt .gt. 3) return
c
      if (ifmt .eq. 0) then
        call ivalin (' How many peaks to list ?',1,mxpk)
        mxpk = max (10, mxpk)
        filenm = 'direct.pc'
        call textin (' New PC-refinement input file ?',filenm)
        close (jfile)
        call xopxua (jfile,filenm,.true.,ierr)
        if (ierr .ne. 0) return
        write (jfile,'(a)')
     +    '! index, theta1, theta2, theta3, PC-value'
      else
        call fvalin (' PC cut-off for use in Rotman ?',1,pcut)
c
        filenm = 'rotman.awk'
        close (jfile)
        call xopxua (jfile,filenm,.true.,ierr)
        if (ierr .eq. 0) then
          write (jfile,'(a)') '/vs/ {print}'
          write (jfile,'(a)') '/rotman @/ {print}'
          write (jfile,'(a)') '/swap @/ {print}'
          write (jfile,'(a)') '/Eulerian angle/ {print}'
          write (jfile,'(a)') '/spherical polar/ {print}'
          write (jfile,'(a)') '/rotation angle/ {print}'
        end if
        filenm = 'rotman.inp'
        call textin (' New rotman input file ?',filenm)
        close (jfile)
        call xopxua (jfile,filenm,.true.,ierr)
        if (ierr .ne. 0) return
        write (jfile,'(a)')
     +    '{ rotman.inp => relations between the highest }'
        write (jfile,'(a,f6.3,a)')
     +    '{ peaks of the RF after PC refinement (>=',pcut,' }'
        write (jfile,*)
        write (jfile,'(a)')
     +    'xrefin  {insert unit cell and symmops here}'
        write (jfile,'(a)')
     +    'end'
        write (jfile,*)
c
        write (*,*) 'Put cell constants and symmetry operators'
        write (*,*) 'in the Rotman input file yourself.'
        write (*,*) 'Also add a STOP statement at the end.'
        write (*,*) 'Type: awk -f rotman.awk rotman.out > rotman.lis'
        write (*,*) 'after running Rotman to get a summary'
        write (*,*) 'of the results.'
        write (*,*)
        call fvalin (' Tolerance (degrees) ?',1,toler)
        toler = max (0.0,toler)
      end if
c
      ndata = 0
      write (*,*)
c
 6000 format (' Peak ',i4,' Ori = ',3f8.1,' Ref = ',3f10.3/
     +        '       Ind, RF, PC = ',i6,1x,f6.3,1x,f8.4)
 6010 format (' New maximum : ',i6,' PC = ',f8.4,' : ',3f10.3)
 6014 format ('    Nr','   PCmax',' Max/Sig',' (Max-Ave)/Sig',
     +  ' Max/RMS',' Euler angles of this solution')
 6015 format (i6,f8.4,f8.4,f14.4,f8.4,3f10.3)
 6020 format (' ........... : ',i6,' PC = ',f8.4,' : ',3f10.3)
 6030 format (' Number of peaks in cluster : ',i4)
 6040 format (1x,a,3(1x,f10.3))
c
   10 continue
c
      if (ifmt .eq. 0) then
        read (ifile,*,end=100) rx,ry,rz,pc3
      else if (ifmt .eq. 1) then
        read (ifile,*,end=100) ox,oy,oz,rx,ry,rz,ind,rf,pc3
      else if (ifmt .eq. 2) then
        read (ifile,*,end=100) ox,oy,oz,rx,ry,rz,ind,rf,pc1,pc3
      else if (ifmt .eq. 3) then
        read (ifile,*,end=100) ox,oy,oz,rx,ry,rz,ind,rf,pc1,pc2,pc3
      end if
c
      ndata = ndata + 1
      ra (1,ndata) = rx
      ra (2,ndata) = ry
      ra (3,ndata) = rz
      pc (ndata) = pc3
      code (ndata) = 0
      rms = rms + pc3*pc3
c
      do i=1,3
   20   continue
        if (ra(i,ndata) .lt. 0.0) then
          ra(i,ndata) = ra(i,ndata) + period(i)
          goto 20
        else if (ra(i,ndata) .ge. period(i)) then
          ra(i,ndata) = ra(i,ndata) - period(i)
          goto 20
        end if
      end do
c
      if (ifmt.gt.0) write (*,6000)
     +  ndata,ox,oy,oz,rx,ry,rz,ind,rf,pc3
c
      if (ndata .lt. maxdat) goto 10
c
      call errcon ('Max nr of data points; rest skipped')
c
  100 continue
      close (ifile)
      write (*,*)
      call jvalut (' Nr of solutions read :',1,ndata)
c
ccc      if (ifmt .eq. 0) then
        call xstats (pc,ndata,xave,xsdv,xmin,xmax,xtot)
        call fvalut (' Minimum PC-value :',1,xmin)
        call fvalut (' Maximum PC-value :',1,xmax)
        call fvalut (' Average PC-value :',1,xave)
        call fvalut (' Sigma   PC-value :',1,xsdv)
        rms = sqrt (rms/float(ndata))
        call fvalut (' RMS     PC-value :',1,rms)
ccc      end if
c
      nmax = 0
      write (*,*)
c
      if (ifmt .eq. 0) write (*,6014)
c
c ... find maximum & cluster
c
  200 continue
      do i=1,ndata
        if (code(i) .eq. 0) goto 210
      end do
      if (ifmt .eq. 0) return
      goto 900
c
  210 continue
      nmax = nmax + 1
      pcmax = pc(i)
      imax = i
      do j=i+1,ndata
        if (code(j).eq.0 .and. pc(j) .gt. pcmax) then
          pcmax = pc(j)
          imax = j
        end if
      end do
c
      if (ifmt .gt. 0) then
        write (*,*)
        write (*,6010) imax,pcmax,(ra(i,imax),i=1,3)
      else
        write (*,6015) imax,pcmax,(pcmax/xsdv),
     +    ((pcmax-xave)/xsdv),(pcmax/rms),(ra(i,imax),i=1,3)
      end if
c
      code (imax) = nmax
c
      if (ifmt .eq. 0) then
        write (jfile,*) imax,(ra(i,imax),i=1,3),pc(imax)
        if (nmax .ge. mxpk) then
          close (jfile)
          return
        end if
        goto 200
      end if
c
      close (kfile)
      if (pc(imax) .lt. pcut) goto 8362
      write (line,*) 'peak.',nmax
      call remspa (line)
      call xopxua (kfile,line,.false.,ierr)
      if (ierr .eq. 0) then
        write (kfile,'(a,i6,a,f8.5,a)')
     +    '{ Peak nr ',nmax,'; PC = ',pc(imax),' }'
        write (kfile,'(a,3f10.3,a)') '  euler = (',
     +    (ra(i,imax),i=1,3),')'
        if (nmax .gt. 1) then
          do j=1,nmax-1
            write (jfile,'(a,i3,a,i3,a)') '{',j,
     +        ' vs ',nmax,'}'
            write (line,*) 'peak.',nmax
            call remspa (line)
            write (jfile,'(9a)') 'rotman @',
     +        line(1:leng1(line))
            write (line,*) 'peak.',j
            call remspa (line)
            write (jfile,'(9a)') '  swap @',
     +        line(1:leng1(line))
            write (jfile,'(a)') '  dist end'
            write (jfile,*)
          end do
        end if
      end if
c
 8362 continue
      nn = 1
      rx = ra(1,imax)
      ry = ra(2,imax)
      rz = ra(3,imax)
      ox = pc(imax)*rx
      oy = pc(imax)*ry
      oz = pc(imax)*rz
      pc3 = pc(imax)
c
      do i=1,ndata
        if (code(i).eq.0 .and. i.ne.imax) then
          if (abs(ra(1,i)-ra(1,imax)) .le. toler .and.
     +        abs(ra(2,i)-ra(2,imax)) .le. toler .and.
     +        abs(ra(3,i)-ra(3,imax)) .le. toler) then
            code (i) = -imax
            nn = nn + 1
            rx = rx + ra(1,i)
            ry = ry + ra(2,i)
            rz = rz + ra(3,i)
            ox = ox + pc(i)*ra(1,i)
            oy = oy + pc(i)*ra(2,i)
            oz = oz + pc(i)*ra(3,i)
            pc3 = pc3 + pc(i)
            write (*,6020) i,pc(i),(ra(j,i),j=1,3)
          end if
        end if
      end do
c
      write (*,6030) nn
      xn = 1.0/float(nn)
      rx = rx * xn
      ry = ry * xn
      rz = rz * xn
      ox = ox / pc3
      oy = oy / pc3
      oz = oz / pc3
      pc3 = pc3 * xn
c
      write (*,6040) 'Average rotation angles :',rx,ry,rz
      write (*,6040) 'Weighted average angles :',ox,oy,oz
      write (*,6040) 'Average PC value        :',pc3
      write (*,6040) 'Maximum PC value        :',pc(imax)
      write (*,6040) 'Max / RMS               :',pc(imax)/rms
      write (*,6040) 'Max / Sigma             :',pc(imax)/xsdv
      write (*,6040) '(Max - Ave) / Sigma     :',
     +  (pc(imax)-xave)/xsdv
c
      goto 200
c
c ... finished clustering
c
  900 continue
c
      write (*,*)
      call ivalut (' Nr of peak clusters found :',1,nmax)
      call fvalut (' Using tolerance (degrees) :',1,toler)
      write (*,*)
c
      return
      end
c
c
c
      subroutine x1dy
c
      implicit none
c
      integer maxdat
      parameter (maxdat = 16*1024)
c
      real xdata (maxdat),ave,sig,xmi,xma,xlo,xhi,dx,xi,xa
c
      integer ifile,ierr,ndata,length,i,leng1
c
      character filenm*80,filepl*80,line*250,fmt*20,myline*128
c
code ...
c
      write (*,*)
      call jvalut (' Max nr of data points :',1,maxdat)
      write (*,*)
      fmt = '(11f7.4)'
      ifile = 11
c
   10 continue
      filenm = ' '
      write (*,*)
      call textin (' Enter name of xxx.3dmatrix file :',filenm)
      if (length(filenm) .lt. 1) return
c
      close (ifile)
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) return
c
 6000 format (a)
c
      read (ifile,6000) line
      if (line(5:15) .ne. 'Translation') then
        call errcon ('Not a translation function 3dmatrix file')
        return
      end if
c
      read (ifile,6000) line
      read (line(12:18),*) ndata
      call jvalut (' Nr of data points :',1,ndata)
      if (ndata .gt. maxdat) then
        call errcon ('Too many data points in this file')
        return
      end if
c
      read (ifile,6000) line
      read (line(9:15),*) ave
      call fvalut (' Average           :',1,ave)
c
      read (ifile,6000) line
      read (line(11:17),*) sig
      call fvalut (' Sigma             :',1,sig)
c
      read (ifile,6000) line
      read (line(9:15),*) xma
      call fvalut (' Maximum           :',1,xma)
c
      read (ifile,6000) line
      read (line(9:15),*) xmi
      call fvalut (' Minimum           :',1,xmi)
c
      read (ifile,6000) line
      read (ifile,6000) line
      read (line(21:31),*) xlo
      call fvalut (' Search start      :',1,xlo)
c
      read (ifile,6000) line
      read (line(21:31),*) xhi
      call fvalut (' Search end        :',1,xhi)
      if (xlo .ge. xhi) then
        call errcon ('Not a search along the B-axis')
        return
      end if
c
      dx = (xhi - xlo) / (ndata-1)
      call fvalut (' Search step       :',1,dx)
c
      read (ifile,6000) line
c
      do i=1,ndata
        read (ifile,6000) line
        read (ifile,6000) line
        read (ifile,*) xdata(i)
      end do
c
 7000 format (6(f11.4,1x))
c
      close (ifile)
      filepl = ' '
      write (*,*)
      call textin (' Enter name of O2D plot file :',filepl)
      if (length(filepl) .lt. 1) return
c
      call xopxua (ifile,filepl,.true.,ierr)
      if (ierr .ne. 0) return
c
      xi = xmi
      xa = xma
      call rvalin (' Min TF value to display ?',1,xi)
      call rvalin (' Max TF value to display ?',1,xa)
c
      call stamp (myline)
      write (ifile,'(2a)') 'REMARK ','XPLO2D - X1DYTRA option'
      write (ifile,'(2a)') 'REMARK ',myline(1:leng1(myline))
      write (ifile,'(2a)') 'REMARK ',filepl(1:leng1(filepl))
      write (ifile,'(2a)') 'REMARK ',filenm(1:leng1(filenm))
      write (ifile,'(2a,i10)') 'REMARK ',
     +  'Nr of datapoints ',ndata
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Average PC value ',ave
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Minimum PC value ',xmi
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Maximum PC value ',xma
      write (ifile,'(2a,f10.6)') 'REMARK ',
     +  'Sigma ',sig
      write (ifile,'(a,i10)') 'COLOUR ',4
      write (ifile,'(a,i10)') 'NPOINT ',ndata
      write (ifile,'(a,4f10.5)') 'XYVIEW ',xlo,xhi,xi,xa
      write (ifile,'(a,2f10.5)') 'XLIMIT ',xlo,dx
c
      write (line,'(a,f8.4,a,f8.4,a,f8.4,a,a)')
     +  'XLABEL Transl B-axis ',
     +  xlo,' - ',xhi,', ',dx,'; file ',
     +  filenm(1:leng1(filenm))
      call pretty (line)
      write (ifile,6000) (line(1:leng1(line)))
c
      write (line,'(a,f8.4,a,f8.4,a,f8.4,a,f8.4)')
     +  'YLABEL TF values ',xmi,' - ',xma,
     +  '; ave ',ave,'; sigma ',sig
      call pretty (line)
      write (ifile,6000) (line(1:leng1(line)))
c
      write (ifile,'(9a)') 'YVALUE ',fmt
      write (ifile,fmt) (xdata(i),i=1,ndata)
c
      close (ifile)
c
      goto 10
c
      end
c
c
c
      subroutine luzat
c
      implicit none
c
      integer maxdat
      parameter (maxdat = 256)
c
      real xdata(maxdat),ydata(maxdat)
c
      integer ifile,ierr,length,jfile,leng1
c
      logical rfree
c
      character filenm*80,filepl*80,line*250
c
code ...
c
      write (*,*)
      call jvalut (' Max nr of data points :',1,maxdat)
      write (*,*)
      ifile = 11
      jfile = 12
      rfree = .false.
c
      filenm = 'luzzati.list'
      call textin (' X-PLOR Luzzati list file ?',filenm)
      if (length(filenm) .lt. 1) return
c
      close (ifile)
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) return
c
 6000 format (a)
c
      read (ifile,6000) line
c
      rewind (ifile)
c
      if (index(line,'TEST') .gt. 0) then
        rfree = .true.
        write (*,*) 'Dataset partitioned for use of Rfree'
      else
        rfree = .false.
        write (*,*) 'Dataset NOT partitioned for use of Rfree'
        goto 100
      end if
c
c ... if here, do Luzzati plot for free R-factor
c
      filepl = 'luzzati_test.plt'
      call textin (' O2D file for Luzzati plot TEST dataset ?',filepl)
      if (length(filepl) .lt. 1) goto 100
c
      close (jfile)
      call xopxua (jfile,filepl,.true.,ierr)
      if (ierr .ne. 0) return
c
      write (jfile,'(9a)') 'REMARK X-PLOR file ',
     +  filenm(1:leng1(filenm))
      write (jfile,'(9a)') 'REMARK Plot file ',
     +  filepl(1:leng1(filepl))
c
      call luplot (ifile,jfile,'Luzzati plot (TEST data; Rfree)',
     +             xdata,ydata,maxdat)
      close (jfile)
c
c ... find start of R-factor list for WORK set
c
   10 read (ifile,6000,end=20) line
      if (index(line,'WORKING') .gt. 0) goto 100
      goto 10
c
   20 continue
      call errcon ('Data for WORK set not found in input file')
      close (ifile)
      return
c
c ... do Luzzati plot for normal R-factor
c
  100 continue
c
      if (rfree) then
        filepl = 'luzzati_work.plt'
        call textin (' O2D file for Luzzati plot WORK dataset ?',
     +               filepl)
      else
        filepl = 'luzzati.plt'
        call textin (' O2D file for Luzzati plot ?',filepl)
      end if
c
      if (length(filepl) .lt. 1) goto 100
c
      close (jfile)
      call xopxua (jfile,filepl,.true.,ierr)
      if (ierr .ne. 0) return
c
      write (jfile,'(9a)') 'REMARK X-PLOR file ',
     +  filenm(1:leng1(filenm))
      write (jfile,'(9a)') 'REMARK Plot file ',
     +  filepl(1:leng1(filepl))
c
      if (rfree) then
        call luplot (ifile,jfile,'Luzzati plot WORK data (R)',
     +               xdata,ydata,maxdat)
      else
        call luplot (ifile,jfile,'Luzzati plot',xdata,ydata,maxdat)
      end if
c
      close (ifile)
      close (jfile)
c
      return
      end
c
c
c
      subroutine luplot (ifi,jfi,ytext,xdata,ydata,maxdat)
c
      implicit none
c
      integer nref,maxref
      parameter (nref = 25, maxref = 10)
c
      character*1 quote
      parameter (quote = '''')
c
      integer maxdat
c
      real xdata(maxdat),ydata(maxdat),x1,x2,x4,xmin,xmax,ymin,ymax
      real xref(nref),yluz(nref),xluz(nref,maxref),xerr(maxref),xx
c
      integer ifi,length,i,jfi,n3,j,nerr,j1,j2,leng1
c
      character line*250,myline*128,ytext*(*)
c
      data xref /0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
     +           0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,
     +           0.3,0.35,0.4,0.45,0.5/
c
      data yluz /0.0,0.025,0.05,0.074,0.098,0.122,0.145,0.168,0.191,
     +           0.214,0.237,0.281,0.319,0.353,0.385,0.414,0.44,
     +           0.463,0.483,0.502,0.518,0.548,0.564,0.574,0.580/
c
      data nerr /5/
      data xerr /0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.8,1.0/
c
      save nerr
      save xerr
c
code ...
c
      call jvalut (' Max nr of error curves :',1,maxref)
      call jvalin (' Nr of error curves to draw ?',1,nerr)
      if (nerr .lt. 2 .or. nerr .gt. maxref) then
        call errcon ('Invalid nr of error curves; default 5')
        nerr = 5
      end if
c
      call fvalin (' Error levels for curves ?',nerr,xerr)
      do i=1,nerr
        xerr(i) = max(0.001,xerr(i))
      end do
c
      do i=1,nerr
        xx = 1.0/xerr(i)
        do j=1,nref
          xluz (j,i) = xx * xref(j)
        end do
ccc        call fvalut (' LUZ :',nref,xluz(1,i))
      end do
c
 6000 format (a)
 6100 format (a6,1x,10i6)
 6200 format (a6,1x,6f12.5)
 6300 format (a6,1x,a)
c
   10 continue
      read (ifi,6000,err=900,end=900) line
      if (index(line,'resol.-range') .gt. 0) goto 100
      goto 10
c
  100 continue
      i = 0
      xmin = 0.1
      xmax = 0.5
      ymin = 0.0
      ymax = 0.5
c
  110 continue
      read (ifi,6000,err=900,end=120) line
      if (length(line) .lt. 1) goto 120
c
ccc      print *,line(1:leng1(line))
c
      read (line(6:),*,err=900,end=900) x1,x2,n3,x4
      i = i + 1
      if (i .le. maxdat) then
        xdata (i) = 2.0/(x1+x2)
        ydata (i) = x4
        xmin = min (xmin,xdata(i))
        xmax = max (xmax,xdata(i))
        ymin = min (ymin,ydata(i))
        ymax = max (ymax,ydata(i))
      else
        call errcon ('Too many data points; use fewer MBINS')
        i = maxdat
        goto 120
      end if
      goto 110
c
  120 continue
      call jvalut (' Nr of data points :',1,i)
      if (i .lt. 2) then
        call errcon ('Too few data points')
        return
      end if
c
      if (xmin .lt. 0.1) xmin = 0.0
c
      xx = xmax + 0.1
      if (xmax .gt. 0.5) then
        xmax = 0.1 * (int(xmax/0.1) + 1)
      end if
      xmax = xmax + 0.15
c
      if (ymax .gt. 0.5) then
        ymax = 0.1 * (int(ymax/0.1) + 1)
      end if
      ymax = ymax + 0.1
c
      write (jfi,'(2a)') 'REMARK ','XPLO2D - LUZZATI option'
c
      call stamp (myline)
      write (jfi,6000) ('REMARK '//myline(1:leng1(myline)))
c
      write (jfi,'(2a)') 'REMARK ',
     +  'Estimated coordinate errors are shown with the curves'
c
      write (jfi,'(2a,i3)') 'REMARK ',
     +  'Nr of datapoints',i
c
      write (jfi,'(2a)') 'REMARK ',
     +  'Remember: positional error = SQRT(3) * coordinate error'
      write (jfi,'(2a)') 'REMARK ',
     +  'Remember: distance error = SQRT(6) * coordinate error'
c
      write (myline,6100) 'NPOINT',i
      write (jfi,6000) myline(1:leng1(myline))
c
      write (myline,6100) 'COLOUR',4
      write (jfi,6000) myline(1:leng1(myline))
c
      write (myline,6200) 'XYVIEW',xmin,xmax,ymin,ymax
      write (jfi,6000) myline(1:leng1(myline))
c
      myline = '1/Resolution (1/A)'
      write (jfi,6300) 'XLABEL',myline(1:leng1(myline))
c
      myline = ytext
      write (jfi,6300) 'YLABEL',myline(1:leng1(myline))
c
c ... write X data points
c
      write (jfi,6300) 'XVALUE','*'
      write (jfi,'(6f12.6)') (xdata(j),j=1,i)
c
      write (jfi,6300) 'YVALUE','*'
      write (jfi,'(6f12.6)') (ydata(j),j=1,i)
c
c ... write reference curves
c
      do i=1,nerr
c
        write (jfi,6000) '!'
        write (jfi,6000) 'MORE'
c
        j1 = 1
        j2 = nref
        do j=1,nref
          if (xluz(j,i) .lt. xmin) j1=j+1
        end do
        do j=nref,1,-1
          if (xluz(j,i) .gt. xx) j2=j-1
        end do
c
ccc        print *,'I,J1,J2 - ',i,j1,j2
c
        write (myline,6100) 'NPOINT',(j2-j1+1)
        write (jfi,6000) myline(1:leng1(myline))
c
        write (myline,6100) 'COLOUR',1
        write (jfi,6000) myline(1:leng1(myline))
c
        write (jfi,6300) 'XVALUE','*'
        write (jfi,'(6f12.6)') (xluz(j,i),j=j1,j2)
c
        write (jfi,6300) 'YVALUE','*'
        write (jfi,'(6f12.6)') (yluz(j),j=j1,j2)
c
        write (myline,'(a,f12.6,1x,f12.6,1x,i2,1x,a1,f4.2,a1)')
     +    'TEXT   ',(xluz(j2,i)+0.005),(yluz(j2)-0.002),10,
     +    quote,xerr(i),quote
        write (jfi,6000) myline(1:leng1(myline))
c
      end do
c
      write (*,*) 'Plot file written'
c
      return
c
  900 continue
      call errcon ('While reading input file')
c
      return
      end
c
c
c
      subroutine badcon
c
      implicit none
c
      integer maxres
      parameter (maxres = 5000)
c
      integer counts(maxres)
      integer ifile,nres,length,ierr,nline,nbad,ntot,i
c
      character filenm*80,line*250,resnam(0:maxres)*14,chain*4
c
code ...
c
      write (*,*)
      call jvalut (' Max nr of residues :',1,maxres)
      ifile = 11
c
      chain = 'AAAA'
      nres = 0
      do i=1,maxres
        counts(i) = 0
      end do
c
      write (*,*)
      call textin (' Monitor which chain ?',chain)
      call upcase (chain)
c
      write (*,*)
      filenm = ' '
      call textin (' Enter name of X-PLOR output file :',filenm)
      if (length(filenm) .lt. 1) return
c
      close (ifile)
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) return
c
 6000 format (a)
c
c%atoms "AAAA-116 -LYS -HZ3 " and "SOLV-350 -HOH -H2  "(XSYM#  4) only  1.46 A apart
c%atoms "SOLV-283 -HOH -H2  " and "SOLV-284 -HOH -H1  " only  1.49 A apart
c2345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      nline = 0
      nbad = 0
c
   10 continue
      read (ifile,6000,end=9000,err=9900) line
      nline = nline + 1
      if (line(1:7) .ne. ' %atoms') goto 10
      nbad = nbad + 1
      call upcase (line)
c
      if (line(10:13) .eq. chain) then
        ntot = ntot + 1
        resnam (0) = line(10:23)
        if (nres .gt. 0) then
          do i=1,nres
            if (resnam(0) .eq. resnam(i)) then
              counts(i) = counts(i) + 1
              goto 20
            end if
          end do
        end if
c
        nres = nres + 1
        if (nres .gt. maxres) then
          nres = nres - 1
          call errcon (' Too many residues')
          goto 20
        end if
        resnam (nres) = resnam(0)
        counts (nres) = 1
      end if
c
   20 continue
c
      if (line(36:39) .eq. chain) then
        ntot = ntot + 1
        resnam (0) = line(36:49)
        if (nres .gt. 0) then
          do i=1,nres
            if (resnam(0) .eq. resnam(i)) then
              counts(i) = counts(i) + 1
              goto 30
            end if
          end do
        end if
c
        nres = nres + 1
        if (nres .gt. maxres) then
          nres = nres - 1
          call errcon (' Too many residues')
          goto 20
        end if
        resnam (nres) = resnam(0)
        counts (nres) = 1
      end if
c
   30 continue
      goto 10
c
 9000 continue
      close (ifile)
      write (*,*)
      call jvalut (' Nr of lines read from file :',1,nline)
      call jvalut (' Nr of bad contact lines    :',1,nbad)
      call jvalut (' Nr of chain hits           :',1,ntot)
      call jvalut (' Nr of chain residues       :',1,nres)
c
      if (nres .gt. 0) then
        write (*,*)
        write (*,*) 'Residue        Nr of bad contacts'
        do i=1,nres
          write (*,*) resnam(i),counts(i)
        end do
      end if
c
      write (*,*)
      return
c
 9900 continue
      call errcon ('While reading file')
      goto 9000
c
      end
c
c
c
      subroutine autodi (lforce,lnoh)
c
      implicit none
c
      integer maxatm,maxtyp,maxnbr,maxbon,maxang,maximp,maxdih
      parameter (maxatm = 500, maxtyp=500, maxnbr=10)
      parameter (maxbon = maxtyp*maxtyp/2)
      parameter (maxang = maxtyp*maxtyp)
      parameter (maximp = maxatm, maxdih = 2*maxatm)
c
      integer maxelm,maxext,maxfrm
      parameter (maxelm = 150, maxext=100, maxfrm=100)
c
      real twopi,rtodeg,degtor
      parameter (twopi  = 6.2831853071796)
      parameter (rtodeg = 360.0 / twopi)
      parameter (degtor=twopi/360.0)
c
      real extxyz(3,maxatm,maxext)
      real atmxyz(3,maxatm),dismat(maxatm,maxatm),nbrdis(maxnbr,maxatm)
      real sbon(maxbon),sang(maxang),simp(maximp),sdih(maxdih)
      real cimp(maximp),cdih(maxdih),amass(maxatm),xdummy(maxext+1)
      real mibon(maxbon),mabon(maxbon),miang(maxang),maang(maxang)
      real miimp(maximp),maimp(maximp),midih(maxdih),madih(maxdih)
      real batom(maxatm),qatom(maxatm),forcec(4),rtdum(12)
      real atmrad(maxatm),tntwgt(5),onowgt(4),valtor(maxdih)
      real samebo,samean,xang,flat,sixty,wgt,mass
      real dist,angle,tangle,dumcut,sum1,sum2,xdis,qdummy,xave,xsdv
      real xmin,xmax,xtot,rmsd,radius,cutbnd,xx
      real lrbond,lrangl,lrdihe,lrimpr
c
      integer nbr(maxatm),nbrptr(maxnbr,maxatm),deftor(4,maxdih)
      integer nat,length,ierr,i,ifile,j,nbo,ntyp,ki,kj,k,l,ll,m
      integer nbtyp,nbon(maxbon),natyp,nang(maxang)
      integer nityp,nimp(maximp),ndih(maxdih),jfile,numh(maxtyp)
      integer ndtyp,nr,leng1,atomnr(maxatm),nbrtyp(maxnbr,maxatm)
      integer ityp(maxatm),first(maxtyp),i1,i2,i3,nwityp(maxatm)
      integer check1(maxnbr),check2(maxnbr),check3(maxnbr)
      integer elmcnt(maxelm),nextra,ix,nweird,kfile,idum,lfile
      integer nfatal,nform,j1,ncnt,j2,nl
c
      logical conect(maxatm,maxatm),isdihe(maxatm,maxatm)
      logical okbond(maxatm,maxatm),afftor(maxatm,maxdih)
      logical lok(maxatm),lextra(maxatm)
      logical lhydro,loccu,lbfac,lhybr,lthird,lform,leach
      logical ldeffi,ldefpa,ltnt,lono,lforce,lnoh
      logical l0,l60,l180
c
      character filenm*80,line*250,atmnam(maxatm)*4,atmtyp(maxatm)*4
      character tbon(2,maxbon)*4,tang(3,maxang)*4,timp(4,maximp)*4
      character resnam*3,resnum*5,chem(maxatm)*2
      character typnam(maxtyp)*4,tfile*80,pfile*80,mfile*80,answer*1
      character ante*4,segid*4,cfile*80,id*1,ask*1,tdih(4,maxdih)*4
      character fulnam*20,tdum*2,xfile*80,formul*80,myform*80
      character newtyp(maxatm)*4,elmnam(maxelm)*2,nowres*9,extrac*1
      character ubon(2,maxbon)*4,uang(3,maxang)*4,uimp(4,maximp)*4
      character udih(4,maxdih)*4,tntfile*80,tntline*250
      character onofile*80,onoline*250,orform(maxfrm)*80
      character tornam(maxdih)*6
c
      data forcec / 1000.0, 500.0, 750.0, 750.0 /
      data tntwgt / 0.02, 2.0, 2.0, 0.02, 10.0 /
      data onowgt / 0.02, 2.0, 2.0, 20.0 /
c
code ...
c
      ifile = 10
      jfile = 11
      kfile = 12
      lfile = 13
c
      nat = 0
      ntyp = 0
c
c      cutoff = 2.3
c      cutlite = 1.8
c      cdum (1) = cutlite
c      cdum (2) = cutoff
c
      cutbnd = 0.35
c
      samebo = 0.03
      samean = 5.0
c
      flat = 10.0
      sixty = 10.0
c
      lrbond = 0.07
      lrangl = 8.0
      lrdihe = 15.0
      lrimpr = 8.0
c
      id = '_'
      loccu = .false.
      lbfac = .false.
      leach = .false.
      lhybr = .false.
      lthird = .false.
      lform = .false.
      ldeffi = .true.
      ldefpa = .true.
c
      resnam = '???'
      nfatal = 0
      nform = 0
c
      write (*,*)
      call jvalut (' Max nr of atoms       :',1,maxatm)
      call jvalut (' Max nr of atom types  :',1,maxtyp)
      call jvalut (' Max nr of neighbours  :',1,maxnbr)
      call jvalut (' Max nr of bond types  :',1,maxbon)
      call jvalut (' Max nr of angle types :',1,maxang)
      call jvalut (' Max nr of dihedrals   :',1,maxdih)
      call jvalut (' Max nr of impropers   :',1,maximp)
c
      write (*,*)
      write (*,*) 'You can now also get TNT and O 6.x dictionaries.'
      write (*,*) 'This requires that every atom gets its own type.'
      ask = 'N'
      call textin (' Do you want to generate a TNT dictionary ?',ask)
      call upcase (ask)
      ltnt = (ask .eq. 'Y')
      call textut (' Do you want to generate a TNT dictionary :',ask)
c
      call textin (' Do you want to generate an O  dictionary ?',ask)
      call upcase (ask)
      lono = (ask .eq. 'Y')
      call textut (' Do you want to generate an O  dictionary :',ask)
c
      if (ltnt .or. lono) then
        loccu = .false.
        lbfac = .false.
        leach = .true.
        lhybr = .true.
        goto 1400
      end if
c
      write (*,*)
      write (*,*) 'This option has now been improved somewhat in'
      write (*,*) 'that you may help the program figuring out'
      write (*,*) 'which atoms are equivalent etc.'
      write (*,*) 'To this end, the OCCUPANCY column may contain'
      write (*,*) 'the number of hydrogen atoms attached to each'
      write (*,*) 'atom, and the B-FACTOR column may contain a'
      write (*,*) 'number (integer) which is identical for'
      write (*,*) 'identical atom types.'
      write (*,*) 'Only reply Yes to the following questions if'
      write (*,*) 'you have edited the input PDB file appropriately.'
      write (*,*) 'Note that you do not have to use BOTH options.'
c
      ask = 'N'
      call textin (' Do occupancies represent the number of Hs ?',ask)
      call upcase (ask)
      loccu = (ask .eq. 'Y')
      call textut (' Do occupancies represent the number of Hs :',ask)
c
      ask = 'N'
      call textin (' Do B-factors represent atom type flags    ?',ask)
      call upcase (ask)
      lbfac = (ask .eq. 'Y')
      call textut (' Do B-factors represent atom type flags    :',ask)
c
      if (.not. lbfac) then
        write (*,*)
        write (*,*) 'If you like, the program can give every atom'
        write (*,*) 'its own type.  This is useful when you have a'
        write (*,*) 'very high resolution structure of the small'
        write (*,*) 'molecule, or if the program gives warnings'
        write (*,*) 'about large ranges for angles, etc.'
        ask = 'Y'
        call textin (' Give every atom its own type ?',ask)
        call upcase (ask)
        leach = (ask .eq. 'Y')
        call textut (' Give every atom its own type :',ask)
      end if
c
      if (.not. loccu) then
        write (*,*)
        write (*,*) 'The program can also try to figure out the'
        write (*,*) 'nr of hydrogens attached to each C, N or O'
        write (*,*) '(often correctly nowadays !)'
        ask = 'Y'
        call textin (' Figure out nr of Hs automatically ?',ask)
        call upcase (ask)
        lhybr = (ask .eq. 'Y')
        call textut (' Figure out nr of Hs automatically :',ask)
      end if
c
 1400 continue
      write (*,*)
      write (*,*) 'To reduce the input, the program can generate'
      write (*,*) 'the name of all output files itself.'
      ask = 'Y'
      call textin (' Use default names for files ?',ask)
      call upcase (ask)
      ldeffi = (ask .eq. 'Y')
c
      write (*,*)
      write (*,*) 'To reduce the input, the program can use'
      write (*,*) 'sensible defaults for all numerical parameters.'
      ask = 'Y'
      call textin (' Use defaults for numerical parameters ?',ask)
      call upcase (ask)
      ldefpa = (ask .eq. 'Y')
c
      write (*,*)
      write (*,*) 'Supply the name of a PDB file containing your'
      write (*,*) 'ligand or whatever (without hydrogens, please)'
      write (*,*) 'The *first* residue in the file will be processed'
      write (*,*) 'All others will be used as extra example to get'
      write (*,*) 'better bonds length, angle, etc. statistics !'
      filenm = ' '
      call textin (' Name of PDB file ?',filenm)
      if (length(filenm) .lt. 1) return
      call textut (' Name of PDB file :',filenm)
c
      close (ifile)
      call xopxoa (ifile,filenm,.true.,ierr)
      if (ierr .ne. 0) return
c
      cfile = filenm(1:leng1(filenm))//'_clean'
      write (*,*)
      if (.not. ldeffi) then
        write (*,*) 'You will get a new, clean PDB file which can'
        write (*,*) 'be read by X-PLOR without parsing problems'
        call textin (' Name of new, clean PDB file ?',cfile)
        if (length(cfile) .lt. 1) return
      end if
      call textut (' Name of new, clean PDB file :',cfile)
c
      close (jfile)
      call xopxua (jfile,cfile,.true.,ierr)
      if (ierr .ne. 0) return
c
      xfile = filenm(1:leng1(filenm))//'_xplo2d'
      write (*,*)
      if (.not. ldeffi) then
        write (*,*) 'You will also get a new PDB file which has'
        write (*,*) 'XPLO2Ds atom type as the B-factor and the'
        write (*,*) 'number of hydrogens as the occupancy (if'
        write (*,*) 'you asked the program to deduce these).'
        write (*,*) 'If the topology and parameter files contain'
        write (*,*) 'atom-type conflicts (or errors in the number'
        write (*,*) 'of Hs), just edit this PDB file and run it'
        write (*,*) 'through XPLO2D again !'
        call textin (' Name of this PDB file ?',xfile)
        if (length(xfile) .lt. 1) return
      end if
      call textut (' Name of this PDB file :',xfile)
c
c --- GET DOWN TO BUSINESS --------------------------------------------------
c
 4713 format (12x,a4,14x,3f8.3,2f6.2)
 4715 format (' Atom # ',i4,' = ',a4,' @ ',3f10.3,2f6.2)
 4720 format (' Bond ',a4,' - ',a4,' = ',f8.3,' A')
c
c ... read PDB file
c
cATOM      1  C2  LEO     1      18.532   0.254  57.030  0.00 20.00      ALEO
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      write (*,*)
      call prompt (' Reading & writing PDB file ...')
      nowres = '!@#$%^&*('
      nextra = 0
c
c ... process PDB file
c
   10 continue
      read (ifile,'(a)',end=19) line
      if (line(1:6) .eq. 'ATOM  ' .or.
     +    line(1:6) .eq. 'HETATM') then
c
        if (nat .eq. 0) nowres = line(18:26)
c
        if (nat .gt. 0) then
          if (nowres .ne. line(18:26)) goto 110
        else
          resnam = line (18:20)
          resnum = line (22:26)
          segid  = line (73:76)
          write (*,*)
          call textut (' Using residue :',resnam)
          call textut (' Identifier    :',resnum)
c
          if (nform .gt. 0) then
            lform = .false.
            do i=1,nform
              if (.not. lform) then
                if (orform(i)(1:3) .eq. resnam) then
                  formul = orform(i) (5:)
                  lform = (length(formul) .gt. 1)
                end if
              end if
            end do
          end if
c
          if (lform) then
            call textut (' Using formula :',formul)
          else
            call prompt (' No chemical formula found')
          end if
          write (*,*)
c
        end if
c
        if (nat .ge. maxatm) then
          call errcon ('Too many atoms; rest skipped')
          goto 19
        end if
        nat = nat + 1
        read (line,4713) atmnam(nat),(atmxyz(j,nat),j=1,3),
     +    qatom(nat),batom(nat)
        call remspa (atmnam(nat)(2:4))
        write (*,4715) nat,atmnam(nat),(atmxyz(j,nat),j=1,3),
     +    qatom(nat),batom(nat)
c
c ... skip atom if duplicate name
c
        if (nat .gt. 1) then
          do j=1,nat-1
            if (atmnam(nat) .eq. atmnam(j)) then
              call errcon ('Duplicate atom name - skipped')
              call textut (' Duplicate of :',atmnam(nat))
              nat = nat - 1
              goto 10
            end if
          end do
        end if
c
c ... try to figure out chemical element type of this atom
c
        call elinfo (atmnam(nat)(1:2),fulnam,nr,mass,radius,.false.)
        if (radius .le. 0.0) radius = 1.0
c
c ... if single character symbol and found, accept it
c
        if (nr .gt. 0 .and. atmnam(nat)(1:1).eq. ' ') then
          if (lnoh .and. nr .eq. 1) then
            call prompt (' Skipping hydrogen !')
            nat = nat - 1
            goto 10
          end if
          atomnr (nat) = nr
          amass (nat) = mass
          chem (nat) = atmnam(nat)(1:2)
          atmrad (nat) = radius
          goto 11
        end if
c
c ... if two characters, check formula (if read)
c
        if (nr .gt. 0) then
          if (lform) then
            i1 = index (formul,atmnam(nat)(1:2))
            if (i1 .gt. 1) then
              atomnr (nat) = nr
              amass (nat) = mass
              chem (nat) = atmnam(nat)(1:2)
              atmrad (nat) = radius
              goto 11
            end if
c
c ... if no formula, we'll have to accept it
c
          else
            atomnr (nat) = nr
            amass (nat) = mass
            chem (nat) = atmnam(nat)(1:2)
            atmrad (nat) = radius
            goto 11
          end if
        end if
c
c ... if here, still not found
c
c ... try hydrogen
c
        if (atmnam(nat)(1:1) .eq. 'H' .or.
     +      atmnam(nat)(2:2) .eq. 'H' .or.
     +      lhydro(atmnam(nat)) ) then
          if (lnoh) then
            call prompt (' Skipping hydrogen !')
            nat = nat - 1
            goto 10
          end if
          atomnr (nat) = 1
          call errcon ('Unknown chemical - assuming hydrogen !!!')
          amass (nat) = 1.008
          chem (nat) = ' H'
          atmrad (nat) = 0.23
          goto 11
        end if
c
c ... try the second character by itself
c
        tdum = ' '//atmnam(nat)(2:2)
        call elinfo (tdum,fulnam,nr,mass,radius,.false.)
        if (nr .gt. 0) then
          atomnr (nat) = nr
          call textut (' Unknown chemical - assuming :',fulnam)
          amass (nat) = mass
          chem (nat) = tdum
          atmrad (nat) = radius
          goto 11
        end if
c
c ... try the first character by itself
c
        tdum = ' '//atmnam(nat)(1:1)
        call elinfo (tdum,fulnam,nr,mass,radius,.false.)
        if (nr .gt. 0) then
          atomnr (nat) = nr
          call textut (' Unknown chemical - assuming :',fulnam)
          amass (nat) = mass
          chem (nat) = tdum
          atmrad (nat) = radius
          goto 11
        end if
c
c ... if nothing worked, assume carbon
c
        atomnr (nat) = 6
        call errcon ('Unknown chemical - assuming carbon !!!')
        amass (nat) = 12.011
        chem (nat) = ' C'
        atmrad (nat) = 0.68
c
   11   continue
c
        qatom (nat) = float(nint(qatom(nat)))
        batom (nat) = float(nint(batom(nat)))
c
c ... 040901 - update (possible de-spaced) atom name
c
        line (13:16) = atmnam(nat)
        line (73:76) = segid
        line (1:6) = 'ATOM  '
        line (55:66) = '  1.00 20.00'
        write (jfile,'(a)') line(1:leng1(line))
c
      else
c
        call textut (' >',line)
        call upcase (line)
c
cREMARK ADR Formula C15 H23 N5 O14 P2
c12345678901234567890123456789012345678901234567890
c
        if (line(1:6)   .eq. 'REMARK' .and.
     +      line(12:18) .eq. 'FORMULA') then
          nform = min (nform + 1, maxfrm)
          orform (nform) = ' '
          orform (nform)(1:3) = line (8:10)
          orform (nform)(5:)  = line (19:)
          call textut (' Formula :',orform(nform))
c
c            formul = line (19:)
c            lform = (length(formul) .gt. 1)
c            if (lform)
c     +        call textut (' ==> List of elements :',formul)
c
        end if
c
cFORMUL   2  REA    C20 H28 O2
cFORMUL   2  ACO  4(C23 H38 N7 O17 P3 S1)                            
c12345678901234567890123456789012345678901234567890
c
        if (line(1:6) .eq. 'FORMUL') then
          nform = min (nform + 1, maxfrm)
          orform (nform) = ' '
          orform (nform)(1:3) = line (13:15)
          orform (nform)(5:)  = line (19:)
          call textut (' Formula :',orform(nform))
c
c          if (resnam .ne. '???') then
c            if (line(13:15) .eq. resnam) then
c              formul = line (18:)
c              lform = (length(formul) .gt. 1)
c            end if
c          else
c            formul = line (18:)
c            lform = (length(formul) .gt. 1)
c          end if
c            if (lform)
c     +        call textut (' ==> List of elements :',formul)
c
        end if
c
      end if
c
      goto 10
c
c ... handle extra examples of the compound
c
  110 continue
      if (line(18:20) .ne. resnam) then
        read (ifile,'(a)',end=219) line
        goto 110
      end if
      if (nextra .gt. 0) then
        do i=1,nat
          if (.not. lextra(i)) then
            nextra = nextra - 1
            call textut (' Missing :',atmnam(i))
            goto 112
          end if
        end do
        call jvalut (' Example OK #',1,nextra)
      end if
c
  112 continue
      nowres = line(18:26)
      call textut (' Extra example ?',nowres)
      if ( (nextra+1) .gt. maxext) then
        call errcon ('Too many extra examples')
        call jvalut (' Maximum allowed :',1,maxext)
        goto 19
      end if
      nextra = nextra + 1
      do i=1,nat
        lextra (i) = .false.
      end do
      backspace (ifile)
c
c ... process this extra copy of the residue
c
  120 continue
      read (ifile,'(a)',end=219) line
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') goto 120
c
ccc      call textut (' >',line)
c
      if (line(18:26) .ne. nowres) goto 110
c
      do i=1,nat
ccc        print *,'|',atmnam(i),'|',line(13:16),'|'
        if (atmnam(i) .eq. line(13:16)) then
          lextra (i) = .true.
          read (line(31:54),'(3f8.3)') (extxyz(j,i,nextra),j=1,3)
ccc          print *,' OK ',atmnam(i)
          goto 120
        end if
      end do
      call textut (' Atom not found :',line(13:16))
c
      goto 120
c
  219 continue
      if (nextra .gt. 0) then
        do i=1,nat
          if (.not. lextra(i)) then
            nextra = nextra - 1
            goto 19
          end if
        end do
        call jvalut (' Example OK #',1,nextra)
      end if
c
   19 continue
      nextra = max (0, min ( maxext, nextra ) )
      write (jfile,'(a)') 'END'
      close (ifile)
      close (jfile)
      write (*,*)
c
      call jvalut (' Nr of atoms read     :',1,nat)
      call jvalut (' Nr of extra examples :',1,nextra)
c
      write (*,*)
      do i=1,nat
        write (*,6023) i,atmnam(i),atomnr(i),chem(i),
     +    amass(i),atmrad(i)
      end do
c
 6023 format (' # ',i4,' |',a4,'| Elem# ',i3,' |',a2,'| Mass = ',
     +  f8.3,' ... Cov. bond radius = ',f6.2)
c
ccccc      if (nat .lt. 2) return
c
c ... calc RMSDs for extra examples
c
      if (nextra .gt. 0) then
        write (*,*)
        do ix=1,nextra
          call lsqgjk (atmxyz(1,1),extxyz(1,1,ix),nat,rmsd,rtdum,ierr)
          if (ierr .ne. 0) then
            write (*,6860) ix
          else
            write (*,6865) ix,rmsd
          end if
        end do
      end if
c
 6860 format (' ERROR - cannot calculate RMSD to example # ',i3)
 6865 format (' RMSD (all atoms) for example # ',i3,' = ',f8.3,' A')
c
c ... print formula
c
      do i=1,maxelm
        elmcnt (i) = 0
        elmnam (i) = '??'
      end do
c
      do i=1,nat
        call elinfo (chem(i),fulnam,nr,mass,radius,.false.)
        if (nr .gt. 0 .and. nr .le. maxelm) then
          elmcnt (nr) = elmcnt (nr) + 1
          elmnam (nr) = chem (i)
        end if
      end do
c
      myform = ' '
      i1 = 1
      do i=1,maxelm
        if (elmcnt(i) .gt. 0) then
          write (myform(i1:),'(1x,a2,i3)') elmnam(i),elmcnt(i)
          i1 = i1 + 6
        end if
      end do
      call pretty (myform)
c
      write (*,*)
      if (lform) call textut (' ==> List of elements :',formul)
      call textut (' ==> Deduced formula  :',myform)
c
      write (*,*)
      write (*,*) 'Provide a single character that will be used to'
      write (*,*) 'generate atom type names.  For instance, if you'
      write (*,*) 'supply an "X", your carbon types will be CX1,'
      write (*,*) 'CX2, etc.'
      call textin (' Character ?',id)
      if (id .eq. ' ') id = 'X'
      call upcase (id)
      call textut (' Character :',id)
c
c      write (*,*)
c      write (*,*) 'I need to know which atoms are chemically bonded.'
c      write (*,*) 'You must supply two distance cut-offs for bonded'
c      write (*,*) 'atoms.  The first is used if both atoms are light'
c      write (*,*) '(H,He,Li,Be,B,C,N,O,F,Ne), the second if either'
c      write (*,*) 'or both are heavy (Na,Mg,Al,Si,P,S,Cl,Ar,...).'
c      call fvalin (' "Light" and "heavy" cut-offs ?',2,cdum)
c      call fvalut (' "Light" and "heavy" cut-offs :',2,cdum)
c      cutlite = cdum(1)
c      cutoff  = cdum(2)
c
      write (*,*)
      if (.not. ldefpa) then
        write (*,*) 'Atoms are considered bonded if their distance'
        write (*,*) 'is not greater than the sum of their covalent'
        write (*,*) 'bond radii plus a small tolerance.'
        write (*,*) 'What value of this tolerance shall I use ?'
        call fvalin (' Bond tolerance ?',1,cutbnd)
      end if
      call fvalut (' Bond tolerance :',1,cutbnd)
c
      write (*,*)
      if (.not. ldefpa) then
        write (*,*) 'In order to assign atom types, I have to compare'
        write (*,*) 'neighbour atom types, bond lengths and (in the'
        write (*,*) 'second sphere of neighbours) bond angles'
        write (*,*) 'Supply tolerance levels for me to decide if'
        write (*,*) 'two bond lengths or angles are "identical"'
        call fvalin (' Tolerance for "identical" bond lengths ?',
     +    1,samebo)
        call fvalin (' Tolerance for "identical" bond angles  ?',
     +    1,samean)
      end if
      call fvalut (' Tolerance for "identical" bond lengths :',
     +  1,samebo)
      call fvalut (' Tolerance for "identical" bond angles  :',
     +  1,samean)
c
      write (*,*)
      if (.not. ldefpa) then
        write (*,*) 'To help me decide if a DIHEdral is special,'
        write (*,*) 'i.e., close to a multiple of +/- 60, +/- 90,'
        write (*,*) '+/- 120, +/- 180 or 0 degrees, I need a'
        write (*,*) 'tolerance level. This may be greater for'
        write (*,*) '60/90/120 than for 0/180 dihedrals (the latter'
        write (*,*) 'enforce flat rings and double bonds).'
        call fvalin (
     +    ' Tolerance for +-60/90/120 DIHEdrals ?',1,sixty)
        call fvalin (
     +    ' Tolerance for 0/+-180/360 DIHEdrals ?',1,flat)
      end if
      call fvalut (
     +  ' Tolerance for +-60/90/120 DIHEdrals :',1,sixty)
      call fvalut (
     +  ' Tolerance for 0/+-180/360 DIHEdrals :',1,flat)
c
      write (*,*)
      if (.not. ldefpa) then
        write (*,*) 'The program will warn you if the observed'
        write (*,*) 'range of values for boond lengths etc. is'
        write (*,*) 'large.'
        call fvalin (' Large range for bond lengths (A) ?',1,lrbond)
        call fvalin (' Large range for bond angles (d)  ?',1,lrangl)
        call fvalin (' Large range for dihedrals (d)    ?',1,lrdihe)
        call fvalin (' Large range for impropers (d)    ?',1,lrimpr)
      end if
      call fvalut (' Large range for bond lengths (A) :',1,lrbond)
      call fvalut (' Large range for bond angles (d)  :',1,lrangl)
      call fvalut (' Large range for dihedrals (d)    :',1,lrdihe)
      call fvalut (' Large range for impropers (d)    :',1,lrimpr)
c
      write (*,*)
      if (.not. ldefpa) then
        write (*,*) 'Supply values for the force constants for'
        write (*,*) 'bond lengths, bond angles, fixed dihedrals'
        write (*,*) 'and impropers, respectively.  The default'
        write (*,*) 'values are in the Engh & Huber ball-park.'
        call fvalin (' Force constants ?',4,forcec)
      end if
      call fvalut (' Force constants :',4,forcec)
c
      if (ltnt) then
        write (*,*)
        if (.not. ldefpa) then
          write (*,*) 'Supply values for the TNT standard deviations'
          write (*,*) 'for bond lengths, bond angles, torsions,'
          write (*,*) 'planes/trigonal centres, and BCORREL (resp.).'
          call fvalin (' TNT standard deviations ?',5,tntwgt)
        end if
        call fvalut (' TNT standard deviations :',5,tntwgt)
      end if
c
      if (lono) then
        write (*,*)
        if (.not. ldefpa) then
          write (*,*) 'Supply values for the O standard deviations'
          write (*,*) 'for bond lengths, bond angles, fixed'
          write (*,*) 'torsions/impropers, and free torsions.'
          call fvalin (' O standard deviations ?',4,onowgt)
        end if
        call fvalut (' O standard deviations :',4,onowgt)
      end if
c
      write (*,*)
      tfile = resnam//'.top'
      call remspa (tfile)
      call locase (tfile)
      if (.not. ldeffi) then
        write (*,*) 'The first file you get is a TOPOLOGY file'
        write (*,*) 'You must check the MASSes, add CHARges,'
        write (*,*) 'UNcomment any DIHEdrals you want to impose'
        write (*,*) 'and check the IMPRopers, DONOrs and ACCEptors'
        call textin (' Name of X-PLOR topology file ?',tfile)
        if (length(tfile) .lt. 1) return
      end if
      call textut (' Name of X-PLOR topology file :',tfile)
c
      write (*,*)
      pfile = resnam//'.par'
      call remspa (pfile)
      call locase (pfile)
      if (.not. ldeffi) then
        write (*,*) 'The second file you get is a PARAMETER file'
        write (*,*) 'You must check the NONBonded parameters and add'
        write (*,*) 'any DIHEdrals you want to use'
        write (*,*) 'Weights are in the same ball park as those used'
        write (*,*) 'by Engh & Huber; you may want to change them'
        call textin (' Name of X-PLOR parameter file ?',pfile)
        if (length(pfile) .lt. 1) return
      end if
      call textut (' Name of X-PLOR parameter file :',pfile)
c
      write (*,*)
      mfile = resnam//'_min.inp'
      call remspa (mfile)
      call locase (mfile)
      if (.not. ldeffi) then
        write (*,*) 'The third X-PLOR file you get is an INPUT file'
        write (*,*) 'Run this file with X-PLOR to perform pure'
        write (*,*) 'geometric (and van der Waals) energy minimisation'
        write (*,*) 'Then check if the resulting structure looks like'
        write (*,*) 'your input structure (especially if the input was'
        write (*,*) 'a high-resolution small molecule X-ray structure!)'
        write (*,*) 'If they differ, patch up some of the parameters,'
        write (*,*) 'impose more IMPRoper or DIHEdral restraints, etc.'
        write (*,*) 'Remember that the result of the minimisation'
        write (*,*) 'shows the structure that X-PLOR will try to get at'
        write (*,*) 'once you start including experimental data !!!'
        write (*,*)
        write (*,*) '     *** GARBAGE IN ==> GARBAGE OUT ***'
        write (*,*)
c
        call textin (' Name of X-PLOR minimisation file ?',mfile)
        if (length(mfile) .lt. 1) return
      end if
      call textut (' Name of X-PLOR minimisation file :',mfile)
c
      if (ltnt) then
        write (*,*)
        tntfile = resnam//'.tnt'
        call remspa (tntfile)
        call locase (tntfile)
        if (.not. ldeffi) then
ccc          write (*,*) 'You will get a TNT dictionary file for free'
          write (*,*)
c
          call textin (' Name of TNT dictionary file ?',tntfile)
          if (length(tntfile) .lt. 1) return
        end if
        call textut (' Name of TNT dictionary file :',tntfile)
      end if
c
      if (lono) then
        write (*,*)
        onofile = resnam//'.odb'
        call remspa (onofile)
        call locase (onofile)
        if (.not. ldeffi) then
ccc          write (*,*) 'You will get an O dictionary file for free'
          write (*,*)
c
          call textin (' Name of O dictionary file ?',onofile)
          if (length(onofile) .lt. 1) return
        end if
        call textut (' Name of O dictionary file :',onofile)
      end if
c
      write (*,*)
      call prompt (' Looking for bonded atoms ...')
c
 2239 continue
c
      nbo = 0
      if (nat .lt. 2) goto 6500
c
      do i=1,nat-1
        do j=i+1,nat
          if (atmnam(i)(1:2) .ne. ' H' .or.
     +        atmnam(j)(1:2) .ne. ' H') then
c
            dismat(i,j) = dist(i,j,atmxyz)
c
            if (nextra .gt. 0) then
              xdis = dismat(i,j)
              xdummy (1) = xdis
              do ix=1,nextra
                qdummy = dist (i,j,extxyz(1,1,ix))
                xdis = xdis + qdummy
                xdummy (ix+1) = qdummy
              end do
              dismat(i,j) = xdis / float (nextra + 1)
            end if
c
            if (dismat(i,j) .lt. 0.8) then
              write (*,*)
              call errcon ('Distance shorter than 0.8 A !')
              call textut (' Between :',atmnam(i))
              call textut ('     and :',atmnam(j))
              nfatal = nfatal + 1
            end if
c
            dismat(j,i) = dismat(i,j)
            conect (i,j) = .false.
            conect (j,i) = .false.
            isdihe (i,j) = .false.
            isdihe (j,i) = .false.
c
            if (lhydro(atmnam(i)) .and.
     +          lhydro(atmnam(j))) goto 2205
c
c            if (atomnr(i).le.10 .and. atomnr(j).le.10) then
c              dumcut = cutlite
c            else
c              dumcut = cutoff
c            end if
c
            dumcut = atmrad(i) + atmrad(j) + cutbnd
c
            if (dismat(i,j) .le. dumcut) then
              nbo = nbo + 1
              conect (i,j) = .true.
              conect (j,i) = .true.
              if (nextra .le. 0) then
                write (*,4720) atmnam(i),atmnam(j),dismat(i,j)
              else
                call xstats (xdummy,nextra+1,xave,xsdv,xmin,xmax,xtot)
                write (*,4723) atmnam(i),atmnam(j),dismat(i,j),
     +            xsdv,xmin,xmax
                if ( (xmax-xmin) .gt. 0.05 ) then
                  call fvalut (
     +    ' WARNING - Large variation in previous bond length :',1,
     +    (xmax-xmin))
                end if
              end if
            end if
 2205       continue
          end if
        end do
        dismat(i,i) = 0.0
        conect(i,i) = .false.
      end do
c
      if (nfatal .gt. 0 .and. (.not. lforce)) then
        call errcon ('Invalid bonded distances - aborting !')
        return
      end if
      nfatal = 0
c
 4723 format (' Bond ',a4,' - ',a4,' = ',f8.3,
     +        ' A; ESD,Min,Max : ',3f8.3)
c
 6500 continue
      write (*,*)
      call jvalut (' Nr of bonds found :',1,nbo)
c
c ... enough bonds found ?
c
      if (nbo .lt. (nat-1)) then
        if (cutbnd .gt. 0.5) then
          call jvalut (' Nr of atoms :',1,nat)
          call jvalut (' Expected min nr of bonds :',1,(nat-1))
          call jvalut (' Actual nr of bonds found :',1,nbo)
          call errstp ('Too few bonds found !')
        end if
        call errcon ('Too few bonds found')
        cutbnd = cutbnd + 0.1
        call fvalut (' Increase CUTBND to (A):',1,cutbnd)
        goto 2239
      end if
c
ccccc      if (nbo .lt. 1) return
c
      write (*,*)
      call prompt (' Generating neighbour lists ...')
c
      do i=1,nat
        nbr(i) = 0
        do j=1,nat
          if (conect(i,j)) then
            if (nbr(i) .ge. maxnbr) then
              call errcon (' Atom has too many neighbours')
              call textut (' Atom :',atmnam(i))
            else
              nbr(i) = nbr(i) + 1
              nbrptr (nbr(i),i) = j
              nbrtyp (nbr(i),i) = atomnr(j)
              nbrdis (nbr(i),i) = dismat(i,j)
            end if
          end if
        end do
        if (nbr(i) .le. 0) then
          write (*,4725) atmnam(i)
        else
          write (*,4730) atmnam(i),(nbrtyp(j,i),nbrdis(j,i),j=1,nbr(i))
        end if
      end do
c
      if (lhybr) then
        call counth (nat,maxatm,qatom,dismat,atmxyz,atomnr,
     +               atmnam,nbr,nbrptr,maxnbr)
        loccu = .true.
ccc        call fvalut (' Q ',nat,qatom)
      end if
c
      do i=1,nat
        atmtyp (i) = ' '
        ityp (i) = -1
      end do
c
      if (lbfac) then
c
        write (*,*)
        call prompt (' Gathering atom types by B-factor ...')
        do i=1,nat
          if (atmtyp(i) .eq. ' ') then
            ntyp = ntyp + 1
            write (line,'(a2,a1,i4)') chem(i),id,ntyp
            call remspa (line)
            ityp (i) = ntyp
            first (ntyp) = i
            atmtyp (i) = line(1:4)
            typnam (ntyp) = line(1:4)
            write (*,4735) atmnam(i),atmtyp(i)
            numh (ntyp) = 0
            if (loccu) numh (ntyp) = nint(qatom(i))
c
            if (i .ge. nat) goto 129
            do j=i+1,nat
              if (batom(j) .eq. batom(i)) then
                atmtyp (j) = atmtyp (i)
                ityp (j) = ityp (i)
                write (*,4735) atmnam(j),atmtyp(j)
              end if
            end do
          end if
        end do
        goto 129
c
      else if (leach) then
c
        write (*,*)
        call prompt (' Assigning every atom its own type ...')
        do i=1,nat
          ntyp = ntyp + 1
          write (line,'(a2,a1,i4)') chem(i),id,ntyp
          call remspa (line)
          ityp (i) = ntyp
          first (ntyp) = i
          atmtyp (i) = line(1:4)
          typnam (ntyp) = line(1:4)
          write (*,4735) atmnam(i),atmtyp(i)
          numh (ntyp) = 0
          if (loccu) numh (ntyp) = nint(qatom(i))
        end do
        goto 129
c
      else
c
        write (*,*)
        call prompt (' First round of atom typing ...')
c
        do i=1,nat
          if (atmtyp(i) .eq. ' ') then
            ntyp = ntyp + 1
            write (line,'(a2,a1,i4)') chem(i),id,ntyp
            call remspa (line)
            ityp (i) = ntyp
            first (ntyp) = i
            atmtyp (i) = line(1:4)
            typnam (ntyp) = line(1:4)
            write (*,4735) atmnam(i),atmtyp(i)
            numh (ntyp) = 0
            if (loccu) numh (ntyp) = nint(qatom(i))
            if (i .ge. nat) goto 125
c
            do j=i+1,nat
c
              if (atmtyp(j) .ne. ' ') goto 123
              if (atomnr(j) .ne. atomnr(i)) goto 123
c
ccc              if (atmnam(j)(1:2) .ne. atmnam(i)(1:2)) goto 123
c
              if (nbr(i) .ne. nbr(j)) goto 123
c
c ... same number of hydrogens ?
c
              if (loccu) then
                if (abs(qatom(i)-qatom(j)) .gt. 0.01) goto 123
              end if
c
c ... same type of neighbours at similar distances ?
c
c ... (1) sum of bonded distances must be similar
c
              sum1 = 0.0
              sum2 = 0.0
              do ki = 1,nbr(i)
                sum1 = sum1 + nbrdis(ki,i)
                sum2 = sum2 + nbrdis(ki,j)
              end do
              if (abs(sum1-sum2) .gt. float(nbr(i))*samebo) then
                goto 123
              end if
c
c ... (2) each nbr must have a counterpart for the other atom
c         (in fact, this is not strict enough: there should be
c         an exact one-to-one mapping !)
c
              do ki=1,nbr(i)
                do kj=1,nbr(j)
                  if (nbrtyp(ki,i) .eq. nbrtyp(kj,j)) then
                    if (abs(nbrdis(ki,i)-nbrdis(kj,j)) .le. samebo)
     +                goto 113
                  end if
                end do
                goto 123
  113           continue
              end do
c
              atmtyp (j) = atmtyp (i)
              ityp (j) = ityp (i)
              write (*,4735) atmnam(j),atmtyp(j)
c
  123       continue
            end do
          end if
  125     continue
        end do
      end if
c
  129 continue
      write (*,*)
      call jvalut (' Nr of atoms :',1,nat)
      call jvalut (' Nr of types :',1,ntyp)
c
      if (lbfac .or. leach) goto 1234
c
c ... check again
c
      write (*,*)
      call prompt (' Second round of atom typing ...')
c
 888  continue
      i1 = 0
      do i=1,nat
        lok(i) = .false.
        newtyp (i) = atmtyp (i)
        nwityp (i) = ityp (i)
      end do
c
 801  continue
      do i=1,nat-1
c
        do k=1,nbr(i)
          check1 (k) = ityp (nbrptr(k,i))
        end do
ccc        call isort (nbr(i),check1)
        if (nbr(i) .gt. 0) call hsorti (nbr(i),check1)
c
        do j=i+1,nat
          if (lok(j)) goto 800
          if (atmtyp(i) .ne. atmtyp(j)) goto 800
          if (nbr(i) .ne. nbr(j)) goto 899
c
          do k=1,nbr(j)
            check2 (k) = ityp (nbrptr(k,j))
          end do
ccc          call isort (nbr(j),check2)
          if (nbr(i) .gt. 0) call hsorti (nbr(j),check2)
c
          do k=1,nbr(i)
            if (check1(k) .ne. check2(k)) goto 899
          end do
c
c ... more checks ?
c
          lok (j) = .true.
          goto 800
c
c ... if here, make atom J a new type
c
  899     continue
          ntyp = ntyp + 1
          write (line,'(a2,a1,i4)') chem(j),id,ntyp
          call remspa (line)
          nwityp (j) = ntyp
          first (ntyp) = j
          newtyp (j) = line(1:4)
          typnam (ntyp) = line(1:4)
          write (*,4735) atmnam(j),newtyp(j)
          numh (ntyp) = 0
          if (loccu) numh (ntyp) = nint(qatom(j))
          i1 = i1 + 1
          lok (j) = .true.
c
          if (j .eq. nat) goto 800
c
c ... if it is a new type, check if any of the others should get
c     the same type
c
          do l=j+1,nat
            if (atmtyp(l) .ne. atmtyp(i)) goto 850
            if (nbr(l) .ne. nbr(j)) goto 850
c
            do k=1,nbr(l)
              check3 (k) = ityp (nbrptr(k,l))
            end do
ccc            call isort (nbr(l),check3)
            if (nbr(i) .gt. 0) call hsorti (nbr(l),check3)
c
            do k=1,nbr(j)
              if (check2(k) .ne. check3(k)) goto 850
            end do
c
c ... yes, this one looks similar to the new one
c
            newtyp (l) = newtyp (j)
            nwityp (l) = nwityp (j)
            write (*,4735) atmnam(l),newtyp(l)
            lok (l) = .true.
c
  850       continue
          end do
c
c ... now that (at least) one atom's type has changed, re-check all
c
  800     continue
        end do
      end do
      call ivalut (' Nr of new types :',1,i1)
c
      do i=1,nat
        atmtyp (i) = newtyp (i)
        ityp (i)   = nwityp (i)
      end do
c
      write (*,*)
      call jvalut (' Nr of atoms :',1,nat)
      call jvalut (' Nr of types :',1,ntyp)
c
      if (lthird) goto 1234
c
      lthird = .true.
c
      write (*,*)
      write (*,*) 'At present, atoms are classified *only* based on'
      write (*,*) 'the number & types of their neighbours and the'
      write (*,*) 'bond lengths'
      write (*,*) 'If you want, I can make the typing more strict,'
      write (*,*) 'by also comparing neighbours of neighbours'

      answer = 'N'
      call textin (' Do a third round of typing ?',answer)
      call upcase (answer)
      call textut (' Do a third round of typing :',answer)
      if (answer .ne. 'Y') goto 1234
c
      write (*,*)
      call prompt (' Third round of atom typing ...')
c
      goto 888
c
c      do i=1,nat
c        if (nbr(i) .gt. 0) then
c          do j=1,nbr(i)
c            nbrtyp(j,i) = ityp(nbrptr(j,i))
c          end do
c        end if
c      end do
c
c      do i=1,nat-1
c
c	print *,' >>> ',atmnam(i),' ',atmtyp(i)
c
c        do j=i+1,nat
c          if (atmtyp(i) .ne. atmtyp(j)) goto 225
c          if (nbr(i) .le. 0 .and. nbr(j) .le. 0) goto 225
c
c	print *,' >>>>>> ',atmnam(j),' ',atmtyp(j)
c
c          do k=1,nbr(i)
c            kk = nbrptr(k,i)
c            if (nbr(kk) .le. 0) goto 224
c            do l=1,nbr(kk)
c              ll = nbrptr(l,kk)
c              if (ll .eq. i) goto 223
c              xang = angle (ll,kk,i,atmxyz)
c              xtyp = atmtyp (ll)
c
c ... check if similar type neighbour's neighbour with similar angle
c
c              do k1=1,nbr(j)
c                kk1 = nbrptr(k1,j)
c                if (nbr(kk1) .le. 0) goto 223
c                do l1=1,nbr(kk1)
c                  ll1 = nbrptr(l1,kk1)
c                  if (ll1 .eq. j) goto 223
c                  if (atmtyp(ll1).eq.xtyp .and.
c     +    abs(angle(ll1,kk1,j,atmxyz)-xang) .le. samean) goto 221
c                end do
c                goto 222
c  221           continue
c              end do
c
c              goto 223
c
c 	print *,atmnam(i),'|',atmtyp(i)
c 	print *,'not matched :'
c 	print *,xang,' ',xtyp
c
c ... found a discrepancy
c
c  222         continue
c              ntyp = ntyp + 1
c              write (line,'(a2,a1,i4)') chem(j),id,ntyp
c              call remspa (line)
c              ityp (j) = ntyp
c              first (ntyp) = j
c              atmtyp (j) = line(1:4)
c              typnam (ntyp) = line(1:4)
c              write (*,4735) atmnam(j),atmtyp(j)
c              goto 225
c
c  223         continue
c            end do
c  224       continue
c          end do
c  225     continue
c        end do
c      end do
c
c ... finished atom typing
c
 1234 continue
      write (*,*)
      write (*,*)
     +   'Making sure all atom types have unique names ...'
      call uniqat (ntyp,typnam)
c
      write (*,*)
      write (*,*) 'Finished classifying atoms'
      call prompt (' Current atom types :')
      do i=1,nat
c
        atmtyp (i) = typnam (ityp(i))
c
        write (*,4740) i,atmnam(i),atmtyp(i)
      end do
c
      write (*,*)
      call jvalut (' Nr of atoms :',1,nat)
      call jvalut (' Nr of types :',1,ntyp)
c
c --- TOPOLOGY FILE
c
      close (ifile)
      call xopxua (ifile,tfile,.true.,ierr)
      if (ierr .ne. 0) return
c
      line = 'Remarks '//tfile
      write (ifile,'(a)') line(1:leng1(line))
c
      call stamp (line)
      line = 'Remarks '//line
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks Auto-generated by XPLO2D from file '//filenm
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks You *MUST* check/edit MASSes and CHARges !!!'
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks Check DONOrs and ACCEptors'
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks Verify IMPRopers yourself'
      write (ifile,'(a)') line(1:leng1(line))
c
      line =
     +  'Remarks DIHEdrals which are not flat are commented out'
      write (ifile,'(a)') line(1:leng1(line))
c
      write (ifile,*)
      write (ifile,'(a)') ' set echo=false end'
      write (ifile,*)
      write (ifile,'(a)') ' { Note: edit masses if necessary }'
      do i=1,ntyp
        j = first(i)
c        if (atomnr(j) .eq. 6) then
c          write (line,'(a,a4,1x,f10.5,a,i2,a)')
c     +      ' MASS ',typnam(i),12.01100+numh(i)*1.008,
c     +      ' ! 12.01100 + 1.008 * ',numh(i),' (Hs)'
c        else if (atomnr(j) .eq. 8) then
c          write (line,'(a,a4,1x,f10.5,a,i2,a)')
c     +      ' MASS ',typnam(i),15.99940+numh(i)*1.008,
c     +      ' ! 15.99940 + 1.008 * ',numh(i),' (Hs)'
c        else if (atomnr(j) .eq. 7) then
c          write (line,'(a,a4,1x,f10.5,a,i2,a)')
c     +      ' MASS ',typnam(i),14.00670+numh(i)*1.008,
c     +      ' ! 14.00670 + 1.008 * ',numh(i),' (Hs)'
c        else if (atomnr(j) .eq. 16) then
c          write (line,'(a,a4,1x,f10.5,a,i2,a)')
c     +      ' MASS ',typnam(i),32.06000+numh(i)*1.008,
c     +      ' ! 32.06000 + 1.008 * ',numh(i),' (Hs)'
c        else if (atomnr(j) .eq. 1) then
c          write (line,'(a,a4,1x,f10.5,a)')
c     +      ' MASS ',typnam(i),'  1.00800 ! Hydrogen'
c        else
c
c .. write the mass of this element
c
          write (line,'(a,a4,1x,f10.5,a,a2,a,1x,f10.5,a,i2,a)')
     +      ' MASS ',typnam(i),amass(j)+numh(i)*1.008,
     +      ' ! assuming ',chem(j),
     +      ' -> ',amass(j),' + 1.008 * ',numh(i),' (Hs)'
          call pretty (line(23:))
c        end if
        write (ifile,'(a)') line(1:leng1(line))
      end do
c
c MASS SMX1  150.36000 ! assuming this is Samarium ???   150.36000 + 1.008 for each H
c12345678901234567890123456789012345678901234567890123456789012345678901234567890123456
c
      write (ifile,*)
      write (ifile,'(a)') ' autogenerate angles=true end'
      write (ifile,*)
c
      line = 'RESIdue '//resnam
      write (ifile,'(a)') line(1:leng1(line))
      write (ifile,*)
      write (ifile,'(a)')
     +  ' { Note: electrostatics should normally not be used in }'
      write (ifile,'(a)')
     +  ' { crystallographic refinement since it can produce }'
      write (ifile,'(a)')
     +  ' { artefacts. For this reason, all charges are set to }'
      write (ifile,'(a)')
     +  ' { zero by default. Edit them if necessary }'
      write (ifile,'(a)') 'GROUp'
      do i=1,nat
        line = ' ATOM '//atmnam(i)//'  TYPE '//atmtyp(i)//
     +         '  CHARge  0.0  END'
        if (lhybr .or. loccu) then
          write (line(leng1(line)+1:),'(a,i2)')
     +      ' ! Nr of Hs = ',nint(qatom(i))
        end if
        write (ifile,'(a)') line(1:leng1(line))
      end do
c
      if (nat .lt. 2) goto 6510
c
c ... TO DO: add hydrogens if LOCCU and polar ?
c     lot of book-keeping ! (need to add bonds and bond/angle params
c     as well and donors ...)
c
      write (ifile,*)
      line = ' '
      i1 = 2
      do i=1,nat-1
        do j=i+1,nat
          if (conect(i,j)) then
            line(i1:) = 'BOND ' // atmnam(i) //
     +                  ' ' // atmnam(j)
            i1 = i1 + 20
            if (i1 .ge. 65) then
              write (ifile,'(a)') line(1:leng1(line))
              line = ' '
              i1 = 2
            end if
          end if
        end do
      end do
      if (i1 .gt. 2) write (ifile,'(a)') line(1:leng1(line))
c
      ndtyp = 0
      write (ifile,*)
      write (ifile,'(a)')
     +  ' { Note: edit these DIHEdrals if necessary }'
      do i=1,nat-1
        do j=i+1,nat
          if (.not.conect(i,j) .or. i.eq.j) goto 300
          if (isdihe(i,j) .or. isdihe(j,i)) goto 300
          do k=1,nat
            if (.not.conect(i,k) .or. j.eq.k .or. i.eq.k) goto 301
            do l=1,nat
              if (.not.conect(j,l) .or. l.eq.i .or.
     +            l.eq.j .or. l.eq.k) goto 302
c
              xang = tangle(k,i,j,l,atmxyz)
              call fixang (xang)
              if (xang .lt. -150.) xang = xang + 360.
c
              l0 = (abs(xang) .le. flat)
              l180 = (abs(xang+180.0) .le. flat .or.
     +                abs(xang-180.0) .le. flat)
              l60 = (abs(xang-60.0)  .le. sixty .or.
     +               abs(xang-90.0)  .le. sixty .or.
     +               abs(xang-120.0) .le. sixty .or.
     +               abs(xang+60.0)  .le. sixty .or.
     +               abs(xang+90.0)  .le. sixty .or.
     +               abs(xang+120.0) .le. sixty)
c
              if (l60 .or. l0 .or. l180) then
c
                if (l180) then
                  line = '  DIHEdral '//atmnam(k)//' '//atmnam(i)//
     +                   ' '//atmnam(j)//' '//atmnam(l)//
     +                   ' ! flat ? (180 degrees = trans)'
                else if (l0) then
                  line = '  DIHEdral '//atmnam(k)//' '//atmnam(i)//
     +                   ' '//atmnam(j)//' '//atmnam(l)//
     +                   ' ! flat ? (0 degrees = cis)'
                else
                  line = '! DIHEdral '//atmnam(k)//' '//atmnam(i)//
     +                   ' '//atmnam(j)//' '//atmnam(l)//
     +                   ' ! flexible dihedral ???'
                end if
c
                write (line(length(line)+2:),'(f8.2)') xang
                write (ifile,'(a)') line(1:leng1(line))
c
              end if
c
c ... flag this bond as done
c
                isdihe (i,j) = .true.
                isdihe (j,i) = .true.
c
c ... keep track of parameters
c
                if (ndtyp .gt. 0) then
                  do m=1,ndtyp
                    if ( (tdih(1,m) .eq. atmtyp(k) .and.
     +                    tdih(2,m) .eq. atmtyp(i) .and.
     +                    tdih(3,m) .eq. atmtyp(j) .and.
     +                    tdih(4,m) .eq. atmtyp(l)) ) then
                      ndih (m) = ndih (m) + 1
                      sdih (m) = sdih (m) + sin (xang*degtor)
                      cdih (m) = cdih (m) + cos (xang*degtor)
                      midih (m) = min (midih(m), xang)
                      madih (m) = max (madih(m), xang)
                      if (nextra .gt. 0) then
                        do ix=1,nextra
                          xang = tangle(k,i,j,l,extxyz(1,1,ix))
                          call fixang (xang)
                          if (xang .lt. -150.) xang = xang + 360.
                          ndih (m) = ndih (m) + 1
                          sdih (m) = sdih (m) + sin (xang*degtor)
                          cdih (m) = cdih (m) + cos (xang*degtor)
                          midih (m) = min (midih(m), xang)
                          madih (m) = max (madih(m), xang)
                        end do
                      end if
                      goto 310
                    else if ( (tdih(1,m) .eq. atmtyp(l) .and.
     +                    tdih(2,m) .eq. atmtyp(j) .and.
     +                    tdih(3,m) .eq. atmtyp(i) .and.
     +                    tdih(4,m) .eq. atmtyp(k)) ) then
c
ccc                      xang = -xang
ccc                      if (xang .lt. -150.) xang = xang + 360.
c
                      ndih (m) = ndih (m) + 1
                      sdih (m) = sdih (m) + sin (xang*degtor)
                      cdih (m) = cdih (m) + cos (xang*degtor)
                      midih (m) = min (midih(m), xang)
                      madih (m) = max (madih(m), xang)
                      if (nextra .gt. 0) then
                        do ix=1,nextra
                          xang = tangle(k,i,j,l,extxyz(1,1,ix))
                          call fixang (xang)
                          if (xang .lt. -150.) xang = xang + 360.
                          ndih (m) = ndih (m) + 1
                          sdih (m) = sdih (m) + sin (xang*degtor)
                          cdih (m) = cdih (m) + cos (xang*degtor)
                          midih (m) = min (midih(m), xang)
                          madih (m) = max (madih(m), xang)
                        end do
                      end if
                      goto 310
                    end if
                  end do
                end if
c
                ndtyp = ndtyp + 1
                ndih (ndtyp) = 1
                sdih (ndtyp) = sin (xang*degtor)
                cdih (ndtyp) = cos (xang*degtor)
                midih (ndtyp) = xang
                madih (ndtyp) = xang
                tdih (1,ndtyp) = atmtyp(k)
                tdih (2,ndtyp) = atmtyp(i)
                tdih (3,ndtyp) = atmtyp(j)
                tdih (4,ndtyp) = atmtyp(l)
                udih (1,ndtyp) = atmnam(k)
                udih (2,ndtyp) = atmnam(i)
                udih (3,ndtyp) = atmnam(j)
                udih (4,ndtyp) = atmnam(l)
                if (nextra .gt. 0) then
                  m = ndtyp
                  do ix=1,nextra
                    xang = tangle(k,i,j,l,extxyz(1,1,ix))
                    call fixang (xang)
                    if (xang .lt. -150.) xang = xang + 360.
                    ndih (m) = ndih (m) + 1
                    sdih (m) = sdih (m) + sin (xang*degtor)
                    cdih (m) = cdih (m) + cos (xang*degtor)
                    midih (m) = min (midih(m), xang)
                    madih (m) = max (madih(m), xang)
                  end do
                end if
c
  310           continue
c
cxyz              end if
  302         continue
            end do
  301       continue
          end do
  300     continue
        end do
      end do
c
      nityp = 0
      write (ifile,*)
      write (ifile,'(a)')
     +  ' { Note: edit these IMPRopers if necessary }'
c
      do i=1,nat
c
ccc        if (atomnr(i) .ne. 6 .and. atomnr(i) .ne. 7) goto 400
c
c ... generate impropers for all 3 and 4 nbr Cs and Ns
c
        if (nbr(i) .lt. 3) goto 400
c
c ... looks like a tetrahedral carbon (or an aromatic branch which must
c     be flat, or a carboxylate C, ditto)
c
        i1 = 1
        i2 = 2
        i3 = 3
c
c        if (xang .lt. 0.0) then
c          i1 = 3
c          i2 = 2
c          i3 = 1
c          xang =
c     +      tangle (i,nbrptr(i1,i),nbrptr(i2,i),nbrptr(i3,i),atmxyz)
c        end if
c
c ... make sure the three first atoms are not co-linear !!!
c
        qdummy = angle (i,nbrptr(i1,i),nbrptr(i2,i),atmxyz)
        qdummy = min (abs(qdummy),abs(qdummy-180.0),abs(qdummy+180.0))
        if (qdummy .le. 5.0) then
ccc          print *,i,nbrptr(i1,i),nbrptr(i2,i),qdummy
          i2 = 3
          i3 = 2
          qdummy = angle (i,nbrptr(i1,i),nbrptr(i2,i),atmxyz)
          qdummy = min (abs(qdummy),abs(qdummy-180.0),abs(qdummy+180.0))
ccc          print *,i,nbrptr(i1,i),nbrptr(i2,i),qdummy
          if (qdummy .le. 5.0) then
            call errcon ('Four co-linear atoms ?????')
          end if
        end if
c
c ... NOTE: should 4 nbrs give rise to 4 IMPR perhaps ???
c           or do the ANGLes keep the others in check ???
c
        xang =
     +    tangle (i,nbrptr(i1,i),nbrptr(i2,i),nbrptr(i3,i),atmxyz)
c
        line = ' IMPRoper '//atmnam(i)//' '//atmnam(nbrptr(i1,i))//
     +       ' '//atmnam(nbrptr(i2,i))//' '//atmnam(nbrptr(i3,i))//
     +       ' ! chirality or flatness improper '
        write (line(length(line)+2:),'(f8.2)') xang
        write (ifile,'(a)') line(1:leng1(line))
c
c ... keep track of parameters
c
        if (nityp .gt. 0) then
          do j=1,nityp
            if (timp(1,j) .eq. atmtyp(i)) then
              if ( timp (2,j) .eq. atmtyp(nbrptr(i1,i)) .and.
     +             timp (3,j) .eq. atmtyp(nbrptr(i2,i)) .and.
     +             timp (4,j) .eq. atmtyp(nbrptr(i3,i))) then
                nimp (j) = nimp (j) + 1
                simp (j) = simp (j) + sin (xang*degtor)
                cimp (j) = cimp (j) + cos (xang*degtor)
                miimp (j) = min (miimp(j), xang)
                maimp (j) = max (maimp(j), xang)
                if (nextra .gt. 0) then
                  do ix=1,nextra
                    xang = tangle (i,nbrptr(i1,i),nbrptr(i2,i),
     +                             nbrptr(i3,i),extxyz(1,1,ix))
                    nimp (m) = nimp (m) + 1
                    simp (m) = simp (m) + sin (xang*degtor)
                    cimp (m) = cimp (m) + cos (xang*degtor)
                    miimp (m) = min (miimp(m), xang)
                    maimp (m) = max (maimp(m), xang)
                  end do
                end if
                goto 400
              end if
            end if
          end do
        end if
c
        nityp = nityp + 1
        nimp (nityp) = 1
        simp (nityp) = sin (xang*degtor)
        cimp (nityp) = cos (xang*degtor)
        miimp (nityp) = xang
        maimp (nityp) = xang
        timp (1,nityp) = atmtyp(i)
        timp (2,nityp) = atmtyp(nbrptr(i1,i))
        timp (3,nityp) = atmtyp(nbrptr(i2,i))
        timp (4,nityp) = atmtyp(nbrptr(i3,i))
        uimp (1,nityp) = atmnam(i)
        uimp (2,nityp) = atmnam(nbrptr(i1,i))
        uimp (3,nityp) = atmnam(nbrptr(i2,i))
        uimp (4,nityp) = atmnam(nbrptr(i3,i))
        if (nextra .gt. 0) then
          m = nityp
          do ix=1,nextra
            xang = tangle (i,nbrptr(i1,i),nbrptr(i2,i),
     +                     nbrptr(i3,i),extxyz(1,1,ix))
            nimp (m) = nimp (m) + 1
            simp (m) = simp (m) + sin (xang*degtor)
            cimp (m) = cimp (m) + cos (xang*degtor)
            miimp (m) = min (miimp(m), xang)
            maimp (m) = max (maimp(m), xang)
          end do
        end if
c
  400   continue
      end do
c
 6510 continue
      write (ifile,*)
      write (ifile,'(a)')
     +  ' { Note: edit any DONOrs and ACCEptors if necessary }'
c
      do i=1,nat
c
c ... OXYGEN
c
        if (atomnr(i) .eq. 8) then
c
c ... donor ?
c
          if ( (loccu .or. lhybr) .and. qatom(i).gt.0.1) then
            i1 = nint (qatom(i))
            do i2=1,i1
              write (line,'(a,i1,1x,a)') '! DONOr H?',i2,atmnam(i)
            end do
          else if (.not. (loccu .or. lhybr) ) then
            if (nbr(i) .le. 1) then
              line = '! DONOr H? '//atmnam(i)//
     +               ' ! only true if -OHx (x>0)'
              write (ifile,'(a)') line(1:leng1(line))
            end if
          end if
c
c ... acceptor
c
          if (nbr(i) .le. 0) then
            ante = '"  "'
          else
            ante = atmnam(nbrptr(1,i))
          end if
          line = ' ACCEptor '//atmnam(i)//' '//ante
          write (ifile,'(a)') line(1:leng1(line))
c
c ... NITROGEN
c
        else if (atomnr(i) .eq. 7) then
c
c ... donor ?
c
          if ( (loccu .or. lhybr) .and. qatom(i).gt.0.1) then
            i1 = nint (qatom(i))
            do i2=1,i1
              write (line,'(a,i1,1x,a)') '! DONOr H?',i2,atmnam(i)
              write (ifile,'(a)') line(1:leng1(line))
            end do
          else if (.not. (loccu .or. lhybr) ) then
            if (nbr(i) .le. 2) then
              line = '! DONOr H? '//atmnam(i)//
     +               ' ! only true if -NHx (x>0)'
              write (ifile,'(a)') line(1:leng1(line))
            end if
          end if
        end if
c
      end do
c
      write (ifile,*)
      write (ifile,'(3a)') 'END { RESIdue ',resnam,' }'
      write (ifile,*)
c
      close (ifile)
c
c --- TNT FILE
c
      if (ltnt) then
c
ccc        close (kfile)
        call xopxua (kfile,tntfile,.true.,ierr)
        if (ierr .ne. 0) return
c
        line = '! Filename = '//tntfile
        write (kfile,'(a)') line(1:leng1(line))
c
        call stamp (line)
        line = '! '//line
        write (kfile,'(a)') line(1:leng1(line))
c
        line = '! Auto-generated by XPLO2D from file '//filenm
        write (kfile,'(a)') line(1:leng1(line))
c
        write (kfile,'(20(a/))')
     +  '!',
     +  '! Note: impropers are mapped to TRIGONAL (-17 to +17),',
     +  '!       or CHIRAL (+-55 to +-17) restraints. Values',
     +  '!       outside these ranges yield an error message',
     +  '!',
     +  '! Note: PLANES are derived from 0 or +- 180 torsions',
     +  '!       and are implemented as sets of 4-atom planes.',
     +  '!       You can merge these into larger planes manually.',
     +  '!       Also note that cis and trans can NOT be',
     +  '!       distinguished by a PLANE restraint !  If this',
     +  '!       is important, use a TORSION 1000 or 1180',
     +  '!       restraint instead of the corresponding PLANE !',
     +  '!',
     +  '! Note: non-flat torsion angles lead to a TORSION',
     +  '!       restraint which is commented out by a # sign.',
     +  '!       To include such restraints in your refinement,',
     +  '!       remove the # sign(s). Also note that all such',
     +  '!       torsions by default have a periodicity of 1 !',
     +  '!'
c
      end if
c
c --- ONO FILE
c
      if (lono) then
c
ccc        close (lfile)
        call xopxua (lfile,onofile,.true.,ierr)
        if (ierr .ne. 0) return
c
        line = '! Filename = '//onofile
        write (lfile,'(a)') line(1:leng1(line))
c
        call stamp (line)
        line = '! '//line
        write (lfile,'(a)') line(1:leng1(line))
c
        line = '! Auto-generated by XPLO2D from file '//filenm
        write (lfile,'(a)') line(1:leng1(line))
c
        write (lfile,'(20(a:/))')
     +  '!',
     +  '!!!!! ADD THIS FILE TO YOUR .bonds_angles DATABLOCK  !!!!!',
     +  '!!!!! THEN UPDATE THE .bonds_angles DATABLOCK HEADER !!!!!',
     +  '!',
     +  '! --------------------------------------------------------',
     +  '!'
c
        write (lfile,'(2a)') 'residue ',resnam
c
c ... CENTRE = central atom:
c
c     - if ' CA ' present, use it
c     - if ' P  ' present, use it
c     - if ' C1 ' present, use it
c     - otherwise, use first atom
c
        do i=1,nat
          if (atmnam (i) .eq. ' CA ') then
            write (lfile,'(a/a,a))')
     +        '!','centre ',atmnam(i)
            goto 4110
          end if
        end do
c
        do i=1,nat
          if (atmnam (i) .eq. ' P  ') then
            write (lfile,'(a/a,a))')
     +        '!','centre ',atmnam(i)
            goto 4110
          end if
        end do
c
        do i=1,nat
          if (atmnam (i) .eq. ' C1 ') then
            write (lfile,'(a/a,a))')
     +        '!','centre ',atmnam(i)
            goto 4110
          end if
        end do
c
        write (lfile,'(a/a,a))')
     +    '!','centre ',atmnam(1)
c
 4110   continue
c
c ... ATOM = list of atoms
c
        write (lfile,'(a)') '!'
c
        do i=1,nat,10
          if (i .eq. 1) then
            line = 'atom'
          else
            line = ' '
          end if
          j = min (i+9, nat)
          if (j .lt. nat) then
            extrac = '\\'
          else
            extrac = ' '
          end if
          ll = max(4,length(line)+2)
          write (line(ll:),'(11(1x,a4))') (atmnam(k),k=i,j),extrac
          call pretty (line(ll:))
          write (lfile,'(a)') line(1:length(line))
        end do
c
c ... FRAGMENT_ALL = list of atoms
c
        write (lfile,'(a)') '!'
c
        do i=1,nat,10
          if (i .eq. 1) then
            line = 'fragment_all'
          else
            line = ' '
          end if
          j = min (i+9, nat)
          if (j .lt. nat) then
            extrac = '\\'
          else
            extrac = ' '
          end if
          ll = max(4,length(line)+2)
          write (line(ll:),'(11(1x,a4))') (atmnam(k),k=i,j),extrac
          call pretty (line(ll:))
          write (lfile,'(a)') line(1:length(line))
        end do
c
c ... FRAGMENT_SC = list of atoms
c
        write (lfile,'(a)') '!'
c
        do i=1,nat,10
          if (i .eq. 1) then
            line = 'fragment_sc'
          else
            line = ' '
          end if
          j = min (i+9, nat)
          if (j .lt. nat) then
            extrac = '\\'
          else
            extrac = ' '
          end if
          ll = max(4,length(line)+2)
          write (line(ll:),'(11(1x,a4))') (atmnam(k),k=i,j),extrac
          call pretty (line(ll:))
          write (lfile,'(a)') line(1:length(line))
        end do
c
c ... SIDE-CHAIN = list of atoms
c
        write (lfile,'(a)') '!'
c
        do i=1,nat,10
          if (i .eq. 1) then
            line = 'side-chain'
          else
            line = ' '
          end if
          j = min (i+9, nat)
          if (j .lt. nat) then
            extrac = '\\'
          else
            extrac = ' '
          end if
          ll = max(4,length(line)+2)
          write (line(ll:),'(11(1x,a4))') (atmnam(k),k=i,j),extrac
          call pretty (line(ll:))
          write (lfile,'(a)') line(1:length(line))
        end do
c
c ... CONNECT_ALL = lists of connected atoms
c
        write (lfile,'(a/a/a)') '!','! connectivity','!'
c
        do i=1,nat
          do j=1,nat
            okbond (i,j) = .false.
          end do
        end do
c
        j1 = 0
c
 4120   continue
c
        i1 = 0
c
        do i=1,nat
          do j=1,nat
            if ( conect(i,j) .and. (.not. okbond(i,j)) ) then
              j1 = j1 + 1
              line = 'connect_all '//atmnam(i)//' '//atmnam(j)
              i1 = 2
              k = j
              okbond (i,j) = .true.
              okbond (j,i) = .true.
              goto 4122
            end if
          end do
        end do
        goto 4126
c
 4122   continue
        if (i1 .ge. 12) goto 4124
        do i=1,nat
          if (conect(k,i) .and. (.not. okbond(k,i)) ) then
            line = line(1:leng1(line))//' '//atmnam(i)
            i1 = i1 + 1
            okbond (i,k) = .true.
            okbond (k,i) = .true.
            k = i
            goto 4122
          end if
        end do
c
 4124   continue
        call pretty (line)
        write (lfile,'(a)') line(1:leng1(line))
        goto 4120
c
 4126   continue
c
c ... check orphan atoms without any bonds
c
        do i=1,nat
          do j=1,nat
            if (conect(i,j)) goto 4128
          end do
          write (lfile,'(a7,1x,a,1x,a)') 'connect_all',
     +      atmnam(i),atmnam(i)
          j1 = j1 + 1
 4128     continue
        end do
c
ccc        call jvalut (' Nr of CONNECT records :',1,j1)
c
c ... TORSION = flexible torsions
c
        write (lfile,'(20(a:/))')
     +    '!',
     +    '! torsion-angle definitions',
     +    '! note: only the first 12 torsions can be used in O',
     +    '!'
c
c ... torsion datablock
c
 4130 continue
c
      ncnt = 0
      do i=1,nat
        do j=1,nat
          if (i.eq.j) goto 4140
          if (.not. conect(i,j)) goto 4140
          do k=1,nat
            if (k.eq.j .or. k.eq.i) goto 4138
            if (.not. conect(i,k)) goto 4138
            do l=1,nat
              if (l.eq.j .or. l.eq.i .or. l.eq.k) goto 4136
              if (.not. conect(l,j)) goto 4136
c
c ... we have found a dihedral
c
              xx = tangle(k,i,j,l,atmxyz)
              call fixang (xx)
              if (xx .lt. -150.) xx = xx + 360.
ccc              write (*,'(/1x,a8,4(1x,a4),1x,f8.2)')
ccc     +          'DIHEDRAL',atmnam(k),atmnam(i),
ccc     +          atmnam(j),atmnam(l),xx
c
c ... skip if torsion near 0 or +/- 180 degrees
c
              if (abs (xx) .le. flat) then
ccc                call prompt (' Skip -> torsion ~ 0')
                goto 4136
              else if (abs(xx-180) .le. flat .or.
     +                 abs(xx+180) .le. flat) then
ccc                call prompt (' Skip -> torsion ~ 180')
                goto 4136
              end if
c
c ... get all atoms affected by the torsion
c
              do i1=1,nat
                do j1=1,nat
                  okbond (i1,j1) = .false.
                end do
              end do
c
c ... flag all atoms connected to atom "J", except "I", as affected
c
              do i1=1,nat
                if (i1.ne.i) then
                  if (conect(j,i1)) okbond (i1,i1) = .true.
                end if
              end do
c
 4132         continue
              j2 = 0
              do i1=1,nat
                if (okbond(i1,i1)) then
                  do j1=1,nat
                    if (j1.ne.j) then
                      if (conect(i1,j1) .and.
     +                    (.not. okbond(j1,j1))) then
                         okbond (j1,j1) = .true.
                         j2 = j2 + 1
                      end if
                    end if
                  end do
                end if
              end do
c
c ... if I or K are flagged, we probably have a ring torsion
c
              if (okbond(i,i) .or. okbond(k,k)) then
ccc                call prompt (' Skip -> ring torsion')
                goto 4136
              end if
c
              if (j2 .gt. 0) goto 4132
c
              j2 = 1
              line = atmnam (l)
              do i1=1,nat
                if (i1 .ne. l .and. okbond(i1,i1)) then
                  line = line(1:leng1(line))//' '//
     +                    atmnam(i1)
                  j2 = j2 + 1
                end if
              end do
c
c ... if almost all atoms affected, better to use it the other
c     way around
c
ccc              call textut (' Affected atoms :',line)
              if (j2 .gt. ((nat-2)/2) ) then
ccc                call prompt (' Skip -> too many affected atoms')
                goto 4136
              end if
c
c ... check if it is a permutation of a previous torsion
c
              if (ncnt .ge. 1) then
                do i1=1,ncnt
                  if (
     +         (i.eq.deftor(2,i1) .and. j.eq.deftor(3,i1))) then
                    do j1=1,nat
                      if (okbond(j1,j1) .neqv. afftor(j1,i1))
     +                  goto 4134
                    end do
ccc                    call textut (' Permutation of :',tornam(i1))
ccc                    call prompt (' Skip -> permutation')
                    goto 4136
 4134               continue
                  end if
                end do
              end if
c
              if (ncnt .ge. maxdih) then
                call jvalut (' Max nr of torsions :',1,maxdih)
                call errcon (' Too many torsion angles !!!')
                goto 4136
              end if
c
c ... we have a new torsion !!!
c
              ncnt = ncnt + 1
              write (tornam(ncnt),'(a3,i3)') 'tor',ncnt
              call remspa (tornam(ncnt))
ccc              call textut (' OKAY ==> new torsion :',tornam(ncnt))
c
              deftor (1,ncnt) = k
              deftor (2,ncnt) = i
              deftor (3,ncnt) = j
              deftor (4,ncnt) = l
c
              valtor (ncnt) = float(nint(xx))
c
              do i1=1,nat
                afftor(i1,ncnt) = okbond(i1,i1)
              end do
c
 4136         continue
            end do
 4138       continue
          end do
 4140     continue
        end do
      end do
c
ccc      write (*,*)
ccc      call ivalut (' Nr of unique rotatable torsions :',1,ncnt)
c
c ... now write the file
c
      nl = 1
c
      do i=1,ncnt
c
ccc        write (line,'(a,1x,a,1x,f8.0,4(1x,a4))') 'torsion',
ccc     +    tornam(i),valtor(i),atmnam(deftor(1,i)),
ccc     +    atmnam(deftor(2,i)),atmnam(deftor(3,i)),
ccc     +    atmnam(deftor(4,i))
c
c ... 20040616 - target values are no longer wanted by alwyn
c
        write (line,'(a,1x,a,4(1x,a4))') 'torsion',
     +    tornam(i),atmnam(deftor(1,i)),atmnam(deftor(2,i)),
     +    atmnam(deftor(3,i)),atmnam(deftor(4,i))
c
        call pretty (line)
c
        do j=1,nat
          if (afftor(j,i)) then
            if (length(line) .ge. 60) then
              write (lfile,'(a,1x,a1)') line(1:leng1(line)),'\\'
              nl = nl + 1
              line = '   '//atmnam(j)
            else
              line = line(1:leng1(line))//' '//atmnam(j)
            end if
            call pretty (line(4:))
          end if
        end do
        write (lfile,'(a)') line(1:leng1(line))
        nl = nl + 1
      end do
c
ccc      call jvalut (' Nr of lines written :',1,nl)
c
      end if
c
c --- PARAMETER FILE
c
ccc      close (ifile)
      call xopxua (ifile,pfile,.true.,ierr)
      if (ierr .ne. 0) return
c
      line = 'Remarks '//pfile
      write (ifile,'(a)') line(1:leng1(line))
c
      call stamp (line)
      line = 'Remarks '//line
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks Auto-generated by XPLO2D from file '//filenm
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks Parameters for residue type '//resnam
      write (ifile,'(a)') line(1:leng1(line))
c
      write (ifile,*)
      write (ifile,'(a)') ' set echo=false end'
c
      if (nat .lt. 2) goto 6520
c
c ... collect all bond types
c
      nbtyp = 0
      do i=1,nat-1
        do j=i+1,nat
          if (conect(i,j)) then
            if (nbtyp .gt. 0) then
              do k=1,nbtyp
                if ( ( atmtyp(i).eq.tbon(1,k) .and.
     +                 atmtyp(j).eq.tbon(2,k) ) .or.
     +               ( atmtyp(i).eq.tbon(2,k) .and.
     +                 atmtyp(j).eq.tbon(1,k) ) ) then
                  nbon (k) = nbon (k) + 1
                  sbon (k) = sbon (k) + dismat(i,j)
                  mibon (k) = min (mibon(k), dismat(i,j))
                  mabon (k) = max (mabon(k), dismat(i,j))
                  if (nextra .gt. 0) then
                    do ix=1,nextra
                      xdis = dist (i,j,extxyz(1,1,ix))
                      nbon (k) = nbon (k) + 1
                      sbon (k) = sbon (k) + xdis
                      mibon (k) = min (mibon(k), xdis)
                      mabon (k) = max (mabon(k), xdis)
                    end do
                  end if
                  goto 410
                end if
              end do
              nbtyp = nbtyp + 1
              nbon (nbtyp) = 1
              sbon (nbtyp) = dismat(i,j)
              mibon (nbtyp) = dismat(i,j)
              mabon (nbtyp) = dismat(i,j)
              tbon (1,nbtyp) = atmtyp(i)
              tbon (2,nbtyp) = atmtyp(j)
              ubon (1,nbtyp) = atmnam(i)
              ubon (2,nbtyp) = atmnam(j)
              if (nextra .gt. 0) then
                do ix=1,nextra
                  xdis = dist (i,j,extxyz(1,1,ix))
                  nbon (nbtyp) = nbon (nbtyp) + 1
                  sbon (nbtyp) = sbon (nbtyp) + xdis
                  mibon (nbtyp) = min (mibon(nbtyp), xdis)
                  mabon (nbtyp) = max (mabon(nbtyp), xdis)
                end do
              end if
            else
              nbtyp = nbtyp + 1
              nbon (nbtyp) = 1
              sbon (nbtyp) = dismat(i,j)
              mibon (nbtyp) = dismat(i,j)
              mabon (nbtyp) = dismat(i,j)
              tbon (1,nbtyp) = atmtyp(i)
              tbon (2,nbtyp) = atmtyp(j)
              ubon (1,nbtyp) = atmnam(i)
              ubon (2,nbtyp) = atmnam(j)
              if (nextra .gt. 0) then
                do ix=1,nextra
                  xdis = dist (i,j,extxyz(1,1,ix))
                  nbon (nbtyp) = nbon (nbtyp) + 1
                  sbon (nbtyp) = sbon (nbtyp) + xdis
                  mibon (nbtyp) = min (mibon(nbtyp), xdis)
                  mabon (nbtyp) = max (mabon(nbtyp), xdis)
                end do
              end if
            end if
          end if
  410     continue
        end do
      end do
c
c ... average bond lengths and write to parameter file
c
      write (*,*)
      call jvalut (' Nr of bond types found :',1,nbtyp)
      if (nbtyp .ge. 1) then
        wgt = forcec(1)
        i1 = 0
        i2 = 0
        write (ifile,*)
        write (ifile,'(a)') ' { Note: edit if necessary }'
        if (ltnt) write (kfile,'(a/a/a)')
     +    '!','! BOND restraints','!'
        if (lono) write (lfile,'(a/a/a)')
     +    '!','! bond lengths','!'
c
        do i=1,nbtyp
          sbon (i) = sbon (i) / float(nbon(i))
c
          if (sbon(i) .le. 0.5) then
            write (line,'(a,f8.3)')
     +        '! >>> WARNING - next bond shorter than 0.5 A : ',sbon(i)
            write (ifile,'(a)') line(1:leng1(line))
            if (ltnt) write (kfile,'(a)') line(1:leng1(line))
            if (lono) write (lfile,'(a)') line(1:leng1(line))
            i2 = i2 + 1
          end if
c
          if (nbon(i) .gt. 1) then
            write (line,6110)
     +        ' BOND ',tbon(1,i),' ',tbon(2,i),' ',wgt,' ',
     +        sbon(i),' ! Nobs = ',nbon(i),' ... Range = ',
     +        mibon(i),mabon(i)
            write (tntline,6210)
     +        resnam,sbon(i),tntwgt(1),
     +        ubon(1,i),ubon(2,i),' ! Nobs = ',nbon(i),
     +        ' ... Range = ',mibon(i),mabon(i)
            write (onoline,6310)
     +        ubon(1,i),ubon(2,i),sbon(i),onowgt(1)
          else
            write (line,6110)
     +        ' BOND ',tbon(1,i),' ',tbon(2,i),' ',wgt,' ',
     +        sbon(i),' ! Nobs = ',nbon(i)
            write (tntline,6210)
     +        resnam,sbon(i),tntwgt(1),
     +        ubon(1,i),ubon(2,i),' ! Nobs = ',nbon(i)
            write (onoline,6310)
     +        ubon(1,i),ubon(2,i),sbon(i),onowgt(1)
          end if
c
          write (ifile,'(a)') line(1:leng1(line))
          if (ltnt) then
            write (kfile,'(a)') tntline(1:leng1(tntline))
            write (tntline,6215)
     +        resnam,0.0,tntwgt(5),
     +        ubon(1,i),ubon(2,i)
            write (kfile,'(a)') tntline(1:leng1(tntline))
          end if
          if (lono) then
            call pretty (onoline)
            write (lfile,'(a)') onoline(1:leng1(onoline))
          end if
c
 6110     format (a,a4,a,a4,a,f8.1,a,f6.3,a,i4,a,2f6.3)
 6210     format ('GEOMETRY ',a,' BOND ',f8.3,1x,f8.3,2x,
     +      a,' ',a,a,i4,a,2f6.3)
 6215     format ('GEOMETRY ',a,' BCORREL ',f8.3,1x,f8.3,2x,
     +      a,' ',a)
 6310     format ('bond_distance ',a,1x,a,1x,f8.3,1x,f8.3)
c
          if ((mabon(i)-mibon(i)) .ge. lrbond) then
            write (line,'(a,f8.3)')
     +        '! >>> WARNING - large range for previous bond : ',
     +        (mabon(i)-mibon(i))
            write (ifile,'(a)') line(1:leng1(line))
            if (ltnt) write (kfile,'(a)') line(1:leng1(line))
            if (lono) write (lfile,'(a)') line(1:leng1(line))
            i1 = i1 + 1
          end if
c
        end do
        if (i1 .gt. 0) then
          call ivalut (' WARNING - Bonds with large ranges :',1,i1)
        end if
        if (i2 .gt. 0) then
          call ivalut (' WARNING - Bonds shorter than 0.5 A :',1,i2)
        end if
      end if
c
c ... find all angle types
c
      natyp = 0
      do i=1,nat
        do j=1,nat-1
          if (i.eq.j) goto 500
          if (.not. conect(i,j)) goto 500
          do k=j+1,nat
            if (k.eq.i) goto 510
            if (k.eq.j) goto 510
            if (.not. conect(i,k)) goto 510
c
            if (natyp .gt. 0) then
              do l=1,natyp
                if (atmtyp(i).eq.tang(2,l)) then
                  if ( (atmtyp(j).eq.tang(1,l) .and.
     +                  atmtyp(k).eq.tang(3,l) ) .or.
     +                 (atmtyp(j).eq.tang(3,l) .and.
     +                  atmtyp(k).eq.tang(1,l) ) ) then
                    nang (l) = nang (l) + 1
                    xang = angle(j,i,k,atmxyz)
                    sang (l) = sang (l) + xang
                    miang (l) = min (miang(l), xang)
                    maang (l) = max (maang(l), xang)
                    if (nextra .gt. 0) then
                      do ix=1,nextra
                        xang = angle (j,i,k,extxyz(1,1,ix))
                        nang (l) = nang (l) + 1
                        sang (l) = sang (l) + xang
                        miang (l) = min (miang(l), xang)
                        maang (l) = max (maang(l), xang)
                      end do
                    end if
                    goto 510
                  end if
                end if
              end do
              natyp = natyp + 1
              nang (natyp) = 1
              xang = angle(j,i,k,atmxyz)
              sang (natyp) = xang
              miang (natyp) = xang
              maang (natyp) = xang
              tang (1,natyp) = atmtyp(j)
              tang (2,natyp) = atmtyp(i)
              tang (3,natyp) = atmtyp(k)
              uang (1,natyp) = atmnam(j)
              uang (2,natyp) = atmnam(i)
              uang (3,natyp) = atmnam(k)
              if (nextra .gt. 0) then
                do ix=1,nextra
                  xang = angle (j,i,k,extxyz(1,1,ix))
                  nang (natyp) = nang (natyp) + 1
                  sang (natyp) = sang (natyp) + xang
                  miang (natyp) = min (miang(natyp), xang)
                  maang (natyp) = max (maang(natyp), xang)
                end do
              end if
            else
              natyp = natyp + 1
              nang (natyp) = 1
              xang = angle(j,i,k,atmxyz)
              sang (natyp) = xang
              miang (natyp) = xang
              maang (natyp) = xang
              tang (1,natyp) = atmtyp(j)
              tang (2,natyp) = atmtyp(i)
              tang (3,natyp) = atmtyp(k)
              uang (1,natyp) = atmnam(j)
              uang (2,natyp) = atmnam(i)
              uang (3,natyp) = atmnam(k)
              if (nextra .gt. 0) then
                do ix=1,nextra
                  xang = angle (j,i,k,extxyz(1,1,ix))
                  nang (natyp) = nang (natyp) + 1
                  sang (natyp) = sang (natyp) + xang
                  miang (natyp) = min (miang(natyp), xang)
                  maang (natyp) = max (maang(natyp), xang)
                end do
              end if
            end if
  510     continue
          end do
  500   continue
        end do
      end do
c
c ... average bond angles and write to parameter file
c
      write (*,*)
      call jvalut (' Nr of angle types found :',1,natyp)
      if (natyp .ge. 1) then
        wgt = forcec(2)
        i1 = 0
        write (ifile,*)
        write (ifile,'(a)') ' { Note: edit if necessary }'
        if (ltnt) write (kfile,'(a/a/a)')
     +    '!','! ANGLE restraints','!'
        if (lono) write (lfile,'(a/a/a)')
     +    '!','! bond angles','!'
c
        do i=1,natyp
          sang (i) = sang (i) / float(nang(i))
c
          if (sang(i) .lt. 80.0) then
            write (*,*)
            call errcon ('Angle of less than 80 degrees !')
            call textut (' Between :',uang(1,i))
            call textut ('     and :',uang(2,i))
            call textut ('     and :',uang(3,i))
            nfatal = nfatal + 1
c
            write (line,'(a,f8.3)')
     +        '! >>> WARNING - next angle less than 80 degrees : ',
     +        sang(i)
            write (ifile,'(a)') line(1:leng1(line))
            if (ltnt) write (kfile,'(a)') line(1:leng1(line))
            if (lono) write (lfile,'(a)') line(1:leng1(line))
            i2 = i2 + 1
          end if
c
          if (nang(i) .gt. 1) then
            write (line,6120)
     +        ' ANGLe ',tang(1,i),' ',tang(2,i),' ',tang(3,i),' ',
     +        wgt,' ',sang(i),' ! Nobs = ',nang(i),' ... Range = ',
     +        miang(i),maang(i)
            write (tntline,6220)
     +        resnam,sang(i),tntwgt(2),
     +        uang(1,i),uang(2,i),uang(3,i),
     +        ' ! Nobs = ',nang(i),' ... Range = ',miang(i),maang(i)
            write (onoline,6320)
     +        uang(1,i),uang(2,i),uang(3,i),sang(i),onowgt(2)
          else
            write (line,6120)
     +        ' ANGLe ',tang(1,i),' ',tang(2,i),' ',tang(3,i),' ',
     +        wgt,' ',sang(i),' ! Nobs = ',nang(i)
            write (tntline,6220)
     +        resnam,sang(i),tntwgt(2),
     +        uang(1,i),uang(2,i),uang(3,i),
     +        ' ! Nobs = ',nang(i)
            write (onoline,6320)
     +        uang(1,i),uang(2,i),uang(3,i),sang(i),onowgt(2)
          end if
c
 6120     format (a,a4,a,a4,a,a4,a,f8.1,a,f8.2,a,i4,a,2f8.2)
 6220     format ('GEOMETRY ',a,' ANGLE ',f8.3,1x,f8.3,2x,
     +      a,' ',a,' ',a,a,i4,a,2f8.2)
 6320     format ('bond_angle ',a,1x,a,1x,a,1x,f8.2,1x,f8.2)
c
          write (ifile,'(a)') line(1:leng1(line))
          if (ltnt) write (kfile,'(a)') tntline(1:leng1(tntline))
          if (lono) then
            call pretty (onoline)
            write (lfile,'(a)') onoline(1:leng1(onoline))
          end if
c
          if ((maang(i)-miang(i)) .ge. lrangl) then
            write (line,'(a,f8.3)')
     +        '! >>> WARNING - large range for previous angle : ',
     +        (maang(i)-miang(i))
            write (ifile,'(a)') line(1:leng1(line))
            if (ltnt) write (kfile,'(a)') line(1:leng1(line))
            if (lono) write (lfile,'(a)') line(1:leng1(line))
            i1 = i1 + 1
          end if
c
        end do
        if (i1 .gt. 0) then
          call ivalut (' WARNING - Angles with large ranges :',1,i1)
        end if
        if (i2 .gt. 0) then
          call ivalut (' WARNING - Angles less than 80 degrees :',1,i2)
        end if
      end if
c
      if (nfatal .gt. 0 .and. (.not. lforce)) then
        call errcon ('Invalid angles - aborting !')
        return
      end if
      nfatal = 0
c
c ... average dihedrals and write to parameter file
c
      write (*,*)
      call jvalut (' Nr of dihedral types found :',1,ndtyp)
      if (ndtyp .ge. 1) then
        wgt = forcec(3)
        i1 = 0
        write (ifile,*)
        write (ifile,'(a)') ' { Note: edit if necessary }'
        if (ltnt) write (kfile,'(a/a/a)')
     +    '!','! TORSION and PLANE restraints','!'
        if (lono) write (lfile,'(a/a/a)')
     +    '!','! dihedral angles (fixed and flexible)','!'
c
        do i=1,ndtyp
          sdih (i) = sdih (i) / float(ndih(i))
          cdih (i) = cdih (i) / float(ndih(i))
          sdih (i) = rtodeg * atan2 (sdih(i), cdih(i))
          xang = sdih(i)
          sdih (i) = 30.0 * float(nint(sdih(i)/30.0))
          if (sdih(i) .lt. -179.9) sdih(i) = sdih(i) + 360.0
c
          l0 = (abs(xang) .le. flat)
          l180 = (abs(xang+180.0) .le. flat .or.
     +            abs(xang-180.0) .le. flat)
          l60 = (abs(xang-60.0)  .le. sixty .or.
     +           abs(xang-90.0)  .le. sixty .or.
     +           abs(xang-120.0) .le. sixty .or.
     +           abs(xang+60.0)  .le. sixty .or.
     +           abs(xang+90.0)  .le. sixty .or.
     +           abs(xang+120.0) .le. sixty)
c
          idum = 1000 + nint(abs(sdih(i)))
          if (sdih(i) .lt. 0) idum = - idum
c
          if (ndih(i) .gt. 1) then
            write (line,6130)
     +        ' DIHEdral ',tdih(1,i),' ',tdih(2,i),' ',tdih(3,i),' ',
     +        tdih(4,i),' ',wgt,' 0 ',sdih(i),' ! Nobs = ',ndih(i),
     +        ' ... Range = ',midih(i),madih(i)
          else
            write (line,6130)
     +        ' DIHEdral ',tdih(1,i),' ',tdih(2,i),' ',tdih(3,i),' ',
     +        tdih(4,i),' ',wgt,' 0 ',sdih(i),' ! Nobs = ',ndih(i),
     +        ' ... Value = ',xang
          end if
c
          if (ltnt) then
            if (l180 .or. l0) then
              if (ndih(i) .gt. 1) then
                write (tntline,6235)
     +            resnam,4,tntwgt(4),udih(1,i),udih(2,i),
     +            udih(3,i),udih(4,i),
     +            ' ! Torsion Nobs = ',ndih(i),
     +            ' ... Range = ',midih(i),madih(i)
              else
                write (tntline,6235)
     +            resnam,4,tntwgt(4),udih(1,i),udih(2,i),
     +            udih(3,i),udih(4,i),
     +            ' ! Torsion Nobs = ',ndih(i),
     +            ' ... Value = ',xang
              end if
            else if (l60) then
              if (ndih(i) .gt. 1) then
                write (tntline,6230)
     +            resnam,idum,tntwgt(3),udih(1,i),udih(2,i),
     +            udih(3,i),udih(4,i),' ! Nobs = ',ndih(i),
     +            ' ... Range = ',midih(i),madih(i)
              else
                write (tntline,6230)
     +            resnam,idum,tntwgt(3),udih(1,i),udih(2,i),
     +            udih(3,i),udih(4,i),' ! Nobs = ',ndih(i),
     +            ' ... Value = ',xang
              end if
              tntline = '# '//tntline
            end if
          end if
c
          if (lono) then
            if (l180 .or. l0) then
              write (onoline,6335)
     +            udih(1,i),udih(2,i),udih(3,i),
     +            udih(4,i),sdih(i),onowgt(3)
            else if (l60) then
              write (onoline,6330)
     +            udih(1,i),udih(2,i),udih(3,i),
     +            udih(4,i),sdih(i),onowgt(4)
            else
              write (onoline,6330)
     +            udih(1,i),udih(2,i),udih(3,i),
     +            udih(4,i),float(nint(xang)),onowgt(4)
            end if
          end if
c
          if (l0 .or. l180 .or. l60) then
            write (ifile,'(a)') line(1:leng1(line))
            if (ltnt) write (kfile,'(a)') tntline(1:leng1(tntline))
          end if
c
c ... O needs *all* torsions to be able to generate coordinates
c     with the build_residue command
c
          if (lono) then
            call pretty (onoline)
            write (lfile,'(a)') onoline(1:leng1(onoline))
            if (l180) then
              write (lfile,'(a,f8.2)')
     +          '! TRANS torsion - Actual value : ',xang
            else if (l0) then
              write (lfile,'(a,f8.2)')
     +          '! CIS torsion - Actual value : ',xang
            else if (l60) then
              write (lfile,'(a,f8.2)')
     +          '! 60/90/120 torsion - Actual value : ',xang
            else
              write (lfile,'(a,f8.2)')
     +          '! Descriptive torsion - Actual value : ',xang
            end if
          end if
c
 6130     format (a,a4,a,a4,a,a4,a,a4,a,f8.1,a,f8.2,a,i4,a,2f8.2)
 6230     format ('GEOMETRY ',a,' TORSION ',i5,1x,f8.2,2x,
     +      a,' ',a,' ',a,' ',a,a,i4,a,2f8.2)
 6235     format ('GEOMETRY ',a,' PLANE ',i2,1x,f8.2,2x,
     +      a,' ',a,' ',a,' ',a,a,i4,a,2f8.2)
 6330     format ('torsion_flexible',4(1x,a),2(1x,f8.2))
 6335     format ('torsion_fixed',4(1x,a),2(1x,f8.2))
c
          call fixdif (madih(i),midih(i),qdummy)
          if (qdummy .ge. lrdihe) then
            write (line,'(a,f8.3,a)')
     +        '! >>> WARNING - large range for previous dihedral : ',
     +        qdummy,' - dihedral may be flexible !!!'
            if (l0 .or. l180 .or. l60) then
              write (ifile,'(a)') line(1:leng1(line))
              if (ltnt) write (kfile,'(a)') line(1:leng1(line))
              i1 = i1 + 1
            end if
            if (lono) write (lfile,'(a)') line(1:leng1(line))
          end if
c
        end do
        if (i1 .gt. 0) then
          call ivalut (' WARNING - Dihedrals with large ranges :',1,i1)
        end if
      end if
c
c ... average impropers and write to parameter file
c
      write (*,*)
      call jvalut (' Nr of improper types found :',1,nityp)
      if (nityp .ge. 1) then
        wgt = forcec(4)
        i1 = 0
        nweird = 0
        write (ifile,*)
        write (ifile,'(a)') ' { Note: edit if necessary }'
        if (ltnt) write (kfile,'(a/a/a)')
     +    '!','! CHIRAL and TRIGONAL restraints','!'
        if (lono) write (lfile,'(a/a/a)')
     +    '!','! improper torsion angles','!'
c
        do i=1,nityp
          simp (i) = simp (i) / float(nimp(i))
          cimp (i) = cimp (i) / float(nimp(i))
          simp (i) = rtodeg * atan2 (simp(i), cimp(i))
          xang = simp(i)
c
          if (abs(simp(i)) .le. 5.0) then
            simp (i) = 0.0
c          else if (abs(simp(i)-180.0) .le. 5.0 .or.
c     +             abs(simp(i)+180.0) .le. 5.0) then
c            simp (i) = 180.0
          else if (abs(simp(i)-35.0) .le. 5.0) then
            simp (i) = 35.0
          else if (abs(simp(i)+35.0) .le. 5.0) then
            simp (i) = -35.0
          else
            if (abs(simp(i)) .le. 15.0) then
              write (line,'(a,f8.2,a)')
     +          '! >>> NOTE - unusual value for following improper : ',
     +          simp(i),' reset to 0.0'
              write (ifile,'(a)') line(1:leng1(line))
              simp (i) = 0.0
              nweird = nweird + 1
            else if (simp(i) .gt. 15.0) then
              write (line,'(a,f8.2,a)')
     +          '! >>> NOTE - unusual value for following improper : ',
     +          simp(i),' reset to +35.0'
              write (ifile,'(a)') line(1:leng1(line))
              simp (i) = 35.0
              nweird = nweird + 1
            else if (simp(i) .lt. -15.0) then
              write (line,'(a,f8.2,a)')
     +          '! >>> NOTE - unusual value for following improper : ',
     +          simp(i),' reset to -35.0'
              write (ifile,'(a)') line(1:leng1(line))
              simp (i) = -35.0
              nweird = nweird + 1
            end if
c
c            write (line,'(9a)')
c     +        '! >>> WARNING - unusual value for following improper',
c     +        ' (not near zero or +/- 35 degrees)'
c            write (ifile,'(a)') line(1:leng1(line))
c
          end if
c
          if (nimp(i) .gt. 1) then
            write (line,6140)
     +        ' IMPRoper ',timp(1,i),' ',timp(2,i),' ',timp(3,i),' ',
     +        timp(4,i),' ',wgt,' 0 ',simp(i),' ! Nobs = ',nimp(i),
     +        ' ... Range = ',miimp(i),maimp(i)
          else
            write (line,6140)
     +        ' IMPRoper ',timp(1,i),' ',timp(2,i),' ',timp(3,i),' ',
     +        timp(4,i),' ',wgt,' 0 ',simp(i),' ! Nobs = ',nimp(i),
     +        ' ... Value = ',xang
          end if
          write (ifile,'(a)') line(1:leng1(line))
 6140     format (a,a4,a,a4,a,a4,a,a4,a,f8.1,a,f8.3,a,i4,a,2f8.3)
c
          if (ltnt) then
            if (simp(i) .ge. -55.0 .and. simp(i) .le. -17.0) then
              if (nimp(i) .gt. 1) then
                write (tntline,6240)
     +            resnam,uimp(1,i),uimp(2,i),
     +            uimp(3,i),uimp(4,i),' ! Nobs = ',nimp(i),
     +          ' ... Range = ',-maimp(i),-miimp(i)
              else
                write (tntline,6240)
     +            resnam,uimp(1,i),uimp(2,i),
     +            uimp(3,i),uimp(4,i),' ! Nobs = ',nimp(i),
     +          ' ... Value = ',-xang
              end if
            else if (simp(i) .le. 55.0 .and. simp(i) .ge. 17.0) then
              if (nimp(i) .gt. 1) then
                write (tntline,6240)
     +            resnam,uimp(1,i),uimp(4,i),
     +            uimp(3,i),uimp(2,i),' ! Nobs = ',nimp(i),
     +          ' ... Range = ',miimp(i),maimp(i)
              else
                write (tntline,6240)
     +            resnam,uimp(1,i),uimp(4,i),
     +            uimp(3,i),uimp(2,i),' ! Nobs = ',nimp(i),
     +          ' ... Value = ',xang
              end if
            else if (simp(i) .le. 17.0 .and. simp(i) .ge. -17.0) then
              if (nimp(i) .gt. 1) then
                write (tntline,6245)
     +            resnam,tntwgt(4),uimp(1,i),uimp(2,i),
     +            uimp(3,i),uimp(4,i),' ! Nobs = ',nimp(i),
     +          ' ... Range = ',miimp(i),maimp(i)
              else
                write (tntline,6245)
     +            resnam,tntwgt(4),uimp(1,i),uimp(2,i),
     +            uimp(3,i),uimp(4,i),' ! Nobs = ',nimp(i),
     +          ' ... Value = ',xang
              end if
            else
              if (nimp(i) .gt. 1) then
                write (tntline,6250)
     +            '! Unusual improper value for ',resnam,
     +            uimp(1,i),uimp(2,i),uimp(3,i),uimp(4,i),
     +            ' Nobs = ',nimp(i),' ... Range = ',miimp(i),maimp(i)
              else
                write (tntline,6250)
     +            '! Unusual improper value for ',resnam,
     +            uimp(1,i),uimp(2,i),uimp(3,i),uimp(4,i),
     +            ' Nobs = ',nimp(i),' ... Value = ',xang
              end if
            end if
            write (kfile,'(a)') tntline(1:leng1(tntline))
          end if
c
          if (lono) then
            write (onoline,6335)
     +            uimp(1,i),uimp(2,i),uimp(3,i),
     +            uimp(4,i),simp(i),onowgt(3)
            call pretty (onoline)
            write (lfile,'(a)') onoline(1:leng1(onoline))
            write (lfile,'(a,f8.2)') '! Actual value : ',xang
          end if
c
 6240 format ('GEOMETRY ',a,' CHIRAL 1 1 ',
     +        a,' ',a,' ',a,' ',a,a,i4,a,2f8.2)
 6245 format ('GEOMETRY ',a,' TRIGONAL 0 ',f8.2,1x,
     +        a,' ',a,' ',a,' ',a,a,i4,a,2f8.2)
 6250 format (a,6(1x,a),i4,a,2f8.2)
c
          call fixdif (maimp(i),miimp(i),qdummy)
          if (qdummy .ge. lrimpr) then
            write (line,'(a,f8.3,a)')
     +        '! >>> WARNING - large range for previous improper : ',
     +        qdummy,' - check carefully !!!'
            write (ifile,'(a)') line(1:leng1(line))
            if (ltnt) write (kfile,'(a)') line(1:leng1(line))
            if (lono) write (lfile,'(a)') line(1:leng1(line))
            i1 = i1 + 1
          end if
c
        end do
        if (i1 .gt. 0) then
          call ivalut (' WARNING - Impropers with large ranges :',1,i1)
        end if
        if (nweird .gt. 0) then
          call ivalut (' WARNING - Unusual impropers reset :',1,nweird)
        end if
      end if
c
c ... non-bonded parameters
c
 6520 continue
c
 6969 format (a,a4,a,a,a)
c
      write (ifile,*)
      write (ifile,'(a)') ' { Note: edit if necessary }'
      do i=1,ntyp
        j = first(i)
        call elinfo (chem(j),fulnam,nr,mass,radius,.false.)
        if (chem(j) .eq. ' C') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' 0.1200  3.7418    0.1000  3.3854',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' O') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' 0.1591  2.8509    0.1591  2.8509',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' N') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' 0.2384  2.8509    0.2384  2.8509',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' S') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' 0.0430  3.3676    0.0430  3.3676',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' P') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' 0.5849  3.3854    0.5849  3.3854',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' H') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' 0.0498  1.4254    0.0498  1.4254',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'SE') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0430 3.510       .0430   3.510 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'BR') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .3200 3.884       .3200   3.884 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'CL') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .2600 3.671       .2600   3.671 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' F') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .1000 3.029       .1500   2.76  ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' I') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .8000 3.635       .8000   3.635 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'AL') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .2500 4.740       .2500   4.740 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'BA') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .1600 3.385       .1600   3.385 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'CA') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0600 3.107       .0600   3.107 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'CS') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0045 5.183       .0045   5.183 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'CU') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0600 2.459       .0600   2.459 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'FE') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0200 2.610       .0200   2.610 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. ' K') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0100 4.187       .0100   4.187 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'LI') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0070 2.202       .0070   2.202 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'MG') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0450 2.564       .0450   2.564 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'MN') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .7000 2.851       .7000   2.851 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'NA') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0260 2.940       .0260   2.940 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'NI') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0500 1.782       .0500   1.782 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'PD') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0100 2.388       .0100   2.388 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'PT') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0100 2.370       .0100   2.370 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'RB') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0065 4.741       .0065   4.741 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'RH') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0250 2.797       .0250   2.797 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'RU') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0250 2.940       .0250   2.940 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'SI') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .3100 4.455       .2000   3.92  ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'SN') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0450 3.385       .0450   3.385 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'SR') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .1719 3.523       .1719   3.523 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'YB') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0600 3.127       .0600   3.127 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'ZN') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .2500 1.942       .2500   1.942 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else if (chem(j) .eq. 'ZR') then
          write (line,6969)
     +      ' NONBonded ',typnam(i),' .0100 3.617       .0100   3.617 ',
     +      ' ! assuming ',fulnam(1:leng1(fulnam))
        else
          if (nr .gt. 0) then
            if (lforce) then
              write (line,'(a,a4,a,a,a)')
     +          ' NONBonded ',typnam(i),
     +          ' 0.1 2.0 0.1 2.0 ! >>> WARNING - No data for ',
     +          fulnam(1:leng1(fulnam)),'; values are dummy numbers !'
            else
              write (line,'(a,a4,a,a,a)')
     +          '! NONBonded ',typnam(i),
     +          ' 0.0 0.0 0.0 0.0 ! >>> WARNING - No data for ',
     +          fulnam(1:leng1(fulnam)),'; set yourself'
            end if
          else
            if (lforce) then
              write (line,'(a,a4,a,a)')
     +          ' NONBonded ',typnam(i),
     +          ' 0.1 2.0 0.1 2.0 ! >>> WARNING - Unknown element;',
     +          ' values are dummy numbers !'
            else
              write (line,'(a,a4,a,a)')
     +          '! NONBonded ',typnam(i),
     +          ' 0.0 0.0 0.0 0.0 ! >>> WARNING - Unknown element;',
     +          ' set yourself'
            end if
          end if
        end if
ccc        call pretty (line(3:))
        write (ifile,'(a)') line(1:leng1(line))
      end do
c
      write (ifile,*)
      write (ifile,'(a)') ' set echo=true end'
      write (ifile,*)
c
      close (ifile)
c
      if (ltnt) then
        if (nat .lt. 2) then
          call prompt ('0Monoatomic compound => No TNT file')
          close (unit=kfile,status='DELETE')
        else
          write (kfile,'(a)') '!'
          close (kfile)
        end if
      end if
c
      if (lono) then
        if (nat .lt. 2) then
          call prompt ('0Monoatomic compound => No O file')
          close (unit=lfile,status='DELETE')
        else
          write (lfile,'(a)') '!'
          close (lfile)
        end if
      end if
c
c --- XPLO2D PDB FILE WITH B=TYPE AND Q=#HYDROGENS
c
      call xopxua (ifile,xfile,.true.,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      write (ifile,'(9a)') 'REMARK ',line(1:leng1(line))
      write (ifile,'(a)') 'REMARK XPLO2D pseudo-PDB file'
      write (ifile,'(a)') 'REMARK B-factors   <===> atom types'
      write (ifile,'(a)') 'REMARK Occupancies <===> nr of hydrogens'
c
      do i=1,nat
c
        write (line,4713) atmnam(i),(atmxyz(j,i),j=1,3),
     +    qatom(i),float(ityp(i))
        line (73:76) = segid
        line (1:6) = 'ATOM  '
        write (line(7:11),'(i5)') i
        line (18:20) = resnam
        line (22:26) = resnum
        line (77:78) = chem(i)
        write (ifile,'(a)') line(1:leng1(line))
c
      end do
c
      write (ifile,'(a)') 'END   '
c
      close (ifile)
c
c --- MINIMISATION INPUT FILE
c
      close (ifile)
      call xopxua (ifile,mfile,.true.,ierr)
      if (ierr .ne. 0) return
c
      line = 'Remarks '//mfile
      write (ifile,'(a)') line(1:leng1(line))
c
      call stamp (line)
      line = 'Remarks '//line
      write (ifile,'(a)') line(1:leng1(line))
c
      line = 'Remarks Auto-generated by XPLO2D from file '//filenm
      write (ifile,'(a)') line(1:leng1(line))
c
      line =
     +  'Remarks Energy-minimisation input file for residue type '//
     +  resnam
      write (ifile,'(a)') line(1:leng1(line))
c
      write (ifile,*)
      write (ifile,'(a)') 'topology'
      write (ifile,'(a)') '!  @tophcsdx.pro'
      write (ifile,'(a)') '  @'//tfile(1:leng1(tfile))
      write (ifile,'(a)') 'end'
c
      write (ifile,*)
      write (ifile,'(a)') 'parameters'
      write (ifile,'(a)') '!  @parhcsdx.pro'
      write (ifile,'(a)') '  @'//pfile(1:leng1(pfile))
      write (ifile,'(a)') '  nbonds'
      write (ifile,'(a)') '    atom cdie shift eps=8.0  e14fac=0.4'
      write (ifile,'(a)') '    cutnb=7.5 ctonnb=6.0 ctofnb=6.5'
      write (ifile,'(a)') '    nbxmod=5 vswitch'
      write (ifile,'(a)') '  end'
      write (ifile,'(a)') '  remark dielectric constant eps set to 8.0'
      write (ifile,'(a)') 'end'
c
      write (ifile,*)
      write (ifile,'(a)') 'flags exclude elec ? end'
c
      write (ifile,*)
      write (ifile,'(a)') 'segment name="'//segid//'"'
      write (ifile,'(a)') '  chain'
      write (ifile,'(a)') '   coordinates @'//cfile(1:leng1(cfile))
      write (ifile,'(a)') '  end'
      write (ifile,'(a)') 'end'
      write (ifile,'(a)') 'coordinates @'//cfile(1:leng1(cfile))
c
      write (ifile,*)
      write (ifile,'(a)') '! Remarks Build hydrogens (if any)'
      write (ifile,'(a)') '! hbuild'
      write (ifile,'(a)') '!   selection = (hydrogen and not known)'
      write (ifile,'(a)') '!   phistep=45'
      write (ifile,'(a)') '! end'
c
      write (ifile,*)
      write (ifile,'(a)') 'Remarks DELETE hydrogens and unknown atoms'
      write (ifile,'(a)') 'delete'
      write (ifile,'(a)') '  selection = (hydrogen or (not known))'
      write (ifile,'(a)') 'end'
c
      write (ifile,*)
      write (ifile,'(a)')
     +  '! Remarks If you want to shake up the coordinates a bit ...'
      write (ifile,'(a)') '! do (x=x+rand(0.1)-0.05) (all)'
      write (ifile,'(a)') '! do (y=y+rand(0.1)-0.05) (all)'
      write (ifile,'(a)') '! do (z=z+rand(0.1)-0.05) (all)'
c
      write (ifile,*)
      write (ifile,'(a)')
     +  'print threshold=0.02 bonds'
      write (ifile,'(a)')
     +  'print threshold=3.0 angles'
      write (ifile,'(a)')
     +  'print threshold=3.0 dihedrals'
      write (ifile,'(a)')
     +  'print threshold=3.0 impropers'
c
      if (nat .lt. 2) then
        write (ifile,*)
        write (ifile,'(a)') 'Remarks No Powell (monoatomic compound)'
        write (ifile,'(a)') '! minimise powell'
        write (ifile,'(a)') '!   nstep=250 drop=40.0'
        write (ifile,'(a)') '! end'
      else
        write (ifile,*)
        write (ifile,'(a)') 'Remarks Do Powell energy minimisation'
        write (ifile,'(a)') 'minimise powell'
        write (ifile,'(a)') '  nstep=250 drop=40.0'
        write (ifile,'(a)') 'end'
      end if
c
      write (ifile,*)
      call locase (resnam)
      write (ifile,'(a)')
     +  'write coordinates output='//resnam//'_min.pdb end'
      write (ifile,'(a)')
     +  'write structure   output='//resnam//'.psf end'
c
      write (ifile,*)
      write (ifile,'(a)')
     +  '! constraints interaction (not hydro) (not hydro) end'
c
      write (ifile,*)
      write (ifile,'(a)')
     +  'print threshold=0.02 bonds'
      write (ifile,'(a)')
     +  'print threshold=3.0 angles'
      write (ifile,'(a)')
     +  'print threshold=3.0 dihedrals'
      write (ifile,'(a)')
     +  'print threshold=3.0 impropers'
c
      write (ifile,*)
      write (ifile,'(a)')
     +  'flags exclude * include vdw end energy end'
      write (ifile,'(a)')
     +  'distance from=(not hydro) to=(not hydro) cutoff=2.6 end'
c
      write (ifile,*)
      write (ifile,'(a)') 'stop'
      write (ifile,*)
c
      close (ifile)
c
      write (*,*)
      call prompt (' All done ...')
c
 4725 format (' Atom ',a4,' has NO neighbours !')
 4730 format (' Atom ',a4,' >',(6(' [#',i3,'] ',f5.3,:)))
 4735 format (' Atom ',a4,' gets type ',a4)
 4740 format (' Atom # ',i4,' = ',a4,' --> ',a4)
c
      return
      end
c
c
c
      subroutine phipsi ()
c
c ... generate an include file with phi,psi restraints
c
      integer maxtyp,maxidl
      parameter (maxtyp=20, maxidl=5)
c
      real ideal (2,maxidl,0:maxtyp)
      real tol,phi,psi,qphi,qpsi,qdis,xnear,stol
c
      integer nideal(0:maxtyp)
      integer i,ierr,leng1,j,k,myres,ityp,ibest,jbest,kbest,nres
c
      logical linter,xinter,lall,lgen
c
      character restyp(0:maxtyp)*3
      character reslin(0:maxtyp)*80
      character file1*80,file2*80,line*256,myseg*4,mytyp*3
c
code ...
c
      linter = xinter()
c
      tol = 60.0
c
c ... set up ideal stuff (residue type 0 = all non-gly/pro residues)
c
      reslin (0) = 'XXX 4 -67 -25 -114 144 -79 141 62 43'
      reslin (1) = 'ALA 4 -62 -29 -132 151 -76 145 61 43'
      reslin (2) = 'LEU 3 -66 -28  -97 139         66 29'
      reslin (3) = 'PRO 3 -57 -22  -61 151                -74 72'
      reslin (4) = 'THR 2 -74 -21 -106 149'
      reslin (5) = 'CYS 3 -70 -24 -107 143         64 40'
      reslin (6) = 'HIS 5 -73 -21 -128 147 -76 141 62 42 -125 71'
      reslin (7) = 'ILE 2 -65 -32 -105 135'
      reslin (8) = 'MET 4 -64 -28 -125 144 -85 137 64 38'
      reslin (9) = 'SER 4 -69 -20 -135 155 -78 151 60 48'
      reslin (10)= 'VAL 2 -66 -32 -108 139'
      reslin (11)= 'PHE 4 -68 -28 -126 149 -87 139 67 33'
      reslin (12)= 'ARG 4 -66 -27 -125 148 -83 142 63 42'
      reslin (13)= 'TYR 4 -71 -24 -121 148 -79 142 66 33'
      reslin (14)= 'TRP 4 -67 -26 -121 146 -78 140 68 24'
      reslin (15)= 'ASP 3 -71 -19          -95 133 62 44'
      reslin (16)= 'ASN 3 -76 -14         -102 132 63 43'
      reslin (17)= 'GLU 4 -64 -28 -122 143 -82 140 63 44'
      reslin (18)= 'GLN 3 -66 -27 -102 143         62 46'
      reslin (19)= 'LYS 3 -66 -27 -100 143         62 45'
      reslin (20)= 'GLY 0'
c
      do i=0,maxtyp
ccc      call textut (' Now :',reslin(i))
        restyp(i) = reslin(i)(1:3)
        read (reslin(i)(4:5),'(i2)') nideal(i)
ccc      print *,restyp(i),nideal(i)
        if (nideal(i) .gt. 0) then
          read (reslin(i)(6:),*) ((ideal(k,j,i),k=1,2),j=1,nideal(i))
        end if
      end do
c
      file1 = 'phipsi.list'
      call textin (' Name of MOLEMAN2 Phi/Psi list file ?',file1)
      call xopxoa (1,file1,linter,ierr)
      if (ierr .ne. 0) return
c
      file2 = 'phipsi.xplor'
      call textin (' Name of X-PLOR Phi/Psi restraint file ?',file2)
      call xopxua (2,file2,linter,ierr)
      if (ierr .ne. 0) return
c
      line = 'N'
      write (*,*) 'You may choose to generate restraints for All'
      write (*,*) 'suitable residues, or only for those that already'
      write (*,*) 'lie Near a favoured region.'
      call textin (' Option (A/N) ?',line)
      call upcase (line)
      call remspa (line)
      lall = (line(1:1) .eq. 'A')
      if (.not. lall) then
        call fvalin (' Allowed distance (1-180 degrees) ?',1,tol)
        tol = max (1.0, min (180., tol))
      end if
c
      line = 'R'
      write (*,*) 'You may choose to use General restraint values or'
      write (*,*) 'Residue-specific ones for the ideal phi,psi values.'
      call textin (' Option (G/R) ?',line)
      call upcase (line)
      call remspa (line)
      lgen = (line(1:1) .eq. 'G')
c
      call stamp (line)
      write (2,6010,err=9000) ' REMARK',line(1:leng1(line))
      write (2,6010,err=9000) ' REMARK Filename:',file2(1:leng1(file2))
c
      write (*,*)
      write (2,6010,err=9000)
      if (lall) then
        write (*,*) 'Generate restraints for All residues'
        write (2,6010,err=9000)
     +    ' { Generate restraints for All residues }'
      else
        write (*,*) 'Generate restraints for Near-ideal residues'
        write (2,6010,err=9000)
     +    ' { Generate restraints for Near-ideal residues }'
        call fvalut (' Tolerance (degrees) :',1,tol)
        write (2,'(a,f8.2,a)',err=9000) 
     + ' { Max distance = ',tol,' degrees }'
      end if
      stol = tol*tol
c
      if (lgen) then
        write (*,*) 'Use General ideal phi,psi values'
        write (2,6010,err=9000)
     +    ' { Use General ideal phi,psi values }'
        call fvalut (' Ideals :',2*nideal(0),ideal(1,1,0))
      else
        write (*,*) 'Use Residue-specific ideal phi,psi values'
        write (2,6010,err=9000)
     +    ' { Use Residue-specific ideal phi,psi values }'
        do i=0,maxtyp
          write (*,6000) restyp(i),
     +      ((ideal(k,j,i),k=1,2),j=1,nideal(i))
        end do
      end if
c
 6000 format (' For ',a3,' : ',15f6.0)
 6010 format (10(a,1x))
 6020 format (a)
c
      write (2,6010,err=9000)
      write (2,6010,err=9000) 
     +  ' evaluate ($phi_wt = 20.0)'
      write (2,6010,err=9000) 
     +  ' evaluate ($psi_wt = 20.0)'
c
      write (2,6010,err=9000)
      write (2,6010,err=9000) 
     +  ' set echo=off end'
c
      write (2,6010,err=9000)
      write (2,6010,err=9000) 
     +  ' { define phi,psi restraints }'
      write (2,6010,err=9000) 
     +  ' parameter'
c
c         -------
c12345678901234567890
c
      nres = 0
c
   11 continue
      read (1,6020,err=9100,end=9100) line
      if (line(1:16) .eq. '         -------') goto 10
      goto 11
c
c PRO-    1 -FCAA     -999.9     -171.3      173.0        2.1
c123456789012345678901234567890123456789012345678901234567890
c
   10 continue
      read (1,6020,err=9100,end=8000) line
c
ccc      call textut (' ...',line)
c
      if (line(2:4) .eq. '   ') goto 10
      mytyp = line(2:4)
      read (line(7:10),'(i4)',err=10,end=10) myres
      myseg = line(13:16)
      read (line(20:),*,err=10,end=10) phi,psi
c
      call textut (' Residue :',line)
c
c ... skip if no valid phi or psi
c
      if (abs(phi) .gt. 180.0 .or. abs(psi) .gt. 180.0) then
        call prompt (' SKIP - invalid PHI and/or PSI')
        goto 10
      end if
c
c ... find residue type if necessary
c
      ityp = 0
      if (.not. lgen) then
        do i=1,maxtyp
          if (mytyp .eq. restyp(i)) ityp = i
        end do
      end if
c
c ... if this residue type has no preferred values, skip it
c
      if (nideal(ityp) .le. 0) then
        call prompt (' SKIP - no preferred PHI,PSI values')
        goto 10
      end if
c
c ... find nearest ideal value
c
      xnear = 999.0 * 999.0
      ibest = -999
      jbest = -999
      kbest = -999
      do i=1,nideal(ityp)
        do j=-1,1
          qphi = phi + float(j)*360.
          do k=-1,1
            qpsi = psi + float(k)*360.
            qdis = (ideal(1,i,ityp)-qphi)**2 +
     +             (ideal(2,i,ityp)-qpsi)**2
            if (qdis .lt. xnear) then
              ibest = i
              jbest = j
              kbest = k
              xnear = qdis
            end if
          end do
        end do
      end do
c
ccc      print *,ibest,ideal(1,ibest,ityp),ideal(2,ibest,ityp),sqrt(xnear)
c
      if (ibest .le. 0 .or. jbest .lt. -1 .or.
     +    kbest .lt. -1) then
        call prompt (' SKIP - no nearest neighbour ??? HJAELP !!!')
        goto 10
      end if
c
c ... skip if not nearby (option N)
c
      if (.not. lall) then
        if (xnear .gt. stol) then
          write (*,6810) ibest,ideal(1,ibest,ityp),
     +      ideal(2,ibest,ityp),sqrt(xnear)
          call prompt (' SKIP - too far from ideal values')
          goto 10
        end if
      end if
c
      write (*,6800) ibest,ideal(1,ibest,ityp),ideal(2,ibest,ityp),
     +  sqrt(xnear)
c
      write (2,6600,err=9000) myseg,myres-1,myseg,myres,myseg,myres,
     +  myseg,myres,ideal(1,ibest,ityp)
      write (2,6700,err=9000) myseg,myres,myseg,myres,myseg,myres,
     +  myseg,myres+1,ideal(2,ibest,ityp)
c
      nres = nres + 1
c
      goto 10
c
 6600 format (/
     +  '   dihedral'/
     +  '     (name C  and segid="',a4,'" and resi ',i4,')'/
     +  '     (name N  and segid="',a4,'" and resi ',i4,')'/
     +  '     (name CA and segid="',a4,'" and resi ',i4,')'/
     +  '     (name C  and segid="',a4,'" and resi ',i4,')'/
     +  '     $phi_wt   0 ',f8.1)
c
 6700 format (/
     +  '   dihedral'/
     +  '     (name N  and segid="',a4,'" and resi ',i4,')'/
     +  '     (name CA and segid="',a4,'" and resi ',i4,')'/
     +  '     (name C  and segid="',a4,'" and resi ',i4,')'/
     +  '     (name N  and segid="',a4,'" and resi ',i4,')'/
     +  '     $psi_wt   0 ',f8.1)
c
 6800 format (' RESTRAIN to #',i1,' = ',2f8.0,' (dist=',f8.1,')')
 6810 format (' NEAREST  to #',i1,' = ',2f8.0,' (dist=',f8.1,')')
c
c ... end-of-file
c
 8000 continue
c
      write (2,6010,err=9000)
      write (2,6010,err=9000) 
     +  ' end'
      write (2,6010,err=9000)
      write (2,6010,err=9000) 
     +  ' set echo=on end'
      write (2,6010,err=9000)
      write (2,'(a,i8,a)',err=9000) 
     +  ' { Nr of phi,psi restraint pairs :',nres,' }'
      write (2,6010,err=9000)
c
      call jvalut (' Nr of restrained residues :',1,nres)
      goto 9900
c
 9000 continue
      call errcon ('While writing file')
      goto 9900
c
 9100 continue
      call errcon ('While reading file')
c
 9900 continue
      close (1)
      close (2)
c
      return
      end
c
c
c
      subroutine counth (nat,maxatm,qatom,dismat,
     +                   atmxyz,atomnr,atmnam,nbr,nbrptr,maxnbr)
c
c ... try to figure out nr of hydrogens
c
      implicit none
c
c ... single C-C = 1.541 A
c     aromatic   = 1.395 A
c     double     = 1.337 A
c     triple     = 1.204 A
c
c ... single C-N = ~1.475 A
c     double     = ~1.35
c     triple     =  1.158
c
c ... single C-O = ~1.45 A
c     double     =  1.23
c
c ... single C-S = ~1.82 A
c     double     =  1.72
c
      real ccsingle,ccdouble,cnsingle,cndouble,cosingle,cssingle
      parameter (ccsingle=1.48, ccdouble=1.28, cnsingle=1.40)
      parameter (cndouble=1.20, cosingle=1.30, cssingle=1.77)
c
      integer nat,maxatm,maxnbr
c
      real dismat(maxatm,maxatm),atmxyz(3,maxatm),qatom(maxatm)
      real imp,tangle
c
      integer nbrptr(maxnbr,maxatm),nbr(maxatm),atomnr(maxatm)
      integer i,j,nhs,j1,j2,j3,ntoth
c
      character atmnam(maxatm)*4
c
code ...
c
      write (*,*)
      call prompt (' Trying to guess number of hydrogens ...')
c
      call fvalut (' Minimum C-C single bond (A) :',1,ccsingle)
      call fvalut (' Minimum C-C double bond (A) :',1,ccdouble)
      call fvalut (' Minimum C-N single bond (A) :',1,cnsingle)
      call fvalut (' Minimum C-N double bond (A) :',1,cndouble)
      call fvalut (' Minimum C-O single bond (A) :',1,cosingle)
      call fvalut (' Minimum C-S single bond (A) :',1,cssingle)
c
      ntoth = 0
c
      do i=1,nat
        nhs = 0
        if (atomnr(i) .eq. 6) then
c
c ... CARBON
c
          if (nbr(i) .eq. 4) then
            nhs = 0
          else if (nbr(i) .eq. 0) then
            nhs = 4
          else if (nbr(i) .eq. 1) then
c
c ... Carbon with one neighbour
c
            j = nbrptr(1,i)
            nhs = 3
            if (atomnr (j) .eq. 6) then
              if (dismat(i,j) .ge. ccsingle) then
                nhs = 3
              else if (dismat(i,j) .ge. ccdouble) then
                nhs = 2
              else
                nhs = 1
              end if
            else if (atomnr (j) .eq. 7) then
              if (dismat(i,j) .ge. cnsingle) then
                nhs = 3
              else if (dismat(i,j) .ge. cndouble) then
                nhs = 2
              else
                nhs = 1
              end if
            else if (atomnr (j) .eq. 8) then
              if (dismat(i,j) .ge. cosingle) then
                nhs = 3
              else
                nhs = 2
              end if
            else if (atomnr (j) .eq. 16) then
              if (dismat(i,j) .ge. cssingle) then
                nhs = 3
              else
                nhs = 2
              end if
            end if
c
c ... Carbon with two neighbours
c
          else if (nbr(i) .eq. 2) then
            j1 = nbrptr(1,i)
            j2 = nbrptr(2,i)
            nhs = 2
            if (atomnr(j1) .eq. 6) then
              if (dismat(i,j1) .le. ccsingle) nhs = nhs - 1
            end if
            if (atomnr(j2) .eq. 6) then
              if (dismat(i,j2) .le. ccsingle) nhs = nhs - 1
            end if
            if (atomnr(j1) .eq. 7) then
              if (dismat(i,j1) .le. cnsingle) nhs = nhs - 1
            end if
            if (atomnr(j2) .eq. 7) then
              if (dismat(i,j2) .le. cnsingle) nhs = nhs - 1
            end if
            if (atomnr(j1) .eq. 8) then
              if (dismat(i,j1) .le. cosingle) nhs = nhs - 1
            end if
            if (atomnr(j2) .eq. 8) then
              if (dismat(i,j2) .le. cosingle) nhs = nhs - 1
            end if
            if (atomnr(j1) .eq. 16) then
              if (dismat(i,j1) .le. cssingle) nhs = nhs - 1
            end if
            if (atomnr(j2) .eq. 16) then
              if (dismat(i,j2) .le. cssingle) nhs = nhs - 1
            end if
            nhs = max (1,nhs)
c
c ... Carbon with three neighbours
c
          else if (nbr(i) .eq. 3) then
            j1 = nbrptr(1,i)
            j2 = nbrptr(2,i)
            j3 = nbrptr(3,i)
            imp = tangle(i,j1,j2,j3,atmxyz)
            if (abs(imp) .le. 15.0) then
              nhs = 0
            else
              nhs = 1
            end if
          end if
          write (*,6000) 'Carbon',atmnam(i),nbr(i),nhs
c
        else if (atomnr(i) .eq. 8) then
c
c ... OXYGEN
c
          if (nbr(i) .eq. 2) then
            nhs = 0
          else if (nbr(i) .eq. 0) then
            nhs = 2
          else if (nbr(i) .eq. 1) then
            j = nbrptr(1,i)
            nhs = 1
            if (atomnr(j) .eq. 6) then
              if (dismat(i,j) .le. cosingle) then
                nhs = 0
              end if
            end if
          end if
          write (*,6000) 'Oxygen',atmnam(i),nbr(i),nhs
c
        else if (atomnr(i) .eq. 16) then
c
c ... SULPHUR
c
          if (nbr(i) .eq. 2) then
            nhs = 0
          else if (nbr(i) .eq. 0) then
            nhs = 2
          else if (nbr(i) .eq. 1) then
            j = nbrptr(1,i)
            nhs = 1
            if (atomnr(j) .eq. 6) then
              if (dismat(i,j) .le. cssingle) then
                nhs = 0
              end if
            end if
          end if
          write (*,6000) 'Sulphur',atmnam(i),nbr(i),nhs
c
        else if (atomnr(i) .eq. 7) then
c
c ... NITROGEN
c
          if (nbr(i) .eq. 3) then
            nhs = 0
          else if (nbr(i) .eq. 0) then
            nhs = 3
          else if (nbr(i) .eq. 1) then
            j = nbrptr(1,i)
            nhs = 1
            if (atomnr(j) .eq. 6) then
              if (dismat(i,j) .le. cnsingle) then
                nhs = 0
              end if
            end if
          else if (nbr(i) .eq. 2) then
            j1 = nbrptr(1,i)
            j2 = nbrptr(2,i)
            nhs = 1
            if (atomnr(j1) .eq. 6) then
              if (dismat(i,j1) .le. cndouble) nhs = nhs - 1
            end if
            if (atomnr(j2) .eq. 6) then
              if (dismat(i,j2) .le. cndouble) nhs = nhs - 1
            end if
            nhs = max (0, nhs)
          end if
          write (*,6000) 'Nitrogen',atmnam(i),nbr(i),nhs
c
        end if
c
        qatom (i) = float(nhs)
        ntoth = ntoth + nhs
c
      end do
c
 6000 format (1x,a10,1x,a4,' Nbrs: ',i2,' ~Hs:',i2)
c
      call ivalut (' Est. total nr of hydrogens :',1,ntoth)
c
      return
      end
c
c
c
      subroutine uniqat (ntyp,typnam)
c
      implicit none
c
      integer ntyp,i,j,k,l,ia,iz
c
      character typnam(ntyp)*4
c
code ...
c
      ia = ichar('A')
      iz = ichar('Z')
c
      do i=2,ntyp
c
  111   continue
c
        do j=1,i-1
          if (typnam(j) .eq. typnam(i)) goto 100
        end do
        goto 999
c
c ... not unique
c
  100   continue
        call textut (' ERROR - Non-unique type name :',typnam(j))
c
        do k=4,1,-1
          do l=ia,iz
            typnam(i)(k:k) = char(l)
            do j=1,i-1
              if (typnam(j) .eq. typnam(i)) goto 200
            end do
            goto 998
  200       continue
          end do
        end do
c
        call errcon ('Could NOT fix the problem !!!')
        goto 999
c
c ... corrected the error
c
  998   continue
        call textut ('         Changed to :',typnam(i))
c
c ... check again !
c
ccc        goto 111
c
c ... done
c
  999   continue
      end do
c
      return
      end
c
c
c
      subroutine dihere
c
c ... generate DIHEdral restraints to restrain a (low resolution)
c     model to be similar to another (high resolution) one
c
c ... pointers for residue atoms:
c     1 = N
c     2 = CA
c     3 = C
c     4 = CB
c     5 = ?G?
c     6 = ?D?
c
      implicit none
c
      integer maxres, maxatm
      parameter (maxres = 10000, maxatm=6*maxres)
c
      real atmxyz(3,maxatm)
      real weight(4)
      real dist,tangle,dum,phi,psi,chi1,chi2
c
      integer iatptr(6,maxres)
      integer iresid(maxres)
      integer iunit,i,j,ierr,junit,leng1,iat,jat,nphi,npsi,nchi1
      integer nres,natoms,j1,j2,i1,i2,i3,i4,ires,nchi2
c
ccc      logical okay(maxres)
      logical xinter,linter,lhydro
c
      character segid(maxres)*4,atmnam(maxatm)*4,resnam(maxres)*10
      character filpdb*80,filepp*80
      character line*128,nowseg*4,nowres*10
c
      data iunit,junit /10, 11/
c
      data weight /20.0, 20.0, 15.0, 10.0/
c
code ...
c
      linter = xinter ()
c
      call jvalut (' Max nr of residues :',1,maxres)
c
      write (*,*)
      filpdb = 'in.pdb'
      call textin (' Input PDB file ?',filpdb)
      call xopxoa (iunit,filpdb,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('Could not open PDB file')
        return
      end if
c
      write (*,*)
      filepp = 'dihe_restr.xplor'
      call textin (' Phi/Psi/Chi-1/Chi-2 restraints file ?',filepp)
c
      write (*,*)
      call fvalin (
     +  ' Restraint weights for Phi/Psi/Chi-1/Chi-2 ?',4,weight)
c
c ... read the PDB file
c
      natoms = 0
      nres = 0
      nowseg = '?!$#'
      nowres = '??????????'
c
   10 continue
      read (iunit,'(a)',end=100,err=999) line
      if (line(1:6) .ne. 'ATOM  ') goto 10
c
      if (nowres .ne. line(18:27) .or.
     +    nowseg .ne. line(73:76)) then
c
        if (nres .gt. 1) then
          if (iat .lt. 3) nres = nres - 1
        end if
c
        nowres = line(18:27)
        nowseg = line(73:76)
c
        nres = nres + 1
        iat = 0
        do i=1,6
          iatptr (i,nres) = -1
        end do
        read (line(23:26),'(i4)') iresid(nres)
        segid (nres) = nowseg
        resnam (nres) = nowres
c
      end if
c
      jat = -1
      if (lhydro(line(13:16))) goto 50
      if (line(13:16) .eq. ' N  ') then
        jat = 1
      else if (line(13:16) .eq. ' CA ') then
        jat = 2
      else if (line(13:16) .eq. ' C  ') then
        jat = 3
      else if (line(13:16) .eq. ' CB ') then
        jat = 4
      else if (line(15:16) .eq. 'G ' .or.
     +         line(15:16) .eq. 'G1') then
        jat = 5
      else if (line(15:16) .eq. 'D ' .or.
     +         line(15:16) .eq. 'D1') then
        jat = 6
      end if
c
   50 continue
      if (jat .gt. 0) then
        iat = iat + 1
        natoms = natoms + 1
        iatptr (jat,nres) = natoms
        read (line(31:54),'(3f8.3)') (atmxyz(j,natoms),j=1,3)
        atmnam (natoms) = line(13:16)
      end if
c
cATOM      2  CA  PRO A   1      18.150  13.525  43.680  1.00 28.82      AAAA
c123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c
      goto 10
c
  100 continue
      close (iunit)
      call prompt (' PDB file read OK')
      call jvalut (' Nr of residues :',1,nres)
      call jvalut (' Nr of atoms    :',1,natoms)
c
c ... open the three output files
c
      call xopxua (junit,filepp,linter,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      write (junit,6010,err=9000)
     +  ' REMARK',line(1:leng1(line))
c
      write (junit,6010,err=9000)
     +  ' REMARK Filename:',filepp(1:leng1(filepp))
c
 6010 format (/10(a,1x))
 6020 format (/a)
 6030 format (/a,f8.1,a)
c
      write (junit,6030,err=9000) 
     +  ' evaluate ($phi_wt  = ',weight(1),
     +  ') { Use 0.0 to switch PHI   restraints off }'
      write (junit,6030,err=9000) 
     +  ' evaluate ($psi_wt  = ',weight(2),
     +  ') { Use 0.0 to switch PSI   restraints off }'
      write (junit,6030,err=9000) 
     +  ' evaluate ($chi1_wt = ',weight(3),
     +  ') { Use 0.0 to switch CHI-1 restraints off }'
      write (junit,6030,err=9000) 
     +  ' evaluate ($chi2_wt = ',weight(4),
     +  ') { Use 0.0 to switch CHI-2 restraints off }'
      write (junit,6010,err=9000) 
     +  ' set echo=off end'
      write (junit,6010,err=9000) 
     +  ' { define phi, psi, chi-1, chi-2 restraints }'
      write (junit,6010,err=9000) 
     +  ' parameter'
c
 6666 format (/
     +  '   dihedral'/
     +  '     (name ',a4,' and segid="',a4,'" and resi ',i6,')'/
     +  '     (name ',a4,' and segid="',a4,'" and resi ',i6,')'/
     +  '     (name ',a4,' and segid="',a4,'" and resi ',i6,')'/
     +  '     (name ',a4,' and segid="',a4,'" and resi ',i6,')'/
     +  '     ',a,f8.1,1x,a,a,a)
c
      nphi = 0
      npsi = 0
      nchi1 = 0
      nchi2 = 0
      call prompt (' Generating DIHEdral restraints ...')
c
c ... loop over residues
c
      do ires=1,nres
c
c ... phi ?
c
        if (ires .le. 1) goto 2000
        j1 = iatptr (2,ires-1)
        j2 = iatptr (2,ires)
        if (j1 .le. 0 .or. j2 .le. 0) goto 2000
        dum = dist(j1,j2,atmxyz)
        if (dum .gt. 4.5) goto 2000
        i1 = iatptr(3,ires-1)
        i2 = iatptr(1,ires)
        i3 = iatptr(2,ires)
        i4 = iatptr(3,ires)
        if (i1.le.0.or.i2.le.0.or.i3.le.0.or.i4.le.0) goto 2000
        phi = tangle (i1,i2,i3,i4,atmxyz)
        write (junit,6666,err=9000)
     +    atmnam(i1),segid(ires-1),iresid(ires-1),
     +    atmnam(i2),segid(ires),iresid(ires),
     +    atmnam(i3),segid(ires),iresid(ires),
     +    atmnam(i4),segid(ires),iresid(ires),
     +    '$phi_wt  0 ',phi,' { PHI   ',resnam(ires),' }'
        nphi = nphi + 1
c
c ... psi ?
c
 2000   continue
        if (ires .ge. nres) goto 3000
        j1 = iatptr (2,ires+1)
        j2 = iatptr (2,ires)
        if (j1 .le. 0 .or. j2 .le. 0) goto 3000
        dum = dist(j1,j2,atmxyz)
        if (dum .gt. 4.5) goto 3000
        i1 = iatptr(1,ires)
        i2 = iatptr(2,ires)
        i3 = iatptr(3,ires)
        i4 = iatptr(1,ires+1)
        if (i1.le.0.or.i2.le.0.or.i3.le.0.or.i4.le.0) goto 3000
        psi = tangle (i1,i2,i3,i4,atmxyz)
        write (junit,6666,err=9000)
     +    atmnam(i1),segid(ires),iresid(ires),
     +    atmnam(i2),segid(ires),iresid(ires),
     +    atmnam(i3),segid(ires),iresid(ires),
     +    atmnam(i4),segid(ires+1),iresid(ires+1),
     +    '$psi_wt  0 ',psi,' { PSI   ',resnam(ires),' }'
        npsi = npsi + 1
c
c ... chi1 ?
c
 3000   continue
        if (iatptr(5,ires) .le. 0) goto 5000
        i1 = iatptr(1,ires)
        i2 = iatptr(2,ires)
        i3 = iatptr(4,ires)
        i4 = iatptr(5,ires)
        if (i1.le.0.or.i2.le.0.or.i3.le.0.or.i4.le.0) goto 5000
        chi1 = tangle (i1,i2,i3,i4,atmxyz)
        write (junit,6666,err=9000)
     +    atmnam(i1),segid(ires),iresid(ires),
     +    atmnam(i2),segid(ires),iresid(ires),
     +    atmnam(i3),segid(ires),iresid(ires),
     +    atmnam(i4),segid(ires),iresid(ires),
     +    '$chi1_wt 0 ',chi1,' { CHI-1 ',resnam(ires),' }'
        nchi1 = nchi1 + 1
c
c ... chi2 ?
c
        i1 = iatptr(2,ires)
        i2 = iatptr(4,ires)
        i3 = iatptr(5,ires)
        i4 = iatptr(6,ires)
        if (i1.le.0.or.i2.le.0.or.i3.le.0.or.i4.le.0) goto 5000
        chi2 = tangle (i1,i2,i3,i4,atmxyz)
        write (junit,6666,err=9000)
     +    atmnam(i1),segid(ires),iresid(ires),
     +    atmnam(i2),segid(ires),iresid(ires),
     +    atmnam(i3),segid(ires),iresid(ires),
     +    atmnam(i4),segid(ires),iresid(ires),
     +    '$chi2_wt 0 ',chi2,' { CHI-2 ',resnam(ires),' }'
        nchi2 = nchi2 + 1
c
 5000   continue
c
      end do
c
      call jvalut (' Nr of PHI   restraints :',1,nphi)
      call jvalut (' Nr of PSI   restraints :',1,npsi)
      call jvalut (' Nr of CHI-1 restraints :',1,nchi1)
      call jvalut (' Nr of CHI-2 restraints :',1,nchi2)
c
      write (junit,6010,err=9000) ' end'
      write (junit,6010,err=9000) 
     +  ' set echo=off end'
      write (junit,*)
c
      call prompt (' Restraints file written')
c
 1000 continue
      close (junit)
c
c ... all done
c
      return
c
c ... read error
c
  999 continue
      call errcon ('While reading PDB file')
      close (iunit)
      return
c
c ... write error
 9000 continue
      call errcon ('While writing file')
      goto 1000
c
      end
