      program flood
c
c ... FLOOD - try & find how many water molecules fit
c             inside a cavity
c
c ... Gerard Kleywegt @ 930930
c
      implicit none
c
      character*12 prognm,vers
c
      parameter (prognm = 'FLOOD', vers = '981109/1.1.1')
c
      integer maxsiz,maxrho,maxwat
      parameter (maxsiz = 6*1024*1024)
      parameter (maxrho = 128*1024)
      parameter (maxwat = 50000)
c
      real shadow(maxsiz),rho(maxrho),cell(6),watxyz(3,maxwat)
      real cvol,vvol
c
      integer cavity(maxsiz),grid(3),origin(3),extent(3),uvw(3)
      integer buffer(maxsiz)
      integer watijk(3,maxwat),iunit,ispace,ifmt,i,npnt,n1
c
      logical incore
c
      character single*1,filnam*80
c
      equivalence (shadow(1),cavity(1))
c
      data iunit /11/

c
code ...
c
      call gkinit (prognm,vers)
c
      call jvalut (' Max size of cavity map      :',1,maxsiz)
      call jvalut (' Max nr of solvent molecules :',1,maxwat)
c
c ... what type of file ?
c
      write (*,*)
      single = 'E'
      call textin (' Read E(zd), N(ew-ezd) or M(ask) ?',single)
      call upcase (single)
      call textut (' File type :',single)
      if (single .eq. 'E') then
        ifmt = 6
        filnam = 'cavity.ezd'
      else if (single .eq. 'M') then
        ifmt = 7
        filnam = 'cavity.mask'
      else if (single .eq. 'N') then
        ifmt = 8
        filnam = 'cavity.nezd'
      else
        call errstp ('Unsupported format type')
      end if
c
c ... file name ?
c
      write (*,*)
      call textin (' Name of cavity file             ?',filnam)
      call textut (' Name of cavity file             :',filnam)
c
      call prompt (' Read map header')
      incore = .false.
      call maphdr (filnam,iunit,ifmt,origin,extent,grid,
     +  uvw,cell,ispace,rho,maxrho,shadow,maxsiz,incore)
c
      if (ispace.le.0) then
        call errstp ('While reading header')
      end if
c
      if (.not. incore) then
        call errstp ('Map not in core')
      end if
c
      call prompt (' Map read okay')
c
c ... close map properly
c
      call mapclo (iunit)
c
      if (grid(1) .ne. grid(2) .or.
     +    grid(1) .ne. grid(3) .or.
     +    grid(2) .ne. grid(3)) then
        call errcon ('Grid not the same for x,y,z')
      end if
c
      if (cell(1) .ne. cell(2) .or.
     +    cell(1) .ne. cell(3) .or.
     +    cell(2) .ne. cell(3)) then
        call errcon ('Cell axes not the same for x,y,z')
      end if
c
      if (cell(4) .ne. 90.0 .or.
     +    cell(5) .ne. 90.0 .or.
     +    cell(6) .ne. 90.0) then
        call errstp ('Cell angles must be 90.0')
      end if
c
      npnt = extent(1)*extent(2)*extent(3)
      call jvalut (' Nr of points in map :',1,npnt)
c
      call fixer (cavity,shadow,npnt)
      do i=1,npnt
        if (cavity(i) .ne. 0 .and. cavity(i) .ne. 1) then
          call jvalut (' Map point nr :',1,i)
          call jvalut (' Has value    :',1,cavity(i))
          call errstp ('Cavity MUST contain ONLY 0 and 1')
        end if
      end do
c
      n1 = 0
      do i=1,npnt
        if (cavity(i) .eq. 1) n1 = n1 + 1
      end do
      call voxvol (cell,grid,cvol,vvol)
      call jvalut (' Nr of grid points in cavity :',1,n1)
      if (n1 .lt. 1) then
        call errstp ('No points in cavity')
      end if
c
      call rvalut (' Cell volume  :',1,cvol)
      call rvalut (' Voxel volume :',1,vvol)
c
      call drown (iunit,maxwat,cavity,buffer,grid,origin,
     +  extent(1),extent(2),extent(3),cell,watxyz,watijk,
     +  n1,vvol)
c
      call gkquit ()
c
      end
c
c
c
      subroutine fixer (cavity,shadow,npnt)
c
      implicit none
c
      integer npnt,i
c
      real shadow(npnt)
      integer cavity(npnt)
c
      do i=1,npnt
        cavity(i) = nint(shadow(i))
      end do
c
      return
      end
c
c
c
      subroutine drown (iunit,maxwat,cavity,buffer,grid,origin,
     +  ext1,ext2,ext3,cell,watxyz,watijk,n1,vvol)
c
      implicit none
c
      real twopi
      parameter (twopi=6.2831853071796)
c
      integer maxwat,ext1,ext2,ext3
c
      integer cavity(ext1,ext2,ext3),buffer(ext1,ext2,ext3)
      integer watijk(3,maxwat)
      real    watxyz(3,maxwat)
c
      real cell(6),watrad,f2c(3,3),x1(3),x2(3),togrid(3),rad2
      real rad22,vvol,svol
c
      integer iunit,grid(3),origin(3),nwat,off,n1,n2
      integer i,j,k,nfirst,ierr,nmax,leng1
c
      logical okay,xinter
c
      character solfil*80,solres*3,solatm*4,macfil*80,line*80
      character molnam*6
c
code ...
c
      call orthog (cell,f2c,0)
c
      do i=1,3
        togrid(i) = 1.0 / float(grid(i))
        x1(i) = float(origin(i)-1)*togrid(i)
        togrid (i) = cell(i)*togrid(i)
      end do
      call mulmtx (f2c,x1,x2,3,3,1)
c
      write (*,*)
      watrad = 1.4
      call fvalin (' Radius for solvent molecule ?',1,watrad)
      call fvalut (' Radius for solvent molecule :',1,watrad)
      if (watrad .le. 0.01) then
        call errstp (' Radius must be > 0.01 A')
      end if
      rad2 = watrad*watrad
      rad22 = 4.0*rad2
      svol = 2.0 * twopi * (watrad**3) / 3.0
      call fvalut (' Volume of one solvent mol   :',1,svol)
c
      write (*,*)
      solfil = 'sol.pdb'
      call textin (' Output solvent PDB file     ?',solfil)
      call textut (' Output solvent PDB file     :',solfil)
      call xopxua (iunit,solfil,xinter,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening solvent PDB file')
      end if
      close (iunit)
c
      write (*,*)
      molnam = 'SOLV'
      call textin (' Solvent molecule name in O  ?',molnam)
      call upcase (molnam)
      call textut (' Solvent molecule name in O  :',molnam)
c
      write (*,*)
      solres = 'HOH'
      call textin (' Solvent residue name        ?',solres)
      call upcase (solres)
      call textut (' Solvent residue name        :',solres)
c
      write (*,*)
      solatm = ' O  '
      call textin (' Solvent atom name           ?',solatm)
      call upcase (solatm)
      call textut (' Solvent atom name           :',solatm)
c
      write (*,*)
      nfirst = 500
      call jvalin (' Nr of first solvent residue ?',1,nfirst)
      call jvalut (' Nr of first solvent residue :',1,nfirst)
      if (nfirst .lt. 1 .or. nfirst .gt. 900) then
        call errstp ('Value must be in range 1 - 900')
      end if
c
      write (*,*)
      macfil = 'sol.omac'
      call textin (' Output O macro file         ?',macfil)
      call textut (' Output O macro file         :',macfil)
      call xopxua (iunit,macfil,xinter,ierr)
      if (ierr .ne. 0) then
        call errstp ('While opening O macro file')
      end if
      close (iunit)
      write (*,*)
c
      off = 1 + nint ( watrad /
     +      min ( togrid(1), togrid(2), togrid(3)) )
c
 6010 format ('ATOM  ',i5,1x,a4,1x,a3,2x,i4,4x,3f8.3,2f6.2)
c
c ... do it
c
      call prompt (' Trying order XYZ ')
      call copmij (buffer,cavity,ext1,ext2,ext3)
      nwat = 0
      okay = .false.
c
 1000 continue
      okay = .false.
      do i=1,ext1
        do j=1,ext2
          do k=1,ext3
            if (buffer(i,j,k) .eq. 1) then
              call checkr (i,j,k,x1,x2,togrid,nwat,
     +          watxyz,watijk,rad22,rad2,off,buffer,okay,
     +          maxwat,ext1,ext2,ext3)
            end if
          end do
        end do
      end do
c
      if (okay) goto 1000
      call jvalut (' Nr of solvent molecules :',1,nwat)
c
      call solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +  svol,vvol)
c
      call xopxua (iunit,solfil,xinter,ierr)
      call stamp (line)
      write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
      j = nfirst
      do i=1,nwat
        write (iunit,6010) i,solatm,solres,j,watxyz(1,i),
     +    watxyz(2,i),watxyz(3,i),1.0,20.0
        j = j + 1
      end do
      write (iunit,'(a)') 'END'
      close (iunit)
      nmax = nwat
c
c ... do it
c
      call prompt (' Trying order XZY ')
      call copmij (buffer,cavity,ext1,ext2,ext3)
      nwat = 0
      okay = .false.
c
 1100 continue
      okay = .false.
      do i=1,ext1
        do k=1,ext3
          do j=1,ext2
            if (buffer(i,j,k) .eq. 1) then
              call checkr (i,j,k,x1,x2,togrid,nwat,
     +          watxyz,watijk,rad22,rad2,off,buffer,okay,
     +          maxwat,ext1,ext2,ext3)
            end if
          end do
        end do
      end do
c
      if (okay) goto 1100
      call jvalut (' Nr of solvent molecules :',1,nwat)
c
      if (nwat .gt. nmax) then
        call solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +    svol,vvol)
        call xopxua (iunit,solfil,xinter,ierr)
        call stamp (line)
        write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        j = nfirst
        do i=1,nwat
          write (iunit,6010) i,solatm,solres,j,watxyz(1,i),
     +      watxyz(2,i),watxyz(3,i),1.0,20.0
          j = j + 1
        end do
        write (iunit,'(a)') 'END'
        close (iunit)
        nmax = nwat
      end if
c
c ... do it
c
      call prompt (' Trying order YZX ')
      call copmij (buffer,cavity,ext1,ext2,ext3)
      nwat = 0
      okay = .false.
c
 1200 continue
      okay = .false.
      do j=1,ext2
        do k=1,ext3
          do i=1,ext1
            if (buffer(i,j,k) .eq. 1) then
              call checkr (i,j,k,x1,x2,togrid,nwat,
     +          watxyz,watijk,rad22,rad2,off,buffer,okay,
     +          maxwat,ext1,ext2,ext3)
            end if
          end do
        end do
      end do
c
      if (okay) goto 1200
      call jvalut (' Nr of solvent molecules :',1,nwat)
c
      if (nwat .gt. nmax) then
        call solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +    svol,vvol)
        call xopxua (iunit,solfil,xinter,ierr)
        call stamp (line)
        write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        j = nfirst
        do i=1,nwat
          write (iunit,6010) i,solatm,solres,j,watxyz(1,i),
     +      watxyz(2,i),watxyz(3,i),1.0,20.0
          j = j + 1
        end do
        write (iunit,'(a)') 'END'
        close (iunit)
        nmax = nwat
      end if
c
c ... do it
c
      call prompt (' Trying order YXZ ')
      call copmij (buffer,cavity,ext1,ext2,ext3)
      nwat = 0
      okay = .false.
c
 1300 continue
      okay = .false.
      do j=1,ext2
        do i=1,ext1
          do k=1,ext3
            if (buffer(i,j,k) .eq. 1) then
              call checkr (i,j,k,x1,x2,togrid,nwat,
     +          watxyz,watijk,rad22,rad2,off,buffer,okay,
     +          maxwat,ext1,ext2,ext3)
            end if
          end do
        end do
      end do
c
      if (okay) goto 1300
      call jvalut (' Nr of solvent molecules :',1,nwat)
c
      if (nwat .gt. nmax) then
        call solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +    svol,vvol)
        call xopxua (iunit,solfil,xinter,ierr)
        call stamp (line)
        write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        j = nfirst
        do i=1,nwat
          write (iunit,6010) i,solatm,solres,j,watxyz(1,i),
     +      watxyz(2,i),watxyz(3,i),1.0,20.0
          j = j + 1
        end do
        write (iunit,'(a)') 'END'
        close (iunit)
        nmax = nwat
      end if
c
c ... do it
c
      call prompt (' Trying order ZXY ')
      call copmij (buffer,cavity,ext1,ext2,ext3)
      nwat = 0
      okay = .false.
c
 1400 continue
      okay = .false.
      do k=1,ext3
        do i=1,ext1
          do j=1,ext2
            if (buffer(i,j,k) .eq. 1) then
              call checkr (i,j,k,x1,x2,togrid,nwat,
     +          watxyz,watijk,rad22,rad2,off,buffer,okay,
     +          maxwat,ext1,ext2,ext3)
            end if
          end do
        end do
      end do
c
      if (okay) goto 1400
      call jvalut (' Nr of solvent molecules :',1,nwat)
c
      if (nwat .gt. nmax) then
        call solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +    svol,vvol)
        call xopxua (iunit,solfil,xinter,ierr)
        call stamp (line)
        write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        j = nfirst
        do i=1,nwat
          write (iunit,6010) i,solatm,solres,j,watxyz(1,i),
     +      watxyz(2,i),watxyz(3,i),1.0,20.0
          j = j + 1
        end do
        write (iunit,'(a)') 'END'
        close (iunit)
        nmax = nwat
      end if
c
c ... do it
c
      call prompt (' Trying order ZYX ')
      call copmij (buffer,cavity,ext1,ext2,ext3)
      nwat = 0
      okay = .false.
c
 1500 continue
      okay = .false.
      do k=1,ext3
        do j=1,ext2
          do i=1,ext1
            if (buffer(i,j,k) .eq. 1) then
              call checkr (i,j,k,x1,x2,togrid,nwat,
     +          watxyz,watijk,rad22,rad2,off,buffer,okay,
     +          maxwat,ext1,ext2,ext3)
            end if
          end do
        end do
      end do
c
      if (okay) goto 1500
      call jvalut (' Nr of solvent molecules :',1,nwat)
c
      if (nwat .gt. nmax) then
        call solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +    svol,vvol)
        call xopxua (iunit,solfil,xinter,ierr)
        call stamp (line)
        write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
        j = nfirst
        do i=1,nwat
          write (iunit,6010) i,solatm,solres,j,watxyz(1,i),
     +      watxyz(2,i),watxyz(3,i),1.0,20.0
          j = j + 1
        end do
        write (iunit,'(a)') 'END'
        close (iunit)
        nmax = nwat
      end if
c
      write (*,*)
      call jvalut (' Max nr of solvent molecules :',1,nmax)
c
c ... write O macro
c
      call xopxua (iunit,macfil,xinter,ierr)
      call stamp (line)
      write (iunit,'(a1,1x,a)') '!',line(1:leng1(line))
c
      write (line,*) 's_a_i ',solfil(1:leng1(solfil)),' ',molnam
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      write (line,*) 'mol ',molnam,' zo ; end'
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      write (line,*) 'centre_zone ',molnam,' ',
     +  nfirst,(nfirst+nmax-1)
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      if (solatm(1:2) .eq. ' C') then
        write (line,*) 'db_set_data .cpk_radii 6 6 ',watrad
      else if (solatm(1:2) .eq. ' N') then
        write (line,*) 'db_set_data .cpk_radii 7 7 ',watrad
      else if (solatm(1:2) .eq. ' O') then
        write (line,*) 'db_set_data .cpk_radii 8 8 ',watrad
      else if (solatm(1:2) .eq. ' S') then
        write (line,*) 'db_set_data .cpk_radii 16 16 ',watrad
      else
        write (line,*) 'db_set_data .cpk_radii ; ',watrad
      end if
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      write (line,*) 'sketch_cpk ',molnam,' bell message Done'
      call pretty (line)
      write (iunit,'(a)') line(1:leng1(line))
c
      close (iunit)
c
      return
      end
c
c
c
      subroutine checkr (i,j,k,x1,x2,togrid,nwat,watxyz,
     +                   watijk,rad22,rad2,off,buffer,okay,
     +                   maxwat,ext1,ext2,ext3)
c
      implicit none
c
      integer maxwat,ext1,ext2,ext3
c
      integer buffer(ext1,ext2,ext3)
      integer watijk(3,maxwat)
      real    watxyz(3,maxwat)
c
      real x1(3),x2(3),togrid(3),rad2,dist,rad22,d1,d2
c
      integer nwat,off,i1,j1,k1,i,j,k,l
c
      logical okay
c
code ...
c
      x1(1) = x2(1) + float(i)*togrid(1)
      x1(2) = x2(2) + float(j)*togrid(2)
      x1(3) = x2(3) + float(k)*togrid(3)
      if (nwat .gt. 0) then
        do l=1,nwat
          dist = (x1(1)-watxyz(1,l))**2 +
     +           (x1(2)-watxyz(2,l))**2 +
     +           (x1(3)-watxyz(3,l))**2
ccc          print *,l,dist,rad22
          if (dist .lt. rad22) goto 1100
        end do
      end if
      nwat = nwat + 1
      watijk (1,nwat) = i
      watijk (2,nwat) = j
      watijk (3,nwat) = k
      do l=1,3
        watxyz(l,nwat) = x1(l)
      end do
c
ccc      write (*,6000) nwat,i,j,k,x1(1),x1(2),x1(3)
 6000 format (' SOL # ',i6,' IJK= ',3i5,' XYZ = ',3f8.3)
c
      do i1=i-off,i+off
        if (i1 .lt. 1 .or. i1 .gt. ext1) goto 1210
        x1(1) = x2(1) + float(i1)*togrid(1)
        d1 = (x1(1)-watxyz(1,nwat))**2
        do j1=j-off,j+off
          if (j1 .lt. 1 .or. j1 .gt. ext2) goto 1220
          x1(2) = x2(2) + float(j1)*togrid(2)
          d2 = (x1(2)-watxyz(2,nwat))**2
          do k1=k-off,k+off
            if (k1 .lt. 1 .or. k1 .gt. ext3) goto 1230
            x1(3) = x2(3) + float(k1)*togrid(3)
            dist = d1 + d2 + (x1(3)-watxyz(3,nwat))**2
            if (dist .le. rad2) buffer (i1,j1,k1) = -1
 1230       continue
          end do
 1220     continue
        end do
 1210   continue
      end do
c               
      okay = .true.
c
 1100 continue
c
      return
      end
c
c
c
c
      subroutine solly (cavity,buffer,ext1,ext2,ext3,n1,n2,nwat,
     +                  svol,vvol)
c
      implicit none
c
      integer ext1,ext2,ext3
c
      integer cavity(ext1,ext2,ext3),buffer(ext1,ext2,ext3)
c
      real svol,vvol,tvol
c
      integer n1,n2,nwat,i,j,k
c
code ...
c
      n2 = 0
      do i=1,ext1
        do j=1,ext2
          do k=1,ext3
            if (cavity(i,j,k) .eq. 1) then
              if (buffer(i,j,k) .eq. -1) n2 = n2 + 1
            end if
          end do
        end do
      end do
c
      tvol = float(nwat)*svol
      call gvalut (' Total volume of solvent molecules :',1,tvol)
      call jvalut (' Cavity points available           :',1,n1)
      call gvalut (' Volume (A3)   :',1,(float(n1)*vvol))
      call jvalut (' Cavity points occupied by solvent :',1,n2)
      call gvalut (' Volume (A3)   :',1,(float(n2)*vvol))
      call fvalut (' Utilisation % :',1,(100.0*float(n2)/float(n1)))
      write (*,*)
c
      return
      end
