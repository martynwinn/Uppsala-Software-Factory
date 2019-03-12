c
c
c
      subroutine maphdr (file, in, flag, 
     $                   ioxyz, imxyz, mygrid, iuvw, cell, spgrp,
     $                   rho, size, savrho, sizsav, incore)
c
c ---	A routine to read in 'header' information from the
c	density maps
c ---	For CCP4 maps, if order is wrong, it gets stored away in-core
c ---	Alwyn Jones 3-Aug-91
c
c ---	Alterations
c
c ---	Phil Evans   MRC LMB, Cambridge                August 1991
c	Added option to read formatted CCP4 maps and 
c	added common block /cmscal/
c
c ---   Gerard Kleywegt March/April 1993
c       - added MASK format
c       - improved error trapping
c       - added NEWEZD format
c       - added binary XPLOR format
c
c --- gjk @ 940212 - interface with new MAPCLO routine
c                    to properly close read maps
c
      implicit none
c
      character file*(*)
      integer size, sizsav
      integer nfast,nmed,nslow
      integer cmnd, imxyz(3), in, ioxyz(3), iuvw(3), mygrid(3), 
     $  spgrp, flag
      real cell(6), rho(size), savrho(sizsav)
      real*8 drcell(6)
      logical incore,xinter
c
c --- For CCP4 maps
c
      integer nu, nv
      logical lform
      real cscale,cplus
      common /cmscal/ cscale,cplus,nu,nv,lform
      save /cmscal/
c
c ---	For Protein maps
c
      integer title(20),nrho
      character*4 title_c(20) ! mrh   

      real rhorng(6)
      integer ibp0, iperm(3), igrid(3), mifd(9)
c
      character tit*80, fmt*80, line*80
      character qxyz(3)
      integer titnum
      integer lmode, i, ierror, ipt, ip1, ip2, ir1, ir2, isec, 
     $  j, nsym, numbay
      integer ct, ct1,ierr
      real scale, avea,sdva,vara,xmin,xmax,xtot
c
      integer*2 jj1,jj2,jj3,jj4,jj5
c
ccc      data  title /20*'    '/
      data  title_c /20*'    '/ ! mrh   
      equivalence (title, title_c) ! mrh

      data qxyz /'X','Y','Z'/
c
code ...
c
      call textut (' Input map :',file)
c
      cmnd = iabs(flag)
      incore = .false.
      spgrp = 1
c
      call setmgt (in,cmnd)
c
c --- "PROTEIN" style map
c
      if (cmnd .eq. 1) then
c ---	Rewind the electron density file.
        close (in)
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        rhorng(1) = 10000000.0
        rhorng(2) = -10000000.0
        rhorng(3) = 0.0
        rhorng(4) = 0.0
        rhorng(5) = 0.0
        rhorng(6) = 0.0
        nrho	= 0
        ibp0	= 0
c ---	  Read the type 0 header.
        call getr0 (in, title, ierror)
        if (ierror .ne. 0) goto 8900
        write(6,10) title
10      format (/,' Map title is :',20A4,/)
c ---	  Read the type 1 header and symmetry card.
        call getr1 (in, cell, nsym, ierror)
        if (ierror .ne. 0) goto 8900
c ---	  Read the type 2 map index parameter record
        call getr2 (in, iperm, igrid, mifd, ierror)
        if (ierror .ne. 0) goto 8900
c ---	  Load common blocks
        do 100 i=1,3
          ipt = (i-1)*3+1
          ioxyz(iperm(i)) = mifd(ipt)
          imxyz(iperm(i)) = mifd(ipt+1)- mifd(ipt)+1
100       mygrid(iperm(i)) = igrid(i)
        iuvw(1)=iperm(2)
        iuvw(2)=iperm(1)
        iuvw(3)=iperm(3)
c
	else if (cmnd .eq. 2 .or. cmnd .eq. 3) then
c
c ---   Ten-eycke style maps: cmnd = 2 for FFT-Y, 3 for TENEYCK2
c
        close(in)
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
c ---	  read title and as much information as possible
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        call ivalin (' Fast, medium, slow axes ?',3,iuvw)
cc        write (6,20)
cc        read (5,*) iuvw
        read (in,err=8900) title
        if (cmnd .eq. 2) then
          read (in,err=8900) isec, ip1, ip2, ir1, ir2	! SECN, FAST, MEDIUM AXES
        else					!TENEYCK2
          read (in,err=8900) jj1, jj2, jj3, jj4, jj5
          isec = jj1
          ir1 = jj2
          ir2 = jj3
          ip1 = jj4
          ip2 = jj5
          rewind (in)
c ---       Skip over title
          read (in,err=8900)
        end if
        mygrid (1) = 10
        mygrid (2) = 10
        mygrid (3) = 10
        call ivalin (' Grid ?',3,mygrid)
cc        write (6,30)
cc        read (5,*) mygrid
        numbay = 10
        call ivalin (' Number of Y-sections ?',1,numbay)
cc	  write (6,40)
cc        read (5,*) numbay
        cell (1) = 100.0
        cell (2) = 100.0
        cell (3) = 100.0
        cell (4) =  90.0
        cell (5) =  90.0
        cell (6) =  90.0
        call fvalin (' Unit cell constants ?',6,cell)
cc        write (6,80)
cc        read (5, *) cell
c
        ioxyz(iuvw(1))=ip1
        ioxyz(iuvw(2))=ir1
        ioxyz(iuvw(3))=isec
        imxyz(iuvw(1))=ip2-ip1+1
        imxyz(iuvw(2))=ir2-ir1+1
        imxyz(iuvw(3))=numbay
c
      else if (cmnd .eq. 4) then
c
c ---    CCP4 maps, formatted or binary
c
c Is it formatted or binary? lform = .true. if formatted
c
        call opnmfl(file,lform,in,ierr)
        if (ierr .ne. 0) goto 9000
c
        if (lform) then
c Formatted, so read formatted header
          call rdfhdr (in, file, tit, numbay, iuvw, mygrid,
     $        isec,ip1,ip2,ir1,ir2,cell,spgrp,lmode,
     $        rhorng(1),rhorng(2),rhorng(3),rhorng(4),
     $        cscale, cplus)
c
c ... flag file type as FORMATTED CCP4
c
          call setmgt (in,-4)
c
        else
c Binary
          call mrdhdr (in, file, tit, numbay, iuvw, mygrid,
     $        isec,ip1,ip2,ir1,ir2,cell,spgrp,lmode,
     $        rhorng(1),rhorng(2),rhorng(3),rhorng(4))
        endif
cxyz
c        call prompt ('MRDHDR done')
c
c ---	 Switch these arrays around, so ordered according to x,y,z.
c
        ioxyz(iuvw(1)) = ip1
        ioxyz(iuvw(2)) = ir1
        ioxyz(iuvw(3)) = isec
        imxyz(iuvw(1)) = ip2-ip1+1
        imxyz(iuvw(2)) = ir2-ir1+1
        imxyz(iuvw(3)) = numbay
        nu = ip2-ip1+1
        nv = ir2-ir1+1
c
        if (iuvw(1) .eq. 2 .and. iuvw(2) .eq. 1 .and. iuvw(3) .eq. 3)
     $    goto 130
        if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) goto 130
        do 110 i=1,numbay
          if (lform) then
            call rdfsec(in, ierror, rho, nu, nv, cscale, cplus)
          else
c
c ... 940225 - replace MGULP by MGULPR so it can read CCP4 masks
c
            call mgulpr (in, rho, ierror)
cxyz
c        call prompt ('MGULPR done')
c
            if (ierror .ne. 0) goto 8900
          endif
cxyz
c        call prompt ('PCKRHO call next')
          call pckrho (savrho, imxyz(1), imxyz(2), imxyz(3), i,
     $          rho, imxyz(iuvw(1)), imxyz(iuvw(2)), iuvw)
cxyz
c        call prompt ('PCKRHO done')
c
110     continue
        incore = .true.
c
c ---	X-plor style giant maps. Assume sorted ZYX
c     010513 - also AMBER-style xplor maps (type=23)
c
      else if (cmnd .eq. 5 .or. cmnd .eq. 23) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        read (in, 1020,err=8900) i
        do 140 j=1,i
          read (in,1030,err=8900) tit
          call textut (' Title :',tit)
140     continue
        read (in, 1050,err=8900) (mygrid(i), ioxyz(i), imxyz(i), i=1,3)
        do 150 i=1,3
150       imxyz(i) = imxyz(i)- ioxyz(i) + 1
        read (in, 1060,err=8900) cell
        read (in, 1030,err=8900) tit
        if (tit(1:3) .ne. 'ZYX') then
          call errcon (' Map should be sorted ZYX !')
          goto 9000
        end if
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
c
c ---	 OLDEZDensity style giant maps
c
      else if (cmnd .eq. 6) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        read (in, 2010,err=8900) ioxyz
        read (in, 2010,err=8900) imxyz
        read (in, 2010,err=8900) mygrid
        read (in, 2020,err=8900) cell
        read (in, 2030,err=8900) line
        if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) then
          call errcon (' Map too big !')
          call jvalut (' Requested size :',1,
     +                 (imxyz(1)*imxyz(2)*imxyz(3)))
          call jvalut (' Available size :',1,sizsav)
          goto 9000
        end if
c
        if (line(1:4) .ne. 'MAP ') then
          call errcon ('Fifth line of OLDEZD header corrupted')
          goto 9000
        end if
        ct = 4
2100    if (line(ct:ct) .eq. ' ') then
          ct = ct+1
          if (ct .gt. 80) then
            call errcon ('Fifth line of OLDEZD header corrupted')
            goto 9000
          end if
          goto 2100
        end if
        ct1 = ct
2110    if (line(ct:ct) .ne. ')') then
          ct = ct+1
          if (ct .gt. 80) then
            call errcon ('Fifth line of OLDEZD header corrupted')
            goto 9000
          end if
          goto 2110
        end if
        fmt = line(ct1:ct)
        line(1:ct) = ' '
        call remspa (line)
cc        call alignl (line)
        read (line, 2060,err=8900) scale
        call rvalut (' Scale constant :',1,scale)
cc        type*,'scale constant=', scale
        ct = imxyz(1)*imxyz(2)*imxyz(3)
        read (in, fmt, end=8900, err=8900) (savrho(i), i=1,ct)
        do 2120 i=1, ct
2120      savrho(i) = savrho(i)/scale
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        incore = .true.
        call xstats (savrho,ct,avea,sdva,xmin,xmax,xtot)
        vara = sdva*sdva
        write (*,6010) ct,xmin,xmax,xtot,avea,vara,sdva
c
        goto 130
c
c ---	 MASKs
c
      else if (cmnd .eq. 7) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        call maskin (in,savrho,ioxyz,imxyz,mygrid,cell,sizsav,ierr)
        if (ierr .ne. 0) goto 9000
        ct = imxyz(1)*imxyz(2)*imxyz(3)
c
c ... MASKIN has written integers into array SAVRHO
c     convert them back into "real integers"
c
        do i=1,ct
          call r2r (savrho(i),savrho(i))
        end do
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        incore = .true.
        goto 130
c
c ---	 (NEW)EZDensity style giant maps
c
      else if (cmnd .eq. 8) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        do i=1,3
          ioxyz(i) = 0
          imxyz(i) = 0
          mygrid(i) = 0
          cell(i) = 0.0
          cell(i+3) = 0.0
        end do
        scale = 0.0
c
        read (in,'(a)',err=8900) line
        call upcase (line)
        if (line (1:7) .ne. 'EZD_MAP') then
          call errcon ('Not a (NEW)EZD file !')
          goto 8900
        end if
c
 2873   continue
        read (in,'(a)',err=8900) line
c
        if (line(1:1) .eq. '!') then
          call textut (' >',line)
          goto 2873
        end if
c
        call upcase (line)
        if (line(1:4).eq.'CELL') then
          read (line(5:),*,err=8900) cell
          goto 2873
        end if
c
        if (line(1:6).eq.'ORIGIN') then
          read (line(7:),*,err=8900) ioxyz
          goto 2873
        end if
c
        if (line(1:6).eq.'EXTENT') then
          read (line(7:),*,err=8900) imxyz
          goto 2873
        end if
c
        if (line(1:4).eq.'GRID') then
          read (line(5:),*,err=8900) mygrid
          goto 2873
        end if
c
        if (line(1:5).eq.'SCALE') then
          read (line(6:),*,err=8900) scale
          call rvalut (' Scale constant :',1,scale)
          goto 2873
        end if
c
        if (line(1:3) .eq. 'MAP') goto 2876
c
        call errcon ('Invalid line in (NEW)EZD file')
        call textut (' >',line)
        goto 2873
c
 2876   continue
        if (scale.eq.0.0 .or.
     +      imxyz(3).eq.0 .or. cell(6).eq.0.0) then
          call errcon ('Not all data in header')
          goto 8900
        end if
c
        if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) then
          call errcon (' Map too big !')
          call jvalut (' Requested size :',1,
     +                 (imxyz(1)*imxyz(2)*imxyz(3)))
          call jvalut (' Available size :',1,sizsav)
          goto 9000
        end if
c
        ct = imxyz(1)*imxyz(2)*imxyz(3)
        read (in, *, end=8900, err=8900) (savrho(i), i=1,ct)
        do 2920 i=1, ct
2920      savrho(i) = savrho(i)/scale
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        incore = .true.
        call xstats (savrho,ct,avea,sdva,xmin,xmax,xtot)
        vara = sdva*sdva
        write (*,6010) ct,xmin,xmax,xtot,avea,vara,sdva
c
        goto 130
c
c ---	BINARY X-plor style giant maps. Assume sorted ZYX
c
      else if (cmnd .eq. 9) then
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        read (in,err=8900) i,(tit(1:80),j=1,i)
        call textut (' Title :',tit)
c
        read (in,err=8900)
     +    mygrid(1),ioxyz(1),imxyz(1),
     +    mygrid(2),ioxyz(2),imxyz(2),
     +    mygrid(3),ioxyz(3),imxyz(3)
cc      print *,(mygrid(i), ioxyz(i), imxyz(i), i=1,3)
        do 139 i=1,3
139       imxyz(i) = imxyz(i) - ioxyz(i) + 1
cc        read (in) i
cc      print *,' I = ',i
        read (in,err=8900) drcell(1),drcell(2),drcell(3),
     +    drcell(4),drcell(5),drcell(6)
cc      print *,'... CELL ...',(drcell(i),i=1,6)
        do i=1,6
          cell(i)=drcell(i)
        end do
        read (in,err=8900) (tit(i:i),i=1,3)
cc      print *,tit(1:3)
        if (tit(1:3) .ne. 'ZYX') then
          call errcon (' Map should be sorted ZYX !')
          goto 9000
        end if
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
c
c ... check if not created on a 64-bit machine
c
        do i=1,3
          if (mygrid(i).le.1 .or. imxyz(i).le.1 .or.
     +        cell(i).le.1.0 .or. cell(i+3).le.1.0) then
            write (6,70) ioxyz,imxyz,mygrid,
     +        cell,(qxyz(iuvw(j)),j=1,3)
            call errcon ('Map probably created on a 64-bit machine')
            call errcon ('Cannot handle this !')
            goto 9000
          end if
        end do
c
        else if (cmnd .eq. 13) then
c
c ---   TNT style maps: cmnd = 13
c
        close(in)
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
c ---     read title and as much information as possible
c         Title records start with a '*'
c         except the last one...
c      
        titnum = 0
141     read(in) tit
        titnum = titnum +1
        if (tit(1:1).eq.'*') goto 141
        read(in) isec, ir1, ir2, ip1, ip2, nfast, nmed, nslow
c
c  Now we reposition the map file at the first record that contains
c  electron density information
c
        rewind(in)
        do 142 i=1,titnum
           read(in) tit
           call textut (' Title :',tit)
142     continue
        iuvw(1) = 1
        iuvw(2) = 2
        iuvw(3) = 3
        call ivalin (' Fast, medium, slow axes ?',3,iuvw)
        numbay = 10
        call ivalin (' Number of sections ?',1,numbay)
        cell (1) = 100.0
        cell (2) = 100.0
        cell (3) = 100.0
        cell (4) =  90.0
        cell (5) =  90.0
        cell (6) =  90.0
        call fvalin (' Unit cell constants ?',6,cell)
c
        mygrid(1) = nfast
        mygrid(2) = nmed
        mygrid(3) = nslow
        ioxyz(iuvw(1))=ir1
        ioxyz(iuvw(2))=ip1
        ioxyz(iuvw(3))=isec
        imxyz(iuvw(1))=ir2-ir1+1
        imxyz(iuvw(2))=ip2-ip1+1
        imxyz(iuvw(3))=numbay
c
c --- Command not found...
c
      else
        call errcon ('Unknown MAP format')
        call ivalut (' Format :',1,cmnd)
cc        write(6,60) cmnd
        goto 9000
      endif
c
c --- Come here to output some facts about map...
c
 130  continue
      write (6,70)ioxyz,imxyz,mygrid,cell,(qxyz(iuvw(i)),i=1,3)
c
      if (imxyz(iuvw(1))*imxyz(iuvw(2)) .gt. size) then
        call errcon (' Levels contain too many points')
        call jvalut (' Available :',1,size)
        call jvalut (' Requested :',1,
     +               (imxyz(iuvw(1))*imxyz(iuvw(2))))
        goto 9000
      end if
c
      if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) then
        call errcon (' Map too big !')
        call jvalut (' Requested size :',1,
     +               (imxyz(1)*imxyz(2)*imxyz(3)))
        call jvalut (' Available size :',1,sizsav)
        goto 9000
      end if
c
      return
c
c ... error traps
c
 8800 call errcon ('While opening input file')
      call mapclo (in)
      goto 9000
c
 8900 call errcon ('While reading input file')
      call mapclo (in)
      goto 9000
c
 9000 continue
      if (flag .gt. 0) then
        call errstp ('MAPHDR - Sorry !')
      end if
      call errcon ('While reading map header. Sorry !')
      spgrp = -1
c
      return
c
70    format(/,' Parameters as read from the map file',/,
     $ ' Origin ...................... ',3I10,/,
     $ ' Extent ...................... ',3I10,/,
     $ ' Grid ........................ ',3I10,/,
     + ' Cell axes ................... ',3f10.2/
     + ' Cell angles ................. ',3f10.2/
     $ ' UVW (fast, medium, slow) .... ',3a10,/)
1020  format (/,i8)
1030  format (a)
1040  format (1x,a)
1050  format (9i8)
1060  format (6e12.5)
2010  format(3i5)
2020  format(6f10.0)
2030  format (a)
2040  format (' No space to store EZD map')
2050  format (' Fifth line wrong')
2060  format (e20.6)
2070  format (' Error reading map file')
c
 6010 format (' Size ',i10/
     +        ' ED min, max, total ',1p,3e12.4/
     +        ' ED ave, var, stdev ',3e12.4/)
c
      end
