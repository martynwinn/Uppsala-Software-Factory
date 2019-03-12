c
c
c
      subroutine mapin  (file, in, flag, 
     $            ioxyz, imxyz, mygrid, iuvw, cell, 
     $            rho, size, savrho, sizsav, incore, level)
c
c ---	Read one level of a map section into array
c ---	Cleaned up 3-Aug-91 Alwyn Jones
c
      implicit none
c
      character file*80
      integer cmnd, imxyz(3), in, ioxyz(3), iuvw(3), level, mygrid(3), 
     $ size, sizsav
      real cell(6), rho(size), savrho(sizsav)
c
c ... patch for binary xplor files
c
      integer maxpbf
      parameter (maxpbf=512000)
      real*8 xpbuf(maxpbf)
c
c --- For CCP4 maps
      integer nu, nv
      logical incore, lform
      real cscale,cplus
      common /cmscal/ cscale,cplus,nu,nv,lform
      save /cmscal/
c
      integer i, ibp0, ipt, j, junk, nlev, npoints, nrows, nrho,k
      real rhorng(6)
      integer iperm(3), mifd(9)
c
      integer ierror,flag
c
      logical nd6
c
code ...
c
      cmnd = abs (flag)
c
cc      print *,' Cmnd = ',cmnd
c
c --- "PROTEIN" style map
c
      if (cmnd .eq. 1) then
c
        iperm(1) = iuvw(2)
        iperm(2) = iuvw(1)
        iperm(3) = iuvw(3)
        do i=1,3
          ipt = (i-1)*3+1
          mifd(ipt) = ioxyz (iperm(i))
          mifd(ipt+1) = imxyz (iperm(i))+ mifd(ipt)-1
        end do
        call getr5 (in, mifd, size, rho, rhorng, nrho, ibp0, ierror)
        if (ierror .ne. 0) goto 1000
c
c --- Ten Eyck style maps: icommand = 2 for FFT-Y, 3 for TENEYCK2
c --- TNT maps also...
c
      else if ((cmnd.eq.2).or.(cmnd.eq.3).or.(cmnd.eq.13)) then
c
        nlev = imxyz(iuvw(1))*imxyz(iuvw(2))
        if (cmnd .eq. 2) then	   !FFT-Y
          read (in, end=1000) (rho(i), i=1,nlev)
          read (in) junk
        else if (cmnd .eq. 13) then
          read (in, end=1000) junk
          read (in) (rho(i), i=1,nlev)
        else                       !TENEYCK2 
          read (in, end=1000) junk
          nrows = imxyz(iuvw(2))	   !IUVW(2) -> ROWS
          npoints = imxyz(iuvw(1))	 !IUVW(1) -> POINTS
          do j = 1,nrows
            read (in) (rho(i), i=(j-1)*npoints+1,j*npoints)
          end do
        end if
c
c ---    CCP4 maps.
c
      else if (cmnd .eq. 4) then
        if (incore) then
          call unprho (savrho, imxyz(1), imxyz(2), imxyz(3), level+1,
     $      rho, imxyz(2), imxyz(1), iuvw)
        else
c ---     In each subsequent pass, read one layer of the map
          if (lform) then
            call rdfsec(in, ierror, rho, nu, nv, cscale, cplus)
          else
c
c ... 940225 - replace MGULP by MGULPR so it can read CCP4 masks
c              BUT IT DOESN'T WORK AS SIMPLE AS THAT OF COURSE !
c
            call mgulpr (in, rho, ierror)
c
            if (ierror .ne. 0) goto 1000
          endif
        end if
c
      else if (cmnd .eq. 5) then
c
c --- X-plor maps. these are calculated x-fast,y-medium, z-sections
c
        read (in, 20, err=1000) i
        j = imxyz(1) * imxyz(2)
ccc        print *,' LEVEL ',i,' # ',j
        read (in, 40, err=1000) (savrho(i), i=1,j)
        call swpord (savrho, rho, imxyz(1), imxyz(2))
c
      else if (cmnd .eq. 23) then
c
c --- AMBER-style X-plor maps
c
        read (in, 20, err=1000) i
        k = 1
        j = 6 * nint(float(imxyz(1))/6.0)
        nd6 = (imxyz(1) .eq. j)
        do j=1,imxyz(2)
          read (in, 40,err=1000) (savrho(k+i), i=1,imxyz(1))
          k = k + imxyz(1)
          if (nd6) read (in,*)
        end do
        call swpord (savrho, rho, imxyz(1), imxyz(2))
c
c --- OLDEZD
c
      else if (cmnd .eq. 6) then
        call unprho (savrho, imxyz(1), imxyz(2), imxyz(3), level+1,
     $               rho, imxyz(2), imxyz(1), iuvw)
c
c --- MASK (OLD, NEW or COMPRESSED format)
c
      else if (cmnd .eq. 7) then
        call unprho (savrho, imxyz(1), imxyz(2), imxyz(3), level+1,
     $		     rho, imxyz(2), imxyz(1), iuvw)
c
c --- (NEW)EZD
c
      else if (cmnd .eq. 8) then
        call unprho (savrho, imxyz(1), imxyz(2), imxyz(3), level+1,
     $               rho, imxyz(2), imxyz(1), iuvw)
c
      else if (cmnd .eq. 9) then
c
c --- BINARY X-plor maps. these are calculated x-fast,y-medium, z-sections
c
        read (in,err=1000) i
cc        print *,' sector ',i
        j = imxyz(1)* imxyz(2)
        if (j .gt. maxpbf) then
          call errcon ('MAPIN - Buffer allocation too small')
          call jvalut (' Available :',1,maxpbf)
          call jvalut (' Required  :',1,j)
          goto 1000
        end if
        read (in,err=1000) (xpbuf(i), i=1,j)
        do i=1,j
          savrho(i)=xpbuf(i)
        end do
        call swpord (savrho, rho, imxyz(1), imxyz(2))
c
c --- Command not found...
c
      else
        call errcon ('Unknown MAP format')
        call ivalut (' Format :',1,cmnd)
        goto 1000
      endif
c
      level = level+1
      return
c
1000  continue
      if (flag .gt. 0) then
        call errstp ('MAPIN - While reading map')
      end if
      call errcon ('While reading map')
      level = -1
c
      return
c
20    format (i8)
40    format (6e12.5)
c
      end
