c
c ... subroutines for MAPMAN
c
c *** 961125 - this file now contains all subroutines which DO NOT
c              use the MAPMAN.INCL file (this was easier while
c              implementing dynamic memory allocation)
c
c
c
      subroutine gjk001 (mapje,i1,i2,i3,ipl,ext1,ext2,ext3,iunit,ierr)
c
      implicit none
c
      integer ext1,ext2,ext3
c
      real mapje(ext1,ext2,ext3)
c
      integer i1,i2,i3,ipl,iunit,ierr,j1,j2
c
code ...
c
      write (iunit,'(a)',err=9000) 'ZVALUE (1p,6e13.4)'
c
      if (i3 .eq. 1) then
        write (iunit,'(1p,6e13.4)',err=9000)
     +    ((mapje(ipl,j1,j2),j1=1,ext2),j2=1,ext3)
      else if (i3 .eq. 2) then
        write (iunit,'(1p,6e13.4)',err=9000)
     +    ((mapje(j1,ipl,j2),j1=1,ext1),j2=1,ext3)
      else if (i3 .eq. 3) then
        write (iunit,'(1p,6e13.4)',err=9000)
     +    ((mapje(j1,j2,ipl),j1=1,ext1),j2=1,ext2)
      end if
c
      ierr = 0
      return
c
 9000 continue
      ierr = -1
      return
c
      end
c
c
c
      subroutine mpstat (n,map)
c
c ... extensive map statistics
c
c ... adapted from "Numerical Recipes", pp. 454 - 459
c
      implicit none
c
      integer n
c
      real map(n)
c
      real xn,s,p,ave,adev,var,skew,curt
c
      integer i,leng1
c
      character line*80
c
code ...
c
      write (*,10) n
      if (n .lt. 2) then
        call errcon ('Too few data points')
        return
      end if
c
      xn = float (n)
      s = 0.0
      p = 0.0
      do i=1,n
        s=s+map(i)
        p=p+map(i)*map(i)
      end do
c
      ave=s/xn
      p=sqrt(p/xn)
      write (*,20) ave,p
c
      adev = 0.0
      var = 0.0
      skew = 0.0
      curt = 0.0
      do i=1,n
        s=map(i)-ave
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
      end do
c
      adev=adev/xn
      var=var/float(n-1)
      p=sqrt(var)
      write (*,30) adev,var,p
c
      if (var .ne. 0.0) then
        skew=skew/(xn*p*p*p)
        curt = curt/(xn*var*var) - 3.0
c
        p=sqrt(6.0/xn)
        s=skew/p
        if (s .ge. 3.0) then
          line = 'Asymmetric tail extending towards positive density'
        else if (s .le. -3.0) then
          line = 'Asymmetric tail extending towards negative density'
        else
          line = 'Skewness not significant; symmetric distribution'
        end if
        write (*,40) skew,p,s,line(1:leng1(line))
c
        p=sqrt(24.0/xn)
        s=curt/p
        if (s .ge. 3.0) then
          line = 'Leptokurtic (sharp) distribution'
        else if (s .le. -3.0) then
          line = 'Platykurtic (flat) distribution'
        else
          line = 'Mesokurtic (normal-like) distribution'
        end if
        write (*,50) curt,p,s,line(1:leng1(line))
c
      else
        call errcon ('Zero variance; cannot compute higher moments')
      end if
c
      write (*,*)
c
   10 format (/' Extensive map statistics'/
     +  ' Number of data points ................... ',i15)
   20 format (
     +  ' Average electron density ................ ',1pe15.5/
     +  ' RMS density ("sigma") ................... ',e15.5)
   30 format (
     +  ' Mean Absolute Deviation from the mean ... ',1pe15.5/
     +  ' Variance ................................ ',e15.5/
     +  ' Standard deviation ...................... ',e15.5)
   40 format (
     +  ' Skewness (third moment) ................. ',1pe15.5/
     +  ' Sigma(skewness) estimate ................ ',e15.5/
     +  ' Skewness / Sigma(skewness) .............. ',e15.5/
     +  ' ',a)
   50 format (
     +  ' Kurtosis (fourth moment) ................ ',1pe15.5/
     +  ' Sigma(kurtosis) estimate ................ ',e15.5/
     +  ' Kurtosis / Sigma(kurtosis) .............. ',e15.5/
     +  ' ',a)
c
      return
      end
c
c
c
      subroutine setpkm (buf,map,na,nb,nc,pklim,pklev,ns,icnt)
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,pklim(2,3),icnt
      real map(na,nb,nc),buf(na,nb,nc),ns,pklev
c
code ...
c
      icnt = 1
      do i1=1,3
        icnt = icnt * (pklim(2,i1)-pklim(1,i1)+1)
      end do
      call jvalut (' Nr of points in pick volume :',1,icnt)
      if (icnt .lt. 1) return
c
      icnt = 0
      do i1=pklim(1,1),pklim(2,1)
        do i2=pklim(1,2),pklim(2,2)
          do i3=pklim(1,3),pklim(2,3)
            if (map (i1,i2,i3) .ge. pklev) then
              buf (i1,i2,i3) = ns
              icnt = icnt + 1
            end if
          end do
        end do
      end do
c
      call jvalut (' Nr of points >= pick level  :',1,icnt)
c
      return
      end
c
c
c
      subroutine pick3d (buf,map,na,nb,nc,
     +                   pklim,pklev,np,maxpk,peaks,ihow)
c
c ... ihow = 0 -> interpolated intensity
c     ihow = 1 -> sum of positive values in 2*2*2 box
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,pklim(2,3),np,j1,j2,j3
      integer maxpk,ii,ihow
c
      real map(na,nb,nc),buf(na,nb,nc),pklev,peaks(4,maxpk)
      real xdum,values(27),xnew(3),xval,sum
c
      logical init
c
      data init /.false./
c
      save init
c
code ...
c
ccc      init = .false.
c
cc      print *,' IHOW = ',ihow
c
      np = 0
      do i1=pklim(1,1),pklim(2,1)
        do i2=pklim(1,2),pklim(2,2)
          do i3=pklim(1,3),pklim(2,3)
            if (buf (i1,i2,i3) .gt. 0.0) then
              xdum = map(i1,i2,i3)
c
c ... if only high values need to be picked, skip maximum test
c
              if (ihow .eq. 2) goto 12
c
              if (xdum .lt. map(i1-1,i2-1,i3-1)) goto 10
              if (xdum .lt. map(i1-1,i2-1,i3  )) goto 10
              if (xdum .lt. map(i1-1,i2-1,i3+1)) goto 10
              if (xdum .lt. map(i1-1,i2  ,i3-1)) goto 10
              if (xdum .lt. map(i1-1,i2  ,i3  )) goto 10
              if (xdum .lt. map(i1-1,i2  ,i3+1)) goto 10
              if (xdum .lt. map(i1-1,i2+1,i3-1)) goto 10
              if (xdum .lt. map(i1-1,i2+1,i3  )) goto 10
              if (xdum .lt. map(i1-1,i2+1,i3+1)) goto 10
c
              if (xdum .lt. map(i1  ,i2-1,i3-1)) goto 10
              if (xdum .lt. map(i1  ,i2-1,i3  )) goto 10
              if (xdum .lt. map(i1  ,i2-1,i3+1)) goto 10
              if (xdum .lt. map(i1  ,i2  ,i3-1)) goto 10
              if (xdum .lt. map(i1  ,i2  ,i3+1)) goto 10
              if (xdum .lt. map(i1  ,i2+1,i3-1)) goto 10
              if (xdum .lt. map(i1  ,i2+1,i3  )) goto 10
              if (xdum .lt. map(i1  ,i2+1,i3+1)) goto 10
c
              if (xdum .lt. map(i1+1,i2-1,i3-1)) goto 10
              if (xdum .lt. map(i1+1,i2-1,i3  )) goto 10
              if (xdum .lt. map(i1+1,i2-1,i3+1)) goto 10
              if (xdum .lt. map(i1+1,i2  ,i3-1)) goto 10
              if (xdum .lt. map(i1+1,i2  ,i3  )) goto 10
              if (xdum .lt. map(i1+1,i2  ,i3+1)) goto 10
              if (xdum .lt. map(i1+1,i2+1,i3-1)) goto 10
              if (xdum .lt. map(i1+1,i2+1,i3  )) goto 10
              if (xdum .lt. map(i1+1,i2+1,i3+1)) goto 10
c
c ... okay, we found a local maximum
c
   12         continue
              if (np .ge. maxpk) then
                call errcon ('Too many peaks')
                goto 20
              end if
c
              np = np + 1
              peaks (1,np) = float(i1)
              peaks (2,np) = float(i2)
              peaks (3,np) = float(i3)
              peaks (4,np) = xdum
c
c ... interpolation
c
              if (ihow .eq. 2) goto 14
c
              ii = 0
              do j1=i1-1,i1+1
                do j2=i2-1,i2+1
                  do j3=i3-1,i3+1
                    ii = ii + 1
                    values (ii) = map(j1,j2,j3)
                  end do
                end do
              end do
c
cc              call fvalut (' Before :',4,peaks(1,np))
c
              call svdmx3 (values,peaks(1,np),peaks(4,np),
     +                     xnew,xval,init)
c
C     SUBROUTINE SVDMX3 (F,O,OEXTR,S,FEXTR,INIT)
C     ================= get precise 3D peak position through interpolation
C
C F (27)      function values of the 27 points
C O (3)       coordinates of the local maximum
C OEXTR       function value at the local maximum
C S (3)       coordinates of the interpolated maximum
C FEXTR       function value at this point
C INIT        logical, .FALSE. on first call, .TRUE. thereafter
C
c
              peaks (1,np) = xnew(1)
              peaks (2,np) = xnew(2)
              peaks (3,np) = xnew(3)
              peaks (4,np) = xval
c
cc              call fvalut (' After  :',4,peaks(1,np))
c
c ... zap block around this peak (& "integrate")
c
              sum = 0.0
              do j1=i1-1,i1+1
                do j2=i2-1,i2+1
                  do j3=i3-1,i3+1
                    if (map (j1,j2,j3) .gt. 0.0) then
                      sum = sum + map (j1,j2,j3)
                    end if
                    buf (j1,j2,j3) = -1.0
                  end do
                end do
              end do
c
c ... if ihow = 1 -> use "integral"
c
              if (ihow .eq. 1) peaks (4,np) = sum
c
   14         continue
c
c ... done with this peak
c
   10         continue
            end if
          end do
        end do
      end do
c
   20 continue
      call jvalut (' Nr of peaks found :',1,np)
c
      return
      end
c
c
c
	subroutine omapin (file,lun,ed,orgn,ext,grid,cell,space,
     +                     krecl,lbrix,ierr)
c
c ... OMAPIN - read DSN6 or BRIX style O maps
c
c ... file = char = name of input file
c     lun  = int  = unit nr
c     ed   = real = array to store density
c     orgn = int  = origin of map
c     ext  = int  = extent of map
c     grid = int  = grid of map
c     cell = real = cell constants
c     space = int = max size of ED array
c     krecl = int = size of records (must be 128)
c     lbrix = log = if true, read BRIX; if false, read DNS6
c     ierr  = int = err flag; if =0 then all is well
c
	implicit none
c
	character file*(*)
	integer lun, orgn(3), ext(3), grid(3), space, krecl, ierr
	real ed(space), cell(6)
c
	integer ct, i, j, k, l, m, n, pt, pt1, pt2, pt3, plus
	integer iindex, iv, inx, iny, inz, errcod, idum
	integer*2 header(256)
	byte box(8,8,8), v(4)
	real prod,sigma
        logical little,litend,lbrix
        character str*512
	equivalence (iv, v(1))
c
code ...
c
        ierr = -1
c
        little = litend()
c
c ... set all four bytes of IV to zero
c
        iv = 0
c
c	type*,'hi'
c      print *,'Open file with record length ',krecl
c        close (lun)
	open (unit=lun, file=file, status='old', err=999,
     $	      recl=krecl, access='direct')
c      print *,'Opened okay'
c	type*,' 0'
c        type*,'i/o error = ',errcod
c
        if (lbrix) then
c
          if (little) then
            call prompt (' Little-endian machine')
          else
            call prompt (' Big-endian machine')
          end if
c
          read (lun, rec=1, iostat=ierr) str
          if (ierr .ne. 0) goto 9999
          if (str(1:3) .ne. ':-)') then
            call errcon ('Not a BRIX file')
            goto 9999
          end if
c
          call textut (' Header :',str)
c
c          read (str, 10) (orgn(i),i=1,3), 
c     +      (ext(i),i=1,3), (grid(i),i=1,3),
c     +      (cell(i),i=1,6), prod, plus, sigma
c
c  10      format (10x,3i5,7x,3i5,5x,3i5,6x,6f10.3,
c     +            5x,f12.5,5x,i8,7x,f12.5)
c
          call upcase (str)
c
          j = index (str,'ORIGIN')
          if (j .le. 0) then
            call errcon ('Could not find ORIGIN in header')
            goto 9999
          end if
          read (str(j+6:),*) (orgn(i),i=1,3)
c
          j = index (str,'EXTENT')
          if (j .le. 0) then
            call errcon ('Could not find EXTENT in header')
            goto 9999
          end if
          read (str(j+6:),*) (ext(i),i=1,3)
c
          j = index (str,'GRID')
          if (j .le. 0) then
            call errcon ('Could not find GRID in header')
            goto 9999
          end if
          read (str(j+4:),*) (grid(i),i=1,3)
c
          j = index (str,'CELL')
          if (j .le. 0) then
            call errcon ('Could not find CELL in header')
            goto 9999
          end if
          read (str(j+4:),*) (cell(i),i=1,6)
c
          j = index (str,'PROD')
          if (j .le. 0) then
            call errcon ('Could not find PROD in header')
            goto 9999
          end if
          read (str(j+4:),*) prod
c
          j = index (str,'PLUS')
          if (j .le. 0) then
            call errcon ('Could not find PLUS in header')
            goto 9999
          end if
          read (str(j+4:),*) plus
c
          j = index (str,'SIGMA')
          if (j .le. 0) then
            call errcon ('Could not find SIGMA in header')
            goto 9999
          end if
          read (str(j+5:),*) sigma
c
          call prompt (' Header read OK')
          call rvalut (' Prod  :',1,prod)
          call jvalut (' Plus  :',1,plus)
          call rvalut (' Sigma :',1,sigma)
c
c 10     format (":-) Origin", 3i5," Extent", 3i5, " Grid", 3i5,
c     $        " Cell ", 6f10.3, " Prod", f12.5, " Plus",i8, 
c     $        " Sigma ", f12.5)
c
        else
    	    call maprdr (lun, 1, header, 256, lbrix, errcod)
          if (errcod .ne. 0) then
            call prompt (' Error while reading header record')
            call ivalut (' IOSTAT code          :',1,errcod)
            goto 9999
          end if
c
          if (little) then
            call prompt (' Little-endian machine; will swap header')
cccc          call swpbyt (header,256)
            call bytswp (header)
          else
            call prompt (' Big-endian machine; will NOT swap header')
          end if
c
c        type *,(header(i),i=1,20)
  	  j = header(18)
	  k = header(19)
	  if (j .eq. 0) j = 100
	  if (k .eq. 0) k = 100
c	type*,' 2'
	  do 100 i=1,3
	    orgn(i) = header(i)
	    ext(i) = header(i+3)
	    grid(i) = header(i+6)
	    cell(i) = float(header(i+9))/float(j)
100	    cell(i+3) = float(header(i+12))/float(j)
c
c        type *,orgn
c        type *,ext
c        type *,grid
c        type *,cell
c	type*,' 3'
c
	  prod = float(header(16))/float(k)
	  plus = header(17)
          call rvalut (' Prod :',1,prod)
          call jvalut (' Plus :',1,plus)
        end if
c
c        call ivalut (' Cell grid   :',3,grid)
c        call ivalut (' Cell origin :',3,orgn)
c        call ivalut (' Cell extent :',3,ext)
c        call fvalut (' Cell axes   :',3,cell(1))
c        call fvalut (' Cell angles :',3,cell(4))
c
	inx = ext(1)/8
	iny = ext(2)/8
	inz = ext(3)/8
	if (mod(ext(1),8) .ge. 1) inx = inx+1
	if (mod(ext(2),8) .ge. 1) iny = iny+1
	if (mod(ext(3),8) .ge. 1) inz = inz+1
        i = inx * iny * inz
        call jvalut (' Nr of pages :',1,i)
c	type*,inx, iny, inz
	ct = 1
c	type*,'ext', ext
c
        idum = ext(1)*ext(2)*ext(3)
        call jvalut (' Map size :',1,idum)
        if (idum .gt. space) then
          call errcon ('Map too big')
          call jvalut (' Available memory :',1,space)
          call jvalut (' Required  memory :',1,idum)
          return
        end if
c
        call prompt (' Reading map ...')
c
	do 200 k=1,inz
	do 200 j=1,iny
	do 200 i=1,inx
	  pt = ct
	  ct = ct+1
c
	  call maprdr (lun, ct, box, 256, lbrix, errcod)
          if (errcod .ne. 0) then
            call ivalut (' Read error at record :',1,ct)
            call ivalut (' IOSTAT code          :',1,errcod)
            goto 9999
          end if
c
c ... swap bytes ???
c
        if (lbrix) then
c
c ... 980922 - do NOT swap bytes of BRIX file
c
ccc          if (little) call bytswp (box)
c
        else
          call bytswp (box)
        end if
c
	  do 300 n=1,8
	  do 300 m=1,8
	  do 300 l=1,8
	    pt3 = (k-1)*8+ n
	    pt2 = (j-1)*8+ m
	    pt1 = (i-1)*8+ l
	    if (pt1 .gt. ext(1)) goto 300
	    if (pt2 .gt. ext(2)) goto 300
	    if (pt3 .gt. ext(3)) goto 300
	    iindex= (pt3-1)*ext(1)*ext(2)+ (pt2-1)*ext(1)+ pt1
c
c	type99,i,j,k,l,m,n,pt1,pt2,pt3,iindex
c99	format (10i6)
c
c ... the opposite of MAKBYT
c     the following does *NOT* work : iv = box(l,m,n)
c
c            if (lbrix) then
c              ed (iindex) = float ( box(l,m,n) - plus ) / prod
c            else
c
c ... set proper byte, depending on ended-ness
c
               if (little) then
                 v(1) = box (l,m,n)
               else
                 v(4) = box (l,m,n)
               end if
c
ccc            if (iv .lt. 0 .or. iv .gt. 255) then
ccc              print *,' OOPS - ',iv,' ~ ',box(l,m,n)
ccc            end if
c
	         ed(iindex) = float(iv-plus)/prod
c            end if
c
300	  continue
200	continue
c
        ierr = 0
ccc        call jvalut (' Nr of records read :',1,ct)
        close (lun)
        return
c
  999   continue
        call errcon ('While opening file')
        ierr = -2
        close (lun)
        return
c
 9999   continue
        call errcon ('While reading file')
        ierr = -3
        close (lun)
        return
c
	end
c
c
c
       subroutine maprdr(lun, ipt, aray, mysize, lbrix, errcod)
c
c ---   Read a random access record of 'size' integer*2 words into 'aray'
c       from 'lun' at record 'ipt'
c
        implicit none
c
        integer size, lun, ipt, errcod, mysize
        integer*2 aray(mysize)
        integer i
c
        logical lbrix
c
c        character*1 record(512)
c
code ...
c
        size = mysize
c
c        if (lbrix) then
c          read (lun,rec=ipt,iostat=errcod) (record(i),i=1,size)
c          do i=1,size
c            aray(i) = ichar (record (i))
c          end do
c        else
c        print *,'About to read record'
  100   continue
c        print *,'size = ',size
c   	  read (lun,rec=ipt,iostat=errcod) aray
   	  read (lun,rec=ipt,iostat=errcod) (aray(i),i=1,size)
c        if (errcod .ne. 0) then
c          size = size - 1
c          goto 100
c        end if
c        print *,'Read header with size=',size
c        print *,'Done reading record; iostat =',errcod
c        end if
c
        return
        end
c
c
c
       subroutine swpbyt (a,index)
c
c ---   swap the bytes within an integer*2 array
c
        implicit none
c
        integer*2 a(*)
        integer index
        integer i
        integer*2 ivalue
        byte b, value(2)
        equivalence (ivalue,value(1))
c
code ...
c
        do 100 i=1,index
          ivalue = a(i)
          b = value(1)
          value(1) = value(2)
          value(2) = b
100       a(i) = ivalue
c
        return
        end
c
c
c
      subroutine envset (env,ed,e1,e2,e3,step,base,max,ierr)
c
      implicit none
c
      integer e1, e2, e3, max,ierr
      real ed(e1, e2, e3), step, base, qqq
      integer env(e1, e2, e3)
c
      integer i, j, k, l
c
code ...
c
      max = -9999999
      qqq = 1.0 / step
c
      do 100 k=1,e3
        do 100 j=1,e2
          do 100 i=1,e1
            l = int ((ed(i,j,k)-base)*qqq)
            if (ed(i,j,k) .ge. base) l = l+1
            if (l .le. 0) l = 0
            if (l .gt. max) max= l
            env(i,j,k) = l
100   continue
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine edchck(envlop,nx,ny,nz,npnt,ierr)
c
      implicit none
c
      integer nx,ny,nz,npnt,ierr
      integer envlop(nx,ny,nz)
      integer ix,iy,iz
c
code ...
c
      do iy = 1,ny
        do iz = 1,nz
          envlop(1,iy,iz) = 0
          envlop(nx,iy,iz) = 0
        enddo
      enddo
c
      do ix = 1,nx
        do iz = 1,nz
          envlop(ix,1,iz) = 0
          envlop(ix,ny,iz) = 0
        enddo
      enddo
c
      do ix = 1,nx
        do iy = 1,ny
          envlop(ix,iy,1) = 0
          envlop(ix,iy,nz) = 0
        enddo
      enddo
c
      npnt = 0
      do ix = 1,nx
        do iy = 1,ny
          do iz = 1,nz
            if(envlop(ix,iy,iz).gt.0)npnt = npnt+1
          enddo
        enddo
      enddo
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine sklton(envlop,nx,ny,nz,ierr)
c
c ---      given the density array,make a skeleton of it.
c
      implicit none
c
      integer nx,ny,nz,ierr
      integer envlop(nx,ny,nz),naybox(3,3,3)
      integer nxm1,nym1,nzm1,lev,i,j,k,ix,iy,iz,jx,jy,jz
      integer ict,nx1,nx2,ny1,ny2,nz1,nz2,nx3,nx4,ny3,ny4
      integer nz3,nz4,kx,ky,kz,lx,ly,lz
c
      logical flag
c
code ...
c
c --- remove all points below current value 
c      do k=1,21
c        write(1,15) k
c        write(1,10) ((envlop(i,j,k),j=1,11),i=1,16)
c      enddo
c10    format( 16(/11i4))
c15    format( //' this is the masked density',i4,//)
c
      nxm1 = nx-1
      nym1 = ny-1
      nzm1 = nz-1
      do 1000 lev=1,11
        do 100 i=2,nzm1
        do 100 j=2,nym1
        do 100 k=2,nxm1
          if(envlop(k,j,i) .le. 0)envlop(k,j,i) = 0
          if(lev .eq. 11)then
c ---   in the 11'th pass all levels are checked
            if(envlop(k,j,i) .eq. 0)goto 100
          else
c ---      greer's condition 1,if not equal current level then keep it
            if(envlop(k,j,i) .ne. lev)goto 100
          endif
c ---      greer's condition 2
          flag = .false.
          if(envlop(k,j-1,i) .eq. 0)flag = .true.
          if(envlop(k,j+1,i) .eq. 0)flag = .true.
          if(envlop(k-1,j,i) .eq. 0)flag = .true.
          if(envlop(k+1,j,i) .eq. 0)flag = .true.
          if( .not. flag) goto 100
c
          flag = .false.
          if(envlop(k,j,i-1) .eq. 0)flag = .true.
          if(envlop(k,j,i+1) .eq. 0)flag = .true.
          if(envlop(k,j-1,i) .eq. 0)flag = .true.
          if(envlop(k,j+1,i) .eq. 0)flag = .true.
          if( .not. flag) goto 100
c
          flag = .false.
          if(envlop(k,j,i-1) .eq. 0)flag = .true.
          if(envlop(k,j,i+1) .eq. 0)flag = .true.
          if(envlop(k-1,j,i) .eq. 0)flag = .true.
          if(envlop(k+1,j,i) .eq. 0)flag = .true.
          if( .not. flag) goto 100
c ---      greer's condition 3
          ict = 0
          nx1 = k-1
          nx2 = k+1
          ny1 = j-1
          ny2 = j+1
          nz1 = i-1
          nz2 = i+1
          do 110 iz=nz1,nz2
          do 110 iy=ny1,ny2
          do 110 ix=nx1,nx2
            if(envlop(ix,iy,iz) .gt. 0)ict = ict+1
  110     continue
c ---      if it is a tip (ict .le. 2) we cannot delete it
          if(ict .le. 2) goto 100
c ---      greer's condition 4 , does it maintain continuity ?
          do 120 iz=1,3
          do 120 iy=1,3
          do 120 ix=1,3
  120       naybox(ix,iy,iz) = 0
c ---      search the 27 naybours for a hit
          do 130 iz=nz1,nz2
          do 130 iy=ny1,ny2
          do 130 ix=nx1,nx2
            if(ix.eq.k.and.iy.eq.j.and.iz.eq.i)goto 130
            if(envlop(ix,iy,iz) .gt. 0)goto 140
130       continue
          goto 100
c ---      mark nearest neighbours to first point in 3*3*3 list
140       call naydx(ix,iy,iz,nx,ny,nz,nx3,nx4,ny3,ny4,nz3,nz4)
          do 143 jz=nz3,nz4
          do 143 jy=ny3,ny4
          do 143 jx=nx3,nx4
            if(jx.eq.k.and.jy.eq.j.and.jz.eq.i)goto 143
            if(envlop(jx,jy,jz) .gt. 0)then
              call markit(jx,jy,jz,k,j,i,naybox,ict)
            end if
143       continue
c ---      go over 3*3*3,check if marked,if it is mark it's neighbours.
c ---      ict is only updated if a new mark is made.if a pass has no new
c ---      marks then we have finished.
150       ict = 0
          do 160 iz=nz1,nz2
          do 160 iy=ny1,ny2
          do 160 ix=nx1,nx2
            if(naybox(ix-nx1+1,iy-ny1+1,iz-nz1+1) .ne. 0)then
              call naydx(ix,iy,iz,nx,ny,nz,nx3,nx4,ny3,ny4,nz3,nz4)
              do 163 jz=nz3,nz4
              do 163 jy=ny3,ny4
              do 163 jx=nx3,nx4
                if(jx.eq.k.and.jy.eq.j.and.jz.eq.i)goto 163
                if(envlop(jx,jy,jz) .gt. 0)then
                  call markit(jx,jy,jz,k,j,i,naybox,ict)
                end if
163           continue
            end if
160       continue
          if(ict .gt. 0)goto 150
c ---      no new neighbours.
c ---      all marked,now see if any data points not in naybox
          do 170 kz=nz1,nz2
          do 170 ky=ny1,ny2
          do 170 kx=nx1,nx2
            if (kx .eq. k .and.
     $          ky .eq. j .and.
     $          kz .eq. i ) goto 170
            if(envlop(kx,ky,kz) .le. 0) goto 170
c ---      are there any points which are not also in the box ?
            lx = kx-nx1+1
            ly = ky-ny1+1
            lz = kz-nz1+1
            if(naybox(lx,ly,lz) .eq. 1)goto 170
c ---      if we come through,then the point is needed to conserve connectivity
c      and cannot be removed
            goto 100
  170    continue
c ---      all satisified,so transfer it to subset r
         envlop(k,j,i) = -envlop(k,j,i)
         goto 100
  100  continue
c ---      end of this pass, so tranfer subset r to subset n
       do 200 i=2,nzm1
       do 200 j=2,nym1
       do 200 k=2,nxm1
         if(envlop(k,j,i).lt.0)envlop(k,j,i) = 0
 200   continue
 1000  continue
c
c      do k=1,21
c            write(1,152) k
c            write(1,102) ((envlop(i,j,k),j=1,11),i=1,16)
c      enddo
c102      format( 16(/11i4))
c152      format( //' this is the skeletonized density',i4,//)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine naydx(ix,iy,iz,nx,ny,nz,nx3,nx4,ny3,ny4,nz3,nz4)
c
      implicit none
c
      integer ix,iy,iz,nx,ny,nz,nx3,nx4,ny3,ny4,nz3,nz4
c
code ...
c
      nz3 = iz-1
      nz4 = iz+1
      if(nz3 .lt. 1)nz3 = 1
      if(nz4 .gt.nz)nz4 = nz
      ny3 = iy-1
      ny4 = iy+1
      if(ny3 .lt. 1)ny3 = 1
      if(ny4 .gt.ny)ny4 = ny
      nx3 = ix-1
      nx4 = ix+1
      if(nx3 .lt. 1)nx3 = 1
      if(nx4 .gt.nx)nx4 = nx
c
      return
      end
c
c
c

      subroutine markit(ix,iy,iz,k,j,i,naybox,ict)
c
      implicit none
c
      integer naybox(3,3,3),ix,iy,iz,k,j,i,ict,i1,i2,i3
c
code ...
c
      i1 = ix-k+2
      i2 = iy-j+2
      i3 = iz-i+2
      if(i1 .le. 0)return
      if(i1 .gt. 3)return
      if(i2 .le. 0)return
      if(i2 .gt. 3)return
      if(i3 .le. 0)return
      if(i3 .gt. 3)return
      if(naybox(i1,i2,i3) .gt. 0)return
      naybox(i1,i2,i3) = 1
      ict = ict+1
c
      return
      end
c
c
c
      subroutine chncon (envlop,nx,ny,nz,ipext,cell,grid,origin,
     +  maxbop,maxbob,maxboc,pntlis,a,levlis,brnlis,
     +  conlis,ierr)
c
c      find a point in the bones and use it as the starting point of a
c      chain. every branch point is used as starting points for new chains.
c
      implicit none
c
      integer nx,ny,nz,ipext,ierr
      integer envlop(nx,ny,nz)
c
      real cell(6)
      integer grid(3),origin(3)
c
ccc      common/page/orgxyz(3),extxyz(3),uvw(3),grid(3),cell(6)
ccc      integer orgxyz,extxyz,uvw,grid
c
      integer maxbop,maxbob,maxboc
      real pntlis(maxbop,3),a(3,3)
      integer levlis(maxbop),brnlis(maxbob,6),conlis(3,maxboc)
c
ccc      common/point/pntlis(16000,3),a(3,3),levlis(16000)
ccc      common/lists/conlis(3,20000),brnlis(10000,6)
ccc      integer conlis,brnlis,levlis
c
      real unknx(3),dx
c
      integer npcn,nccn,nlab,nbran,irand
      integer pcount,bcount,ccount,fcount,i,j
      integer finlis,even,odd,unkn,side,main,unus
      integer lno,ix,iy,iz,maxp,maxb,maxc
      integer maxx,maxy,maxz,iscnx,iscny,iscnz,index,ival,iblng
      integer ixm1,ixp1,iym1,iyp1,izm1,izp1,nflg,lev,jx,jy,jz
      integer lx,ly,lz,lim1,inlow,iext,maxlen
      integer frgsz(15),nfrg(16)
c
c      integer icoo(3)
c      equivalence(ix,icoo(1)),(iy,icoo(2)),(iz,icoo(3))
c
      data finlis/3/,even/0/,odd/1/,unkn/1/,side/2/,main/3/,unus/4/
      data unknx/1500.,1500.,1500./
      data frgsz/1,2,3,4,5,6,7,8,9,10,15,20,25,50,100/
c
code ...
c
      ierr = -1
c
      maxp = maxbop
      maxb = maxbob
      maxc = maxboc
      maxlen = 0
c
c      set the pointers of the lists :
c
      nbran = 0
      pcount=0
      bcount=0
      ccount=0
      fcount=0
      lno=1
      do i = 1,16
        nfrg(i) = 0
      enddo
c
c      get matrix for converting grid coord. to orthog. angstroms :
c
      call orthog (cell,a,0)
      do i=1,3
      do j=1,3
        a(i,j)=a(i,j)/(float(grid(j)))
      enddo
      enddo
c
c      scan the envelope :
c
      maxx=nx-1
      maxy=ny-1
      maxz=nz-1
c    points of the edges are not included in the
c    scan but if there exists points on the edges
c    with a value > 0 these will be detected and
c    connected during the neighbour search.
      do 100 iscnz=2,maxz
      do 100 iscny=2,maxy
      do 100 iscnx=2,maxx
        ix=iscnx
        iy=iscny
        iz=iscnz      
        if(envlop(ix,iy,iz).eq.0)goto 100
c            -- point in the bones found
c            -- store it in the point list :
        pcount=pcount+1
        if(pcount.gt.maxp) then
          call errcon ('CHNCON ** 1 ** too many points/connections')
          call jvalut (' Required  :',1,pcount)
          call jvalut (' Allocated :',1,maxp)
          goto 9999
        end if
c
        index=pcount
        ival=envlop(ix,iy,iz)
        envlop(ix,iy,iz)=0
        call points(index,ix,iy,iz,ival,maxbop,pntlis,levlis,a,origin)
        iblng=0
120     continue
c            -- and store a "move to" in the connectivity list :
        ccount=ccount+1
        if(ccount.gt.maxc) then
          call errcon ('CHNCON ** 2 ** too many points/connections')
          call jvalut (' Required  :',1,ccount)
          call jvalut (' Allocated :',1,maxc)
          goto 9999
        end if
c
        conlis(lno,ccount)=index*2
110     continue
c            scan the neighbours :
        ixm1=ix-1
        ixp1=ix+1
        iym1=iy-1
        iyp1=iy+1
        izm1=iz-1
        izp1=iz+1
        nflg=0
        lev = 0
c
        do 151 jx=ixm1,ixp1
        do 152 jy=iym1,iyp1
        do 153 jz=izm1,izp1
          if(jx.lt.1 .or. jx.gt.nx)goto 151
          if(jy.lt.1 .or. jy.gt.ny)goto 152
          if(jz.lt.1 .or. jz.gt.nz)goto 153
          if(envlop(jx,jy,jz).ne.0)then
c                  found a neighbour :
            nflg=nflg+1
            if(envlop(jx,jy,jz).gt.lev)then
              lev = envlop(jx,jy,jz)
c                        the first neighbour to be used:
              lx=jx
              ly=jy
              lz=jz
            endif
            if(nflg.eq.2)then
c                        a second neighbour i.e. it's a branchpoint
c                        -- store in branch list :
              bcount = bcount+1
              if(bcount.gt.maxb) then
                call errcon ('CHNCON ** 3 ** too many branches')
          call jvalut (' Required  :',1,bcount)
          call jvalut (' Allocated :',1,maxb)
                goto 9999
              end if
c
              nbran = nbran+1
              brnlis(bcount,1) = index
              brnlis(bcount,2) = ix
              brnlis(bcount,3) = iy
              brnlis(bcount,4) = iz
              brnlis(bcount,5) = 1
              brnlis(bcount,6) = iblng
            endif
          endif
153     continue
152     continue
151     continue
c
125     if(nflg.ne.0)then
c            continue the chain
c            -- store in point list 
          ix=lx
          iy=ly
          iz=lz
          pcount=pcount+1
          if(pcount.gt.maxp) then
            call errcon ('CHNCON ** 4 ** too many points/connections')
          call jvalut (' Required  :',1,pcount)
          call jvalut (' Allocated :',1,maxp)
            goto 9999
          end if
c
          index=pcount
          ival=envlop(ix,iy,iz)
          envlop(ix,iy,iz)=0
          call points(index,ix,iy,iz,ival,maxbop,pntlis,levlis,a,origin)
c            -- store a "line to" in connectivity list
          ccount=ccount+1
          if(ccount.gt.maxc) then
            call errcon ('CHNCON ** 5 ** too many points/connections')
          call jvalut (' Required  :',1,ccount)
          call jvalut (' Allocated :',1,maxc)
            goto 9999
          end if
c
          conlis(lno,ccount)=index*2+1
          iblng=iblng+1
          goto 110
       else
c            tip is reached
c            -- go back to the latest branchpoint
130      continue
         if(bcount.eq.0)goto 140
         index=brnlis(bcount,1)
         ix=brnlis(bcount,2)
         iy=brnlis(bcount,3)
         iz=brnlis(bcount,4)
         iblng=brnlis(bcount,6)
c                  -- check if there is still a neighbour
         ixm1=ix-1
         ixp1=ix+1
         iym1=iy-1
         iyp1=iy+1
         izm1=iz-1
         izp1=iz+1
c
          do 161 jx=ixm1,ixp1
          do 162 jy=iym1,iyp1
          do 163 jz=izm1,izp1
            if(jx.lt.1 .or. jx.gt.nx)goto 161
            if(jy.lt.1 .or. jy.gt.ny)goto 162
            if(jz.lt.1 .or. jz.gt.nz)goto 163
            if(envlop(jx,jy,jz).ne.0) then
              lx=jx
              ly=jy
              lz=jz
              brnlis(bcount,5)=brnlis(bcount,5)+1
              nflg=1
c                              -- store a "move to" in list :
              ccount=ccount+1
              if(ccount.gt.maxc) then
                call errcon (
     +            'CHNCON ** 6 ** too many points/connections')
          call jvalut (' Required  :',1,ccount)
          call jvalut (' Allocated :',1,maxc)
                goto 9999
              end if
c
              conlis(lno,ccount)=index*2
              goto 125
            endif
163       continue
162       continue
161       continue
c                  no more branches at this point
c                  
c                  -- so take the longest branch as the 'mainchain'
          call branch(lno,ccount,bcount,fcount,maxbop,maxbob,
     +      maxboc,brnlis,conlis,ierr)
          if (ierr .ne. 0) return
          goto 130
        endif
c      --finished a whole chain 
c      --transfer it to the final list
140     if(ccount.ge.1)then
          fcount = fcount+1
          if(fcount.gt.maxc) then
            call errcon (
     +        'CHNCON ** 7 ** too many points/connections')
          call jvalut (' Required  :',1,fcount)
          call jvalut (' Allocated :',1,maxc)
            goto 9999
          end if
c
          lim1 = 1
          call movlis(lno,finlis,fcount,lim1,ccount,even,
     +      maxboc,conlis,ierr)
          if(fcount.gt.maxc) then
            call errcon (
     +        'CHNCON ** 8 ** too many points/connections')
          call jvalut (' Required  :',1,fcount)
          call jvalut (' Allocated :',1,maxc)
            goto 9999
          end if
c
          do i = 1,15
            if(ccount.le.frgsz(i))then
              nfrg(i) = nfrg(i)+1
              goto 200
            endif
          enddo
200       if(ccount.gt.frgsz(15))nfrg(16) = nfrg(16)+1
          maxlen = max (maxlen, ccount)
          ccount = 0
        endif
100   continue
c      the connectivity is now defined
c      redefine the levels
c      all points are given the level 'unknown'
      do 400 i = 1,pcount
        levlis(i) = unkn
400   continue
c      insert extra points
      inlow = pcount+1
      pcount = pcount+ipext
      iext = ipext
      if(pcount.gt.maxp)then
        write(*,1200)
        pcount = maxp
        iext = maxp-pcount+1
      endif
c
c ... instead of storing the spares at (1500,1500,1500)
c     store them near randomly selected other bones
c
c ... initialise random number generator with MCLOCK()
c
      call gkrand (dx,0.0,0.0,-1)
c
      do 300 i = inlow,pcount
c
c ... find a random bones atom
c
        call gkrand (dx,1.0,float(inlow-1),0)
        irand = max(1,min(inlow,nint(dx)))
c
c ... put it (1,1,1) A away from that atom
c
        do 310 j = 1,3
ccc          pntlis(i,j) = unknx(j)
          pntlis(i,j) = pntlis(irand,j) + 1.0
310     continue
        levlis(i) = unus
300   continue
c      write it out
      npcn = pcount
      nccn = fcount
      nlab = 0
      call envsav (npcn, nccn, nlab, finlis, conlis)
      ierr = 0
c
ccc      write(3) npcn,nccn,nlab
ccc      write(3) ((pntlis(index,i),i=1,3),index=1,npcn)
ccc      write(3) (levlis(index),index=1,npcn)
ccc      write(3) (conlis(finlis,i),i=1,nccn)
c
      write(*,1000) pcount,iext,fcount,nbran
      write(*,1100) (frgsz(i),nfrg(i),i=1,15),frgsz(15),
     +  nfrg(16),maxlen
c
 1000 format (/
     +  ' Nr of points in the bones : ',i8/
     +  '            (includes spare points !)'/
     +  ' Nr of spare points        : ',i8/
     +  ' Nr of connections         : ',i8/
     +  ' Nr of branch-points       : ',i8)
c
 1100 format (/
     +  ' Nr of fragments'/
     +  15(' Size less than or equal to ',i3,' = ',i8/),
     +    ' Size greater than          ',i3,' = ',i8/
     +    ' Maximum fragment size       = ',i8/)
c
 1200 format(/' ***** Not enough space for extra points !')
c
      return
c
c ... an error occured
c
 9999 continue
      ierr = -2
c
      return
      end                  
c
c
c
      subroutine branch(lno,ccount,bcount,fcount,
     +                  maxbop,maxbob,maxc,brnlis,conlis,ierr)
c
      implicit none
c
      integer maxbop,maxbob,maxc,ierr
      integer brnlis(maxbob,6),conlis(3,maxc)
c
ccc      common/lists/conlis(3,20000),brnlis(10000,6)
ccc      integer conlis,brnlis
c
      integer ccount,bcount,fcount,lim1,lim2,i,j,k
      integer icc,lng,ibeg,lno,finlis,even,odd,lsfrom,lsto
      integer ioflw,in,nbrn,bln,in1,inb,max1,max2,icon,mi1,mi2
      integer maxboc
c
      integer point(6)
c
      data finlis/3/,even/0/,odd/1/,ioflw/32000/
c
code ...
c
      ierr = 0
      maxboc = maxc
c
      in=brnlis(bcount,1)
c     index of the branchpoint
      nbrn=brnlis(bcount,5)
c     no. of branches at this point
      bln=brnlis(bcount,6)
c     distance to the beginning of the chain
c
      in = in*2
c     the 'move to' to search for
      in1 = in+1
c     the 'line to' to search for
c
      icc = ccount+1
      inb = 0
      lng = 0
      max1 = 0
      max2 = 0
      lsto = 1
      if(lno.eq.1)lsto = 2
      lsfrom = lno
      if(nbrn.le.1)then
        bcount = bcount-1
        return
      endif
c      -- find the two longest branches and their lengths
100   icc = icc-1
      if(icc.lt.0) then
        write(*,101) icc
        ierr = -1
        return
      endif
101   format(' ***** Error in branch-search; icc = ',i8)
c
      icon = conlis(lno,icc)
      if(icon.ne.in .and. icon.ne.in1) then
        lng = lng+1
        goto 100
      endif
c      -- found a branchpoint
      inb = inb+1
      if(lng.gt.max1) then
c            -- the branch is the longest so far
        if(max1.ne.0) then
          if(max2.ne.0) then
c                  -- transfer the shortest branch to the final list
            lim1 = mi2
            lim2 = mi2+max2
            call squeze(lsfrom,finlis,fcount,ccount,lim1,lim2,
     +        maxboc,conlis,ierr)
c
            if(fcount.gt.maxc) then
              call errcon (
     +          'BRANCH ** br1 ** too many points/connections')
          call jvalut (' Required  :',1,fcount)
          call jvalut (' Allocated :',1,maxc)
              ierr = -1
              return
            end if
c
            brnlis(bcount,5) = brnlis(bcount,5)-1
            if(lim2.eq.1) then
              mi1 = mi1-(max2+1)
            endif
          endif
c                  -- make the previously longest branch no. 2
          max2 = max1
          mi2 = mi1
       endif
c            -- and store lenght and index of the new no. 1
       max1 = lng
       mi1 = icc
      else
       if(lng.ge.max2) then
         if(max2.ne.0) then
c                  -- the branch is long enough to be second
c                  -- transfer the previous no. 2 to the final list
           lim1 = mi2
           lim2 = mi2+max2
           call squeze(lsfrom,finlis,fcount,ccount,lim1,lim2,
     +        maxboc,conlis,ierr)
c
           if(fcount.gt.maxc) then
              call errcon (
     +          'BRANCH ** br2 ** too many points/connections')
          call jvalut (' Required  :',1,fcount)
          call jvalut (' Allocated :',1,maxc)
              ierr = -1
              return
            end if
c
            brnlis(bcount,5) = brnlis(bcount,5)-1
            if(lim2.eq.1) then
              mi1 = mi1-(lng+1)
            endif
          endif
          max2 = lng
          mi2 = icc
        else
c            -- this one didn't make it -- transfer it to the final list
          lim1 = icc
          lim2 = icc+lng
          call squeze(lsfrom,finlis,fcount,ccount,lim1,lim2,
     +        maxboc,conlis,ierr)
c
          if(fcount.gt.maxc) then
            call errcon (
     +        'BRANCH ** br3 ** too many points/connections')
          call jvalut (' Required  :',1,fcount)
          call jvalut (' Allocated :',1,maxc)
            ierr = -1
            return
          end if
c
          brnlis(bcount,5) = brnlis(bcount,5)-1
          if(lim2.eq.1) then
            mi1 = mi1-(lng+1)
            mi2 = mi2-(lng+1)
          endif
        endif
      endif
c
      lng = 0
      if(inb.ne.nbrn) goto 100
c      -- choose between the two branches at this point and the connection
c         from the starting point
      if(bln.ge.max2) then
c            -- the shortest of the two branches didn't make it
c            -- transfer it to the final list
        lim1 = mi2
        lim2 = mi2+max2
        call squeze(lsfrom,finlis,fcount,ccount,lim1,lim2,
     +        maxboc,conlis,ierr)
c
        if(fcount.gt.maxc) then
          call errcon (
     +      'BRANCH ** br4 ** too many points/connections')
          call jvalut (' Required  :',1,fcount)
          call jvalut (' Allocated :',1,maxc)
          ierr = -1
          return
        end if
c
        bcount = bcount-1
        return
      else
        if(bln.eq.0) then
c                  -- this is a special case : the starting point happens
c                     to be the first branchpoint as well.
c                  -- so we have to invert one of the chains and connect
c                     the other to it.
          ibeg = 1
          lim1 = mi2
          lim2 = mi2+max2
          call invlis(lsfrom,lsto,ibeg,lim1,lim2,even,maxbop,maxbob,
     +      maxboc,brnlis,conlis,ierr)
c                  -- the last element in the new list is overwritten by
c                     the first element from the longest chain 
          lim1 = mi1
          lim2 = mi1+max1
          call movlis(lsfrom,lsto,ibeg,lim1,lim2,odd,
     +      maxboc,conlis,ierr)
          bcount = bcount-1
          ccount = ibeg
          if(ccount.gt.maxc) then
            call errcon (
     +        'BRANCH ** br5 ** too many points/connections')
          call jvalut (' Required  :',1,ccount)
          call jvalut (' Allocated :',1,maxc)
            ierr = -1
            return
          end if
c
          lno = lsto
          return
        else
c            -- the connection from the starting point is shortest at the 
c               the moment
c            -- invert the tree structure to investigate it further
          ibeg = 1
          lim1 = mi1
          lim2 = mi1+max1
          call invlis(lsfrom,lsto,ibeg,lim1,lim2,even,maxbop,maxbob,
     +      maxboc,brnlis,conlis,ierr)
c
          lim1 = mi2
          lim2 = mi2+max2
          call movlis(lsfrom,lsto,ibeg,lim1,lim2,odd,
     +      maxboc,conlis,ierr)
c
          ibeg = ibeg+1
          lim1 = 1
          lim2 = mi1
          if(mi1.gt.mi2)lim2 = mi2
          call invlis(lsfrom,lsto,ibeg,lim1,lim2,even,maxbop,maxbob,
     +      maxboc,brnlis,conlis,ierr)
        endif
c
        ccount = ibeg
        if(ccount.gt.maxc) then
          call errcon (
     +      'BRANCH ** br6 ** too many points/connections')
          call jvalut (' Required  :',1,ccount)
          call jvalut (' Allocated :',1,maxc)
          ierr = -1
          return
        end if
c
        lno = lsto
c            --and invert the branch-list as well
        j = 1
        k = bcount
        lng = max1+brnlis(k,6)
10      do i=1,6
          point(i) = brnlis(k-j+1,i)
        enddo
        do i=1,5
          brnlis(k-j+1,i) = brnlis(j,i)
        enddo
        brnlis(k-j+1,6) = lng-brnlis(j,i)
        do i=1,5
          brnlis(j,i) = point(i)
        enddo
        brnlis(j,6) = lng-point(6)
        j = j+1
        if(j.lt.k-j+1)goto 10
        if(j.eq.k-j+1)then
          brnlis(j,6) = lng-brnlis(j,6)
        endif
        return
      endif
c
      end
c
c
c
      subroutine movlis (lsfrom,lsto,ibeg,lim1,lim2,md2,
     +  maxboc,conlis,ierr)
c
      implicit none
c
      integer maxboc,ierr
      integer conlis(3,maxboc)
c
ccc      common/lists/conlis(3,20000),brnlis(5000,6)
ccc      integer conlis,brnlis
c
      integer ibeg,lim1,lim2,lsfrom,lsto,md2,m2
c
code ...
c
ccc      print *,'MOVLIS ',maxboc,lsfrom,lim1
c
      m2 = mod(conlis(lsfrom,lim1),2)
      conlis(lsto,ibeg) = conlis(lsfrom,lim1)+md2-m2
10    if(lim1.eq.lim2) goto 999
      ibeg = ibeg+1
      lim1 = lim1+1
      conlis(lsto,ibeg) = conlis(lsfrom,lim1)
      goto 10
c
999   continue
      ierr = 0
c
      return
      end
c
c
c
      subroutine invlis (lsfrom,lsto,ibeg,lim1,lim2,md2,
     +  maxbop,maxbob,maxboc,brnlis,conlis,ierr)
c
      implicit none
c
      integer maxbop,maxbob,maxboc,ierr
      integer brnlis(maxbob,6),conlis(3,maxboc)
c
ccc      common/lists/conlis(3,20000),brnlis(5000,6)
ccc      integer conlis,brnlis
c
      integer ibeg,lim1,lim2,lsfrom,lsto,iend,jend,md2,m2,icon
c
code ...
c
      ierr = 0
c
      iend = lim2-lim1+ibeg+1
      jend = iend
      m2 = mod(conlis(lsfrom,lim2),2)
      conlis(lsto,ibeg) = conlis(lsfrom,lim2)-m2+md2
c
10    if(lim1.eq.lim2) return
      ibeg = ibeg+1
      lim2 = lim2-1
      m2 = mod(conlis(lsfrom,lim2),2)
      if(m2.eq.1)then
        conlis(lsto,ibeg) = conlis(lsfrom,lim2)
        goto 10
      else
        icon = conlis(lsfrom,lim2)
        conlis(lsto,ibeg) = icon+1
        if(lim1.eq.lim2) return
20      lim2 = lim2-1
        iend = iend-1
        m2 = mod(conlis(lsfrom,lim2),2)
        if(m2.eq.1)then
          if(icon.eq.conlis(lsfrom,lim2)-m2) goto 50
        endif
        conlis(lsto,iend) = conlis(lsfrom,lim2)
        if(lim1.ne.lim2) goto 20
        goto 51

50      iend = iend+1
51      ibeg = ibeg+1
        if(iend.ne.jend)then
          conlis(lsto,ibeg) = conlis(lsto,iend)
          goto 50
        else
          conlis(lsto,ibeg) = conlis(lsfrom,lim2)-m2
          if(lim1.ne.lim2)then
            goto 10
          else
            ibeg = ibeg-1
            brnlis(1,5) = brnlis(1,5)-1
            return
          endif
        endif
      endif
c
      end
c
c
c
      subroutine squeze (lsfrom,lsto,fcount,ccount,lim1,lim2,
     +                   maxboc,conlis,ierr)
c
      implicit none
c
      integer maxboc,ierr
      integer conlis(3,maxboc)
c
ccc      common/lists/conlis(3,20000),brnlis(5000,6)
ccc      integer conlis,brnlis
c
      integer fcount,ccount,ccn,lim1,lim2,lsfrom,lsto
      integer even,odd,lim,md2,lng
c
      data even/0/,odd/1/
c
code ...
c
      ierr = 0
c
      lim = lim1
      ccn = ccount
      md2 = even
      lng = lim2-lim1+1
      fcount = fcount+1
      call movlis(lsfrom,lsto,fcount,lim1,lim2,md2,
     +      maxboc,conlis,ierr)
c
      if(lim2.lt.ccn) then
        lim1 = lim
        lim2 = lim2+1
        if(conlis(lsfrom,lim1).eq.conlis(lsfrom,lim2)+1) md2 = odd
        call movlis(lsfrom,lsfrom,lim1,lim2,ccn,md2,
     +      maxboc,conlis,ierr)
        lim2 = 1
      else
        lim2 = 0
      endif
      ccount = ccount-lng
c
      return
      end
c
c
c
      subroutine points(index,ix,iy,iz,ival,
     +                  maxbop,pntlis,levlis,a,ioxyz)
c
c      store point (orthogonal angstroms coordinates) and ival in
c      the point list.
c
      implicit none
c
      integer maxbop
      real pntlis(maxbop,3),a(3,3)
      integer levlis(maxbop)
c
ccc      common/point/pntlis(16000,3),a(3,3),levlis(16000)
c
ccc      common/lists/conlis(3,20000),brnlis(10000,6)
ccc      integer conlis,brnlis,levlis
c
ccc      common/page/ioxyz(3),iexxyz(3),iuvw(3),igrid(3),cell(6)
c
      integer ixyz(3),ioxyz(3)
      integer index,ix,iy,iz,ival,i,j
c
      real xyz(3)
c
code ...
c
      ixyz(1) = ix+ioxyz(1)-1
      ixyz(2) = iy+ioxyz(2)-1
      ixyz(3) = iz+ioxyz(3)-1
c
      do i=1,3
        xyz(i)=0
        do j=1,3
          xyz(i)=a(i,j)*ixyz(j)+xyz(i) 
        enddo
        pntlis(index,i)=xyz(i)
      enddo
      levlis(index) = ival
c
      return
      end
c
c
c
      subroutine mainlv (maxbop,maxbob,maxboc,npcn,nccn,
     +  mlng,levlis,conlis,ierr)
c
c  ---  85-07-12 s.thirup, dept. of molecular biology, uppsala.
c
c  ---  programme for changing the levels of bones points.
c      points in chains longer than mlng are assigned the value 3 (mainchain)
c      the rest of the points are assigned the value 2 (sidechain)
c
      implicit none
c
      integer maxbop,maxbob,maxboc,npcn,nccn,ierr
      integer levlis(maxbop),conlis(maxboc)
c
ccc      common/point/pntlis(16000,3),a(3,3),levlis(16000)
ccc      common/cons/npcn,nccn,conlis(20000)
ccc      integer npcn,nccn,conlis,levlis,nlab
c
      integer even,odd,unkn,side,main,nside,nmain
      integer pcount,ccount,i,icon,md2,inew,lng,nlev,mlng
c
      logical change
c
      data even/0/,odd/1/,unkn/1/,side/2/,main/3/
c
code ...
c
      nmain = 0
      nside = 0
c points are assigned level according to the length of the chain
      pcount = npcn
      ccount = nccn
      change = .true.
c initialization changed from f to t. mk 851021.
      i = 0
c
  100 i = i+1
      icon = conlis(i)
      md2 = mod(icon,2)
      if(md2.eq.even)then
        if(change)then
          inew = i
          lng = 0
          change = .false.
        else
          i = inew
          nlev = side
          change = .true.
          md2 = even
          call newlev(i,md2,nlev,nmain,nside,
     +      maxbop,maxboc,levlis,conlis)
        endif
      else
        if(change)then
          call newlev(i,md2,nlev,nmain,nside,
     +      maxbop,maxboc,levlis,conlis)
        else
          lng = lng+1
          if(lng.ge.mlng)then
            nlev = main
            i = inew
            change = .true.
            md2 = even
            call newlev(i,md2,nlev,nmain,nside,
     +        maxbop,maxboc,levlis,conlis)
          endif
        endif
      endif
      if(i.lt.ccount)goto 100
c
      if(i.eq.ccount)then
        if(.not.change)then
          nlev = side
          i = inew
          change = .true.
          md2 = even
          call newlev(i,md2,nlev,nmain,nside,
     +      maxbop,maxboc,levlis,conlis)
        endif
      endif
c
c write it out
c
cccc      write(4) npcn,nccn,nlab
cccc      write(4) ((pntlis(index,i),i=1,3),index=1,npcn)
cccc      write(4) (levlis(index),index=1,npcn)
cccc      write(4) (conlis(i),i=1,nccn)
c
      write(*,1000) mlng,nmain,nside
      ierr = 0
c
 1000 format (' Min length for main-chain fragments : ',i8/
     +        ' Nr of main-chain connections        : ',i8/
     +        ' Nr of side-chain connections        : ',i8)
c
      return
      end
c
c
c
      subroutine newlev(i,md2,nlev,nmain,nside,
     +      maxbop,maxboc,levlis,conlis)
c
      implicit none
c
      integer maxbop,maxboc
      integer levlis(maxbop),conlis(maxboc)
c
ccc      common/point/pntlis(16000,3),a(3,3),levlis(16000)
ccc      common/cons/pcount,ccount,conlis(20000)
ccc      integer pcount,ccount,conlis,levlis
c
      integer side,main,nside,nmain
      integer i,md2,nlev,icon,inx
c
      data side/2/,main/3/      
c
code ...
c
      icon = conlis(i)
      inx = nint((icon-md2)/2.0)
ccc      if(levlis(inx).lt.nlev)levlis(inx) = nlev
      levlis(inx) = nlev
      if(nlev.eq.side)nside = nside+1
      if(nlev.eq.main)nmain = nmain+1
c
      return
      end
c
c
c
      subroutine obones (maxbop,maxbob,maxboc,bonatm,bonlnk,
     +  bonxyz,oht,oib,out,mol,junk,ierr)
c
c ---      Write it out as a formatted o_db set of vectors
c
      implicit none
c
      integer maxbop,maxbob,maxboc,bonatm,bonlnk,out,ierr
      real bonxyz(maxbop,3)
      integer oht(maxbop),oib(maxboc)
c
ccc      common/point/bonxyz(16000,3),a(3,3),oht(16000)
ccc      common/cons/bonatm,bonlnk,oib(20000)
ccc      integer bonatm,bonlnk,oib,oht
ccc      real a, bonxyz
c
      character line*128, c6*6, fmt*15, mol*(*), param*25, type*1
      integer color(10), ct, ctskl(10), i, j, pt(2), junk(*)
c
      data color/6010, 12010, 18010, 24010, 30010, 36010, 4*0/
c
code ...
c
      ierr = 0
c
      do i=1,10
        ctskl (i) = 0
      end do
c
      do 100 i=1,bonatm
        if (oht(i) .le. 0) then
ccc          type*,' le 0 ', oht(i), i
        end if
        if (oht(i) .eq. 1) then
ccc          type*,' eq 1 ', oht(i), i
        end if
c
        if (oht(i) .le. 10) then
          ctskl(oht(i)) = ctskl(oht(i))+1
        end if
100   continue
c
      call prompt (' Histogram of bones types 1,2,3....')
      write (line, 20) ctskl
      call pretty (line)
      line = ' '//line
      call prompt (line)
      call prompt (' Type 2 = side-chain; 3 = main-chain')
      call prompt (' Type 1 = left-overs; 4 = spare parts')
c
      call remspa (mol)
      call upcase (mol)
c
      param = mol(1:5)//'_ATOM_XYZ'
      call remspa (param)
      type = 'R'
      ct = bonatm*3
      fmt = '(3f10.3)'
      write (out,30,err=9) param, type, ct, fmt
      write (out,fmt=fmt,err=9) ((bonxyz(j,i),i=1,3),j=1,bonatm)
c
      param = mol(1:5)//'_ATOM_BONE'
      call remspa (param)
      type = 'I'
      ct = bonatm
      fmt = '(35i2)'
      write (out,30,err=9) param, type, ct, fmt
      write (out,fmt=fmt,err=9) (oht(j),j=1,bonatm)
c
      param = mol(1:5)//'_CONNECTIVITY'
      call remspa (param)
      type = 'I'
      ct = bonlnk
      fmt = '(10i6)'
      write (out,30,err=9) param, type, ct, fmt
      write (out,fmt=fmt,err=9) (oib(j),j=1,bonlnk)
c
      param = mol(1:5)//'_ATOM_COLOUR'
      call remspa (param)
      type = 'I'
      ct = bonatm
      fmt = '(10i6)'
c
c --- get colour
c
      do 110 i=1,bonatm
        if (oht(i) .le. 10) then
          junk(i) = color(oht(i))
        else
          junk(i) = 0
        end if
110   continue
      write (out,30,err=9) param, type, ct, fmt
      write (out,fmt=fmt,err=9) (junk(j),j=1,bonatm)
c
      param = mol(1:5)//'_RESIDUE_NAME'
      call remspa (param)
      c6 = 'BONES'
      type = 'C'
      ct = 1
      fmt = '(5a)'
      write (out,30,err=9) param, type, ct, fmt
      write (out,60,err=9) c6      
c
      param = mol(1:5)//'_RESIDUE_TYPE'
      call remspa (param)
      c6 = 'SKL'
      type = 'C'
      ct = 1
      fmt = '(5a)'
      write (out,30,err=9) param, type, ct, fmt
      write (out,60,err=9) c6      
c
      param = mol(1:5)//'_RESIDUE_POINTERS'
      call remspa (param)
      pt(1) = 1
      pt(2) = bonatm
      type = 'I'
      ct = 2
      fmt = '(5i10)'
      write (out,30,err=9) param, type, ct, fmt
      write (out,50,err=9) pt
c
      write (out,*)
c
      return
c
c ... write error
c
    9 continue
      call errcon ('While writing BONES file')
      ierr = -1
      return
c
10    format (a)
20    format (' Bone type : ',10i8)
30    format (a,1x,a,i10,1x,a)
40    format (3f10.3)
50    format (5i10)
60    format (5a)
c
      end
c
c
c
      subroutine swappy (iunit,krecl,brick,ierr)
c
      integer iunit,krecl,ierr
      integer*2 record(256)
      integer extent(3),grid(3),origin(3),plus
      integer isca1,isca2,i,ir,nw,inx,iny,inz,nr
c
      real cell(6),prod
c
      character brick*(*)
c
code ...
c
      ierr = 0
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
      ir=0
10    continue
c
      ir=ir+1
      read(iunit,rec=ir,err=14)record
      call bytswp(record)
      write(iunit,rec=ir,err=15)record
      nw=ir
c
      if(ir.ne.1)go to 10
c
c ... header record
c
      isca1 = record(18)
      isca2 = record(19)
      if (isca1 .eq. 0) isca1 = 100
      if (isca2 .eq. 0) isca2 = 100
      do i=1,3
        origin(i) = record(i)
        extent(i) = record(i+3)
        grid  (i) = record(i+6)
        cell  (i) = float(record(i+9))/float(isca1)
        cell(i+3) = float(record(i+12))/float(isca1)
      end do
      prod = float(record(16))/float(isca2)
      plus = record(17)
c
      inx = extent(1)/8
      if(mod(extent(1),8) .ge. 1)inx = inx+1
      iny = extent(2)/8
      if(mod(extent(2),8) .ge. 1)iny = iny+1
      inz = extent(3)/8
      if(mod(extent(3),8) .ge. 1)inz = inz+1
c
      nr=(inx*iny*inz)+1
c
      write(*,30)origin,extent,grid,cell,prod,plus,
     +  isca1,isca2,inx,iny,inz,nr
c
30    format(' Origin          ',3i5,/,
     +       ' Extent          ',3i5,/,
     +       ' Grid            ',3i5,/,
     +       ' Cell    ',6f8.2,/,
     +       ' Prod/Plus       ',f7.2,1x,i5,/,
     +       ' Scale1/2        ',2i5/
     +       ' Bricks in x,y,z ',3i5/
     +       ' Nr of records   ',i10)
c
      go to 10
c
80    continue
      call jvalut (' Total number of records written :',1,nw)
      close (iunit)
c
      return
c
   13 continue
      call errcon ('While opening BRICK file')
      close (iunit)
      ierr = -1
      return
c
   14 continue
      if (nw .eq. nr) goto 80
      call errcon ('While reading BRICK file')
      call jvalut (' Record nr :',1,ir)
      close (iunit)
      ierr = -2
      return
c
   15 continue
      call errcon ('While writing BRICK file')
      call jvalut (' Record nr :',1,ir)
      close (iunit)
      ierr = -3
      return
c
      end
c
c
c
      subroutine pbones (maxbop,maxbob,maxboc,bonatm,bonlnk,
     +  bonxyz,olev,ocon,iunit,onbr,cntr,bonlen,bobfac,cell,
     +  ierr)
c
c --- prune BONES and write PDB file
c
      implicit none
c
      integer maxlev
      parameter (maxlev=1000)
c
      integer maxbop,maxbob,maxboc,bonatm,bonlnk,iunit,bonlen,ierr
      real bonxyz(maxbop,3)
      integer olev(maxbop),ocon(maxboc),onbr(maxbob,6)
      integer done(maxlev),cntr(*)
c
      real cell(6),ca2fa(12),bobfac
c
      integer i,j,k,l,prev,st(0:6),ct(2,10),leng1
c
      logical succes
c
      character line*128
c
code ...
c
      ierr = 0
c
c ... set up neighbour arrays
c
      do i=1,bonatm
        do j=1,6
          onbr (i,j) = 0
        end do
        cntr (i) = 0
      end do
c
      prev = -1
      do i=1,bonlnk
        j = ocon(i)
        k = j/2
        l = 2*k
        if (l .eq. j) then
c
c ... start of a new fragment
c
          prev = k
c
        else
c
c ... continuation of a fragment
c
          cntr(prev)=cntr(prev)+1
          onbr (prev,cntr(prev)) = k
          cntr(k)=cntr(k)+1
          onbr (k,cntr(k)) = prev
          prev = k
        end if
      end do
c
      do i=0,6
        st (i) = 0
      end do
c
      do i=1,bonatm
        st (cntr(i)) = st(cntr(i)) + 1
      end do
c
      do i=0,6
        write (*,1000) i,st(i)
      end do
 1000 format (' Nr of BONES with ',i1,' neighbours = ',i8)
c
c ... keep only side-chain atoms connected to the main chain
c
      call prompt (' Pruning BONES ...')
      write (*,*)
c
      do i=1,bonatm
        if (olev(i) .eq. 2) then
ccc          print *,i
          succes = .false.
          j = i
          k = 1
          done (1) = i
          call iscona (j,k,bonlen,succes,cntr,onbr,olev,
     +                 done,maxbob,maxbop)
c
          if (succes) olev(i) = -2
ccc          if (succes) print *,' OKAY ! - ',i
        end if
      end do
c
c ... select all main chain atoms (code=3)
c
      do i=1,bonatm
        if (olev(i) .eq. 3) olev(i) = -3
      end do
c
c ... write pruned BONES to a PDB file
c
      do i=1,10
        ct (1,i) = 0
        ct (2,i) = 0
      end do
c
      call stamp (line)
      write (iunit,'(a6,1x,a)') 'REMARK',line(1:leng1(line))
c
      write (iunit,'(a,3f9.3,3f7.2)') 'CRYST1',(cell(i),i=1,6)
      write (iunit,'(a)')
     +  'ORIGX1      1.000000  0.000000  0.000000        0.00000'
      write (iunit,'(a)')
     +  'ORIGX2      0.000000  1.000000  0.000000        0.00000'
      write (iunit,'(a)')
     +  'ORIGX3      0.000000  0.000000  1.000000        0.00000'
c
      call orthog (cell,ca2fa,1)
      write (iunit,'(a6,4x,3f10.6,f15.5)')
     +  'SCALE1',(ca2fa(i),i=1,10,3)
      write (iunit,'(a6,4x,3f10.6,f15.5)')
     +  'SCALE2',(ca2fa(i),i=2,11,3)
      write (iunit,'(a6,4x,3f10.6,f15.5)')
     +  'SCALE3',(ca2fa(i),i=3,12,3)
c
      l = 0
      do i=1,bonatm
        k = olev(i)
        if (k .gt. 0) then
          ct (1,k) = ct (1,k) + 1
        else
          k = -k
          ct (2,k) = ct (2,k) + 1
          olev (i) = k
c
          l = l + 1
          write (iunit,6010) l,' CA ','ALA',l,bonxyz(i,1),
     +      bonxyz(i,2),bonxyz(i,3),1.0,bobfac
        end if
      end do
c
      call jvalut (' Nr of BONES written to PDB file :',1,l)
c
      call prompt (' Nr of BONES with code 1,2,... written:')
      write (line,'(10i8)') (ct(2,i),i=1,10)
      call pretty (line)
      call textut (' Nr : ',line)
c
      call prompt (' Nr of BONES with code 1,2,... NOT written:')
      write (line,'(10i8)') (ct(1,i),i=1,10)
      call pretty (line)
      call textut (' Nr : ',line)
c
 6010 format ('ATOM  ',i5,1x,a4,1x,a3,1x,i5,4x,3f8.3,2f6.2)
c
      return
      end
c
c
c
      subroutine iscona (inow,idepth,mdepth,succes,
     +                   cntr,onbr,olev,done,maxbob,maxbop)
c
c ... recursive procedure to check if a side-chain
c     BONES atom is somehow connected to the main-chain
c
      integer inow,maxbob,maxbop,idepth,mdepth
      integer olev(maxbop),onbr(maxbob,6),done(*),cntr(*)
c
      integer i,j,inbr,jdepth,iindex
c
      logical succes
c
code ...
c
c ... if SUCCES, unwind recursion
c
      if (succes) return
c
c ... don't go too deep (avoid loops)
c
      if (idepth .ge. mdepth) return
c
c ... if this BONE has no neighbours, forget it
c
      if (cntr(inow) .le. 0) return
c
c ... see if any of its neighbours is a main-chain atom
c
ccc      print *,' INOW = ',inow
      do i=1,cntr(inow)
c
        if (i .le. 6) then
          inbr = onbr(inow,i)
c
c ... the ALPHA version may crash without the following WRITE
c     statement ....
c
ccc          write (*,6996) inow,idepth,cntr(inow),inbr,olev(inbr)
 6996     format ('+',5i10)
c
          if (inbr .gt. 0) then
            if (olev(inbr) .eq. 3) then
               succes = .true.
               return
            end if
          end if
        else
          print *,'OOPS - I,INOW,CNTR,IDEPTH = ',
     +      i,inow,cntr(inow),idepth
        end if
c
      end do
c
c ... if no succes, try the neighbours' neighbours
c
      jdepth = idepth + 1
ccc      print *,' RECURSION ',jdepth
      do i=1,cntr(inow)
        inbr = onbr(inow,i)
c
c ... the ALPHA version may crash without the following WRITE
c     statement ....
c
ccc          write (*,6996) inow,idepth,cntr(inow),inbr,olev(inbr)
c
c ... check if this neighbour not already visited (to avoid
c     jumping back and forth between neighbouring bones atoms)
c
        if (inbr .gt. 0) then
          j = iindex (inbr,0,idepth,done)
c
          if (j .le. 0) then
            done (jdepth) = inbr
cc          print *,' TRY - ',inbr
            call isconb (inbr,jdepth,mdepth,succes,cntr,onbr,
     +                   olev,done,maxbob,maxbop)
            if (succes) return
          end if
        end if
c
      end do
c
      return
      end
c
c
c
      subroutine isconb (inow,idepth,mdepth,succes,
     +                   cntr,onbr,olev,done,maxbob,maxbop)
c
c ... fake routine to enforce recursion
c
      integer inow,maxbob,maxbop,idepth,mdepth
      integer olev(maxbop),onbr(maxbob,6),done(*),cntr(*)
c
      logical succes
c
code ...
c
      call iscona (inow,idepth,mdepth,succes,
     +             cntr,onbr,olev,done,maxbob,maxbop)
c
      return
      end
c
c
c
      subroutine x3dmin (iunit,map,ext1,ext2,ext3,ierr)
c
      implicit none
c
      integer ext1,ext2,ext3,iunit,ierr,i,j,k
c
      real map(ext1,ext2,ext3)
c
code ...
c
      ierr = -1
c
      read (iunit,*,err=10,end=20)
      do i=1,ext1
        read (iunit,*,err=10,end=20)
        do j=1,ext2
          read (iunit,*,err=10,end=20)
          read (iunit,*,err=10,end=20) (map(i,j,k),k=1,ext3)
          read (iunit,*,err=10,end=20)
        end do
        read (iunit,*,err=10,end=20)
      end do
      read (iunit,*,err=10,end=20)
c
      ierr = 0
c
      return
c
   10 continue
      call errcon ('While reading matrix')
      ierr = -2
      print *,i,j,k
      return
c
   20 continue
      call errcon ('End-of-file while reading matrix')
      ierr = -3
      print *,i,j,k
      return
c
      end
c
c
c
      subroutine gethis (map,np,dmin,step,cnt,nbins,ierr)
c
      implicit none
c
      integer np,nbins
c
      real map (np)
      real dmin,step
c
      integer cnt(nbins)
      integer i,idum,ierr
c
code ...
c
      ierr = 0
c
      do i=1,nbins
        cnt (i) = 0
      end do
c
      do i=1,np
        idum = 1 + int ( (map(i)-dmin)/step)
        if (idum .le. 0 .or. idum .gt. nbins) then
          call errcon ('Programming error; wrong bin')
          call jvalut (' Calculated bin :',1,idum)
          ierr = -1
          return
        end if
        cnt(idum) = cnt(idum) + 1
      end do
c
      return
      end
c
c
c
      subroutine rfur1 (map,iunit,nx,ny,nz,xdum,maxlen,ierr)
c
c ... read PHASES map
c
      implicit none
c
      integer nx,ny,nz,maxlen
      real map(nx,ny,nz),xdum(maxlen)
c
      integer iunit,ierr,iy,iz,i
c
code ...
c
      ierr = -1
c
      write (*,6000) ny*nz,nx
 6000 format (' Reading ',i8,' records of ',i5,' points ...')
c
      do iy=1,ny
        do iz=1,nz
          read (iunit,end=9000,err=9000) (xdum(i),i=1,nx)
          do i=1,nx
            map (i,iy,iz) = xdum (i)
          end do
        end do
      end do
c
      ierr = 0
      return
c
 9000 continue
      ierr = -1
      call errcon ('While reading')
      call jvalut (' Approximately at record :',1,(iy*nz + iz))
      return
c
      end
c
c
c
      subroutine rfur2 (map,iunit,nx,ny,nz,mybyte,maxlen,ierr)
c
c ... read PHASES mask
c
      implicit none
c
      integer nx,ny,nz,maxlen
      real map(nx,ny,nz)
      byte mybyte(maxlen)
c
      integer iunit,ierr,iy,iz,i,ix
c
code ...
c
      ierr = -1
c
      write (*,6000) ny*nz,nx
 6000 format (' Reading ',i8,' records of ',i5,' points ...')
c
      do iy=1,ny
        do iz=1,nz
          read (iunit,end=9000,err=9000) (mybyte(i),i=1,nx)
          do i=1,nx
            ix = mybyte(i)
            map (i,iy,iz) = float(ix)
          end do
        end do
      end do
c
      ierr = 0
      return
c
 9000 continue
      ierr = -1
      call errcon ('While reading')
      call jvalut (' Approximately at record :',1,(iy*nz + iz))
      return
c
      end
c
c
c
      subroutine dopro1 (map,ex1,ex2,ex3,buffer,ityp) 
c
      implicit none
c
      integer ex1,ex2,ex3
c
      real map(ex1,ex2,ex3),buffer(*)
c
      integer ityp,i,j,k
c
code ...
c
      if (ityp .eq. 1) then
        do i=1,ex1
          buffer (i) = 0.0
          do j=1,ex2
            do k=1,ex3
              buffer (i) = buffer (i) + map(i,j,k)
            end do
          end do
        end do
c
      else if (ityp .eq. 2) then
        do i=1,ex2
          buffer (i) = 0.0
          do j=1,ex1
            do k=1,ex3
              buffer (i) = buffer (i) + map(j,i,k)
            end do
          end do
        end do
c
      else if (ityp .eq. 3) then
        do i=1,ex3
          buffer (i) = 0.0
          do j=1,ex1
            do k=1,ex2
              buffer (i) = buffer (i) + map(j,k,i)
            end do
          end do
        end do
c
      else
        call errcon ('DOPRO1 - Invalid projection axis')
      end if
c
      return
      end
c
c
c
      subroutine dopro2 (map,ex1,ex2,ex3,buffer,ityp,i1,i2) 
c
      implicit none
c
      integer ex1,ex2,ex3
c
      real map(ex1,ex2,ex3),buffer(*)
c
      integer ityp,i,j,k,i1,i2,ib
c
code ...
c
      if (ityp .eq. 1) then
        do i=1,ex2
          do j=1,ex3
            ib = i+ex2*(j-1)
            buffer (ib) = 0.0
            do k=i1,i2
              buffer (ib) = buffer (ib) + map(k,i,j)
            end do
          end do
        end do
c
      else if (ityp .eq. 2) then
        do i=1,ex1
          do j=1,ex3
            ib = i+ex1*(j-1)
            buffer (ib) = 0.0
            do k=i1,i2
              buffer (ib) = buffer (ib) + map(i,k,j)
            end do
          end do
        end do
c
      else if (ityp .eq. 3) then
        do i=1,ex1
          do j=1,ex2
            ib = i+ex1*(j-1)
            buffer (ib) = 0.0
            do k=i1,i2
              buffer (ib) = buffer (ib) + map(i,j,k)
            end do
          end do
        end do
c
      else
        call errcon ('DOPRO2 - Invalid projection direction')
      end if
c
      return
      end
c
c
c
      subroutine filter (how,map,buf,nc,nx,ny,nz,alpha,beta,gamma,
     +                   kuse,ierr)
c
c ... apply filters to a map
c
      implicit none
c
      integer maxbuf
      parameter (maxbuf=10000)
c
      integer nx,ny,nz,nc,ierr,kuse
c
      real map(nx,ny,nz),buf(nx,ny,nz)
      real tmp(maxbuf),conv(-1:1,-1:1,-1:1)
      real sumx,sumxsq,f,q,alpha,beta,gamma,bg
c
      integer inx(maxbuf)
      integer i,j,k,l,m,n,knt
c
      character how*(*)
c
code ...
c
      ierr = -1
      if (nc .le. 0) then
        call errcon ('Cube size must be positive')
        return
      else if (nc*nc*nc .gt. maxbuf) then
        call errcon ('Cube too big for buffer')
        return
      else if (2*nc .gt. nx) then
        call errcon ('Cube too large for X direction')
        return
      else if (2*nc .gt. ny) then
        call errcon ('Cube too large for Y direction')
        return
      else if (2*nc .gt. nz) then
        call errcon ('Cube too large for Z direction')
        return
      end if
      if (kuse .le. 1) then
        call errcon ('K should greater than 1')
        return
      else if (kuse .gt. (2*nc+1)**3) then
        call errcon ('Value of K larger than nr of points in cube')
        return
      end if
c
      call inimap (buf,nx,ny,nz,0.0)
c
c ... MINIMUM
c
      if (how (1:2) .eq. 'MI') then
        call prompt (' MInimum - working ...')
        call jvalut (' Cube size :',1,nc)
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              q = map(i,j,k)
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    if (map(l,m,n).lt.q) q = map(l,m,n)
                  end do
                end do
              end do
              buf (i,j,k) = q
            end do
          end do
        end do
c
c ... MAXIMUM
c
      else if (how (1:2) .eq. 'MA') then
        call prompt (' MAximum - working ...')
        call jvalut (' Cube size :',1,nc)
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              q = map(i,j,k)
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    if (map(l,m,n).gt.q) q = map(l,m,n)
                  end do
                end do
              end do
              buf (i,j,k) = q
            end do
          end do
        end do
c
c ... MEDIAN
c
      else if (how (1:2) .eq. 'ME') then
        call prompt (' MEdian - working ...')
        call jvalut (' Cube size :',1,nc)
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    tmp (knt) = map(l,m,n)
                    inx (knt) = knt
                  end do
                end do
              end do
c
c ... sort density values
c
ccc      if (tmp(1) .ne. 0.0) then
ccc      print *,' I,J,K = ',i,j,k
ccc      call fvalut (' Dens :',(2*nc+1)**3,tmp)
ccc      call ivalut (' Indx :',(2*nc+1)**3,inx)
ccc      end if
c
              call shellr (tmp,inx,knt)
c
ccc      if (tmp(1) .ne. 0.0) then
ccc      call fvalut (' Aftr :',(2*nc+1)**3,tmp)
ccc      call ivalut (' Sort :',(2*nc+1)**3,inx)
ccc      l = knt/2
ccc      call fvalut (' Median :',1,tmp(l))
ccc      end if
c
c ... take median
c
              knt = knt/2
              buf (i,j,k) = tmp(knt)
c
            end do
          end do
        end do
c
c ... K LOWEST
c
      else if (how (1:2) .eq. 'KL') then
        call prompt (' K Lowest - working ...')
        call jvalut (' Cube size  :',1,nc)
        call jvalut (' Value of K :',1,kuse)
        f = 1.0 / float (kuse)
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    tmp (knt) = map(l,m,n)
                    inx (knt) = knt
                  end do
                end do
              end do
c
c ... sort density values
c
              call shellr (tmp,inx,knt)
c
c ... average KUSE lowest
c
              q = 0.0
              do l=1,kuse
                q =q + tmp(l)
              end do
c
              buf (i,j,k) = q * f
c
            end do
          end do
        end do
c
c ... K HIGHEST
c
      else if (how (1:2) .eq. 'KH') then
        call prompt (' K Highest - working ...')
        call jvalut (' Cube size  :',1,nc)
        call jvalut (' Value of K :',1,kuse)
        f = 1.0 / float (kuse)
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    tmp (knt) = map(l,m,n)
                    inx (knt) = knt
                  end do
                end do
              end do
c
c ... sort density values
c
              call shellr (tmp,inx,knt)
c
c ... average KUSE highest
c
              q = 0.0
              do l=knt,knt-kuse+1,-1
                q =q + tmp(l)
              end do
c
              buf (i,j,k) = q * f
c
            end do
          end do
        end do
c
c ... GRADIENT
c
      else if (how (1:2) .eq. 'GR') then
        nc = 1
        call prompt (' GRadient - working ...')
c
c ... X
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              buf (i,j,k) = map(i+1,j,k) - map(i-1,j,k)
            end do
          end do
        end do
c
c ... Y
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              map (i,j,k) = buf(i,j+1,k) - buf(i,j-1,k)
            end do
          end do
        end do
c
c ... Z
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              buf (i,j,k) = map(i,j,k+1) - map(i,j,k-1)
            end do
          end do
        end do
c
c ... LAPLACE
c
      else if (how (1:2) .eq. 'LA') then
        nc = 1
        call prompt (' LAplace - working ...')
c
c ... X
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              buf (i,j,k) = 2.0*map(i,j,k) - map(i+1,j,k) - map(i-1,j,k)
            end do
          end do
        end do
c
c ... Y
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              map (i,j,k) = 2.0*buf(i,j,k) - buf(i,j+1,k) - buf(i,j-1,k)
            end do
          end do
        end do
c
c ... Z
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              buf (i,j,k) = 2.0*map(i,j,k) - map(i,j,k+1) - map(i,j,k-1)
            end do
          end do
        end do
c
c ... AVERAGE
c
      else if (how (1:2) .eq. 'AV') then
        call prompt (' AVerage - working ...')
        call jvalut (' Cube size :',1,nc)
        f = 1.0 / float( (2*nc+1)**3 )
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              sumx = 0.0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    sumx = sumx + map(l,m,n)
                  end do
                end do
              end do
              buf (i,j,k) = sumx * f
            end do
          end do
        end do
c
c ... SIGNAL/NOISE
c
      else if (how (1:2) .eq. 'SI') then
        call prompt (' SIgnal/noise - working ...')
        call jvalut (' Cube size :',1,nc)
        f = 1.0 / float( (2*nc+1)**3 )
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              sumx = 0.0
              sumxsq = 0.0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    sumx = sumx + map(l,m,n)
                    sumxsq = sumxsq + map(l,m,n)*map(l,m,n)
                  end do
                end do
              end do
              q = f * (sumxsq - sumx*sumx*f)
              if (q .le. 0.0) then
                buf (i,j,k) = 0.0
              else
                buf (i,j,k) = sumx*f / sqrt(q)
              end if
            end do
          end do
        end do
c
c ... STATISTICAL DIFFERENCING
c
      else if (how (1:2) .eq. 'ST') then
        call prompt (' STatistical differencing - working ...')
        call jvalut (' Cube size :',1,nc)
        f = 1.0 / float( (2*nc+1)**3 )
        bg = beta*gamma
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              sumx = 0.0
              sumxsq = 0.0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    sumx = sumx + map(l,m,n)
                    sumxsq = sumxsq + map(l,m,n)*map(l,m,n)
                  end do
                end do
              end do
              q = f * (sumxsq - sumx*sumx*f)
              if (q .le. 0.0) then
                buf (i,j,k) = 0.0
              else
                sumx = sumx * f
                buf (i,j,k) = (1.0-alpha)*sumx +
     +            (map(i,j,k)-sumx)*bg/(gamma+beta*sqrt(q))
              end if
            end do
          end do
        end do
c
c ... VARIABLE THRESHOLD
c
      else if (how (1:2) .eq. 'VA') then
        call prompt (' VAriable threshold - working ...')
        call jvalut (' Cube size :',1,nc)
        f = 1.0 / float( (2*nc+1)**3 )
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              knt = 0
              sumx = 0.0
              sumxsq = 0.0
              do n=k-nc,k+nc
                do m=j-nc,j+nc
                  do l=i-nc,i+nc
                    knt = knt + 1
                    sumx = sumx + map(l,m,n)
                    sumxsq = sumxsq + map(l,m,n)*map(l,m,n)
                  end do
                end do
              end do
              q = f * (sumxsq - sumx*sumx*f)
              if (q .le. 0.0) then
                buf (i,j,k) = 0.0
              else
                sumx = sumx * f
                buf (i,j,k) = map(i,j,k) -
     +            (alpha*sqrt(q) + beta*sumx)
              end if
            end do
          end do
        end do
c
c ... SMOOTH
c
      else if (how (1:2) .eq. 'SM') then
        nc = 1
        call prompt (' SMooth - working ...')
        f = 0.0
        do i=-1,1
          do j=-1,1
            do k=-1,1
              conv (i,j,k) = 1.0
              if (i .eq. 0) conv(i,j,k) = conv(i,j,k) + 1.0
              if (j .eq. 0) conv(i,j,k) = conv(i,j,k) + 1.0
              if (k .eq. 0) conv(i,j,k) = conv(i,j,k) + 1.0
              if (i.eq.0 .and. j.eq.0 .and. k.eq.0) conv(i,j,k) = 5.0
              f = f + conv(i,j,k)
            end do
          end do
        end do
        call fvalut (' Filter :',27,conv)
        call fvalut (' Sum    :',1,f)
        f = 1.0 / f
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              q = 0.0
              do n=-nc,nc
                do m=-nc,nc
                  do l=-nc,nc
                    q = q + conv (l,m,n) * map(i+l,j+m,k+n)
                  end do
                end do
              end do
              buf (i,j,k) = q * f
            end do
          end do
        end do
c
c ... EDGE
c
      else if (how (1:2) .eq. 'ED') then
        nc = 1
        call prompt (' EDge - working ...')
        call fvalut (' Value for A :',1,alpha)
        call fvalut (' Value for B :',1,beta)
        f = 0.0
        do i=-1,1
          do j=-1,1
            do k=-1,1
              conv (i,j,k) = alpha
              f = f + conv(i,j,k)
            end do
          end do
        end do
        f = f - alpha
        conv (0,0,0) = beta - 26.0*alpha
        f = f + conv (0,0,0)
        do i=-1,1
          do j=-1,1
            do k=-1,1
              conv (i,j,k) = conv (i,j,k) / f
            end do
          end do
        end do
        call fvalut (' Filter :',27,conv)
c
        do k=nc+1,nz-nc
          do j=nc+1,ny-nc
            do i=nc+1,nx-nc
              q = 0.0
              do n=-nc,nc
                do m=-nc,nc
                  do l=-nc,nc
                    q = q + conv (l,m,n) * map(i+l,j+m,k+n)
                  end do
                end do
              end do
              buf (i,j,k) = q
            end do
          end do
        end do
c
c ... invalid option
c
      else
        call errcon ('Invalid filter option')
        return
      end if
c
c ... copy back to map
c
      call copmap (map,buf,nx,ny,nz)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine eramap (map,ori,nx,ny,nz,lim)
c
      implicit none
c
      integer nx,ny,nz,ori(3),lim(2,3),i,j,k,nnn
c
      real map(nx,ny,nz),sum,asm
c
code ...
c
      sum = 0.0
      asm = 0.0
      nnn = 0
c
      do i = lim(1,1)-ori(1)+1,lim(2,1)-ori(1)+1
        do j = lim(1,2)-ori(2)+1,lim(2,2)-ori(2)+1
          do k = lim(1,3)-ori(3)+1,lim(2,3)-ori(3)+1
            sum = sum + map(i,j,k)
            asm = asm + abs(map(i,j,k))
            nnn = nnn + 1
            map (i,j,k) = 0.0
          end do
        end do
      end do
c
      call jvalut (' Nr of points set to zero :',1,nnn)
      call rvalut (' Sum of zeroed density    :',1,sum)
      call rvalut (' Sum of | zeroed dens |   :',1,asm)
      if (nnn .gt. 0) then
        call rvalut (' Average zeroed density   :',1,sum/float(nnn))
        call rvalut (' Average | zeroed dens |  :',1,asm/float(nnn))
      end if
c
      return
      end
c
c
c
      subroutine intmap (map,ori,nx,ny,nz,lim)
c
      implicit none
c
      integer nx,ny,nz,ori(3),lim(2,3),i,j,k,nnn
c
      real map(nx,ny,nz),sum,asm
c
code ...
c
      sum = 0.0
      asm = 0.0
      nnn = 0
c
      do i = lim(1,1)-ori(1)+1,lim(2,1)-ori(1)+1
        do j = lim(1,2)-ori(2)+1,lim(2,2)-ori(2)+1
          do k = lim(1,3)-ori(3)+1,lim(2,3)-ori(3)+1
            sum = sum + map(i,j,k)
            asm = asm + abs(map(i,j,k))
            nnn = nnn + 1
          end do
        end do
      end do
c
      call jvalut (' Nr of points summed :',1,nnn)
      call rvalut (' Sum of density      :',1,sum)
      call rvalut (' Sum of | density |  :',1,asm)
      if (nnn .gt. 0) then
        call rvalut (' Average zeroed density   :',1,sum/float(nnn))
        call rvalut (' Average | zeroed dens |  :',1,asm/float(nnn))
      end if
c
      return
      end
c
c
c
      subroutine values (map,ori,nx,ny,nz,lim)
c
      implicit none
c
      integer nx,ny,nz,ori(3),lim(2,3),i,j,k,nnn
c
      real map(nx,ny,nz),sum,asm
c
      character answer*1
c
code ...
c
      nnn = (lim(2,1)-lim(1,1)+1)*(lim(2,2)-lim(1,2)+1)*
     +      (lim(2,3)-lim(1,3)+1)
      call jvalut (' Nr of map points to list :',1,nnn)
c
      if (nnn .gt. 50) then
        answer = 'N'
        call textin (' Are you SURE (Y/N) ?',answer)
        call upcase (answer)
        if (answer .ne. 'Y') return
      end if
c
      sum = 0.0
      asm = 0.0
      nnn = 0
c
 6000 format (/1x,'    X','    Y','    Z','      Density')
 6010 format (1x,i5,i5,i5,1x,1pe12.4,3i5)
c
ccc      print *,' I limits ',lim(1,1)-ori(1)+1,lim(2,1)-ori(1)+1
ccc      print *,' J limits ',lim(1,2)-ori(2)+1,lim(2,2)-ori(2)+1
ccc      print *,' K limits ',lim(1,3)-ori(3)+1,lim(2,3)-ori(3)+1
c
      do i = lim(1,1)-ori(1)+1,lim(2,1)-ori(1)+1
        write (*,6000)
        do j = lim(1,2)-ori(2)+1,lim(2,2)-ori(2)+1
          do k = lim(1,3)-ori(3)+1,lim(2,3)-ori(3)+1
            sum = sum + map(i,j,k)
            asm = asm + abs(map(i,j,k))
            nnn = nnn + 1
            write (*,6010) i+ori(1)-1,j+ori(2)-1,
     +        k+ori(3)-1,map(i,j,k)
          end do
        end do
      end do
c
      write (*,*)
      call jvalut (' Nr of points listed :',1,nnn)
      call rvalut (' Sum of density      :',1,sum)
      call rvalut (' Sum of | density |  :',1,asm)
      if (nnn .gt. 0) then
        call rvalut (' Average density     :',1,sum/float(nnn))
        call rvalut (' Average | dens |    :',1,asm/float(nnn))
      end if
c
      return
      end
c
c
c
      subroutine pastem (map1,nx,ny,nz,ori1,lim1,
     +                   map2,mx,my,mz,ori2,lim2,ierr)
c
      implicit none
c
      integer nx,ny,nz,mx,my,mz,ori1(3),ori2(3),lim1(2,3),lim2(3)
      integer i,j,k,ext1(3),ext2(3),ierr,nnn,i1,j1,k1,i2,j2,k2
c
      real map1(nx,ny,nz),map2(mx,my,mz),sumold,sumnew
c
      logical sorry
c
      character xyz(3)*1
c
      data xyz /'X','Y','Z'/
c
code ...
c
      ierr = -1
c
      ext1(1)=nx
      ext1(2)=ny
      ext1(3)=nz
c
      ext2(1)=mx
      ext2(2)=my
      ext2(3)=mz
c
 6000 format (' Invalid ',a1,' limits - FROM = ',2i6,' TO = ',2i6)
c
      sorry = .false.
      do i=1,3
        if (lim1(1,i) .lt. ori1(i) .or.
     +      lim1(2,i) .gt. (ori1(i)+ext1(i)-1) .or.
     +      lim2(i) .lt. ori2(i) .or.
     +      (lim2(i)+lim1(2,i)-lim1(1,i)) .gt. 
     +           (ori2(i)+ext2(i)-1)) then
          write (*,6000) xyz(i),lim1(1,i),lim1(2,i),
     +      lim2(i),(lim2(i)+lim1(2,i)-lim1(1,i))
          sorry = .true.
        end if
      end do
      if (sorry) then
        ierr = -2
        call errcon ('Invalid limits - Sorry')
        return
      end if
c
      nnn = 0
      sumold = 0.0
      sumnew = 0.0
      do i = lim1(1,1),lim1(2,1)
        i1 = i-ori1(1)+1
        i2 = lim2(1) + (i-lim1(1,1)) - ori2(1) + 1
ccc		print *,' I, I1, I2 = ',i,i1,i2
        do j = lim1(1,2),lim1(2,2)
          j1 = j-ori1(2)+1
          j2 = lim2(2) + (j-lim1(1,2)) - ori2(2) + 1
ccc            if (i .eq. lim1(1,1))
ccc     +        print *,' J, J1, J2 = ',j,j1,j2
          do k = lim1(1,3),lim1(2,3)
            k1 = k-ori1(3)+1
            k2 = lim2(3) + (k-lim1(1,3)) - ori2(3) + 1
ccc            if (i .eq. lim1(1,1))
ccc     +        print *,' K, K1, K2 = ',k,k1,k2
            nnn = nnn + 1
            sumold = sumold + map2(i2,j2,k2)
            sumnew = sumnew + map1(i1,j1,k1)
            map2(i2,j2,k2) =  map1(i1,j1,k1)
          end do
        end do
      end do
c
      call jvalut (' Nr of points copied :',1,nnn)
      call rvalut (' Sum of OLD density  :',1,sumold)
      call rvalut (' Sum of NEW density  :',1,sumnew)
c
      if (nnn .gt. 0) ierr = 0
c
      return
      end
c
c
c
      subroutine combim (exoper,map1,na1,nb1,nc1,
     +                   map2,na2,nb2,nc2,lim1,lim2)
c
c ... combine maps
c
      implicit none
c
      integer na1,nb1,nc1,na2,nb2,nc2,lim1(2,3),lim2(2,3)
      integer i1,i2,i3,io1,io2,io3,ntot
c
      real map1(na1,nb1,nc1),map2(na2,nb2,nc2)
c
      real avx,avy,avxy,avxsq,avysq,sumfo,sumfc,fofc,rmsd
      real avalue,bvalue,aabs,babs,f,ccoef,shape,q1,q2,gamma
c
      character exoper*(*),oper*(10)
c
code ...
c
      oper = exoper
c
      io1 = lim2(1,1) - lim1(1,1)
      io2 = lim2(1,2) - lim1(1,2)
      io3 = lim2(1,3) - lim1(1,3)
c
c ... SIMILARITY M1 - M2 (note: should be done in double precision !)
c
      if (oper(1:3) .eq. 'SIM') then
c
        ntot = 0
        avx = 0.
        avy = 0.
        avxy = 0.
        avxsq = 0.
        avysq = 0.
        sumfo = 0.
        sumfc = 0.
        fofc = 0.
        rmsd = 0.
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              ntot = ntot + 1
              avalue = map1(i1,i2,i3)
              bvalue = map2(i1+io1,i2+io2,i3+io3)
              aabs = abs(avalue)
              babs = abs(bvalue)
              avx = avx+ avalue
              avy = avy+ bvalue
              avxsq = avxsq+ avalue**2
              avysq = avysq+ bvalue**2
              avxy = avxy+ avalue*bvalue
              fofc = fofc + abs(aabs-babs)
              sumfo = sumfo + aabs
              sumfc = sumfc + babs
              rmsd = rmsd + (avalue-bvalue)**2
            end do
          end do
        end do
c
        call jvalut (' Nr of points in common grid :',1,ntot)
        if (ntot .le. 0) return
c
        f = float(ntot)
        gamma = 
     +    (sqrt(avxsq/f- (avx/f)**2)* sqrt(avysq/f- (avy/f)**2))
        if (gamma .ne. 0.0) then
          ccoef = (avxy/f- avx*avy/(f*f))/ gamma
          call fvalut (' Correlation coefficient :',1,ccoef)
        else
          if ( (avxy/f- avx*avy/(f*f)) .eq. 0.0) then
            ccoef = 1.0
            call fvalut (' Correlation coefficient :',1,ccoef)
          else
            call errcon ('Cannot compute correlation coefficient')
          end if
        end if
c
        if (sumfo .ne. 0.0) then
          sumfo = fofc/sumfo
          call fvalut (' R-factor w.r.t. map 1   :',1,sumfo)
        else
          call errcon ('Cannot compute R-factor w.r.t. map 1')
        end if
c
        if (sumfc .ne. 0.0) then
          sumfc = fofc/sumfc
          call fvalut (' R-factor w.r.t. map 2   :',1,sumfc)
        else
          call errcon ('Cannot compute R-factor w.r.t. map 2')
        end if
c
        rmsd = sqrt (rmsd/f)
        call rvalut (' RMS difference          :',1,rmsd)
c
        gamma = (sqrt(avxsq)*sqrt(avysq))
        if (gamma .ne. 0.0) then
          shape = avxy / gamma
          call fvalut (' Shape similarity index  :',1,shape)
        else
          call errcon ('Cannot compute shape similarity index')
        end if
c
        call prompt (' R-factors based on UNSCALED data !')
c
c ... ADD M2 to M1
c
      else if (oper(1:3) .eq. 'ADD' .or.
     +         oper(1:3) .eq. 'GK+') then
c
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              map1(i1,i2,i3) = map1(i1,i2,i3) + 
     +                         map2(i1+io1,i2+io2,i3+io3)
            end do
          end do
        end do
c
c ... SUBTRACT M2 from M1
c
      else if (oper(1:3) .eq. 'GK-') then
c
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              map1(i1,i2,i3) = map1(i1,i2,i3) - 
     +                         map2(i1+io1,i2+io2,i3+io3)
            end do
          end do
        end do
c
c ... MULTIPLY M2 and M1
c
      else if (oper(1:3) .eq. 'GK*') then
c
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              map1(i1,i2,i3) = map1(i1,i2,i3) * 
     +                         map2(i1+io1,i2+io2,i3+io3)
            end do
          end do
        end do
c
c ... M1 = MIN(M1,M2) & M2 = MAX(M1,M2)
c
      else if (oper(1:3) .eq. 'MIN') then
c
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              if (map1(i1,i2,i3) .gt.
     +            map2(i1+io1,i2+io2,i3+io3)) then
                q1 = map1(i1,i2,i3)
                map1(i1,i2,i3) = 
     +                    map2(i1+io1,i2+io2,i3+io3)
                map2(i1+io1,i2+io2,i3+io3) = q1
              end if
            end do
          end do
        end do
c
c ... CORRELATE M1 & M2
c     GAMMA = 2*M1*M2/(M1**2+M2**2)
c     M1 = GAMMA*M1
c     M2 = GAMMA*M2
c
      else if (oper(1:3) .eq. 'COR') then
c
        do i1=lim1(1,1),lim1(2,1)
          do i2=lim1(1,2),lim1(2,2)
            do i3=lim1(1,3),lim1(2,3)
              q1 = map1(i1,i2,i3)
              q2 = map2(i1+io1,i2+io2,i3+io3)
              gamma = (q1*q1 + q2*q2)
              if (gamma .ne. 0.0) then
                gamma = 2.0*q1*q2 / gamma
                map1(i1,i2,i3) = q1*gamma
                map2(i1+io1,i2+io2,i3+io3) = q2*gamma
              end if
            end do
          end do
        end do
c
c ... INVALID OPERAND
c
      else
        call errcon ('Invalid operation :'//oper)
      end if
c
      return
      end
c
c
c
      subroutine mapmat (map,npnt,prod,plus)
c
      implicit none
c
      integer npnt,i
c
      real map(npnt),prod,plus
c
code ...
c
      do i=1,npnt
        map(i) = prod*map(i) + plus
      end do
c
      return
      end
c
c
c
      subroutine proplu (x1,x2,prod,plus)
c
      implicit none
c
      real xlo,xhi,prod,xplus,x1,x2
c
      integer plus
c
code ...
c
      xlo = x1
      xhi = x2
      call rlohi (xlo,xhi)
      if (xlo .eq. xhi) then
        call prompt (
     +    ' WARNING --- All map points have identical values !')
        if (xlo .eq. 0.0) then
          prod = 1.0
        else
          prod = 1.0 / xlo
        end if
        plus = 0
        goto 100
      end if
      prod = 255.0 / (xhi-xlo)
      if (abs(prod). le. 1.0e-2) prod = 0.01
      xplus = -xlo*prod
      if (prod .gt. 327.66) then
        write (*,*)
        call errcon ('PROD > 327.66; reset prod/plus yourself !')
        call prompt (' Or multiply your map with 10 or 100 or ...')
        call prompt (' This error ONLY affects the MAppage command,')
        call prompt (' not the BRix command !')
        write (*,*)
        xplus = xplus*(327.66/prod)
        prod = 327.66
      end if
      plus = nint (xplus)
      plus = max (-32765, min (plus, 32766))
c
c ... list current values
c
  100 continue
      write (*,'(1x,a,1p,2e12.4)')
     +  'Requested dynamic range :',xlo,xhi
      write (*,'(1x,a,1pe12.4,i8)')
     +  'Value of Prod and Plus  :',prod,plus
      call dynran (prod,plus,xlo,xhi)
      write (*,'(1x,a,1p,2e12.4)')
     +  'Actual dynamic range    :',xlo,xhi
c
      return
      end
c
c
c
      subroutine dynran (prod,plus,xlo,xhi)
c
      implicit none
c
      real prod,xplus,xlo,xhi
c
      integer plus
c
code ...
c
      xplus = plus
      xlo = -xplus/prod
      xhi = (255.0-xplus)/prod
      call rlohi (xlo,xhi)
c
      return
      end
c
c
c
      subroutine storem (map,ix,iy,iz,ilev,rho)
c
c ... store a plane the way Alwyn's routines expect them
c
      implicit none
c
      integer ix,iy,iz,ilev,i1,i2
c
      real map(ix,iy,iz),rho(iy,ix)
c
code ...
c
      do i1=1,ix
        do i2=1,iy
          rho (i2,i1) = map (i1,i2,ilev)
        end do
      end do
c
      return
      end
c
c
c
      subroutine paged (rholev, ix, iy, ilev, prod, plus, slice,
     +                  maxrho,iunit,indxra,nx,ny,nz,ierr)
c
c --- A new level of density, store it and if necessary
c     write it out as 3-d non-overlapping boxes of 8*8*8 values
c
      implicit none
c
      integer ix, iy, maxrho
      real rholev(iy,ix)
      integer slice (8,maxrho)
c
      real prod
c
      integer plus,iunit,indxra,ierr,nx,ny,nz
c
      integer i, ict, i1, ilev, j, jct, j1, j2, j3
      integer k, k1, k2, k3, value
c
      integer*2 irecrd(256)
      byte record(512)
      equivalence (record, irecrd)
c
code ...
c
      ierr = 0
c
      i1 = mod (ilev,8)
      if (i1 .eq. 0) i1 = 8
      ict = 0
      do i=1,iy
        do j=1,ix
          ict = ict+1
          value = rholev(i,j) * prod + plus
          if (value .gt. 255) value = 255
          if (value .lt.   0) value =   0
          slice (i1,ict) = value
        end do
      end do
c
c ---	Now pick out our non-overlapping bricks ?
c
      if (i1 .ne. 8) return
c
      value = 0
      ierr = -1
c
c ---	Loop over possible y-pages
c
      do 110 j=1,ny
        j1 = (j-1)*8+1
        j2 =   j  *8
c
c ---     Loop over possible x-pages
c
        do 110 k=1,nx
c
c ---	  Now get our loop parameters
c
          k1 = (k-1)*8+1
          k2 =   k  *8
          ict = 0
c ---	    Loop over z-levels
          do 120 i=1,8
c ---	    Loop over y-indeces of current page
            do 120 j3=j1,j2
              jct = (j3-1)*ix+k1-1
c ---	      Loop over x-indeces of current page
              do 120 k3=k1,k2
                ict = ict+1
                jct = jct+1
c ---	        If either direction over edge,pack record
                if(j3 .gt. iy  .or.  k3 .gt. ix)  then
                  call makbyt (record(ict), value)
                else
                  call makbyt (record(ict), slice(i,jct))
                end if
120       continue
          indxra = indxra+1
          call bytswp (record)
          write (iunit, rec=indxra, err=1300) irecrd
110   continue
c
      ierr = 0
      return
c
 1300 return
      end
c
c
c
      subroutine rest (ix, iy, ilev, prod, plus, slice,
     +                 maxrho,iunit,indxra,nx,ny,nz,ierr)
c
c ---	Write out rest of density
c ---	write it out as 3-d elements
c
      implicit none
c
      integer ix, iy, ilev, maxrho
      integer slice(8,maxrho)
c
      real prod
c
      integer plus,indxra,ierr,iunit,nx,ny,nz
c
      integer i, ict, i1, j, jct, j1, j2, j3 
      integer k, k1, k2, k3
c
      byte record(512)
      integer*2 irecrd(256)
      equivalence (record(1),irecrd(1))
c
code ...
c
      ierr = -1
c
      i1 = mod(ilev,8)
      if(i1 .eq. 0)i1= 8
      ict = 0
c
c ---	Now pick out our overlapping bricks
c ---	Loop over possible y-pages
c
      do 100 j=1,ny
        j1 = (j-1)*8+1
        j2 =   j  *8
c ---     Loop over possible x-pages
        do 100 k=1,nx
c ---	   Now get our loop parameters
          k1 = (k-1)*8+1
          k2 =   k  *8
          ict = 0
c ---	     Loop over z-levels
          do 110 i=1,8
c ---	     Loop over y-indeces of current page
            do 110 j3 = j1,j2
              jct = (j3-1)*ix+k1-1
c ---	       Loop over x-indeces of current page
              do 110 k3 = k1,k2
                ict = ict+1
                jct = jct+1
c ---	         If either direction over edge,pack record
c	         or if z - direction over edge,pack record
                if (j3 .gt. iy  .or.  k3 .gt. ix .or.
     $              i .gt. i1 ) then
                  call makbyt (record(ict), plus)
                else
                  call makbyt (record(ict), slice(i,jct))
                end if
110       continue
          indxra = indxra+1
          call bytswp (record)
          write (iunit, rec=indxra,err=1300) irecrd
100   continue
c
      ierr = 0
c
      return
c
 1300 return
      end
c
c
c
      logical function litend ()
C
C Check endedness, return true if little-ended (Vax-like),
C false if big-ended (IBM-like)
C
      integer i
c
      byte b(4)
c
      equivalence (i,b(1))
c
code ...
c
      i=1
      litend = (b(1) .ne. 0)
c
      return
      end
c
c
c
      subroutine bytswp (rec)
c
      implicit none
c
      byte rec(2,256), one
c
      integer i
c
code ...
c
      do i=1,256
        one = rec (1,i)
        rec(1,i) = rec(2,i)
        rec(2,i) = one
      end do
c
      return
      end
c
c
c
      subroutine makbyt (value, ivalue)
c
c ---	This is a machine specific routine.
c	Depending on whether the machine is big or little endian
c	the last assignment statement will have to be changed.
c ---	On Vax,IBM  PC's need one assignment
c	Other machines use the other
c
      implicit none
c
      byte value
      integer ivalue
c
c      integer i
c      byte j(4)
c      equivalence (i,j)
c
code ...
c
C+PRE This should work on all machines
      value = ivalue
c
c	i = ivalue
c ---	Most unix machines and IBM mainframes
c	value = j(4)
c ---	VAX/IBM PC
c	value = j(1)
c
	return
	end
c
c
c
      subroutine zeromp (map,npnt,xlo,xhi)
c
      implicit none
c
      integer npnt,i,n1,n2
c
      real map(npnt),xlo,xhi
c
code ...
c
      n1 = 0
      n2 = 0
      do i=1,npnt
        if (map(i) .lt. xlo) then
          map(i) = 0.0
          n1 = n1 + 1
        else if (map(i) .gt. xhi) then
          map(i) = 0.0
          n2 = n2 + 1
        end if
      end do
c
      call jvalut (' Less    than lower bound :',1,n1)
      call jvalut (' Greater than upper bound :',1,n2)
c
      return
      end
c
c
c
      subroutine chaval (map,npnt,xlo,xhi,xnew)
c
      implicit none
c
      integer npnt,i,n1
c
      real map(npnt),xlo,xhi,xnew
c
code ...
c
      n1 = 0
      do i=1,npnt
        if (map(i) .gt. xlo .and. map(i) .lt. xhi) then
          map(i) = xnew
          n1 = n1 + 1
        end if
      end do
c
      call jvalut (' Nr of points changed :',1,n1)
c
      return
      end
c
c
c
      subroutine pageb (olun, rhosec, slices, ix, iy, nx, ny,
     +                  ilev, icnt, prod, plus, err)
c
c ... Morten's BRIX routine
c     950512 - hacked by GJK to fit in with MAPMAN
c
      implicit none
c
      integer olun, ix, iy, nx, ny, err, icnt, plus
      real rhosec(iy, ix), prod
      integer slices(8,ix*iy)
c
c A new level of density, store it and if necessary, write it out as 
c 3-d non-overlapping boxes of 8*8*8 values
c Original logic by Alwyn Jones, a long time ago.
c Modified as library routine, new brick format, Morten Kjeldgaard, Nov 1993
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer i, ict, i1, ilev, j, jct, j1, j2, j3, k, k1, k2, k3, value
      character*1 record(512) 
c
code ...
c
      err = 0
c
      i1 = mod (ilev,8)
      if (i1 .eq. 0) i1 = 8
      ict = 0

      if (ilev .ne. 0) then
         do i=1,iy
            do j=1,ix
              ict = ict+1
              value = rhosec(i,j) * prod + plus
              if (value .gt. 255) value = 255
              if (value .lt.   0) value =   0
              slices (i1,ict) = value
            enddo
         enddo
      endif
c
c Pick out our non-overlapping bricks?
      if (i1 .ne. 8 .and. ilev .ne. 0) return
c
      value = 0
c     Loop over possible y-pages
      do j=1,ny
         j1 = (j-1)*8+1
         j2 =   j  *8
c
c        Loop over possible x-pages
         do k=1,nx
	    k1 = (k-1)*8+1
	    k2 =   k  *8
	    ict = 0
c
c           Loop over z-levels
	    do i=1,8
c
c              Loop over y-indeces of current page
               do j3=j1,j2
                  jct = (j3-1)*ix+k1-1
c
c                 Loop over x-indeces of current page
                  do k3=k1,k2
                     ict = ict+1
                     jct = jct+1
c
c                    If either direction over edge, pack record
                     if (j3 .gt. iy  .or.  k3 .gt. ix 
     $                    .or. i .gt. i1)  then
                        record(ict) = char(0)
                     else
ccc                        record(ict) = char(int(slices(i,jct)))
                        record(ict) = char(slices(i,jct))
                     end if
                  enddo
               enddo
            enddo
c
            icnt = icnt + 1
            write (olun, rec=icnt, iostat=err) record
c
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine getmst (map1,na1,nb1,nc1,
     +                   map2,na2,nb2,nc2,lim1,lim2,
     +                   cutoff,buffer,maxpnt)
c
c ... combine maps
c
      implicit none
c
      integer na1,nb1,nc1,na2,nb2,nc2,lim1(2,3),lim2(2,3)
      integer i1,i2,i3,io1,io2,io3,ntot,maxpnt,ipass
c
      real map1(na1,nb1,nc1),map2(na2,nb2,nc2)
      real buffer(maxpnt)
c
      real cutoff,q1,q2,xave,xsdv,xmin,xmax,xtot
c
      logical lmask
c
code ...
c
      io1 = lim2(1,1) - lim1(1,1)
      io2 = lim2(1,2) - lim1(1,2)
      io3 = lim2(1,3) - lim1(1,3)
c
c ... IPASS = 1 -> inside mask
c     IPASS = 2 -> outside mask
c
      do ipass=1,2
        ntot = 0
        do i1=1,na1
          do i2=1,nb1
            do i3=1,nc1
              lmask = .false.
              q1 = map1(i1,i2,i3)
c
c ... inside mask grid ?
c
              if (i1 .lt. lim1(1,1)) goto 100
              if (i1 .gt. lim1(2,1)) goto 100
              if (i2 .lt. lim1(1,2)) goto 100
              if (i2 .gt. lim1(2,2)) goto 100
              if (i3 .lt. lim1(1,3)) goto 100
              if (i3 .gt. lim1(2,3)) goto 100
              q2 = map2(i1+io1,i2+io2,i3+io3)
              lmask = (q2 .ge. cutoff)
c
  100         continue
              if (ipass .eq. 1) then
                if (lmask) then
                  ntot = ntot + 1
                  buffer (ntot) = q1
                end if
              else if (ipass .eq. 2) then
                if (.not. lmask) then
                  ntot = ntot + 1
                  buffer (ntot) = q1
                end if
              end if
c
            end do
          end do
        end do
c
        if (ipass .eq. 1) then
          call prompt ('0Density inside the mask :')
        else
          call prompt ('0Density outside the mask :')
        end if
c
        call jvalut (' Nr of points :',1,ntot)
c
        if (ntot .gt. 0) then
          call xstats (buffer,ntot,xave,xsdv,xmin,xmax,xtot)
          call rvalut (' Average density :',1,xave)
          call rvalut (' St. deviation   :',1,xsdv)
          call rvalut (' Minimum value   :',1,xmin)
          call rvalut (' Maximum value   :',1,xmax)
          call rvalut (' Sum of values   :',1,xtot)
        else
          call errcon ('No points found !')
        end if
c
      end do
c
      return
      end
c
c
c
      subroutine peekit (map,grid,origin,ext1,ext2,ext3,cell,
     +  option,iunit,junit,mymode,npt,rad)
c
      implicit none
c
      integer ext1,ext2,ext3
      real map(ext1,ext2,ext3)
c
      integer grid(3),origin(3)
      integer iunit,junit,npt,nlines,natoms,i,off,ierr,j,k,l
      integer ii,jj,kk,leng1,i1,i2,i3,nsum,nout
c
      real cell(6),rad,a(3,3),b(3,3),c(3,3),g(3),x(3),fake(6)
      real x1(3),x2(3),xp(3),sum,value,bnew,distce
c
      character option*2,mymode*1
      character line*128,atmnam*14
c
code ...
c
      nlines = 0
      natoms = 0
      nout = 0
c
c ... foreplay
c
      do i=1,3
        fake(i) = 1.0
        fake(i+3) = cell(i+3)
      end do
      call orthog (fake, a, 1)
      call orthog (fake, c, 0)
      call orthog (cell, b, 1)
      do i=1,3
        g(i) = cell(i)/float(grid(i))
      end do
      off = 1 + int (0.1 + (rad/(min(g(1),g(2),g(3)))))
      call fvalut (' Grid spacing (A) :',3,g)
      call prompt (' Peeking ...')
c
   10 continue
      read (iunit,'(a)',end=100,err=9999) line
      nlines = nlines + 1
      if (line(1:6) .ne. 'ATOM  ' .and.
     +    line(1:6) .ne. 'HETATM') then
        if (junit.gt.0) write (junit,'(a)') line(1:leng1(line))
        goto 10
      end if
      natoms = natoms + 1
      bnew = -9999.999
c
      call upcase (line)
      read (line,6000,err=9999) atmnam,(x(j),j=1,3)
c
 6000 format (12x,a14,4x,3f8.3)
c
      call mulmtx (a, x, x1, 3, 3, 1)
c
      x2(1) = (x1(1)/g(1)) - float(origin(1)) + 1.0
      x2(2) = (x1(2)/g(2)) - float(origin(2)) + 1.0
      x2(3) = (x1(3)/g(3)) - float(origin(3)) + 1.0
c
      i1 = nint(x2(1))
      i2 = nint(x2(2))
      i3 = nint(x2(3))
c
      if (i1 .le. 0) goto 20
      if (i2 .le. 0) goto 20
      if (i3 .le. 0) goto 20
c
      if (i1 .gt. ext1) goto 20
      if (i2 .gt. ext2) goto 20
      if (i3 .gt. ext3) goto 20
c
c ... option = VAlue
c
      if (option .eq. 'VA') then
        if (mymode .eq. 'N') then
          bnew = map (i1,i2,i3)
        else if (mymode .eq. 'S') then
          call intrpl (map, ext1, ext2, ext3, x2, value, ierr)
          if (ierr .ne. 0) then
            call errcon ('Spline interpolation error')
            call textut (' For atom :',atmnam)
          else
            bnew = value
          end if
        else if (mymode .eq. 'I') then
          call linint (map, ext1, ext2, ext3, x2, value, ierr)
          if (ierr .ne. 0) then
            call errcon ('Linear interpolation error')
            call textut (' For atom :',atmnam)
          else
            bnew = value
          end if
        end if
      else if (option .eq. 'CU') then
        sum = 0.0
        nsum = 0
        do ii=i1-npt,i1+npt
          do jj=i2-npt,i2+npt
            do kk=i3-npt,i3+npt
              if (ii .le. 0) goto 30
              if (jj .le. 0) goto 30
              if (kk .le. 0) goto 30
              if (ii .gt. ext1) goto 30
              if (jj .gt. ext2) goto 30
              if (kk .gt. ext3) goto 30
              nsum = nsum + 1
              if (mymode .eq. 'I') then
                sum = sum + map(ii,jj,kk)
              else if (mymode .eq. 'A') then
                sum = sum + abs(map(ii,jj,kk))
              else if (mymode .eq. 'R') then
                sum = sum + map(ii,jj,kk)*map(ii,jj,kk)
              else if (mymode .eq. 'M') then
                sum = sum + map(ii,jj,kk)
              end if
   30         continue
            end do
          end do
        end do
        if (nsum .eq. 0) then
          call errcon ('No points in cube !')
          call textut (' For atom :',atmnam)
        else
          if (mymode .eq. 'I') then
            bnew = sum
          else if (mymode .eq. 'A') then
            bnew = sum
          else if (mymode .eq. 'R') then
            bnew = sqrt ( sum/float(nsum) )
          else if (mymode .eq. 'M') then
            bnew = sum / float(nsum)
          end if
        end if
      else if (option .eq. 'SP') then
        sum = 0.0
        nsum = 0
        do j=-off,off
          xp(3) = x1(3)+ j*g(3)
          do k=-off,off
            xp(2) = x1(2)+ k*g(2)
            do l=-off,off
              xp(1) = x1(1)+ l*g(1)
              call mulmtx (c, xp, x2, 3, 3, 1)
              if (distce(x,x2) .gt. rad) goto 40
              ii = nint(xp(1)/g(1))- origin(1) +1
              jj = nint(xp(2)/g(2))- origin(2) +1
              kk = nint(xp(3)/g(3))- origin(3) +1
              if (ii .le. 0) goto 40
              if (jj .le. 0) goto 40
              if (kk .le. 0) goto 40
              if (ii .gt. ext1) goto 40
              if (jj .gt. ext2) goto 40
              if (kk .gt. ext3) goto 40
              nsum = nsum + 1
              if (mymode .eq. 'I') then
                sum = sum + map(ii,jj,kk)
              else if (mymode .eq. 'A') then
                sum = sum + abs(map(ii,jj,kk))
              else if (mymode .eq. 'R') then
                sum = sum + map(ii,jj,kk)*map(ii,jj,kk)
              else if (mymode .eq. 'M') then
                sum = sum + map(ii,jj,kk)
              end if
   40         continue
            end do
          end do
        end do
        if (nsum .eq. 0) then
          call errcon ('No points in sphere !')
          call textut (' For atom :',atmnam)
        else
          if (mymode .eq. 'I') then
            bnew = sum
          else if (mymode .eq. 'A') then
            bnew = sum
          else if (mymode .eq. 'R') then
            bnew = sqrt ( sum/float(nsum) )
          else if (mymode .eq. 'M') then
            bnew = sum / float(nsum)
          end if
        end if
      end if
c
      goto 23
c
c1234567890123456789012345678901234567890123456789012345678901234567890
cATOM     59  N   LYS     8      30.052  19.368  32.702  1.00  9.39      1CBS 281
c
   20 continue
      call textut (' Warning - atom outside map :',line(1:30))
      nout = nout + 1
c
   23 continue
      if (abs(bnew) .ge. 1000.0) then
        write (line(61:66),'(f6.0)') bnew
      else if (abs(bnew) .ge. 100.0) then
        write (line(61:66),'(f6.1)') bnew
      else if (abs(bnew) .ge. 10.0) then
        write (line(61:66),'(f6.2)') bnew
      else if (abs(bnew) .ge. 1.0) then
        write (line(61:66),'(f6.3)') bnew
      else if (abs(bnew) .ge. 0.1) then
        write (line(61:66),'(f6.4)') bnew
      else
        write (line(61:66),'(f6.4)') bnew
      end if
c
      write (junit,'(a)') line(1:leng1(line))
c
      write (*,'(1x,a,a,1pe12.4)') line(1:30),
     +  ' --> PEEK-A-BOO --> ',bnew
c
      goto 10
c
  100 continue
      ierr = 0
      call prompt (' New PDB file created')
      call jvalut (' Nr of atoms read and processed   :',1,natoms)
      call jvalut (' Nr outside the limits of the map :',1,nout)
c
      return
c
 9999 continue
      call errcon ('While reading PDB file')
      return
c
      end
c
c
c
      subroutine dominv (map1,map2,na,nb,nc)
c
c ... copy map2 to map1
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3
      real map1(na,nb,nc),map2(na,nb,nc)
c
code ...
c
      do i1=1,na
        do i2=1,nb
          do i3=1,nc
            map1 (i1,i2,i3) = map2 (na+1-i1,nb+1-i2,nc+1-i3)
          end do
        end do
      end do
c
      return
      end
c
c
c
      subroutine rem08 (map,iunit,nx,ny,nz,ierr)
c
c ... read EM08 map
c
      implicit none
c
      integer maxbuf
      parameter (maxbuf=10000)
c
      integer nx,ny,nz
      real map(nx,ny,nz)
c
      real*4 xdata(maxbuf)
c
      integer iunit,ierr,i,j,k
c
code ...
c
      ierr = -1
c
      if (nz .gt. maxbuf) then
        call errcon ('Records too long; change MAXBUF')
        return
      end if
c
c ... why won't it read ???
c
c      do i=1,10
c        read (iunit) xdata(i)
c      end do
c      rewind (iunit)
c
      do k=1,nx
        do j=1,ny
          read (iunit,end=9000,err=9000) (xdata(i),i=1,nz)
          do i=1,nz
            map (k,j,i) = xdata (i)
          end do
cccc     +      (map(k,j,i),i=1,nz)
cccc     +      ((map(k,j,i),i=1,nz),j=1,ny)
cccc     +      (((map(k,j,i),i=1,nz),j=1,ny),k=1,nx)
        end do
      end do
c
      ierr = 0
      return
c
 9000 continue
      ierr = -1
      call errcon ('While reading EM08 map')
c
      print *,'XYZ = ',k,j,i
      print *,'LEN = ',nx,ny,nz
c
      return
      end
c
c
c
      subroutine irough (mapn,mapo,masky,nx,ny,nz,npos,nneg,ibox,
     +  qdmin,qdmax,qdave,qdsdv)
c
      implicit none
c
      integer nx,ny,nz
      real mapn(nx,ny,nz),mapo(nx,ny,nz)
      integer masky(nx,ny,nz)
c
      real*8 dmin,dmax,dave,dsdv,rms,xdum,lrms
c
      real qdmin,qdmax,qdave,qdsdv,xmin,xmax
c
      integer npos,nneg,ibox
      integer nxyz,i1,i2,i3,np,nn,j1,j2,j3,npt,nmax,nl
c
code ...
c
      dmin = qdmin
      dmax = qdmax
c
      call prompt (' Initialising map and mask ...')
      nxyz = nx * ny * nz
      call inimsk (masky,nx,ny,nz,0)
      call inimap (mapn,nx,ny,nz,0.0)
c
      do i1=1+ibox,nx-ibox
        do i2=1+ibox,ny-ibox
          do i3=1+ibox,nz-ibox
            masky (i1,i2,i3) = 1
          end do
        end do
      end do
c
      if (npos .gt. 0) then
        call prompt (' Removing highest peaks ...')
        np = 0
  100   continue
        xmax = qdmin
        j1 = 0
        j2 = 0
        j3 = 0
        do i1=1+ibox,nx-ibox
          do i2=1+ibox,ny-ibox
            do i3=1+ibox,nz-ibox
              if (masky(i1,i2,i3) .ne. 0) then
                if (mapo(i1,i2,i3) .gt. xmax) then
                  xmax = mapo(i1,i2,i3)
                  j1 = i1
                  j2 = i2
                  j3 = i3
                end if
              end if
            end do
          end do
        end do
c
        call rvalut (' Peak value :',1,xmax)
        np = np + 1
        do i1=j1-ibox,j1+ibox
          if (i1 .lt. 1 .or. i1 .gt. nx) goto 110
          do i2=j2-ibox,j2+ibox
            if (i2 .lt. 1 .or. i2 .gt. ny) goto 112
            do i3=j3-ibox,j3+ibox
              if (i3 .lt. 1 .or. i3 .gt. nz) goto 114
              masky (i1,i2,i3) = 0
  114         continue
            end do
  112       continue
          end do
  110     continue
        end do
c
        if (np .lt. npos) goto 100
      end if
c
      if (nneg .gt. 0) then
        call prompt (' Removing lowest peaks ...')
        nn = 0
  200   continue
        xmin = qdmax
        j1 = 0
        j2 = 0
        j3 = 0
        do i1=1+ibox,nx-ibox
          do i2=1+ibox,ny-ibox
            do i3=1+ibox,nz-ibox
              if (masky(i1,i2,i3) .ne. 0) then
                if (mapo(i1,i2,i3) .lt. xmin) then
                  xmin = mapo(i1,i2,i3)
                  j1 = i1
                  j2 = i2
                  j3 = i3
                end if
              end if
            end do
          end do
        end do
c
        call rvalut (' Peak value :',1,xmin)
        nn = nn + 1
        do i1=j1-ibox,j1+ibox
          if (i1 .lt. 1 .or. i1 .gt. nx) goto 210
          do i2=j2-ibox,j2+ibox
            if (i2 .lt. 1 .or. i2 .gt. ny) goto 212
            do i3=j3-ibox,j3+ibox
              if (i3 .lt. 1 .or. i3 .gt. nz) goto 214
              masky (i1,i2,i3) = 0
  214         continue
            end do
  212       continue
          end do
  210     continue
        end do
c
        if (nn .lt. nneg) goto 200
      end if
c
      call prompt (' Calculating local RMS density ...')
      npt = 0
      rms = 0.0
      nmax = (2*ibox+1)*(2*ibox+1)*(2*ibox+1)
      call ivalut (' Box size :',1,nmax)
c
      dmin = 9999.0
      dmax = -9999.0
      dave = 0.0
      dsdv = 0.0
c
      do i1=1+ibox,nx-ibox
        do i2=1+ibox,ny-ibox
          do i3=1+ibox,nz-ibox
            if (masky(i1,i2,i3) .ne. 0) then
              nl = 0
              lrms = 0.0
              do j1=i1-ibox,i1+ibox
                do j2=i2-ibox,i2+ibox
                  do j3=i3-ibox,i3+ibox
                    if (masky(j1,j2,j3) .ne. 0) then
                      nl = nl + 1
                      lrms = lrms + mapo(j1,j2,j3)*mapo(j1,j2,j3)
                    end if
                  end do
                end do
              end do
              if ( (2*nl) .lt. nmax ) then
                masky (i1,i2,i3) = 0
              else
                npt = npt + 1
                lrms = dsqrt (lrms/dble(nl))
                mapn (i1,i2,i3) = lrms
                rms = rms + lrms*lrms
                if (lrms .lt. dmin) dmin = lrms
                if (lrms .gt. dmax) dmax = lrms
                dave = dave + lrms
              end if
            end if
          end do
        end do
      end do
c
      call prompt (' Finishing touch ...')
      dsdv = rms
      rms = sqrt (rms/dble(npt))
      lrms = 1.0D0 / rms
      call mapmat (mapn,nxyz,real(lrms),0.0)
c
      dmin = dmin * lrms
      dmax = dmax * lrms
      dave = dave * lrms / dble(npt)
      dsdv = dsdv * lrms * lrms
      xdum = (dsdv/dble(npt)) - (dave*dave)
      dsdv = 0.0D0
      if (xdum .ge. 0.0D0) dsdv = dsqrt(xdum)
c
      call jvalut (' Nr of map points used :',1,npt)
      call rvalut (' RMS (local RMS dens)  :',1,real(rms))
      call prompt (' Statistics after dividing by RMS(RMS) :')
      call rvalut (' Minimum :',1,real(dmin))
      call rvalut (' Maximum :',1,real(dmax))
      call rvalut (' Average :',1,real(dave))
      call rvalut (' St.dev. :',1,real(dsdv))
c
      if (dmin .gt. 0.0D0) dmin = 0.0D0
c
      qdmin = dmin
      qdmax = dmax
      qdave = dave
      qdsdv = dsdv
c
      return
      end
c
c
c
      subroutine wherev (map,ori,ext1,ext2,ext3,val,tol)
c
      implicit none
c
      integer ext1,ext2,ext3,ori(3),i,j,k,n,ntot
c
      real map(ext1,ext2,ext3),val,tol,dlo,dhi,perc
c
code ...
c
      dlo = val - tol
      dhi = val + tol
      call rlohi (dlo,dhi)
      n = 0
      ntot = ext1*ext2*ext3
c
      do i=1,ext1
        do j=1,ext2
          do k=1,ext3
            if (map(i,j,k).ge.dlo) then
              if (map(i,j,k).le.dhi) then
                write (*,6000) ori(1)+i-1,ori(2)+j-1,
     +            ori(3)+k-1,map(i,j,k)
                n = n + 1
              end if
            end if
          end do
        end do
      end do
c
 6000 format (' Value @ X, Y, Z ',3i6,' = ',1pe12.4)
c
      call jvalut (' Nr of points listed :',1,n)
      call jvalut (' Nr of points in map :',1,ntot)
      perc = 100.0 * float(n) / float(ntot)
      call fvalut (' Percentage listed   :',1,perc)
c
      return
      end
c
c
c
      subroutine vrdots (map,na,nb,nc,grid,origin,cell,
     +                   vlo,vhi,iunit,ierr)
c
      implicit none
c
      integer na,nb,nc,grid(3),origin(3),iunit,ierr
c
      real map(na,nb,nc)
      real cell(6),a(3,3),x(3),y(3),vlo,vhi
c
      integer i,j,k,l,leng1,nat
c
      character line*80
c
code ...
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  PointSet { startIndex 0 '/
     +  '    numPoints ',i8,' } }')
c
      call orthog (cell,a,0)
c
      write (iunit,6000,err=9999)
c
      nat = 0
      do i=1,na
        x(1) = float(origin(1)+i-1) *
     +         (cell(1)/float(grid(1))) / cell(1)
        do j=1,nb
          x(2) = float(origin(2)+j-1) *
     +           (cell(2)/float(grid(2))) / cell(2)
          do k=1,nc
            if (map(i,j,k) .ge. vlo) then
              if (map(i,j,k) .le. vhi) then
                x(3) = float(origin(3)+k-1) *
     +                 (cell(3)/float(grid(3))) / cell(3)
                call mulmtx (a,x,y,3,3,1)
                write (line,6010,err=9999) (y(l),l=1,3)
                call pretty (line)
                write (iunit,'(a)',err=9999) line(1:leng1(line))
                nat = nat + 1
              end if
            end if
          end do
        end do
      end do
c
      write (iunit,6020,err=9999)
      write (iunit,6030,err=9999) nat
c
      call jvalut (' Nr of points written :',1,nat)
      ierr = 0
      return
c
 9999 continue
      call errcon ('While writing VRML file')
      ierr = -1
      return
c
      end
c
c
c
      subroutine vrdraw (map,na,nb,nc,grid,origin,cell,
     +                   vlevel,iunit,buffer,ierr)
c
      implicit none
c
      integer maxbuf
      parameter (maxbuf=200000)
c
      integer na,nb,nc,grid(3),origin(3),ext(3),iunit,ierr
c
      real map(na,nb,nc)
      real buffer(*)
      real cell(6),a(3,3),vlevel,rpl,zc(6)
c
      integer i,i1,i2,ipl,itype,ivrml
      integer zo(3),zg(3),nlev,icol,npt,ncont
c
      character plan(3)*1
      character bits(2*maxbuf)*1
c
      common /LPOST/ a,zc,zo,zg,rpl,ipl,i1,i2,itype,ivrml,npt,ncont
c
code ...
c
      nlev = 1
      icol = 1
c
      call orthog (cell,a,0)
c
      do i=1,3
        zc(i) = cell(i)
        zc(i+3) = cell(i+3)
        zo(i) = origin(i)
        zg(i) = grid(i)
      end do
c
      ext(1) = na
      ext(2) = nb
      ext(3) = nc
c
      plan(1) = 'X'
      plan(2) = 'Y'
      plan(3) = 'Z'
c
      ivrml = iunit
c
      do itype=1,3
c
        ncont = 0
        npt = 0
c
        call textut (
     +    ' Contour planes perpendicular to :',plan(itype))
        call ivalut (' Number of planes :',1,ext(itype))
c
ccc        call rvalut (' Contour level :',1,vlevel)
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
        if (ext(i1)*ext(i2) .gt. maxbuf) then
          call errcon ('CONTPL - Buffer too small; tell Gerard !')
          goto 1000
        end if
c
ccc        print *,'ITYPE,I1,I2 = ',itype,i1,i2
c
        do ipl=1,ext(itype)
c
          rpl = float(origin(itype)+ipl-1) *
     +          (cell(itype)/float(grid(itype))) / cell(itype)
c
ccc          print *,'PLANE ',ipl,rpl
c
          call cutpln (itype,ipl,map,na,nb,nc,buffer,
     +               ext(i1),ext(i2))
          call hpcntr (buffer,bits,ext(i1),ext(i2),vlevel,icol,nlev)
        end do
c
 1000   continue
        call jvalut (' Nr of contours drawn  :',1,ncont)
        call jvalut (' Nr of points included :',1,npt)
c
      end do
c
      call prompt (' Contours written')
      ierr = 0
c
      return
      end
c
c
c
      subroutine aplot(x,y,ipt)
c     =========================
c
c ... RB's magic bit
c
c ... join IPT coordinate pairs X,Y of contour crossings, in grid units
c     for m*n grid , coordinates are in range 0,0 to m-1,n-1
c
      implicit none
c
      integer ipt,i,ipl,i1,i2,itype,iunit,imax,grid(3),origin(3),j
      integer leng1,npt,ncont
c
      real x(ipt),y(ipt),x1,y1,xl,yl,a(3,3),xx(3),yy(3),rpl,cell(6)
c
      character line*256
c
      common /LPOST/ a,cell,origin,grid,rpl,ipl,i1,i2,itype,
     +               iunit,npt,ncont
c
code ...
c
c ... save 1st and last points before transformation
c
      x1=x(1)
      y1=y(1)
      xl=x(ipt)
      yl=y(ipt)
c
      ncont = ncont + 1
      npt = npt + ipt
c
c ... work out coordinates in Cartesian space
c
      write (iunit,6000,err=9999)
      do i=1,ipt
        xx(itype) = rpl
        xx(i1) = (origin(i1)+x(i)) *
     +          (cell(i1)/float(grid(i1))) / cell(i1)
        xx(i2) = (origin(i2)+y(i)) *
     +          (cell(i2)/float(grid(i2))) / cell(i2)
        call mulmtx (a,xx,yy,3,3,1)
        write (line,6010) yy(1),yy(2),yy(3)
        call pretty (line)
        write (iunit,'(a)',err=9999) line(1:leng1(line))
      end do
      write (iunit,6020,err=9999)
      write (iunit,6030,err=9999)
      do i=1,ipt,20
        imax = min(ipt,i+19)
        write (line,6040) (j-1,j=i,imax)
        call remspa (line)
        write (iunit,'(a)',err=9999) line(1:leng1(line))
      end do
      write (iunit,6050,err=9999)
c
c ... restore 1st & last points
c
      x(1)=x1
      y(1)=y1
      x(ipt)=xl
      y(ipt)=yl
c
      return
c
 6000 format ('Separator { Coordinate3 { point [')
 6010 format (1x,3(f8.3,1x),',')
 6020 format (' ] }')
 6030 format ('  IndexedLineSet { coordIndex [ ')
 6040 format (20(i8,','))
 6050 format ('-1 ] } }')
c
 9999 continue
      call errcon ('While writing VRML file')
      call xvrml_close ()
c
      return
      end
c
c
c
      subroutine acolor (icol)
c     ========================
c
c ... set colour index
c
      implicit none
c
      integer icol
c
code ...
c
c ... do nothing for VRML
c
      return
      end
c
c
c
      subroutine cutpln (itype,ipl,mapje,ext1,ext2,ext3,
     +                   buffer,jext1,jext2)
c
      implicit none
c
      integer ext1,ext2,ext3,jext1,jext2
c
      real mapje(ext1,ext2,ext3)
      real buffer(jext1,jext2)
c
      integer itype,ipl,j1,j2
c
code ...
c
      if (itype .eq. 1) then
        do j1=1,ext2
          do j2=1,ext3
            buffer(j1,j2) = mapje (ipl,j1,j2)
          end do
        end do
      else if (itype .eq. 2) then
        do j1=1,ext1
          do j2=1,ext3
            buffer(j1,j2) = mapje (j1,ipl,j2)
          end do
        end do
      else
        do j1=1,ext1
          do j2=1,ext2
            buffer(j1,j2) = mapje (j1,j2,ipl)
          end do
        end do
      end if
c
      return
      end
C
C
C
      SUBROUTINE hpcntr (A,BITS,M,N,CONT,LCOL,NCONT)
C     ==============================================
C
      PARAMETER(IPTMAX=1024)
C
C***************************************
C
C                  contour routine
C
C***************************************
C
C  a is function to be contoured, dimensioned a(m,n) or a(m*n)
C  m = number of rows of a
C  n = number of columns of a
C  bits = logical*1 array dimensioned at least 2*m*n in calling routine
C  cont = array of ncont contour levels
C  lcol = array of ncont contour level colours
c  ncont = number of contour levels
c
C  x,y  arrays of dimension iptmax used in subroutine to store contour l
C
C  idiv = number of line segments between contour level crossings in uni
C     cell.  if >= 1, then routine finds idiv-1 coordinates between cros
C     by linear interpolation on the unit cell.
C  contour output via calls to subroutine aplot (x,y,ipt)
C     x,y are coordinates with 0<=x<m-1 and 0<=y<n-1.
C  subroutine constl(cl) is called to set up plotting style for contour
C
C
C
      DIMENSION IDIR(4), ISIDE(5), IVCT(4), JVCT(4), KVCT(4), LVCT(4),
     *SIDE(4), NVCT(4)
      CHARACTER*1 BITS(*)
      DIMENSION A(*),C(4),CONT(*)
      INTEGER LCOL(*)
      DIMENSION X(IPTMAX),Y(IPTMAX)
      DATA  IDIR/3,4,1,2/, ISIDE/1,2,3,4,1/, IVCT/0,-1,0,1/,
     *JVCT/-1,0,1,0/, SIDE/0.,0.,1.0,1.0/, KVCT/4,1,2,3/, LVCT/2,3,4,1/
C
Cc
Cc
      NX=M-1
      NY=N-1
      NCOUNT=0
      CLOLD=1.E28
      IDIV=1
      MN=M*N
      MN2=2*MN
      DIV=IDIV
C
C loop contour levels cl
      DO 190 ICONT=1,NCONT
      CL=CONT(ICONT)
C
C set colour index for this level
      CALL ACOLOR(LCOL(ICONT))
C
C set up plotting style for this contour level
C     call constl(cl)
C
      DO 10 I=1,MN2
   10 BITS(I)='f'
      NVCT(1)=0
      NVCT(2)=MN
      NVCT(3)=M
      NVCT(4)=MN+1
      IPT=1
      MM=M-1
      NN=N-1
C     search for contour crossing between adjacent column of array a(i,j
      I=0
      J=1
      ISUB=0
      JSUB=0
      IRTN=1
  100 IF (J .GT. N) GO TO 140
  110 IF (I .GE. MM) GO TO 130
      I=I+1
      ISUB=ISUB+1
      JSUB=JSUB+1
      IF (A(ISUB)-CL) 115,600,120
  115 IF (A(ISUB+1)-CL) 110,110,125
  120 IF (A(ISUB+1)-CL) 125,110,110
  125 IF (BITS(JSUB+NVCT(1)).EQ.'t')  GO TO 110
      XSTART=(CL-A(ISUB))/(A(ISUB+1)-A(ISUB))
      YSTART=0
      GO TO 200
  130 I=0
      ISUB=ISUB+1
      JSUB=JSUB+1
      J=J+1
      GO TO 100
C     search for contour crossing between adjacent rows of array a(i,j)
  140 I=0
      J=1
      JSUB=0
      ISUB=0
      IRTN=2
  150 IF (J .GT. NN) GO TO 190
  160 IF (I .GE. M) GO TO 180
      I=I+1
      ISUB=ISUB+1
      JSUB=JSUB+1
      IF (A(ISUB)-CL) 165,160,170
  165 IF (A(ISUB+M)-CL) 160,160,175
  170 IF (A(ISUB+M)-CL) 175,160,160
  175 IF (BITS(JSUB+NVCT(2)).EQ.'t')  GO TO 160
      YSTART=(CL-A(ISUB))/(A(ISUB+M)-A(ISUB))
      XSTART=0
      GO TO 200
  180 I=0
      J=J+1
      GO TO 150
C
C     begin following contour line... save indices for return to search
  200 ISAVE=I
      JSAVE=J
      ISUBSV=ISUB
      JSUBSV=JSUB
      XSAVE=XSTART
      YSAVE=YSTART
      X(1)=XSTART+FLOAT(I-1)
      Y(1)=YSTART+FLOAT(J-1)
      IENT=IRTN
      IRS=0
      GO TO 250
C     dump line and follow contour line on opposite side of starting pio
C     when used a second time this entry returns to search
  205 IRS=1
  210 IF (IPT .GT. 1) CALL  APLOT (X,Y,IPT)
      IPT=1
      I=ISAVE
      J=JSAVE
      ISUB=ISUBSV
      JSUB=JSUBSV
      XSTART=XSAVE
      YSTART=YSAVE
      X(1)=XSTART+FLOAT(I-1)
      Y(1)=YSTART+FLOAT(J-1)
      IF (IRS.NE.0) GO TO (110,160), IRTN
      IEXIT=IRTN
      IRS=1
      GO TO 240
C     return from following contour line through a cell
  230 IF (BITS(JSUB+NVCT(IEXIT)).EQ.'t')  GO TO 205
  240 I=I+IVCT(IEXIT)
      J=J+JVCT(IEXIT)
      JSUB=I+(J-1)*M
      ISUB=JSUB
      IENT=IDIR(IEXIT)
  250 BITS(JSUB+NVCT(IENT))='t'
      IF (I.LT.1 .OR. I.GT.MM .OR. J.LT.1 .OR. J.GT.NN)  GO TO 210
C     find contour crossing in new cell
  260 IF (ISUB+1.GT.MN .OR. ISUB+M.GT.MN
     1     .OR. ISUB+1+M.GT.MN)  GO TO 210
      C(1)=A(ISUB+1)
      C(2)=A(ISUB)
      C(3)=A(ISUB+M)
      C(4)=A(ISUB+1+M)
      JRTN=1
      ICNT=1
      JCNT=1
      DO 290 IROUND=1,4
      IF (IROUND .EQ. IENT) GO TO 290
      I1=ISIDE(IROUND)
      I2=ISIDE(IROUND+1)
      IF (C(I1)-CL) 270,285,275
  270 IF (C(I2)-CL) 290,290,280
  275 IF (C(I2)-CL) 280,290,290
  280 IEXIT=IROUND
      ICNT=ICNT+1
      GO TO 290
  285 JEXIT=IROUND
      JCNT=JCNT+1
  290 CONTINUE
      GO TO (300,310,700,210), JCNT
  300 GO TO (210,320,210,800), ICNT
  310 GO TO (710,320,210,210), ICNT
  320 GO TO (330,340,350,360), IENT
  330 GO TO (210,410,500,410), IEXIT
  340 GO TO (510,210,510,400), IEXIT
  350 GO TO (500,410,210,410), IEXIT
  360 GO TO (510,400,510,210), IEXIT
C     follow contour line across a cell to a side
  400 XSTART=SIDE(IENT)
  410 XFIN=SIDE(IEXIT)
      XINC=(XFIN-XSTART)/DIV
      XBASE=FLOAT(I-1)
      YBASE=FLOAT(J-1)
      A1=CL-C(2)
      A2=C(1)-C(2)
      A3=C(3)-C(2)
      A4=C(2)-C(1)+C(4)-C(3)
      DO 440 INTERP=1,IDIV
      XSTART=XSTART+XINC
      YSTART=(A1-A2*XSTART)/(A3+A4*XSTART)
      IF (IPT.LT.IPTMAX)  GO TO 430
      CALL APLOT(X, Y, IPT)
      X(1)=X(IPT)
      Y(1)=Y(IPT)
      IPT=1
  430 IPT=IPT+1
      X(IPT)=XBASE+XSTART
      Y(IPT)=YBASE+YSTART
  440 CONTINUE
      GO TO (230,210,615,635), JRTN
  500 YSTART=SIDE(IENT)
C     follow contour line across a cell to a top or bottom
  510 YFIN=SIDE(IEXIT)
      XBASE=FLOAT(I-1)
      YINC=(YFIN-YSTART)/DIV
      YBASE=FLOAT(J-1)
      A1=CL-C(2)
      A2=C(3)-C(2)
      A3=C(1)-C(2)
      A4=C(2)-C(1)+C(4)-C(3)
      DO 540 INTERP=1,IDIV
      YSTART=YSTART+YINC
      XSTART=(A1-A2*YSTART)/(A3+A4*YSTART)
      IF (IPT.LT.IPTMAX) GO TO 530
      CALL APLOT(X, Y, IPT)
      X(1)=X(IPT)
      Y(1)=Y(IPT)
      IPT=1
  530 IPT=IPT+1
      X(IPT)=XBASE+XSTART
      Y(IPT)=YBASE+YSTART
  540 CONTINUE
      GO TO (230,210,615,635), JRTN
C     follow contour line from corner to corner
  600 K1=ISUB-M
      K2=ISUB+1-M
      K3=ISUB+1
      K4=ISUB+1+M
      K5=ISUB+M
      K6=ISUB-1+M
      K7=ISUB-1
      C1=A(K1)
      C2=A(K2)
      C3=A(K3)
      C4=A(K4)
      C5=A(K5)
      C6=A(K6)
      C7=A(K7)
      C8=A(ISUB)
      IF (ISUB.LT.1 .OR. ISUB.GT.MN) GO TO 640
      X(1)=FLOAT(I-1)
      Y(1)=FLOAT(J-1)
      IF (J .EQ. 1 .OR. J .EQ. M)  GO TO 610
      IF (K1.LT.1 .OR. K1.GT.MN)  GO TO 610
      IF (K2.LT.1 .OR. K2.GT.MN)  GO TO 610
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 610
      IF (K4.LT.1 .OR. K4.GT.MN)  GO TO 610
      IF (K5.LT.1 .OR. K5.GT.MN)  GO TO 610
      IF (C3.NE.CL)  GO TO 610
      IF (C1 .EQ. CL .AND. C2 .EQ. CL .AND.
     *    C4 .EQ. CL .AND. C5 .EQ. CL)  GO TO 610
      X(2)=X(1)+1.
      Y(2)=Y(1)
      CALL APLOT (X, Y, 2)
      GO TO 620
  610 IF (J .EQ. 1)  GO TO 620
      IF (K1.LT.1 .OR. K1.GT.MN)  GO TO 620
      IF (K2.LT.1 .OR. K2.GT.MN)  GO TO 620
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 620
      IF (C2 .NE. CL)  GO TO 620
      IF (C1 .EQ. CL .OR. C3 .EQ. CL)  GO TO 620
      IF (C1 .GT. CL .AND. C3 .GT. CL .OR.
     *    C1 .LT. CL .AND. C3 .LT. CL)  GO TO 620
      C(1)=C2
      C(2)=C1
      C(3)=C8
      C(4)=C3
      J=J-1
      JRTN=3
      IENT=3
      IEXIT=1
      GO TO 500
  615 IF (IPT .GT. 1)  CALL APLOT (X, Y, IPT)
      IPT=1
      J=J+1
      X(1)=FLOAT(I-1)
      Y(1)=FLOAT(J-1)
  620 IF (J .EQ. M .OR. I .EQ. 1)  GO TO 630
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 630
      IF (K4.LT.1 .OR. K4.GT.MN)  GO TO 630
      IF (K5.LT.1 .OR. K5.GT.MN)  GO TO 630
      IF (K6.LT.1 .OR. K6.GT.MN)  GO TO 630
      IF (K7.LT.1 .OR. K7.GT.MN)  GO TO 630
      IF (C5 .NE. CL)  GO TO 630
      IF (C3 .EQ. CL .AND. C4 .EQ. CL .AND.
     *    C6 .EQ. CL .AND. C7 .EQ. CL)   GO TO 630
      X(2)=X(1)
      Y(2)=Y(1)+1.
      CALL APLOT (X, Y, 2)
      GO TO 640
  630 IF (J .EQ. M)  GO TO 640
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 640
      IF (K4.LT.1 .OR. K4.GT.MN)  GO TO 640
      IF (K5.LT.1 .OR. K5.GT.MN)  GO TO 640
      IF (C4 .NE. CL)  GO TO 640
      IF (C3 .EQ. CL .OR. C5 .EQ. CL)  GO TO 640
      IF (C3 .GT. CL .AND. C5 .GT. CL .OR.
     *    C3 .LT. CL .AND. C5 .LT. CL)  GO TO 640
      C(1)=C3
      C(2)=C8
      C(3)=C5
      C(4)=C4
      JRTN=4
      IENT=1
      IEXIT=3
      GO TO 500
  635 IF (IPT .GT. 1) CALL APLOT (X, Y, IPT)
      IPT=1
      X(1)=FLOAT(I-1)
      Y(1)=FLOAT(J-1)
  640 GO TO (110,160), IRTN
C    follow contour line from side to corner or corners
  700 JRTN=2
      IOPP=IDIR(IENT)
      I1=ISIDE(IOPP)
      I2=ISIDE(IOPP+1)
      IEXIT=IOPP
      C(I1)=C(KVCT(I1))
      C(I2)=C(LVCT(I2))
      GO TO 320
  710 JRTN=2
      IEXIT=JEXIT
      GO TO 320
C     follow contour line through saddle point
  800 IOPP=IDIR(IENT)
      I1=ISIDE(IENT)
      C1=C(I1)
      I2=ISIDE(IENT+1)
      C2=C(I2)
      I3=ISIDE(IOPP)
      C3=C(I3)
      I4=ISIDE(IOPP+1)
      C4=C(I4)
      IF ((C1-CL)/(C1-C2) .EQ. (C4-CL)/(C4-C3))  GO TO 820
      IF ((C1-CL)/(C1-C4) .GT. (C2-CL)/(C2-C3))  GO TO 810
      IEXIT=I4
      GO TO 320
  810 IEXIT=I2
      GO TO 320
  820 C(I3)=C(I2)
      C(I4)=C(I1)
      IEXIT=I3
      GO TO 320
C
C
  190 CONTINUE
      RETURN
      END
c
c
c
      subroutine rempi (map,iunit,nx,ny,nz,ierr)
c
c ... read MPI EM map
c
      implicit none
c
      integer nx,ny,nz
      real map(nx,ny,nz),dummy(64)
c
      integer iunit,ierr,i,j,k,l
c
code ...
c
      ierr = -1
c
      call prompt (' Reading data ...')
      rewind (iunit)
      read (iunit,end=9000,err=9000) (dummy(l),l=1,64),
     +  (((map(k,j,i),k=1,nx),j=1,ny),i=1,nz)
c
ccc     +  (((map(k,j,i),i=1,nz),j=1,ny),k=1,nx)
c
ccc      call rvalut (' MAP :',10,map(1,1,1))
c
      ierr = 0
      return
c
 9000 continue
      ierr = -1
      call errcon ('While reading MPI map')
c
      return
      end
c
c
c
      subroutine setval (map,nx,ny,nz,first,dx,dy,dz)
c
      implicit none
c
      integer i,j,k,nx,ny,nz
c
      real map(nx,ny,nz)
      real first,dx,dy,dz,xf,yf,zf
c
code ...
c
      do i=1,nx
        xf = first + (i-1)*dx
        do j=1,ny
          yf = xf + (j-1)*dy
          do k=1,nz
            map (i,j,k) = yf + (k-1)*dz
          end do
        end do
      end do
c
      return
      end
c
c
c
      subroutine grsfit (
     +  mapobs,npnt,no1,no2,no3,
     +  origin,grid,cell,
     +  mapcal,dumask,
     +  maskit,mappy1,mappy2,minmap,
     +  store1,store2,nstore,
     +  nowatc,nowatr,nowatb,mnowat,
     +  userad,iunit,junit,proscl,ierr)
c
      implicit none
c
c ... MAXAPR = max nr of atoms per residue
c
      integer maxapr
      parameter (maxapr = 1000)
c
      integer npnt,no1,no2,no3,nstore,minmap,mnowat
      real mapobs(npnt),mapcal(npnt)
      integer dumask(npnt)
      real mappy1(minmap),mappy2(minmap)
      real nowatc(mnowat),nowatr(mnowat),nowatb(mnowat)
      real store1(nstore),store2(nstore),cell(6)
      integer maskit(minmap)
c
      real atmxyz(3,maxapr)
      real occup(maxapr),bfact(maxapr)
      real g(3),rtunit(12),cog(3),cello(6),c2f(12),f2c(12)
      real c2x(12),glocog(3)
      real rsfrad,sumocc,srmsd,sshap,scorr,bave,xocc,resocc
      real srf1,srf2,srf3,srf4,ssf1,ssf2,srsr,sumdif,sumsum
      real sumfo,sumfc,scalef,gmin,userad,avy,avysq,avxy,avxmy
      real xpdb,xsdv,xmin,xmax,xtot
c
      integer origin(3),grid(3),orgnz(3),extz(3),extent(3),ctkind(4)
      integer iunit,junit,ierr,nres,natom,nline,i,j,nhydro
      integer natot,ngrid,nalt,nsum,nofit,nzeroc,nerror
      integer nhet,numhet,idum,nfunny,nzeres,ii,itype,k,length
      integer nnwc,nnwr,nnwb,n3sig,ncog,nrgt5,nrgt1,nuniat
      integer nbatch,ncount
c
      logical lalt,lhydro,leof,lerror,hasc,hasn,hasca,lwater
      logical lcrys1,lscal1,lscal2,lscal3,ldosca,proscl,doinit
c      
      character atmlin(maxapr)*80,line*80,nowres*10,prevat*4
      character iskind(4)*1,altids*30
c
      data rtunit / 1.0,0.0,0.0, 0.0,1.0,0.0,
     +              0.0,0.0,1.0, 0.0,0.0,0.0 /
c
      data iskind /'A','N','W','H'/
c
      data lcrys1,lscal1,lscal2,lscal3 /4*.false./
c
      data c2f /12*0.0/, f2c /12*0.0/, c2x /12*0.0/
c
code ...
c
      call prompt (' Initialising ...')
      ierr = 0
      do i=1,3
        g(i) = cell(i)/float(grid(i))
      end do
      gmin = min(g(1),g(2),g(3))
      rsfrad = userad
      call fvalut (' Grid spacing (A)      :',3,g)
      call jvalut (' Max size of mini-maps :',1,minmap)
      extent (1) = no1
      extent (2) = no2
      extent (3) = no3
c
      if (rsfrad .le. 0.6*gmin) then
        call errcon (' Masking radius too small - aborting')
        call fvalut (' Value (A) :',1,userad)
        rsfrad = 0.6 * gmin
        call fvalut (' Suggested minimum (A) :',1,rsfrad)
        ierr = -1
        return
      end if
c
      nres = 0
      natom = 0
      nuniat = 0
      prevat = '!@#$'
      nline = 0
      nhydro = 0
      natot = 0
      nalt = 0
      nsum = 0
      nofit = 0
      numhet = 0
      nzeroc = 0
      nzeres = 0
      nfunny = 0
      nerror = 0
      nowres = '!@#$%^&&*('
      leof = .false.
      altids = 'A'
      nnwc = 0
      nnwr = 0
      nnwb = 0
      ldosca = .false.
c
      do i=1,4
        ctkind (i) = 0
      end do
c
      ncog = 0
      do i=1,3
        glocog (i) = 0.0
      end do
c
      if (proscl) then
c
c ... masked scaling
c
        call prompt (' Doing masked scaling ...')
c
        ncount = 0
        nbatch = 0
c
 1810   continue
        read (iunit,6000,end=1815,err=9900) line
c
        if (line(1:6) .eq. 'ATOM  ' .or.
     +      line(1:6) .eq. 'HETATM') goto 1910
c
c ... CRYST1 etc.
c
        if (line(1:6) .eq. 'CRYST1') then
          read (line,'(6x,3f9.3,3f7.2)',err=1901) (cello(i),i=1,6)
          lcrys1 = .true.
          call orthog (cello,f2c,0)
          call orthog (cello,c2x,1)
          goto 1905
 1901     continue
          goto 1810
        end if
c
        if (line(1:6) .eq. 'SCALE1') then
          read (line,'(10x,3f10.6,f15.5)',err=1902) (c2f(i),i=1,10,3)
          lscal1 = .true.
          goto 1905
 1902     continue
          goto 1810
        end if
c
        if (line(1:6) .eq. 'SCALE2') then
          read (line,'(10x,3f10.6,f15.5)',err=1903) (c2f(i),i=2,11,3)
          lscal2 = .true.
          goto 1905
 1903     continue
          goto 1810
        end if
c
        if (line(1:6) .eq. 'SCALE3') then
          read (line,'(10x,3f10.6,f15.5)',err=1904) (c2f(i),i=3,12,3)
          lscal3 = .true.
          goto 1905
 1904     continue
          goto 1810
        end if
c
        goto 1810
c
c ... work out SCALE stuff once all info is there
c
 1905   continue
c
        if (lcrys1 .and. lscal1 .and. lscal2 .and. lscal3) then
          ldosca = .false.
          do i=1,9
            if (abs(c2f(i)-c2x(i)) .ge. 0.00001) then
              ldosca = .true.
              goto 1907
            end if
          end do
          do i=10,12
            if (abs(c2f(i)-c2x(i)) .ge. 0.0001) then
              ldosca = .true.
              goto 1907
            end if
          end do
 1907     continue
        end if
c
        goto 1810
c
c ... process ATOM and HETATM cards ...
c
 1910   continue
c
        if (lhydro(line(13:16))) goto 1810
        read (line(55:60),'(f6.2)') xocc
        if (xocc .gt. 0.0) then
          ncount = ncount + 1
          nbatch = nbatch + 1
          read (line(31:54),'(3f8.3)') (atmxyz(j,nbatch),j=1,3)
        end if
c
        if (nbatch .ge. maxapr) then
          if (ldosca) then
            call vecrtv (atmxyz,atmxyz,nbatch,c2f(1),c2f(10))
            call vecrtv (atmxyz,atmxyz,nbatch,f2c(1),f2c(10))
          end if
          doinit = (nbatch .eq. ncount)
          call grpmsk (cell,grid,nbatch,atmxyz,rsfrad,origin,
     +                 extent,dumask,npnt,doinit,ierr)
          call jvalut (' Atoms added to mask :',1,nbatch)
          nbatch = 0
        end if
c
        goto 1810
c
c ... end of file
c
 1815   continue
c
        if (nbatch .gt. 0) then
          if (ldosca) then
            call vecrtv (atmxyz,atmxyz,nbatch,c2f(1),c2f(10))
            call vecrtv (atmxyz,atmxyz,nbatch,f2c(1),f2c(10))
          end if
          doinit = (nbatch .eq. ncount)
          call grpmsk (cell,grid,nbatch,atmxyz,rsfrad,origin,
     +                 extent,dumask,npnt,doinit,ierr)
          call jvalut (' Atoms added to mask :',1,nbatch)
          nbatch = 0
        end if
c
        call jvalut (' Total atoms in mask :',1,ncount)
        if (ncount .lt. 1) then
          call errcon ('No atoms found - aborting')
          ierr = -1
          return
        end if
c
        sumfo = 0.0
        sumfc = 0.0
        do i=1,npnt
          if (dumask(i) .eq. 1) then
            sumfo = sumfo + abs(mapobs(i))
            sumfc = sumfc + abs(mapcal(i))
          end if
        end do
        call rvalut (' Sum |obs_map| in mask :',1,sumfo)
        call rvalut (' Sum |cal_map| in mask :',1,sumfc)
        if (sumfo .le. 0.0 .or. sumfc .le. 0.0) then
          call errcon ('Map is zero throughout - aborting')
          ierr = -1
          return
        end if
        scalef = sumfo/sumfc
        call rvalut (' Scale factor for cal_map :',1,scalef)
c
        write (junit,6010) '!'
        write (nowres,'(1pe10.3)') sumfo
        write (junit,6010) '! Sum |obs_map| in mask = ',nowres
        write (nowres,'(1pe10.3)') sumfc
        write (junit,6010) '! Sum |cal_map| in mask = ',nowres
        write (nowres,'(1pe10.3)') scalef
        write (junit,6010) '! Scale factor for cal_map = ',nowres
c
c ... rewind PDB file
c
        rewind (iunit)
c
        lcrys1 = .false.
        lscal1 = .false.
        lscal2 = .false.
        lscal3 = .false.
c
      else
c
c ... Quick-n-dirty scaling
c
        call prompt (' Doing quick-n-dirty scaling ...')
c
        sumfo = 0.0
        sumfc = 0.0
        do i=1,npnt
          sumfo = sumfo + abs(mapobs(i))
          sumfc = sumfc + abs(mapcal(i))
        end do
        call rvalut (' Sum |obs_map| :',1,sumfo)
        call rvalut (' Sum |cal_map| :',1,sumfc)
        if (sumfo .le. 0.0 .or. sumfc .le. 0.0) then
          call errcon ('Map is zero throughout - aborting')
          ierr = -1
          return
        end if
        scalef = sumfo/sumfc
        call rvalut (' Scale factor for cal_map :',1,scalef)
c
        write (junit,6010) '!'
        write (nowres,'(1pe10.3)') sumfo
        write (junit,6010) '! Sum |obs_map| = ',nowres
        write (nowres,'(1pe10.3)') sumfc
        write (junit,6010) '! Sum |cal_map| = ',nowres
        write (nowres,'(1pe10.3)') scalef
        write (junit,6010) '! Scale factor for cal_map = ',nowres
c
      end if
c
cATOM   1254  O   HIS   165      14.687  19.539  60.412  1.00 19.71      AAAA
cATOM   1255  N   ARG   166      16.090  18.542  58.968  1.00 19.08      AAAA
cATOM   2121  CB ATHR D  32      15.033  24.886   8.363  0.64 21.00           C  
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      call prompt (' Processing PDB file ...')
c
      write (junit,6010) '!'
      write (junit,6010)
     +  '! # = sequential residue index'
      write (junit,6010)
     +  '! residue_id = columns 18:27 of the ATOM/HETATM card'
      write (junit,6010)
     +  '! RSCC = residue-based real-space correlation coefficient'
      write (junit,6010)
     +  '! RSR = residue-based real-space R-value'
      write (junit,6010)
     +  '! OWAB = occupancy-weighted average B = (SUM B*Q)/(SUM Q)'
      write (junit,6010)
     +  '! Natom = nr of atoms encountered in residues'
      write (junit,6010)
     +  '! S_occ = sum of their occupancies (normally same as Natom)'
      write (junit,6010)
     +  '! Nhet = number of HETATM cards encountered in residue'
      write (junit,6010)
     +  '! AC = alternative conformation flag (T or F)'
      write (junit,6010)
     +  '! Scale_RHOc = residue-based density scale (not used)'
      write (junit,6010)
     +  '! Ngridpts = nr of grid points summed for RSCC and RSR'
      write (junit,6010)
     +  '! R? = Error? = computation-error flag (T or F)'
      write (junit,6010)
     +  '! K = Kind = kind of residue (Amino, Nucleic, Water, Hetero)'
      write (junit,6010)
     +  '! Nuniq = number of unique atom types in the residue'
      write (junit,6010)
     +  '! <Occ> = aver. occup. of the unique atoms = S_occ / Nuniq'
      write (junit,6010)
     +  '! <X>,<Y>,<Z> = occup.-weighted centre-of-gravity coords'
      write (junit,6010) '!'
      write (junit,6010)    '!        1         2         3',
     +  '         4         5         6         7         8',
     +  '         9        10        11        12        13'
      write (junit,6010) '!23456789+123456789+123456789+',
     +  '123456789+123456789+123456789+123456789+123456789+',
     +  '123456789+123456789+123456789+123456789+123456789+'
c
      write (junit,6010) '!    # [','residue_id',']   RSCC',
     +  '   RSR','    OWAB',' Natom',' S_occ','  Nhet',' AC',
     +  '  Scale_RHOc','  Ngridpts',' R?',' K',' Nuniq',
     +  ' <Occ>','       <X>','       <Y>','       <Z>'
c
 6010 format (20a)
c
   10 continue
      read (iunit,6000,end=15,err=9900) line
      nline = nline + 1
c
      if (line(1:6) .eq. 'ATOM  ' .or.
     +    line(1:6) .eq. 'HETATM') goto 910
c
c ... CRYST1
c
      if (line(1:6) .eq. 'CRYST1') then
        call textut (' ... ',line)
        read (line,'(6x,3f9.3,3f7.2)',err=901) (cello(i),i=1,6)
        lcrys1 = .true.
        call orthog (cello,f2c,0)
        call orthog (cello,c2x,1)
ccc        call fvalut (' F2C matrix :',12,f2c)
        goto 905
  901   continue
        goto 10
      end if
c
c ... SCALE1
c
      if (line(1:6) .eq. 'SCALE1') then
        call textut (' ... ',line)
        read (line,'(10x,3f10.6,f15.5)',err=902) (c2f(i),i=1,10,3)
        lscal1 = .true.
        goto 905
  902   continue
        goto 10
      end if
c
c ... SCALE2
c
      if (line(1:6) .eq. 'SCALE2') then
        call textut (' ... ',line)
        read (line,'(10x,3f10.6,f15.5)',err=903) (c2f(i),i=2,11,3)
        lscal2 = .true.
        goto 905
  903   continue
        goto 10
      end if
c
c ... SCALE3
c
      if (line(1:6) .eq. 'SCALE3') then
        call textut (' ... ',line)
        read (line,'(10x,3f10.6,f15.5)',err=904) (c2f(i),i=3,12,3)
        lscal3 = .true.
        goto 905
  904   continue
        goto 10
      end if
c
      goto 10
c
c ... work out SCALE stuff once all info is there
c
  905 continue
c
      if (lcrys1 .and. lscal1 .and. lscal2 .and. lscal3) then
ccc        call fvalut (' C2F matrix :',12,c2f)
ccc        call fvalut (' C2X matrix :',12,c2x)
        ldosca = .false.
        do i=1,9
          if (abs(c2f(i)-c2x(i)) .ge. 0.00001) then
            ldosca = .true.
            goto 907
          end if
        end do
        do i=10,12
          if (abs(c2f(i)-c2x(i)) .ge. 0.0001) then
            ldosca = .true.
            goto 907
          end if
        end do
  907   continue
        if (ldosca) then
          call prompt (
     +      ' Fractionalisation matrix from cell constants')
          call prompt (
     +      ' differs from PDB SCALEn matrix -> will apply')
          call prompt (
     +      ' fractionalisation and orthogonalisation !')
        else
          call prompt (
     +      ' No fractionalisation and orthogonalisation needed')
        end if
      end if
c
      goto 10
c
c ... process ATOM and HETATM cards ...
c
  910 continue
c
      if (lhydro(line(13:16))) then
        nhydro = nhydro + 1
        goto 10
      end if
c
ccc      read (line(55:60),*) xocc
ccc      if (xocc .lt. 0.01) then
ccc        nzeroc = nzeroc + 1
ccc        goto 10
ccc      end if
c
      natot = natot + 1
c
      if (line(17:17) .ne. ' ') then
        i = length(altids)
        if (index(altids(1:i),line(17:17)) .le. 0) then
          altids (i+1:i+1) = line(17:17)
        end if
      end if
c
      if (natom .eq. 0) then
        natom = 1
        nhet = 0
        nuniat = 1
        if (line(1:6) .eq. 'HETATM') nhet = nhet + 1
        nowres = line (18:27)
        prevat = line (13:16)
        atmlin (1) = line
        goto 10
      else
        if (nowres .ne. line(18:27)) then
          backspace (iunit)
          nline = nline - 1
          natot = natot - 1
          prevat = '!@#$'
          goto 20
        else
          natom = natom + 1
          if (prevat .ne. line(13:16)) nuniat = nuniat + 1
          if (line(1:6) .eq. 'HETATM') nhet = nhet + 1
          if (natom .le. maxapr) atmlin (natom) = line
          prevat = line (13:16)
          goto 10
        end if
      end if
c
c ... shouldn't get here
c
      call errcon ('Huh ? Programming error ? Tell Gerard !')
      goto 10
c
c ... end of file
c
   15 continue
      leof = .true.
      if (natom .le. 0) goto 9000
c
c ... end of a residue
c
   20 continue
      nres = nres + 1
      if (nhet .gt. 0) numhet = numhet + 1
      sumocc = 0.0
      lalt = .false.
      bave = 0.0
c
      if (natom .gt. maxapr) then
        call errcon ('Too many atoms in residue - skipped')
        nofit = nofit + 1
        goto 1234
      end if
c
c ... figure out type of residue: 1=Amino acid, 2=Nucleic acid,
c     3=Water, 4=Hetero-compound
c
      itype = 4
      if (lwater(nowres(1:3))) then
        itype =3
      else
        call nuctyp (nowres(1:3),k)
        if (k .gt. 0) then
          itype = 2
        else
          hasn = .false.
          hasc = .false.
          hasca = .false.
          do i=1,natom
            if (atmlin(i)(13:16) .eq. ' CA ') then
              hasca = .true.
            else if (atmlin(i)(13:16) .eq. ' C  ') then
              hasc = .true.
            else if (atmlin(i)(13:16) .eq. ' N  ') then
              hasn = .true.
            end if
          end do
          if (natom .eq. 1 .and. hasca) itype = 1
          if (hasca .and. hasc .and. hasn) itype = 1
        end if
      end if
      ctkind (itype) = ctkind (itype) + 1
c
      ii = 0
c
      do j=1,3
        cog(j) = 0.0
      end do
c
      do i=1,natom
ccc        write (*,6000) atmlin(i)
        read (atmlin(i)(55:60),'(f6.2)') xocc
        if (xocc .gt. 0.0) then
          ii = ii + 1
          occup (ii) = xocc
          read (atmlin(i)(31:54),'(3f8.3)') (atmxyz(j,ii),j=1,3)
          read (atmlin(i)(61:66),'(f6.2)') bfact(ii)
c
          cog(1) = cog(1) + xocc*atmxyz(1,ii)
          cog(2) = cog(2) + xocc*atmxyz(2,ii)
          cog(3) = cog(3) + xocc*atmxyz(3,ii)
c
          sumocc = sumocc + occup(ii)
          bave = bave + occup(ii)*bfact(ii)
          if (.not. lalt) then
            lalt = (atmlin(ii)(17:17) .ne. ' ')
          end if
          if (itype .ne. 3) then
            ncog = ncog + 1
            do j=1,3
              glocog (j) = glocog (j) + atmxyz(j,ii)
            end do
          end if
        else
          nzeroc = nzeroc + 1
        end if
      end do
c
      if (lalt) nalt = nalt + 1
      if (abs(float(natom)-sumocc) .ge. 0.001) nsum = nsum + 1
c
c ... calculate residue centre-of-gravity weighted by occupancy
c
      if (sumocc .gt. 0.0) then
        do j=1,3
          cog (j) = cog(j)/sumocc
        end do
        if (ldosca) then
          call vecrtv (cog,cog,3,c2f(1),c2f(10))
          call vecrtv (cog,cog,3,f2c(1),f2c(10))
        end if
      else
        do j=1,3
          cog (j) = 0.0
        end do
      end if
c
      natom = ii
c
ccc      print *,' NCOG ',ncog
c
c ... need to fract/orthog ???
c
      if (ldosca) then
        if (natom .gt. 0) then
ccc          call fvalut (' Atom1 #1 :',3,atmxyz)
          call vecrtv (atmxyz,atmxyz,natom,c2f(1),c2f(10))
ccc          call fvalut (' Atom1 #2 :',3,atmxyz)
          call vecrtv (atmxyz,atmxyz,natom,f2c(1),f2c(10))
ccc          call fvalut (' Atom1 #3 :',3,atmxyz)
        end if
      end if
c
      if (sumocc .gt. 0.0) then
        bave = bave / sumocc
      end if
c
      scorr = -9.99
      srsr  = -9.99
      ngrid = 0
      ssf1 = 1.0
      lerror = .false.
c
c ... skip if sum of occupancies is zero
c
      if (natom .lt. 1) then
        call prompt (
     +    ' WARNING - Only zero occupancy atoms in residue - skipped')
        nzeres = nzeres + 1
ccc        lerror = .true.
        goto 1234
      end if
c
c ... do RS-fit calculations for this residue
c
c ... calculate size of mask for this residue
c
      call sizmsk (cell,grid,natom,atmxyz,rsfrad,orgnz,extz,ierr)
c
ccc      call jvalut (' Origin :',3,orgnz)
ccc      call jvalut (' Extent :',3,extz)
c
      idum = extz(1)*extz(2)*extz(3)
c
ccc      call jvalut (' Calculated mini-map/masksize :',1,idum)
c
      if (idum .gt. minmap .or. idum .lt. 1) then
        call errcon ('Invalid mini-map/masksize - skipped')
        nofit = nofit + 1
        lerror = .true.
        goto 1234
      end if
c
c ... initialising not necessary since it is done inside GRPMSK and SUB02X
c
c      do i=1,idum
c        mappy1 (i) = 0.0
c        mappy2 (i) = 0.0
c        maskit (i) = 0
c      end do
c
c ... actually generate the mask
c
      doinit = .true.
      call grpmsk (cell,grid,natom,atmxyz,rsfrad,orgnz,extz,
     +             maskit,minmap,doinit,ierr)
c
c ... extract MapObs around this residue
c
      call sub02x (
     +    mapobs,no1,no2,no3,origin,g,cell,
     +    mappy1,extz(1),extz(2),extz(3),orgnz,g,cell,
     +    maskit, rtunit, 1, rtunit, 1,
     +    rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
c
c ... extract MapCal around this residue
c
      call sub02x (
     +    mapcal,no1,no2,no3,origin,g,cell,
     +    mappy2,extz(1),extz(2),extz(3),orgnz,g,cell,
     +    maskit, rtunit, 1, rtunit, 1,
     +    rtunit,ierr,avy,avysq,avxy,avxmy,.false.)
c
      ngrid = 0
      idum = extz(1)*extz(2)*extz(3)
      do i=1,idum
        if (maskit(i) .eq. 1) then
          ngrid = ngrid + 1
          store1 (ngrid) = mappy1 (i)
          store2 (ngrid) = mappy2 (i)
          maskit (i) = 0
        end if
      end do
c
c ... RF1 = R-factor with respect to XDATA
c ... RF2 = R-factor with respect to YDATA
c ... RF3 = scaled R-factor wrt XDATA, scale SF1 for YDATA
c ... RF4 = scaled R-factor wrt YDATA, scale SF2 for XDATA
c ... CORR = correlation coefficient
c
      if (ngrid .ge. 3) then
        call xystat (store1,store2,ngrid,
     +    srmsd,sshap,scorr,srf1,srf2,srf3,srf4,ssf1,ssf2)
        sumdif = 0.0
        sumsum = 0.0
        do i=1,ngrid
c
c ... use global scale factor for MapCal
c
ccc          store2 (i) = ssf1 * store2 (i)
          store2 (i) = scalef * store2 (i)
          sumdif = sumdif + abs (store1(i)-store2(i))
          sumsum = sumsum + abs (store1(i)+store2(i))
        end do
        if (sumsum .gt. 0) srsr = sumdif / sumsum
      else
        call errcon ('Fewer than 3 points in residue mask - skipped')
        nofit = nofit + 1
        lerror = .true.
        goto 1234
      end if
c
      if (srsr .lt. 0.0) then
        call errcon ('Negative real-space R-value - funny')
        nfunny = nfunny + 1
        lerror = .true.
        goto 1234
      end if
c
c ... write result to screen and list file
c
 1234 continue
c
      if (lerror) nerror = nerror + 1
c
c ... average residue occupancy = SUM occ / Nr of unique atoms
c
      resocc = 0.0
      if (nuniat .gt. 0) resocc = sumocc / float(nuniat)
c
      write (junit,6200,err=25) nres,nowres,scorr,srsr,bave,
     +  natom,sumocc,nhet,lalt,ssf1,ngrid,lerror,iskind(itype),
     +  nuniat,resocc,cog(1),cog(2),cog(3)
c
      write (*,6100) nres,nowres,natom,scorr,srsr,bave,resocc
c
c ... store non-water values
c
      if (itype .ne. 3) then
        if (scorr .ge. -1.0 .and. scorr .le. 1.0) then
          nnwc = nnwc + 1
          nowatc (nnwc) = scorr
        end if
        if (srsr .ge. 0.0) then
          nnwr = nnwr + 1
          nowatr (nnwr) = srsr
        end if
        if (bave .ge. 0.0) then
          nnwb = nnwb + 1
          nowatb (nnwb) = bave
        end if
      end if
c
 6100 format (' # ',i6,' [',a,'] Nat ',i6,
     +  ' RSCC ',f6.3,' RSR ',f6.3,' <B> ',f8.2,' <OCC> ',f6.3)
c
 6200 format (i6,' [',a10,'] ',f6.3,f6.3,f8.2,i6,f6.2,i6,1x,l1,
     +  1x,1pe12.4,0p,i10,1x,l1,2x,a1,i6,f6.3,3f10.3)
c
   25 continue
c
      natom = 0
      nowres = '!@#$%^&&*('
c
      if (.not. leof) goto 10
c
      goto 9000
c
c ... read error
c
 9900 continue
      call errcon ('While reading PDB file - aborting')
c
c ... all done
c
 9000 continue
      write (*,*)
      if (ldosca) then
        call prompt (' CRYST1 and SCALEn information WAS used')
      else
        call prompt (' CRYST1 and SCALEn information was NOT used')
      end if
      call jvalut (' Nr of lines read from file   :',1,nline)
      call jvalut (' Nr of hydrogen atoms skipped :',1,nhydro)
      call jvalut (' Nr of zero-occ atoms         :',1,nzeroc)
      call jvalut (' Nr of non-H atoms read       :',1,natot)
      call jvalut (' Nr of residues encountered   :',1,nres)
      call jvalut (' ... of which amino acids     :',1,ctkind(1))
      call jvalut (' ... of which nucleic acids   :',1,ctkind(2))
      call jvalut (' ... of which waters          :',1,ctkind(3))
      call jvalut (' ... of which heteros         :',1,ctkind(4))
      call jvalut (' Nr with alt. conf.           :',1,nalt)
      i = length(altids)
      if (i .gt. 1) then
        do j=2,i
          call textut (' Alt. loc. identifier         :',altids(j:j))
        end do
      end if
      call jvalut (' Nr with HETERO atoms         :',1,numhet)
      call jvalut (' Nr with zero summed occup.   :',1,nzeres)
      call jvalut (' Nr with odd summed occup.    :',1,nsum)
      call jvalut (' Nr not calculated at all     :',1,nofit)
      call jvalut (' Nr with computation errors   :',1,nfunny)
      call jvalut (' Total nr with errors         :',1,nerror)
c
      write (junit,6010) '!'
      if (ldosca) then
        write (junit,6010)
     +    '! CRYST1 and SCALEn information WAS used'
      else
        write (junit,6010)
     +    '! CRYST1 and SCALEn information was NOT used'
      end if
c
      write (nowres,'(i10)') nline
      write (junit,6010) '! Nr of lines read from file   : ',nowres
      write (nowres,'(i10)') nhydro
      write (junit,6010) '! Nr of hydrogen atoms skipped : ',nowres
      write (nowres,'(i10)') nzeroc
      write (junit,6010) '! Nr of zero-occ atoms         : ',nowres
      write (nowres,'(i10)') natot
      write (junit,6010) '! Nr of non-H atoms read       : ',nowres
      write (nowres,'(i10)') nres
      write (junit,6010) '! Nr of residues encountered   : ',nowres
c
      write (nowres,'(i10)') ctkind(1)
      write (junit,6010) '! ... of which amino acids     : ',nowres
      write (nowres,'(i10)') ctkind(2)
      write (junit,6010) '! ... of which nucleic acids   : ',nowres
      write (nowres,'(i10)') ctkind(3)
      write (junit,6010) '! ... of which waters          : ',nowres
      write (nowres,'(i10)') ctkind(4)
      write (junit,6010) '! ... of which heteros         : ',nowres
c
      write (nowres,'(i10)') nalt
      write (junit,6010) '! Nr with alt. conf.           : ',nowres
      i = length(altids)
      if (i .gt. 1) then
        do j=2,i
          write (junit,6010)
     +      '! Alt. loc. identifier         : ',altids(j:j)
        end do
      end if
      write (nowres,'(i10)') numhet
      write (junit,6010) '! Nr with HETERO atoms         : ',nowres
      write (nowres,'(i10)') nzeres
      write (junit,6010) '! Nr with zero summed occup.   : ',nowres
      write (nowres,'(i10)') nsum
      write (junit,6010) '! Nr with odd summed occup.    : ',nowres
      write (nowres,'(i10)') nofit
      write (junit,6010) '! Nr not calculated at all     : ',nowres
      write (nowres,'(i10)') nfunny
      write (junit,6010) '! Nr with computation errors   : ',nowres
      write (nowres,'(i10)') nerror
      write (junit,6010) '! Total nr with errors         : ',nowres
c
      if (ncog .gt. 0) then
        write (junit,6010) '!'
        write (*,*)
        call jvalut (' Non-water atoms in CoG calcn :',1,ncog)
        write (nowres,'(i10)') ncog
        write (junit,6010) '! Non-water atoms in CoG calcn : ',nowres
ccc        print *,' NCOG 1 ',ncog
c
        do i=1,3
          glocog (i) = glocog(i) / float (ncog)
        end do
ccc        print *,' NCOG 2 ',ncog
c
c ... on SGI, this IF block destroys the value of NCOG !!!???
c
        if (ldosca) then
ccc          call fvalut (' C-of-G #1 :',3,glocog)
          call vecrtv (glocog,glocog,3,c2f(1),c2f(10))
ccc          call fvalut (' C-of-G #2 :',3,glocog)
          call vecrtv (glocog,glocog,3,f2c(1),f2c(10))
ccc          call fvalut (' C-of-G #3 :',3,glocog)
        end if
ccc        print *,' NCOG 3 ',ncog
c
        call fvalut (' ... Centre-of-Gravity X      :',1,glocog(1))
        call fvalut (' ... Centre-of-Gravity Y      :',1,glocog(2))
        call fvalut (' ... Centre-of-Gravity Z      :',1,glocog(3))
        write (nowres,'(f10.3)') glocog(1)
        write (junit,6010) '! ... Centre-of-Gravity X      : ',nowres
        write (nowres,'(f10.3)') glocog(2)
        write (junit,6010) '! ... Centre-of-Gravity Y      : ',nowres
        write (nowres,'(f10.3)') glocog(3)
        write (junit,6010) '! ... Centre-of-Gravity Z      : ',nowres
      end if
c
      write (junit,6010) '!'
      write (*,*)
      call jvalut (' Non-water residues with RSCC :',1,nnwc)
      write (nowres,'(i10)') nnwc
      write (junit,6010) '! Non-water residues with RSCC : ',nowres
      if (nnwc .gt. 1) then
        call xstats (nowatc,nnwc,xpdb,xsdv,xmin,xmax,xtot)
        call fvalut (' ... Average RSCC for these   :',1,xpdb)
        call fvalut (' ... St. dev. RSCC            :',1,xsdv)
        write (nowres,'(f10.4)') xpdb
        write (junit,6010) '! ... Average RSCC for these   : ',nowres
        write (nowres,'(f10.4)') xsdv
        write (junit,6010) '! ... St. dev. RSCC            : ',nowres
        n3sig = 0
        xtot = xpdb - 3.0*xsdv
        call fvalut (' ... 3 Sigma cut-off          :',1,xsdv)
        write (nowres,'(f10.4)') xtot
        write (junit,6010) '! ... 3 Sigma cut-off          : ',nowres
        do i=1,nnwc
          if (nowatc(i) .lt. xtot) n3sig = n3sig + 1
        end do
        call jvalut (' ... 3 Sigma outliers (nr)    :',1,n3sig)
        write (nowres,'(i10)') n3sig
        write (junit,6010) '! ... 3 Sigma outliers         : ',nowres
        xmin = 100.0 * float (n3sig) / float (nnwc)
        call fvalut (' ... 3 Sigma outliers (%-age) :',1,xmin)
        write (nowres,'(f10.1)') xmin
        write (junit,6010) '! ... 3 Sigma outliers (%-age) : ',nowres
      end if
c
      write (junit,6010) '!'
      write (*,*)
      call jvalut (' Non-water residues with RSR  :',1,nnwr)
      write (nowres,'(i10)') nnwr
      write (junit,6010) '! Non-water residues with RSR  : ',nowres
      if (nnwr .gt. 1) then
        call xstats (nowatr,nnwr,xpdb,xsdv,xmin,xmax,xtot)
        call fvalut (' ... Average RSR for these    :',1,xpdb)
        call fvalut (' ... St. dev. RSR             :',1,xsdv)
        write (nowres,'(f10.4)') xpdb
        write (junit,6010) '! ... Average RSR for these    : ',nowres
        write (nowres,'(f10.4)') xsdv
        write (junit,6010) '! ... St. dev. RSR             : ',nowres
        n3sig = 0
        xtot = xpdb + 3.0*xsdv
        call fvalut (' ... 3 Sigma cut-off          :',1,xsdv)
        write (nowres,'(f10.4)') xtot
        write (junit,6010) '! ... 3 Sigma cut-off          : ',nowres
        do i=1,nnwr
          if (nowatr(i) .gt. xtot) n3sig = n3sig + 1
        end do
        call jvalut (' ... 3 Sigma outliers (nr)    :',1,n3sig)
        write (nowres,'(i10)') n3sig
        write (junit,6010) '! ... 3 Sigma outliers         : ',nowres
        xmin = 100.0 * float (n3sig) / float (nnwr)
        call fvalut (' ... 3 Sigma outliers (%-age) :',1,xmin)
        write (nowres,'(f10.1)') xmin
        write (junit,6010) '! ... 3 Sigma outliers (%-age) : ',nowres
c
        nrgt1 = 0
        nrgt5 = 0
        do i=1,nnwr
          if (nowatr(i) .gt. 0.5) then
            nrgt5 = nrgt5 + 1
            if (nowatr(i) .gt. 1.0) then
              nrgt1 = nrgt1 + 1
            end if
          end if
        end do
        call jvalut (' Non-water residues RSR > 0.5 :',1,nrgt5)
        write (nowres,'(i10)') nrgt5
        write (junit,6010) '! Non-water residues RSR > 0.5 : ',nowres
        call jvalut (' Non-water residues RSR > 1.0 :',1,nrgt1)
        write (nowres,'(i10)') nrgt1
        write (junit,6010) '! Non-water residues RSR > 1.0 : ',nowres
      end if
c
      write (junit,6010) '!'
      write (*,*)
      call jvalut (' Non-water residues with OWAB :',1,nnwb)
      write (nowres,'(i10)') nnwb
      write (junit,6010) '! Non-water residues with OWAB : ',nowres
      if (nnwb .gt. 1) then
        call xstats (nowatb,nnwb,xpdb,xsdv,xmin,xmax,xtot)
        call fvalut (' ... Average OWAB for these   :',1,xpdb)
        call fvalut (' ... St. dev. OWAB            :',1,xsdv)
        write (nowres,'(f10.2)') xpdb
        write (junit,6010) '! ... Average OWAB for these   : ',nowres
        write (nowres,'(f10.2)') xsdv
        write (junit,6010) '! ... St. dev. OWAB            : ',nowres
        n3sig = 0
        xtot = xpdb + 3.0*xsdv
        call fvalut (' ... 3 Sigma cut-off          :',1,xsdv)
        write (nowres,'(f10.2)') xtot
        write (junit,6010) '! ... 3 Sigma cut-off          : ',nowres
        do i=1,nnwb
          if (nowatb(i) .gt. xtot) n3sig = n3sig + 1
        end do
        call jvalut (' ... 3 Sigma outliers (nr)    :',1,n3sig)
        write (nowres,'(i10)') n3sig
        write (junit,6010) '! ... 3 Sigma outliers         : ',nowres
        xmin = 100.0 * float (n3sig) / float (nnwb)
        call fvalut (' ... 3 Sigma outliers (%-age) :',1,xmin)
        write (nowres,'(f10.1)') xmin
        write (junit,6010) '! ... 3 Sigma outliers (%-age) : ',nowres
      end if
c
      write (junit,6010) '!'
c
c ... don't remember why I do the following ...
c
      do i=1,minmap
        if (maskit(i) .ne. 0) maskit (i) = 1
      end do
c
 6000 format (a)
c
      return
      end
c
c
c
      subroutine gmfmap (
     +  fitmap,no1,no2,no3,
     +  origin,grid,cell,
     +  mapobs,mapcal,
     +  maskit,minmap,
     +  store1,store2,nstore,
     +  what,how,mfrad,ierr)
c
      implicit none
c
      integer no1,no2,no3,nstore,minmap
      real fitmap(no1,no2,no3),mapobs(no1,no2,no3),mapcal(no1,no2,no3)
      real store1(nstore),store2(nstore)
      integer maskit(minmap)
c
      real cell(6),g(3),xyz(3)
      real mfrad,srmsd,sshap,scorr
      real srf1,srf2,srf3,srf4,ssf1,ssf2,sumdif,sumsum
      real sumfo,sumfc,scalef,setto
c
      integer origin(3),grid(3),orgnz(3),extz(3)
      integer ierr,i,j,k,l,nokay,ngrid
      integer io1,io2,io3,npnt,ii,jj,kk,idum
c
      logical lcc,lpar,doinit
c      
      character how*(*),what*(*)
c
      data xyz /0.0,0.0,0.0/
c
code ...
c
      write (*,*)
      call prompt (' Initialising ...')
      ierr = 0
c
      lcc = (what(1:1) .ne. 'R')
      if (lcc) then
        call prompt (
     +    ' Calculating local correlation-coefficient map')
      else
        call prompt (
     +    ' Calculating local R-value map')
      end if
c
      lpar = (how(1:1) .ne. 'S')
      if (lpar) then
        call prompt (
     +    ' using a parallellopiped around each grid point')
      else
        call prompt (
     +    ' using a sphere around each grid point')
      end if
c
c      if (.not. lpar) then
c        call errcon ('Spheres option not implemented yet !')
c        call prompt (' Using parallellopipeds instead')
c        lpar = .true.
c      end if
c
      do i=1,3
        g(i) = cell(i)/float(grid(i))
      end do
      call fvalut (' Grid spacing (A) :',3,g)
c
      io1 = max (1, int(0.1 + (mfrad/g(1))))
      io2 = max (1, int(0.1 + (mfrad/g(2))))
      io3 = max (1, int(0.1 + (mfrad/g(3))))
      call ivalut (' Range +/- I :',1,io1)
      call ivalut (' Range +/- J :',1,io2)
      call ivalut (' Range +/- K :',1,io3)
c
      if (lpar) then
c
        ngrid = (2*io1+1) * (2*io2+1) * (2*io3+1)
        call ivalut (' Grid points in parallellopiped :',1,ngrid)
        if (ngrid .gt. nstore) then
          call errcon ('Too many points - aborting')
          call ivalut (' Max allowed :',1,nstore)
          ierr = -1
          return
        end if
c
      else
c
        call sizmsk (cell,grid,1,xyz,mfrad,orgnz,extz,ierr)
        call jvalut (' Origin :',3,orgnz)
        call jvalut (' Extent :',3,extz)
        idum = extz(1)*extz(2)*extz(3)
        call jvalut (' Calculated mini-mask size :',1,idum)
        if (idum .gt. minmap .or. idum .lt. 1) then
          call errcon ('Mini-mask too large - aborting')
          call ivalut (' Max allowed :',1,minmap)
          ierr = -1
          return
        end if
        doinit = .true.
        call grpmsk (cell,grid,1,xyz,mfrad,orgnz,extz,
     +               maskit,minmap,doinit,ierr)
        nokay = 0
        do i=1,idum
          if (maskit(i) .eq. 1) nokay = nokay + 1
        end do
        call ivalut (' Grid points in sphere :',1,nokay)
        if (nokay .gt. nstore) then
          call errcon ('Too many points - aborting')
          call ivalut (' Max allowed :',1,nstore)
          ierr = -1
          return
        end if
c
      end if
c
      npnt = no1 * no2 * no3
      setto = 2.0
      if (lcc) setto = -2.0
      call mfinit (fitmap,mapobs,mapcal,npnt,sumfo,sumfc,setto)
      call rvalut (' Sum |map 1| :',1,sumfo)
      call rvalut (' Sum |map 2| :',1,sumfc)
      if (sumfo .le. 0.0 .or. sumfc .le. 0.0) then
        call errcon ('Map is zero throughout - aborting')
        ierr = -1
        return
      end if
      scalef = sumfo/sumfc
      call rvalut (' Scale factor for map 2 :',1,scalef)
c
      write (*,*)
      call prompt (' Processing ...')
c
      if (lpar) then
c
c ... parallellopiped
c
        do i=1+io1,no1-io1
          do j=1+io2,no2-io2
            do k=1+io3,no3-io3
c
              nokay = 0
              do ii=i-io1,i+io1
                do jj=j-io2,j+io2
                  do kk=k-io3,k+io3
                    nokay = nokay + 1
                    store1 (nokay) = mapobs(ii,jj,kk)
                    store2 (nokay) = mapcal(ii,jj,kk)
                  end do
                end do
              end do
c
              if (lcc) then
                call xystat (store1,store2,nokay,
     +            srmsd,sshap,scorr,srf1,srf2,srf3,srf4,ssf1,ssf2)
                fitmap (i,j,k) = scorr
              else
                sumdif = 0.0
                sumsum = 0.0
                do l=1,ngrid
                  store2 (l) = scalef * store2 (l)
                  sumdif = sumdif + abs (store1(l)-store2(l))
                  sumsum = sumsum + abs (store1(l)+store2(l))
                end do
                if (sumsum .gt. 0) then
                  fitmap (i,j,k) = sumdif / sumsum
                end if
              end if
c
            end do
          end do
        end do
c
c ... sphere
c
      else
c
        do i=1+io1,no1-io1
          do j=1+io2,no2-io2
            do k=1+io3,no3-io3
c
              nokay = 0
              do ii=i-io1,i+io1
                do jj=j-io2,j+io2
                  do kk=k-io3,k+io3
ccc       print *,ii-i-orgnz(1)+1,jj-j-orgnz(2)+1,kk-k-orgnz(3)+1
                    call getmsk(maskit,
     +                extz(1),extz(2),extz(3),
     +                ii-i-orgnz(1)+1,
     +                jj-j-orgnz(2)+1,
     +                kk-k-orgnz(3)+1,
     +                idum,ierr)
                    if (ierr .eq. 0 .and. idum .eq. 1) then
                      nokay = nokay + 1
                      store1 (nokay) = mapobs(ii,jj,kk)
                      store2 (nokay) = mapcal(ii,jj,kk)
                    end if
                  end do
                end do
              end do
c
              if (nokay .gt. 0) then
                if (lcc) then
                  call xystat (store1,store2,nokay,
     +              srmsd,sshap,scorr,srf1,srf2,srf3,srf4,ssf1,ssf2)
                  fitmap (i,j,k) = scorr
                else
                  sumdif = 0.0
                  sumsum = 0.0
                  do l=1,nokay
                    store2 (l) = scalef * store2 (l)
                    sumdif = sumdif + abs (store1(l)-store2(l))
                    sumsum = sumsum + abs (store1(l)+store2(l))
                  end do
                  if (sumsum .gt. 0) then
                    fitmap (i,j,k) = sumdif / sumsum
                  end if
                end if
              end if
c
ccc              write (*,'(3i6,i8,1pe12.4)') i,j,k,nokay,fitmap(i,j,k)
c
            end do
          end do
        end do
c
      end if
c
      call prompt (' Finished !')
c
 1234 continue
c
      return
      end
c
c
c
      subroutine mfinit (fitmap,mapobs,mapcal,npnt,sumfo,sumfc,setto)
c
      implicit none
c
      integer npnt
c
      real fitmap(npnt),mapobs(npnt),mapcal(npnt)
      real sumfo,sumfc,setto
c
      integer i
c
code ...
c
      sumfo = 0.0
      sumfc = 0.0
      do i=1,npnt
        fitmap (i) = setto
        sumfo = sumfo + abs(mapobs(i))
        sumfc = sumfc + abs(mapcal(i))
      end do
c
      return
      end
c
c
c
      subroutine sizmsk (cell,grid,maxatm,xyzr,radius,
     +                   origin,extent,ierr)
c
c ... calculate origin and extent of mini-mask around a set of atoms
c
      implicit none
c
      integer maxatm
      integer i,j,extent(3),grid(3),origin(3),ct,ierr,nextra
c
      real cell(6),g(3),xyzr(3,maxatm)
      real b(3,3),ming(3),minr,maxg(3),maxr,radius
      real x1(3),x2(3),minc(3),maxc(3)
c
code ...
c
      ierr = -1
c
      minr = max (0.1, radius)
      maxr = minr
      nextra = 3
c
c --- Check the coords to see the required origin and extent
c
      call orthog (cell, b, 1)
c
      do i=1,3
        ming(i) =  9999999.0
        maxg(i) = -9999999.0
        minc(i) =  9999999.0
        maxc(i) = -9999999.0
      end do
c
      do 150 i=1,maxatm
        do 160 j=1,3
160       x1(j) = xyzr(j,i)
        call mulmtx (b, x1, x2, 3, 3, 1)
        do 170 j=1,3
c
          minc (j) = min (minc(j),xyzr(j,i))
          maxc (j) = max (minc(j),xyzr(j,i))
c
          x1(j) = x2(j)*float(grid(j))
c
          ming (j) = min (ming(j),x1(j))
          maxg (j) = max (maxg(j),x1(j))
c
170     continue
150   continue
c
c      call fvalut (' Lower bounds (coordinates) :',3,minc)
c      call fvalut (' Upper bounds (coordinates) :',3,maxc)
c      call fvalut (' Lower bounds (grid points) :',3,ming)
c      call fvalut (' Upper bounds (grid points) :',3,maxg)
c
      do i=1,3
        g(i) = cell(i)/float(grid(i))
        origin(i) = nint (ming(i)- minr/g(i))-1
        j = nint (maxg(i)+ maxr/g(i))+1
        extent(i) = j- origin(i)+1
      end do
c
c --- Now extend by some grid points around this origin and extent
c
      j = nextra
      ct = 1
      do i=1,3
        extent(i) = extent(i) + 2*j
        origin(i) = origin(i) - j
        ct = ct * extent(i)
      end do
c
c --- Report
c
c      call jvalut (' Mask origin :',3,origin)
c      call jvalut (' Mask extent :',3,extent)
c      call jvalut (' Grid points :',1,ct)
c      call jvalut (' Mask grid   :',3,grid)
c      call fvalut (' Mask cell   :',6,cell)
c
c --- DONE
c
      ierr = 0
c
      return
c
      end
c
c
c
      subroutine grpmsk (cell,grid,maxatm,xyzr,radius,
     +                   origin,extent,mask,size,doinit,ierr)
c
      implicit none
c
      integer maxatm
c
      integer i,j,k,l,extent(3),grid(3),origin(3)
      integer ct,i1,i2,i3,off,ierr,size,inx,nxy
      integer mask(size)
c
      real cell(6),g(3),distce,x(3),xp(3),xyzr(3,maxatm)
      real a(3,3),b(3,3),c(3,3),fake(6)
      real x1(3),x2(3),gmin
      real radius
c
      logical doinit
c
code ...
c
      ierr = -1
c
      do i=1,3
        g(i) = cell(i)/float(grid(i))
      end do
c
c --- Blank mask
c
      if (doinit) then
        ct = extent(3)*extent(2)*extent(1)
        do k=1,ct
          mask(k) = 0
        end do
      end if
c
      do i=1,3
        fake(i) = 1.0
        fake(i+3) = cell(i+3)
      end do
      call orthog (fake, a, 1)
      call orthog (fake, c, 0)
      call orthog (cell, b, 1)
c
c --- Do the actual mask job
c
      nxy = extent(1)*extent(2)
      gmin = min (g(1),g(2),g(3))
c
      do 120 i=1,maxatm
        do j=1,3
          x(j) = xyzr(j,i)
        end do
        call mulmtx (a, x, x1, 3, 3, 1)
        off = 1 + radius/gmin
        do 130 j=-off,off	
          do 130 k=-off,off	
            do 130 l=-off,off
              xp(1) = x1(1)+ l*g(1)
              xp(2) = x1(2)+ k*g(2)
              xp(3) = x1(3)+ j*g(3)
              call mulmtx (c, xp, x2, 3, 3, 1)
              if (distce(x,x2) .gt. radius) goto 130
              i1 = nint(xp(1)/g(1))- origin(1) +1
              i2 = nint(xp(2)/g(2))- origin(2) +1
              i3 = nint(xp(3)/g(3))- origin(3) +1
              if (i1 .le. 0) goto 130
              if (i2 .le. 0) goto 130
              if (i3 .le. 0) goto 130
              if (i1 .gt. extent(1)) goto 130
              if (i2 .gt. extent(2)) goto 130
              if (i3 .gt. extent(3)) goto 130
              inx = (i3-1)*nxy + (i2-1)*extent(1) + i1
              mask(inx) = 1
 130    continue
 125    continue
 120  continue
c
c --- DONE
c
      ierr = 0
c
      return
c
      end
c
c
c
      subroutine sub02x (
     +    mapa,exta1,exta2,exta3,orgna,spaca,cella,
     +    mapb,extb1,extb2,extb3,orgnb,spacb,cellb,
     +    maskb, rtbtoa, ctrt, rtsym, ctsym, or2or, ierr,
     +    avy,avysq,avxy,avxmy,lposit)
c
c ... average into a (possibly) different spacegroup
c
c ... Gerard Kleywegt, 930317,18
c     based on code by Alwyn Jones, 1991
c
      implicit none
c
      real small
      parameter (small=1.0e-9)
c
      integer orgna(3), exta1, exta2, exta3
      real mapa(exta1, exta2, exta3),spaca(3),cella(6)
c
      integer orgnb(3), extb1, extb2, extb3
      real mapb(extb1, extb2, extb3),spacb(3),cellb(6)
      integer maskb(extb1, extb2, extb3)
c
      real rtbtoa(12,*),rtsym(12,*),or2or(12),val1,f,cc
      real avy(*),avysq(*),avxy(*),avxmy(*),avx,avxsq,avax
      integer ctrt,ctsym,ierr,nbit,ngjk
c
      integer errcod, i, j, k, l, loop, m, ext(3), ijk(3)
c      integer bobo1,bobo2,bobo3
      integer nerr1,nerr2
c
      real forgna(3), fexta(3), mapbit(4,4,4)
      real value, x(3), gexta(3)
c
      real f2ca(3,3),c2fa(3,3),f2cb(3,3),c2fb(3,3)
      real fmsk(3),cmsk(3),cmska(3),ncmska(3),fcmska(3),xdum
c
      logical lposit
c
      equivalence (ijk(1),i)
      equivalence (ijk(2),j)
      equivalence (ijk(3),k)
c
code ...
c
      ierr = -1
      nbit = 0
c
      call orthog (cella, f2ca, 0)
      call orthog (cella, c2fa, 1)
      call orthog (cellb, f2cb, 0)
      call orthog (cellb, c2fb, 1)
c
      ext(1) = exta1
      ext(2) = exta2
      ext(3) = exta3
c
c ... get edges of map A in fractional coordinates
c
      do i=1,3
        forgna(i) = float(orgna(i))         *spaca(i)/cella(i)
        fexta(i)  = float(ext(i)+orgna(i)-1)*spaca(i)/cella(i)
        gexta(i)  = float(ext(i)+orgna(i)-2)*spaca(i)/cella(i)
      end do
c
c      call fvalut (' FORGNA :',3,forgna)
c      call fvalut (' FEXTA  :',3,fexta)
c      call fvalut (' GEXTA  :',3,gexta)
c
      avx   = 0.0
      avax  = 0.0
      avxsq = 0.0
      do loop=1,ctrt
        avy(loop)    = 0.0
        avysq(loop)  = 0.0
        avxy (loop)  = 0.0
        avxmy (loop) = 0.0
      end do
      ngjk = 0
      nerr1 = 0
      nerr2 = 0
c
c      call cntmsk (maskb,extb1,extb2,extb3,bobo1)
c      call jvalut (' Points in mask :',1,bobo1)
c      bobo2 = bobo1/10
c      bobo3 = bobo2
c
c ... loop over the mask points
c
      do 100 k=1,extb3
      do 100 j=1,extb2
      do 100 i=1,extb1
        mapb (i,j,k) = 0.0
        if (maskb(i,j,k) .ne. 1) goto 100
c
        ngjk = ngjk + 1
c
c        if (ngjk .eq. bobo3) then
c          xdum = 100.0 * float(ngjk-1)/float(bobo1)
c          call fvalut (' Progress (% mask) :',1,xdum)
c          bobo3 = bobo3 + bobo2
c        end if
c
c ... FMSK = fractional coords of mask point in B
c
        do l=1,3
          fmsk (l) = float(ijk(l)-1+orgnb(l))*spacb(l)/cellb(l)
        end do
c
c ... CMSK = cartesian coords of mask point in B
c     CMSKA = ditto, but now in A
c
        call mulmtx (f2cb,fmsk,cmsk,3,3,1)
        call vecrtv (cmsk,cmska,1,or2or(1),or2or(10))
c
c ... loop over the NCS operators in B
c
        do 120 loop=1,ctrt
c
c ... NCMSKA = ditto, in A, but in NCS mate
c     FCMSKA = ditto, but in fractional coords
c
          call vecrtv (cmska,ncmska,1,rtbtoa(1,loop),rtbtoa(10,loop))
          call mulmtx (c2fa,ncmska,fcmska,3,3,1)
c
c ... force this point into the asymmetric unit of map A
c
          call frcsym (fcmska,forgna,gexta,rtsym,ctsym,errcod)
c
c ... if not possible, it is close to the extent; build a 4x4x4 map
c
          if (errcod .ne. 0) then
c
            nbit = nbit + 1
            call bldbit (mapa,exta1,exta2,exta3,orgna,spaca,cella,
     +                   forgna,fexta,rtsym,ctsym,mapbit,fcmska,errcod)
c
            if (errcod .ne. 0) then
              if (nerr1 .lt. 10) then
                call errcon ('Serious FRCSYM error')
                call ivalut (' Mask point in reference  :',3,ijk)
                call fvalut (' Fractional               :',3,fmsk)
                call fvalut (' Cartesian                :',3,cmsk)
                call fvalut (' Cartesian in second xtal :',3,cmska)
                call fvalut (' Cartesian of NCS mate    :',3,ncmska)
                call fvalut (' Fractional               :',3,fcmska)
              else if (nerr1 .eq. 10) then
                call prompt (
     +            ' NOTE: further FRCSYM errors but not listed!!!')
              end if
              nerr1 = nerr1 + 1
              goto 120
            end if
c
c ... interpolate in 4x4x4 map
c
            do l=1,3
              x(l) = fcmska(l)*cella(l)/spaca(l)- float(orgna(l))+ 1
              m = x(l)
              x(l) = x(l)- float(m-1)
            end do
            call intrpl (mapbit, 4, 4, 4, x, value, errcod)
c
c ... if first FRCSYM worked okay, do straightforward interpolation
c
          else
c
            do l=1,3
              x(l) = fcmska(l)*cella(l)/spaca(l)- float(orgna(l)) + 1
            end do
c
            call intrpl (mapa, exta1, exta2, exta3, x, value, errcod)
c
          end if
c
          if (errcod .eq. 0) then
            mapb (i,j,k) = mapb (i,j,k) + value
c
              if (loop .eq. 1) then
                val1 = value
                avx = avx + value
                avxsq = avxsq + value*value
                avax  = avax + abs(value)
              end if
c
              avy (loop)   = avy(loop)   + value
              avysq (loop) = avysq(loop) + value*value
              avxy (loop)  = avxy(loop)  + value*val1
              avxmy (loop) = avxmy(loop) + abs(value-val1)
c
          else
            if (nerr2 .lt. 10) then
              call errcon ('Interpolation error')
              call ivalut (' Mask point in reference  :',3,ijk)
              call fvalut (' Fractional               :',3,fmsk)
              call fvalut (' Cartesian                :',3,cmsk)
              call fvalut (' Cartesian in second xtal :',3,cmska)
              call fvalut (' Cartesian of NCS mate    :',3,ncmska)
              call fvalut (' Fractional               :',3,fcmska)
            else if (nerr2 .eq. 10) then
              call prompt (
     +      ' NOTE: further interpolation errors but not listed!!!')
            end if
            nerr2 = nerr2 + 1
          end if
  120   continue
  100 continue
c
      ierr = 0
c
c      call jvalut (' Calls to BLDBIT      :',1,nbit)
c      call jvalut (' Severe FRCSYM errors :',1,nerr1)
c      call jvalut (' Interpolation errors :',1,nerr2)
c
c ... print correlation coefficients for the various operators
c
c      call jvalut (' Nr of mask points :',1,ngjk)
      f = float (ngjk)
      do loop=1,ctrt
        cc = (avxy(loop)/f- avx*avy(loop)/(f*f))/ 
     +    (sqrt(avxsq/f- (avx/f)**2) *
     +     sqrt(avysq(loop)/f- (avy(loop)/f)**2))
c        write (*,'(a,i3,a,f8.5)') ' Corr. coeff. for operator ',
c     +    loop,' = ',cc
        cc = avxmy(loop) / avax
c        write (*,'(a,i3,a,f8.5)') ' R-factor for operator ',
c     +    loop,' w.r.t. operator 1 = ',cc
      end do
c
c ... average (if necessary)
c
      if (lposit) then
        call rvalut (' Positivity; set <= 0 to:',1,small)
        xdum = 1.0 / float(ctrt)
        do 300 k=1,extb3
        do 300 j=1,extb2
        do 300 i=1,extb1
c ... changed next line for HP
          if (maskb(i,j,k) .eq. 1) then
            if (mapb(i,j,k) .gt. 0.0) then
              mapb(i,j,k) = mapb(i,j,k) * xdum
            else
              mapb(i,j,k) = small
            end if
          end if
  300   continue
      else
c        call prompt (' Averaging without positivity constraint')
        if (ctrt .gt. 1) then
          xdum = 1.0 / float(ctrt)
          do 200 k=1,extb3
          do 200 j=1,extb2
          do 200 i=1,extb1
200         mapb(i,j,k) = mapb(i,j,k) * xdum
        end if
      end if
c
      return
      end
c
c
c
      subroutine getmsk (maskit,no1,no2,no3,i1,i2,i3,iget,ierr)
c
      implicit none
c
      integer no1,no2,no3,i1,i2,i3,iget,ierr
      integer maskit(no1,no2,no3)
c
code ...
c
      ierr = -1
      if (i1 .lt. 1) return
      if (i2 .lt. 1) return
      if (i3 .lt. 1) return
c
      ierr = -2
      if (i1 .gt. no1) return
      if (i2 .gt. no2) return
      if (i3 .gt. no3) return
c
      iget = maskit (i1,i2,i3)
      ierr = 0
c
      return
      end
c
c
c
      subroutine wrimpi (file,iunit,map,nx,ny,nz,ierr)
c
c ... write MPI-style EM map
c
      implicit none
c
      integer nx,ny,nz
      real map(nx,ny,nz)
c
      real*4 emdata(40)
c
      integer iunit,ierr,i,j,k,ii,jj,kk,nb
c
      logical xinter
c
      character file*(*),coment*80
c
code ...
c
      ierr = -1
c
      call textut (' Opening MPI map :',file)
c
      call xopxub (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening file')
        return
      end if
c
      call prompt (' Writing data ...')
      call prompt (
     +  ' Using 4 bytes per pixel, and all EM data set to zero')
c
      nb = 4
      call stamp (coment)
      do i=1,40
        emdata (i) = 0.0
      end do
c
      write (iunit,err=9000) nb,nx,ny,nz,
     +  (coment(j:j),j=1,80),(emdata(k),k=1,40),
     +  (((map(kk,jj,ii),kk=1,nx),jj=1,ny),ii=1,nz)
c
      ierr = 0
      return
c
 9000 continue
      ierr = -1
      call errcon ('While writing MPI map')
c
      return
      end
