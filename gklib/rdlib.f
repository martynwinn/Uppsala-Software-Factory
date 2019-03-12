c
c
c
      subroutine rdlib (iunit,maxtyp,typ1lc,cmpmat,ierr)
c
c ... RDLIB = read substitution matrix file
c
      implicit none
c
      integer maxaat
      parameter (maxaat=30)
c
      integer maxtyp,iunit,ierr,ntyp,i,j,length
c
      real cmpmat(maxtyp,maxtyp),rmat(maxaat,maxaat)
      real avemat
c
      integer imat(maxaat,maxaat)
      integer iptr(maxaat)
c
      logical lmatok
c
      character typ1lc(maxtyp)*1,mytype(maxaat)*1
      character line*256,myform*40
c
code ...
c
      ierr = -1
      ntyp = -1
      lmatok = .false.
c
   10 continue
      read (iunit,'(a)',err=100,end=200) line
      if (line(1:1) .eq. '!' .or. line(1:1) .eq. '#') then
        if (length(line) .gt. 1) call textut (' Comment :',line)
        goto 10
      end if
c
      call upcase (line)
ccc      print*,line(1:leng1(line))
      if (line(1:4) .eq. 'TYPE') then
        read (line(5:),'(i3,1x,a)',err=100,end=100) ntyp,myform
        call remspa (myform)
ccc        print *,'NTYP, MYFORM ',ntyp,' |',myform,'|'
        read (iunit,myform,err=100,end=100) (mytype(i),i=1,ntyp)
        goto 10
      end if
c
      if (line(1:4) .eq. 'MATI') then
        if (ntyp .le. 0) then
          call errcon ('MATI before TYPE')
          return
        end if
        read (line(5:),'(a)',err=100,end=100) myform
        call remspa (myform)
        call textut (' Read INTR matrix with format :',myform)
ccc        print *,'MYFORM |',myform,'|'
        do i=1,ntyp
          read (iunit,myform,err=100,end=100) (imat(j,i),j=1,ntyp)
          do j=1,ntyp
            rmat (j,i) = float(imat(j,i))
          end do
        end do
        lmatok = .true.
        goto 10
      end if
c
      if (line(1:4) .eq. 'MATR') then
        if (ntyp .le. 0) then
          call errcon ('MATR before TYPE')
          return
        end if
        read (line(5:),'(a)',err=100,end=100) myform
        call remspa (myform)
        call textut (' Read REAL matrix with format :',myform)
ccc        print *,'MYFORM |',myform,'|'
        do i=1,ntyp
          read (iunit,myform,err=100,end=100) (rmat(j,i),j=1,ntyp)
        end do
        lmatok = .true.
        goto 10
      end if
c
      goto 10
c
  100 continue
      call errcon ('While reading library file')
      return
c
  200 continue
c
      if (.not .lmatok) then
        call errcon ('No matrix found')
        return
      end if
c
      do i=1,maxtyp
        do j=1,maxtyp
          cmpmat (j,i) = -9.99
        end do
        iptr (i) = -1
        do j=1,ntyp
          if (mytype(j) .eq. typ1lc (i)) then
            iptr (i) = j
            goto 300
          end if
        end do
        call textut (' Residue type not found :',typ1lc(i))
  300   continue
      end do
c
      do i=1,maxtyp
        if (iptr(i) .gt. 0) then
          do j=1,maxtyp
            if (iptr(j) .gt. 0) then
              cmpmat (j,i) = rmat (j,i)
            end if
          end do
ccc          call fvalut (typ1lc(i),maxtyp,cmpmat(1,i))
        end if
      end do
c
      avemat = 0.0
      do i=1,maxtyp
        do j=1,maxtyp
          avemat = avemat + cmpmat (i,j)
          if (abs(cmpmat(i,j)-cmpmat(j,i)) .gt. 0.01) then
            call errcon ('Matrix is asymmetric')
            call textut (' Amino acid I :',typ1lc(i))
            call textut (' Amino acid J :',typ1lc(j))
            call fvalut (' Matrix (I,J) :',1,cmpmat(i,j))
            call fvalut (' Matrix (J,I) :',1,cmpmat(j,i))
            return
          end if
        end do
      end do
c
      avemat = avemat / float(maxtyp*maxtyp)
      call fvalut (' Average matrix value :',1,avemat)
c
      ierr = 0
c
      return
      end
