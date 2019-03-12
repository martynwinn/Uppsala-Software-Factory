c
c
c
      subroutine xopxoa (iunit,file,lretry,ierr)
c
c ... try to open existing ASCII file using GKPATH env. var.
c
      implicit none
c
      integer iunit,ierr,length,i,l1,l2
c
      logical lretry
c
      character file*(*),myfile*500
c
      integer ndir,maxdir,dirlen
      parameter (maxdir=50,dirlen=100)
      character patdir(maxdir)*(dirlen)
c
      common /gpath1/ ndir
      common /gpath2/ patdir
c
      character separ*1
      data separ /'/'/
c
code ...
c
      ierr = -1
      if (length(file) .lt. 1) then
        call errcon ('XOPXOA - No file name provided')
        return
      end if
c
      if (iunit .le. 0) then
        call errcon ('XOPXOA - Zero or negative unit number')
        return
      end if
c
c ... (1) - if file name contains directory separator (/), or if
c           GKPATH empty, then don't use GKPATH
c
      if (index(file,separ) .gt. 0 .or. ndir .le. 0) then
ccc      print *,' CALL XOLASC'
        call xolasc (iunit,file,lretry,ierr)
        return
      end if
c
c ... (2) - see if normal file exists
c
      call xexist (iunit,file,lretry,ierr)
      if (ierr .eq. 0) return
c
c ... (3) - try GKPATH directories
c
      l1 = length (file)
      do i=1,ndir
        l2 = length(patdir(i))
        if (l1+l2 .le. len(myfile)) then
          myfile = patdir(i)(1:l2)//file(1:l1)
          call xexist (iunit,myfile,lretry,ierr)
          if (ierr .eq. 0) then
            call textut (' ==> Found file in GKPATH :',myfile)
            return
          end if
        end if
      end do
c
c ... (4) - nothing worked
c
      call errcon ('XOPXOA - File not found in GKPATH')
      call xolasc (iunit,file,lretry,ierr)
c
      return
      end
