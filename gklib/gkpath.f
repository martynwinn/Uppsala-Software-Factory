c
c
c
      subroutine gkpath ()
c
c ... check if environment variable GKPATH is defined; if so, store
c
      implicit none
c
      integer ierr,length,i,l
c
      character value*1000
c
      integer ndir,maxdir,dirlen
      parameter (maxdir=50,dirlen=100)
      character patdir(maxdir)*(dirlen)
c
      common /gpath1/ ndir
      common /gpath2/ patdir
c
      character separ*1,dirsep
      data separ /':'/, dirsep /'/'/
c
code ...
c
      ndir = 0
c
      call gknval ('GKPATH',value,ierr)
      if (ierr .ne. 0) return
      if (length(value) .lt. 1) return
      if (value(1:1) .eq. separ) value = value(2:)
      i = length(value)
      if (value(i:i) .eq. separ) value(i:) = ' '
c
ccc      call textut (' GKPATH :',value)
c
   10 continue
      if (ndir .eq. maxdir) return
c
      i = index (value,separ)
c
      if (i .lt. 1) then
        if (ndir .lt. maxdir) then
          if (length(value) .le. dirlen) then
            ndir = ndir + 1
            patdir (ndir) = value
            call remspa (patdir(ndir))
            l = length(patdir(ndir))
            if (patdir(ndir)(l:l) .ne. dirsep) then
              patdir(ndir)(l+1:l+1) = dirsep
            end if
c
ccc      call textut (' PATH :',patdir(ndir))
c
          end if
        end if
        return
      end if
c
      if (ndir .lt. maxdir) then
        if (length(value) .le. dirlen) then
          ndir = ndir + 1
          patdir (ndir) = value (1:i-1)
          call remspa (patdir(ndir))
          l = length(patdir(ndir))
          if (patdir(ndir)(l:l) .ne. dirsep) then
            patdir(ndir)(l+1:l+1) = dirsep
          end if
c
ccc      call textut (' PATH :',patdir(ndir))
c
        end if
      end if
c
      value = value (i+1:)
      goto 10
c
      end
