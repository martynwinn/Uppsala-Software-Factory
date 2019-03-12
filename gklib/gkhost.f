c
c ==========================================================================
c
      subroutine gkhost (string)
c
c --- GKHOST (string) - get hostname
c
c --- G J Kleywegt @ 930320
c
      implicit none
c
      integer namlen,i,gethostname
c
      character string*(*)
      character*1 name(64)
c
code ...
c
      string = 'localhost (OSX/Apple)'
c
c      string = ' '
c      namlen = 64
c      i = gethostname (name,namlen)
c      if (i .eq. 0) then
c        do i=1,namlen
c          string (i:i) = name(i)
c        end do
c      end if
c
      return
      end
