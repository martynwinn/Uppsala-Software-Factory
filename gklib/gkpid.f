c
c ==========================================================================
c
      subroutine gkpid (id)
c
c --- GKPID (id) - get process id
c
c --- G J Kleywegt @ 930320
c
      implicit none
c
      integer id,getpid
c
code ...
c
      id = getpid ()
c
      return
      end
