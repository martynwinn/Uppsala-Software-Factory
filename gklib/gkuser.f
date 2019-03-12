c
c ===========================================================================
c
      subroutine gkuser (usernm)
c
      implicit none
c
c --- GKUSER (...) => returns user's login name
c
c --- G J Kleywegt @ 920403
c
      character usernm*(*)
c
code ...
c
      call getlog (usernm)
c
      return
      end
