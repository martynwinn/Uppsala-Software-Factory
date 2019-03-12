c
c ===========================================================================
c
      subroutine gkecpu (total,user,sys)
c
      implicit none
c
c --- GKECPU (...) => returns elapsed cpu time
c
c --- G J Kleywegt @ 920403
c
      real etime,total,user,sys,dum(2)
c
code ...
c
      total = etime (dum)
      user = dum(1)
      sys  = dum(2)
c
      return
      end
