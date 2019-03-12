c
c ===========================================================================
c
      subroutine gkdcpu (total,user,sys)
c
      implicit none
c
c --- GKDCPU (...) => returns DELTA cpu time (since last call)
c
c --- G J Kleywegt @ 920403
c
      real dtime,total,user,sys,dum(2)
c
code ...
c
      total = dtime (dum)
      user = dum(1)
      sys  = dum(2)
c
      return
      end
