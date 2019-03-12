c
c ===========================================================================
c
c OSX_SYSUBS.F - system-dependent routines - Gerard Kleywegt - Uppsala - 2002
c             - this version for OSX (Apple) - but identical to LINUX version
c
c ===========================================================================
c
c subroutine flusho (iunit)
c subroutine spawno (command)
c subroutine gkdate (str24)
c subroutine gkecpu (total,user,sys)
c subroutine gkdcpu (total,user,sys)
c subroutine gknval (nam,val,ierr)
c subroutine gkuser (usernm)
c subroutine gkatty (iunit,ltty,ntty)
c subroutine defina (myname,fpath,f1,fname,f2,fext,f3,fvers)
c subroutine gkrand (xrand,xlo,xhi)
c subroutine gkargs ()
c subroutine gkerr  (string)
c subroutine gkhost (string)
c subroutine gkpid  (id)
c subroutine gksys  (string)
c subroutine gklibf (file)
c
c ===========================================================================
c
      subroutine flusho (iunit)
c
      implicit none
c
c === FLUSHO (IUNIT) => flush output buffer on IUNIT
c
c !!! If not required or not available, comment out the
c !!! call to FLUSH (Unix routine)
c
c === G J Kleywegt @ 920331
c
      integer iunit
c
code ...
c
      call flush (iunit)
c
      return
      end
