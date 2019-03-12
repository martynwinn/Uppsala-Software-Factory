c
c ==========================================================================
c
      logical function xbatch ()
c
c --- LOG FUNC XBATCH () - returns .TRUE. if process is batch
c
c --- G J Kleywegt @ 920915
c
      implicit none
c
      logical linter,lbatch,lsocket
c
      common /promod/ linter,lbatch,lsocket
c
code ...
c
      xbatch = lbatch
c
      return
      end
