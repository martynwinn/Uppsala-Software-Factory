c
c ==========================================================================
c
      logical function xinter ()
c
c --- LOG FUNC XINTER () - returns .TRUE. if process is interactive
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
      xinter = linter
c
      return
      end
