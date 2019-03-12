c
c ==========================================================================
c
      logical function xsocket ()
c
c --- LOG FUNC XSOCKET () - returns .TRUE. if process is socket
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
      xsocket = lsocket
c
      return
      end
