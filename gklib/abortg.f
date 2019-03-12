c
c ===========================================================================
c
      subroutine abortg (why)
c
c ... fatal error occured; print message and bail out
c
      implicit none
c
      integer lr,length
c
      logical xsocket
c
      character why*(*),esline*128
c
      entry endit  (why)
c
      entry fatal  (why)
c
      entry errstp (why)
c
code ...
c
      lr = length (why)
      if (lr .gt. 0) then
        if (xsocket()) then
          esline = (' FATAL ERROR --- '//why(1:lr)//
     +      ' - execution aborted !')
          call oprint (6,esline)
        else
          write (*,1000) why(1:lr)
        end if
      else
        if (xsocket()) then
          call oprint (6,'***UNKNOWN FATAL ERROR***')
        else
          write (*,1000) '***UNKNOWN FATAL ERROR***'
        end if
      end if
c
      call gkquit
c
 1000 format (' FATAL ERROR --- ',a,' - execution aborted !')
c
      return
      end
