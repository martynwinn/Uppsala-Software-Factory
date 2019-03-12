c
c ===========================================================================
c
      subroutine errcon (why)
c
      implicit none
c
      integer lr,length
c
      logical xsocket
c
      character why*(*),esline*128
c
code ...
c
      lr = length (why)
      if (lr .gt. 0) then
        if (xsocket()) then
          esline = (' ERROR --- '//why(1:lr))
          call oprint (6,esline)
        else
          write (*,1000) why(1:lr)
        end if
      else
        if (xsocket()) then
          call oprint (6,'***UNKNOWN ERROR***')
        else
          write (*,1000) '***UNKNOWN ERROR***'
        end if
      end if
c
 1000 format (' ERROR --- ',a)
c
      return
      end
