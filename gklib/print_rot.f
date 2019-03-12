c
c ==========================================================================
c
      subroutine print_rot (mat)
c
      implicit NONE
c
      real mat(9)
      real det3
c
code ...
c
      if (abs(det3(mat)) .le. 2.0) then
        write (*,6000) 'X',mat(1),mat(4),mat(7)
        write (*,6000) 'Y',mat(2),mat(5),mat(8)
        write (*,6000) 'Z',mat(3),mat(6),mat(9)
      else
        write (*,7000) 'X',mat(1),mat(4),mat(7)
        write (*,7000) 'Y',mat(2),mat(5),mat(8)
        write (*,7000) 'Z',mat(3),mat(6),mat(9)
      end if
c
 6000 format (' New ',a1,' = ',f10.7,' * X + ',f10.7,' * Y + ',
     +  f10.7,' * Z')
 7000 format (' New ',a1,' = ',f10.4,' * X + ',f10.4,' * Y + ',
     +  f10.4,' * Z')
c
      return
      end
