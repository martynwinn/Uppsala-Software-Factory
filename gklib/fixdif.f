c
c=====================================================
c
      subroutine fixdif (ang1,ang2,dif)
c
c ... fix difference of two angles in range <-180.0,+180.0]
c
      implicit none
c
      real ang1,ang2,dif
c
code ...
c
      dif = ang1 - ang2
c
  210 if (dif .le. -180.0) then
        dif = dif + 360.0
        goto 210
      end if
c
  220 if (dif .gt. 180.0) then
        dif = dif - 360.0
        goto 220
      end if
c
      return
      end
