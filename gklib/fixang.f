c
c=====================================================
c
      subroutine fixang (ang)
c
c ... fix angle in range <-180.0,+180.0]
c
      implicit none
c
      real ang
c
code ...
c
  210 if (ang .le. -180.0) then
        ang = ang + 360.0
        goto 210
      end if
c
  220 if (ang .gt. 180.0) then
        ang = ang - 360.0
        goto 220
      end if
c
      return
      end
