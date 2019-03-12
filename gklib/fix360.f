c
c=====================================================
c
      subroutine fix360 (ang)
c
c ... fix angle in range [0.0,+360.0>
c
      implicit none
c
      real ang
c
code ...
c
  210 if (ang .lt. 0.0) then
        ang = ang + 360.0
        goto 210
      end if
c
  220 if (ang .ge. 360.0) then
        ang = ang - 360.0
        goto 220
      end if
c
      return
      end
