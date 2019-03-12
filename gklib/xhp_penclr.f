c
      SUBROUTINE xhp_penclr(J)
c
C ---	CHANGE PEN COLOUR,DO NOTHING IF INVALID COLOUR
C ---	J=0,NO PEN.J=1,PEN 1ETC.J>4INVALID
c
      implicit none
c
      integer j
c
      IF(J.LT.0)RETURN
      IF(J.GT.4)RETURN
      CALL xhp_PEN(3,J)
c
      RETURN
      END
