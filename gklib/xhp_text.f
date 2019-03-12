c
      SUBROUTINE xhp_text(WHAT,NUMBA)
c
C ---	DRAW TEXT AT CURRENT PEN POSITION
C ---	WHAT CONTAINS NUMBA ASCII CHARACTERS
c
      include 'xhp_graf.incl'
c
      integer numba
c
      BYTE WHAT(NUMBA)
      BYTE P(3)
C
! mrh
      character*1 P_c(3)
      equivalence(P,P_c)
c
      DATA P_c/'L','B','3'/
C

      WRITE(LUN,10)P(1),P(2),WHAT,P(3)
      RETURN
c
  10  FORMAT(1X,80A1)
c
      END
