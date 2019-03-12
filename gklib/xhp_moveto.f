c
      SUBROUTINE xhp_moveto(X,Y)
c
C ---	MOVE TO X,Y ABSOLUT,BUT PICK UP PEN FIRST
C
C ---	X,Y ARE IN USER SPACE
C
      include 'xhp_graf.incl'
c
      real x,y
c
      integer i,j
c
      BYTE P(3)
! mrh
      character*1 P_c(3)
      equivalence(P,P_c)
c
      DATA P_c/'P','A',';'/
c
      CALL xhp_PEN(1)
      I=X*XSCALE(1)+XSCALE(2)
      J=Y*YSCALE(1)+YSCALE(2)
      WRITE(LUN,10)P(1),P(2),I,J,P(3)
      RETURN
c
  10  FORMAT(1X,2A1,2I6,A1)
c
      END
