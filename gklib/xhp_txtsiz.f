c
      SUBROUTINE xhp_txtsiz(ISIZE)
c
C ---	DEFINE TEXTSIZE AS MULTIPLE OF 1/4 STANDARD
c
      include 'xhp_graf.incl'
c
      real x,y
c
      integer isize
c
      BYTE P(3)
c
! mrh
      character*1 P_c(3)
      equivalence(P,P_c)
c
      DATA P_c/'S','I',';'/
C

C      DATA P/'S','I',';'/
c
      X = 0.19*FLOAT(ISIZE)*0.25
      Y = 0.27*FLOAT(ISIZE)*0.25
      WRITE(LUN,10)P(1),P(2),X,Y,P(3)
c
      RETURN
c
  10  FORMAT(1X,2A1,2F5.2,A1)
c
      END
