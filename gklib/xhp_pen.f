c
      SUBROUTINE xhp_pen(I,J)
c
C ---	I=1,PEN UP
C ---	I=2 PEN DOWN
C ---	I=3 CHANGE COLOUR TO J
c
      include 'xhp_graf.incl'
c
      integer i,j
c
! mrh
      BYTE P(3),Q(3),R(3)
      character*1 P_c(3), Q_c(3),R_c(3)
      equivalence(P,P_c), (Q,Q_c),(R, R_c)
c
      DATA P_c/'P','U',';'/,
     *     Q_c/'P','D',';'/,
     *     R_c/'S','P',';'/
C
      GOTO (1,2,3),I
c
1     WRITE(LUN,100)P
      RETURN
c
2     WRITE(LUN,100)Q
      RETURN
c
3     WRITE(LUN,200)R(1),R(2),J,R(3)
      RETURN
c
  100 FORMAT(1X,3A1)
  200 FORMAT(1X,2A1,I1,A1)
c
      END
