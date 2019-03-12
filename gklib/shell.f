c
c ===========================================================================
c
      subroutine shell (A,INX,N)
c
C ... SHELL SORT ALGORITHM
c
c ... A contains N INTEGER values by which the index array INX
c     should be sorted; array A itself is also sorted, so if you don't
c     this to happen, use a scratch array !!!
c
      implicit none
c
      INTEGER n,A(N),inx(n),W,i,m,k,j,ipm
c
code ...
c
      I=1
  10  CONTINUE
      IF (I.GT.N) GOTO 20
      M=2*I-1
      I=2*I
      GOTO 10
  20  M=M/2
      IF (M.EQ.0) RETURN
      K=N-M
      DO 40 J=1,K
        I=J
  50    CONTINUE
        IPM=I+M
        IF (A(IPM).GE.A(I)) GOTO 40
c
c ... swap As and INXs
c
        W=A(I)
        A(I)=A(IPM)
        A(IPM)=W
c
        w=inx(i)
        inx(i)=inx(ipm)
        inx(ipm)=w
c
        I=I-M
        IF (I.GE.1) GOTO 50
  40  CONTINUE
      GOTO 20
c
      END
