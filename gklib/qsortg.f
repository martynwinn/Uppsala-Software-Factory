c
c ========================================================================
c
      subroutine qsortg (x,n)
C
C =========================================================================
C
C SUBROUTINE : QSORTG
C -------------------
C
C PURPOSE : sorting an array X
C
C CALL    : QSORTG (X,N)
C
C REFER.  : ???
C
C =========================================================================
C
      implicit NONE
c
      REAL X(*),xmax
c
      integer n,i,j,imax,j2
c
code ...
c
      IF (N.LE.1) RETURN
C
      J = 0
    2 J = J+1
C GJK@12JAN89 : the following line included !!! (namely for xmax = x(j) !!)
      IMAX = J
      XMAX = X(J)
      J2 = J+1
C GJK@19DEC88 : replaced .EQ. by .GT. in the following line !!!
      IF (J2.GT.N) goto 3
      DO 10 I = J2,N
        IF (X(I).LT.XMAX) GOTO 10
        IMAX = I
        XMAX = X(I)
   10 CONTINUE
      X(IMAX) = X(J)
      X(J) = XMAX
      GOTO 2
c
c ... sort from low to high
c
    3 continue
      do i=1,n/2
        xmax = x(i)
        x(i) = x(n-i+1)
        x(n-i+1) = xmax
      end do
c
      return
      end
