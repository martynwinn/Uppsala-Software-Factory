c
c ========================================================================
c
      subroutine isort (n,iar)
C
C =========================================================================
C
C SUBROUTINE : ISORT
C ------------------
C
C PURPOSE : sorting an array IAR (using straight insertion); this routine
C           is only efficient for small N (say, N < 50)
C
C CALL    : ISORT (N,IAR)
C
C REFER.  : W H Press, B P Flannery, S A Teukolsky & W T Vetterling,
C           "Numerical Recipes - The Art of Scientific Computing",
C           Cambridge University Press, Cambridge, 1986, pp. 227-228
C
C =========================================================================
C
      implicit none
c
      INTEGER n,IAR(N),i,j,ia
C
      IF (N.le.1) GOTO 12
C
      DO J = 2,N
        IA = IAR (J)
        DO 11 I = J-1,1,-1
          IF (IAR(I).LE.IA) GOTO 10
          IAR (I+1) = IAR (I)
   11   CONTINUE
        I = 0
   10   IAR (I+1) = IA
      end do
c
   12 CONTINUE
C
      RETURN
      END
