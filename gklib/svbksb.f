c
c ========================================================================
c
      subroutine svbksb (U,W,V,M,N,MP,NP,B,X)
C
C =========================================================================
C
C SUBROUTINE : SVBKSB
C -------------------
C
C PURPOSE : solves AX = B for a vector X, where a is specified by the arrays
C           U, W, V as returned by SVDCMP. M and N are the logical dimensions
C           of A and will be equal for square matrices. MP and NP are the
C           physical dimensions of A. B is the input right-hand side. X is
C           the output solution vector. No input quantities are destroyed,
C           so the routine may be called repeatedly with different B's
C
C CALL    : SVBKSB (U,W,V,M,N,MP,NP,B,X)
C
C REFER.  : W H Press, B P Flannery, S A Teukolsky & W T Vetterling,
C           "Numerical Recipes - The Art of Scientific Computing",
C           Cambridge University Press, Cambridge, 1986, pp. 52-64
C
C =========================================================================
C
      implicit none
c
      integer nmax
      PARAMETER (NMAX = 100)
c
      real zero
      PARAMETER (ZERO = 0.0)
c
      integer mp,np,m,n,i,j,jj
c
      REAL U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      real s
c
code ...
c
      DO 12 J = 1,N
        S = ZERO
        IF(W(J).NE.ZERO)THEN
          DO 11 I = 1,M
            S = S+U(I,J)*B(I)
11        CONTINUE
          S = S/W(J)
        ENDIF
        TMP(J) = S
12    CONTINUE
      DO 14 J = 1,N
        S = ZERO
        DO 13 JJ = 1,N
          S = S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J) = S
14    CONTINUE
      RETURN
      END
