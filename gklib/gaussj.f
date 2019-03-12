c
c ========================================================================
c
      subroutine  gaussj (A,N,NP,B,M,MP)
C
C =========================================================================
C
C SUBROUTINE : GAUSSJ
C -------------------
C
C PURPOSE : solving a set of equations AND inverting the matrix A
C
C CALL    : GAUSSJ (A,N,NP,B,M,MP)
C           A  = matrix of equations (N*N); physical dim NP*NP
C           B  = matrix of right-hand side vectors (M of them, of length N);
C                physical dim NP*MP
C           On output, A is replaced by its inverse and B by the matrix
C           of corresponding solution vectors
C
C REFER.  : W H Press, B P Flannery, S A Teukolsky & W T Vetterling,
C           "Numerical Recipes - The Art of Scientific Computing",
C           Cambridge University Press, Cambridge, 1986, pp. 24-29
C
C =========================================================================
C
      implicit none
c
      integer nmax
      PARAMETER (NMAX = 50)
c
      integer n,np,m,mp,i,j,k,l,icol,irow,ll
      INTEGER IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
c
      REAL A(NP,NP),B(NP,MP)
      real big,dum,pivinv
c
code ...
c
      DO 11 J = 1,N
        IPIV(J) = 0
11    CONTINUE
      DO 22 I = 1,N
        BIG = 0.
        DO 13 J = 1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K = 1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG = ABS(A(J,K))
                  IROW = J
                  ICOL = K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                call errcon ('GAUSSJ - Singular matrix')
                return
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL) = IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L = 1,N
            DUM = A(IROW,L)
            A(IROW,L) = A(ICOL,L)
            A(ICOL,L) = DUM
14        CONTINUE
          DO 15 L = 1,M
            DUM = B(IROW,L)
            B(IROW,L) = B(ICOL,L)
            B(ICOL,L) = DUM
15        CONTINUE
        ENDIF
        INDXR(I) = IROW
        INDXC(I) = ICOL
        IF (A(ICOL,ICOL).EQ.0.) THEN
          call errcon ('GAUSSJ - Singular matrix')
          return
        ENDIF
        PIVINV = 1./A(ICOL,ICOL)
        A(ICOL,ICOL) = 1.
        DO 16 L = 1,N
          A(ICOL,L) = A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L = 1,M
          B(ICOL,L) = B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL = 1,N
          IF(LL.NE.ICOL)THEN
            DUM = A(LL,ICOL)
            A(LL,ICOL) = 0.
            DO 18 L = 1,N
              A(LL,L) = A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L = 1,M
              B(LL,L) = B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L = N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K = 1,N
            DUM = A(K,INDXR(L))
            A(K,INDXR(L)) = A(K,INDXC(L))
            A(K,INDXC(L)) = DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
