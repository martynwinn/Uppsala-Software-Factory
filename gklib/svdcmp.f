c
c ========================================================================
c
      subroutine svdcmp (A,M,N,MP,NP,W,V)
C
C =========================================================================
C
C SUBROUTINE : SVDCMP
C -------------------
C
C PURPOSE : given matrix A (log dim M by N; phys dim MP by NP) this
C           routine computes its Singular Value Decomposition,
C           A = U W V(Transpose). U replaces A on output. The diagonal
C           matrix of singular values W is output as a vector W. The
C           matrix V (NOT its transpose !) is output as V. M > =  N; if
C           not, then A should be filled up with zeros until M = N.
C
C CALL    : SVDCMP (A,M,N,MP,NP,W,V)
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
      real zero,one
      PARAMETER (ZERO = 0.0,ONE = 1.0)
c
      integer m,n,mp,np,i,l,k,j,nm,its
c
      REAL A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      real g,scale,anorm,s,f,h,c,x,y,z
c
code ...
c
      G = ZERO
      SCALE = ZERO
      ANORM = ZERO
      DO 25 I = 1,N
        L = I+1
        RV1(I) = SCALE*G
        G = ZERO
        S = ZERO
        SCALE = ZERO
        IF (I.LE.M) THEN
          DO 11 K = I,M
            SCALE = SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.ZERO) THEN
            DO 12 K = I,M
              A(K,I) = A(K,I)/SCALE
              S = S+A(K,I)*A(K,I)
12          CONTINUE
            F = A(I,I)
            G = -SIGN(SQRT(S),F)
            H = F*G-S
            A(I,I) = F-G
            IF (I.NE.N) THEN
              DO 15 J = L,N
                S = ZERO
                DO 13 K = I,M
                  S = S+A(K,I)*A(K,J)
13              CONTINUE
                F = S/H
                DO 14 K = I,M
                  A(K,J) = A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K =  I,M
              A(K,I) = SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I) = SCALE *G
        G = ZERO
        S = ZERO
        SCALE = ZERO
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K = L,N
            SCALE = SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.ZERO) THEN
            DO 18 K = L,N
              A(I,K) = A(I,K)/SCALE
              S = S+A(I,K)*A(I,K)
18          CONTINUE
            F = A(I,L)
            G = -SIGN(SQRT(S),F)
            H = F*G-S
            A(I,L) = F-G
            DO 19 K = L,N
              RV1(K) = A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J = L,M
                S = ZERO
                DO 21 K = L,N
                  S = S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K = L,N
                  A(J,K) = A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K = L,N
              A(I,K) = SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM = MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I = N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.ZERO) THEN
            DO 26 J = L,N
              V(J,I) = (A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J = L,N
              S = ZERO
              DO 27 K = L,N
                S = S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K = L,N
                V(K,J) = V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J = L,N
            V(I,J) = ZERO
            V(J,I) = ZERO
31        CONTINUE
        ENDIF
        V(I,I) = ONE
        G = RV1(I)
        L = I
32    CONTINUE
      DO 39 I = N,1,-1
        L = I+1
        G = W(I)
        IF (I.LT.N) THEN
          DO 33 J = L,N
            A(I,J) = ZERO
33        CONTINUE
        ENDIF
        IF (G.NE.ZERO) THEN
          G = ONE/G
          IF (I.NE.N) THEN
            DO 36 J = L,N
              S = ZERO
              DO 34 K = L,M
                S = S+A(K,I)*A(K,J)
34            CONTINUE
              F = (S/A(I,I))*G
              DO 35 K = I,M
                A(K,J) = A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J = I,M
            A(J,I) = A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J =  I,M
            A(J,I) = ZERO
38        CONTINUE
        ENDIF
        A(I,I) = A(I,I)+ONE
39    CONTINUE
      DO 49 K = N,1,-1
        DO 48 ITS = 1,50
          DO 41 L = K,1,-1
            NM = L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C = ZERO
          S = ONE
          DO 43 I = L,K
            F = S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G = W(I)
              H = SQRT(F*F+G*G)
              W(I) = H
              H = ONE/H
              C =  (G*H)
              S = -(F*H)
              DO 42 J = 1,M
                Y = A(J,NM)
                Z = A(J,I)
                A(J,NM) = (Y*C)+(Z*S)
                A(J,I) = -(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z = W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.ZERO) THEN
              W(K) = -Z
              DO 44 J = 1,N
                V(J,K) = -V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.50) then
            call errcon ('SVDCMP - No convergence (50 iters)')
            return
          end if
          X = W(L)
          NM = K-1
          Y = W(NM)
          G = RV1(NM)
          H = RV1(K)
          F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G = SQRT(F*F+ONE)
          F = ((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C = ONE
          S = ONE
          DO 47 J = L,NM
            I = J+1
            G = RV1(I)
            Y = W(I)
            H = S*G
            G = C*G
            Z = SQRT(F*F+H*H)
            RV1(J) = Z
            C = F/Z
            S = H/Z
            F =  (X*C)+(G*S)
            G = -(X*S)+(G*C)
            H = Y*S
            Y = Y*C
            DO 45 NM = 1,N
              X = V(NM,J)
              Z = V(NM,I)
              V(NM,J) =  (X*C)+(Z*S)
              V(NM,I) = -(X*S)+(Z*C)
45          CONTINUE
            Z = SQRT(F*F+H*H)
            W(J) = Z
            IF (Z.NE.ZERO) THEN
              Z = ONE/Z
              C = F*Z
              S = H*Z
            ENDIF
            F =  (C*G)+(S*Y)
            X = -(S*G)+(C*Y)
            DO 46 NM = 1,M
              Y = A(NM,J)
              Z = A(NM,I)
              A(NM,J) =  (Y*C)+(Z*S)
              A(NM,I) = -(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L) = ZERO
          RV1(K) = F
          W(K) = X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
