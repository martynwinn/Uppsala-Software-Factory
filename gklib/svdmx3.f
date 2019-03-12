c
c ========================================================================
c
      subroutine svdmx3 (F,O,OEXTR,S,FEXTR,init)
c
ccc,NOK,NERR1,NERR2,NUND)
C
C =========================================================================
C
C     SUBROUTINE SVDMX3 (F,O,OEXTR,S,FEXTR,INIT,NOK,NERR1,NERR2,NUND)
C     ================= get precise 3D peak position through interpolation
C
C =========================================================================
C
C Gerard J Kleywegt - Dept of Organic Chemistry - University of Utrecht
C
C 12AUG88      Latest update @ 09JAN89/900510/930303
C
C ======================================================================
C
C F (27)      function values of the 27 points
C O (3)       coordinates of the local maximum
C OEXTR       function value at the local maximum
C S (3)       coordinates of the interpolated maximum
C FEXTR       function value at this point
C INIT        .FALSE. on first call, .TRUE. thereafter
C
C ======================================================================
C
C fit a curve :
C ... f =     c1*x    + c2*y    + c3*z +
C ...         c4*x**2 + c5*y**2 + c6*z**2 +
C ...         c7*x*y  + c8*x*z  + c9*y*z
C through 27 points (-1,-1,-1), (-1,-1,0), ..., (1,1,1) whose function
C values are stored in vector F (the value at the central point 0,0,0
C will be subtracted from this in order to reduce the number of unknowns
C by one (9 instead of 10) and to ensure that the hypersurface will
C pass through this point)
C
C the coordinates of the maximum are returned in vector S
C the maximum function value is returned in FEXTR
C
C uses SVD-routines from WH Press et al., "Numerical Recipes"
C uses Gauss-Jordan elimination routines from the same source
C
C a second order polynomial f(x,y,z) is fitted through 27
C points with respective (differential) x,y,z-coordinates of:
C (-1,-1,-1), (-1,-1,0), (-1,-1,1), (-1,0,-1), ...., (1,1,1)
C array F contains the function values at these points
C matrix C contains the coefficients (fixed for fixed x,y,z-values) :
C   C(m,1) = x
C   C(m,2) = y
C   C(m,3) = z
C   C(m,4) = x**2
C   C(m,5) = y**2
C   C(m,6) = z**2
C   C(m,7) = x*y
C   C(m,8) = x*z
C   C(m,9) = y*z
C (m = 1,..,27, with m = 9*x + 3*y + z + 14)
C
C array X will finally contain the coefficients
C the fitting polynomial is then : f(x,y,z)-f(0,0,0) = X(1)*x + X(2)*y +
C ...    X(3)*z + X(4)*x**2 + X(5)*y**2 + X(6)*z**2 + X(7)*x*y +
C ...    X(8)*x*z + X(9)*y*z
C
C by taking the derivatives of F with reSPECt to x,y and z and setting
C dF/dx = dF/dy = dF/dz = 0, one obtains 3 equations with 3 unknowns (the
C x,y,z-coordinates of the extremum)
C the ensuing equations are solved by means of Gauss-Jordan elimination
C (see also: "Numerical Recipes")
C
C MP      number of equations
C NP      number of unknowns
C NDIM    dimension of space (do not alter !!!)
C CSTSVD  constant used in SVD process (see "Numerical Recipes")
C
      implicit none
c
      integer mp,np,ndim
      PARAMETER (MP = 27, NP = 9, NDIM = 3)
c
      real cstsvd
      PARAMETER (CSTSVD = 1.0E-6)
c
      integer nill
      real zero,half,one
      PARAMETER (NILL=0, ZERO=0.0, HALF=0.5, ONE=1.0)
C
      REAL F(MP), C(MP,NP), W(MP), V(MP,NP), X(MP), S(NDIM), O(NDIM)
      REAL A(NDIM,NDIM),FOFF,WMIN,WMAX,XMIN,XMAX,XEXTR,YEXTR,ZEXTR
      REAL FEXTR,OEXTR
C
      INTEGER M,N,I,J
c
cc      integer NOK,NERR1,NERR2,NUND
C
      LOGICAL INIT
C
C ... the equation matrix (fixed since differential coords are used)
C
      DATA C /9*-1,9*0,9*1,
     +        3*-1,3*0,3*1,3*-1,3*0,3*1,3*-1,3*0,3*1,
     +        -1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,
     +        -1,0,1,-1,0,1,-1,0,1,
     +        9*1,9*0,9*1,
     +        3*1,3*0,3*1,3*1,3*0,3*1,3*1,3*0,3*1,
     +        1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,
     +        1,0,1,1,0,1,1,0,1,
     +        3*1,3*0,3*-1,9*0,3*-1,3*0,3*1,
     +        1,0,-1,1,0,-1,1,0,-1,9*0,-1,0,1,-1,0,1,-1,0,1,
     +        1,0,-1,3*0,-1,0,1,1,0,-1,3*0,-1,0,1,
     +        1,0,-1,3*0,-1,0,1/
C
C ... SAVE statements for Unix
C
      SAVE C,W,V
C
Code ...
C
C ... there are 27 equations with 9 unknowns :
      M = MP
      N = NP
C ... subtract f(0,0,0) from the other function values
      FOFF = F(14)
      DO 2 I=1,M
        F(I) = F(I) - FOFF
    2 CONTINUE
C
C ... do singular value decomposition (but only ONCE !!!)
      IF (.NOT.INIT) THEN
        CALL SVDCMP (C,M,N,MP,NP,W,V)
C ... blank out W if necessary
        WMAX = W(1)
        DO 10 J=2,N
          WMAX = MAX (WMAX,W(J))
   10   CONTINUE
        WMIN = WMAX * CSTSVD
        DO 20 J=1,N
          IF (W(J).LT.WMIN) W(J) = ZERO
   20   CONTINUE
        INIT = .TRUE.
      ENDIF
C
C ... compute solution vector X by back-substitution
      CALL SVBKSB (C,W,V,M,N,MP,NP,F,X)
C ... set very small coefficients to ZERO
      XMAX = ABS(X(1))
      DO 11 J=2,N
        XMAX = MAX (XMAX,ABS(X(J)))
   11 CONTINUE
      XMIN = ABS (XMAX * CSTSVD)
      DO 21 J=1,N
        IF (ABS(X(J)).LT.XMIN) X(J) = ZERO
   21 CONTINUE
C
C ... determine position of extremum
      A (1,1) = 2.0*X(4)
      A (1,2) = X(7)
      A (1,3) = X(8)
      A (2,1) = A (1,2)
      A (2,2) = 2.0*X(5)
      A (2,3) = X(9)
      A (3,1) = A (1,3)
      A (3,2) = A (2,3)
      A (3,3) = 2.0*X(6)
C
      S (1) = - ONE * X(1)
      S (2) = - ONE * X(2)
      S (3) = - ONE * X(3)
C
C ... solve A.u = s (with u = vector(x,y,z) and s = vector(xextr,yextr,zextr))
      CALL GAUSSJ (A,NDIM,NDIM,S,1,1)
      XEXTR = S(1)
      YEXTR = S(2)
      ZEXTR = S(3)
C ... compute function value at extremum
      FEXTR = X(1)*XEXTR + X(2)*YEXTR + X(3)*ZEXTR +
     +        X(4)*XEXTR**2 + X(5)*YEXTR**2 + X(6)*ZEXTR**2 +
     +        X(7)*XEXTR*YEXTR + X(8)*XEXTR*ZEXTR + X(9)*XEXTR*ZEXTR
C ... check whether extremum is really a maximum !!!
      IF (FEXTR.LT.F(14)) THEN
        XEXTR = ZERO
        YEXTR = ZERO
        ZEXTR = ZERO
        FEXTR = F(14)
cc        call errcon ('SVDMX3 - Interpolated extremum NOT a maximum')
cc        NERR1 = NERR1 + 1
C ... check whether extremum falls within x=-1,0,1 ; y=-1,0,1 ; z=-1,0,1 !!!
      ELSE IF (ABS(XEXTR).GT.ONE.OR.ABS(YEXTR).GT.ONE.OR.
     +               ABS(ZEXTR).GT.ONE) THEN
        XEXTR = ZERO
        YEXTR = ZERO
        ZEXTR = ZERO
        FEXTR = F(14)
cc        call errcon ('SVDMX3 - Interpolated extremum OUTSIDE range')
cc        NERR2 = NERR2 + 1
      ELSE
cc        NOK = NOK + 1
      ENDIF
C
      S(1)  = O(1) + XEXTR
      S(2)  = O(2) + YEXTR
      S(3)  = O(3) + ZEXTR
      FEXTR = OEXTR + FEXTR
C
      RETURN
C
      END
