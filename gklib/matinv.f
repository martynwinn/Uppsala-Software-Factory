c
c ========================================================================
c
	subroutine matinv(a,n,d,l,m)
c
      implicit none
c
      real A(*),L(*),M(*),d,biga,hold
c
      integer n,nk,k,kk,j,iz,i,ij,ki,ji,jp,jk,ik
      integer kj,jq,jr
C
C        PURPOSE
C           INVERT A MATRIX
C
C        USAGE
C           CALL MATINV(A,N,D,L,M)
C
C        DESCRIPTION OF PARAMETERS
C           A - INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY
C               RESULTANT INVERSE.
C           N - ORDER OF MATRIX A
C           D - RESULTANT DETERMINANT
C           L - WORK VECTOR OF LENGTH N
C           M - WORK VECTOR OF LENGTH N
C
C        REMARKS
C           MATRIX A MUST BE A GENERAL MATRIX
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT
C           IS ALSO CALCULATED. A DETERMINANT OF ZERO INDICATES THAT
C           THE MATRIX IS SINGULAR.
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A,D,BIGA,HOLD,DABS
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  ABS IN STATEMENT
C        10 MUST BE CHANGED TO DABS.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C        SEARCH FOR LARGEST ELEMENT
C
      D=1.0
      NK=-N
      DO 190 K=1,N
         NK=NK+N
         L(K)=K
         M(K)=K
         KK=NK+K
         BIGA=A(KK)
         DO 30 J=K,N
            IZ=N*(J-1)
            DO 30 I=K,N
               IJ=IZ+I
   10          IF( ABS(BIGA)- ABS(A(IJ))) 20,30,30
   20          BIGA=A(IJ)
               L(K)=I
               M(K)=J
   30    CONTINUE
C
C        INTERCHANGE ROWS
C
         J=L(K)
         IF(J-K) 60,60,40
   40    KI=K-N
         DO 50 I=1,N
            KI=KI+N
            HOLD=-A(KI)
            JI=KI-K+J
            A(KI)=A(JI)
   50    A(JI) =HOLD
C
C        INTERCHANGE COLUMNS
C
   60    I=M(K)
         IF(I-K) 90,90,70
   70    JP=N*(I-1)
         DO 80 J=1,N
            JK=NK+J
            JI=JP+J
            HOLD=-A(JK)
            A(JK)=A(JI)
   80    A(JI) =HOLD
C
C        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
C        CONTAINED IN BIGA)
C
   90    IF(BIGA) 110,100,110
  100    D=0.0
         RETURN
  110    DO 130 I=1,N
            IF(I-K) 120,130,120
  120       IK=NK+I
            A(IK)=A(IK)/(-BIGA)
  130    CONTINUE
C
C        REDUCE MATRIX
C
         DO 160 I=1,N
            IK=NK+I
            HOLD=A(IK)
            IJ=I-N
            DO 160 J=1,N
               IJ=IJ+N
               IF(I-K) 140,160,140
  140          IF(J-K) 150,160,150
  150          KJ=IJ-I+K
               A(IJ)=HOLD*A(KJ)+A(IJ)
  160    CONTINUE
C
C        DIVIDE ROW BY PIVOT
C
         KJ=K-N
         DO 180 J=1,N
            KJ=KJ+N
            IF(J-K) 170,180,170
  170       A(KJ)=A(KJ)/BIGA
  180    CONTINUE
C
C        PRODUCT OF PIVOTS
C
         D=D*BIGA
C
C        REPLACE PIVOT BY RECIPROCAL
C
         A(KK)=1.0/BIGA
  190 CONTINUE
C
C        FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  200 K=(K-1)
      IF(K) 270,270,210
  210 I=L(K)
      IF(I-K) 240,240,220
  220 JQ=N*(K-1)
      JR=N*(I-1)
      DO 230 J=1,N
         JK=JQ+J
         HOLD=A(JK)
         JI=JR+J
         A(JK)=-A(JI)
  230 A(JI) =HOLD
  240 J=M(K)
      IF(J-K) 200,200,250
  250 KI=K-N
      DO 260 I=1,N
         KI=KI+N
         HOLD=A(KI)
         JI=KI-K+J
         A(KI)=-A(JI)
  260 A(JI) =HOLD
      GO TO 200
  270 RETURN
      END
