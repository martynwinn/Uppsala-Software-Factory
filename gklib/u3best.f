c
c ===========================================================================
c
      subroutine u3best (x,y,n,mode,rms,u,t,ier)
c
C**** CALCULATES BEST ROTATION & TRANSLATION BETWEEN TWO VECTOR SETS
C**** SUCH THAT U*X+T IS THE BEST APPROXIMATION TO Y.
C**** THIS VERSION OF THE ALGORITHM IS OPTIMIZED FOR THREE-DIMENSIONAL
C**** REAL VECTOR SPACE.
C**** USE OF THIS ROUTINE IS RESTRICTED TO NON-PROFIT ACADEMIC
C**** APPLICATIONS.
C**** PLEASE REPORT ERRORS TO
C**** PROGRAMMER:  W.KABSCH  MAX-PLANCK-INSTITUTE FOR MEDICAL RESEARCH
C                            JAHNSTRASSE 29, 6900 HEIDELBERG, FRG.
C**** REFERENCES:  W.KABSCH  ACTA CRYST.(1978).A34,827-828
C                  W.KABSCH  ACTA CRYST.(1976).A32,922-923
C
C  X     - X(I,M) ARE COORDINATES OF ATOM # M IN SET X           (GIVEN)
C  Y     - Y(I,M) ARE COORDINATES OF ATOM # M IN SET Y           (GIVEN)
C  N     - N IS NUMBER OF ATOM PAIRS                             (GIVEN)
C  MODE  - 0:CALCULATE RMS ONLY                                  (GIVEN)
C          1:CALCULATE RMS,U,T   (TAKES LONGER)
C  RMS   - SUM OF W*(UX+T-Y)**2 OVER ALL ATOM PAIRS             (RESULT)
C  U     - U(I,J) IS   ROTATION  MATRIX FOR BEST SUPERPOSITION  (RESULT)
C  T     - T(I)   IS TRANSLATION VECTOR FOR BEST SUPERPOSITION  (RESULT)
C  IER   - 0:NO ERROR; -1:N WAS <2;                             (RESULT)
C
C-----------------------------------------------------------------------
c
      implicit none
c
      real  X(3,*),Y(3,*),U(3,3),T(3),rms
c
      integer ip(9),mode,n,ier,i,m,j,m1,m2,m3,l,k
c
      DOUBLE PRECISION R(3,3),XC(3),YC(3),WC,A(3,3),B(3,3),E0,
     1 E(3),E1,E2,E3,D,SPUR,DET,COF,H,G,CTH,STH,SQRTH,SQRT3,P,
     2 RR(6),RR1,RR2,RR3,RR4,RR5,RR6,SS(6),SS1,SS2,SS3,SS4,SS5,SS6,
     $	smax,sigma
c
      EQUIVALENCE (RR1,RR(1)),(RR2,RR(2)),(RR3,RR(3)),
     1            (RR4,RR(4)),(RR5,RR(5)),(RR6,RR(6)),
     2            (SS1,SS(1)),(SS2,SS(2)),(SS3,SS(3)),
     3            (SS4,SS(4)),(SS5,SS(5)),(SS6,SS(6)),
     4            (E1,E(1)),(E2,E(2)),(E3,E(3))
c
      DATA SQRT3/1.73205080756888D+00/
      DATA IP/1,2,4,  2,3,5,  4,5,6/
c
cc:	call PRmpti ('(''mode='',i5)', 1, mode)
c
      IER=-1
      IF (N.LT.2)RETURN
      IER=-2
      WC=0.0D+00
      DO 1 I=1,3
      XC(I)=0.0D+00
1     YC(I)=0.0D+00
      DO 2 M=1,N
      WC=WC+1.0
      DO 2 I=1,3
      XC(I)=XC(I)+X(I,M)
2     YC(I)=YC(I)+Y(I,M)
      IF (WC.LE.0.0D+00)RETURN
      DO 3 I=1,3
      XC(I)=XC(I)/WC
      YC(I)=YC(I)/WC
      DO 3 J=1,3
3     R(I,J)=0.0D+00
      E0=0.0D+00
      DO 4 M=1,N
      DO 4 I=1,3
      E0=E0+((X(I,M)-XC(I))*(X(I,M)-XC(I))+
     $	          (Y(I,M)-YC(I))*(Y(I,M)-YC(I)))
      D=(Y(I,M)-YC(I))
      DO 4 J=1,3
4     R(I,J)=R(I,J)+D*(X(J,M)-XC(J))
C**** CALCULATE DETERMINANT OF R(I,J)
      DET=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))
     1   -R(1,2)*(R(2,1)*R(3,3)-R(2,3)*R(3,1))
     2   +R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
      SIGMA=DET
C**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
      M=0
      DO 5 J=1,3
      DO 5 I=1,J
      M=M+1
5     RR(M)=R(1,I)*R(1,J)+R(2,I)*R(2,J)+R(3,I)*R(3,J)
C***************** EIGENVALUES *****************************************
C**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
      SPUR=(RR1+RR3+RR6)/3.0D+00
      COF=(RR3*RR6-RR5*RR5+RR1*RR6-RR4*RR4+RR1*RR3-RR2*RR2)/3.0D+00
	if (dabs(det) .lt. 1.d-19) then
	  if (det .ne. 0.) det = dsign(1.d-19, det)
	end if
      DET=DET*DET
C**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y-SPUR
      D=SPUR*SPUR
      H=D-COF
      G=SPUR*(COF*1.5D+00-D)-DET*0.5D+00
C**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
      IF (H.LE.D*1.0D-9)GO TO 8
      SQRTH=DSQRT(H)
      D=-G/(H*SQRTH)
      IF (D.GT. 0.9999999D+00)GO TO 6
      IF (D.LT.-0.9999999D+00)GO TO 7
C.....HANDLE CASE OF THREE DISTINCT ROOTS
      D=DACOS(D)/3.0D+00
      CTH=SQRTH*DCOS(D)
      STH=SQRTH*SQRT3*DSIN(D)
      E1=SPUR+CTH+CTH
      E2=SPUR-CTH+STH
      E3=SPUR-CTH-STH
      IF (E3.LT.0.0D+00)E3=0.0D+00
      IF (MODE.EQ.0)GO TO 35
      M1=3
      M2=1
      M3=2
      GO TO 10
C.....HANDLE SPECIAL CASE OF TWO IDENTICAL ROOTS
6     E1=SPUR+SQRTH+SQRTH
      E2=SPUR-SQRTH
      E3=E2
      IF (MODE.EQ.0)GO TO 35
      M=1
      M1=3
      M2=1
      M3=2
      GO TO 20
7     E1=SPUR+SQRTH
      E2=E1
      E3=SPUR-SQRTH-SQRTH
      IF (E3.LT.0.0D+00)E3=0.0D+00
      IF (MODE.EQ.0)GO TO 35
      M=3
      M1=1
      M2=2
      M3=3
      GO TO 20
C.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
8     E1=SPUR
      E2=SPUR
      E3=SPUR
      IF (MODE)26,35,26
C**************** EIGENVECTORS *****************************************
C.....EIGENVECTORS IN CASE OF THREE DISTINCT ROOTS
10    DO 15 L=1,2
      D=E(L)
      SS1=(D-RR3)*(D-RR6)-RR5*RR5
      SS2=(D-RR6)*RR2+RR4*RR5
      SS3=(D-RR1)*(D-RR6)-RR4*RR4
      SS4=(D-RR3)*RR4+RR2*RR5
      SS5=(D-RR1)*RR5+RR2*RR4
      SS6=(D-RR1)*(D-RR3)-RR2*RR2
	smax = ss1
c ---	Start changes
	do 1000 i=1,6
	  if (ss(i) .gt. smax) smax = ss(i)
1000	continue
	do 1100 i=1,6
1100	  ss(i) = ss(i)/smax
c ---	End changes
      J=1
      IF (DABS(SS1).GE.DABS(SS3))GO TO 12
      J=2
      IF (DABS(SS3).GE.DABS(SS6))GO TO 13
11    J=3
      GO TO 13
12    IF (DABS(SS1).LT.DABS(SS6))GO TO 11
13    D=0.0D+00
      J=3*(J-1)
      DO 14 I=1,3
      K=IP(I+J)
ccc      A(I,L)=SS(K)
      A(I,L)=SS(K)*smax
14    D=D+SS(K)*SS(K)
c ---	Start changes
ccc      D=DSQRT(D)
	D=DSQRT(D)*smax
	do 1200 i=1,6
1200	  ss(i) = ss(i)*smax
c ---	End changes
      DO 15 I=1,3
15    A(I,L)=A(I,L)/D
16    A(1,M1)=A(2,M2)*A(3,M3)-A(2,M3)*A(3,M2)
      A(2,M1)=A(3,M2)*A(1,M3)-A(3,M3)*A(1,M2)
      A(3,M1)=A(1,M2)*A(2,M3)-A(1,M3)*A(2,M2)
      GO TO 30
C.....EIGENVECTORS IN CASE OF TWO DISTINCT ROOTS
20    P=0.0D+00
cccccccccccccccc      H=SPUR+G
      H=e2
      DO 21 I=1,3
      K=(I*I+I)/2
      D=DABS(RR(K)-H)
      IF (D.LT.P)GO TO 21
      J=I
      P=D
21    CONTINUE
      P=0.0D+00
      D=0.0D+00
      L=3*(J-1)
      DO 23 I=1,3
      K=IP(I+L)
      A(I,2)=1.0D+00
      IF (I.NE.J)GO TO 22
      A(I,M)=RR(K)-H
      GO TO 23
22    A(I,M)=RR(K)
      P=P-A(I,M)
23    D=D+A(I,M)**2
      A(J,2)=P/A(J,M)
      D=DSQRT(D)
      P=DSQRT(A(1,2)**2+A(2,2)**2+A(3,2)**2)
      DO 24 I=1,3
      A(I,2)=A(I,2)/P
24    A(I,M)=A(I,M)/D
      GO TO 16
C.....EIGENVECTORS IN CASE OF THREE IDENTICAL ROOTS
26    DO 27 I=1,3
      DO 27 J=1,3
      D=0.0D+00
      IF (I.EQ.J)D=1.0D+00
27    A(I,J)=D
cc	write (6,99) a
C**** CALCULATE ROTATION MATRIX
30    DO 32 L=1,2
      D=0.0D+00
      DO 31 I=1,3
      B(I,L)=R(I,1)*A(1,L)+R(I,2)*A(2,L)+R(I,3)*A(3,L)
31    D=D+B(I,L)**2 
cc	write (6,99) b
      D=DSQRT(D)
      DO 32 I=1,3
32    B(I,L)=B(I,L)/D
      B(1,3)=B(2,1)*B(3,2)-B(2,2)*B(3,1)
      B(2,3)=B(3,1)*B(1,2)-B(3,2)*B(1,1)
      B(3,3)=B(1,1)*B(2,2)-B(1,2)*B(2,1)
      DO 33 I=1,3
      DO 33 J=1,3
33    U(I,J)=B(I,1)*A(J,1)+B(I,2)*A(J,2)+B(I,3)*A(J,3)
C**** CALCULATE TRANSLATION VECTOR
      DO 34 I=1,3
34    T(I)=YC(I)-U(I,1)*XC(1)-U(I,2)*XC(2)-U(I,3)*XC(3)
cc99	format (3d15.3)
cc	call prompt (' t done')
cc	call PRmptf('(3f10.3)', 3, t)
C**** CALCULATE RMS ERROR
35	if (e3 .le. 0.0 d 00) then
	  d = 0.
	else
	  d=dsqrt(e3)
	end if
      IF (SIGMA.LT.0.0)D=-D
	if (e2 .gt. 0.0 d 00) then
	  d = d+ dsqrt (e2)
	end if
	if (e1 .gt. 0.0 d 00) then
	  d = d+ dsqrt (e1)
	end if
ccccc      D=D+DSQRT(E2)+DSQRT(E1)
      RMS=ABS(E0-D-D)
C**** NORMAL EXIT
      IER=0
c
      RETURN
      END
