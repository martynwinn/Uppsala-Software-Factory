c
c ===========================================================================
c
      subroutine merpol (ang,amat)
c..
c..  As in the earlier cases this matrix was also pinched from
c..  the manual.
c.
c.  rotation 1 (Phi ) about Z     -----   called F
c.. rotation 2 (Psi )  about the new Y  ------ called P
c.   rotation 3 (Kappa ) about the new Z ------ called K
cc.  rotation 4 (Psi -1) about the new Y
c...  Rotation 5 (Phi-1)  about the new Z
c.. the MATRIX
c..
c   +CK+(CF*CF*SP*SP(1-CK))  +CP*SK+(CF*SF*SP*SP(1-CK) -SF*SP*SK+(CF*CP*SP(1-CK))
C.. -CP*SK+(CF*SF*SP*SP(1-CK)  +CK+(SF*SF*SP*SP(1-CK)) +CF*SP*SK+(SF*CP*SP(1-CK))
C.. +SF*SP*SK+(CF*CP*SP(1-CK)) -CF*SP*SK+(SF*CP*SP(1-CK) +CK+(CP*CP(1-CK))
c..    now make the matrices and from the angles
cc.. remember ang(1) is  Phi ang(2) is Psi and ang(3) is Kappa
cc.
c
      implicit NONE
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real ang(3),amat(3,3)
      double precision cf,cf2,sf,sf2,cp,cp2,sp,sp2,ck,sk,cck,phi,psi,kap
c
code ...
c
      phi = ang(1)*degtor
      psi = ang(2)*degtor
      kap = ang(3)*degtor
c
	cf=cos(phi)
	cf2=cf*cf
	sf=sin(phi)
	sf2=sf*sf
	cp=cos(psi)
      cp2=cp*cp
	sp=sin(psi)
	sp2=sp*sp
	ck=cos(kap) 
	sk=sin(kap)
	cck=1.0D0-ck
cc...
c..
c.. the bad code for making the matrix.  Note many things are calculated many times.
c.. a wonderful example of a bad code which works any way.
c...
	amat(1,1)=ck+(cf2*sp2*cck)
	amat(1,2)=cp*sk+(cf*sf*sp2*cck)
	amat(1,3)=-sf*sp*sk+(cf*sp*cp*cck)
	amat(2,1)=-cp*sk+(cf*sf*sp2*cck)
	amat(2,2)=ck+(sf2*sp2*cck)
	amat(2,3)=cf*sp*sk+(sf*cp*sp*cck)
	amat(3,1)=sf*sp*sk+(cf*cp*sp*cck)
	amat(3,2)=-cf*sp*sk+(sf*cp*sp*cck)
	amat(3,3)=ck+(cp2*cck)
c
	return
	end
