c
c ========================================================================
c
      subroutine matejd (trmat)
c
c ... Eleanor Dodson's routine
c
c     this s/r works out various rotation angles corresponding
c     to a given rotation matrix.  it could easily be extended
c     to give other sets.
c
      implicit none
c
      double precision twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      real trmat(3,3)
      double precision dc(5)
      double precision conv,phi,omega,chi,theta1,theta2,theta3
      double precision trace,det,r13,r23,r31,r32,coschi,sinchi
      double precision chk,somega,beta,r33,r11,alpha,rchk,r22,r21
      double precision r12,gamma,qqq,xplphi,xplpsi,xplkap,thetap,thetam
c
      integer i
c
      character val1*12,val2*12,val3*12
c
code ...
c
c     check matrix determinant eq 1
c
      det = trmat(1,1 )*
     1(trmat(2,2 )*trmat(3,3 )-trmat(3,2 )*trmat(2,3 ))
     1 - trmat(1,2 )*
     1(trmat(2,1 )*trmat(3,3 )-trmat(3,1 )*trmat(2,3 ))
     1 + trmat(1,3 )*
     1(trmat(2,1 )*trmat(3,2)-trmat(3,1 )*trmat(2,2 ))
c
c      write(6,9876) det,((trmat(k,j),j=1,3),k=1,3)
c9876  format(' Matrix & determinant ..',4x,f10.6,3(/,3x,3f10.5))
c
      if(abs(det-1.0).gt.0.001) then
        call errcon ('MATEJD - Determinant differs from 1.000 !')
        write (*,'(1x,a,f15.6)') 'Value :',det
        return
      end if
c
      conv=rtodeg
c
      phi=-999.99
      omega=0.
      chi=0.
      dc(1)=0.0
      dc(2)=0.0
      dc(3)=0.0
      theta1=0.
      theta2=0.
      theta3=0.
      xplphi=0.
      xplpsi=0.
      xplkap=0.
c
      trace=trmat(1,1)+trmat(2,2)+trmat(3,3)
      if(trace.ge.3.0) go to 100
      if(trmat(3,3).ge.0.99999) go to 13
      if(trmat(3,3).le.-0.999999) go to 13
      theta2=conv*acos(trmat(3,3))
c
      r13=trmat(1,3)
      r23=trmat(2,3)
      r31=trmat(3,1)
      r32=trmat(3,2)
      if(r23.ge.0.) go to 9
      r13=-r13
      r23=-r23
      r31=-r31
      r32=-r32
      theta2=-theta2
9     theta3=conv*atan2(r32,-r31)
      theta1=conv*atan2(r23,r13)
      go to 18
c
c     beta = 0    can only find alpha +  gamma.
c     beta = 180  can only find alpha - gamma.
c
13    theta2=0.
      if(trmat(3,3).lt.0.) theta2=180.
      theta3=conv*atan2(trmat(2,1),trmat(2,2))
      theta1=0.
      if(theta2.eq.0)then
        chi=theta3
        omega=0.
        phi=-999.99
      endif
c
cxyz
c prob wrong...      if(theta2.eq.180.) chi=180.
c
      go to 10
c
c     spherical polars ref patterson int tab vol2 p59
c
18    chi=0.0
      dc(1)=0.0
      dc(2)=0.0
      dc(3)=0.0
      phi=-999.99
      omega=-999.99
      coschi=(trace/2.0D0) -0.5D0
      sinchi=sqrt(abs(1.-coschi*coschi))
      chi=conv*atan2(sinchi,coschi)
      chk=abs(chi)
      if(chk.gt.179.0) go to 10
c
      if(chk.lt.1.0) then
c        write(6,*)'  rotation less than 1 degree * uninteresting !!'
        go to 100
      endif
c
      dc(1)=(trmat(3,2)-trmat(2,3))/(2.*sinchi)
      dc(2)=(trmat(1,3)-trmat(3,1))/(2.*sinchi)
      dc(3)=(trmat(2,1)-trmat(1,2))/(2.*sinchi)
c      print *,' 1 ',dc
      chk=dc(1)**2 +dc(2)**2 +dc(3)**2
      if(chk.gt.0.8) go to 11
c
c     if chi =180  coschi=-1,sinchi=0
c
10    dc(1)=sqrt(abs(0.5D0*trmat(1,1)+0.5D0))
      if(dc(1).lt.0.001)dc(1)=0.0
      dc(2)=sqrt(abs(0.5D0*trmat(2,2)+0.5D0))
      if(dc(2).lt.0.001)dc(2)=0.0
      dc(3)=sqrt(abs(0.5D0*trmat(3,3)+0.5D0))
      if(dc(3).lt.0.001)dc(3)=0.0
c      print *,' 2 ',dc
c
c     assume dc(1) positive
c
      if(dc(1).eq.0.0)goto 17
      if(dc(2).ne.0.0)dc(2)=dc(2)*sign(1.0,trmat(1,2))
      if(dc(3).ne.0.0)dc(3)=dc(3)*sign(1.0,trmat(1,3))
c      print *,' 3 ',dc
      goto 11
c
17    if(dc(2).eq.0.0)goto 11
      if(dc(3).eq.0.0)goto 11
      dc(3)=dc(3)*sign(1.0,trmat(2,3))
11    dc(4)=dc(1)
      dc(5)=dc(2)
c      print *,' 4 ',dc
c
12    somega=sqrt(abs(1.0D0-dc(3)*dc(3)))
      omega=conv*atan2(somega,dc(3))
      phi=-999.99
      if(omega.gt.1.0)  phi=conv*atan2(dc(2),dc(1))
c
c  end of copy
c
 100  continue
c      print *,' 5 ',dc
c
      val1 = ' *undefined*'
      if (abs(theta1) .lt. 999.9) write (val1,'(f12.3)') theta1
      val2 = ' *undefined*'
      if (abs(theta2) .lt. 999.9) write (val2,'(f12.3)') theta2
      val3 = ' *undefined*'
      if (abs(theta3) .lt. 999.9) write (val3,'(f12.3)') theta3
      write (6,102) ' Crowther Alpha Beta Gamma          ',
     +  val1,val2,val3
c
      val1 = ' *undefined*'
      if (abs(omega) .lt. 999.9) write (val1,'(f12.3)') omega
      val2 = ' *undefined*'
      if (abs(phi) .lt. 999.9) write (val2,'(f12.3)') phi
      val3 = ' *undefined*'
      if (abs(chi) .lt. 999.9) write (val3,'(f12.3)') chi
      write (6,102) ' Spherical polars Omega Phi Chi     ',
     +  val1,val2,val3
c
 102  format (a,3a12)
 103  format (a,3f12.3)
      write (6,101) (dc(i),i=1,3)
 101  format (' Direction cosines of rotation axis ',3f12.6)
c
c     dave smiths  now
c
      qqq = max (-1.0, min (1.0, trmat(3,1)))
      beta=conv*acos( qqq )
c
      r33=trmat(3,3)
      r23=trmat(2,3)
      r12=trmat(1,2)
      r11=trmat(1,1)
      if(r23.ge.0.) go to 3
      r23=-r23
      r33=-r33
      r12=-r12
      r11=-r11
      beta=-beta
3     continue
      alpha=999.99
      rchk=r23*r23+ r33*r33
      if(rchk.gt.0.)alpha=conv*atan2(-r23,r33)
      rchk=r12*r12+ r11*r11
      gamma=999.99
      if(rchk.gt.0.)gamma=conv*atan2(-r12,r11)
c
      val1 = ' *undefined*'
      if (abs(alpha) .lt. 999.9) write (val1,'(f12.3)') alpha
      val2 = ' *undefined*'
      if (abs(beta) .lt. 999.9) write (val2,'(f12.3)') beta
      val3 = ' *undefined*'
      if (abs(gamma) .lt. 999.9) write (val3,'(f12.3)') gamma
ccc      write (6,102) ' Dave Smith                         ',
ccc     +  val1,val2,val3
c
c ... xplor (rossmann 1962) polar
c
      r11 = trmat(1,1)
      r12 = trmat(1,2)
      r13 = trmat(1,3)
      r21 = trmat(2,1)
      r22 = trmat(2,2)
      r23 = trmat(2,3)
      r31 = trmat(3,1)
      r32 = trmat(3,2)
      r33 = trmat(3,3)
c
      qqq = max (-1.0D0, min (1.0D0, 0.5D0*(r11+r22+r33-1.0D0)))
      xplkap = acos (qqq) * rtodeg
      write (val3,'(f12.3)') xplkap
      val1 = ' *undefined*'
      if (abs(r23-r32).gt.1.0e-10) then
        xplphi = atan2(r21-r12,r23-r32) * rtodeg
        write (val1,'(f12.3)') xplphi
      end if
      val2 = ' *undefined*'
      qqq = (r13-r31)*cos(xplphi*degtor)
      if (abs(qqq).gt.1.0e-10) then
        xplpsi = atan2(r32-r23,qqq) * rtodeg
        write (val2,'(f12.3)') xplpsi
      end if
c
      write (6,102) ' X-PLOR polars Phi Psi Kappa        ',
     +  val1,val2,val3
c
c ... lattman theta+,2,-
c
      val1 = ' *undefined*'
      val3 = ' *undefined*'
      if (abs(r32).gt.1.0e-10 .and. abs(r23).gt.1.0e-10) then
        thetap = rtodeg * (atan2(r31,-r32) + atan2(r13,r23))
        write (val1,'(f12.3)') thetap
        thetam = rtodeg * (atan2(r31,-r32) - atan2(r13,r23))
        write (val3,'(f12.3)') thetam
      end if
      qqq = max (-1.0D0, min (1.0D0, r33))
      theta2 = rtodeg * acos (qqq)
      write (val2,'(f12.3)') theta2
c
      write (6,102) ' Lattmann Theta+ Theta2 Theta-      ',
     +  val1,val2,val3
c
c      write(6,102)alpha,beta,gamma
c102   format(' Dave Smith                         ',3f12.3)
c
      return
      end
