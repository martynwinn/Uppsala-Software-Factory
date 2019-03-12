c
c ===========================================================================
c
      subroutine u3qion (xyzfix,xyzmov,natoms,mode,rms,rot,tra,ierr) 
c
c ... U3QION - get least - squares RT operator using quaternion method
c
c ... See: SK Kearsley, Acta Cryst A45, 208 - 210 (1989) 
c
      implicit none
c
      integer natoms,ierr,mode
      real xyzfix(3,natoms),xyzmov(3,natoms),rot(3,3),tra(3),rms
c
      real*8 lsqmat(4,4),xyzsub(3),xyzadd(3),cogmov(3),cogfix(3) 
      real sngmat(4,4),eigvec(4,4) 
      real eigval(4),q1,q2,q3,q4
c
      integer i,j,nrot
c
code ...
c
      ierr = 0
c
c ... mode = 0 not implemented
c
      if (mode .ne. 1)  then
        ierr = - 1
        call errcon ('U3QION - MODE = 0 not implemented') 
        return
      end if
c
c ... initialise
c
      do i = 1,3
        cogfix (i) = 0.0
        cogmov (i) = 0.0
      end do
c
      do i = 1,4
        do j = 1,4
          lsqmat (i,j) = 0.0D0
        end do
      end do
c
c ... sum coordinates
c
      do i = 1,natoms
        do j = 1,3
          cogfix (j) = cogfix (j) + dble (xyzfix (j,i)) 
          cogmov (j) = cogmov (j) + dble (xyzmov (j,i)) 
        end do
      end do
c
c ... get centres-of-gravity and translation vector
c
      do i = 1,3
        cogmov (i) = cogmov (i) / dble (natoms) 
        cogfix (i) = cogfix (i) / dble (natoms) 
      end do
c
ccc      print *,' C - of - G moving :',cogmov
ccc      print *,' C - of - G fixed  :',cogfix
c
c ... generate matrix of coordinate differential sums
c
      do i=1,natoms
c
c ... coordinate differentials
c
        do j = 1,3
          xyzsub(j) = dble (xyzmov(j,i)) - cogmov (j) - 
     +                (dble (xyzfix(j,i)) - cogfix (j)) 
          xyzadd(j) = dble (xyzmov(j,i)) - cogmov (j) + 
     +                (dble (xyzfix(j,i)) - cogfix (j)) 
        end do
c
c ... diagonal matrix elements
c
        lsqmat (1,1) = lsqmat (1,1) + 
     +     xyzsub(1)**2 + xyzsub(2)**2 + xyzsub(3)**2
        lsqmat (2,2) = lsqmat (2,2) + 
     +     xyzadd(2)**2 + xyzadd(3)**2 + xyzsub(1)**2
        lsqmat (3,3) = lsqmat (3,3) + 
     +     xyzadd(1)**2 + xyzadd(3)**2 + xyzsub(2)**2
        lsqmat (4,4) = lsqmat (4,4) + 
     +     xyzadd(1)**2 + xyzadd(2)**2 + xyzsub(3)**2
c
c ... off-diagonal matrix elements (only half; matrix is symmetric) 
c
        lsqmat (1,2) = lsqmat (1,2) + 
     +     xyzadd(2)*xyzsub(3) - xyzsub(2)*xyzadd(3) 
        lsqmat (1,3) = lsqmat (1,3) + 
     +     xyzsub(1)*xyzadd(3) - xyzadd(1)*xyzsub(3) 
        lsqmat (1,4) = lsqmat (1,4) + 
     +     xyzadd(1)*xyzsub(2) - xyzsub(1)*xyzadd(2) 
        lsqmat (2,3) = lsqmat (2,3) + 
     +     xyzsub(1)*xyzsub(2) - xyzadd(1)*xyzadd(2) 
        lsqmat (2,4) = lsqmat (2,4) + 
     +     xyzsub(1)*xyzsub(3) - xyzadd(1)*xyzadd(3) 
        lsqmat (3,4) = lsqmat (3,4) + 
     +     xyzsub(2)*xyzsub(3) - xyzadd(2)*xyzadd(3) 
c
      end do
c
c ... other half of the off - diagonal matrix elements
c
      lsqmat (2,1) = lsqmat (1,2) 
      lsqmat (3,1) = lsqmat (1,3) 
      lsqmat (4,1) = lsqmat (1,4) 
      lsqmat (3,2) = lsqmat (2,3) 
      lsqmat (4,2) = lsqmat (2,4) 
      lsqmat (4,3) = lsqmat (3,4) 
c
      do i = 1,4
        do j = 1,4
          sngmat (i,j) = lsqmat (i,j) 
        end do
      end do
c
c ... get eigenvalues and eigenvectors
c
      call jacobi (sngmat,4,4,eigval,eigvec,nrot,ierr) 
c
      if (ierr .ne. 0)  then
        call errcon ('In JACOBI - cannot superimpose coordinates') 
        return
      end if
c
c ... sort them
c
ccc      call jvalut (' NROT Jacobi   :',1,nrot) 
ccc      call rvalut (' EIGVAL before :',4,eigval) 
c
      call eigsrt (eigval,eigvec,4,4) 
c
ccc      call rvalut (' EIGVAL after  :',4,eigval) 
ccc      call rvalut (' EIGVEC # 4    :',4,eigvec (1,4) ) 
c
c ... sum of squared errors = eigenvalue nr 4
c
      rms = eigval (4) 
c
c ... the 4-th (and smallest) eigenvalue is the one we want; now generate
c     the rotation matrix from it
c
      q1 = eigvec (1,4) 
      q2 = eigvec (2,4) 
      q3 = eigvec (3,4) 
      q4 = eigvec (4,4) 
c
      rot (1,1) = q1*q1 + q2*q2 - q3*q3 - q4*q4
      rot (2,1) = 2.0* (q2*q3 - q1*q4) 
      rot (3,1) = 2.0* (q2*q4 + q1*q3) 
      rot (1,2) = 2.0* (q2*q3 + q1*q4) 
      rot (2,2) = q1*q1 + q3*q3 - q2*q2 - q4*q4
      rot (3,2) = 2.0* (q3*q4 - q1*q2) 
      rot (1,3) = 2.0* (q2*q4 - q1*q3) 
      rot (2,3) = 2.0* (q3*q4 + q1*q2) 
      rot (3,3) = q1*q1 + q4*q4 - q2*q2 - q3*q3
c
c ... calculate translation vector
c
      do i = 1,3
        tra (i) = sngl(cogmov(i)) - rot(i,1)*sngl(cogfix(1)) - 
     +     rot(i,2)*sngl(cogfix(2)) - rot(i,3)*sngl(cogfix(3)) 
      end do
c
      return
      end
