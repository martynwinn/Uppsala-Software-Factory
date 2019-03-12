c
c ... rt_xtal.f - routines for dealing with rotations, translations
c     and crystallographic matrices etc.
c
c ===========================================================================
c
      subroutine ranrot (rotmat)
c
c ... generate random rotation matrix
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real rotmat(3,3),dx,dy,dz,det3,sum,detr
      real rx(3,3),ry(3,3),rz(3,3),dum(3,3)
c
      integer i,j
c
code ...
c
 6536 continue
c
      do i=1,3
        do j=1,3
          rotmat (i,j) = 0.0
          rx (i,j) = 0.0
          ry (i,j) = 0.0
          rz (i,j) = 0.0
          dum (i,j) = 0.0
        end do
        rotmat (i,i) = 1.0
        rx (i,i) = 1.0
        ry (i,i) = 1.0
        rz (i,i) = 1.0
        dum (i,i) = 1.0
      end do
c
      call gkrand (dx,0.0,360.0,0)
      call gkrand (dy,0.0,360.0,0)
      call gkrand (dz,0.0,360.0,0)
c
      write (*,'(a,3f8.2/)') ' Rotations around X,Y,Z (deg) :',
     +  dx,dy,dz
c
      dx = degtor*dx
      dy = degtor*dy
      dz = degtor*dz
c
      rx (2,2) = cos(dx)
      rx (3,3) = rx (2,2)
      rx (3,2) = sin(dx)
      rx (2,3) = -rx (3,2)
      dx = det3 (rx)
      write (*,'(a12,3(3f13.7,:,/,12x))')
     +   ' X Matrix : ',((rx(i,j),j=1,3),i=1,3)
      call rvalut (' Determinant of X rotation matrix :',1,dx)
c
      ry (1,1) = cos(dy)
      ry (3,3) = ry (1,1)
      ry (3,1) = sin(dy)
      ry (1,3) = -ry (3,1)
      dy = det3 (ry)
      write (*,'(a12,3(3f13.7,:,/,12x))')
     +   ' Y Matrix : ',((ry(i,j),j=1,3),i=1,3)
      call rvalut (' Determinant of Y rotation matrix :',1,dy)
c
      rz (2,2) = cos(dz)
      rz (1,1) = rz (2,2)
      rz (1,2) = sin(dz)
      rz (2,1) = -rz (1,2)
      dz = det3 (rz)
      write (*,'(a12,3(3f13.7,:,/,12x))')
     +   ' Z Matrix : ',((rz(i,j),j=1,3),i=1,3)
      call rvalut (' Determinant of Z rotation matrix :',1,dz)
c
      call mulmat (rx,ry,dum)
      call mulmat (dum,rz,rotmat)
c
      detr = det3 (rotmat)
      write (*,'(a12,3(3f13.7,:,/,12x))')
     +   '   Matrix : ',((rotmat(i,j),j=1,3),i=1,3)
      call rvalut (' Determinant of random rotation matrix :',1,detr)
c
c ... make absolutely sure this is a proper rotation matrix !
c
      if (abs(detr-1.0) .gt. 1.0e-5) goto 6536
c
      do i=1,3
        sum = rotmat(i,1)**2 + rotmat(i,2)**2 + rotmat(i,3)**2
        if (abs(sum - 1.0) .gt. 1.0e-5 ) goto 6536
        sum = rotmat(1,i)**2 + rotmat(2,i)**2 + rotmat(3,i)**2
        if (abs(sum - 1.0) .gt. 1.0e-5 ) goto 6536
      end do
c
      return
      end
