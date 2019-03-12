c
c ===========================================================================
c
      subroutine anancs (ctsym,rtsym,lprint,ierr)
c
c ... check NCS operators
c
      implicit none
c
      double precision twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      integer ctsym
      real rtsym(12,ctsym)
c
      real det,theta,det3,arg,ip1,ip2,ip3,iptot
c
      integer ierr,i
c
      logical lprint,linvert
c
code ...
c
      if (lprint) then
        write (*,6000) ctsym
      end if
 6000 format (/' Nr of RT operators : ',i3)
c
      if (ctsym .lt. 1) then
        call errcon ('No operators')
        ierr = -1
        return
      end if
c
      do i=1,ctsym
        if (lprint)
     +    write (*,'(/a7,i2,a3,3(3f13.7,10x,f13.3,:,/12x))')
     +    ' RT-OP ',i,' = ',
     +    rtsym(1,i),rtsym(4,i),rtsym(7,i),rtsym(10,i),
     +    rtsym(2,i),rtsym(5,i),rtsym(8,i),rtsym(11,i),
     +    rtsym(3,i),rtsym(6,i),rtsym(9,i),rtsym(12,i)
c
c ... check determinant (must be 1.000000)
c     950512 - also allow for determinant -1 !
c
        det = det3(rtsym(1,i))
        linvert = (det .lt. 0.0)
        if (lprint) write (*,'(a,f17.6)')
     +    ' Determinant of rotation matrix',det
        if (abs(abs(det)-1.000) .gt. 0.01) then
          if (.not. lprint) then
            call ivalut (' Operator nr:',1,i)
          end if
          call errcon (
     +      'Determinant differs from +-ONE by more than 0.01')
          ierr = -2
          return
        else if (abs(abs(det)-1.000) .gt. 0.00001) then
          if (.not. lprint) then
            call ivalut (' Operator nr:',1,i)
          end if
          call errcon (
     +      'Determinant differs from +-ONE by more than 0.00001')
        end if
c
        if (linvert) then
          call prompt (' WARNING - Determinant < 0 -> inversion !')
        end if
c
c ... check if column vectors are orthogonal
c
        ip1 = rtsym(1,i)*rtsym(4,i) + rtsym(2,i)*rtsym(5,i) +
     +        rtsym(3,i)*rtsym(6,i)
        ip2 = rtsym(1,i)*rtsym(7,i) + rtsym(2,i)*rtsym(8,i) +
     +        rtsym(3,i)*rtsym(9,i)
        ip3 = rtsym(4,i)*rtsym(7,i) + rtsym(5,i)*rtsym(8,i) +
     +        rtsym(6,i)*rtsym(9,i)
        iptot = abs(ip1) + abs(ip2) + abs(ip3)
        if (lprint) write (*,'(a,2x,3f12.6)')
     +    ' Column-vector products (12,13,23)',ip1,ip2,ip3
        if (iptot .gt. 0.001) then
          call errcon ('Seriously non-orthogonal column vectors')
          ierr = -3
          return
        else if (iptot .gt. 0.00001) then
          call errcon ('WARNING - Non-orthogonal column vectors')
        end if
c
        if (lprint) then
          call matejd (rtsym(1,i))
        end if
c
        arg = (rtsym(1,i) + rtsym(5,i) +
     +         rtsym(9,i) - 1) * 0.5
        if (arg .lt. -1.0 .or. arg .gt. 1.0) then
          call fvalut (
     +      ' WARNING - cosine of rot angle not in [-1,+1] :',
     +      1,arg)
          arg = max (-1.000, min (arg, 1.000))
          call fvalut (' Reset to :',1,arg)
        end if
        theta = acos (arg) * rtodeg
        if (lprint) write (*,'(a,f13.3)')
     +    ' Rotation angle                       ',theta
c
      end do
c
      ierr = 0
c
      return
      end
