c
c ===========================================================================
c
      subroutine anasgs (ctsym,rtsym,lprint,ierr)
c
c ... check spacegroup symmetry operators
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
      real det,theta,det3,arg
c
      integer ierr,i,j
c
      logical lprint,linvert
c
code ...
c
      if (lprint) then
        write (*,6000) ctsym
      end if
 6000 format (/' Nr of spacegroup symmetry operators : ',i3)
c
      if (ctsym .lt. 1) then
        call errcon ('No operators')
        ierr = -1
        return
      end if
c
      do i=1,ctsym
        if (lprint)
     +    write (*,'(/a7,i2,a3,3(3f13.4,10x,f13.3,:,/12x))')
     +    ' SYMOP ',i,' = ',
     +    rtsym(1,i),rtsym(4,i),rtsym(7,i),rtsym(10,i),
     +    rtsym(2,i),rtsym(5,i),rtsym(8,i),rtsym(11,i),
     +    rtsym(3,i),rtsym(6,i),rtsym(9,i),rtsym(12,i)
c
        do j=10,12
          if (rtsym(j,i) .lt. -1.0 .or.
     +        rtsym(j,i) .gt.  1.0) then
            if (.not. lprint) then
              call ivalut (' Operator nr:',1,i)
            end if
            call errcon ('Translation appears non-fractional !')
          end if
        end do
c
c ... check determinant (must be 1.000000)
c     950512 - also allow for determinant -1 !
c
        det = det3(rtsym(1,i))
        linvert = (det .lt. 0.0)
        if (lprint) write (*,'(a,f10.6)')
     +    ' Determinant of rotation matrix = ',det
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
        if (.not. linvert) then
          arg = (rtsym(1,i) + rtsym(5,i) +
     +           rtsym(9,i) - 1) * 0.5
          if (arg .lt. -1.0 .or. arg .gt. 1.0) then
            call fvalut (
     +        ' WARNING - cosine of rot angle not in [-1,+1] :',
     +        1,arg)
            arg = max (-1.000, min (arg, 1.000))
            call fvalut (' Reset to :',1,arg)
          end if
          theta = acos (arg) * rtodeg
          if (lprint) write (*,'(a,f10.6)')
     +      ' Rotation angle                 = ',theta
        end if
c
      end do
c
      ierr = 0
c
      return
      end
