c
c ===========================================================================
c
      subroutine anamat (mat,det)
c
      implicit none
c
      double precision twopi,pi,degtor,rtodeg
      parameter (twopi=6.2831853071796, pi=0.50*twopi)
      parameter (degtor=twopi/360.0,    rtodeg=360.0/twopi)
c
      real det3,mat(3,3),det,theta,dtheta,v(3),vlen,qqq
c
      integer i
c
code ...
c
      det = det3(mat)
      write (*,'(a,f10.6)')
     +  ' Determinant of rotation matrix = ',det
c
      qqq = (mat(1,1) + mat(2,2) + mat(3,3) - 1) * 0.5
      if (qqq .lt. -1.000) then
        call prompt (' WARNING - (M11+M22+M33-1)/2 < -1')
        write (*,'(a,f15.10)') ' Value = ',qqq
        call prompt (' Cannot take ACOS of this number')
        call prompt (' Check your rotation matrix !')
        if (qqq .lt. -1.0001) return
        qqq = -1.000
        call prompt (' Reset to -1.0000')
      else if (qqq .gt. 1.000) then
        call prompt (' WARNING - (M11+M22+M33-1)/2 > +1')
        write (*,'(a,f15.10)') ' Value = ',qqq
        call prompt (' Cannot take ACOS of this number')
        call prompt (' Check your rotation matrix !')
        if (qqq .gt. 1.0001) return
        qqq = 1.000
        call prompt (' Reset to +1.0000')
      end if
c
      qqq = max (-1.000, min(1.000, qqq))
      theta = acos (qqq)
      dtheta = theta*rtodeg
      write (*,'(a,f10.6)')
     +  ' Rotation angle                 = ',dtheta
c
      if (1 .gt. 0) return
c
c ... the following gives sign differences with XPLOR
c
      if (dtheta .eq. 0.0) then
        write (*,'(a)')
     +    ' ERROR - indeterminate direction cosines'
      else if (dtheta .eq. 180.0) then
        v(1) = mat(1,1) + 1.0
        v(2) = mat(2,1)
        v(3) = mat(3,1)
        vlen = sqrt (v(1)**2 + v(2)**2 + v(3)**2)
        write (*,'(a,3f10.4)')
     +    ' Direction cosines of rotation axis : ',
     +    (v(i)/vlen,i=1,3)
      else
        vlen = -0.5 / sin(theta)
        v(1) = mat(2,3)-mat(3,2)
        v(2) = mat(1,3)-mat(3,1)
        v(3) = mat(2,1)-mat(1,2)
        write (*,'(a,3f10.4)')
     +    ' Direction cosines of rotation axis : ',
     +    (v(i)*vlen,i=1,3)
      end if
c
      return
      end
