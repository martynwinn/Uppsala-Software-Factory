c
c ===========================================================================
c
	subroutine lsqgjk (x1, x2, n, rms, rtx, ier)
c
c ... modified from alwyn's lsqslw
c
c ... 2004-09-28 - this now calls a routine that uses quaternions
c                  after i found a case in which U3BEST fails miserably
c
c ---	Find the least squares best bit of the set of points X1,X2
c ---	There are N points, with rotation U and translation T
c ---	The transformations are made on the set X2 to bring them on X1
c ---	This is a fake to call Kabsch's U3BEST until I write my own.
c
	implicit none
c
	real x1(3,*),x2(3,*),u(3,3),t(3),rt(12),rtx(12),rms
c
	integer i,n,mode,ier
c
        equivalence ( u(1,1), rt(1)  )
        equivalence ( t(1),   rt(10) )
c
        data u / 1.0,0.,0., 0.,1.,0., 0.,0.,1. /
        data t / 0.,0.,0. /
c
code ...
c
      mode = 1
      rms = 9999999.999
      ier = 0
c
	if (n .le. 2) then
	  call errcon 
     +	    ('LSQGJK - Not enough atoms for least squares fit')
        ier = -1
        goto 9999
	end if
c
c	call u3best (x2, x1, n, mode, rms, u, t, ier)
c	if(ier .ne. 0)then
c	  call errcon ('In LSQGJK/U3BEST')
c         goto 9999
c	end if
c
c ... call quaternion routine (note: parameter MODE is a dummy !!!)
c
	call u3qion (x2, x1, n, mode, rms, u, t, ier)
	if (ier .ne. 0) then
	  call errcon ('In LSQGJK/U3QION')
        ier = -2
        goto 9999
	end if
c
      if (rms .ge. 0.0) then
        rms = sqrt(rms/float(n))
      else
        call rvalut (
     +    ' WARNING - Reset negative RMSD in LSQGJK :',1,rms)
        rms = 0.0
      end if
c
 9999 continue
      do i=1,12
        rtx(i) = rt(i)
      end do
c
	return
	end
