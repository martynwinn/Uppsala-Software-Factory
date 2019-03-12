c
c ===========================================================================
c
      subroutine mereul (ang,amat)
c
c.  the rotation matrices have been pinched from the MERLOT manual.
c.  This is used  once you solve the rotation function using
c.  merlot , you can rotate your molecule using merlot itself and 
c.  write out a merlot coordinate file - which is not a pdb file
c.  then write a code to rewrite it to the pdb format.  But this is 
c.  way of just extracting the angels is better as you can rotate
c.   other molecules too without bothering about converting them to merlot
c.. format too.
c...
c...........................................................................
c.. the matrix  is the transpose of the ccp4 euler angle matrix
c.. hence call that first and then do the transpose
c
      implicit NONE
c
      real ang(3),amat(3,3),bmat(3,3)
c
      integer i,j
c
code ...
c
      call ccpeul(ang,bmat)
      do i = 1,3
        do j=1,3
          amat(i,j)=bmat(j,i)
        enddo
      enddo
c
      return
      end
