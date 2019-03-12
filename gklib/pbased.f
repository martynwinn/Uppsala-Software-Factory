c
c
c
      subroutine pbased (nres,nchain,jaa,weight)
c
c ... PBASED = Henikoff-Henikoff position-based sequence weights
c              (JMB 243, pp. 574-578, 1994)
c
      implicit none
c
      integer maxchn,maxtyp
      parameter (maxchn=1000,maxtyp=25)
c
      integer nres,nchain
c
      integer jaa(nres,nchain),many(0:maxtyp)
      integer i,j,ii,ndiff
c
      real weight(nchain)
      real sumx
c
code ...
c
      call prompt (' Calculating sequence distances ...')
c
      if (nchain .gt. maxchn) then
        call errstp ('PBASED - Too many molecules !!!')
      end if
c
      if (nchain .eq. 1) then
        weight (1) = 1.0
        return
      end if
c
      if (nchain .eq. 2) then
        weight (1) = 0.5
        weight (2) = 0.5
        return
      end if
c
      do i=1,nchain
        weight (i) = 0.0
      end do
c
c ... count nr of types at each position and their identities
c
      do i=1,nres
        do j=0,maxtyp
          many (j) = 0
        end do
        ndiff = 0
        do j=1,nchain
          ii = jaa(i,j)
          if (many(ii) .eq. 0) ndiff = ndiff + 1
          many (ii) = many (ii) + 1
        end do
        do j = 1,nchain
          ii = jaa(i,j)
          weight (j) = weight (j) +
     +      (1.0 / (float(ndiff) * float(many(ii))))
        end do
      end do
c
      sumx = 0.0
      do i=1,nchain
        weight (i) = weight (i) / float(nres)
        sumx = sumx + weight (i)
      end do
      do i=1,nchain
        weight (i) = weight (i) / sumx
      end do
c
      return
      end
