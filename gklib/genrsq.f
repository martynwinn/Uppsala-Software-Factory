c
c
c
      subroutine genrsq (maxrnd,rndseq,maxtyp,occur)
c
c ... GENRSQ = generate random sequence
c
      implicit none
c
      integer maxaat
      parameter (maxaat=30)
c
      integer maxrnd,maxtyp
c
      real occur(maxtyp),prob(maxaat),rcnt(maxaat)
c
      integer rndseq(maxrnd),counts(maxaat)
c
      real dx
c
      integer i,j,k
c
code ...
c
      write (*,*)
      call prompt (' Generating random sequence ...')
c
      dx = 0.0
      do i=1,maxtyp
        dx = dx + occur(i)
      end do
      do i=1,maxtyp
        occur(i) = occur(i)/dx
      end do
      dx = 0.0
      do i=1,maxtyp
        dx = dx + occur(i)
        prob(i) = dx
        counts(i) = 0
      end do
c
      call fvalut (' Target composition    :',maxtyp,occur)
      call prompt (' Working ...')
c
      do i=1,maxrnd
        call gkrand (dx,0.0,1.0,0)
        k = 1
        do j=2,maxtyp
          if (dx.ge.prob(j-1) .and.
     +        dx.lt.prob(j)) then
            k = j
            goto 10
          end if
        end do
   10   continue
        rndseq (i) = k
        counts (k) = counts (k) + 1
      end do
c
      do i=1,maxtyp
        rcnt(i) = counts(i)/float(maxrnd)
      end do
c
      call fvalut (' Actual composition    :',maxtyp,rcnt)
c
      return
      end
