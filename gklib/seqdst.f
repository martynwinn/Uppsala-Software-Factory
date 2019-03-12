c
c
c
      subroutine seqdst (nres,nchain,jaa,cts,ich,weight)
c
c ... SEQDST = calc sequence weights from "sequence distances"
c              (Sibbald-Argos method)
c
      implicit none
c
      integer maxchn
      parameter (maxchn=1000)
c
      integer nres,nchain
c
      integer jaa(nres,nchain),cts(nres,*),ich(nchain)
      integer i,j,k,n,ii,iloop,nok,ncycle,nloop
c
      real weight(nchain),normwt(maxchn),prevwt(maxchn)
      real dx,tol,xlarge,diff
c
code ...
c
      call prompt (' Calculating sequence distances ...')
c
      if (nchain .gt. maxchn) then
        call errstp ('SEQDST - Too many molecules !!!')
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
      tol = 0.01
      ncycle = 10000
      nloop = max (10,min(100,nchain*nres))
ccc      print *,nloop,ncycle,nloop*ncycle
c
      do i=1,nchain
        weight (i) = 0.0
        normwt (i) = 1.0/float(nchain)
        prevwt (i) = normwt (i)
      end do
c
c ... count nr of types at each position and their identities
c
ccc      do i=1,nchain
ccc        call ivalut (' JAA :',nres,jaa(1,i))
ccc      end do
c
      do i=1,nres
        cts (i,1) = 0
        do j=1,nchain
          ii = jaa(i,j)
          if (cts(i,1) .gt. 0) then
            do k=2,cts(i,1)+1
              if (ii .eq. cts(i,k)) goto 10
            end do
          end if
          cts (i,1) = cts (i,1) + 1
          cts (i,1+cts(i,1)) = ii
   10     continue
        end do
ccc        write (*,6010) i,cts(i,1),(cts(i,j),j=2,1+cts(i,1))
      end do
ccc 6010 format (' Res ',i3,' # ',i2,' -> ',20i3)
c
c ... do Ncycle*Nloop iterations
c
      do iloop=1,nloop
      do i=1,ncycle
c
        do j=1,nchain
          ich (j) = 0
        end do
c
c ... construct random sequence
c
        do j=1,nres
          call gkrand (dx,0.0,float(cts(j,1)),0)
          ii = cts(j,2+int(dx))
ccc          print *,' I J ii ',i,j,ii
c
c ... accumulate sequence distances for all chains
c
          do k=1,nchain
            if (jaa(j,k) .ne. ii) ich(k) = ich(k) + 1
          end do
        end do
c
c ... find minimal distance (II)
c
        ii = nres+1
        do j=1,nchain
          if (ich(j) .lt. ii) ii = ich(j)
        end do
c
c ... count sequences with minimal distance (N)
c
        n = 0
        do j=1,nchain
          if (ich(j) .eq. ii) n = n + 1
        end do
c
c ... update weights for nearest neighbour(s) of random sequence
c
        dx = 1.0/float(n)
        do j=1,nchain
          if (ich(j) .eq. ii) weight(j) = weight(j) + dx
        end do
c
      end do
c
c ... check for convergence every NCYCLE cycles
c
      dx = 0.0
      do k=1,nchain
        dx = dx + weight(k)
      end do
      nok = 0
      xlarge = 0.0
      do k=1,nchain
        normwt(k) = weight(k)/dx
        diff = abs(normwt(k)-prevwt(k))/prevwt(k)
        if (diff .le. tol) nok = nok + 1
        xlarge = max (xlarge,diff)
        if (iloop .lt. nres) prevwt (k) = normwt (k)
      end do
      if (nok .eq. nchain) then
        call jvalut (' Weights converged :',1,ncycle*iloop)
        call fvalut (' Largest shift (%) :',1,100.0*xlarge)
        goto 9999
      end if
c
      end do
c
      call jvalut (' WARNING - Weights did not converge :',
     +  1,ncycle*nloop)
ccc      call fvalut (' Weights  :',nchain,normwt)
ccc      call fvalut (' Previous :',nchain,prevwt)
      call fvalut (' Largest shift (%) :',1,100.0*xlarge)
c
c ... normalise
c
 9999 continue
ccc      call fvalut (' Raw sums :',nchain,weight)
c
      dx = 0.0
      do k=1,nchain
        dx = dx + weight(k)
      end do
      do k=1,nchain
        weight(k) = weight(k)/dx
      end do
c
      return
      end
