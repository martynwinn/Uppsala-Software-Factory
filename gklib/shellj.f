c
c ==========================================================================
c
      subroutine shellj (a,inx,n)
c
c ... shell sort algorithm - without GOTOs !
c     about the same speed as the version with GOTOs
c
c ... a contains n real values by which the index array inx
c     should be sorted; array a itself is also sorted, so if you don't
c     this to happen, use a scratch array !!!
c
      implicit none
c
      integer n,inx(n),i,m,k,j,ipm,iw
c
      real a(n),w
c
code ...
c
      i=1
c
      do while (i .le. n)
        m=2*i-1
        i=2*i
      end do
c
      do while (.true.)
        m=m/2
        if (m.eq.0) return
c
        k=n-m
        do j=1,k
          i=j
          do while (i .ge. 1)
            ipm=i+m
            if (a(ipm).ge.a(i)) then
              i = -1
            else
c
c ... swap as and inxs
c
              w=a(i)
              a(i)=a(ipm)
              a(ipm)=w
c
              iw=inx(i)
              inx(i)=inx(ipm)
              inx(ipm)=iw
c
              i=i-m
            end if
          end do
        end do
      end do
c
      end
