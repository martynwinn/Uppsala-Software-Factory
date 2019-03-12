c
c
c
      subroutine crank (n,w,s)
c
c ... CRANK - used by NONPAR to return ranks of array elements
c
      implicit none
c
      integer n
      real w(n)
c
      real s,rank,t
c
      integer j,jt,ji
c
code ...
c
      s = 0.0
      j = 1
    1 continue
      if (j .lt. n) then
        if (w(j+1) .ne. w(j)) then
          w(j) = j
          j = j + 1
        else
          do jt=j+1,n
            if (w(jt).ne.w(j)) goto 2
          end do
          jt = n+1
    2     continue
          rank = 0.5*(j+jt-1)
          do ji=j,jt-1
            w(ji) = rank
          end do
          t = jt-j
          s = s + t**3 - t
          j = jt
        end if
        goto 1
      end if
c
      if (j .eq. n) w(n) = n
c
      return
      end
