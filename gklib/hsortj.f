c
c
c
      subroutine hsortj (n,ra)
c
c ... HSORTR - heapsort of real array - without GOTOs !
c              (~25% slower than the version with GOTOs)
c
      implicit none
c
      integer n
      real ra(n)
c
      real rra
c
      integer l,ir,i,j
c
code ...
c
      l=n/2+1
      ir=n
c
      do while (.true.)
        if (l.gt.1) then
          l=l-1
          rra=ra(l)
        else
          rra = ra(ir)
          ra(ir) = ra(1)
          ir=ir-1
          if (ir.eq.1) then
            ra(1)=rra
            return
          end if
        end if
        i = l
        j = l + l
c
        do while (j .le. ir)
          if (j.lt.ir) then
            if (ra(j).lt.ra(j+1)) j=j+1
          end if
          if (rra.lt.ra(j)) then
            ra(i)=ra(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if
        end do
        ra(i)=rra
      end do
c
      end
