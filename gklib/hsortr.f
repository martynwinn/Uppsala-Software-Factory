c
c
c
      subroutine hsortr (n,ra)
c
c ... HSORTR - heapsort of real array
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
   10 continue
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
   20 continue
      if (j.le.ir) then
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
        goto 20
      end if
      ra(i)=rra
      goto 10
c
      end
