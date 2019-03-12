c
c
c
      subroutine sort2 (n,ra,rb)
c
c ... SORT2 - sorts RA and simultaneously RB in the same order
c             Numerical Recipes, p. 231
c
      implicit none
c
      integer n
      real ra(n),rb(n)
c
      real rra,rrb
c
      integer l,ir,i,j
c
code ...
c
      l = n/2 + 1
      ir = n
c
   10 continue
      if (l .gt. 1) then
        l = l - 1
        rra = ra(l)
        rrb = rb(l)
      else
        rra = ra(ir)
        rrb = rb(ir)
        ra(ir) = ra(1)
        rb(ir) = rb(1)
        ir = ir - 1
        if (ir .eq. 1) then
          ra(1) = rra
          rb(1) = rrb
          return
        end if
      end if
c
      i = l
      j = l + l
   20 continue
      if (j .le. ir) then
        if (j .lt. ir) then
          if (ra(j) .lt. ra(j+1)) j = j + 1
        end if
        if (rra .lt. ra(j)) then
          ra(i) = ra(j)
          rb(i) = rb(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if
        goto 20
      end if
c
      ra(i) = rra
      rb(i) = rrb
c
      goto 10
c
      end
