c
c ================================================================
c
      subroutine xps_long_legend (text)
c
c ... split a long legend up into smaller ones (text must be a variable !)
c
      integer i,j,k,ll,length
c
      character text*(*)
c
code ...
c
      ll = length(text)
      if (ll .le. 0) return
c
      call pretty (text)
      i = 1
c
   30 continue
      if (i .gt. ll) goto 40
      if (i+70 .ge. ll) then
        k = min (i + 70, ll)
        call xps_legend (text(i:k))
        goto 40
      end if
      k = i + 70
      do j=k,i,-1
        if (text(j:j) .eq. ' ') then
          k = j
          goto 20
        end if
      end do
   20 continue
      call xps_legend (text(i:k))
      i = k + 1
      goto 30
c
   40 continue
c
      return
      end
