c
c
c
      subroutine frcsym (x, min, max, rt, ct, errcod)
c
c ---	FoRCe_SYM
c ---	Given the fractional coords x, force via symops
c	into the envelope. Errcod .ne. 0 if impossible
c ---	Alwyn Jones
c
      implicit none
c
      integer ct, errcod
      real x(3), min(3), max(3), rt(12,*)
c
      integer i, j
      real x1(3)
c
code ...
c
      do 100 i=1,3
        if (x(i) .lt. min(i)) goto 110
        if (x(i) .gt. max(i)) goto 110
100   continue
      errcod = 0
      return
c
110   continue
      do 200 i=1,ct
        errcod = 0
        call vecrtv (x, x1, 1, rt(1,i), rt(10,i))
        do 210 j=1,3
220       if (x1(j) .lt. min(j)) then
            x1(j) = x1(j)+ 1.
            goto 220
          end if
c
230       if (x1(j) .gt. max(j)) then
            x1(j) = x1(j)- 1.
            goto 230
          end if
c
          if (x1(j) .lt. min(j)) then
            errcod = 1
          end if
210     continue
        if (errcod .eq. 0) then
          do 240 j=1,3
240         x(j) = x1(j)
            goto 300
        end if
200   continue
c
c ---	Unable to do it but force positive value and lt 1.
c
      do 400 i=1,3
410     if (x(i) .lt. 0.) then
          x(i) = x(i)+1.
          goto 410
        end if
420     if (x(i) .ge. 1.) then
          x(i) = x(i)-1.
          goto 420
        end if
400   continue
c
300   continue
c
      return
      end
