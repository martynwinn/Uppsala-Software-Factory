c
c
c
      subroutine frctrn (i1, i2, orgn, grid, min, max, 
     $                   rt, ct, errcod)
c
c ---	FoRCe_TRaNslation into the asymmetric unit
c ---	Given the grid index, force via symops
c	into the envelope. Errcod .ne. 0 if impossible
c ---	Alwyn Jones
c
      implicit none
c
      integer ct, errcod, i1(3), i2(3), grid(3), orgn(3)
      real min(3), max(3), rt(12,*)
c
      integer i, j
      real x(3), x1(3)
c
code ...
c
      do 50 i=1,3
50      x(i) = float(i1(i)+ orgn(i)- 1)/float(grid(i))
      do 100 i=1,3
        if (x(i) .lt. min(i)) goto 110
        if (x(i) .gt. max(i)) goto 110
100   continue
      do 60 i=1,3
60      i2(i) = i1(i)
      errcod = 0
      return
c
110   continue
      do 200 i=1,ct
        errcod = 0
        call vecrtv (x, x1, 1, rt(1,i), rt(10,i))
        do 210 j=1,3
220       if (x1(j) .lt. min(j)) then
            if (abs(x1(j)-min(j)) .gt. 0.0001) then
              x1(j) = x1(j)+ 1.
              goto 220
            end if
          end if
c
230       if (x1(j) .gt. max(j)) then
            if (abs(x1(j)-max(j)) .gt. 0.0001) then
              x1(j) = x1(j)- 1.
              goto 230
            end if
          end if
c
          if (x1(j) .lt. min(j)) then
            if (abs(x1(j)-min(j)) .gt. 0.0001) then
              errcod = 1
            end if
          end if
210     continue
c
        if (errcod .eq. 0) then
          do 240 j=1,3
240         x(j) = x1(j)
          goto 300
        end if
c
200   continue
      call errcon ('FRCTRN error')
      call ivalut (' At :',3,i1)
      call fvalut (' X  :',3,x)
      call fvalut (' Min:',3,min)
      call fvalut (' Max:',3,max)
      call ivalut (' NSY:',1,ct)
      do i=1,ct
        call vecrtv (x, x1, 1, rt(1,i), rt(10,i))
        call fvalut (' X1 :',3,x1)
      end do
c
      return
c
300   continue
      do 310 i=1,3
310     i2(i) = nint(x(i)*float(grid(i)))- orgn(i)+ 1
c
      return
      end
