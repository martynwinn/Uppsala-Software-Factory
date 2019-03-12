c
c=====================================================
c
      subroutine subchr (line,char1,char2,nocc)
c
c ... SUBCHR - substitute (and/or count) CHAR1 to CHAR2 in string LINE
c
      implicit none
c
      integer length,ii,jj,nocc
c
      character line*(*),char1*1,char2*1
c
code ...
c
      ii = length (line)
      nocc = 0
c
      if (ii .gt. 0) then
        do jj=1,ii
          if (line(jj:jj) .eq. char1) then
            line(jj:jj) = char2
            nocc = nocc + 1
          end if
        end do
      end if
c
      return
      end
