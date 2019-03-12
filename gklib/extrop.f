c
c ========================================================================
c
      subroutine extrop (option,nopt,maxopt,optpar,ierr)
c
c ... EXTROP - extract option & parameters from a string and put
c              them into an array of strings
c
c ... Gerard Kleywegt @ 920511/940103,0831
c
      implicit none
c
      integer nopt,maxopt,ierr,length,lopt,iold,inow,i
c
      character option*(*),optpar(maxopt)*(*),tab*1
c
code ...
c
      ierr = 0
c
      lopt = length (option)
      nopt = 0
      tab = char(9)
c
      if (lopt .gt. 0) then
        inow = 0
c
   10   continue
        inow = inow + 1
c
        if (inow .gt. lopt) goto 20
c
        if (option (inow:inow) .eq. '"') then
c
          nopt = nopt + 1
          if (nopt .gt. maxopt) goto 99
          iold = inow + 1
          do i=inow+1,lopt
            if (option(i:i) .eq. '"') then
              inow = i-1
              goto 15
            end if
          end do
          inow = lopt
   15     optpar (nopt) = option (iold:inow)
          inow = inow + 1
          goto 10
c
        else if (option (inow:inow) .ne. ' ' .and.
     +           option (inow:inow) .ne. tab) then
c
          nopt = nopt + 1
          if (nopt .gt. maxopt) goto 99
          iold = inow
          do i=inow+1,lopt
            if (option(i:i) .eq. ' ' .or.
     +          option(i:i) .eq. tab) then
              inow = i-1
              goto 25
            end if
          end do
          inow = lopt
   25     optpar (nopt) = option (iold:inow)
          goto 10
c
        end if
c
        goto 10
c
      end if
c
   20 continue
      return
c
c ... too many items on input line
c
   99 continue
      call errcon ('EXTROP - too many items on input line')
      call jvalut (' Maximum allowed :',1,maxopt)
      call prompt (' Rest of line not parsed !!!')
      nopt = nopt - 1
c
      return
      end
