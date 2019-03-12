c
c ===========================================================================
c
      subroutine pretty (str)
c
c  Make a string pretty by reducing the number of spaces between
c  words to 1, i.e. "Dont    worry,     be        happy"
c  will be changed to "Dont worry, be happy".
c
c When ------- Who ---------------- What -------------------------------
c 27-Nov-1990  Morten Kjeldgaard  Aarhus
c 27-Aug-1991  Morten Kjeldgaard  Bug fix
c 12-Aug-1992  Gerard Kleywegt   Adapted
c
	implicit none
c
	character*(*) str
	integer i, j, m, iword, jword, length
c
code ...
c
	m = length (str)
	if (m .le. 2) return
c
	i = 0
	iword = 1
10	i = i+1
	if (i .gt. m) then
	   if (iword .le. m)  str(iword:) = ' '
	   return
	endif
	if (str(i:i) .eq. ' ') goto 10
c
c gjk @ 930607 - don't go beyond end of string
c
	j = i
        if (j .ge. m) goto 21
c
20	j = j+1
	if (str(j:j) .ne. ' ' .and. j.lt.m) goto 20
c
21	if (str(j:j) .eq. ' ') j = j-1
c
	jword = iword+j-i
	str(iword:jword) = str(i:j)
        if ( (jword+1) .le. m) then
          str(jword+1:jword+1) = ' '
        end if
	iword = jword+2
	i = j
	goto 10
c
	end
