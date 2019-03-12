c
c
c
      subroutine pdbinfo (key,line)
c
c ... PDBINFO (KEY,LINE) ... if KEY is interesting, print the line
c
c ... GJK @ 990301
c
      implicit none
c
      integer ll,length
c
      character key*6,line*(*)
c
code ...
c
      if (key .eq. ' ') return
      ll = length (line)
      if (ll .lt. 1) return
c
      if (key .eq. 'HEADER') goto 100
      if (key .eq. 'TITLE ') goto 100
      if (key .eq. 'KEYWDS') goto 100
      if (key .eq. 'EXPDTA') goto 100
      if (key .eq. 'AUTHOR') goto 100
      if (key .eq. 'REVDAT') goto 100
      if (key .eq. 'CRYST1') goto 100
c
      return
c
  100 continue
      write (*,6000) key,line(1:ll)
 6000 format (1x,a6,' : ',a)
c
      return
      end
