c
c ========================================================================
c
      integer function leng1 (string)
c
c ... on ALPHAs, writing string(1:0) gives a core dump
c     therefore, when writing a string (even an empty one)
c     this should be string(1:1)
c
c ... gjk @ 960415
c
      implicit none
c
      integer length,idum
c
      character string*(*)
c
code ...
c
      idum = length (string)
      if (idum .lt. 1) idum = 1
      leng1 = idum
c
      return
      end
