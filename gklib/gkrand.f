c
c ==========================================================================
c
      subroutine gkrand (xrand,xlo,xhi,iseed)
c
c --- GKRAND (...) - pseudo-random number generator
c
c     ISEED = 0 : generate next number in interval [XLO,XHI]
c                 (initialises with seed MCLOCK() if not init'ed)
c     ISEED > 0 : initialise with seed ISEED & generate XRAND
c     ISEED < 0 : initialise with seed MCLOCK() & generate XRAND
c
c --- G J Kleywegt @ 920702
c
      implicit none
c
      integer iseed,nseed,mclock
c
      logical inited
c
      real xlo,xhi,xrand
      real xval,rand
c
      data inited /.false./
c
      save inited
c
code ...
c
      if ( iseed .lt. 0 .or.
     +    (iseed .eq. 0 .and. (.not.inited)) ) then
        nseed = mclock()
        call srand (nseed)
        write (*,6000) nseed
        inited = .true.
      else if (iseed .gt. 0) then
        nseed = iseed
        call srand (nseed)
        write (*,6000) nseed
        inited = .true.
      end if
c
      xval = rand()
      xrand = xlo + xval*(xhi-xlo)
c
      return
c
 6000 format (' => Random number generator initialised with ',
     +  'seed : ',i10)
c
      end
