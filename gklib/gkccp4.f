c
c ==========================================================================
c
      subroutine gkccp4 (fatal)
c
c --- GKCCP4 (log) - check if CCP4_OPEN is defined
c
c --- G J Kleywegt @ 950118
c
      implicit none
c
      integer ierr,length
c
      logical fatal
c
      character name*80,value*256
c
code ...
c
      name = 'CCP4_OPEN'
      call gknval (name,value,ierr)
      if (ierr .ne. 0) goto 900
      if (length(value) .lt. 1) goto 900
c
      return
c
  900 continue
      call errcon ('Environment variable CCP4_OPEN undefined')
      call errcon ('You will *not* be able to write CCP4 maps')
      call errcon ('To fix this, use: setenv CCP4_OPEN UNKNOWN')
      if (fatal) call errstp ('GKCCP4 - Aborting program')
c
      return
      end
