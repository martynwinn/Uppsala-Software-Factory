c
c ==========================================================================
c
      subroutine gklibf (file)
c
c --- GKLIBF (file) - check if GKLIB is defined; if so, prepend to FILE
c
c --- G J Kleywegt @ 950118
c
      implicit none
c
      integer ierr,length,lv,lf
c
      character value*256,file*(*),separator*1,dummy*512
c
c ... SEPARATOR = '/' for UNIX, ']' for VMS
c
      data separator / '/' /
c
code ...
c
c ... check if environment variable GKLIB has been defined
c
      call gknval ('GKLIB',value,ierr)
      if (ierr .ne. 0) return
      if (length(value) .lt. 1) return
c
      lv = length(value)
      if (value(lv:lv) .ne. separator) then
        if (lv .eq. len(value)) return
        lv = lv + 1
        value (lv:lv) = separator
      end if
c
      lf = length(file)
      if ( (lf + lv) .le. len(file) ) then
        dummy = value(1:lv) // file(1:lf)
        file = dummy (1:(lv+lf))
      end if
c
      return
      end
