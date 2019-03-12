c
c
c
      subroutine extint (name,ivalue)
c
c ... EXTINT - check if an INTeger has been defined EXTernally
c
c ... check first if NAME is an environment variable; if so, read
c     the value of IVALUE
c     then check if it is a command-line argument; if so, read the
c     value of IVALUE
c
      implicit none
c
      integer ivalue,ierr,idum,ll,length,i
c
      character name*(*),line*256
c
code ...
c
      ll = length(name)
      if (ll .lt. 1) return
c
c ... (1) check environment variable NAME
c
      call gknval (name,line,ierr)
      if (ierr .eq. 0) then
        call str2i (line,idum,ierr)
        if (ierr .eq. 0) ivalue = idum
      end if
c
c ... (2) check if NAME was a command-line argument
c
      call gknarg (idum)
      if (idum .gt. 1) then
        do i=1,idum-1
          call gkgarg (i,line,ierr)
          if (ierr .eq. 0) then
            call upcase (line)
            if (line(1:ll) .eq. name(1:ll)) then
              call gkgarg (i+1,line,ierr)
              if (ierr .eq. 0) then
                call str2i (line,idum,ierr)
                if (ierr .eq. 0) ivalue = idum
              end if
              goto 999
            end if
          end if
        end do
      end if
c
  999 continue
c
      return
      end
