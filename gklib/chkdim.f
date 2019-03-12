c
c ===========================================================================
c
      subroutine chkdim (ival,mival,maval,rout,text)
c
      implicit none
c
      integer ival,mival,maval,leng1
c
      character text*(*),rout*(*)
c
code ...
c
      if (ival .lt. mival .or. ival .gt. maval) then
        write (*,6000) text(1:leng1(text)),
     +    rout(1:leng1(rout)),ival,mival,maval
        call endit ('Invalid input')
      end if
c
 6000 format (//
     + ' *** ERROR'/
     + ' *** Invalid value entered for ',a/
     + ' *** In routine ',a/
     + ' *** Value ',i10,' not in range ',i10,' - ',i10/
     + ' *** Execution aborted'/
     + ' *** ERROR'//)
c
      return
      end
