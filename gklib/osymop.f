c
c ==========================================================================
c
      subroutine osymop (iunit,file,ierr)
c
      implicit none
c
      integer iunit,ierr,lv,lf,length,leng1
c
      character file*(*),line*256,value*256,separator*1
c
c ... SEPARATOR = '/' for UNIX, ']' for VMS
c
      data separator / '/' /
c
code ...
c
      line = file
c
c ... (1) - try to open using the string as the filename
c
ccc      call textut (' Try to open as :',line)
c
      call xexist (iunit,line,.false.,ierr)
      if (ierr. eq. 0) then
        call textut (' Opened file :',line)
        return
      end if
c
c ... (2) try to prefix with $OSYM if defined
c
      call remspa (line)
c
      call gknval ('OSYM',value,ierr)
      if (ierr .ne. 0) return
c
c ... make sure it's a directory name including trailing '/'
c
      lv = length(value)
      if (value(lv:lv) .ne. separator) then
        if (lv .eq. len(value)) return
        lv = lv + 1
        value (lv:lv) = separator
      end if
c
      call textut (' OSYM :',value)
c
      lf = length(line)
      if ( (lf + lv) .gt. len(line) ) return
      line = value(1:lv) // line(1:lf)
c
ccc      call textut (' Try to open as :',line)
c
      call xexist (iunit,line,.false.,ierr)
      if (ierr .eq. 0) then
        call textut (' Opened file :',line)
        file = line
        return
      end if
c
c ... (3) try to suffix with ".sym"
c
      call locase (line(lv+1:lv+lf))
      line = line(1:leng1(line))//'.sym'
c
ccc      call textut (' Try to open as :',line)
c
      call xexist (iunit,line,.false.,ierr)
      if (ierr .eq. 0) then
        call textut (' Opened file :',line)
        file = line
        return
      end if
c
c ... (4) try to suffix with ".o"
c
      line = line(1:leng1(line)-4)//'.o'
c
ccc      call textut (' Try to open as :',line)
c
      call xexist (iunit,line,.false.,ierr)
      if (ierr .eq. 0) then
        call textut (' Opened file :',line)
        file = line
        return
      end if
c
      call errcon ('Cannot find O symm-op file anywhere')
c
      ierr = -1
c
      return
      end
