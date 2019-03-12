c
c=====================================================
c
      subroutine opoodb (iunit,file,par,partyp,j,fmt,errcod)
c
c ... OPOODB - open an existing O datablock file
c
c ... Gerard Kleywegt @ 930318, 930813
c
      implicit none
c
      integer iunit,j,errcod,nopt,length
c
      logical xinter
c
      character file*(*),par*(*),partyp*(*),fmt*(*)
      character optpar(4)*80,line*120
c
code ...
c
      errcod = -1
      call textut (' Opening O datablock :',file)
c
      if (length(file) .lt. 1) then
        call errcon ('Empty file name')
        return
      end if
c
      call xopxoa (iunit,file,xinter(),errcod)
      if (errcod .ne. 0) then
        call errcon ('While opening datablock file')
        return
      end if
c
   10 continue
      read (iunit,'(a)',err=999,end=999) line
c
      if (line(1:1) .eq. '!') then
        call prompt (line(2:))
        goto 10
      end if
c
      call upcase (line)
      call extrop (line,nopt,4,optpar,errcod)
c
c ... 930813 - datablocks may now be preceded by comment lines
c              starting with '!' in column 1
c
      if (optpar(1)(1:1) .eq. '!') then
        call prompt (line(2:))
        goto 10
      end if
c
      if (nopt .ne. 4 .or. errcod .ne. 0) then
        errcod = -1
        call errcon ('While reading datablock header')
        call textut (' Header :',line)
        return
      end if
c
      par = optpar (1)
      partyp = optpar (2)
      call str2i (optpar(3),j,errcod)
      if (errcod .ne. 0) goto 999
      fmt = optpar (4)
c
      call textut (' Datablock :',par)
      call textut (' Data type :',partyp)
      call textut (' Number    :',optpar(3))
      call textut (' Format    :',fmt)
c
      errcod = 0
      return
c
  999 continue
      call errcon ('While reading datablock file')
      errcod = -2
c
      return
      end
