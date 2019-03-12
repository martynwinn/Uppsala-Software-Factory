c
c ==========================================================================
c
      subroutine oprint (unit,text)
c
c --- OPRINT (UNIT,TEXT) - write TEXT to UNIT and flush
c
c --- G J Kleywegt @ 920916
c
      implicit none
c
      integer unit,length,leng1
c
      character text*(*),line*1024
c
code ...
c
      if (length(text) .gt. 0 .and. unit .gt. 0) then
        line = text
        call pretty (line)
        write (unit,'(9a)',err=999)
     +    'PRINT "',line(1:leng1(line)),'"'
        call flusho (unit)
      end if
c
      return
c
  999 continue
      call errint ('OPRINT - could not send text to O')
c
      return
      end
