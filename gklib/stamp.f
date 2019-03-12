c
c ========================================================================
c
      subroutine stamp (line)
c
c ... return a line 'Created by XXX V. YYY at ZZZ for NN'
c
      implicit none
c
      integer leng1
c
      character line*(*),nowish*24
      character myprog*20,myvers*20,mydate*24,myuser*20
c
      common /prginf/ myprog,myvers,mydate,myuser
c
code ...
c
      call gkdate (nowish)
      write (line,*,err=10)
     +  'Created by ',myprog(1:leng1(myprog)),
     +  ' V. ',myvers(1:leng1(myvers)),' at ',
     +  nowish(1:leng1(nowish)),' for ',
     +  myuser(1:leng1(myuser))
      call pretty (line)
      return
c
   10 continue
      line = 'Created at '//nowish
      return
c
      end
