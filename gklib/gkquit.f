c
c ===========================================================================
c
      subroutine gkquit
c
      implicit none
c
c --- GKQUIT => print goodbye mesage etc.
c
c --- G J Kleywegt @ 920403
c
      real x,y,z
c
      integer lp,length,np,i,leng1
c
      logical xsocket,alwyn
c
      character myprog*20,myvers*20,mydate*24,myuser*20,line*80
      character esline*200
c
      common /prginf/ myprog,myvers,mydate,myuser
      common /gjktaj/ alwyn
c
code ...
c
      call gkecpu (x,y,z)
c
      lp = max (1, length(myprog))
      np = 75/(lp+5)
c
      if (xsocket()) then
        write (line,'(3(1x,a,f8.1))') 'User - ',y,'Sys - ',
     +    z,'Total - ',x
        esline = ('CPU time (sec) : '//line(1:leng1(line)))
        call oprint (6,esline)
        esline = ( 'Closed socket '//(myprog(1:lp))//' V-'//
     +    myvers(1:leng1(myvers)))
        call oprint (6,esline)
        call oclose (6)
        goto 9999
      end if
c
      write (*,*)
      write (*,*)
      write (*,*) '*** ',((myprog(1:lp)//' *** '),i=1,np)
      write (*,*)
      write (*,*) 'Version - ',myvers(1:leng1(myvers))
      write (*,*) 'Started - ',mydate(1:leng1(mydate))
      call gkdate (mydate)
      call pretty (mydate)
      write (*,*) 'Stopped - ',mydate(1:leng1(mydate))
      write (*,*)
      write (*,*) 'CPU-time taken :'
      write (*,'(3(1x,a,f8.1))') 'User    - ',y,'Sys    - ',z,
     +  'Total   - ',x
      write (*,*)
      write (*,*) '*** ',((myprog(1:lp)//' *** '),i=1,np)
      write (*,*)
c
      write (*,*)
     +'>>>>>>>>>>>>>> USF .... Uppsala Software Factory <<<<<<<<<<<<<<'
c
      if (alwyn) then
        write (*,*)
     +'>>>> This program: (c) 1992-2009, G J Kleywegt & T A Jones <<<<'
        write (*,*)
     +'>>>> E-mail: gerard@xray.bmc.uu.se or alwyn@xray.bmc.uu.se <<<<'
      else
        write (*,*)
     +'>>>>>>>>>> This program: (c) 1992-2009, G J Kleywegt <<<<<<<<<<'
        write (*,*)
     +'>>>>>>>>>>>>>>>> E-mail: gerard@xray.bmc.uu.se <<<<<<<<<<<<<<<<'
      end if
c
      write (*,*)
     +'>>>>>>>>>>>>>>>>>> http://xray.bmc.uu.se/usf <<<<<<<<<<<<<<<<<<'
      write (*,*)
      write (*,*) '*** ',((myprog(1:lp)//' *** '),i=1,np)
      write (*,*)
c
      if (index(mydate,'Jun 5') .gt. 0) then
        i = length(mydate)
        read (mydate(i-3:i),'(i4)') np
        np = np - 1962
        write (esline,*)
     +    '... Vicious rumours that Gerard turns ',
     +    np,' today are greatly exaggerated ...'
        call pretty (esline)
        write (*,*) esline(1:leng1(esline))
        write (*,*)
      end if
c
      if (index(mydate,'Aug 30') .gt. 0) then
        i = length(mydate)
        read (mydate(i-3:i),'(i4)') np
        np = np - 1947
        write (esline,*)
     +    '... Libellous rumour-mongers insist that Alwyn turns ',
     +    np,' today ...'
        call pretty (esline)
        write (*,*) esline(1:leng1(esline))
        write (*,*)
     +    '... Why not surprise him with an E-mail today ?'
        write (*,*)
     +    '... Just click on the following link:'
        write (*,*)
     +    '... mailto:alwyn@xray.bmc.uu.se?subject=Happy_Birthday'
        write (*,*)
      end if
c
 9999 continue
      stop '... Toodle pip ...'
c
      end
