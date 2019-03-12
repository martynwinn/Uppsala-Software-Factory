c
c ===========================================================================
c
c GKSUBS.F - utility routines - Gerard Kleywegt - Uppsala - 1992-2009
c
c ===========================================================================
c
      subroutine gkinit (prognm,version)
c
      implicit none
c
c --- GKINIT (...) => print header etc.
c
c --- G J Kleywegt @ 920403/920915
c
      real x,y,z
c
      integer lp,length,np,i,iproc,leng1
c
      logical ltty,xsocket,alwyn
c
      character prognm*(*),version*(*),myprog*20,myvers*20
      character mydate*24,myuser*20,ntty*20,mymode*20
      character myhost*64,mypid*20,esline*128
c
      common /prginf/ myprog,myvers,mydate,myuser
      common /gjktaj/ alwyn
c
      data myprog,myvers,mydate,myuser,mymode /5*'???'/
c
code ...
c
      alwyn = .false.
      goto 10
c
      entry gainit (prognm,version)
c
      alwyn = .true.
c
   10 continue
      call gkdcpu (x,y,z)
      call gkdate (mydate)
      call pretty (mydate)
      myprog = prognm
      myvers = version
      call gkuser (myuser)
      if (length(myuser) .lt. 1) myuser = 'A. Nonymous'
      call gkatty (5,ltty,ntty)
      call gkargs ()
      call gkpath ()
      call gkmode (mymode)
      call gkhost (myhost)
      call gkpid  (iproc)
      write (mypid,*) iproc
      call remspa (mypid)
c
      lp = max (1, length(myprog))
      np = 75/(lp+5)
c
      if (xsocket()) then
        esline = ('Opened socket '//(myprog(1:lp))//' V-'//
     +    myvers(1:leng1(myvers))) 
        call oprint (6,esline)
        return
      end if
c
      write (*,*)
      write (*,*) '*** ',((myprog(1:lp)//' *** '),i=1,np)
      write (*,*)
      write (*,*) 'Version  - ',myvers(1:leng1(myvers))
      if (alwyn) then
        write (*,*)
     +    '(c) 1992-2009 Gerard J. Kleywegt & T. Alwyn Jones,',
     +    ' BMC, Uppsala (SE)'
      else
        write (*,*)
     +    '(c) 1992-2009 Gerard J. Kleywegt, Dept. Cell Mol. Biol.,',
     +    ' Uppsala (SE)'
      end if
      write (*,*) 
     +  'User I/O - routines courtesy of Rolf Boelens, Univ.',
     +  ' of Utrecht (NL)'
      write (*,*) 
     +  'Others   - T.A. Jones, G. Bricogne, Rams, W.A. Hendrickson'
      write (*,*) 
     +  'Others   - W. Kabsch, CCP4, PROTEIN, E. Dodson, etc. etc.'
      write (*,*)
      write (*,*) 'Started  - ',mydate(1:leng1(mydate))
      write (*,*) 'User     - ',myuser(1:leng1(myuser))
      write (*,*) 'Mode     - ',mymode(1:leng1(mymode))
      if (length(myhost) .gt. 0)
     +  write (*,*) 'Host     - ',myhost(1:leng1(myhost))
      if (iproc .gt. 0)
     +  write (*,*) 'ProcID   - ',mypid(1:leng1(mypid))
c
      if (ltty) then
        write (*,*) 'Tty      - ',ntty(1:leng1(ntty))
      else
        write (*,*) 'Not using a tty as input device'
      end if
      write (*,*)
      write (*,*) '*** ',((myprog(1:lp)//' *** '),i=1,np)
      write (*,*)
c
      call gkrefs (myprog,i)
      if (i .gt. 0) then
        write (*,*)
        write (*,*) '*** ',((myprog(1:lp)//' *** '),i=1,np)
        write (*,*)
      end if
c
      return
      end
