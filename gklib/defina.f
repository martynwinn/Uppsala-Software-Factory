c
c ==========================================================================
c
      subroutine defina (myname,fpath,f1,fname,f2,fext,f3,fvers)
c
c --- DEFINA (...) - deduce file name, extension etc. from path name
c
c --- G J Kleywegt @ 911030/920518
c
      implicit none
c
      integer length,i,lp,i1,i2,i3
c
      character*(*) fpath,fname,myname,fext,fvers
      character*1   separator,extension,version,f1,f2,f3
      character     fullnm*255
c
c ... UNIX: separator = /  extension = .  version = NONE
c     VMS : separator = ]  extension = .  version = ;
c
      data separator /'/'/, extension /'.'/, version/';'/
c
code ...
c
      fullnm = myname
      fname = fullnm
      lp = length (fullnm)
      fpath = ' '
      fext = ' '
      fvers = ' '
      f1 = separator
      f2 = extension
      f3 = version
c
c ... first chop off path name (if any)
c
      i1 = 0
      do i=lp,1,-1
        if (fullnm(i:i) .eq. separator) then
          fpath = fullnm (1:i-1)
          i1 = i
          goto 10
        end if
      end do
c
c ... find base name
c
   10 continue
      fname = fullnm (i1+1:lp)
      lp = length(fname)
      i2 = index (fname,extension)
      i3 = index (fname,version)
c
      if (i2 .le. 0 .and. i3 .le. 0) then
c ...   nothing
      else if (i2 .gt. 0 .and. i3 .le. 0) then
        fext  = fname (i2+1:lp)
        fname = fname (1:i2-1)
      else if (i2 .le. 0 .and. i3 .gt. 0) then
        fvers = fname (i3+1:lp)
        fname = fname (1:i3-1)
      else
        fext  = fname (i2+1:i3-1)
        fvers = fname (i3+1:lp)
        fname = fname (1:i2-1)
      end if
c
cc      print *
cc      print *,'FULL ',fullnm(1:length(fullnm))
cc      print *,'PATH ',fpath(1:length(fpath))
cc      print *,'BASE ',fname(1:length(fname))
cc      print *,'EXTN ',fext(1:length(fext))
cc      print *,'VERS ',fvers(1:length(fvers))
cc      print *,'I123 ',i1,i2,i3,lp
c
      return
      end
