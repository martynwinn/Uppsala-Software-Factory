c
c ==========================================================================
c
      subroutine ezdput (unit,
     +                   a,b,c,alpha,beta,gamma,
     +                   nxlo,nylo,nzlo,
     +                   ngx,ngy,ngz,
     +                   nia,nib,nic,
     +                   ndata,xdata,scale,ifmt1,ifmt2)
c
c --- EZDPUT (UNIT,...) - write EZD_MAP stream or file
c
c unit = logical I/O unit
c a,b,c,alpha,beta,gamma = cell constants
c nxlo,nylo,nzlo = grid origin
c ngx,ngy,ngz = grid extent
c nia,nib,nic = nrs of grid points
c ndata = nr of data points (=ngx*ngy*ngz)
c xdata = electron density array
c scale = factor by which density HAS BEEN MULTIPLIED
c ifmt1,ifmt2 = used for FORMAT (e.g., for F12.5,
c               ifmt1 = 12 and ifmt2 = 5)
c
c --- G J Kleywegt @ 920917
c
      implicit none
c
      integer unit,nxlo,nylo,nzlo,ngx,ngy,ngz,nia,nib,nic
      integer ndata,fmt1,fmt2,length,i,i2,j,nfmt,ifmt1,ifmt2
c
      real a,b,c,alpha,beta,gamma,xdata(*),scale
c
      character fmtstr*80,myline*250
      character mydate*24,myuser*20,myprog*20,myvers*20
c
      common /prginf/ myprog,myvers,mydate,myuser
c
code ...
c
c ... the following MUST be the first line in the file
c
      write (unit,6000) 'EZD_MAP'
c
c ... write a comment (which program created this file,
c     when and for which user)
c
      call stamp (myline)
      myline = '! '//myline
      write (unit,6000) myline (1:length(myline))
c
      write (myline,'(a,6(1x,f15.3))')
     +  'CELL',a,b,c,alpha,beta,gamma
      call pretty (myline)
      write (unit,6000) myline (1:length(myline))
c
      write (myline,'(a,3(1x,i15))')
     +  'ORIGIN',nxlo,nylo,nzlo
      call pretty (myline)
      write (unit,6000) myline (1:length(myline))
c
      write (myline,'(a,3(1x,i15))')
     +  'EXTENT',ngx,ngy,ngz
      call pretty (myline)
      write (unit,6000) myline (1:length(myline))
c
      write (myline,'(a,3(1x,i15))')
     +  'GRID',nia,nib,nic
      call pretty (myline)
      write (unit,6000) myline (1:length(myline))
c
      write (myline,'(a,1x,1pe15.8)')
     +  'SCALE',scale
      call pretty (myline)
      write (unit,6000) myline (1:length(myline))
c
      write (unit,6000) 'MAP'
c
      fmt1 = max (15, ifmt1)
      if (7*(fmt1+1) .gt. 250) fmt1 = 34
c
      fmt2 = max (1, min (ifmt2, fmt1 - 2))
cxyz
c ... O 5.9.1 -> fix at SEVEN items per line ...
c
c      nfmt = 80/(fmt1+1)
      nfmt = 7
cxyz
c      myline = ' '
      write (myline,*) '(0P,',nfmt,'(f',fmt1,'.',fmt2,',1x))'
      call remspa (myline)
c      fmtstr = ' '
      fmtstr = myline(1:length(myline))
      call textut (' (NEW)EZD format :',fmtstr)
c
cc      print *,'EZDPUT'
cc      print *,fmt1,fmt2,nfmt,fmtstr
c
      do i=1,ndata,nfmt
        myline = ' '
        i2 = min (ndata, i+nfmt-1)
        write (myline,fmtstr) (xdata(j),j=i,i2)
        call pretty (myline)
        write (unit,6000) myline (1:length(myline))
      end do
c
c ... NOTE: subroutine PRETTY replaces MULTIPLE SPACES by
c           a single space so as to save disk space
c
c ... the last record must be END
c
      write (unit,6000) 'END'
c
      return
c
 6000 format (A)
c
      end
