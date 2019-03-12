c
c
c
      subroutine rdfhdr(lun,fname,title,nsec,iuvw,mxyz,nw1,nu1,nu2,
     $   nv1,nv2,cell,lspgrp,lmode,rhmin,rhmax,rhmean,rhrms,scale,plus)
C     ===============================================================
c
C  Formatted, read keyworded header, fixed format
c
      implicit none
c
      integer lun
      character*(*) fname
      character*80 title
      integer nsec,iuvw(3),mxyz(3),nw1,nu1,nu2,nv1,nv2,lspgrp,lmode
      real cell(6),rhmin,rhmax,rhmean,rhrms,scale,plus
c
      integer i,j,limits(2,3),kdummy,kfail
      character*1 xyz(3),axes(3)
      data xyz/'X','Y','Z'/
c
code ...
c
      kdummy = 0
      kfail = 0
      call ccpdpn(lun,fname,'READONLY','F',kdummy,kfail)
c
      read(lun,6000) 
 6000 format('MAPEXCHANGE HEADER')
      write(6,6100) 
 6100 format(' MAPEXCHANGE HEADER')
      read(lun,6010) title
 6010 format(/A)
      write(6,6110) title
 6110 format(' TITLE'/1X,A)
      read(lun,6020) axes
 6020 format('AXIS    ',3(7X,A1))
      write(6,6120) axes
 6120 format(' AXIS    ',3(7X,A1))
      read(lun,6030) mxyz
 6030 format('GRID    ',3I8)
      write(6,6130) mxyz
 6130 format(' GRID    ',3I8)
      read(lun,6040) limits
 6040 format('XYZLIM  ',6I8)
      write(6,6140) limits
 6140 format(' XYZLIM  ',6I8)
      read(lun,6050) lspgrp
 6050 format('SPACEGROUP',6X,I8)
      write(6,6150) lspgrp
 6150 format(' SPACEGROUP',6X,I8)
      read(lun,6060) lmode
 6060 format('MODE    ',I8)
      write(6,6160) lmode
 6160 format(' MODE    ',I8)
      read(lun,6070) cell
 6070 format('CELL    ',6F10.3)
      write(6,6170) cell
 6170 format(' CELL    ',6F10.3)
      read(lun,6080) rhmin,rhmax,rhmean,rhrms
 6080 format('RHOLIM  ',4G16.6)
      write(6,6180) rhmin,rhmax,rhmean,rhrms
 6180 format(' RHOLIM  ',4G16.6)
      read(lun,6090) scale,plus
 6090 format('PACK    ',2G16.6)
      write(6,6190) scale,plus
 6190 format(' PACK    ',2G16.6)
      read(lun,6095) 
 6095 format('END HEADER')
      write(6,6195) 
 6195 format(' END HEADER')
c
C Get axis order
      do 10, i=1,3
         do 20, j=1,3
            if (axes(i) .eq. xyz(j)) then
               iuvw(i) = j
            endif
 20      continue
 10   continue
c
C Get limits 
      nu1 = limits(1,iuvw(1))
      nu2 = limits(2,iuvw(1))
      nv1 = limits(1,iuvw(2))
      nv2 = limits(2,iuvw(2))
      nw1 = limits(1,iuvw(3))
      nsec= limits(2,iuvw(3))-nw1+1
c
      return
      end
