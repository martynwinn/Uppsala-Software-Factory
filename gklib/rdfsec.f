c
c
c
      subroutine rdfsec(lun,isec,rho,nu,nv,scale,plus)
C     ================================================
c
C Read formatted section from unit lun
c
      implicit none
c
      integer lun,isec,nu,nv
      real rho(nu,nv)
      real scale,plus
c
      integer maxlin
      parameter (maxlin=1000)
      integer i,j,line(maxlin)
c
code ...
c
      read (lun,6001) isec
 6001 format(/7X,I8/)
c
      do 10, j=1,nv
         read(lun,6010) (line(i),i=1,nu)
 6010    format(20I4)
         do 20, i=1,nu
            rho(i,j) = (float(line(i))-plus)/scale
 20      continue
 10   continue
c
      return
      end
