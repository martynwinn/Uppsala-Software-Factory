c
c ========================================================================
c
      subroutine xystat (xdata,ydata,ndata,
     +  srmsd,sshap,scorr,srf1,srf2,srf3,srf4,ssf1,ssf2)
c
c ... XYSTAT - statistics pertaining to two real arrays
c
c ... Gerard Kleywegt @ 930420
c
      implicit none
c
      real xdata(*),ydata(*),srmsd,sshap,srf1,srf2,srf3,srf4,scorr,
     +     ssf1,ssf2
      double precision rmsd,shap,rf1,rf2,rf3,rf4,corr,sf1,sf2,avysq
      double precision avalue,bvalue,aabs,babs,avx,avy,avxy,avxsq
      double precision sumfo,sumfc,fofc,fofc1,fofc2,f,x1,x2,x3,dzero
c
      integer i,ndata
c
code ...
c
      dzero = 0.0D0
c
      srmsd = 0.
      sshap = 0.
      scorr = 0.
      srf1 = 0.
      srf2 = 0.
      srf3 = 0.
      srf4 = 0.
      ssf1 = 0.
      ssf2 = 0.
c
      if (ndata .lt. 2) then
        call errcon ('XYSTAT - fewer than 2 data points')
        return
      end if
c
        avx = dzero
        avy = dzero
        avxy = dzero
        avxsq = dzero
        avysq = dzero
        sumfo = dzero
        sumfc = dzero
        fofc = dzero
        fofc1 = dzero
        fofc2 = dzero
        rmsd = dzero
c
        do i=1,ndata
          sumfo = sumfo + abs(xdata(i))
        end do
c
        do i=1,ndata
          sumfc = sumfc + abs(ydata(i))
        end do
c
        if (sumfc .ne. 0.0) then
          sf1 = dabs(sumfo/sumfc)
        else
          sf1 = 1.0
        end if
c
        if (sumfo .ne. 0.0) then
          sf2 = dabs(sumfc/sumfo)
        else
          sf2 = 1.0
        end if
c
        do i=1,ndata
          avalue = xdata(i)
          bvalue = ydata(i)
          aabs = dabs(avalue)
          babs = dabs(bvalue)
          avx = avx+ avalue
          avy = avy+ bvalue
          avxsq = avxsq+ avalue**2
          avysq = avysq+ bvalue**2
          avxy = avxy+ avalue*bvalue
          fofc = fofc + dabs(aabs-babs)
          fofc1 = fofc1 + dabs(aabs-sf1*babs)
          fofc2 = fofc2 + dabs(sf2*aabs-babs)
          rmsd = rmsd + (avalue-bvalue)**2
        end do
c
c ... CORR = correlation coefficient
c
        f = float(ndata)
        x1 = dsqrt(avxsq/f- (avx/f)**2)
        x2 = dsqrt(avysq/f- (avy/f)**2)
        x3 = avxy/f- avx*avy/(f*f)
        if (x1.eq.dzero .or. x2.eq.dzero) then
ccc          if (x3 .eq. dzero) then
ccc            corr = 1.0D0
ccc          else
            corr = dzero
ccc          end if
        else
          corr = (x3)/(x1 * x2)
        end if
c
c ... RF1 = R-factor with respect to XDATA
c
        rf1 = 999.99
        if (sumfo .ne. 0.0) rf1 = fofc/dabs(sumfo)
c
c ... RF2 = R-factor with respect to YDATA
c
        rf2 = 999.99
        if (sumfc .ne. 0.0) rf2 = fofc/dabs(sumfc)
c
c ... RF3 = scaled R-factor wrt XDATA, scale SF1 for YDATA
c
        rf3 = 999.99
        if (sumfo .ne. 0.0) rf3 = fofc1/dabs(sumfo)
c
c ... RF4 = scaled R-factor wrt YDATA, scale SF2 for XDATA
c
        rf4 = 999.99
        if (sumfc .ne. 0.0) rf4 = fofc2/dabs(sumfc)
c
c ... RMSD = Root-Mean-Square Difference
c
        rmsd = dsqrt (rmsd/f)
c
c ... SHAP = shape similarity
c
        if (avxsq .eq. 0.0 .or. avysq .eq. 0.0) then
          shap = 0.0
        else
          shap = avxy / (dsqrt(avxsq)*dsqrt(avysq))
        end if
c
      srmsd = rmsd
      sshap = shap
      scorr = corr
      srf1 = rf1
      srf2 = rf2
      srf3 = rf3
      srf4 = rf4
      ssf1 = sf1
      ssf2 = sf2
c
      return
      end
