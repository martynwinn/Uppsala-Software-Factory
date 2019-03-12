c
c
c
      subroutine c2car4 (nat,xyz,int,ref)
c
c ... C2CAR4 - convert internal to Cartesian coordinates
c
c ... RETAINS THE COORDINATES OF THE FIRST THREE ATOMS !!!
c
      implicit none
c
      real twopi,pi,degtor,rtodeg,tiny
c
      parameter (twopi=6.2831853071796)
      parameter (pi=0.50*twopi)
      parameter (degtor=twopi/360.0)
      parameter (rtodeg=360.0/twopi)
      parameter (tiny=1.0E-9)
c
      real xyz(3,*),int(3,*)
      real xa,xb,ya,yb,za,zb,xyb,xpa,zpa,xpb,zpb,costh,sinth
      real rbc,cosph,sinph,xqa,zqa,a,d,h,coskh,sinkh,sina,cosa
      real sind,cosd,xd,yd,zd,xpd,ypd,zpd,xqd,yqd,zqd,xrd,zrd
      real ypa,yza
c
      integer ref(3,*),nat,i,n1,n2,n3,k,nskip,nzero
c
code ...
c
      if (nat .lt. 4) return
c
ccc      print *,ref(1,1),ref(1,2),ref(1,3)
c
      ref (1,1) = 1
      ref (1,2) = 1
      ref (1,3) = 2
c
      do i=4,nat
        ref (1,i) = - ref(1,i)
      end do
c
      nzero = 0
c
  100 continue
c
      nskip = 0
c
      do i=4,nat
c
        if (ref(1,i) .gt. 0) goto 110
c
ccc        if (i .lt. 10) print *,i,ref(1,i),ref(2,i),ref(3,i)
ccc        if (i .lt. 10) print *,i,int(1,i),int(2,i),int(3,i)
c
        n3 = -ref(1,i)
        n2 = ref(2,i)
        n1 = ref(3,i)
c
        if (ref(1,n1) .lt. 0 .or. ref(1,n2) .lt. 0 .or.
     +      ref(1,n3) .lt. 0) then
          nskip = nskip + 1
          goto 110
        end if
c
        d =  int(1,i)
        a =  int(2,i)
        h = -int(3,i)
c
        xa = xyz(1,n1) - xyz(1,n3)
        ya = xyz(2,n1) - xyz(2,n3)
        za = xyz(3,n1) - xyz(3,n3)
c
        xb = xyz(1,n2) - xyz(1,n3)
        yb = xyz(2,n2) - xyz(2,n3)
        zb = xyz(3,n2) - xyz(3,n3)
c
        xyb = sqrt(xb**2 + yb**2)
        k = 1
c
        if (xyb .lt. 0.1) then
          k = 0
c
          xpa = za
          zpa = -xa
          xa = xpa
          za = zpa
c
          xpb = zb
          zpb = -xb
          xb = xpb
          zb = zpb
c
          xyb = sqrt(xb**2 + yb**2)
        end if
c
        if (abs(xyb) .lt. tiny) then
          call errcon ('C2CAR4 - Prevent division by zero !')
          call rvalut (' Divide by XYB :',1,xyb)
          call rvalut (' Instead use   :',1,tiny)
          xyb = tiny
          nzero = nzero + 1
        end if
c
        costh = xb/xyb
        sinth = yb/xyb
        xpa = xa*costh + ya*sinth
        ypa = ya*costh - xa*sinth
c
        rbc = sqrt (xb**2 + yb**2 + zb**2)
c
        if (abs(rbc) .lt. tiny) then
          call errcon ('C2CAR4 - Prevent division by zero !')
          call rvalut (' Divide by RBC :',1,rbc)
          call rvalut (' Instead use   :',1,tiny)
          rbc = tiny
          nzero = nzero + 1
        end if
c
        sinph = zb/rbc
        cosph = sqrt (1.0 - sinph**2)
        xqa = xpa*cosph + za*sinph
        zqa = za*cosph - xpa*sinph
c
        yza = sqrt (ypa**2 + zqa**2)
c
        if (abs(yza) .lt. tiny) then
          call errcon ('C2CAR4 - Prevent division by zero !')
          call rvalut (' Divide by YZA :',1,yza)
          call rvalut (' Instead use   :',1,tiny)
          yza = tiny
          nzero = nzero + 1
        end if
c
        coskh = ypa/yza
        sinkh = zqa/yza
c
        sina = sin(degtor*a)
        cosa = cos(degtor*a)
        sind = sin(degtor*h)
        cosd = cos(degtor*h)
c
        xd = d*cosa
        yd = d*sina*cosd
        zd = d*sina*sind
c
        ypd = yd*coskh - zd*sinkh
        zpd = zd*coskh + yd*sinkh
        xpd = xd*cosph - zpd*sinph
c
        zqd = zpd*cosph + xd*sinph
        xqd = xpd*costh - ypd*sinth
        yqd = ypd*costh + xpd*sinth
c
        if (k .ne. 1) then
          xrd = -zqd
          zrd = xqd
          xqd = xrd
          zqd = zrd
        end if
c
        xyz (1,i) = xyz(1,n3) + xqd
        xyz (2,i) = xyz(2,n3) + yqd
        xyz (3,i) = xyz(3,n3) + zqd
c
        ref (1,i) = - ref(1,i)
c
ccc        if (i .lt. 10) print *,i,xyz(1,i),xyz(2,i),xyz(3,i)
ccc        if (i .lt. 10) print *
c
  110   continue
c
      end do
c
      if (nskip .gt. 0) then
        call jvalut (' C2CAR4 - Iterate - Nskip :',1,nskip)
        goto 100
      end if
c
      if (nzero .gt. 0) then
        call errcon ('C2CAR4 - Division-by-zero errors !')
        call jvalut (' >>> Nr of errors :',1,nzero)
        call prompt (' >>> Do not trust the results !')
        call prompt (' >>> Please inform Gerard !')
      end if
c
      return
      end
