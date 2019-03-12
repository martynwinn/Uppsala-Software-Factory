c
c ... subroutines for DATAMAN
c
c ... Gerard Kleywegt
c
c=======================================================================--------
c
      SUBROUTINE   W_SCALE_2 ( FOBS 1, DST SQ1, NUM ORB 1, N OBS 1,
     .                         FOBS 2, DST SQ2, NUM ORB 2, N OBS 2,
     .                DLO SCL, DHI SCL, STEP, N BIN MX, nbinmn, NBINS,
     .                         NREFL 1, NFULL 1, NREFL 2, NFULL 2,
     .                         SUM SQ1, SUM SQ2, SCALE, BTEMP )
c
      implicit NONE
c
      integer nobs1,nobs2,nbinmx,nbinmn
C
      INTEGER   NUM ORB 1 (N OBS 1), NUM ORB 2 (N OBS 2)
      INTEGER   N REFL 1 (N BIN MX), N REFL 2 (N BIN MX)
      INTEGER   N FULL 1 (N BIN MX), N FULL 2 (N BIN MX)
C
      REAL      F OBS 1  (N OBS 1) , F OBS 2  (N OBS 2)
      REAL      DST SQ 1 (N OBS 1) , DST SQ 2 (N OBS 2)
      REAL      SUM SQ1 (N BIN MX ), SUM SQ2 (N BIN MX)
c
      real dloscl,dhiscl,step,scale,btemp
c
      integer nbins,ismin,ismax
C
C
C         Determine a crude scale factor SCALE and temperature factor
C     BTEMP such that  FOBS1  and  FOBS2 * SCALE * EXP(-BTEMP*STOLSQ) 
C     be equivalent. The two sets of amplitudes do not have to belong
C     to the same crystal form.
C
C         Written by G. Bricogne, 18 March 1993, at BMC, Uppsala.
C
C-----------------------------------------------------------------------
C
C     Arguments.
C
C         F OBS  i : list of amplitudes for reflexions of crystal form i
C         DST SQ i : d* squared ( = 4 * (sin(theta)/lambda)**2 for 
C                    reflexions of crystal form i.
C         NUMORB i : orbital multiplicities of reflexions of crystal 
C                    form i. This is defined as the number of reflexions
C                    equivalent to a given unique reflexion under 
C                    symmetry expansion and Friedel expansion. Note that
C                    NUM ORB i = 2 * #Gi / eps_hkl_i for acentric hkl_i
C                              =     #Gi / eps_hkl_i for centric  hkl_i
C                    where #Gi is the number of symmetry operations in 
C                    the space group Gi of form i, and eps_hkl_i is the 
C                    number of operations of Gi under which hkl_i is 
C                    invariant.
C         N OBS i  : number of observed amplitudes in the list for form
C                    no. i .
C
C         DLO SCL  : minimum value of d* squared for resolution range.
C         DHI SCL  : maximum value of d* squared for resolution range.
C
C         STEP     : increment in d* squared values used to define bins.
C                    This is chosen by the routine but its value is 
C                    made available (e.g. for plotting).
C         N BIN MX : maximum number of resolution bins allowed by the 
C                    array dimensions in the calling program.
C         N BINS   : number of bins actually used in the accumulation
C                    of intensity statistics.
C
C         N REFL i : number of unique reflexions for crystal form i 
C                    found within each resolution bin.
C         N FULL i : sum of orbital multiplicities for the reflexions
C                    of crystal form no. i within each bin.
C         SUM SQ i : sum or orbital-multiplicity-weighted intensities
C                    within each bin.
C
C         SCALE    : final scale factor.
C         TEMP     : final temperature factor.
C
C
C-----------------------------------------------------------------------
C         Note: it may be tempting to fudge the orbital multiplicities
C     NUM ORB i  (e.g. by using a value of 2 for acentrics and a value
C     of 1 for centrics) but this temptation should be resisted if at 
C     least one of the forms has a high-multiplicity space group.
C-----------------------------------------------------------------------
C
C
C         Define the resolution bins so that their number does not 
C     exceed  N BIN MX .
C
C
      if (step .lt. 0.0001) then
        STEP  =  0.0025
      end if
c
      N BINS  =  INT ( 0.25 * DHI scl / STEP ) + 1
ccc      print *,'nbins ',nbins,step
      IF ( N BINS .GT. N BIN MX )      THEN
ccc         STEP  =  ( FLOAT ( N BIN MX ) - 1.01 ) / ( 0.25 * DHI scl )
         N BINS  =  N BIN MX
         step = ( 0.25 * DHI scl ) / ( FLOAT ( N BIN MX ) - 1.01 )
      END IF
ccc      print *,'nbins ',nbins,step
c
      if (nbins .lt. nbinmn) then
ccc        STEP  =  ( FLOAT ( N BIN Mn ) - 1.01 ) / ( 0.25 * DHI scl )
        N BINS  =  N BIN Mn
        step = ( 0.25 * DHI scl ) / ( FLOAT ( N BIN Mn ) - 1.01 )
      end if
ccc      print *,'nbins ',nbins,step
c
cc      call jvalut (' Bins :',1,nbins)
cc      call fvalut (' Step :',1,step)
C
C         Accumulate the intensity statistics in these bins for the two 
C     sets of amplitudes.
C
      CALL   DST BIN ( F OBS 1, DST SQ1, NUM ORB 1, N OBS 1,
     .                 DLO SCL, DHI SCL, STEP, N BINS,
     .                 N REFL1, NFULL 1, SUM SQ1 )
C
      CALL   DST BIN ( F OBS 2, DST SQ2, NUM ORB 2, N OBS 2,
     .                 DLO SCL, DHI SCL, STEP, N BINS,
     .                 N REFL2, NFULL 2, SUM SQ2 )
c
ccc      do i=1,nbins
ccc        print *,i,nrefl1(i),nrefl2(i)
ccc      end do
C
C         Determine relative scale and temperature factor by a Wilson
C     plot.
C
      IS MIN  =  INT ( 0.25 * DLO scl / STEP ) + 1
      IS MAX  =  INT ( 0.25 * DHI scl / STEP ) + 1
C
C
      CALL   WILSON2 ( N BINS, STEP, IS MIN, IS MAX, 
     .                 SUM SQ1, SUM SQ2, NFULL 1, NFULL 2, 
     .                 SCALE, BTEMP )
C
C
      RETURN
      END
c
C=======================================================================--------
c
      SUBROUTINE   DST BIN ( F OBS  , DST SQ , NUM ORB, N OBS,
     .                       DLO SCL, DHI SCL, STEP, N BINS,
     .                       N REFL , N FULL , SUM SQ )
c
      implicit NONE
c
      integer nobs,nbins
C
      INTEGER   NUM ORB (N OBS), N REFL (N BINS), N FULL (N BINS)
C
      REAL      F OBS (N OBS), DST SQ (N OBS), SUM SQ (N BINS)
c
      real dloscl,dhiscl,step,stolsq
c
      integer is,i,m
C
C
C         Compiles intensity statistics in resolution bins.
C
C
      DO IS = 1, N BINS
         N REFL (IS) =  0
         N FULL (IS) =  0
         SUM SQ (IS) = 0.0
      END DO
C
C
      DO I = 1, N OBS
C
         IF ( ( DSTSQ (I) .GE. DLO SCL ) .AND. 
     .        ( DSTSQ (I) .LE. DHI SCL )       )      THEN
C
            STOLSQ  =  0.25 * DSTSQ (I)
C
            M  =  NUMORB (I)
            IS =  INT ( STOLSQ / STEP ) + 1
            N REFL (IS)  =  N REFL (IS) + 1
            N FULL (IS)  =  N FULL (IS) + M
C
            SUM SQ (IS)  =  SUM SQ (IS)  +  M * FOBS(I)**2
C
         END IF
C
      END DO
C
      RETURN
      END
c
C=======================================================================--------
c
      SUBROUTINE   WILSON2  ( N BINS, STEP, IS MIN, IS MAX,
     .                       SUM SQ1, SUM SQ2, N COUNT 1, N COUNT 2,
     .                       SCALE, BTEMP )
C
C         Determine a crude scale factor SCALE and temperature factor
C     BTEMP such that  F1  and  F2 * SCALE * EXP(-BTEMP*STOLSQ) be
C     equivalent. The two data sets need not belong to the same crystal
C     form.
C
C
      implicit NONE
c
      integer nbins
c
      INTEGER   N COUNT 1 ( N BINS ), N COUNT 2 ( N BINS )
c
      REAL      SUM SQ1   ( N BINS ), SUM SQ2   ( N BINS )
c
      real step,scale,btemp,sumx,sumy,sumxx,sumxy,sumw,x,y,w
      real xw,yw,slope,ordint
c
      integer ismin,ismax,is,n1,n2,icnt
C
C
      SUM X  = 0.0
      SUM Y  = 0.0
      SUM XX = 0.0
      SUM XY = 0.0
      SUM W  = 0.0
      icnt = 0
C
C
      DO IS = IS MIN, IS MAX
         N1 = N COUNT 1 (IS)
         N2 = N COUNT 2 (IS)
         IF ( ( N1 .GT. 0 ) .AND. ( N2 .GT. 0 ) )      THEN
c
c      print *,'IS ',is,n1,n2
c
            icnt = icnt + 1
            X   =  ( IS - 0.5 ) * STEP
            Y   =  ALOG ( ( SUM SQ2 (IS) / FLOAT (N2) )
     .                  / ( SUM SQ1 (IS) / FLOAT (N1) ) )
            W   =  FLOAT ( N1 * N2 ) / FLOAT ( N1 + N2 )
            XW  =  X * W
            YW  =  Y * W
            SUM X   =  SUM X  + XW
            SUM Y   =  SUM Y  + YW
            SUM XX  =  SUM XX + X*XW
            SUM XY  =  SUM XY + X*YW
            SUM W   =  SUM W  + W
         END IF
      END DO
C
C
      if (icnt .lt. 2) then
        call errcon ('WILSON2 - Not enough filled bins !')
        btemp = 0.0
        scale = 1.0
      else
        if (icnt .lt. 10) then
          call prompt (' WARNING - Fewer than 10 bins used')
        end if
        SLOPE  = ( sum W * SUM XY  -  SUM X * SUM Y )
     .         / ( sum W * SUM XX  -  SUM X * SUM X )
        ORDINT = ( SUM Y - SLOPE * SUM X ) / sum W
        BTEMP  =   0.5 * SLOPE
        SCALE  =   EXP ( - 0.5 * ORDINT )
      end if
C
C
      WRITE(6,1000)   icnt,SCALE, BTEMP
C
C
      RETURN
C
C
 1000 FORMAT(1X,'F''s put on same scale+temp by Wilson plot:',
     +      /1x,' Nr of filled bins = ',i15,
     .      /1X,'          W SCALE  = ',1PE15.5,0P,
     .      /1X,'          W BTEMP  = ',F15.3)
C
C
      END
c
c
c
      logical function oklaue (hkl,laue)
c
      implicit none
c
      integer hkl(3),laue,h,k,l
c
code ...
c
      oklaue = .false.
c
      h = hkl(1)
      k = hkl(2)
      l = hkl(3)
c
      if (laue .eq. 1) then
c
c - LAUE = 1, 1bar,   hkl:h>=0   0kl:k>=0   00l:l>=0
c
        if (h .eq. 0 .and. k .eq. 0) then
          if (l .lt. 0) return
        else if (h .eq. 0) then
          if (k .lt. 0) return
        else
          if (h .lt. 0) return
        end if
        goto 10
c
      else if (laue .eq. 2) then
c
c - LAUE = 2, 1bar,   hkl:k>=0   h0l:l>=0   h00:h>=0
c
        if (k .eq. 0 .and. l .eq. 0) then
          if (h .lt. 0) return
        else if (k .eq. 0) then
          if (l .lt. 0) return
        else
          if (k .lt. 0) return
        end if
        goto 10
c
      else if (laue .eq. 3) then
c
c - LAUE = 3, 1bar,   hkl:l>=0   hk0:h>=0   0k0:k>=0
c
        if (h .eq. 0 .and. l .eq. 0) then
          if (k .lt. 0) return
        else if (l .eq. 0) then
          if (h .lt. 0) return
        else
          if (l .lt. 0) return
        end if
        goto 10
c
      else if (laue .eq. 4) then
c
c 940711 - added k>=0 to hk0
c - LAUE = 4, 2/m,    hkl:k>=0, l>=0     hk0:h>=0, k>=0
c
        if (l .eq. 0) then
          if (h .lt. 0) return
          if (k .lt. 0) return
        else
          if (k .lt. 0) return
          if (l .lt. 0) return
        end if
        goto 10
c
      else if (laue .eq. 5) then
c
c 940711 - added l>=0 to 0kl
c - LAUE = 5, 2/m,    hkl:h>=0, l>=0     0kl:k>=0, l>=0   (2-nd sett.)
c
        if (h .eq. 0) then
          if (k .lt. 0) return
          if (l .lt. 0) return
        else
          if (h .lt. 0) return
          if (l .lt. 0) return
        end if
        goto 10
c
      else if (laue .eq. 6) then
c
c - LAUE = 6, mmm,    hkl:h>=0, k>=0, l>=0
c
        if (h .lt. 0) return
        if (k .lt. 0) return
        if (l .lt. 0) return
        goto 10
c
      else if (laue .eq. 7) then
c
c - LAUE = 7, 4/m,    hkl:h>=0, k>0, l>=0 with  k>=0 for h=0
c
        if (h .lt. 0) return
        if (l .lt. 0) return
        if (h .eq. 0) then
          if (k .lt. 0) return
        else
          if (k .le. 0) return
        end if
        goto 10
c
      else if (laue .eq. 8) then
c
c - LAUE = 8, 4/mmm,  hkl:h>=0, h>=k>=0, l>=0
c
        if (l .lt. 0) return
        if (k .lt. 0) return
        if (h .lt. k) return
        goto 10
c
      else if (laue .eq. 9) then
c
c - LAUE = 9, 3bar,   hkl:h>=0, k<0, l>=0 including 00l
c
        if (h .eq. 0 .and. k .eq. 0) then
c
        else
          if (l .lt. 0) return
          if (k .ge. 0) return
          if (h .lt. 0) return
        end if
        goto 10
c
      else if (laue .eq. 10) then
c
c - LAUE = 10, 3bar,  hkl:h>=0, k>0  including  00l:l>0
c
        if (h .eq. 0 .and. k .eq. 0) then
          if (l .le. 0) return
        else
          if (h .lt. 0) return
          if (k .le. 0) return
        end if
        goto 10
c
      else if (laue .eq. 11) then
c
c - LAUE = 11, 3barm, hkl:h>=0, k>=0 with k<=h; if h=k l>=0
c
        if (h .lt. 0) return
        if (k .lt. 0) return
c
c ... bug fix here, 930726
c
        if (k .gt. h) return
        if (k .eq. h) then
          if (l .lt. 0) return
        end if
c
        goto 10
c
      else if (laue .eq. 12) then
c
c - LAUE = 12, 6/m,   hkl:h>=0, k>0, l>=0  with  k>=0 for h=0
c
        if (l .lt. 0) return
        if (h .lt. 0) return
        if (h .eq. 0) then
          if (k .lt. 0) return
        else
          if (k .le. 0) return
        end if
        goto 10
c
      else if (laue .eq. 13) then
c
c - LAUE = 13, 6/mmm, hkl:h>=0, h>=k>=0, l>=0
c
        if (l .lt. 0) return
        if (k .lt. 0) return
        if (h .lt. k) return
        goto 10
c
      else if (laue .eq. 14) then
c
c - LAUE = 14, m3,    hkl:h>=0, k>=0, l>=0 with l>=h,
c                         k>=h if l=h,
c                         k>h if l>h
c
        if (h .lt. 0) return
        if (k .lt. 0) return
        if (l .lt. h) return
        if (l .eq. h .and. k .lt. h) return
        if (l .gt. h .and. k .le. h) return
        goto 10
c
      else if (laue .eq. 15) then
c
c - LAUE = 15, m3m,   hkl:k>=l>=h>=0
c
        if (h .lt. 0) return
        if (l .lt. h) return
        if (k .lt. l) return
        goto 10
c
      else
        call errcon ('Invalid Laue group number')
        call jvalut (' Laue group number :',1,laue)
      end if
c
c ... if here, the reflection is okay
c
   10 continue
      oklaue = .true.
c
      return
      end
c
c
c
      subroutine hklmat (hnew,knew,lnew,mathkl,ierr)
c
c ... figure out re-indexing matrix from text expressions
c
      implicit none
c
      integer mathkl(3,3),ierr,i,j,sign,ll,length,mult
c
      character hnew*(*),knew*(*),lnew*(*)
      character hklnew(3)*20,hklnam(3)*1
c
      data hklnam /'H','K','L'/
c
code ...
c
      ierr = -1
c
      hklnew (1) = hnew
      hklnew (2) = knew
      hklnew (3) = lnew
c
   10 format (' New ',a1,' = ',a)
   20 format ('     ',a1,' = ',i2,'*H + ',i2,'*K + ',i2,'*L')

      do i=1,3
c
        call remspa (hklnew(i))
        ll=length(hklnew(i))
        if (ll.lt.1) then
          call errcon ('No expression for new '//hklnam(i))
          return
        end if
c
        call upcase (hklnew(i))
        write (*,10) hklnam(i),(hklnew(i)(1:ll))
c
        do j=1,ll
          if (index('HKL+-123456',hklnew(i)(j:j)) .lt. 1) then
            call errcon ('Unrecognised character '//hklnew(i)(j:j))
            return
          end if
        end do
c
        do j=1,3
          mathkl(j,i)=0
        end do
        mathkl(i,i) = 0
c
        sign = 1
        mult = 1
        do j=1,ll
          if (hklnew(i)(j:j) .eq. '+') then
            sign = 1
          else if (hklnew(i)(j:j) .eq. '-') then
            sign = -1
          else if (hklnew(i)(j:j) .eq. '2') then
            mult = 2
          else if (hklnew(i)(j:j) .eq. '3') then
            mult = 3
          else if (hklnew(i)(j:j) .eq. '4') then
            mult = 4
          else if (hklnew(i)(j:j) .eq. '5') then
            mult = 5
          else if (hklnew(i)(j:j) .eq. '6') then
            mult = 6
          else if (hklnew(i)(j:j) .eq. 'H') then
            mathkl (1,i) = sign*mult
          else if (hklnew(i)(j:j) .eq. 'K') then
            mathkl (2,i) = sign*mult
          else if (hklnew(i)(j:j) .eq. 'L') then
            mathkl (3,i) = sign*mult
          end if
        end do
c
        write (*,20) hklnam(i),(mathkl(j,i),j=1,3)
c
      end do
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine reindx (nhkl,hkl,m,ierr)
c
      implicit none
c
      integer nhkl,hkl(3,*),m(3,3),ierr,h,k,l,i
c
code ...
c
      do i=1,nhkl
        h=hkl(1,i)
        k=hkl(2,i)
        l=hkl(3,i)
        hkl(1,i)=m(1,1)*h+m(2,1)*k+m(3,1)*l
        hkl(2,i)=m(1,2)*h+m(2,2)*k+m(3,2)*l
        hkl(3,i)=m(1,3)*h+m(2,3)*k+m(3,3)*l
        if (i.eq.1) write (*,10) h,k,l,hkl(1,i),hkl(2,i),hkl(3,i)
      end do
c
   10 format (' First reflection: ',3i6,' => ',3i6)
c
      call jvalut (' Nr of reflections re-indexed :',1,nhkl)
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine rxplor (line,irf,sf,ierr)
c
c ... read Rfree flag from a line of an XPLOR reflection file
c     (also the SIGMa if encountered)
c
      implicit NONE
c
      real sf
c
      integer ll,i,length,ierr,irf,j1,j2,j3,j4,j
c
      character line*(*)
c
code ...
c
c ... initialise
c
      irf = 0
c
      ierr = -1
c
c ... if line too short, return
c
      ll = length(line)
      if (ll .lt. 5) return
c
c ... replace all '=' by spaces
c
      do i=1,ll
        if (line(i:i) .eq. '=') line(i:i) = ' '
      end do
c
c ... replace multiple spaces by a single space and convert to uppercase
c
      call pretty (line)
      call upcase (line)
      ll = length(line)
c
c ... get TEST flag
c
      i = index (line, 'TEST')
      if (i .lt. 1) return
c
      read (line(i+4:),*,err=9900) irf
c
c ... SIGMa optional, so no error if not encountered
c
      ierr = 0
c
      j1 = index (line, 'SIGM')
      if (j1 .lt. 1) return
c
c ... isolate SIGM part
c
      j2 = 0
      j3 = 0
      j4 = 0
      do j=j1+3,ll
        if (line(j:j) .eq. ' ') then
          j3 = j3 + 1
          if (j3 .eq. 1) j2 = j+1
          if (j3 .eq. 2) then
            j4 = j
            goto 300
          end if
        end if
      end do
      if (j3 .eq. 1) then
        j4 = ll
        goto 300
      end if
      return
c
  300 continue
      read (line(j2:j4),*,err=9900) sf
c
      return
c
c ... read from string error
c
 9900 continue
      return
c
      end
c
c
c
      subroutine estuni (cell,res,lat,nasu,nr)
c
c ... calculate nr of unique reflections using the volume of
c     reciprocal space:
c
c            4/3 * PI * (A*B*C)/(D**3)
c     Nref = -------------------------
c                2 * (1 + F) * N
c
c     PI = 3.14...; A/B/C = cell axes (A); D = resolution limit (A)
c     "2" = Friedel mates; F = centering flag (0 for P/R/C/F, 1 for I);
c     N = nr of asymmetric units of the spacegroup
c
c ... rearranging this expression in terms of D gives a formula
c     to estimate the resolution depending on the number of
c     reflections; this is used in the companion subroutine ESTRES
c
      implicit none
c
      integer nasu,nr
c
      real cell(3),res,pi,dc,x
c
      character lat*1
c
code ...
c
      pi = 0.50000 * 6.2831853071796
c
      dc = 0.0
      if (lat .eq. 'I') dc = 1.0
c
      x = 2.0 * pi * cell(1) * cell(2) * cell(3)
      x = x / (3.0 * (1.0 + dc) * res * res * res * float(nasu))
c
      nr = nint (x)
c
      return
      end
c
c
c
      subroutine estres (cell,res,lat,nasu,nr)
c
c ... see comments in subroutine ESTUNI
c
      implicit none
c
      integer nasu,nr
c
      real cell(3),res,pi,dc,x,third
c
      character lat*1
c
code ...
c
      pi = 0.50000 * 6.2831853071796
      third = 1.000D0 / 3.0000D0
c
      dc = 0.0
      if (lat .eq. 'I') dc = 1.0
c
      x = 2.0 * pi * cell(1) * cell(2) * cell(3)
      x = x / (3.0 * (1.0 + dc) * float(nr) * float(nasu))
c
      res = (x)**third
c
      return
      end
c
c
c
      subroutine readxp (iunit,maxhkl,nhkl,hkl,fo,sig,rfl,ierr)
c
      implicit none
c
      integer maxhkl
c
      real fo(maxhkl),sig(maxhkl)
c
      integer hkl(3,maxhkl),rfl(maxhkl),ierr,iunit,nhkl
      integer nh,nf,ns,nr,nl,ll,i,j1,j2,length
c
      character line*256,line2*256,now*4
c
code ...
c
      nh = nhkl
      nf = nhkl
      ns = nhkl
      nr = nhkl
      nl = 0
c
      ierr = -1
c
c ... read header for X-PLOR/CNS files
c
  111 continue
      read (iunit,'(a)',end=200,err=9000) line
      nl = nl + 1
      if (index(line,'INDE') .gt. 0) then
        call prompt (' Found start of reflection list')
        call prompt (
     +    ' Only INDE, FOBS, SIGM and TEST will be read !')
        goto 11
      end if
      call textut (' >',line)
      goto 111
c
   10 continue
      read (iunit,'(a)',end=200,err=9000) line2
      line = line2
      nl = nl + 1
   11 continue
c
      ll = length (line)
      if (ll .lt. 5) goto 19
c
c ... replace all '=' by spaces
c
      do i=1,ll
        if (line(i:i) .eq. '=') line(i:i) = ' '
      end do
c
c ... replace multiple spaces by a single space and convert to uppercase
c
      call pretty (line)
      call upcase (line)
      ll = length(line)
c
c ... remark ?
c
      if (line(1:1) .eq. '!') goto 19
c
c ... check for header lines
c
cc      j1 = index (line,'DECL')
cc      if (j1 .gt. 0) goto 19
c
c ... INDEx ?
c
      now = 'INDE'
      j1 = index (line,now)
      if (j1 .gt. 0) then
c
        if (nh .ge. maxhkl) then
          call errcon ('Too many reflections !!!')
          call jvalut (' Max nr of reflections :',1,maxhkl)
          call prompt (
     +      ' Use ZP_restart command to allocate more memory !')
          ierr = -1
          return
        end if
c
        nh = nh + 1
c
c        if (nh .gt. maxhkl) then
c          nh = maxhkl
c          call errcon ('Too many reflections; rest skipped')
c          goto 200
c        end if
c
        hkl(1,nh) = -9999
        hkl(2,nh) = -9999
        hkl(3,nh) = -9999
        fo (nh) = -9999.0
        sig (nh) = 0.0
        rfl (nh) = 0
c
        j2 = index (line(j1+1:),' ')
        read (line(j2+j1+1:),*,end=9010,err=9010) hkl(1,nh),
     +    hkl(2,nh),hkl(3,nh)
      end if
c
c ... FOBS ?
c
      now = 'FOBS'
      j1 = index (line,now)
      if (j1 .gt. 0) then
        nf = nf + 1
        if (nf .ne. nh) goto 9100
        j2 = index (line(j1+1:),' ')
        read (line(j2+j1+1:),*,end=9010,err=9010) fo(nf)
      end if
c
c ... SIGM ?
c
      now = 'SIGM'
      j1 = index (line,now)
      if (j1 .gt. 0) then
        ns = ns + 1
        if (ns .ne. nh) goto 9100
        j2 = index (line(j1+1:),' ')
        read (line(j2+j1+1:),*,end=9010,err=9010) sig(ns)
      end if
c
c ... TEST ?
c
      now = 'TEST'
      j1 = index (line,now)
      if (j1 .gt. 0) then
        nr = nr + 1
        if (nr .ne. nh) goto 9100
        j2 = index (line(j1+1:),' ')
        read (line(j2+j1+1:),*,end=9010,err=9010) rfl(nr)
      end if
c
c ... done
c
      goto 10
c
c ... echo this line
c
   19 continue
      call textut (' >>>',line)
      goto 10
c
c ... end of file encountered
c
  200 continue
c
      call jvalut (' Nr of lines read       :',1,nl)
      call jvalut (' Nr of reflections read :',1,(nh-nhkl))
c
      if (ns .le. nhkl) then
        call prompt (' WARNING - no SIGMas found')
      end if
c
      if (nr .le. nhkl) then
        call prompt (' No Rfree TEST flags found')
      end if
c
      if (nh .ne. nf .or.
     +    (ns .gt. 0 .and. ns .ne. nh) .or.
     +    (nr .gt. 0 .and. nr .ne. nh) ) then
        call errcon ('Inconsistent file contents')
        call jvalut (' Total INDE read :',1,nh-nhkl)
        call jvalut (' Total FOBS read :',1,nf-nhkl)
        call jvalut (' Total SIGM read :',1,ns-nhkl)
        call jvalut (' Total TEST read :',1,nr-nhkl)
        return
      end if
c
      ierr = 0
      nhkl = nh
c
      return
c
c ... read error from file
c
 9000 continue
      call errcon ('While reading from file')
      call textut (' Previous line :',line)
      return
c
c ... read error from line
c
 9010 continue
      call errcon ('While reading from line')
      call textut (' Attempting to read :',now)
      call textut (' Offending line :',line2)
ccc      print *,ll,j1,j2
      return
c
c ... inconsistency in file
c
 9100 continue
      call errcon ('Inconsistent file contents')
      call jvalut (' Current INDE #:',1,nh-nhkl)
      call jvalut (' Current FOBS #:',1,nf-nhkl)
      call jvalut (' Current SIGM #:',1,ns-nhkl)
      call jvalut (' Current TEST #:',1,nr-nhkl)
      call textut (' Error near line :',line2)
      return
c
      end
c
c
c
      subroutine datain (iunit,file,type,form,maxhkl,nstart,
     +                   hkl,nhkl,fo,sig,rfl,reso,morbit,centri,
     +                   lcell,cell,ierr)
c
c ... DATAIN - read ASCII HKL file
c
      implicit NONE
c
      integer maxopt
      parameter (maxopt=10)
c
      real fo(*),sig(*),reso(*),cell(6),f,sf,rf
c
      integer maxhkl,nhkl,iunit,length,ierr,i1,i2,i3,if,is,irf,i
      integer hkl(3,maxhkl),rfl(*),morbit(*),leng1,ncol,nopt,nfail
      integer nstart
c
      logical xinter,lfree,lshel,lcell,dumpok
c
      character*1 centri(*)
      character file*(*),type*(*),form*(*),line*120
      character optpar(maxopt)*40
c
code ...
c
      ierr = -1
c
      lcell = .false.
c
      ncol = 5
c
      if (type(1:1) .eq. '?') then
        write (*,'(1x,a)')
     +  ' ',
     +  'Supported READ formats:',
     +  '-----------------------',
     +  'SHELXS  -> sets type HKLFS, format (3i4,2f8.2)',
     +  'PROTEIN -> sets type PROTEIN, user format or (*)',
     +  'MKLCF   -> sets type MKLCF, user format or (*)',
     +  'HKLFS   -> sets type HKLFS, user format or (*)',
     +  '2HKLFS  -> sets type HKLFS, 2 per line, user format or (*)',
     +  'RFREE   -> sets type RFREE, user format or (*)',
     +  'ELEANOR -> sets type ELEANOR, user format or (*)',
     +  'XPLOR   -> sets type XPLOR, format (*)',
     +  'X-PLOR  -> sets type XPLOR, format (*)',
     +  'RXPLOR  -> sets type RXPLOR, format (*)',
     +  'RX-PLOR -> sets type RXPLOR, format (*)',
     +  'CNS     -> sets type CNS, format (*)',
     +  'TNT     -> sets type TNT, format (4x,3i4,2f8.1)',
     +  'MTZDUMP -> MTZDUMP output file, user format or (*)',
     +  '*       -> sets type HKLFS, user format or (*)',
     +  ' '
        return
      end if
c
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
      call xopxoa (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening HKL file')
        return
      end if
c
      if (nstart .gt. 0) then
        call prompt (' Appending reflections to an existing dataset')
        call jvalut (' Nr of reflections before append :',1,nstart)
      else
        call prompt (' Reading a new dataset into memory')
      end if
c
      if (length(type) .lt. 1) type = 'HKLFS'
      if (length(form) .lt. 1) form = '*'
c
      if (type(1:1) .eq. '*') type = 'HKLFS'
      call upcase (type)
      call upcase (form)
c
      lshel = .false.
c
      if (type(1:6) .eq. 'SHELXS') then
        type = 'HKLFS'
        form = '(3i4,2f8.2)'
        lshel = .true.
      else if (type(1:3) .eq. 'TNT') then
        type = 'TNT'
        form = '*'
      else if (type(1:5) .eq. 'XPLOR') then
        type = 'XPLOR'
        form = '*'
      else if (type(1:6) .eq. 'X-PLOR') then
        type = 'XPLOR'
        form = '*'
      else if (type(1:3) .eq. 'CNS') then
        type = 'CNS'
        form = '*'
      else if (type(1:6) .eq. 'RXPLOR') then
        type = 'RXPLOR'
        form = '*'
      else if (type(1:7) .eq. 'RX-PLOR') then
        type = 'RXPLOR'
        form = '*'
      else if (type(1:5) .eq. 'RFREE') then
        type = 'RFREE'
      else if (type(1:7) .eq. 'ELEANOR') then
        type = 'ELEANOR'
      else if (type(1:7) .eq. 'MTZDUMP') then
        type = 'MTZDUMP'
      else if (type(1:7) .eq. 'PROTEIN') then
        type = 'PROTEIN'
        call prompt (
     +    ' WARNING - PROTEIN file; SIGFOs will be zero !')
      else if (type(1:5) .eq. 'MKLCF') then
        type = 'MKLCF'
      else
        type = 'HKLFS'
      end if
c
      call textut (' File   :',file)
      call textut (' Type   :',type)
      call textut (' Format :',form)
      lfree = (form(1:1) .eq. '*')
c
      if (type(1:3) .eq. 'TNT') then
        call prompt (' TNT phases and FOMs ignored !')
      end if
c
      if (type(1:7) .eq. 'MTZDUMP') then
        call prompt (' Scanning MTZDUMP file')
  103   continue
        read (iunit,'(a)',end=200,err=9010) line
c
        if (length(line) .le. 0) goto 103
c
        call pretty (line)
        write (*,'(a3,a)') ' > ',line(1:length(line))
c
        if (line(1:6) .eq. '* Cell') then
          read (iunit,'(a)',end=200,err=9010)
          read (iunit,'(a)',end=200,err=9010) line
          call pretty (line)
          write (*,'(a3,a)') ' > ',line(1:leng1(line))
          read (line,*,err=106) (cell(i),i=1,6)
          call fvalut (' ==> Cell :',6,cell)
          lcell = .true.
  106     continue
        end if
c
        if (line(1:21) .eq. '* Number of Columns =') then
          read (line(23:),*) ncol
          call jvalut (' ==> Nr of columns :',1,ncol)
        end if
c
        if (line(1:6) .eq. 'FORMAT') goto 104
        if (line(1:19) .eq. 'LIST OF REFLECTIONS') then
          read (iunit,'(a)',end=200,err=9010)
          read (iunit,'(a)',end=200,err=9010)
          goto 104
        end if
c
        goto 103
c
c LIST OF REFLECTIONS
c123456789012345678901234567890
c
  104   continue
        write (*,*)
        call prompt (' Found start of reflection list')
c
      else
c
c ... make sure you are at the top of the file
c
        rewind (iunit)
      end if
c
      nhkl = nstart
      nfail = 0
c
      if (type(1:5).eq.'XPLOR' .or. type(1:6).eq.'RXPLOR' .or.
     +    type(1:3).eq.'CNS') then
        call readxp (iunit,maxhkl,nhkl,hkl,fo,sig,rfl,ierr)
        if (ierr .ne. 0) goto 9010
        goto 200
      end if
c
c ... read loop
c
  100 continue
c
      irf = 0
c
      if (type(1:7) .eq. 'PROTEIN') then
        if (lfree) then
          read (iunit,*,end=200,err=9010) i1,i2,i3,f
        else
          read (iunit,form,end=200,err=9010) i1,i2,i3,f
        end if
        sf = 0.0
      else if (type(1:5) .eq. 'MKLCF') then
        if (lfree) then
          read (iunit,*,end=200,err=9010) i1,i2,i3,if,is
        else
          read (iunit,form,end=200,err=9010) i1,i2,i3,if,is
        end if
        f = float (if)
        sf = float (is)
      else if (type(1:5) .eq. 'RFREE') then
        if (lfree) then
          read (iunit,*,end=200,err=9010) i1,i2,i3,if,is,irf
        else
          read (iunit,form,end=200,err=9010) i1,i2,i3,if,is,irf
        end if
        f = float(if)
        sf = float (is)
      else if (type(1:7) .eq. 'ELEANOR') then
        if (lfree) then
          read (iunit,*,end=200,err=9010) i1,i2,i3,if,is,rf
        else
          read (iunit,form,end=200,err=9010) i1,i2,i3,if,is,rf
        end if
        f = float(if)
        sf = float (is)
        irf = nint (1.0 - rf)
      else if (type(1:3) .eq. 'TNT') then
        read (iunit,'(a)',end=200,err=9010) line
        call upcase (line)
        if (line(1:4) .ne. 'HKL ') then
          call textut (' Skipped :',line)
          goto 100
        end if
        read (line,'(4x,3i4,2f8.1)',end=9010,err=9010) i1,i2,i3,f,sf
      else if (type(1:7) .eq. 'MTZDUMP') then
        if (lfree) then
ccc          read (iunit,*,end=200,err=200) i1,i2,i3,f,sf
c
          read (iunit,'(a)',end=200,err=200) line
c
          if (line(1:8) .eq. ' MTZDUMP' .or.
     +        line(1:8) .eq. '<B><FONT') then
            write (*,'(a3,a)') ' > ',line(1:leng1(line))
            goto 200
          end if
c
          call extrop (line,nopt,maxopt,optpar,ierr)
          dumpok = (ierr .eq. 0 .and. nopt .ge. 5)
          if (dumpok) then
            call str2i (optpar(1),i1,ierr)
            dumpok = (dumpok .and. (ierr.eq.0))
            call str2i (optpar(2),i2,ierr)
            dumpok = (dumpok .and. (ierr.eq.0))
            call str2i (optpar(3),i3,ierr)
            dumpok = (dumpok .and. (ierr.eq.0))
            if (dumpok) then
              ierr = -1
              if (optpar(4) .ne. '?') call str2r (optpar(4),f,ierr)
              dumpok = (dumpok .and. (ierr.eq.0))
              ierr = -1
              if (optpar(5) .ne. '?') call str2r (optpar(5),sf,ierr)
              dumpok = (dumpok .and. (ierr.eq.0))
            end if
          end if
c
ccc          read (line,*,end=100,err=100) i1,i2,i3,f,sf
c
          if (ncol .ge. 9) then
            do i=1,1+(ncol-8)/6
              read (iunit,*)
            end do
          end if
c
          if (.not. dumpok) then
            if (nfail .lt. 10) then
              call pretty (line)
              call textut (' Garbage (MNFs ???) :',line)
            else if (nfail .eq. 10) then
              call prompt (
     +          ' Additional garbage line(s) suppressed')
            end if
            nfail = nfail + 1
            goto 100
          end if
c
        else
          read (iunit,form,end=200,err=200) i1,i2,i3,f,sf
        end if
      else
        if (lfree) then
          read (iunit,*,end=200,err=9010) i1,i2,i3,f,sf
        else
          read (iunit,form,end=200,err=9010) i1,i2,i3,f,sf
        end if
      end if
c
      if (nhkl .ge. maxhkl) then
        call errcon ('Too many reflections !!!')
        call jvalut (' Max nr of reflections :',1,maxhkl)
        call prompt (
     +    ' Use ZP_restart command to allocate more memory !')
        goto 9010
      end if
c
      nhkl = nhkl + 1
      hkl (1,nhkl) = i1
      hkl (2,nhkl) = i2
      hkl (3,nhkl) = i3
      fo (nhkl)    = f
      sig (nhkl)   = sf
      rfl (nhkl)   = 0
      if (irf .ne. 0) rfl (nhkl) = 1
      reso (nhkl) = 0.0
      morbit (nhkl) = 0
      centri (nhkl) = '?'
c
      goto 100
c
  200 continue
      close (iunit)
c
c ... SHELXS files contain a dummy HKL at the end
c
      if (lshel) nhkl = nhkl - 1
c
      if (type(1:7) .eq. 'MTZDUMP') then
        if (nfail .gt. 0) then
          write (*,*)
          call errcon ('Could not parse all lines (see above) !!!')
          call jvalut (' Ignored (MNFs ???) :',1,nfail)
          write (*,*)
        end if
      end if
c
      if (nstart .gt. 0) then
        call jvalut (' Nr of reflections read  :',1,nhkl-nstart)
        call jvalut (' Nr of reflections total :',1,nhkl)
      else
        call jvalut (' Nr of reflections read :',1,nhkl)
      end if
      ierr = 0
      if (nhkl .lt. 1) then
        call errcon ('Empty data set')
        ierr = -1
      end if
c
      return
c
 9010 continue
      call errcon ('While reading reflection file')
      call jvalut (' Nr of reflections :',1,nhkl)
      ierr = -1
      close (iunit)
c
      return
c
      end
c
c
c
      subroutine dxplor (line,i1,i2,i3,f,sf,ierr)
c
c ... read h,k,l,fobs,sigma from a line of an XPLOR reflection file
c
      implicit NONE
c
      integer i1,i2,i3,ll,i,j,j1,j2,j3,j4,length,ierr
c
      real f,sf
c
      character line*(*)
c
code ...
c
c ... initialise
c
      i1 = 0
      i2 = 0
      i3 = 0
      f  = 0.0
      sf = 0.0
c
      ierr = -1
c
c ... if line too short, return
c
      ll = length(line)
      if (ll .lt. 5) return
c
c ... replace all '=' by spaces
c
      do i=1,ll
        if (line(i:i) .eq. '=') line(i:i) = ' '
      end do
c
c ... replace multiple spaces by a single space and convert to uppercase
c
      call pretty (line)
      call upcase (line)
      ll = length(line)
c
c ... get H, K and L
c
      j1 = index (line, 'INDE')
      if (j1 .lt. 1) return
c
c ... isolate INDEx part
c
      j2 = 0
      j3 = 0
      j4 = 0
      do j=j1+3,ll
        if (line(j:j) .eq. ' ') then
          j3 = j3 + 1
          if (j3 .eq. 1) j2 = j+1
          if (j3 .eq. 4) then
            j4 = j
            goto 100
          end if
        end if
      end do
      if (j3 .eq. 3) then
        j4 = ll
        goto 100
      end if
      return
c
  100 continue
      read (line(j2:j4),*,err=9900) i1,i2,i3
c
      j1 = index (line, 'FOBS')
      if (j1 .lt. 1) return
c
c ... isolate FOBS part
c
      j2 = 0
      j3 = 0
      j4 = 0
      do j=j1+3,ll
        if (line(j:j) .eq. ' ') then
          j3 = j3 + 1
          if (j3 .eq. 1) j2 = j+1
          if (j3 .eq. 2) then
            j4 = j
            goto 200
          end if
        end if
      end do
      if (j3 .eq. 1) then
        j4 = ll
        goto 200
      end if
      return
c
  200 continue
      read (line(j2:j4),*,err=9900) f
c
c ... SIGMa optional, so no error if not encountered
c
      ierr = 0
c
      j1 = index (line, 'SIGM')
      if (j1 .lt. 1) return
c
c ... isolate SIGM part
c
      j2 = 0
      j3 = 0
      j4 = 0
      do j=j1+3,ll
        if (line(j:j) .eq. ' ') then
          j3 = j3 + 1
          if (j3 .eq. 1) j2 = j+1
          if (j3 .eq. 2) then
            j4 = j
            goto 300
          end if
        end if
      end do
      if (j3 .eq. 1) then
        j4 = ll
        goto 300
      end if
      return
c
  300 continue
      read (line(j2:j4),*,err=9900) sf
c
      return
c
c ... read from string error
c
 9900 continue
      return
c
      end
c
c
c
      subroutine datast (nhkl,hkl,fo,sf,re,lr,om,lm,buffer)
c
      implicit none
c
      real fo(*),sf(*),buffer(*),re(*)
      real ave,sdv,xmin,xmax,xtot,var,x1,x2,x3,x4,x5
      real x6,x7,x8,x9,xdum,cc1,cc2,cc3
c
      integer nhkl,i,j,nzeros,nh0,nk0,nl0,nf0,ns0
      integer hkl(3,nhkl),om(*)
c
      logical lr,lm
c
      character txt(3)*2
c
      data txt /'H ','K ','L '/
c
code ...
c
      if (nhkl .lt. 1) then
        call errcon ('No reflections in this dataset')
        return
      end if
c
      nh0 = 0
      nk0 = 0
      nl0 = 0
      nf0 = 0
      ns0 = 0
c
c ... SANITY CHECKS
c
      write (*,*)
      call jvalut (' Total number of reflections   :',1,nhkl)
      if (nhkl .lt. 250) then
        call jvalut (
     +    ' ERROR --- SANITY CHECK - < 250 reflections !',1,nhkl)
      end if
      do i=1,nhkl
        if (hkl(1,i) .eq. 0) nh0 = nh0 + 1
        if (hkl(2,i) .eq. 0) nk0 = nk0 + 1
        if (hkl(3,i) .eq. 0) nl0 = nl0 + 1
        if (fo(i) .lt. 0.01) nf0 = nf0 + 1
        if (sf(i) .lt. 0.01) ns0 = ns0 + 1
      end do
      xdum = 100.0*float(nh0)/float(nhkl)
      call jvalut (' Reflections with H = 0        :',1,nh0)
      if (xdum .ge. 50.0) then
        call fvalut (
     +    ' ERROR --- SANITY CHECK - >= 50 % has H = 0 !',1,xdum)
      end if
      xdum = 100.0*float(nk0)/float(nhkl)
      call jvalut (' Reflections with K = 0        :',1,nk0)
      if (xdum .ge. 50.0) then
        call fvalut (
     +    ' ERROR --- SANITY CHECK - >= 50 % has K = 0 !',1,xdum)
      end if
      xdum = 100.0*float(nl0)/float(nhkl)
      call jvalut (' Reflections with L = 0        :',1,nl0)
      if (xdum .ge. 50.0) then
        call fvalut (
     +    ' ERROR --- SANITY CHECK - >= 50 % has L = 0 !',1,xdum)
      end if
      xdum = 100.0*float(nf0)/float(nhkl)
      call jvalut (' Reflections with Fobs < 0.01  :',1,nf0)
      if (xdum .gt. 1.0) then
        call fvalut (
     +    ' ERROR --- SANITY CHECK - > 1 % has Fobs < 0.01 !',1,xdum)
      end if
      xdum = 100.0*float(ns0)/float(nhkl)
      call jvalut (' Reflections with SigFo < 0.01 :',1,ns0)
      if (xdum .gt. 1.0) then
        call fvalut (
     +    ' WARNING - > 1 % has SigFo < 0.01 !',1,xdum)
      end if
c
      write (*,*)
      write (*,6000)
     +  'Item','Minimum','Maximum','Average','Sdv ','Var ',
     +  '====','=======','=======','=======','=== ','=== '
c
 6000 format (2(1x,a6,2(2x,a10),3(2x,a10),:,/))
 6010 format (1x,a6,2(2x,i10),3(2x,f10.3))
 6020 format (1x,a6,1p,5(2x,e10.3))
 6030 format (1x,a6,5(2x,f10.3))
c
      do i=1,3
        do j=1,nhkl
          buffer (j) = hkl(i,j)
        end do
        call xstats (buffer,nhkl,ave,sdv,xmin,xmax,xtot)
        var = sdv*sdv
        write (*,6010) txt(i),nint(xmin),nint(xmax),ave,sdv,var
        if (sdv .le. 0.1) then
          if (i .eq. 1) then
            call errcon (
     +        'SANITY CHECK - (almost) all H indices equal !')
          else if (i .eq. 2) then
            call errcon (
     +        'SANITY CHECK - (almost) all K indices equal !')
          else
            call errcon (
     +        'SANITY CHECK - (almost) all L indices equal !')
          end if
        end if
      end do
c
      call xstats (fo,nhkl,ave,sdv,xmin,xmax,xtot)
      var = sdv*sdv
      write (*,6020) 'Fobs',xmin,xmax,ave,sdv,var
      if (sdv .le. 0.1) then
        call errcon (
     +    'SANITY CHECK - (almost) all Fobs equal !')
      end if
c
      call xstats (sf,nhkl,ave,sdv,xmin,xmax,xtot)
      var = sdv*sdv
      write (*,6020) 'SigFo',xmin,xmax,ave,sdv,var
      if (sdv .le. 0.1) then
        call prompt (
     +    ' WARNING - (almost) all SigFo equal !')
      end if
c
      if (lr) then
        call xstats (re,nhkl,ave,sdv,xmin,xmax,xtot)
        var = sdv*sdv
        write (*,6030) 'Reso',xmin,xmax,ave,sdv,var
        if (sdv .le. 0.1) then
          call errcon (
     +      'SANITY CHECK - (almost) all Reso equal !')
        end if
      end if
c
      if (lm) then
        do j=1,nhkl
          buffer (j) = om(j)
        end do
        call xstats (buffer,nhkl,ave,sdv,xmin,xmax,xtot)
        var = sdv*sdv
        write (*,6010) 'O.M.',nint(xmin),nint(xmax),ave,sdv,var
        if (sdv .le. 0.1) then
          call prompt (
     +      ' WARNING - (almost) all O.M. equal !')
        end if
      end if
c
      nzeros = 0
      do j=1,nhkl
        if (abs(sf(j)) .ge. 0.00001) then
          buffer (j) = fo(j)/sf(j)
        else
          buffer (j) = 0
          nzeros = nzeros + 1
        end if
      end do
      call xstats (buffer,(nhkl-nzeros),ave,sdv,xmin,xmax,xtot)
      var = sdv*sdv
      write (*,6020) 'Fo/Sig',xmin,xmax,ave,sdv,var
      if (sdv .le. 0.1) then
        call prompt (
     +    ' WARNING - (almost) all Fo/Sig equal !')
      end if
      if (ave .lt. 1.0) then
        call fvalut (
     +    ' ERROR --- SANITY CHECK - average Fo/Sig < 1 !',1,ave)
      else if (ave .lt. 5.0) then
        call fvalut (
     +    ' WARNING - average Fo/Sig < 5 !',1,ave)
      else if (ave .lt. 10.0) then
        call fvalut (
     +    ' NOTE - average Fo/Sig < 10 !',1,ave)
      end if
c
      if (nzeros .gt. 0) then
        call errcon ('There are reflections with SIGMA = 0 !')
        call jvalut (' Number of SIGMA = 0 :',1,nzeros)
      else
        write (*,*)
c
        call xystat (fo,sf,nhkl,x1,x2,x3,x4,x5,x6,x7,x8,x9)
        call fvalut (' Correlation Fobs-SigFo   :',1,x3)
        cc1 = abs(x3)
c
        call xystat (fo,buffer,nhkl,x1,x2,x3,x4,x5,x6,x7,x8,x9)
        call fvalut (' Correlation Fobs-Fo/Sig  :',1,x3)
        cc2 = abs(x3)
c
        call xystat (sf,buffer,nhkl,x1,x2,x3,x4,x5,x6,x7,x8,x9)
        call fvalut (' Correlation SigFo-Fo/Sig :',1,x3)
        cc3 = abs(x3)
c
c ... anything suspicious about the sigmas ?
c
        if (cc1 .gt. 0.9 .and. cc1 .gt. cc2 .and.
     +      cc1 .gt. cc3) then
          call prompt (
     +      ' WARNING - it looks as if Sigma ~ Cst * Fobs !')
        else if (cc2 .gt. 0.9 .and. cc2 .gt. cc1 .and.
     +           cc2 .gt. cc3) then
          call prompt (
     +      ' WARNING - it looks as if Sigma ~ Cst !')
        else if (cc3 .gt. 0.9 .and. cc3 .gt. cc1 .and.
     +           cc3 .gt. cc2) then
          call prompt (
     +      ' WARNING - it looks as if Sigma ~ Cst * SQRT(Fobs) !')
        end if
c
      end if
c
      write (*,*)
      call jvalut (' Nr of reflections      :',1,nhkl)
c
      return
      end
c
c
c
      subroutine dataut (iunit,file,type,form,maxhkl,
     +                   hkl,nhkl,fo,sig,rfl,which,lcell,cell,
     +                   wcent,cent,ierr)
c
c ... DATAUT - write ASCII HKL file
c
      implicit NONE
c
      real fo(*),sig(*),cell(6),f,sf,rf
c
      integer maxhkl,nhkl,iunit,length,ierr,i1,i2,i3,i,irf,isel
      integer hkl(3,maxhkl),rfl(*),nwrit,leng1,icen,mhkl
c
      logical xinter,lshel,lcell
c
      character file*(*),type*(*),form*(*),which*(*),line*120
      character cifstat(0:1)*1,cent(*)*1,wcent*(*)
c
code ...
c
      ierr = -1
c
c ... CIF reflection status (Observed, FreeR)
c
      cifstat (0) = 'o'
      cifstat (1) = 'f'
c
      if (type(1:1) .eq. '?') then
        write (*,'(1x,a)')
     +  ' ',
     +  'Supported WRITE formats:',
     +  '------------------------',
     +  'SHELXS  -> sets type HKLFS, format (3i4,2f8.2)',
     +  'PROTEIN -> sets type PROTEIN, user format or (*)',
     +  'MKLCF   -> sets type MKLCF, user format or (3i6,2i10)',
     +  'HKLFS   -> sets type HKLFS, user format or (*)',
     +  'RFREE   -> sets type HKLFS+Rfree_01, user format or (*)',
     +  'ELEANOR -> sets type HKLFS+(1.0-Rfree_01), user format or (*)',
     +  'XPLOR   -> sets type XPLOR, format (*)',
     +  'X-PLOR  -> sets type XPLOR, format (*)',
     +  'RXPLOR  -> sets type XPLOR+Rfree_01, format (*)',
     +  'RX-PLOR -> sets type XPLOR+Rfree_01, format (*)',
     +  'CNS     -> sets type CNS+Rfree_01, format (*)',
     +  'XCNS    -> sets type CNS, but no Rfree_01, format (*)',
     +  'TNT     -> sets type TNT, format (*)',
     +  'CIF     -> sets type CIF+Rfree_of, format (*)',
     +  'OHKL    -> sets type OHKL, format (a3,3i5,2f15.5)',
     +  '*       -> sets type HKLFS, user format or (*)',
     +  'Default format HKLFS/PROTEIN is (3i6,2f10.3)',
     +  ' '
        return
      end if
c
      if (length(file) .lt. 1) then
        call errcon ('No file name provided')
        return
      end if
c
      call xopxua (iunit,file,xinter(),ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening HKL file')
        return
      end if
c
      if (length(type) .lt. 1) type = 'HKLFS'
      if (length(form) .lt. 1) form = '*'
c
      if (type(1:1) .eq. '*') type = 'HKLFS'
      call upcase (type)
      call upcase (form)
c
      lshel = .false.
c
      if (type(1:6) .eq. 'SHELXS') then
        type = 'SHELXS'
        form = '(3i4,2f8.2)'
        lshel = .true.
      else if (type(1:4) .eq. 'OHKL') then
        type = 'OHKL'
        form = '(''hkl'',3i5,2f15.5)'
      else if (type(1:3) .eq. 'TNT') then
        type = 'TNT'
        form = '(''HKL '',3i4,2f8.1,''  1000.0  0.0000'')'
      else if (type(1:3) .eq. 'CIF') then
        type = 'CIF'
ccc        form = '(3i10,1p,2e15.6,1x,a1)'
        form = '(3i6,2f15.3,1x,a1)'
      else if (type(1:5) .eq. 'XPLOR') then
        type = 'XPLOR'
        form = 
     +    '('' INDEX='',3i6,'' FOBS='',f12.3,'' SIGMA='',f10.3)'
      else if (type(1:6) .eq. 'X-PLOR') then
        type = 'XPLOR'
        form = 
     +    '('' INDEX='',3i6,'' FOBS='',f12.3,'' SIGMA='',f10.3)'
      else if (type(1:6) .eq. 'RXPLOR') then
        type = 'RXPLOR'
        form = 
     +    '('' INDEX='',3i6,'' FOBS='',f12.3,'' SIGMA='',f10.3'//
     +    ','' TEST='',i3)'
      else if (type(1:7) .eq. 'RX-PLOR') then
        type = 'RXPLOR'
        form = 
     +    '('' INDEX='',3i6,'' FOBS='',f12.3,'' SIGMA='',f10.3'//
     +    ','' TEST='',i3)'
      else if (type(1:3) .eq. 'CNS') then
        type = 'CNS'
        form = 
     +    '('' INDE'',3i6,'' FOBS='',f12.3,'' 0.0 SIGMA='',f10.3'//
     +    ','' TEST='',i3)'
      else if (type(1:4) .eq. 'XCNS') then
        type = 'XCNS'
        form = 
     +    '('' INDE'',3i6,'' FOBS='',f12.3,'' 0.0 SIGMA='',f10.3)'
      else if (type(1:7) .eq. 'PROTEIN') then
        type = 'PROTEIN'
        call prompt (
     +    ' WARNING - PROTEIN file; no SIGFOs !')
      else if (type(1:5) .eq. 'MKLCF') then
        type = 'MKLCF'
        if (form(1:1) .eq. '*') form = '(3i6,2i10)'
      else if (type(1:5) .eq. 'RFREE') then
        type = 'RFREE'
        if (form(1:1) .eq. '*') form = '(3i6,2f10.3,i2)'
      else if (type(1:7) .eq. 'ELEANOR') then
        type = 'ELEANOR'
        if (form(1:1) .eq. '*') form = '(3i6,2f10.3,f4.1)'
      else
        type = 'HKLFS'
      end if
      if (form(1:1) .eq. '*') form = '(3i6,2f10.3)'
c
      call textut (' File   :',file)
      call textut (' Type   :',type)
      call textut (' Format :',form)
c
      isel = -1
      if (which(1:1) .eq. 'T') then
        isel = 1
        call prompt (' Write TEST set only')
      else if (which(1:1) .eq. 'W') then
        isel = 0
        call prompt (' Write WORK set only')
      else
        isel = -1
        call prompt (' Write WORK and TEST set')
      end if
c
      icen = -1
      if (wcent(1:1) .eq. 'A') then
        icen = 1
        call prompt (' Write ACENTRICS only')
      else if (wcent(1:1) .eq. 'C') then
        icen = 0
        call prompt (' Write CENTRICS only')
      else
        icen = -1
        call prompt (' Write BOTH centrics and acentrics')
      end if
c
      if (icen .eq. -1 .and. isel .eq. -1) then
        mhkl = nhkl
      else
        mhkl = 0
        do i=1,nhkl
          if (isel .ge. 0) then
            if (rfl(i) .ne. isel) goto 3697
          end if
          if (icen .eq. 1) then
            if (cent(i) .ne. 'A') goto 3697
          else if (icen .eq. 0) then
            if (cent(i) .ne. 'C') goto 3697
          end if
          mhkl = mhkl + 1
 3697     continue
        end do
      end if
      call jvalut (' Nr of reflections to write :',1,mhkl)
c
c ... write header (if appropriate)
c
      if (type(1:3) .eq. 'TNT') then
        call prompt (' TNT phases set to 1000.0, FOMs to 0.0 !')
        call stamp (line)
        write (iunit,'(a3,1x,a)',err=9010) 'REM',line(1:leng1(line))
c
      else if (type(1:3).eq.'CNS') then
          write (iunit,6000) mhkl
c
      else if (type(1:4).eq.'XCNS') then
          write (iunit,6005) mhkl
c
c         write (iunit,6000)
c         if (type(1:6).eq.'RXPLOR') then
c           write (iunit,6010)
c         end if
c
      else if (type(1:4).eq.'OHKL') then
          call stamp (line)
          write (iunit,6089) line(1:leng1(line))
          write (iunit,6080) 'Ohkl'
          if (lcell) then
            write (iunit,6081) (cell(i),i=1,6)
          else
            write (iunit,6089) 'You must edit the cell constants !'
            write (iunit,6081) 1.0,1.0,1.0,90.0,90.0,90.0
            call prompt (' -> You must edit the cell constants !')
          end if
          write (iunit,6089) 'You must edit the spacegroup name !'
          write (iunit,6083) 'P1'
          call prompt (' -> You must edit the spacegroup name !')
          write (iunit,6082) 'Fo','SigF'
c
      else if (type(1:3) .eq. 'CIF') then
        write (iunit,6020) ';'
        write (iunit,6020) file(1:leng1(file))
        call stamp (line)
        write (iunit,6020) line(1:leng1(line))
        write (iunit,6020) ';'
        write (iunit,6020)
        write (iunit,6020) 'data_r0zzzsf'
        if (lcell) then
          write (iunit,6020)
          write (iunit,6030) '_cell.length_a   ',cell(1)
          write (iunit,6030) '_cell.length_b   ',cell(2)
          write (iunit,6030) '_cell.length_c   ',cell(3)
          write (iunit,6030) '_cell.angle_alpha',cell(4)
          write (iunit,6030) '_cell.angle_beta ',cell(5)
          write (iunit,6030) '_cell.angle_gamma',cell(6)
        end if
        write (iunit,6020)
        write (iunit,6020) 'loop_'
        write (iunit,6020) '  _refln.index_h'
        write (iunit,6020) '  _refln.index_k'
        write (iunit,6020) '  _refln.index_l'
        write (iunit,6020) '  _refln.F_meas_au'
        write (iunit,6020) '  _refln.F_meas_sigma_au'
        write (iunit,6020) '  _refln.status'
        write (iunit,6020)
      end if
c
 6000 format (
     +  ' NREFlection= ',i8/
     +  ' ANOMalous=FALSe { equiv. to HERMitian=TRUE}'/
     +  ' DECLare  NAME=FOBS   DOMAin=RECIprocal  TYPE=COMP  END'/
     +  ' DECLare  NAME=SIGMA  DOMAin=RECIprocal  TYPE=REAL  END'/
     +  ' DECLare  NAME=TEST   DOMAin=RECIprocal  TYPE=INTE  END')
c
 6005 format (
     +  ' NREFlection= ',i8/
     +  ' ANOMalous=FALSe { equiv. to HERMitian=TRUE}'/
     +  ' DECLare  NAME=FOBS   DOMAin=RECIprocal  TYPE=COMP  END'/
     +  ' DECLare  NAME=SIGMA  DOMAin=RECIprocal  TYPE=REAL  END')
c
c 6000 format (
c     +  'DECLare  NAME=FOBS    DOMAin=RECIprocal    TYPE=COMP  END'/
c     +  'DECLare  NAME=SIGMA   DOMAin=RECIprocal    TYPE=REAL  END')
c
c 6010 format (
c     +  'DECLare  NAME=TEST    DOMAin=RECIprocal    TYPE=INTE  END')
c
 6020 format (a)
 6030 format (a,1x,f10.2)
c
 6080 format ('Short_name',1x,a)
 6081 format ('Cell',1x,6f10.3)
 6082 format ('SF_field',9(1x,a))
 6083 format ('Space_group',1x,a)
 6089 format ('! ',a)
c
c ... write loop
c
      call prompt (' Writing reflections ...')
      nwrit = 0
      do i=1,nhkl
c
        i1  = hkl(1,i)
        i2  = hkl(2,i)
        i3  = hkl(3,i)
        f   = fo (i)
        sf  = sig (i)
        irf = rfl(i)
c
c ... test/work/both ?
c
        if (isel .ge. 0) then
          if (irf .ne. isel) goto 3627
        end if
c
c ... centrics/acentrics/both ?
c
        if (icen .eq. 1) then
          if (cent(i) .ne. 'A') goto 3627
        else if (icen .eq. 0) then
          if (cent(i) .ne. 'C') goto 3627
        end if
c
        if (type(1:1) .eq. 'P') then
          write (line,form,err=9010) i1,i2,i3,f
        else if (type(1:1) .eq. 'M') then
          write (line,form,err=9010) i1,i2,i3,nint(f),nint(sf)
        else if (type(1:1) .eq. 'X') then
          write (line,form,err=9010) i1,i2,i3,f,sf
          call pretty (line)
        else if (type(1:3) .eq. 'CIF') then
          if (irf .ne. 1) irf = 0
          write (line,form,err=9010) i1,i2,i3,f,sf,cifstat(irf)
ccc          call pretty (line)
        else if (type(1:6) .eq. 'RXPLOR') then
          if (irf .ne. 1) irf = 0
          write (line,form,err=9010) i1,i2,i3,f,sf,irf
          call pretty (line)
        else if (type(1:3) .eq. 'CNS') then
          write (line,form,err=9010) i1,i2,i3,f,sf,irf
          call pretty (line)
        else if (type(1:4) .eq. 'XCNS') then
          write (line,form,err=9010) i1,i2,i3,f,sf
          call pretty (line)
        else if (type(1:5) .eq. 'RFREE') then
          write (line,form,err=9010) i1,i2,i3,f,sf,irf
        else if (type(1:7) .eq. 'ELEANOR') then
          if (irf .ne. 1) irf = 0
          rf = 1.0 - float(irf)
          write (line,form,err=9010) i1,i2,i3,f,sf,rf
        else if (type(1:4) .eq. 'OHKL') then
          write (line,form,err=9010) i1,i2,i3,f,sf
        else
          write (line,form,err=9010) i1,i2,i3,f,sf
        end if
c
        write (iunit,'(a)',err=9010) line(1:leng1(line))
        nwrit = nwrit + 1
c
 3627   continue
c
      end do
c
c ... SHELXS files contain a dummy HKL at the end
c
      if (lshel) write (iunit,form,err=9010) 0,0,0,0.0,0.0
c
c ... finish CIF file
c
      if (type(1:3) .eq. 'CIF') then
        write (iunit,6020)
        write (iunit,6020) ';'
        write (line,*) 'This file should contain ',nwrit,
     +    ' reflections'
        call pretty (line)
        write (iunit,'(a)',err=9010) line(1:leng1(line))
        write (iunit,6020) ';'
        write (iunit,6020)
      end if
c
c ... close the file
c
      close (iunit)
c
      call jvalut (' Nr of reflections stored   :',1,nhkl)
      call jvalut (' Nr of reflections written  :',1,nwrit)
      if (nwrit .ne. mhkl) then
        call errcon ('Gerard cannot count properly')
        call jvalut (' Expected to write          :',1,mhkl)
      end if
      ierr = 0
      return
c
 9010 continue
      call errcon ('While writing reflection file')
      ierr = -1
      return
c
      end
c
c
c
      subroutine imatmu (a, b, c, dim1, dim2, dim3)
c
c ---	Integer_MATrix_MULtiply
c	Multiply 2 matrices to gether and store in a 3rd matrix.
c	The matrices are C(dim1,dim3) = a(dim1,dim2)*b(dim2,dim3)
c	They must be different storeage locations.
c ---	Written by 
c
      implicit NONE
c
      integer dim1, dim2, dim3
      integer a(dim1,dim2), b(dim2, dim3), c(dim1,dim3)
      integer i, j, k
c
code ...
c
      do 100 i=1,dim1
      do 100 j=1,dim3
        c(i,j) = 0
        do 100 k=1,dim2
100       c(i,j) = c(i,j)+ a(i,k) * b(k,j)
      return
      end
c
c
c
      subroutine mulorb (hkl,ns,sym,maxsym,done,om)
c
c ... calculate orbital multiplicity and generate
c     symmetry-equivalent reflections
c
      implicit none
c
      integer ns,hkl(3),i,j,k,res(3),maxsym,om,nd
      integer done(3,maxsym),sym(3,3,ns)
c
code ...
c
      nd = 0
c
      do i=1,ns
        call imatmu (sym(1,1,i),hkl,res,3,3,1)
c
c ... try HKL and its Friedel mate
c
        do k=1,2
c
c ... generate Friedel mate
c
          if (k.eq.2) then
            res(1) = -res(1)
            res(2) = -res(2)
            res(3) = -res(3)
          end if
c
c ... have we seen this one before ?
c
          if (nd .gt. 0) then
            do j=1,nd
              if (res(1) .eq. done(1,j) .and.
     +            res(2) .eq. done(2,j) .and.
     +            res(3) .eq. done(3,j)) goto 100
            end do
          end if
c
c ... we haven't seen this set of HKL before; add it
c
          nd = nd + 1
          done (1,nd) = res (1)
          done (2,nd) = res (2)
          done (3,nd) = res (3)
c
  100     continue
        end do
      end do
c
      om = nd
c
      return
      end
c
c
c
      logical function iscent (hkl,ns,sym)
c
      integer ns,hkl(3),i,res(3)
c
      integer sym(3,3,ns)
c
code ...
c
      do i=1,ns
        call imatmu (sym(1,1,i),hkl,res,3,3,1)
        if ( hkl(1).eq.-res(1) .and.
     +       hkl(2).eq.-res(2) .and.
     +       hkl(3).eq.-res(3) ) then
          iscent = .true.
          return
        end if
      end do
c
      iscent = .false.
c
      return
      end
c
c
c
      subroutine compar (iset,jset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,buffer,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
c
      real x1,x2,x3,x4,x5,x6,x7,x8,x9,top,bot
c
      integer iset,jset,i,new,icode,ni1,ni2,m
c
code ...
c
      call textut (' Comparing Set 1 =',name(iset))
      call textut ('       and Set 2 =',name(jset))
c
      call prompt (' Encoding reflections of set 2 ...')
      do i=1,numhkl(jset)
        call packin (hkl(1,i,jset),hkl(2,i,jset),
     +    hkl(3,i,jset),0,icode)
        ibuff (3*maxhkl+i) = icode
        ibuff (4*maxhkl+i) = i
      end do
c
      call prompt (' Sorting reflections of set 2 ...')
      call shell (ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),numhkl(jset))
c
      ni1 = maxhkl
      ni2 = 2*maxhkl
      new = 0
      call prompt (' Locating reflections of set 1 in set 2 ...')
c
      do i=1,numhkl(iset)
c
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,icode)
c
        call bindex (icode,ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),
     +    numhkl(jset),m)
c
        if (m .le. 0) goto 190
c
c ... found its mate !
c
        new = new + 1
        ni1 = ni1 + 1
        ni2 = ni2 + 1
        buffer (ni1) = fobs (i,iset)
        buffer (ni2) = fobs (m,jset)
c
  190   continue
      end do
c
c ... the old implementation required Order (N^2) operations
c     e.g., for 100,000 reflections, old = 236 CPU seconds,
c     new = 2.3 CPU seconds (SGI R10k) !
c
c      do i=1,numhkl(jset)
c        call packin (hkl(1,i,jset),hkl(2,i,jset),
c     +    hkl(3,i,jset),0,jcode)
c        icode = iindex (jcode,0,numhkl(iset),ibuff)
c        if (icode .gt. 0) then
c          new = new + 1
c          ni1 = ni1 + 1
c          ni2 = ni2 + 1
c          buffer (ni1) = fobs (icode,iset)
c          buffer (ni2) = fobs (i,jset)
c        end if
c      end do
c
      call jvalut (' HKLs in set 1 :',1,numhkl(iset))
      call jvalut (' HKLs in set 2 :',1,numhkl(jset))
      call jvalut (' HKLs in both  :',1,new)
c
      if (new .le. 0) return
c
      if (new .lt. nint(0.1*numhkl(iset)) .or.
     +    new .lt. nint(0.1*numhkl(jset)) ) then
        call prompt (' WARNING - few reflections in common !')
      end if
c
      call xystat (buffer(maxhkl+1),buffer(2*maxhkl+1),new,
     +             x1,x2,x3,x4,x5,x6,x7,x8,x9)
      call fvalut (' Correlation coeff Fobs :',1,x3)
      call fvalut (' Shape similarity  Fobs :',1,x2)
      call rvalut (' RMS difference Fo1/Fo2 :',1,x1)
      call rvalut (' R=SUM(Fo1-Fo2)/SUM(Fo1):',1,x4)
      call rvalut (' R with (Fo1-S*Fo2)     :',1,x6)
      call rvalut ('          where scale S :',1,x8)
      call rvalut (' R=SUM(Fo1-Fo2)/SUM(Fo2):',1,x5)
      call rvalut (' R with (S*Fo1-Fo2)     :',1,x7)
      call rvalut ('          where scale S :',1,x9)
c
      top = 0.0
      bot = 0.0
c
      do i=1,new
        x1 = buffer(maxhkl+i)
        x2 = x8*buffer(2*maxhkl+i)
        top = top + abs(x1 - x2)
        bot = bot + abs(x1 + x2)
      end do
c
      call prompt (' Rmerge = SUM |F1-S*F2| / SUM |F1+S*F2|')
      call fvalut (' Value of Rmerge :',1,(top/bot))
c
      return
      end
c
c
c
      subroutine yeates (iset,jset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,centri,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      character*1 centri(maxhkl,maxset)
c
      real sumch,sumch2,sumah,sumah2,hhh,ppp,qqq,zeroto
c
      integer iset,jset,i,new,icode,m,na,nc
c
      data zeroto /1.0E-6/
c
code ...
c
      call textut (' Comparing Set 1 =',name(iset))
      call textut ('       and Set 2 =',name(jset))
c
      call prompt (' Encoding reflections of set 2 ...')
      do i=1,numhkl(jset)
        call packin (hkl(1,i,jset),hkl(2,i,jset),
     +    hkl(3,i,jset),0,icode)
        ibuff (3*maxhkl+i) = icode
        ibuff (4*maxhkl+i) = i
      end do
c
      call prompt (' Sorting reflections of set 2 ...')
      call shell (ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),numhkl(jset))
c
      sumch = 0.0
      sumch2 = 0.0
      sumah = 0.0
      sumah2 = 0.0
      nc = 0
      na = 0
      new = 0
c
      call prompt (' Locating reflections of set 1 in set 2 ...')
c
      do i=1,numhkl(iset)
c
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,icode)
c
        call bindex (icode,ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),
     +    numhkl(jset),m)
c
        if (m .le. 0) goto 190
c
c ... found its mate !
c
        new = new + 1
c
        ppp = fobs(i,iset)*fobs(i,iset)
        qqq = fobs(m,jset)*fobs(m,jset)
        if ( (ppp+qqq) .le. zeroto) goto 190
        hhh = abs( (ppp-qqq) / (ppp+qqq) )
c
        if (centri(i,iset) .eq. 'C') then
          nc = nc + 1
          sumch = sumch + hhh
          sumch2 = sumch2 + hhh*hhh
        else
          na = na + 1
          sumah = sumah + hhh
          sumah2 = sumah2 + hhh*hhh
        end if
c
  190   continue
      end do
c
      call jvalut (' HKLs in set 1 :',1,numhkl(iset))
      call jvalut (' HKLs in set 2 :',1,numhkl(jset))
      call jvalut (' HKLs in both  :',1,new)
c
      if (new .le. 0) return
c
      if (new .lt. nint(0.1*numhkl(iset)) .or.
     +    new .lt. nint(0.1*numhkl(jset)) ) then
        call prompt (' WARNING - few reflections in common !')
      end if
c
 6000 format (1x,a,' = ',f6.3,'  Expected identical = ',f6.3,
     +  '  Expected unrelated = ',f6.3)
c
ccc      print *,sumch,sumch2,sumah,sumah2
c
      call jvalut (' Common centrics  :',1,nc)
      sumch = sumch / float(nc)
      sumch2 = sumch2 / float(nc)
      write (*,6000) '<|H|>',sumch,0.0,(2.0/pi)
      write (*,6000) '<H^2>',sumch2,0.0,0.5
c
      call jvalut (' Common acentrics :',1,na)
      sumah = sumah / float(na)
      sumah2 = sumah2 / float(na)
      write (*,6000) '<|H|>',sumah,0.0,0.5
      write (*,6000) '<H^2>',sumah2,0.0,(1.0/3.0)
c
      return
      end
c
c
c
      subroutine kilrog (iset,ih,ik,il,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
c
c ... delete rogue reflections by HKL
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer morbit(maxhkl,maxset),rfree(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      integer iset,i,j,ierr,ih,ik,il,hklbad,ibad,iindex
c
      character line*80
c
code ...
c
      ierr = -1
c
      write (line,*) ih,ik,il
      call pretty (line)
      call textut (' Rogue :',line)
c
      call packin (ih,ik,il,0,hklbad)
      do i=1,numhkl(iset)
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,ibuff(i))
      end do
c
c ... note: we only look up one reflection, so there's no
c     point in using the bisection technique since this
c     requires encoding and sorting all reflections first !
c
      ibad = iindex (hklbad,0,numhkl(iset),ibuff)
c
      if (ibad .le. 0) then
        call errcon ('Rogue hkl not found')
        return
      end if
c
      write (*,6000) ibad,(hkl(j,ibad,iset),j=1,3),
     +      fobs(ibad,iset),sigfob(ibad,iset)
c
 6000 format (' # ',i8,'  HKL ',3i6,'  Fobs & SigFob = ',1p,2e12.4)
c
      do i=ibad,numhkl(iset)-1
        j = i+1
        hkl(1,i,iset) = hkl(1,j,iset)
        hkl(2,i,iset) = hkl(2,j,iset)
        hkl(3,i,iset) = hkl(3,j,iset)
        fobs(i,iset)   = fobs(j,iset)
        sigfob(i,iset) = sigfob(j,iset)
        reso(i,iset)   = reso(j,iset)
        centri(i,iset) = centri(j,iset)
        morbit(i,iset) = morbit(j,iset)
        rfree(i,iset)  = rfree(j,iset)
      end do
      numhkl (iset) = numhkl (iset) - 1
c
      call jvalut (' Nr of reflections now :',1,numhkl(iset))
      ierr = 0
c
      return
      end
c
c
c
      subroutine oddevn (mode,iset,crit,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri,ibuff)
c
c ... delete reflexions with odd or even H or K or L
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer morbit(maxhkl,maxset),rfree(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      integer iset,i,new,ierr,j,k,ikill,nold
c
      logical lodd
c
      character crit*(*),line*80,mode*(*)
c
code ...
c
      ierr = -1
c
      lodd = (mode(1:1) .eq. 'O')
c
      ikill = 0
      if (crit(1:1) .eq. 'H') ikill = 1
      if (crit(1:1) .eq. 'K') ikill = 2
      if (crit(1:1) .eq. 'L') ikill = 3
      if (ikill .eq. 0) then
        call errcon ('Not H, K or L selected')
        return
      end if
c
      nold = numhkl(iset)
      call jvalut (' Nr of reflections before :',1,nold)
c
      if (lodd) then
        write (line,*) crit(1:1),' ODD'
      else
        write (line,*) crit(1:1),' EVEN'
      end if
      call pretty (line)
      call textut (' Kill reflection if :',line)
c
      new = 0
c
      if (lodd) then
c
        do i=1,numhkl(iset)
c
          j = hkl(ikill,i,iset)
          k = 2 * (j/2)
c
c ... EVEN; keep this one
c
          if (j .eq. k) then
            new = new + 1
            ibuff (new) = i
          end if
c
        end do
c
      else
c
        do i=1,numhkl(iset)
c
          j = hkl(ikill,i,iset)
          k = 2 * (j/2)
c
c ... ODD; keep this one
c
          if (j .ne. k) then
            new = new + 1
            ibuff (new) = i
          end if
c
        end do
c
      end if
c
      if (new .eq. 0) then
        call errcon ('This would kill all reflections; aborted')
        return
      end if
c
      if (new .eq. numhkl(iset)) then
        call prompt (' NOTE --- No reflections killed')
        return
      end if
c
      do i=1,new
        j=ibuff(i)
        if (j .ne. i) then
          hkl(1,i,iset) = hkl(1,j,iset)
          hkl(2,i,iset) = hkl(2,j,iset)
          hkl(3,i,iset) = hkl(3,j,iset)
          fobs(i,iset)   = fobs(j,iset)
          sigfob(i,iset) = sigfob(j,iset)
          reso(i,iset)   = reso(j,iset)
          centri(i,iset) = centri(j,iset)
          morbit(i,iset) = morbit(j,iset)
          rfree(i,iset)  = rfree(j,iset)
        end if
      end do
c
      numhkl (iset) = new
      call jvalut (' Nr of reflections after :',1,new)
      ierr = 0
c
      i = nold - new
      if (i .gt. 0) then
        call jvalut (' WARNING - Reflections killed :',1,i)
      end if
c
      return
      end
c
c
c
      subroutine hklspe (iset,condit,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,
     +                   hkl,morbit,rfree,
     +                   centri)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
      integer morbit(maxhkl,maxset),rfree(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      real xdum
c
      integer iset,i,j,n,which(3)
c
      logical iszero(3),other(3)
c
      character condit*3
c
code ...
c
      do i=1,3
        iszero (i) = (condit(i:i).eq.'0')
      end do
      other (1) = (condit(1:1) .eq. 'K' .or.
     +             condit(1:1) .eq. 'L')
      which (1) = 1
      if (condit(1:1) .eq. 'K') which(1)=2
      if (condit(1:1) .eq. 'L') which(1)=3
c
      other (2) = (condit(2:2) .eq. 'H' .or.
     +             condit(2:2) .eq. 'L')
      which (2) = 2
      if (condit(2:2) .eq. 'H') which(2)=1
      if (condit(2:2) .eq. 'L') which(2)=3
c
      other (3) = (condit(3:3) .eq. 'H' .or.
     +             condit(3:3) .eq. 'K')
      which (3) = 3
      if (condit(3:3) .eq. 'H') which(3)=1
      if (condit(3:3) .eq. 'K') which(3)=2
c
      write (*,6901) 'Reflxn','   H','   K','   L',
     +    'Resoln','   Fobs   ',' Sigma(Fobs)',' F/Sig',
     +    'Test','A/C','Orb'
c
      n = 0
      do i=1,numhkl(iset)
        if (iszero(1)) then
          if (hkl(1,i,iset) .ne. 0) goto 10
        end if
        if (iszero(2)) then
          if (hkl(2,i,iset) .ne. 0) goto 10
        end if
        if (iszero(3)) then
          if (hkl(3,i,iset) .ne. 0) goto 10
        end if
c
c ... check HHL etc.
c
        if (other(1)) then
          if (hkl(1,i,iset) .ne. hkl(which(1),i,iset)) goto 10
        end if
c
        if (other(2)) then
          if (hkl(2,i,iset) .ne. hkl(which(2),i,iset)) goto 10
        end if
c
        if (other(3)) then
          if (hkl(3,i,iset) .ne. hkl(which(3),i,iset)) goto 10
        end if
c
c ... write this reflection
c
        n = n + 1
        if (sigfob(i,iset) .ne. 0.0) then
          xdum = fobs(i,iset) / sigfob(i,iset)
        else
          xdum = -1.0
        end if
        write (*,6900) i,(hkl(j,i,iset),j=1,3),reso(i,iset),
     +      fobs(i,iset),sigfob(i,iset),xdum,rfree(i,iset),
     +      centri(i,iset),morbit(i,iset)
c
   10   continue
      end do
c
      call jvalut (' Reflections listed :',1,n)
c
      return
c
 6900 format (1x,i6,1x,3i4,1x,f6.2,1x,1p,2e12.4,1x,0p,f6.2,
     +  1x,i4,3x,a1,1x,i3)
 6901 format (1x,a6,1x,3a4,1x,a6,1x,2a12,1x,a6,1x,a4,1x,a3,1x,a3)
c
      end
c
c
c
      subroutine testrf (iset,maxset,maxhkl,rfree)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      integer rfree(maxhkl,maxset)
c
      integer ntest,nwork,i,iset
c
      real perc
c
code ...
c
      ntest = 0
      nwork = 0
      do i = 1, numhkl (iset)
        if (rfree(i,iset) .eq. 1) then
          ntest = ntest + 1
        else
          nwork = nwork + 1
        end if
      end do
c
      call prompt (' TEST has FLAG=1; WORK has FLAG<>1')
      call jvalut (' Nr of WORK reflections :',1,nwork)
      call jvalut (' Nr of TEST reflections :',1,ntest)
c
      if (numhkl(iset) .gt. 0) then
        perc = 100.0 * float(ntest) / float(numhkl(iset))
      else
        perc = 0.0
      end if
      call fvalut (' Percentage TEST data   :',1,perc)
c
      if (ntest .le. 0) then
        know (kfree,iset) = .false.
        call prompt (' This is NOT a cross-validation dataset')
      else
        know (kfree,iset) = .true.
        call prompt (' This is a cross-validation dataset')
      end if
c
      if (perc .gt. 15.0 .and. perc .le. 25.0) then
        call prompt (' WARNING - more than 15% TEST reflections !')
      else if (perc .gt. 25.0) then
        call prompt (' ERROR - more than 25% TEST reflections !!!')
      end if
c
      if (ntest .lt. 500) then
        call prompt (' WARNING - fewer than 500 TEST reflections !')
      else if (ntest .gt. 2500) then
        call prompt (' WARNING - more than 2500 TEST reflections !')
      end if
c
      return
      end
c
c
c
      subroutine rfadju (iset,maxset,maxhkl,rfree,newper)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      integer rfree(maxhkl,maxset)
c
      integer ntest,nwork,i,iset,iold,inew,nflip
c
      real perc,newper,dx,xcut
c
code ...
c
      ntest = 0
      nwork = 0
      do i = 1, numhkl (iset)
        if (rfree(i,iset) .eq. 0) then
          nwork = nwork + 1
        else
          ntest = ntest + 1
        end if
      end do
c
      call jvalut (' Nr of WORK reflections :',1,nwork)
      call jvalut (' Nr of TEST reflections :',1,ntest)
      if (numhkl(iset) .gt. 0) then
        perc = 100.0 * float(ntest) / float(numhkl(iset))
      else
        perc = 0.0
      end if
      call fvalut (' Percentage TEST data   :',1,perc)
      call fvalut (' Requested percentage   :',1,newper)
c
      if (newper .eq. perc) then
        call errcon ('New percentage equal to current !')
        return
      end if
c
      if (newper .lt. perc) then
        iold = ntest
        inew = nint (0.01 * newper * float(numhkl(iset)))
        xcut = float(inew) / float(iold)
        nflip = 0
        do i = 1, numhkl (iset)
          if (rfree(i,iset) .eq. 1) then
            call gkrand (dx,0.0,1.0,0)
            if (dx .ge. xcut) then
              rfree(i,iset) = 0
              nflip = nflip + 1
            end if
          end if
        end do
        ntest = ntest - nflip
        perc = 100.0 * float(ntest) / float(numhkl(iset))
        call fvalut (' Actual new percentage  :',1,perc)
      else if (newper .gt. perc) then
        iold = ntest
        inew = nint (0.01 * newper * float(numhkl(iset)))
        xcut = float(inew-iold) / float(numhkl(iset))
        nflip = 0
        do i = 1, numhkl (iset)
          if (rfree(i,iset) .eq. 0) then
            call gkrand (dx,0.0,1.0,0)
            if (dx .lt. xcut) then
              rfree(i,iset) = 1
              nflip = nflip + 1
            end if
          end if
        end do
        ntest = ntest + nflip
        perc = 100.0 * float(ntest) / float(numhkl(iset))
        call fvalut (' Actual new percentage  :',1,perc)
      end if
c
      return
      end
c
c
c
      subroutine scattr (iset,pfile,hori,vert,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,rfree,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer rfree(maxhkl,maxset)
c
      integer iset,ierr,i,nd,flag,iunit,leng1
c
      real x,y,z,ave,sdv,xmin,xmax,xdum,ymin,ymax
      real pxmin,pxmax,pymin,pymax
c
      logical xinter
c
      character pfile*(*),hori*(*),vert*(*),line*80
c
code ...
c
      ierr = 0
      flag = 0
      iunit = -1
c
      call textut (' Scatter plot Set =',name(iset))
c
   10 continue
c
c ... collect X and Y values
c
      nd = 0
c
      do i=1,numhkl(iset)
        x = 0.0
        y = 0.0
        if (rfree(i,iset) .eq. flag) then
          nd = nd + 1
          if (hori(1:3) .eq. 'FOB') then
            x = fobs (i,iset)
          else if (hori(1:3) .eq. 'SIG') then
            x = sigfob (i,iset)
          else if (hori(1:3) .eq. 'F/S') then
            x = fobs (i,iset) / sigfob (i,iset)
          else if (hori(1:3) .eq. 'INT') then
            x = fobs (i,iset) ** 2
          else if (hori(1:3) .eq. 'LNI') then
            if (fobs(i,iset).gt.0.0) x = alog(fobs(i,iset)**2)
          else if (hori(1:3) .eq. 'I/S') then
            x = (fobs (i,iset) ** 2) / (sigfob (i,iset) ** 2)
          else if (hori(1:3) .eq. 'RES') then
            x = reso (i,iset)
          else if (hori(1:3) .eq. '1/R') then
            x = 1.0 / reso (i,iset)
          else if (hori(1:3) .eq. 'STL') then
            x = 0.5 / reso (i,iset)
          else if (hori(1:3) .eq. 'DST') then
            z = 0.5 / reso (i,iset)
            x = 4.0 * z * z
          else
            call errcon ('Invalid hori variable')
            ierr = -1
            goto 9000
          end if
          if (vert(1:3) .eq. 'FOB') then
            y = fobs (i,iset)
          else if (vert(1:3) .eq. 'SIG') then
            y = sigfob (i,iset)
          else if (vert(1:3) .eq. 'F/S') then
            y = fobs (i,iset) / sigfob (i,iset)
          else if (vert(1:3) .eq. 'INT') then
            y = fobs (i,iset) ** 2
          else if (vert(1:3) .eq. 'LNI') then
            if (fobs(i,iset).gt.0.0) y = alog(fobs(i,iset)**2)
          else if (vert(1:3) .eq. 'I/S') then
            y = (fobs (i,iset) ** 2) / (sigfob (i,iset) ** 2)
          else if (vert(1:3) .eq. 'RES') then
            y = reso (i,iset)
          else if (vert(1:3) .eq. '1/R') then
            y = 1.0 / reso (i,iset)
          else if (vert(1:3) .eq. 'STL') then
            y = 0.5 / reso (i,iset)
          else if (vert(1:3) .eq. 'DST') then
            z = 0.5 / reso (i,iset)
            y = 4.0 * z * z
          else
            call errcon ('Invalid vert variable')
            ierr = -1
            return
          end if
          buffer (nd) = x
          buffer (maxhkl+nd) = y
        end if
      end do
c
c ... print some info
c
      write (*,*)
      call ivalut (' Rfree flag  :',1,flag)
      call jvalut (' Data points :',1,nd)
c
      if (nd .lt. 1) then
        call errcon ('No data points')
        ierr = -1
        return
      end if
c
      write (*,'(a,a3,a,a3)') ' Plot ',vert,' versus ',hori
      call xstats (buffer(1),nd,ave,sdv,xmin,xmax,xdum)
      write (*,6000) hori,'MIN',xmin,'MAX',xmax,'AVE',ave,'SDV',sdv
      call xstats (buffer(maxhkl+1),nd,ave,sdv,ymin,ymax,xdum)
      write (*,6000) vert,'MIN',ymin,'MAX',ymax,'AVE',ave,'SDV',sdv
c
c ... write plot file
c
      if (iunit .le. 0) then
        iunit = 45
        call xopxua (iunit,pfile,xinter(),ierr)
        if (ierr .ne. 0) goto 9000
c
        z = 0.025 * (xmax-xmin)
        pxmin = xmin - z
        pxmax = xmax + z
        z = 0.025 * (ymax-ymin)
        pymin = ymin - z
        pymax = ymax + z
c
c ... write header
c
        call stamp (line)
        write (iunit,6100) 'REMARK ',line(1:leng1(line))
        write (iunit,6100) 'REMARK DATAMAN scatter plot'
        write (iunit,6100) 'REMARK Filename = ',
     +    pfile(1:leng1(pfile))
        write (iunit,6100) 'REMARK Dataset = ',
     +    name(iset)(1:leng1(name(iset))),' File = ',
     +    file(iset)(1:leng1(file(iset)))
        write (iunit,6100) 'REMARK Comment = ',
     +    coment(iset)(1:leng1(coment(iset)))
c
        write (iunit,6100) 'XLABEL ',hori(1:leng1(hori))
        write (iunit,6100) 'YLABEL ',vert(1:leng1(vert))
        write (iunit,6100) 'LINFIT'
      end if
c
      write (iunit,6100) '!'
      write (iunit,6110) 'NPOINT ',nd
      write (iunit,6110) 'MRKTYP ',(flag+1)
      write (iunit,6110) 'COLOUR ',(4-3*flag)
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (buffer(i),i=1,nd)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (buffer(i),i=maxhkl+1,maxhkl+nd)
      write (iunit,6100) '!'
c
      if (flag .eq. 0 .and. know(kfree,iset)) then
        flag = 1
        write (iunit,6100) 'MORE'
        write (iunit,6100) 'REMARK Work and Test data separated'
        goto 10
      end if
c
      z = 0.025 * (xmax-xmin)
      pxmin = min (pxmin, xmin - z)
      pxmax = max (pxmax, xmax + z)
      z = 0.025 * (ymax-ymin)
      pymin = min (pymin, ymin - z)
      pymax = max (pymax, ymax + z)
c
      write (iunit,6120) 'XYVIEW ',pxmin,pxmax,pymin,pymax
      write (iunit,6100) '!'
      write (iunit,6100) 'END'
c
      call prompt (' Plot file generated')
c
c 'FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|LNI'
c
 9000 continue
      if (iunit .gt. 0) close (iunit)
c
      return
c
 6000 format (1x,a3,1p,4(1x,a3,1x,e12.4))
 6100 format (10a)
 6110 format (a,i10)
 6120 format (a,1p,4e12.4)
c
      end
c
c
c
      subroutine binplt (iset,pfile,hori,vert,binsiz,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),rfree(maxhkl,maxset)
c
      integer maxbin, minbin
      parameter (maxbin = 100, minbin = 5)
c
      integer iset,ierr,i,iunit,nbins,j,leng1
      integer ninbin(maxbin),rinbin(maxbin)
c
      real x,y,z,xmin,xmax,ymin,ymax
      real pxmin,pxmax,pymin,pymax,binsiz
      real xval(maxbin),yval(maxbin),rval(maxbin)
      real rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2
c
      logical xinter
c
      character pfile*(*),hori*(*),vert*(*),line*80
      character esv1*12,esv2*12
c
code ...
c
      ierr = 0
      iunit = -1
c
      call textut (' Bin plot Set =',name(iset))
c
c ... get X range & collect data
c
      do i=1,numhkl(iset)
        buffer (i) = 0.0
        buffer (maxhkl+i) = 0.0
      end do
c
      do i=1,numhkl(iset)
        x = 0.0
        y = 0.0
        if (hori(1:3) .eq. 'FOB') then
          x = fobs (i,iset)
        else if (hori(1:3) .eq. 'SIG') then
          x = sigfob (i,iset)
        else if (hori(1:3) .eq. 'F/S') then
          x = fobs (i,iset) / sigfob (i,iset)
        else if (hori(1:3) .eq. 'INT') then
          x = fobs (i,iset) ** 2
        else if (hori(1:3) .eq. 'LNI') then
          if (fobs(i,iset).gt.0.0) x = alog(fobs(i,iset)**2)
        else if (hori(1:3) .eq. 'I/S') then
          x = (fobs (i,iset) ** 2) / (sigfob (i,iset) ** 2)
        else if (hori(1:3) .eq. 'RES') then
          x = reso (i,iset)
        else if (hori(1:3) .eq. '1/R') then
          x = 1.0 / reso (i,iset)
        else if (hori(1:3) .eq. 'STL') then
          x = 0.5 / reso (i,iset)
        else if (hori(1:3) .eq. 'DST') then
          z = 0.5 / reso (i,iset)
          x = 4.0 * z * z
        else if (hori(1:3) .eq. 'H  ') then
          x = float (hkl(1,i,iset))
        else if (hori(1:3) .eq. 'K  ') then
          x = float (hkl(2,i,iset))
        else if (hori(1:3) .eq. 'L  ') then
          x = float (hkl(3,i,iset))
        else if (hori(1:3) .eq. 'H/A') then
          x = float (hkl(1,i,iset)) / cell (1,iset)
        else if (hori(1:3) .eq. 'K/B') then
          x = float (hkl(2,i,iset)) / cell (2,iset)
        else if (hori(1:3) .eq. 'L/C') then
          x = float (hkl(3,i,iset)) / cell (3,iset)
        else
          call errcon ('Invalid hori variable')
          ierr = -1
          goto 9000
        end if
c
        if (vert(1:3) .eq. 'FOB') then
          y = fobs (i,iset)
        else if (vert(1:3) .eq. 'SIG') then
          y = sigfob (i,iset)
        else if (vert(1:3) .eq. 'F/S') then
          y = fobs (i,iset) / sigfob (i,iset)
        else if (vert(1:3) .eq. 'INT') then
          y = fobs (i,iset) ** 2
        else if (vert(1:3) .eq. 'LNI') then
          if (fobs(i,iset).gt.0.0) y = alog(fobs(i,iset)**2)
        else if (vert(1:3) .eq. 'I/S') then
          y = (fobs (i,iset) ** 2) / (sigfob (i,iset) ** 2)
        else if (vert(1:3) .eq. 'RES') then
          y = reso (i,iset)
        else if (vert(1:3) .eq. '1/R') then
          y = 1.0 / reso (i,iset)
        else if (vert(1:3) .eq. 'STL') then
          y = 0.5 / reso (i,iset)
        else if (vert(1:3) .eq. 'DST') then
          z = 0.5 / reso (i,iset)
          y = 4.0 * z * z
        else if (vert(1:3) .eq. 'NRF') then
          y = 1.0
        else
          call errcon ('Invalid vert variable')
          ierr = -1
          return
        end if
c
        buffer (i) = x
        buffer (maxhkl+i) = y
c
        if (i .eq. 1) then
          xmin = x
          xmax = x
        else
          xmin = min (xmin, x)
          xmax = max (xmax, x)
        end if
      end do
c
      write (*,'(1x,a3,1x,a,1x,1p,2e12.4)') hori,
     +  'Min, Max =',xmin,xmax
c
      write (*,'(1x,a,2i3)') 'Min and Max nr of bins = ',
     +  minbin,maxbin
c
      if (binsiz .gt. 0.0) then
        nbins = (xmax-xmin)/binsiz
      else
        nbins = max (1, nint(-binsiz))
      end if
      nbins = max (minbin, min(nbins,maxbin))
      call ivalut (' Nr of bins :',1,nbins)
      binsiz = (xmax-xmin)/float(nbins)
      call rvalut (' Bin size   :',1,binsiz)
c
      do i=1,nbins
        xval (i) = xmin + float(i-1)*binsiz + 0.5*binsiz
        yval (i) = 0.0
        rval (i) = 0.0
        ninbin (i) = 0
        rinbin (i) = 0
      end do
c
c ... collect data in bins
c
      do i=1,numhkl(iset)
        x = buffer (i)
        y = buffer (maxhkl+i)
        j = max (1, min (1 + int((x-xmin)/binsiz), nbins))
        if (rfree(i,iset) .eq. 0) then
          ninbin (j) = ninbin (j) + 1
          yval (j) = yval (j) + y
        else
          rinbin (j) = rinbin (j) + 1
          rval (j) = rval (j) + y
        end if
      end do
c
      if (vert(1:3) .ne. 'NRF') then
        do i=1,nbins
          if (ninbin(i) .gt. 0) then
            yval (i) = yval (i) / float(ninbin(i))
          end if
          if (rinbin(i) .gt. 0) then
            rval (i) = rval (i) / float(rinbin(i))
          end if
        end do
      end if
c
c ... print table
c
      write (*,'(a,a3,a,a3)') ' Plot ',vert,' versus ',hori
      if (know(kfree,iset)) then
        esv1 = ' <'//vert(1:3)//'> Work'
        esv2 = ' <'//vert(1:3)//'> Test'
        write (*,6300) 'Bin nr','Start value',esv1,
     +    'Nr values',esv2,'Nr values'
      else
        esv1 = ' <'//vert(1:3)//'>'
        write (*,6300) 'Bin nr','Start value',esv1,
     +    'Nr values'
      end if
c
      do i=1,nbins
        if (know(kfree,iset)) then
          write (*,6310) i,xval(i),yval(i),ninbin(i),
     +      rval(i),rinbin(i)
        else
          write (*,6310) i,xval(i),yval(i),ninbin(i)
        end if
      end do
c
 6300 format (1x,6a12)
 6310 format (1x,i12,1p,2e12.4,i12,e12.4,i12)
c
      if (know(kfree,iset)) then
        call xystat (yval,rval,nbins,
     +    rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2)
        write (*,*)
        call prompt (' Comparison for WORK and TEST data :')
        call fvalut (' Correlation coefficient :',1,corr)
        call rvalut (' Scaled R w.r.t. WORK    :',1,rf3)
        call rvalut (' Scaled R w.r.t. TEST    :',1,rf4)
        call rvalut (' RMS difference          :',1,rmsd)
      end if
c
c ... write plot file
c
      iunit = 45
      call xopxua (iunit,pfile,xinter(),ierr)
      if (ierr .ne. 0) goto 9000
c
      z = 0.025 * (xmax-xmin)
      pxmin = xmin - z
      pxmax = xmax + z
      ymin = yval(1)
      ymax = yval(1)
      do i=1,nbins
        ymin = min (ymin,yval(i))
        ymax = max (ymax,yval(i))
        if (know(kfree,iset)) then
          ymin = min (ymin,rval(i))
          ymax = max (ymax,rval(i))
        end if
      end do
      z = 0.025 * (ymax-ymin)
      pymin = ymin - z
      pymax = ymax + z
c
c ... write header
c
      call stamp (line)
      write (iunit,6100) 'REMARK ',line(1:leng1(line))
      write (iunit,6100) 'REMARK DATAMAN bin plot '
      write (iunit,6100) 'REMARK Filename = ',
     +  pfile(1:leng1(pfile))
      write (iunit,6100) 'REMARK Dataset = ',
     +  name(iset)(1:leng1(name(iset))),' File = ',
     +  file(iset)(1:leng1(file(iset)))
      write (iunit,6100) 'REMARK Comment = ',
     +  coment(iset)(1:leng1(coment(iset)))
      write (iunit,6100) 'XLABEL ',hori(1:leng1(hori))
      write (iunit,6100) 'YLABEL ',vert(1:leng1(vert))
      write (iunit,6100) 'LINFIT'
      write (iunit,6120) 'XYVIEW ',pxmin,pxmax,pymin,pymax
      write (iunit,6100) '!'
      write (iunit,6110) 'NPOINT ',nbins
      write (iunit,6110) 'MRKTYP ',1
      write (iunit,6110) 'COLOUR ',4
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (xval(i),i=1,nbins)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (yval(i),i=1,nbins)
      write (iunit,6100) '!'
c
      if (know(kfree,iset)) then
        write (iunit,6100) 'MORE'
        write (iunit,6100) '!'
        write (iunit,6100) 'REMARK Work and Test data separated'
        write (iunit,6110) 'NPOINT ',nbins
        write (iunit,6110) 'MRKTYP ',2
        write (iunit,6110) 'COLOUR ',1
        write (iunit,6100) 'XVALUE *'
        write (iunit,'(1p,6e12.4)') (xval(i),i=1,nbins)
        write (iunit,6100) 'YVALUE *'
        write (iunit,'(1p,6e12.4)') (rval(i),i=1,nbins)
        write (iunit,6100) '!'
      end if
c
      write (iunit,6100) 'END'
c
      call prompt (' Plot file generated')
c
c 'FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|NRF|LNI'
c
 9000 continue
      if (iunit .gt. 0) close (iunit)
c
      return
c
 6000 format (1x,a3,1p,4(1x,a3,1x,e12.4))
 6100 format (10a)
 6110 format (a,i10)
 6120 format (a,1p,4e12.4)
c
      end
c
c
c
      subroutine duoplt (iset,jset,pfile,hori,vert,binsiz,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
c
      integer maxbin, minbin
      parameter (maxbin = 100, minbin = 5)
c
      integer iset,jset,ierr,i,iunit,leng1,nbins,j
      integer ninbin(maxbin),rinbin(maxbin)
c
      real x,y,z,xmin,xmax,ymin,ymax
      real pxmin,pxmax,pymin,pymax,binsiz
      real xval(maxbin),yval(maxbin),rval(maxbin)
      real rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2
c
      logical xinter
c
      character pfile*(*),hori*(*),vert*(*),line*80
      character esv1*12,esv2*12
c
code ...
c
      ierr = 0
      iunit = -1
c
      call textut (' Duo plot Set 1 =',name(iset))
      call textut ('          Set 2 =',name(jset))
c
c ... get X range & collect data
c
      do i=1,numhkl(iset)
        buffer (i) = 0.0
        buffer (maxhkl+i) = 0.0
      end do
c
      do i=1,numhkl(iset)
        x = 0.0
        y = 0.0
        if (hori(1:3) .eq. 'FOB') then
          x = fobs (i,iset)
        else if (hori(1:3) .eq. 'SIG') then
          x = sigfob (i,iset)
        else if (hori(1:3) .eq. 'F/S') then
          x = fobs (i,iset) / sigfob (i,iset)
        else if (hori(1:3) .eq. 'INT') then
          x = fobs (i,iset) ** 2
        else if (hori(1:3) .eq. 'LNI') then
          if (fobs(i,iset).gt.0.0) x = alog(fobs(i,iset)**2)
        else if (hori(1:3) .eq. 'I/S') then
          x = (fobs (i,iset) ** 2) / (sigfob (i,iset) ** 2)
        else if (hori(1:3) .eq. 'RES') then
          x = reso (i,iset)
        else if (hori(1:3) .eq. '1/R') then
          x = 1.0 / reso (i,iset)
        else if (hori(1:3) .eq. 'STL') then
          x = 0.5 / reso (i,iset)
        else if (hori(1:3) .eq. 'DST') then
          z = 0.5 / reso (i,iset)
          x = 4.0 * z * z
        else
          call errcon ('Invalid hori variable')
          ierr = -1
          goto 9000
        end if
c
        if (vert(1:3) .eq. 'FOB') then
          y = fobs (i,iset)
        else if (vert(1:3) .eq. 'SIG') then
          y = sigfob (i,iset)
        else if (vert(1:3) .eq. 'F/S') then
          y = fobs (i,iset) / sigfob (i,iset)
        else if (vert(1:3) .eq. 'INT') then
          y = fobs (i,iset) ** 2
        else if (vert(1:3) .eq. 'LNI') then
          if (fobs(i,iset).gt.0.0) y = alog(fobs(i,iset)**2)
        else if (vert(1:3) .eq. 'I/S') then
          y = (fobs (i,iset) ** 2) / (sigfob (i,iset) ** 2)
        else if (vert(1:3) .eq. 'RES') then
          y = reso (i,iset)
        else if (vert(1:3) .eq. '1/R') then
          y = 1.0 / reso (i,iset)
        else if (vert(1:3) .eq. 'STL') then
          y = 0.5 / reso (i,iset)
        else if (vert(1:3) .eq. 'DST') then
          z = 0.5 / reso (i,iset)
          y = 4.0 * z * z
        else if (vert(1:3) .eq. 'NRF') then
          y = 1.0
        else
          call errcon ('Invalid vert variable')
          ierr = -1
          return
        end if
c
        buffer (i) = x
        buffer (maxhkl+i) = y
c
        if (i .eq. 1) then
          xmin = x
          xmax = x
        else
          xmin = min (xmin, x)
          xmax = max (xmax, x)
        end if
      end do
c
      do i=1,numhkl(jset)
        buffer (2*maxhkl+i) = 0.0
        buffer (3*maxhkl+i) = 0.0
      end do
c
      do i=1,numhkl(jset)
        x = 0.0
        y = 0.0
        if (hori(1:3) .eq. 'FOB') then
          x = fobs (i,jset)
        else if (hori(1:3) .eq. 'SIG') then
          x = sigfob (i,jset)
        else if (hori(1:3) .eq. 'F/S') then
          x = fobs (i,jset) / sigfob (i,jset)
        else if (hori(1:3) .eq. 'INT') then
          x = fobs (i,jset) ** 2
        else if (hori(1:3) .eq. 'I/S') then
          x = (fobs (i,jset) ** 2) / (sigfob (i,jset) ** 2)
        else if (hori(1:3) .eq. 'LNI') then
          if (fobs(i,iset).gt.0.0) x = alog(fobs(i,iset)**2)
        else if (hori(1:3) .eq. 'RES') then
          x = reso (i,jset)
        else if (hori(1:3) .eq. '1/R') then
          x = 1.0 / reso (i,jset)
        else if (hori(1:3) .eq. 'STL') then
          x = 0.5 / reso (i,jset)
        else if (hori(1:3) .eq. 'DST') then
          z = 0.5 / reso (i,jset)
          x = 4.0 * z * z
        end if
c
        if (vert(1:3) .eq. 'FOB') then
          y = fobs (i,jset)
        else if (vert(1:3) .eq. 'SIG') then
          y = sigfob (i,jset)
        else if (vert(1:3) .eq. 'F/S') then
          y = fobs (i,jset) / sigfob (i,jset)
        else if (vert(1:3) .eq. 'INT') then
          y = fobs (i,jset) ** 2
        else if (vert(1:3) .eq. 'LNI') then
          if (fobs(i,iset).gt.0.0) y = alog(fobs(i,iset)**2)
        else if (vert(1:3) .eq. 'I/S') then
          y = (fobs (i,jset) ** 2) / (sigfob (i,jset) ** 2)
        else if (vert(1:3) .eq. 'RES') then
          y = reso (i,jset)
        else if (vert(1:3) .eq. '1/R') then
          y = 1.0 / reso (i,jset)
        else if (vert(1:3) .eq. 'STL') then
          y = 0.5 / reso (i,jset)
        else if (vert(1:3) .eq. 'DST') then
          z = 0.5 / reso (i,jset)
          y = 4.0 * z * z
        else if (vert(1:3) .eq. 'NRF') then
          y = 1.0
        end if
c
        buffer (2*maxhkl+i) = x
        buffer (3*maxhkl+i) = y
c
        xmin = min (xmin, x)
        xmax = max (xmax, x)
      end do
c
      write (*,'(1x,a3,1x,a,1x,1p,2e12.4)') hori,
     +  'Min, Max =',xmin,xmax
c
      write (*,'(1x,a,2i3)') 'Min and Max nr of bins = ',
     +  minbin,maxbin
c
      if (binsiz .gt. 0.0) then
        nbins = (xmax-xmin)/binsiz
      else
        nbins = max (1, nint(-binsiz))
      end if
      nbins = max (minbin, min(nbins,maxbin))
      call ivalut (' Nr of bins :',1,nbins)
      binsiz = (xmax-xmin)/float(nbins)
      call rvalut (' Bin size   :',1,binsiz)
c
      do i=1,nbins
        xval (i) = xmin + float(i-1)*binsiz + 0.5*binsiz
        yval (i) = 0.0
        rval (i) = 0.0
        ninbin (i) = 0
        rinbin (i) = 0
      end do
c
c ... collect data in bins
c
      do i=1,numhkl(iset)
        x = buffer (i)
        y = buffer (maxhkl+i)
        j = max (1, min (1 + int((x-xmin)/binsiz), nbins))
        ninbin (j) = ninbin (j) + 1
        yval (j) = yval (j) + y
      end do
c
      do i=1,numhkl(jset)
        x = buffer (2*maxhkl+i)
        y = buffer (3*maxhkl+i)
        j = max (1, min (1 + int((x-xmin)/binsiz), nbins))
        rinbin (j) = rinbin (j) + 1
        rval (j) = rval (j) + y
      end do
c
      if (vert(1:3) .ne. 'NRF') then
        do i=1,nbins
          if (ninbin(i) .gt. 0) then
            yval (i) = yval (i) / float(ninbin(i))
          end if
          if (rinbin(i) .gt. 0) then
            rval (i) = rval (i) / float(rinbin(i))
          end if
        end do
      end if
c
c ... print table
c
      write (*,'(a,a3,a,a3)') ' Plot ',vert,' versus ',hori
      esv1 = ' <'//vert(1:3)//'> Set 1'
      esv2 = ' <'//vert(1:3)//'> Set 2'
      write (*,6300) 'Bin nr','Start value',esv1,
     +  'Nr values',esv2,'Nr values'
c
      do i=1,nbins
        write (*,6310) i,xval(i),yval(i),ninbin(i),
     +    rval(i),rinbin(i)
      end do
c
 6300 format (1x,6a12)
 6310 format (1x,i12,1p,2e12.4,i12,e12.4,i12)
c
      call xystat (yval,rval,nbins,
     +  rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2)
      write (*,*)
      call prompt (' Comparison for Set 1 and Set 2 data :')
      call fvalut (' Correlation coefficient :',1,corr)
      call rvalut (' Scaled R w.r.t. Set 1   :',1,rf3)
      call rvalut (' Scaled R w.r.t. Set 2   :',1,rf4)
      call rvalut (' RMS difference          :',1,rmsd)
c
c ... write plot file
c
      iunit = 45
      call xopxua (iunit,pfile,xinter(),ierr)
      if (ierr .ne. 0) goto 9000
c
      z = 0.025 * (xmax-xmin)
      pxmin = xmin - z
      pxmax = xmax + z
      ymin = yval(1)
      ymax = yval(1)
      do i=1,nbins
        ymin = min (ymin,yval(i),rval(i))
        ymax = max (ymax,yval(i),rval(i))
      end do
      z = 0.025 * (ymax-ymin)
      pymin = ymin - z
      pymax = ymax + z
c
c ... write header
c
      call stamp (line)
      write (iunit,6100) 'REMARK ',line(1:leng1(line))
      write (iunit,6100) 'REMARK DATAMAN duo plot '
      write (iunit,6100) 'REMARK Filename = ',
     +  pfile(1:leng1(pfile))
      write (iunit,6100) 'REMARK Dataset 1 = ',
     +  name(iset)(1:leng1(name(iset))),' File = ',
     +  file(iset)(1:leng1(file(iset)))
      write (iunit,6100) 'REMARK Comment 1 = ',
     +  coment(iset)(1:leng1(coment(iset)))
      write (iunit,6100) 'REMARK Dataset 2 = ',
     +  name(jset)(1:leng1(name(jset))),' File = ',
     +  file(jset)(1:leng1(file(jset)))
      write (iunit,6100) 'REMARK Comment 2 = ',
     +  coment(jset)(1:leng1(coment(jset)))
      write (iunit,6100) '!'
      write (iunit,6100) 'XLABEL ',hori(1:leng1(hori))
      write (iunit,6100) 'YLABEL ',vert(1:leng1(vert))
      write (iunit,6100) 'LINFIT'
      write (iunit,6120) 'XYVIEW ',pxmin,pxmax,pymin,pymax
      write (iunit,6100) '!'
      write (iunit,6110) 'NPOINT ',nbins
      write (iunit,6110) 'MRKTYP ',1
      write (iunit,6110) 'COLOUR ',4
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (xval(i),i=1,nbins)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (yval(i),i=1,nbins)
      write (iunit,6100) '!'
c
      write (iunit,6100) 'MORE'
      write (iunit,6100) '!'
      write (iunit,6110) 'NPOINT ',nbins
      write (iunit,6110) 'MRKTYP ',2
      write (iunit,6110) 'COLOUR ',1
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (xval(i),i=1,nbins)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (rval(i),i=1,nbins)
      write (iunit,6100) '!'
c
      write (iunit,6100) 'END'
c
      call prompt (' Plot file generated')
c
c 'FOB|SIG|F/S|INT|I/S|RES|1/R|STL|DST|NRF|LNI'
c
 9000 continue
      if (iunit .gt. 0) close (iunit)
c
      return
c
 6000 format (1x,a3,1p,4(1x,a3,1x,e12.4))
 6100 format (99a)
 6110 format (a,i10)
 6120 format (a,1p,4e12.4)
c
      end
c
c
c
      subroutine merger (kset,iset,jset,how,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,rfree,buffer,
     +                   morbit,centri,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer rfree(maxhkl,maxset),morbit(maxhkl,maxset)
      character*1 centri(maxhkl,maxset)
c
      integer iset,jset,kset,ierr,i,new,icode,jcode
      integer n1,n2,n12,mode
c
      real rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2,top,bot,x1,x2
c
      character how*3
c
code ...
c
      ierr = -1
c
      call textut (' Merging Set 1 :',name(iset))
      call textut ('     and Set 2 :',name(jset))
      call textut (' Method        :',how)
      call prompt (' Sets *assumed* to be scaled together')
c
      mode = 3
      if (how .eq. 'AVE') mode = 1
      if (how .eq. 'SIG') mode = 2
c
      call prompt (' Encoding reflections of set 1 ...')
      do i=1,numhkl(iset)
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,icode)
        ibuff (i) = icode
        ibuff (maxhkl+i) = icode
        ibuff (4*maxhkl+i) = i
      end do
c
      call prompt (' Sorting reflections of set 1 ...')
      call shell (ibuff(maxhkl+1),ibuff(4*maxhkl+1),numhkl(iset))
c
      new = 0
      n1 = 0
      n2 = 0
      n12 = 0
c
      call prompt (' Generating merged dataset ...')
c
      do i=1,numhkl(jset)
c
        call packin (hkl(1,i,jset),hkl(2,i,jset),
     +    hkl(3,i,jset),0,jcode)
c
        call bindex (jcode,ibuff(maxhkl+1),
     +    ibuff(4*maxhkl+1),numhkl(iset),icode)
c
        if (icode .gt. 0) then
c
c ... reflection occurs in both datasets
c
          new = new + 1
          n12 = n12 + 1
          hkl (1,new,kset) = hkl(1,i,jset)
          hkl (2,new,kset) = hkl(2,i,jset)
          hkl (3,new,kset) = hkl(3,i,jset)
          if (mode .eq. 1) then
            fobs (new,kset) = 0.5*(fobs(icode,iset)+fobs(i,jset))
            sigfob (new,kset) = 0.5*sqrt(sigfob(icode,iset)**2 +
     +        sigfob(i,jset)**2)
          else if (mode .eq. 2) then
            fobs (new,kset) =
     +        (sigfob(i,jset)*fobs(icode,iset) + 
     +         sigfob(icode,iset)*fobs(i,jset)) /
     +        (sigfob(i,jset) + sigfob(icode,iset))
            sigfob (new,kset) = 2.0 * sigfob(icode,iset) *
     +        sigfob(i,jset) /
     +        (sigfob(icode,iset) + sigfob(i,jset))
          else if (mode .eq. 3) then
            fobs (new,kset) = fobs(icode,iset)
            sigfob (new,kset) = sigfob(icode,iset)
          end if
          reso (new,kset)   = reso (icode,iset)
          morbit (new,kset) = morbit (icode,iset)
          centri (new,kset) = centri (icode,iset)
          rfree (new,kset)  = 0
c
          buffer (2*maxhkl+n12) = fobs(icode,iset)
          buffer (3*maxhkl+n12) = fobs(i,jset)
c
c ... flag this reflection as used
c
          ibuff (icode) = -999
        else
c
c ... reflection only occurs in set 2
c
          new = new + 1
          n2 = n2 + 1
          hkl (1,new,kset) = hkl(1,i,jset)
          hkl (2,new,kset) = hkl(2,i,jset)
          hkl (3,new,kset) = hkl(3,i,jset)
          fobs (new,kset) = fobs(i,jset)
          sigfob (new,kset) = sigfob(i,jset)
          reso (new,kset)   = reso (i,jset)
          rfree (new,kset)  = 0
        end if
      end do
c
      call prompt (' Almost done ...')
c
      do i=1,numhkl(iset)
        if (ibuff(i) .gt. -999) then
c
c ... reflection only occurs in set 1
c
          new = new + 1
          n1 = n1 + 1
          hkl (1,new,kset) = hkl(1,i,iset)
          hkl (2,new,kset) = hkl(2,i,iset)
          hkl (3,new,kset) = hkl(3,i,iset)
          fobs (new,kset) = fobs(i,iset)
          sigfob (new,kset) = sigfob(i,iset)
          reso (new,kset)   = reso (i,iset)
          rfree (new,kset)  = 0
        end if
      end do
c
      write (*,*)
      call jvalut (' HKLs only in set 1 :',1,n1)
      call jvalut (' HKLs only in set 2 :',1,n2)
      call jvalut (' HKLs in both sets  :',1,n12)
      call jvalut (' Total nr of HKLs   :',1,new)
c
      if (n12 .gt. 0) then
        call xystat (buffer(2*maxhkl+1),buffer(3*maxhkl+1),n12,
     +    rmsd,shap,corr,rf1,rf2,rf3,rf4,sf1,sf2)
        write (*,*)
        call prompt (' Comparison for Set 1 and Set 2 Fobs :')
        call fvalut (' Correlation coefficient :',1,corr)
        call fvalut (' Shape similarity        :',1,shap)
        call rvalut (' Unscaled R w.r.t. <F1>  :',1,rf1)
        call rvalut (' Unscaled R w.r.t. <F2>  :',1,rf2)
        call rvalut (' Scaled R w.r.t. <F1>    :',1,rf3)
        call rvalut ('            Scale factor :',1,sf1)
        call rvalut (' Scaled R w.r.t. <F2>    :',1,rf4)
        call rvalut ('            Scale factor :',1,sf2)
        call rvalut (' RMS difference          :',1,rmsd)
      end if
c
      if (new .le. 0) return
c
      if (new .lt. nint(0.1*numhkl(iset)) .or.
     +    new .lt. nint(0.1*numhkl(jset)) ) then
        call prompt (' WARNING - few reflections in common !')
      end if
c
      top = 0.0
      bot = 0.0
c
      if (n12 .gt. 0) then
        do i=1,n12
          x1 = buffer(2*maxhkl+i)
          x2 = buffer(3*maxhkl+i)
          top = top + abs(x1 - x2)
          bot = bot + abs(x1 + x2)
        end do
c
        call prompt (' Rmerge = SUM |F1-F2| / SUM |F1+F2|')
        call fvalut (' Value of Rmerge :',1,(top/bot))
      end if
c
      numhkl (kset) = new
      ierr = 0
c
      return
      end
c
c
c
      subroutine effres (iset,lat,nasu,
     +                   maxset,maxhkl,fobs,sigfob)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
c
      integer nasu,iset,cnt(0:5),i,nz,j,k
c
      real res
c
      character lat*1
c
code ...
c
      do i=0,5
        cnt (i) = 0
      end do
c
      nz = 0
      do i=1,numhkl(iset)
        if (sigfob(i,iset) .le. 0.0) then
          nz = nz + 1
        else if (fobs(i,iset) .le. 0.0) then
          nz = nz + 1
        else
          j = int (fobs(i,iset)/sigfob(i,iset))
          j = max (0, min (5,j))
          do k=0,j
            cnt (k) = cnt (k) + 1
          end do
        end if
      end do
c
      if (nz .gt. 0) then
        call jvalut (
     +    ' Reflections with F <= 0 or Sigma <= 0 omitted :',1,nz)
      end if
c
      do i=0,5
        call estres (cell(1,iset),res,lat,nasu,cnt(i))
        write (*,6000) i,cnt(i),res
      end do
c
 6000 format (' Nr of HKLs with F >= ',i1,' * Sigma = ',i8,
     +  ' ==> Eff. D ~ ',f8.2,' A')
c
      return
      end
c
c
c
      subroutine compl (iset,res1,res2,lat,nasu,
     +                   maxset,maxhkl,reso)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real reso(maxhkl,maxset)
c
      integer iset,nasu,nz,n1,n2,i
c
      real res1,res2,co
c
      character lat*1
c
code ...
c
      call rlohi (res1,res2)
      nz = 0
      do i=1,numhkl(iset)
        if (reso(i,iset) .ge. res1 .and.
     +      reso(i,iset) .le. res2) nz = nz +1
      end do
c
      call estuni (cell(1,iset),res1,lat,nasu,n1)
      call estuni (cell(1,iset),res2,lat,nasu,n2)
      co = 100.0 * nz / float (n1 - n2)
      write (*,6000) nz,res1,res2,(n1-n2),co
c
 6000 format (1x,i8,' HKLs in [',f6.2,'-',f6.2,'] Max ~ ',
     +  i8,' Cmplt ~ ',f6.2,' %')
c
      return
      end
c
c
c
      subroutine minres (iset,lat,nasu,nncs,nres)
c
      include 'dataman.incl'
c
      integer nstrat
      parameter (nstrat = 9)
c
      real rho(nstrat),res
c
      integer iset,nasu,nncs,nres,i,j,nr
      integer npar(nstrat)
c
      character lat*1
      character stratxt(nstrat)*60
c
      data stratxt /
     +  'Rigid-body refinement       (6*Nncs)',
     +  'Torsion  /Grouped Bs/NCS    (~4*Nres)',
     +  'Torsion  /Grouped Bs/No NCS (~4*Nres*Nncs)',
     +  'Cartesian/Grouped Bs/NCS    (~26*Nres)',
     +  'Cartesian/Grouped Bs/No NCS (~26*Nres*Nncs)',
     +  'Cartesian/Isotr Bs  /NCS    (~32*Nres)',
     +  'Cartesian/Isotr Bs  /No NCS (~32*Nres*Nncs)',
     +  'Cartesian/Anisotr Bs/NCS    (~72*Nres)',
     +  'Cartesian/Anisotr Bs/No NCS (~72*Nres*Nncs)'/
c
code ...
c
      call textut (' Refinement strategy for set :',name(iset))
      call fvalut (' Unit cell axes  :',3,cell(1,iset))
      call asciut (' Lattice type    :',1,lat)
      call ivalut (' Nr asymm. units :',1,nasu)
      call ivalut (' Nr of residues  :',1,nres)
      call ivalut (' NCS molecules   :',1,nncs)
c
      npar (1) = 6 * nncs
      npar (2) = 4 * nres
      npar (3) = 4 * nres * nncs
      npar (4) = 26 * nres
      npar (5) = 26 * nres * nncs
      npar (6) = 32 * nres
      npar (7) = 32 * nres * nncs
      npar (8) = 72 * nres
      npar (9) = 72 * nres * nncs
c
      write (*,*)
      call prompt (' The following refinement strategies are used :')
      do i=1,nstrat
        write (*,6000) i,stratxt(i)
      end do
 6000 format (' Nr ',i2,' = ',a60)
c
      write (*,*)
      call prompt (' The optimal strategy depends on the resolution;')
      call prompt (' Nreflections should be > ~1.5 Nparameters !!!!!')
      call prompt (' The following table shows the MINIMUM effective')
      call prompt (' resolution for which this is the case for these')
      call prompt (' refinement strategies:')
      do i=1,nstrat
        nr = nint (1.5 * float(npar(i)))
        call estres (cell(1,iset),res,lat,nasu,nr)
        write (*,6010) i,npar(i),res
      end do
 6010 format (' Nr ',i2,' ~ ',i10,' parameters => Dmin ~ ',f6.2,' A')
c
      write (*,*)
      call prompt (' RHO = Nref / Npar is listed in the following')
      call prompt (' table as a function of EFFECTIVE resolution and')
      call prompt (' refinement strategy:')
c
      write (*,*)
      write (*,6020) 'Res(A)','Nrefl','RHO',(j,j=1,nstrat)
 6020 format (1x,a6,1x,a10,1x,a3,1x,9(i7))
 6030 format (1x,f6.2,1x,i10,5x,9(f7.1))
c
      do i=40,5,-1
        res = 0.1 * float(i)
        call estuni (cell(1,iset),res,lat,nasu,nr)
        do j=1,nstrat
          rho (j) = float (nr) / float (npar(j))
        end do
        write (*,6030) res,nr,(rho(j),j=1,nstrat)
      end do
c
      write (*,*)
c
      return
      end
c
c
c
      subroutine anisop (iset,pfile,ierr,
     +                   maxset,maxhkl,fobs,reso,hkl)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset),reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
c
      integer maxbin
      parameter (maxbin = 1000)
c
      integer iset,ierr,i,iunit,nbins,h,k,l,nh,nk,nl,leng1
      integer hinbin(0:maxbin),kinbin(0:maxbin),linbin(0:maxbin)
c
      real x,z,xmin,xmax,ymin,ymax
      real pxmin,pxmax,pymin,pymax
      real hdval(0:maxbin),kdval(0:maxbin),ldval(0:maxbin)
      real hfval(0:maxbin),kfval(0:maxbin),lfval(0:maxbin)
c
      logical xinter
c
      character pfile*(*),line*80
c
code ...
c
      ierr = 0
      iunit = -1
c
      call textut (' HKL-aniso plot Set =',name(iset))
c
      nbins = maxbin
c
      do i=0,nbins
        hdval (i) = 0.0
        kdval (i) = 0.0
        ldval (i) = 0.0
        hfval (i) = 0.0
        kfval (i) = 0.0
        lfval (i) = 0.0
        hinbin (i) = 0
        kinbin (i) = 0
        linbin (i) = 0
      end do
c
c ... store relevant data
c
      do i=1,numhkl(iset)
c
        h = iabs(hkl(1,i,iset))
        hinbin (h) = hinbin (h) + 1
        hdval (h) = hdval (h) + 1.0/reso (i,iset)
        hfval (h) = hfval (h) + fobs(i,iset)
c
        k = iabs(hkl(2,i,iset))
        kinbin (k) = kinbin (k) + 1
        kdval (k) = kdval (k) + 1.0/reso (i,iset)
        kfval (k) = kfval (k) + fobs(i,iset)
c
        l = iabs(hkl(3,i,iset))
        linbin (l) = linbin (l) + 1
        ldval (l) = ldval (l) + 1.0/reso (i,iset)
        lfval (l) = lfval (l) + fobs(i,iset)
c
      end do
c
      do i=0,nbins
        if (hinbin(i) .gt. 0) then
          x = hfval (i) / float( hinbin(i) )
          hfval (i) = alog (x)
          hdval (i) = hdval (i) / float ( hinbin(i) )
        end if
        if (kinbin(i) .gt. 0) then
          x = kfval (i) / float( kinbin(i) )
          kfval (i) = alog (x)
          kdval (i) = kdval (i) / float ( kinbin(i) )
        end if
        if (linbin(i) .gt. 0) then
ccc       print *,' L ',i,linbin(i),lfval(i),ldval(i)
          x = lfval (i) / float( linbin(i) )
          lfval (i) = alog (x)
          ldval (i) = ldval (i) / float ( linbin(i) )
        end if
      end do
c
c ... print table
c
      write (*,*)
      write (*,'(1x,a,a)')
     +  'HKL |   <1/R>h    LN<F>h       # |   <1/R>k    LN<F>k',
     +  '       # |   <1/R>l    LN<F>l       # |'
      do i=0,nbins
        if ( (hinbin(i)+kinbin(i)+linbin(i)).gt.0) then
          write (*,6310) i,hdval(i),hfval(i),hinbin(i),
     +      kdval(i),kfval(i),kinbin(i),
     +      ldval(i),lfval(i),linbin(i)
        end if
      end do
c
 6310 format (1x,i3,' | ',1p,3(2e10.3,i6,' | '))
c
c ... write plot file
c
      iunit = 45
      call xopxua (iunit,pfile,xinter(),ierr)
      if (ierr .ne. 0) goto 9000
c
      nh = -1
      nk = -1
      nl = -1
      xmin = 999.9
      xmax = -999.9
      ymin = 999.9
      ymax = -999.9
      do i=0,nbins
c
        if (hinbin(i).gt.20) then
          nh = nh + 1
          if (nh.lt.i) then
            hdval (nh) = hdval (i)
            hfval (nh) = hfval (i)
            hinbin (nh) = hinbin (i)
          end if
          xmin = min (xmin,hdval(nh))
          xmax = max (xmax,hdval(nh))
          ymin = min (ymin,hfval(nh))
          ymax = max (ymax,hfval(nh))
        end if
c
        if (kinbin(i).gt.20) then
          nk = nk + 1
          if (nk.lt.i) then
            kdval (nk) = kdval (i)
            kfval (nk) = kfval (i)
            kinbin (nk) = kinbin (i)
          end if
          xmin = min (xmin,kdval(nk))
          xmax = max (xmax,kdval(nk))
          ymin = min (ymin,kfval(nk))
          ymax = max (ymax,kfval(nk))
        end if
c
        if (linbin(i).gt.20) then
          nl = nl + 1
          if (nl.lt.i) then
            ldval (nl) = ldval (i)
            lfval (nl) = lfval (i)
            linbin (nl) = linbin (i)
          end if
          xmin = min (xmin,ldval(nl))
          xmax = max (xmax,ldval(nl))
          ymin = min (ymin,lfval(nl))
          ymax = max (ymax,lfval(nl))
        end if
      end do
c
      z = 0.025 * (xmax-xmin)
      pxmin = xmin - z
      pxmax = xmax + z
      z = 0.025 * (ymax-ymin)
      pymin = ymin - z
      pymax = ymax + z
c
c ... write header
c
      call stamp (line)
      write (iunit,6100) 'REMARK ',line(1:leng1(line))
      write (iunit,6100) 'REMARK DATAMAN HKL anisotropy plot '
      write (iunit,6100) 'REMARK Filename = ',
     +  pfile(1:leng1(pfile))
      write (iunit,6100) 'REMARK Dataset = ',
     +  name(iset)(1:leng1(name(iset))),' File = ',
     +  file(iset)(1:leng1(file(iset)))
      write (iunit,6100) 'REMARK Comment = ',
     +  coment(iset)(1:leng1(coment(iset)))
      write (iunit,6100) 'REMARK Only bins with > 20 reflections ',
     +  'have been included in the plot'
      write (iunit,6100) 'REMARK H=solid blue, K=dashed red ',
     +  'L=dotted green curve'
      write (iunit,6100)
     +  'XLABEL <1/R> in bins of constant |H|, |K|, |L|'
      write (iunit,6100) 'YLABEL LN<F> for these reflections'
      write (iunit,6120) 'XYVIEW ',pxmin,pxmax,pymin,pymax
      write (iunit,6100) '!'
c
      write (iunit,6110) 'NPOINT ',nh+1
      write (iunit,6110) 'MRKTYP ',1
      write (iunit,6110) 'COLOUR ',4
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (hdval(i),i=0,nh)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (hfval(i),i=0,nh)
      write (iunit,6100) '!'
c
      write (iunit,6100) 'MORE'
      write (iunit,6100) '!'
      write (iunit,6110) 'NPOINT ',nk+1
      write (iunit,6110) 'COLOUR ',1
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (kdval(i),i=0,nk)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (kfval(i),i=0,nk)
      write (iunit,6100) '!'
c
      write (iunit,6100) 'MORE'
      write (iunit,6100) '!'
      write (iunit,6110) 'NPOINT ',nl+1
      write (iunit,6110) 'COLOUR ',2
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (ldval(i),i=0,nl)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (lfval(i),i=0,nl)
      write (iunit,6100) '!'
c
      write (iunit,6100) 'END'
c
      call prompt (' Plot file generated')
c
 9000 continue
      if (iunit .gt. 0) close (iunit)
c
      return
c
 6000 format (1x,a3,1p,4(1x,a3,1x,e12.4))
 6100 format (10a)
 6110 format (a,i10)
 6120 format (a,1p,4e12.4)
c
      end
c
c
c
      subroutine rfbins (iset,maxset,maxhkl,rfree,reso,
     +                   buffer,newper,nbins,opt)
c
      include 'dataman.incl'
c
      integer maxbin
      parameter (maxbin = 100)
c
      integer maxset,maxhkl
      real reso(maxhkl,maxset),buffer(maxhkl)
      integer rfree(maxhkl,maxset)
c
      integer nwork(maxbin),ntest(maxbin),nflip(maxbin)
      integer iset,nbins,i,j,k,ibin
c
      real perc(maxbin),xcut(maxbin)
      real newper,xmin,xmax,dx,step
c
      character opt*2
c
code ...
c
      if (nbins .lt. 3 .or. nbins .gt. maxbin) then
        call errcon ('Invalid number of bins')
        return
      end if
c
      xmin = 1.0 / (reso(1,iset)*reso(1,iset))
      xmax = xmin
      do i=1,numhkl(iset)
        buffer (i) = 1.0 / (reso(i,iset)*reso(i,iset))
        if (buffer(i) .lt. xmin) xmin = buffer(i)
        if (buffer(i) .gt. xmax) xmax = buffer(i)
      end do
      call rvalut (' 4*(sin(theta)/lambda)**2 min :',1,xmin)
      call rvalut (' 4*(sin(theta)/lambda)**2 max :',1,xmax)
      call jvalut (' Nr of bins :',1,nbins)
      step = ( xmax - xmin ) / (float(nbins))
      call rvalut (' Bin size   :',1,step)
c
      do i=1,nbins
        nwork (i) = 0
        ntest (i) = 0
        nflip (i) = 0
        xcut  (i) = -1.0
      end do
c
      do i=1,numhkl(iset)
        ibin = 1 + int ( (buffer(i)-xmin)/step )
        ibin = max (1, min (ibin,nbins))
        if (rfree(i,iset) .eq. 0) then
          nwork (ibin) = nwork (ibin) + 1
        else
          ntest (ibin) = ntest (ibin) + 1
        end if
      end do
c
      if (opt .eq. 'FI') then
        call prompt (' Filling up bins ...')
        do i=1,nbins
          k = ntest(i) + nwork(i)
          if (k .gt. 0) then
            perc (i) = 100.0 * float(ntest(i)) / float(k)
            if (perc(i) .lt. newper) then
              j = nint (0.01 * newper * float(k))
              xcut (i) = float(j-ntest(i)) / float(k)
            end if
          else
            perc (i) = 0.0
          end if
        end do
c
        do i=1,numhkl(iset)
          ibin = 1 + int ( (buffer(i)-xmin)/step )
          ibin = max (1, min (ibin,nbins))
          if (xcut(ibin) .gt. 0.0) then
            if (rfree(i,iset) .eq. 0) then
              call gkrand (dx,0.0,1.0,0)
              if (dx .lt. xcut(ibin)) then
                rfree(i,iset) = 1
                nflip (ibin) = nflip (ibin) + 1
              end if
            end if
          end if
        end do
c
      else if (opt .eq. 'CU') then
        call prompt (' Cutting down bins ...')
        do i=1,nbins
          k = ntest(i) + nwork(i)
          if (k .gt. 0) then
            perc (i) = 100.0 * float(ntest(i)) / float(k)
            if (perc(i) .gt. newper) then
              j = nint (0.01 * newper * float(k))
              xcut (i) = float(j) / float(ntest(i))
            end if
          else
            perc (i) = 0.0
          end if
        end do
c
        do i=1,numhkl(iset)
          ibin = 1 + int ( (buffer(i)-xmin)/step )
          ibin = max (1, min (ibin,nbins))
          if (xcut(ibin) .gt. 0.0) then
            if (rfree(i,iset) .eq. 1) then
              call gkrand (dx,0.0,1.0,0)
              if (dx .ge. xcut(ibin)) then
                rfree(i,iset) = 0
                nflip (ibin) = nflip (ibin) + 1
              end if
            end if
          end if
        end do
c
      else if (opt .eq. 'BI') then
        do i=1,nbins
          k = ntest(i) + nwork(i)
          if (k .gt. 0) then
            perc (i) = 100.0 * float(ntest(i)) / float(k)
          else
            perc (i) = 0.0
          end if
        end do
      end if
c
      if (opt .eq. 'BI') then
        write (*,6000) 'Bin','4STOLSQ limits','Resol limits',
     +    'Nrefl','Ntest','%test'
      else
        write (*,6000) 'Bin','4STOLSQ limits','Resol limits',
     +    'Nrefl','Ntest','%test','New Ntest & %test'
      end if
c
 6000 format (/1x,a5,1x,a16,1x,a16,1x,a8,1x,a6,1x,a6,1x,a13)
c
      do i=1,nbins
c
        k = ntest(i) + nwork(i)
        if (opt .eq. 'FI') then
          j = ntest(i) + nflip(i)
        else if (opt .eq. 'CU') then
          j = ntest(i) - nflip(i)
        end if
c
        if (opt .eq. 'BI') then
          write (*,6100) i,
     +      xmin+float(i-1)*step,
     +      xmin+float(i)*step,
     +      1.0/sqrt(xmin+float(i-1)*step),
     +      1.0/sqrt(xmin+float(i)*step),
     +      k,ntest(i),perc(i)
        else
          write (*,6100) i,
     +      xmin+float(i-1)*step,
     +      xmin+float(i)*step,
     +      1.0/sqrt(xmin+float(i-1)*step),
     +      1.0/sqrt(xmin+float(i)*step),
     +      k,ntest(i),perc(i),j,
     +      100.0*float(j)/float(k)
        end if
c
      end do
c
      write (*,*)
c
 6100 format (1x,i5,1x,2f8.4,1x,2f8.3,1x,i8,1x,i6,1x,f6.2,
     +  1x,i6,1x,f6.2)
c
      return
      end
c
c
c
      subroutine rfshel (iset,rf,nb,
     +                   maxset,maxhkl,maxbuf,
     +                   reso,hkl,rfree,buffer,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer rfree(maxhkl,maxset)
c
      integer i,iset,nb,ioff,j,nsh,nrf,i0,i1,i2,k1,k2,k
c
      real rf,hrf
c
code ...
c
      ioff = numhkl(iset)
      call prompt (' Encoding reflections of this set ...')
      do i=1,numhkl(iset)
        buffer (ioff+i) = reso (i,iset)
        ibuff (i) = i
      end do
c
      call prompt (' Sorting reflections by resolution ...')
      call shellr (buffer(ioff+1),ibuff(1),numhkl(iset))
c
      nsh = nint ( float (numhkl(iset)) / float (nb) )
      nrf = int ( 0.01 * rf * float (nsh) )
      hrf = nint ( 0.5 * float (nrf) )
      call jvalut (' Nr of reflections        :',1,numhkl(iset))
      call jvalut (' Nr of resolution shells  :',1,nb)
      call jvalut (' Reflections per shell    :',1,nsh)
      call fvalut (' Percentage TEST reflect. :',1,rf)
      call jvalut (' Test reflections / shell :',1,nrf)
c
      do i=1,nb
c
c ... get reflection at centre of shell
c
        i0 = nint (( float(i - 1) + 0.5 ) *  float (nsh))
c
c ... get reflections at limits of resolution shell & print
c     a message
c
        i1 = max (1, nint (float(i-1)*float(nsh)))
        i2 = min (nint (float(i)*float(nsh)), numhkl(iset))
        k1 = ibuff (i1)
        k2 = ibuff (i2)
        write (*,6000) i,reso(k2,iset),reso(k1,iset)
c
c ... get reflections at limits of shell to be excised & print
c     a message
c
        i1 = max (1, i0 - int (hrf))
        i2 = min (i0 + int (hrf), numhkl(iset))
        k1 = ibuff (i1)
        k2 = ibuff (i2)
        write (*,6010) i,reso(k2,iset),reso(k1,iset),
     +    hkl(1,k2,iset),hkl(2,k2,iset),hkl(3,k2,iset),
     +    hkl(1,k1,iset),hkl(2,k1,iset),hkl(3,k1,iset)
c
c ... flag reflections in this shell as TEST reflections
c
        do j = i1, i2
          k = ibuff (j)
          rfree (k,iset) = 1
        end do
c
      end do
c
      write (*,*)
c
 6000 format (/
     +  ' -> Real shell # ',i4,' Resolution = ',f7.3,' A - ',f7.3,' A')
 6010 format (
     +  '    TEST Shell # ',i4,' Resolution = ',f7.3,' A - ',f7.3,' A'/
     +  '    First HKL = ',3i6/
     +  '    Last  HKL = ',3i6)
c
      return
      end
c
c
c
      subroutine rint (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,buffer,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
c
      integer maxmul
      parameter (maxmul=1000)
c
      real sumtel,sumnoe,ave
c
      integer iptr(maxmul),icnt(maxmul)
      integer iset,i,j,k,l,nn,ncont,nuni,nmult,nonce,nn2
c
code ...
c
      call jvalut (' Maximum multiplicity :',1,maxmul)
      call prompt (' Encoding reflections ...')
      nn = numhkl(iset)
      nn2 = 2*numhkl(iset)
      do i=1,numhkl(iset)
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,ibuff(i))
        ibuff (nn+i) = 0
        buffer (nn2+i) = fobs(i,iset)*fobs(i,iset)
      end do
      sumtel = 0.0
      sumnoe = 0.0
      ncont = 0
      nuni = 0
      nonce = 0
      nmult = 0
c
      do i=1,maxmul
        icnt (i) = 0
      end do
c
      call prompt (' Calculating Rint ...')
      do i=1,numhkl(iset)-1
        if (ibuff(nn+i) .eq. 0) then
          k = 1
          iptr (k) = i
          do j=i+1,numhkl(iset)
            if (ibuff(j).eq.ibuff(i)) then
              if (k .eq. maxmul) then
                call errcon ('Too high multiplicity !')
                goto 100
              else
                k = k + 1
                iptr (k) = j
              end if
            end if
          end do
 100      continue
          ibuff (nn+i) = k
          if (k .le. 1) then
            icnt (1) = icnt (1) + 1
            nonce = nonce + 1
            goto 200
          end if
          nmult = nmult + 1
          l = min (k,maxmul)
          icnt (l) = icnt(l) + 1
          ave = 0.0
          do l=1,k
            j = iptr(l)
            ave = ave + buffer(nn2+j)
          end do
          ave = ave / float(k)
ccc          print *,' *** ',hkl(1,i,iset),hkl(2,i,iset),
ccc     +       hkl(3,i,iset),k
ccc          call ivalut (' IPTR :',k,iptr)
          do l=1,k
            j = iptr (l)
            sumtel = sumtel + abs (buffer(nn2+j) - ave)
            sumnoe = sumnoe + buffer(nn2+j)
            ibuff (nn+j) = -1
ccc            print *,j,ave,buffer(nn2+j),sumtel,sumnoe
          end do
          ncont = ncont + k
 200      continue
          ibuff (nn+i) = k
          nuni = nuni + 1
        end if
      end do
c
      write (*,*)
      do i=1,maxmul
        if (icnt(i) .gt. 0) write (*,6000) i,icnt(i)
      end do
 6000 format (' Nr of reflexions with multiplicity ',i3,' = ',i10)
c
      write (*,*)
      call jvalut (' Nr of reflections :',1,numhkl(iset))
      call jvalut (' Nr of single obs  :',1,nonce)
      call jvalut (' Nr of mult obs    :',1,nmult)
      call jvalut ('   times they occur:',1,ncont)
      call jvalut (' Nr of unique refl :',1,nuni)
ccc      print *,sumtel,sumnoe
      write (*,*)
      if (nmult .gt. 0) then
        call prompt ('            Sum(hkl) Sum(i) | I - <I> |')
        call prompt (' Rint (I) = ---------------------------')
        call prompt ('                Sum(hkl) Sum(i) |I|')
        write (*,*)
        call rvalut (' Sum(hkl) Sum(i) | I - <I> | :',1,sumtel)
        call rvalut (' Sum(hkl) Sum(i) |I|         :',1,sumnoe)
        write (*,*)
        call fvalut (' Value of Rint (I) :',1,sumtel/sumnoe)
      else
        call prompt (' There are NO multiple observations !!!')
      end if
c
      return
      end
c
c
c
      subroutine hklkhl (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,hkl,buffer,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
c
      real x1,x2,x3,x4,x5,x6,x7,x8,x9,top,bot,topi,boti
c
      integer iset,i,new,jcode,icode,ni1,ni2,nhh,nsi
c
code ...
c
      call prompt (' Encoding reflections ...')
      do i=1,numhkl(iset)
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,icode)
        ibuff (i) = icode
        ibuff (3*maxhkl+i) = icode
        ibuff (4*maxhkl+i) = i
      end do
c
      call prompt (' Sorting reflections ...')
      call shell (ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),numhkl(iset))
c
      ni1 = maxhkl
      ni2 = 2*maxhkl
      new = 0
      nhh = 0
      nsi = 0
      call prompt (' Checking reflections ...')
c
      do i=1,numhkl(iset)-1
c
c ... pack KHL
c
        call packin (hkl(2,i,iset),hkl(1,i,iset),
     +    hkl(3,i,iset),0,jcode)
c
c ... HKL = KHL -> HHL
c
        if (jcode .eq. ibuff(i)) then
          nhh = nhh + 1
          goto 1854
        end if
c
        call bindex (jcode,ibuff(3*maxhkl+1),ibuff(4*maxhkl+1),
     +    numhkl(iset),icode)
c
c ... if ICODE < I -> already done
c
        if (icode .gt. 0 .and. icode .le. i) goto 1854
c
c ... if ICODE > 0 -> found symmetry mate
c
        if (icode .gt. 0) then
          new = new + 1
          ni1 = ni1 + 1
          ni2 = ni2 + 1
          buffer (ni1) = fobs (icode,iset)
          buffer (ni2) = fobs (i,iset)
        else
          nsi = nsi + 1
        end if
c
 1854   continue
      end do
c
      write (*,*)
      call jvalut (' Total nr of reflections   :',1,numhkl(iset))
      call jvalut (' Nr of HHL reflections     :',1,nhh)
      call jvalut (' Nr of single observations :',1,nsi)
      call jvalut (' Nr of HKL & KHL observ.   :',1,new)
      call jvalut (' Nr of reduced reflections :',1,(new+nhh+nsi))
c
      if (new .le. 0) return
c
      if (new .lt. nint(0.1*numhkl(iset))) then
        call prompt (' WARNING - few double observations !')
      end if
c
      write (*,*)
      call xystat (buffer(maxhkl+1),buffer(2*maxhkl+1),new,
     +             x1,x2,x3,x4,x5,x6,x7,x8,x9)
      call fvalut (' Correlation coeff Fobs :',1,x3)
      call fvalut (' Shape similarity  Fobs :',1,x2)
      call rvalut (' RMS difference Fhk/Fkh :',1,x1)
      call rvalut (' R=SUM(Fhk-Fkh)/SUM(Fhk):',1,x4)
      call rvalut (' R with (Fhk-S*Fkh)     :',1,x6)
      call rvalut ('          where scale S :',1,x8)
      call rvalut (' R=SUM(Fhk-Fkh)/SUM(Fkh):',1,x5)
      call rvalut (' R with (S*Fhk-Fkh)     :',1,x7)
      call rvalut ('          where scale S :',1,x9)
c
      top = 0.0
      bot = 0.0
      topi = 0.0
      boti = 0.0
c
      do i=1,new
        x1 = buffer(maxhkl+i)
        x2 = buffer(2*maxhkl+i)
        top = top + abs(x1 - x2)
        bot = bot + abs(x1 + x2)
        x1 = x1 * x1
        x2 = x2 * x2
        topi = topi + abs (x1 - x2)
        boti = boti + x1 + x2
      end do
c
      write (*,*)
      call prompt (' Rsym (F) = SUM |Fhkl-Fkhl| / SUM |Fhkl+Fkhl|')
      call fvalut (' Value of Rsym (F) :',1,(top/bot))
c
      write (*,*)
      call prompt (' Rsym (I) = SUM |Ihkl-Ikhl| / SUM |Ihkl+Ikhl|')
      call prompt (' Approximation: I = F*F')
      call fvalut (' Value of Rsym (I) :',1,(topi/boti))
c
      return
      end
c
c
c
      subroutine multip (iset,
     +                   maxset,maxhkl,maxbuf,
     +                   hkl,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
c
      real x1
c
      integer iset,i,icode,ncode,nuni,nmax
c
code ...
c
      call prompt (' Encoding reflections ...')
      do i=1,numhkl(iset)
        call packin (hkl(1,i,iset),hkl(2,i,iset),
     +    hkl(3,i,iset),0,icode)
        ibuff (maxhkl+i) = icode
        ibuff (2*maxhkl+i) = i
      end do
c
      call prompt (' Sorting reflections ...')
      call shell (ibuff(maxhkl+1),ibuff(2*maxhkl+1),numhkl(iset))
c
      do i=1,numhkl(iset)
        ibuff (i) = 0
      end do
c
      call prompt (' Checking reflections ...')
      icode = ibuff(maxhkl+1)
      ncode = 1
      nuni = 1
      nmax = 1
c
      do i=2,numhkl(iset)
        if (ibuff(maxhkl+i) .eq. icode) then
ccc          print *,' EQ ',i,ncode,icode,ibuff(maxhkl+i)
          ncode = ncode + 1
        else
ccc          print *,' NE ',i,ncode,icode,ibuff(maxhkl+i)
          ibuff (ncode) = ibuff(ncode) + 1
          if (ncode .gt. nmax) nmax = ncode
          icode = ibuff(maxhkl+i)
          ncode = 1
          nuni = nuni + 1
        end if
      end do
      ibuff (ncode) = ibuff(ncode) + 1
      if (ncode .gt. nmax) nmax = ncode
c
      write (*,*)
      call jvalut (' Total nr of reflections :',1,numhkl(iset))
      call jvalut (' Unique reflections      :',1,nuni)
      x1 = float(numhkl(iset)) / float (nuni)
      call fvalut (' Average redundancy      :',1,x1)
      call ivalut (' Maximum redundancy      :',1,nmax)
c
      if (nmax .gt. 1) then
        call prompt (' WARNING - Redundant reflections !!!')
      end if
c
      write (*,*)
      write (*,6000) 'Redundancy','Nr of unique reflections'
      write (*,6000) '----------','------------------------'
c
 6000 format (1x,a15,1x,a25)
 6010 format (1x,i10,5x,1x,10x,i10)
c
      do i=1,nmax
        if (i .le. 10 .or.
     +      (i .gt. 10 .and. ibuff(i) .gt. 0)) then
          write (*,6010) i,ibuff(i)
        end if
      end do
c
      return
      end
c
c
c
      subroutine sysabs (iset,mode,
     +                   maxset,maxhkl,maxbuf,fobs,sigfob,hkl,
     +                   rfree,reso,morbit,centri,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset),ibuff(maxbuf)
      integer morbit(maxhkl,maxset),rfree(maxhkl,maxset)
      character centri(maxhkl,maxset)*1
c
      real vecpro,q
c
      integer iset,i,j,k,res(3),ivec,nabs,new,nold
c
      logical lkill
c
      character mode*4
c
code ...
c
      nabs = 0
c
      lkill = (mode(1:1).eq.'k' .or. mode(1:1).eq.'K')
c
      if (numhkl(iset) .lt. 1) then
        call errcon ('No reflections')
        return
      end if
c
      if (nuniq(iset) .lt. 1) then
        call errcon ('No (unique) symmetry operators')
        return
      end if
c
      new = 0
      nold = numhkl(iset)
c
      do i=1,numhkl(iset)
c
        do j=1,nuniq(iset)
c
c ... apply symm-op G
c
          call imatmu (transp(1,j,iset),hkl(1,i,iset),res,3,3,1)
c
c ... check if G(h) = h
c
          if (hkl(1,i,iset) .ne. res(1)) goto 100
          if (hkl(2,i,iset) .ne. res(2)) goto 100
          if (hkl(3,i,iset) .ne. res(3)) goto 100
c
c ... G(h) = h; now test if (h.t) = integer
c
          vecpro = symmop(10,j,iset)*float(hkl(1,i,iset)) +
     +             symmop(11,j,iset)*float(hkl(2,i,iset)) +
     +             symmop(12,j,iset)*float(hkl(3,i,iset))
c
          ivec = nint(vecpro)
          if (abs(vecpro-float(ivec)) .le. 0.001) goto 100
c
c ... (h.t) not integer --> it is a systematic absence
c
          nabs = nabs + 1
          if (nabs .eq. 1) write (*,6980) 'Reflxn','H','K','L',
     +      'Fobs  ','Sigma  ','Test','F/Sigma'
c
          q = -1.0
          if (sigfob(i,iset).ne.0.0) q=fobs(i,iset)/sigfob(i,iset)
          write (*,6990) i,(hkl(k,i,iset),k=1,3),
     +          fobs(i,iset),sigfob(i,iset),rfree(i,iset),q
ccc     +          vecpro,ivec
c
          goto 200
c
 6980 format (a8,3a6,a12,a12,a5,a8)
 6990 format (i8,3i6,1p,2e12.4,i5,0p,f8.2)
c
  100     continue
        end do
c
c ... keep this one
c
        new = new + 1
        ibuff (new) = i
c
  200   continue
      end do
c
      call jvalut (' Nr of systematic absences :',1,nabs)
      if (.not. lkill) return
      if (nabs .eq. 0) return
c
      if (new .eq. 0) then
        call errcon ('This would kill all reflections; aborted')
        return
      end if
c
      do i=1,new
        j=ibuff(i)
        if (j .ne. i) then
          hkl(1,i,iset) = hkl(1,j,iset)
          hkl(2,i,iset) = hkl(2,j,iset)
          hkl(3,i,iset) = hkl(3,j,iset)
          fobs(i,iset)   = fobs(j,iset)
          sigfob(i,iset) = sigfob(j,iset)
          reso(i,iset)   = reso(j,iset)
          centri(i,iset) = centri(j,iset)
          morbit(i,iset) = morbit(j,iset)
          rfree(i,iset)  = rfree(j,iset)
        end if
      end do
c
      numhkl (iset) = new
      call jvalut (' Nr of reflections left :',1,new)
c
      i = nold - new
      if (i .gt. 0) then
        call jvalut (' WARNING - Reflections killed :',1,i)
      end if
c
      return
      end
c
c
c
      subroutine randno (iset,nb,xminoi,xmanoi,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,reso,buffer,ibuff)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer ibuff(maxbuf)
c
      integer i,iset,nb,ioff,j,nsh,i0,i1,i2,k1,k2,nn,ierr,jj
c
      real top,bot,topi,boti,xminoi,xmanoi,qq,ish,dx,dy,sign,old
      real bot1,boti1
c
code ...
c
      ioff = numhkl(iset)
      call prompt (' Copying & encoding reflections ...')
      do i=1,numhkl(iset)
        buffer (ioff+i) = reso (i,iset)
        ibuff (i) = i
      end do
c
      call prompt (' Sorting reflections by resolution ...')
      call shellr (buffer(ioff+1),ibuff(1),numhkl(iset))
c
      nsh = nint ( float (numhkl(iset)) / float (nb) )
      call jvalut (' Nr of reflections        :',1,numhkl(iset))
      call jvalut (' Nr of resolution shells  :',1,nb)
      call jvalut (' Reflections per shell    :',1,nsh)
c
      call fvalut (' Minimum noise %          :',1,xminoi)
      call fvalut (' Maximum noise %          :',1,xmanoi)
c
      top = 0.0
      bot = 0.0
      bot1 = 0.0
      topi = 0.0
      boti = 0.0
      boti1 = 0.0
c
      do i=1,nb
c
c ... get reflection at centre of shell
c
        i0 = nint (( float(i - 1) + 0.5 ) *  float (nsh))
c
c ... get reflections at limits of resolution shell & print
c     a message
c
        i1 = max (1, nint (float(i-1)*float(nsh)))
        i2 = min (nint (float(i)*float(nsh)), numhkl(iset))
        k1 = ibuff (i1)
        k2 = ibuff (i2)
        write (*,6000) i,reso(k2,iset),reso(k1,iset)
c
c ... get average I in shell
c
        ish = 0.0
        nn = 0
        do j=i1,i2
          jj = ibuff(j)
          ish = ish + fobs(jj,iset)*fobs(jj,iset)
          nn = nn + 1
        end do
c
        if (nn .le. 0) then
          call errcon ('No reflections in this shell ???')
          ierr = -1
          return
        end if
c
        ish = ish / float(nn)
        call jvalut (' Nr of reflection in shell:',1,nn)
        call rvalut (' Average intensity        :',1,ish)
c
c ... add random noise
c
        do j=i1,i2
          jj = ibuff (j)
c
          call gkrand (dx,xminoi,xmanoi,0)
c
c ... determine sign
c
          call gkrand (dy,-1.0,1.0,0)
          if (dy .le. 0.0) then
            sign = -1.0
          else
            sign = 1.0
          end if
c
          old = fobs(jj,iset)
          qq = fobs(jj,iset)*fobs(jj,iset) + sign*0.01*dx*ish
c
ccc          if (j.eq.i1) then
ccc            print *,sign,dx,ish,fobs(jj,iset)
ccc          end if
c
          if (qq .gt. 0.0) then
            fobs (jj,iset) = sqrt(qq)
          else
            fobs (jj,iset) = fobs(jj,iset) / 10.0
          end if
          top = top + abs(old-fobs(jj,iset))
          bot = bot + abs(old+fobs(jj,iset))
          bot1 = bot1 + abs(old)
          topi = topi + abs(old*old-fobs(jj,iset)*fobs(jj,iset))
          boti = boti + abs(old*old+fobs(jj,iset)*fobs(jj,iset))
          boti1 = boti1 + abs(old*old)
c
        end do
c
      end do
c
      write (*,*)
      ierr = 0
c
      call prompt (' Rmerge (F) = SUM |Fold-Fnew| / SUM |Fold+Fnew|')
      call fvalut (' Value of Rmerge (F) :',1,(top/bot))
      call prompt (' Rm" (F) = SUM |Fold-Fnew| / SUM |Fold|')
      call fvalut (' Value of Rm" (F) :',1,(top/bot1))
      call prompt (' Rmerge (I) = SUM |Iold-Inew| / SUM |Iold+Inew|')
      call fvalut (' Value of Rmerge (I) :',1,(topi/boti))
      call prompt (' Rm" (I) = SUM |Iold-Inew| / SUM |Iold|')
      call fvalut (' Value of Rm" (I) :',1,(topi/boti1))
c
 6000 format (/
     +  ' -> Real shell # ',i4,' Resolution = ',f7.3,' A - ',f7.3,' A')
c
      return
      end
c
c
c
      subroutine parity (iset,maxset,maxhkl,fobs,hkl)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
c
      real intens,ohk,ehk,ohl,ehl,okl,ekl,ohkl,ehkl,r1,r2,r3
      real oh,eh,ok,ek,ol,el
c
      integer iset,i,j,k,nohk,nohl,nokl,nohkl,noh,nok,nol
c
code ...
c
      if (numhkl(iset) .lt. 1) then
        call errcon ('No reflections')
        return
      end if
c
      oh = 0.0
      eh = 0.0
      noh = 0
      ok = 0.0
      ek = 0.0
      nok = 0
      ol = 0.0
      el = 0.0
      nol = 0
      ohk = 0.0
      ehk = 0.0
      nohk = 0
      ohl = 0.0
      ehl = 0.0
      nohl = 0
      okl = 0.0
      ekl = 0.0
      nokl = 0
      ohkl = 0.0
      ehkl = 0.0
      nohkl = 0
c
      do i=1,numhkl(iset)
c
        intens = fobs(i,iset) ** 2
c
        j = hkl(1,i,iset)
        if (j .ne. (2*(j/2))) then
          noh = noh + 1
          oh = oh + intens
        else
          eh = eh + intens
        end if
c
        j = hkl(2,i,iset)
        if (j .ne. (2*(j/2))) then
          nok = nok + 1
          ok = ok + intens
        else
          ek = ek + intens
        end if
c
        j = hkl(3,i,iset)
        if (j .ne. (2*(j/2))) then
          nol = nol + 1
          ol = ol + intens
        else
          el = el + intens
        end if
c
        j = hkl(1,i,iset)+hkl(2,i,iset)
        if (j .ne. (2*(j/2))) then
          nohk = nohk + 1
          ohk = ohk + intens
        else
          ehk = ehk + intens
        end if
c
        j = hkl(1,i,iset)+hkl(3,i,iset)
        if (j .ne. (2*(j/2))) then
          nohl = nohl + 1
          ohl = ohl + intens
        else
          ehl = ehl + intens
        end if
c
        j = hkl(2,i,iset)+hkl(3,i,iset)
        if (j .ne. (2*(j/2))) then
          nokl = nokl + 1
          okl = okl + intens
        else
          ekl = ekl + intens
        end if
c
        j = hkl(1,i,iset)+hkl(2,i,iset)+hkl(3,i,iset)
        if (j .ne. (2*(j/2))) then
          nohkl = nohkl + 1
          ohkl = ohkl + intens
        else
          ehkl = ehkl + intens
        end if
c
      end do
c
      write (*,*)
      call jvalut (' H odd  :',1,noh)
      call jvalut (' H even :',1,(numhkl(iset)-noh))
      if (noh .gt. 0 .and. noh .lt. numhkl(iset)) then
        r1 = (oh/float(noh))
        call gvalut (' <I(H odd)>  :',1,r1)
        r2 = (eh/float(numhkl(iset)-noh))
        call gvalut (' <I(H even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) A centering ???')
        end if
      end if
c
      write (*,*)
      call jvalut (' K odd  :',1,nok)
      call jvalut (' K even :',1,(numhkl(iset)-nok))
      if (nok .gt. 0 .and. nok .lt. numhkl(iset)) then
        r1 = (ok/float(nok))
        call gvalut (' <I(K odd)>  :',1,r1)
        r2 = (ek/float(numhkl(iset)-nok))
        call gvalut (' <I(K even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) B centering ???')
        end if
      end if
c
      write (*,*)
      call jvalut (' L odd  :',1,nol)
      call jvalut (' L even :',1,(numhkl(iset)-nol))
      if (nol .gt. 0 .and. nol .lt. numhkl(iset)) then
        r1 = (ol/float(nol))
        call gvalut (' <I(L odd)>  :',1,r1)
        r2 = (el/float(numhkl(iset)-nol))
        call gvalut (' <I(L even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) C centering ???')
        end if
      end if
c
      k = 0
c
      write (*,*)
      call jvalut (' H+K odd  :',1,nohk)
      call jvalut (' H+K even :',1,(numhkl(iset)-nohk))
      if (nohk .gt. 0 .and. nohk .lt. numhkl(iset)) then
        r1 = (ohk/float(nohk))
        call gvalut (' <I(H+K odd)>  :',1,r1)
        r2 = (ehk/float(numhkl(iset)-nohk))
        call gvalut (' <I(H+K even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) C-face centering ???')
          k = k + 1
        end if
      end if
c
      write (*,*)
      call jvalut (' H+L odd  :',1,nohl)
      call jvalut (' H+L even :',1,(numhkl(iset)-nohl))
      if (nohl .gt. 0 .and. nohl .lt. numhkl(iset)) then
        r1 = (ohl/float(nohl))
        call gvalut (' <I(H+L odd)>  :',1,r1)
        r2 = (ehl/float(numhkl(iset)-nohl))
        call gvalut (' <I(H+L even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) B-face centering ???')
          k = k + 1
        end if
      end if
c
      write (*,*)
      call jvalut (' K+L odd  :',1,nokl)
      call jvalut (' K+L even :',1,(numhkl(iset)-nokl))
      if (nokl .gt. 0 .and. nokl .lt. numhkl(iset)) then
        r1 = (okl/float(nokl))
        call gvalut (' <I(K+L odd)>  :',1,r1)
        r2 = (ekl/float(numhkl(iset)-nokl))
        call gvalut (' <I(K+L even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) A-face centering ???')
          k = k + 1
        end if
      end if
c
      if (k .eq. 3) then
        call prompt (' (Pseudo) F (all-face) centering ???')
      end if
c
      write (*,*)
      call jvalut (' H+K+L odd  :',1,nohkl)
      call jvalut (' H+K+L even :',1,(numhkl(iset)-nohkl))
      if (nohkl .gt. 0 .and. nohkl .lt. numhkl(iset)) then
        r1 = (ohkl/float(nohkl))
        call gvalut (' <I(H+K+L odd)>  :',1,r1)
        r2 = (ehkl/float(numhkl(iset)-nohkl))
        call gvalut (' <I(H+K+L even)> :',1,r2)
        r3 = r1 / r2
        call fvalut (' Ratio :',1,r3)
        if (r3 .le. 0.8) then
          call prompt (' (Pseudo) I (body) centering ???')
        end if
      end if
c
      return
      end
c
c
c
      subroutine fillin (iset,maxset,maxhkl,fobs,sigfob,
     +                   reso,buffer,nbins,ierr)
c
      include 'dataman.incl'
c
      integer maxbin
      parameter (maxbin = 100)
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxhkl)
c
      real aveint(maxbin)
      real xmin,xmax,dx,step,q1,q2
c
      integer nref(maxbin),nnot(maxbin)
      integer iset,nbins,i,ibin,ierr,nok,nfi
c
      logical lmiss
c
code ...
c
      ierr = -1
c
      if (nbins .lt. 3 .or. nbins .gt. maxbin) then
        call errcon ('Invalid number of bins')
        return
      end if
c
      xmin = 1.0 / (reso(1,iset)*reso(1,iset))
      xmax = xmin
      do i=1,numhkl(iset)
        buffer (i) = 1.0 / (reso(i,iset)*reso(i,iset))
        if (buffer(i) .lt. xmin) xmin = buffer(i)
        if (buffer(i) .gt. xmax) xmax = buffer(i)
      end do
      call rvalut (' 4*(sin(theta)/lambda)**2 min :',1,xmin)
      call rvalut (' 4*(sin(theta)/lambda)**2 max :',1,xmax)
      call jvalut (' Nr of bins :',1,nbins)
      step = ( xmax - xmin ) / (float(nbins))
      call rvalut (' Bin size   :',1,step)
c
      do i=1,nbins
        nref (i) = 0
        nnot (i) = 0
        aveint  (i) = 0.0
      end do
c
      nok = 0
      nfi = 0
      do i=1,numhkl(iset)
        ibin = 1 + int ( (buffer(i)-xmin)/step )
        ibin = max (1, min (ibin,nbins))
        if (sigfob(i,iset) .le. small) then
          nfi = nfi + 1
          nnot (ibin) = nnot (ibin) + 1
        else
          nok = nok + 1
          nref (ibin) = nref (ibin) + 1
          aveint (ibin) = aveint (ibin) +
     +                    fobs(i,iset)*fobs(i,iset)
        end if
      end do
c
      write (*,6000) 'Bin','4STOLSQ limits',
     +  'Nobs','Nfill','<Ibin>','Ffill'
c
      lmiss = .false.
      q1 = 999.9
      q2 = -999.9
      do i=1,nbins
        if (nref(i) .gt. 0) then
          aveint (i) = aveint (i) / float(nref(i))
        end if
        dx = sqrt (aveint(i))
        write (*,6010) i,xmin+float(i-1)*step,
     +    xmin+float(i)*step,nref(i),nnot(i),aveint(i),dx
        aveint (i) = dx
        dx = 0.0
        if (nnot(i).ne.0) dx = float(nref(i))/float(nnot(i))
        lmiss = (nnot(i) .gt. 0 .and. nref(i) .le. 0)
        if (nnot(i) .gt. 0) then
          if (dx .lt. q1) q1 = dx
          if (dx .gt. q2) q2 = dx
        end if
      end do
c
 6000 format (1x,a5,1x,a17,2(1x,a8),2(1x,a12))
 6010 format (1x,i5,2(1x,f8.5),2(1x,i8),2(1x,1pe12.4))
c
      call jvalut (' Nr of measured reflections   :',1,nok)
      call jvalut (' Nr of reflections to fill in :',1,nfi)
      if (nfi .le. 0) return
c
      if (lmiss) then
        call errcon (
     +  'At least one shell without observations to fill !!!')
        return
      end if
c
      call fvalut (' Min Nobs/Nfill ratio :',1,q1)
      call fvalut (' Max Nobs/Nfill ratio :',1,q2)
      if (q1 .lt. 2.0) then
        call prompt (' WARNING - Min Nobs/Nfill < 2 !')
      end if
c
      do i=1,numhkl(iset)
        ibin = 1 + int ( (buffer(i)-xmin)/step )
        ibin = max (1, min (ibin,nbins))
        if (sigfob(i,iset) .le. small) then
          fobs (i,iset) = aveint(ibin)
          sigfob (i,iset) = aveint(ibin)
        end if
      end do
c
      call prompt (' Done')
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine eoplot (iset,pfile,horkorl,ierr,
     +                   maxset,maxhkl,maxbuf,fobs,hkl)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset)
      integer hkl(3,maxhkl,maxset)
c
      integer maxbin
      parameter (maxbin = 512)
c
      integer iset,ierr,i,iunit,ndata,j,leng1,ihkl,imin,imax
      integer ninbin(-maxbin:maxbin)
c
      real xmin,xmax
      real xval(-maxbin:maxbin)
      real xdata(maxbin+1),ydata(maxbin+1)
c
      logical xinter
c
      character pfile*(*),horkorl*(*),line*80
c
code ...
c
      ierr = 0
      iunit = -1
c
      call textut (' Even/Odd plot set =',name(iset))
      call textut (' Miller index type =',horkorl)
      call textut (' Plot file name    =',pfile)
c
      if (horkorl(1:1) .eq. 'H') then
        ihkl = 1
      else if (horkorl(1:1) .eq. 'K') then
        ihkl = 2
      else
        ihkl = 3
      end if
c
      do i=-maxbin,maxbin
        ninbin (i) = 0
        xval (i) = 0
      end do
c
      do i=1,numhkl(iset)
        j = hkl (ihkl,i,iset)
        ninbin (j) = ninbin (j) + 1
        xval (j) = xval (j) + fobs(i,iset)
      end do
c
      imax = 0
      imin = 0
      xmin = 0.0
      xmax = -1.0E6
      do i=-maxbin,maxbin
        if (ninbin(i) .gt. 0) then
          if (imin .eq. 0) imin = i
          imax = i
          xval (i) = xval(i) / float(ninbin(i))
          if (xval(i) .gt. xmax) xmax = xval(i)
        end if
      end do
c
      xmax = 1.05 * xmax
      j = max (-imin,imax)
      imin = - j - 1
      imax = j + 1
c
      ndata = 0
      do i=imin,imax
        if (ninbin(i) .gt. 0) then
          if (i .eq. 2*(i/2)) then
            ndata = ndata + 1
            xdata (ndata) = float(i)
            ydata (ndata) = xval(i)
          end if
        end if
      end do
c
c ... write plot file
c
      iunit = 45
      call xopxua (iunit,pfile,xinter(),ierr)
      if (ierr .ne. 0) goto 9000
c
c ... write header
c
      call stamp (line)
      write (iunit,6100) 'REMARK ',line(1:leng1(line))
      write (iunit,6100) 'REMARK DATAMAN Even/Odd plot '
      write (iunit,6100) 'REMARK Filename = ',
     +  pfile(1:leng1(pfile))
      write (iunit,6100) 'REMARK Dataset = ',
     +  name(iset)(1:leng1(name(iset))),' File = ',
     +  file(iset)(1:leng1(file(iset)))
      write (iunit,6100) 'REMARK Comment = ',
     +  coment(iset)(1:leng1(coment(iset)))
      write (iunit,6100) 'XLABEL Miller index ',horkorl(1:1)
      write (iunit,6100) 'YLABEL <Fobs>'
      write (iunit,6120) 'XYVIEW ',float(imin),float(imax),xmin,xmax
      write (iunit,6100) '!'
c
      write (iunit,6100) 'REMARK EVEN indices: solid blue curve'
      write (iunit,6110) 'NPOINT ',ndata
      write (iunit,6110) 'MRKTYP ',1
      write (iunit,6110) 'COLOUR ',4
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (xdata(i),i=1,ndata)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (ydata(i),i=1,ndata)
      write (iunit,6100) '!'
c
      ndata = 0
      do i=imin,imax
        if (ninbin(i) .gt. 0) then
          if (i .ne. 2*(i/2)) then
            ndata = ndata + 1
            xdata (ndata) = float(i)
            ydata (ndata) = xval(i)
          end if
        end if
      end do
c
      write (iunit,6100) 'MORE'
      write (iunit,6100) '!'
      write (iunit,6100) 'REMARK ODD indices: dashed red curve'
      write (iunit,6110) 'NPOINT ',ndata
      write (iunit,6110) 'MRKTYP ',2
      write (iunit,6110) 'COLOUR ',1
      write (iunit,6100) 'XVALUE *'
      write (iunit,'(1p,6e12.4)') (xdata(i),i=1,ndata)
      write (iunit,6100) 'YVALUE *'
      write (iunit,'(1p,6e12.4)') (ydata(i),i=1,ndata)
      write (iunit,6100) '!'
c
      write (iunit,6100) 'END'
c
      call prompt (' Plot file generated')
c
 9000 continue
      if (iunit .gt. 0) close (iunit)
c
      return
c
 6000 format (1x,a3,1p,4(1x,a3,1x,e12.4))
 6100 format (10a)
 6110 format (a,i10)
 6120 format (a,1p,4e12.4)
c
      end
c
c
c
      subroutine reclat (iset,iunit,filenm,ramper,ierr,
     +                   maxset,maxhkl,maxbuf,
     +                   fobs,sigfob,reso,hkl,buffer)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl,maxbuf
      real fobs(maxhkl,maxset),sigfob(maxhkl,maxset)
      real reso(maxhkl,maxset),buffer(maxbuf)
      integer hkl(3,maxhkl,maxset)
c
      real xmin,xmax,h1,s1,v1,h2,s2,v2,h3,delta,r,g,b,h
c
      integer iset,ierr,iunit,hklmax,i,ncol,j,leng1,k
c
      logical xinter,lramp
c
      character filenm*(*),ramper*(*)
      character line*128
c
code ...
c
      call textut (' Ramp_odl :',name(iset))
      call textut (' ODL file :',filenm)
      call upcase (ramper)
      call remspa (ramper)
      call textut (' Ramp by  :',ramper)
c
      call xopxua (iunit,filenm,xinter(),ierr)
      if (ierr .ne. 0) return
c
 6000 format (a,3(1x,i15),1x,a,1x,i15)
c
      hklmax = 0
      do i=1,numhkl(iset)
        hklmax = max (hklmax, abs(hkl(1,i,iset)),
     +                abs(hkl(2,i,iset)), abs(hkl(3,i,iset)))
      end do
      call jvalut (' Max Abs (HKL) =',1,hklmax)
      hklmax = 10 * (1 + int (0.1 * float(hklmax)))
c
      write (iunit,6000,err=900) 'begin reclat'
      write (iunit,6000,err=900) 'colour cyan'
      write (iunit,6000,err=900) 'text_colour cyan'
c
c ... H axis
c
      write (line,6000,err=900) 'move ',-hklmax,0,0
      call pretty (line)
      write (iunit,6000,err=900) line(1:leng1(line))
      write (line,6000,err=900) 'line ',hklmax,0,0
      call pretty (line)
      write (iunit,6000,err=900) line(1:leng1(line))
      do i=-hklmax,hklmax,10
        write (line,6000,err=900) 'text ',i,0,0,'H = ',i
        call pretty (line)
        write (iunit,6000,err=900) line(1:leng1(line))
      end do
c
c ... K axis
c
      write (line,6000,err=900) 'move ',0,-hklmax,0
      call pretty (line)
      write (iunit,6000,err=900) line(1:leng1(line))
      write (line,6000,err=900) 'line ',0,hklmax,0
      call pretty (line)
      write (iunit,6000,err=900) line(1:leng1(line))
      do i=-hklmax,hklmax,10
        write (line,6000,err=900) 'text ',0,i,0,'K = ',i
        call pretty (line)
        write (iunit,6000,err=900) line(1:leng1(line))
      end do
c
c ... L axis
c
      write (line,6000,err=900) 'move ',0,0,-hklmax
      call pretty (line)
      write (iunit,6000,err=900) line(1:leng1(line))
      write (line,6000,err=900) 'line ',0,0,hklmax
      call pretty (line)
      write (iunit,6000,err=900) line(1:leng1(line))
      do i=-hklmax,hklmax,10
        write (line,6000,err=900) 'text ',0,0,i,'L = ',i
        call pretty (line)
        write (iunit,6000,err=900) line(1:leng1(line))
      end do
c
      write (iunit,6000,err=900) 'colour magenta'
c
      lramp = .false.
c
      if (ramper(1:3) .eq. 'FOB') then
        xmin = fobs(1,iset)
        xmax = xmin
        do i=1,numhkl(iset)
          buffer (i) = fobs (i,iset)
          xmin = min (xmin, buffer(i))
          xmax = max (xmax, buffer(i))
        end do
        call prompt (' Ramp by Fobs')
        call rvalut (' Minimum :',1,xmin)
        call rvalut (' Maximum :',1,xmax)
        if (xmin .lt. xmax) lramp = .true.
c
      else if (ramper(1:3) .eq. 'SIG') then
        xmin = sigfob(1,iset)
        xmax = xmin
        do i=1,numhkl(iset)
          buffer (i) = sigfob (i,iset)
          xmin = min (xmin, buffer(i))
          xmax = max (xmax, buffer(i))
        end do
        call prompt (' Ramp by Sigma(Fobs)')
        call rvalut (' Minimum :',1,xmin)
        call rvalut (' Maximum :',1,xmax)
        if (xmin .lt. xmax) lramp = .true.
c
      else if (ramper(1:3) .eq. 'F/S') then
        xmin = fobs(1,iset)/sigfob(1,iset)
        xmax = xmin
        do i=1,numhkl(iset)
          buffer (i) = fobs (i,iset) / sigfob (i,iset)
          xmin = min (xmin, buffer(i))
          xmax = max (xmax, buffer(i))
        end do
        call prompt (' Ramp by Fobs / Sigma(Fobs)')
        call rvalut (' Minimum :',1,xmin)
        call rvalut (' Maximum :',1,xmax)
        if (xmin .lt. xmax) lramp = .true.
c
      else if (ramper(1:3) .eq. 'RES') then
        if (know(kreso,iset)) then
          xmin = reso(1,iset)
          xmax = xmin
          do i=1,numhkl(iset)
            buffer (i) = reso (i,iset)
            xmin = min (xmin, buffer(i))
            xmax = max (xmax, buffer(i))
          end do
          call prompt (' Ramp by Resolution (A)')
          call rvalut (' Minimum :',1,xmin)
          call rvalut (' Maximum :',1,xmax)
          if (xmin .lt. xmax) lramp = .true.
        else
          call errcon ('Resolution not calculated; cannot ramp')
        end if
      else if (ramper(1:3) .eq. 'NON') then
        call prompt (' No colour ramping')
      else
        call errcon ('Invalid ramping criterion')
        call prompt (' No colour ramping')
      end if
c
      if (lramp) then
        call prompt (' Will do colour ramping:')
        call prompt (' From BLUE for low  values')
        call prompt (' To   RED  for high values')
c
        call rgbhsv (0.0,0.0,1.0,h1,s1,v1)
        call rgbhsv (1.0,0.0,0.0,h2,s2,v2)
        if (h1 .eq. h2) h2 = h2 + 180.0
        h3 = (h2 - h1) / 100.0
        delta = (h2-h1) / (xmax-xmin)
      end if
c
      if (lramp) then
        do i=1,numhkl(iset)
          k = int(float(int((buffer(i)-xmin)*delta/h3))*h3 + h1) 
          h = float(k)
          if (h2 .gt. h1) then
            h = max (h1, min (h2, h))
          else
            h = max (h2, min (h1, h))
          end if
          call fix360 (h)
          call hsvrgb (h,1.0,1.0,r,g,b)
          call rgb2o (r,g,b,ncol)
c
cc      if (buffer(i) .ge. (xmax-0.5) .or.
cc     +    buffer(i) .le. (xmin+0.5)) then
cc        write (*,'(3i5,3f10.4,i15)') (hkl(j,i,iset),j=1,3),
cc     +    buffer(i),blue,red,ncol
cc      end if
c
          write (line,6000,err=900) 'colour',ncol
          call pretty (line)
          write (iunit,6000,err=900) line(1:leng1(line))
c
          write (line,6000,err=900) 'dot ',
     +      (hkl(j,i,iset),j=1,3)
          call pretty (line)
          write (iunit,6000,err=900) line(1:leng1(line))
        end do
      else
        do i=1,numhkl(iset)
          write (line,6000,err=900) 'dot ',
     +      (hkl(j,i,iset),j=1,3)
          call pretty (line)
          write (iunit,6000,err=900) line(1:leng1(line))
        end do
      end if
c
      write (iunit,6000,err=900) 'end_object'
      close (iunit)
c
      call prompt (' ODL file written')
      ierr = 0
c
      return
c
  900 continue
      call errcon ('While writing ODL file')
      close (iunit)
      ierr = -1
c
      return
      end
c
c
c
      subroutine dosigs (mode,iset,ierr,maxset,maxhkl,fobs,sigfob,
     +                   sigmin,sigmax)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset)
      real sigfob(maxhkl,maxset)
      real sigmin,sigmax
c
      integer mode,i,ierr,iset,nmin,nmax
c
code ...
c
      ierr = 1
c
c ... mode = 1 - FAKE
c
      if (mode .eq. 1) then
        call prompt (' Faking sigmas ... Sigma = SQRT( |Fobs| )')
        do i=1,numhkl(iset)
          sigfob (i,iset) = sqrt (abs (fobs(i,iset)))
        end do
c
        ierr = 0
c
c ... mode = 2 - LIMIT
c
      else if (mode .eq. 2) then
        call prompt (' Limiting sigmas ...')
        call rvalut (' Minimum value :',1,sigmin)
        call rvalut (' Maximum value :',1,sigmax)
        nmin = 0
        nmax = 0
        do i=1,numhkl(iset)
          if (sigfob (i,iset) .lt. sigmin) then
            nmin = nmin + 1
            sigfob (i,iset) = sigmin
          else if (sigfob (i,iset) .gt. sigmax) then
            nmax = nmax + 1
            sigfob (i,iset) = sigmax
          end if
        end do
        call jvalut (' Reset < minimum :',1,nmin)
        call jvalut (' Reset > maximum :',1,nmax)
c
        ierr = 0
c
c ... unknown mode
c
      else
        call errcon ('DOSIGS - Unknown mode !?')
        ierr = -1
      end if
c
      return
      end
c
c
c
      subroutine sigace (iset,maxset,maxhkl,fobs,sigfob,centri)
c
      include 'dataman.incl'
c
      integer maxset,maxhkl
      real fobs(maxhkl,maxset)
      real sigfob(maxhkl,maxset)
      real suma,sumc,sumsa,sumsc
      real*8 duma,dumc,dumsa,dumsc
c
      integer i,iset,nc,na,nz
c
      character centri(maxhkl,maxset)*1
c
code ...
c
      nc = 0
      na = 0
      nz = 0
c
      duma = 0.0D0
      dumc = 0.0D0
      dumsa = 0.0D0
      dumsc = 0.0D0
c
      do i=1,numhkl(iset)
        if (sigfob(i,iset) .le. 0.0) then
          nz = nz + 1
        else
          if (centri(i,iset) .eq. 'A') then
            na = na + 1
            duma = duma + dble(abs(fobs(i,iset)/sigfob(i,iset)))
            dumsa = dumsa + dble(sigfob(i,iset))
          else
            nc = nc + 1
            dumc = dumc + dble(abs(fobs(i,iset)/sigfob(i,iset)))
            dumsc = dumsc + dble(sigfob(i,iset))
          end if
        end if
      end do
c
      call jvalut (' Nr of reflections  :',1,numhkl(iset))
      call jvalut (' Nr with Sigma <= 0 :',1,nz)
      call jvalut (' Nr acentric        :',1,na)
      if (na .gt. 0) then
        duma = duma / dble(na)
        suma = duma
        call fvalut (' <Fobs/Sigma(Fobs)> :',1,suma)
        dumsa = dumsa / dble(na)
        sumsa = dumsa
        call fvalut (' <Sigma(Fobs)>      :',1,sumsa)
      end if
      call jvalut (' Nr centric         :',1,nc)
      if (nc .gt. 0) then
        dumc = dumc / dble(nc)
        sumc = dumc
        call fvalut (' <Fobs/Sigma(Fobs)> :',1,sumc)
        dumsc = dumsc / dble(nc)
        sumsc = dumsc
        call fvalut (' <Sigma(Fobs)>      :',1,sumsc)
      end if
c
      return
      end
c
c
c
      subroutine bindex (icode,jcodes,jindex,numj,m)
c
c ... use bisecting to see if ICODE occurs in the sorted array
c     of encoded reflections. this requires only 2LOG(N) steps
c     instead of N/2 steps for a naive loop over all reflections
c
      implicit none
c
      integer icode,numj,m
      integer jcodes(numj),jindex(numj)
c
      integer jlo,klo,jhi,khi,jxx,kxx
c
code...
c
      m = -1
c
      jlo = 1
      klo = jcodes(jlo) - icode
      if (klo .eq. 0) then
        m = jindex(jlo)
        return
      end if
      if (klo .gt. 0) return
c
      jhi = numj
      khi = jcodes(jhi) - icode
      if (khi .eq. 0) then
        m = jindex(jhi)
        return
      end if
      if (khi .lt. 0) return
c
  100 continue
      jxx = (jlo+jhi)/2
      kxx = jcodes(jxx) - icode
c
c ... did we find it ?
c
      if (kxx .eq. 0) then
        m = jindex(jxx)
        return
      end if
c
c ... replace the bracket with the same sign
c
      if (kxx .gt. 0) then
        jhi = jxx
        khi = kxx
      else
        jlo = jxx
        klo = kxx
      end if
c
      if (abs(jhi-jlo) .le. 1) return
c
      goto 100
c
      end
