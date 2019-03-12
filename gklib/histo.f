c
c
c
      subroutine histo (ndata,xdata,nr,rh,nh)
c
c ... list a histogram; real data array XDATA(NDATA)
c     nr = nr of histogram limits
c     RH(NR) = values of these limits
c     NH(NR+1) = counts (note dimensioning + 1 !!!)
c
c ... Gerard Kleywegt, Uppsala, 1993
c
      implicit NONE
c
      real xdata(*),rh(*),perc,xn,cum
c
      integer ndata,nr,nh(*),i,j
c
code ...
c
      do i=1,nr+1
        nh (i) = 0
      end do
c
      do i=1,ndata
        if (xdata(i) .lt. rh(1)) then
          nh (1) = nh (1) + 1
        else if (xdata(i) .ge. rh(nr)) then
          nh (nr+1) = nh (nr+1) + 1
        else
          do j=1,nr-1
            if (xdata(i) .ge. rh(j) .and.
     +          xdata(i) .lt. rh(j+1)) then
              nh (j+1) = nh (j+1) + 1
              goto 113
            end if
          end do
  113     continue
        end if
      end do
c
      xn = 100.0/float(ndata)
      perc = nh(1)*xn
      cum = perc
      write (*,7010) rh(1),nh(1),perc,cum
c
      do i=1,nr-1
        perc = nh(i+1)*xn
        cum = cum + perc
        if (nh(i+1) .gt. 0)
     +    write (*,7020) rh(i),rh(i+1),nh(i+1),perc,cum
      end do
c
      perc = nh(nr+1)*xn
      cum = cum + perc
      write (*,7030) rh(nr),nh(nr+1),perc,cum
c
 7010 format (/
     +  ' Nr <  ',f12.4,19x,' : ',i8,
     +  ' (',f6.2,' %; Cum ',f6.2,' %)')
 7020 format (
     +  ' Nr >= ',f12.4,' and < ',f12.4,' : ',i8,
     +  ' (',f6.2,' %; Cum ',f6.2,' %)')
 7030 format (
     +  ' Nr >= ',f12.4,19x,' : ',i8,
     +  ' (',f6.2,' %; Cum ',f6.2,' %)'/)
c
      return
      end
