c
c ... xrave.f - RAVE-related routines
c
c ... Gerard Kleywegt & Alwyn Jones
c
c
c
      subroutine cntnbr (mask,shadow,na,nb,nc,n1,n2)
c
c ... cnt nr of nbrs which are in mask for non-mask points
c     (n1=0, n2=1) or the nr of nbrs outside the mask for
c     mask points (n1=1, n2=0); NBRS = [0,26]
c
      implicit none
c
      integer na,nb,nc,i1,i2,i3,ncnt,n1,n2
      integer mask(na,nb,nc),shadow(na,nb,nc)
c
code ...
c
      call inimsk (shadow,na,nb,nc,0)
c
      do i1=2,na-1
        do i2=2,nb-1
          do i3=2,nc-1
            ncnt = 0
            if (mask(i1,i2,i3) .eq. n1) then
c
              if (mask(i1-1,i2-1,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2-1,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2-1,i3+1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2  ,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2  ,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2  ,i3+1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2+1,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2+1,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1-1,i2+1,i3+1) .eq. n2) ncnt = ncnt + 1
c
              if (mask(i1  ,i2-1,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1  ,i2-1,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1  ,i2-1,i3+1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1  ,i2  ,i3-1) .eq. n2) ncnt = ncnt + 1
c
              if (mask(i1  ,i2  ,i3+1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1  ,i2+1,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1  ,i2+1,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1  ,i2+1,i3+1) .eq. n2) ncnt = ncnt + 1
c
              if (mask(i1+1,i2-1,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2-1,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2-1,i3+1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2  ,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2  ,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2  ,i3+1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2+1,i3-1) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2+1,i3  ) .eq. n2) ncnt = ncnt + 1
              if (mask(i1+1,i2+1,i3+1) .eq. n2) ncnt = ncnt + 1
c
            end if
            shadow (i1,i2,i3) = ncnt
          end do
        end do
      end do
c
      return
      end
