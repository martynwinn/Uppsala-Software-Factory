c
c ========================================================================
c
      subroutine jacobi (a,n,np,d,v,nrot,ierr)
c
c ... calc eigenvalues and eigenvectors for real symmetric N*N matrix A
c     D returns eigenvalues, V contains column eigen vectors
c
c ... Numerical Recipes, pp. 346
c
      implicit none
c
      integer nmax,itmax
      parameter (nmax=1000, itmax=100)
c
      integer n,np,nrot,ierr
c
      real a(np,np),d(np),v(np,np),b(nmax),z(nmax)
      real sm,g,h,s,t,c,tau,theta,tresh
c
      integer ip,iq,i,j
c
code ...
c
      ierr = -1
c
      if (n .gt. nmax) then
        call errcon ('Increase dimensioning in routine JACOBI')
        return
      end if
c
      do ip=1,n
        do iq=1,n
          v(ip,iq) = 0.
        end do
        v(ip,ip) = 1.
      end do
c
      do ip=1,n
        b(ip) = a(ip,ip)
        d(ip) = b(ip)
        z(ip) = 0.
      end do
c
      nrot = 0
c
      do i=1,itmax
        sm = 0.0
        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
          end do
        end do
ccc        print *,' Iter # ',i,' Sum ',sm
        if (sm .eq. 0.0) then
          ierr = 0
          return
        end if
        if (i.lt.4) then
          tresh = 0.2*sm/(n**2)
        else
          tresh = 0.
        end if
c
        do ip=1,n-1
          do iq=ip+1,n
            g=100.0*abs(a(ip,iq))
            if ( (i.gt.4) .and. (abs(d(ip))+g.eq.abs(d(ip))) .and.
     +           (abs(d(iq))+g.eq.abs(d(iq))) ) then
              a(ip,iq) = 0.0
            else if (abs(a(ip,iq)) .gt. tresh) then
              h=d(iq)-d(ip)
              if (abs(h)+g.eq.abs(h)) then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1.0/(abs(theta)+sqrt(1.0+theta**2))
                if (theta.lt.0.) t=-t
              end if
              c=1.0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq) = 0.0
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              end do
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              end do
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              end do
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              end do
              nrot = nrot + 1
            end if
          end do
        end do
        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.0
        end do
      end do
c
      call errcon ('Too many iterations in JACOBI')
c
      return
      end
