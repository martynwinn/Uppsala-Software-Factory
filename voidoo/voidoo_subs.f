c
c ... cavity_subs.f - subroutines for VOIDOO
c
c ... Gerard Kleywegt
c
c ... NOTE - the old RECURSIVE zapping routines are still here
c            (RECZAP ... RE6ZAP), but they have been renamed to
c            XXRECZAP ... XXRE6ZAP.
c            These recursive routines work fine, but they lead
c            to crashes if the grid gets too big (on a 4D GTX
c            this happend if the grid contained more than about
c            half a million points). The numerical results are
c            the same with both sets of routines !!!
c
c ===========================================================================
c
      subroutine def_grd (gsize,xmin,xmax,ymin,ymax,zmin,zmax,
     +                    ngrid1,ngrid2,ngrid3)
c
c ... def_grd (...) - define a proper grid
c
      implicit none
c
      real gsize,xmin,xmax,ymin,ymax,zmin,zmax
c
      integer ngrid1,ngrid2,ngrid3
c
code ...
c
c ... NOTE: make the grid 2 points wider than strictly necessary !
c           We do this because then we can exclude the
c           border of the box from consideration in the
c           "zapping" routines which makes life a lot easier
c
      xmin  = gsize * float (int (xmin/gsize) - 1)
      xmax  = gsize * float (int (xmax/gsize) + 1)
      ymin  = gsize * float (int (ymin/gsize) - 1)
      ymax  = gsize * float (int (ymax/gsize) + 1)
      zmin  = gsize * float (int (zmin/gsize) - 1)
      zmax  = gsize * float (int (zmax/gsize) + 1)
c
      ngrid1 = 1 + nint ((xmax - xmin)/gsize)
      ngrid2 = 1 + nint ((ymax - ymin)/gsize)
      ngrid3 = 1 + nint ((zmax - zmin)/gsize)
c
      xmax = xmin + gsize * float (ngrid1-1)
      ymax = ymin + gsize * float (ngrid2-1)
      zmax = zmin + gsize * float (ngrid3-1)
c
      return
      end
c
c ===========================================================================
c
      subroutine set_up_grid (mode,tty,fnow,probe,natoms,ngrid,gsize,
     +    xmin,xmax,ymin,ymax,zmin,zmax,xg,yg,zg,nxyz,nxy,
     +    xd,yd,zd,dvdw,ibuff,ierror,prot,notp,cavi,ltrace)
c
c ... fill a grid
c
      implicit none
c
      real dx,dy,dz,fnow,total,user,sys,gsize,radius
      real xmin,xmax,ymin,ymax,zmin,zmax,xg(*),yg(*),zg(*)
      real xd(*),yd(*),zd(*),dvdw(*),probe,sqrrad
c
      integer atom,j0,k0,nhole,istart,jstart,kstart,nzap
      integer i1,i2,j1,j2,k1,k2,tty,natoms,ngrid(3)
      integer nxyz,nxy,ierror,i,j,k,mode
c
      logical ltrace
c
      character*1 ibuff(*),prot,notp,cavi
c
code ...
c
      if (ltrace) then
        write (tty,1200) ' ','Set up grid',
     +                       '-----------',' '
      else
        write (tty,1200) 'Setting up grid ...'
      end if
c
c ... Initialise grid
c
      if (ltrace) write (tty,1200) '(1) Initialise grid'
c
      do i=1,nxyz
        ibuff (i) = notp
      end do
c
      if (ltrace) then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)')
     +    '20 CPU total/user/sys :',total,user,sys
        write (tty,1200) '(2) Find empty space'
      end if
c
c ... Find empty space
c
c
cc      write (*,*) fnow,probe
cc      write (*,*) (ngrid(i),i=1,3)
cc      write (*,*) nxy,nxyz
c
      do atom=1,natoms
c
        radius = fnow*dvdw(atom) + probe
        sqrrad = radius**2
c
        i1 = int ( (xd(atom) - radius - xmin) / gsize ) + 1
        i1 = max (1, min (ngrid(1), i1))
        i2 = int ( (xd(atom) + radius - xmin) / gsize ) + 2
        i2 = max (1, min (ngrid(1), i2))
        j1 = int ( (yd(atom) - radius - ymin) / gsize ) + 1
        j1 = max (1, min (ngrid(2), j1))
        j2 = int ( (yd(atom) + radius - ymin) / gsize ) + 2
        j2 = max (1, min (ngrid(2), j2))
        k1 = int ( (zd(atom) - radius - zmin) / gsize ) + 1
        k1 = max (1, min (ngrid(3), k1))
        k2 = int ( (zd(atom) + radius - zmin) / gsize ) + 2
        k2 = max (1, min (ngrid(3), k2))
c
cc        write (*,'(7i10)') atom,i1,i2,j1,j2,k1,k2
c
        do i=i1,i2
          dx = (xg(i)-xd(atom))**2
          do j=j1,j2
            dy = dx + (yg(j)-yd(atom))**2
            j0 = (j-1)*ngrid(1) + i
            do k=k1,k2
              k0 = (k-1)*nxy + j0
              if (ibuff(k0) .eq. notp) then
                dz = dy + (zg(k)-zd(atom))**2
                if (dz .le. sqrrad) then
                  ibuff (k0) = prot
cc      write (*,'(a,3i4,i10,f10.2)') ' zapped ',i,j,k,k0,dz
                end if
              end if
            end do
          end do
        end do
      end do
c
      nhole = 0
      do i=1,nxyz
        if (ibuff(i) .eq. notp) nhole = nhole + 1
      end do
c
      call jvalut (' Nr of points in grid :',1,nxyz)
      call jvalut (' Not the protein      :',1,nhole)
      call jvalut (' The protein itself   :',1,nxyz-nhole)
c
c ... for cavity refinement, skip the rest of this routine
c
      if (mode .ne. 1) goto 300
c
      if (ltrace) then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '21 CPU total/user/sys :',
     +    total,user,sys
        write (tty,1200) '(3) Blank border'
      end if
c
c ... Blank border
c
      call set_border (0,' ',prot,ibuff,ngrid)
c
      nhole = 0
      do i=1,nxyz
        if (ibuff(i) .eq. notp) nhole = nhole + 1
      end do
c
      call jvalut (' Nr of non-protein points now :',1,nhole)
      if (nhole .lt. 1) then
cc        call errcon ('Is this a block-shaped protein ??!!')
cc        ierror = -1
        return
      end if
c
c ... Blank external space
c
      if (ltrace) then
        call gkdcpu (total,user,sys)
        write (tty,'(1x,a,3f10.1)') '22 CPU total/user/sys :',
     +    total,user,sys
        write (tty,1200) '(4) Blank external space'
      end if
c
c ... find first '1' (assumed to be external to the protein)
c
      do i=2,ngrid(1)-1
        do j=2,ngrid(2)-1
          j0 = (j-1)*ngrid(1) + i
          do k=2,ngrid(3)-1
            k0 = (k-1)*nxy + j0
            if (ibuff(k0) .eq. notp) then
              istart = i
              jstart = j
              kstart = k
              goto 100
            end if
          end do
        end do
      end do
c
      call errcon ('No external points found !!??')
      ierror = -2
      return
c
  100 continue
      if (ltrace) write (tty,'(1x,a,3i6)')
     +  'First point = ',istart,jstart,kstart
      nzap = 0
      call reczap (istart,jstart,kstart,nzap,notp,cavi)
c
      call jvalut (' Nr of points outside protein :',1,nzap)
c
      nhole = 0
      do i=1,nxyz
        if (ibuff(i) .eq. notp) nhole = nhole + 1
        if (ibuff(i) .eq. cavi) ibuff(i) = prot
      end do
c
      call jvalut (' Nr of cavity points          :',1,nhole)
      if (nhole .lt. 1) then
        call errcon ('This is a solid protein without cavities !')
cc        ierror = -1
        return
      end if
c
  300 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '23 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ','Set up grid done',
     +                     '----------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine reczap (istart,jstart,kstart,nzap,from,to)
c
c ... iterative zapping routine
c
      include 'cavity.incl'
c
      integer istart,jstart,kstart,inx,nzap
      integer ii,jj,kk,j0,k0,jj0,kk0,nchange
c
      character*1 from,to
c
code ...
c
      if (istart.le.1 .or. istart.ge.ngrid(1)) return
      if (jstart.le.1 .or. jstart.ge.ngrid(2)) return
      if (kstart.le.1 .or. kstart.ge.ngrid(3)) return
c
      inx = istart+(jstart-1)*ngrid(1)+(kstart-1)*nxy
      if (ibuff(inx) .eq. from) then
        ibuff (inx) = to
        nzap = nzap + 1
c
        nchange = 0
        do i=istart-1,istart+1
          if (i.le.1 .or. i.ge.ngrid(1)) goto 1234
          do j=jstart-1,jstart+1
            if (j.le.1 .or. j.ge.ngrid(2)) goto 1236
            j0 = (j-1)*ngrid(1) + i
            do k=kstart-1,kstart+1
              if (k.le.1 .or. k.ge.ngrid(3)) goto 1238
              k0 = (k-1)*nxy + j0
              if (ibuff (k0) .eq. from) then
                nchange = nchange + 1
                ibuff (k0) = to
              end if
 1238         continue
            end do
 1236       continue
          end do
 1234     continue
        end do
c
        nzap = nzap + nchange
        if (nchange .eq. 0) return
        ibuff (inx) = temp
c
c ... okay, the cavity contained more than one point; now sweep through
c     the matrix until nothing changes anymore
c
   10   continue
        nchange = 0
c
        do i=2,ngrid(1)-1
          do j=2,ngrid(2)-1
            j0 = (j-1)*ngrid(1) + i
            do k=2,ngrid(3)-1
              k0 = (k-1)*nxy + j0
              if (ibuff (k0) .eq. to) then
                do ii=i-1,i+1
                  do jj=j-1,j+1
                    jj0 = (jj-1)*ngrid(1) + ii
                    do kk=k-1,k+1
                      kk0 = (kk-1)*nxy + jj0
                      if (ibuff(kk0) .eq. from) then
                        nchange = nchange + 1
                        ibuff(kk0) = to
                      end if
                    end do
                  end do
                end do
                ibuff (k0) = temp
              end if
            end do
          end do
        end do
c
        nzap = nzap + nchange
        if (nchange .gt. 0) goto 10
c
        call set_border (1,to,from,ibuff,ngrid)
c
        do i=1,nxyz
          if (ibuff(i) .eq. temp) ibuff (i) = to
        end do
c
cc        print *,' *1* NZAP with border    = ',nzap
c
        nzap = 0
        do i=1,nxyz
          if (ibuff(i) .eq. to) nzap = nzap + 1
        end do
c
cc        print *,' *1* NZAP without border = ',nzap
c
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine re3zap (istart,jstart,kstart,nzap,from,to,
     +                   ilo,jlo,klo,ihi,jhi,khi)
c
c ... iterative zapping routine
c
      include 'cavity.incl'
c
      integer istart,jstart,kstart,inx,nzap,ilo,jlo,klo,ihi,jhi,khi
      integer ii,jj,kk,j0,k0,jj0,kk0,nchange
c
      character*1 from,to
c
code ...
c
      if (istart.le.1 .or. istart.ge.ngrid(1)) return
      if (jstart.le.1 .or. jstart.ge.ngrid(2)) return
      if (kstart.le.1 .or. kstart.ge.ngrid(3)) return
c
      inx = istart+(jstart-1)*ngrid(1)+(kstart-1)*nxy
      if (ibuff(inx) .eq. from) then
        ibuff (inx) = to
        nzap = nzap + 1
c
        nchange = 0
        do i=istart-1,istart+1
          if (i.le.1 .or. i.ge.ngrid(1)) goto 1234
          do j=jstart-1,jstart+1
            if (j.le.1 .or. j.ge.ngrid(2)) goto 1236
            j0 = (j-1)*ngrid(1) + i
            do k=kstart-1,kstart+1
              if (k.le.1 .or. k.ge.ngrid(3)) goto 1238
              k0 = (k-1)*nxy + j0
              if (ibuff (k0) .eq. from) then
                nchange = nchange + 1
                ibuff (k0) = to
                ilo = min(ilo,i)
                jlo = min(jlo,j)
                klo = min(klo,k)
                ihi = max(ihi,i)
                jhi = max(jhi,j)
                khi = max(khi,k)
              end if
 1238         continue
            end do
 1236       continue
          end do
 1234     continue
        end do
c
        nzap = nzap + nchange
        if (nchange .eq. 0) return
        ibuff (inx) = temp
c
c ... okay, the cavity contained more than one point; now sweep through
c     the matrix until nothing changes anymore
c
   10   continue
        nchange = 0
c
        do i=2,ngrid(1)-1
          do j=2,ngrid(2)-1
            j0 = (j-1)*ngrid(1) + i
            do k=2,ngrid(3)-1
              k0 = (k-1)*nxy + j0
              if (ibuff (k0) .eq. to) then
                do ii=i-1,i+1
                  do jj=j-1,j+1
                    jj0 = (jj-1)*ngrid(1) + ii
                    do kk=k-1,k+1
                      kk0 = (kk-1)*nxy + jj0
                      if (ibuff(kk0) .eq. from) then
                        nchange = nchange + 1
                        ibuff(kk0) = to
                      end if
                    end do
                  end do
                end do
                ibuff (k0) = temp
              end if
            end do
          end do
        end do
c
        nzap = nzap + nchange
        if (nchange .gt. 0) goto 10
c
        call set_border (1,to,from,ibuff,ngrid)
c
        do i=1,nxyz
          if (ibuff(i) .eq. temp) ibuff (i) = to
        end do
c
cc        print *,' *3* NZAP with border    = ',nzap
c
        nzap = 0
        do i=1,nxyz
          if (ibuff(i) .eq. to) nzap = nzap + 1
        end do
c
cc        print *,' *3* NZAP without border = ',nzap
c
c ... find lowest/highest indices
c
        do i=2,ngrid(1)-1
          do j=2,ngrid(2)-1
            j0 = (j-1)*ngrid(1) + i
            do k=2,ngrid(3)-1
              k0 = (k-1)*nxy + j0
              if (ibuff (k0) .eq. to) then
                ilo = min(ilo,i)
                jlo = min(jlo,j)
                klo = min(klo,k)
                ihi = max(ihi,i)
                jhi = max(jhi,j)
                khi = max(khi,k)
              end if
            end do
          end do
        end do
c
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine re5zap (istart,jstart,kstart,nzap,from,to,
     +                   mgrid1,mgrid2,mgrid3,mxy)
c
c ... iterative zapping routine
c
      include 'cavity.incl'
c
      integer istart,jstart,kstart,inx,nzap,mgrid1,mgrid2,mgrid3,mxy
      integer ii,jj,kk,j0,k0,jj0,kk0,nchange,igrid(3)
c
      character*1 from,to
c
code ...
c
      if (istart.le.1 .or. istart.ge.mgrid1) return
      if (jstart.le.1 .or. jstart.ge.mgrid2) return
      if (kstart.le.1 .or. kstart.ge.mgrid3) return
c
      inx = istart+(jstart-1)*mgrid1+(kstart-1)*mxy
      if (ibuff2 (inx) .eq. from) then
        ibuff2 (inx) = to
        nzap = nzap + 1
c
        nchange = 0
        do i=istart-1,istart+1
          if (i.le.1 .or. i.ge.mgrid1) goto 1234
          do j=jstart-1,jstart+1
            if (j.le.1 .or. j.ge.mgrid2) goto 1236
            j0 = (j-1)*mgrid1 + i
            do k=kstart-1,kstart+1
              if (k.le.1 .or. k.ge.mgrid3) goto 1238
              k0 = (k-1)*mxy + j0
              if (ibuff2 (k0) .eq. from) then
                nchange = nchange + 1
                ibuff2 (k0) = to
              end if
 1238         continue
            end do
 1236       continue
          end do
 1234     continue
        end do
c
        nzap = nzap + nchange
        if (nchange .eq. 0) return
        ibuff2 (inx) = temp
c
c ... okay, the cavity contained more than one point; now sweep through
c     the matrix until nothing changes anymore
c
   10   continue
        nchange = 0
c
        do i=2,mgrid1-1
          do j=2,mgrid2-1
            j0 = (j-1)*mgrid1 + i
            do k=2,mgrid3-1
              k0 = (k-1)*mxy + j0
              if (ibuff2 (k0) .eq. to) then
                do ii=i-1,i+1
                  do jj=j-1,j+1
                    jj0 = (jj-1)*mgrid1 + ii
                    do kk=k-1,k+1
                      kk0 = (kk-1)*mxy + jj0
                      if (ibuff2(kk0) .eq. from) then
                        nchange = nchange + 1
                        ibuff2 (kk0) = to
                      end if
                    end do
                  end do
                end do
                ibuff2 (k0) = temp
              end if
            end do
          end do
        end do
c
        nzap = nzap + nchange
        if (nchange .gt. 0) goto 10
c
        igrid (1) = mgrid1
        igrid (2) = mgrid2
        igrid (3) = mgrid3
c
        call set_border (1,to,from,ibuff2,igrid)
c
        do i=1,mxy*mgrid3
          if (ibuff2(i) .eq. temp) ibuff2 (i) = to
        end do
c
cc        print *,' *5* NZAP with border    = ',nzap
c
        nzap = 0
        do i=1,mxy*mgrid3
          if (ibuff2(i) .eq. to) nzap = nzap + 1
        end do
c
cc        print *,' *5* NZAP without border = ',nzap
c
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine xxreczap (istart,jstart,kstart,nzap,from,to)
c
c ... recursive zapping routine (NEEDS RE2ZAP !!!)
c
      include 'cavity.incl'
c
      integer istart,jstart,kstart,inx,nzap
c
      character*1 from,to
c
code ...
c
      if (istart.le.1 .or. istart.ge.ngrid(1)) return
      if (jstart.le.1 .or. jstart.ge.ngrid(2)) return
      if (kstart.le.1 .or. kstart.ge.ngrid(3)) return
c
      inx = istart+(jstart-1)*ngrid(1)+(kstart-1)*nxy
      if (ibuff(inx) .eq. from) then
        ibuff (inx) = to
        nzap = nzap + 1
        do i=istart-1,istart+1
          do j=jstart-1,jstart+1
            do k=kstart-1,kstart+1
              call xxre2zap (i,j,k,nzap,from,to)
            end do
          end do
        end do
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine xxre2zap (i,j,k,nzap,from,to)
c
c ... recursive zapping routine (NEEDED FOR RECZAP !!!)
c
      include 'cavity.incl'
c
      integer nzap
c
      character*1 from,to
c
code ...
c
      call xxreczap (i,j,k,nzap,from,to)
c
      return
      end
c
c ===========================================================================
c
      subroutine xxre3zap (istart,jstart,kstart,nzap,from,to,
     +                   ilo,jlo,klo,ihi,jhi,khi)
c
c ... recursive zapping routine (NEEDS RE4ZAP !!!)
c
      include 'cavity.incl'
c
      integer istart,jstart,kstart,inx,nzap,ilo,jlo,klo,ihi,jhi,khi
c
      character*1 from,to
c
code ...
c
      if (istart.le.1 .or. istart.ge.ngrid(1)) return
      if (jstart.le.1 .or. jstart.ge.ngrid(2)) return
      if (kstart.le.1 .or. kstart.ge.ngrid(3)) return
c
      inx = istart+(jstart-1)*ngrid(1)+(kstart-1)*nxy
      if (ibuff(inx) .eq. from) then
        ibuff (inx) = to
        nzap = nzap + 1
c
        ilo = min(ilo,istart)
        jlo = min(jlo,jstart)
        klo = min(klo,kstart)
        ihi = max(ihi,istart)
        jhi = max(jhi,jstart)
        khi = max(khi,kstart)
c
        do i=istart-1,istart+1
          do j=jstart-1,jstart+1
            do k=kstart-1,kstart+1
              call xxre4zap (
     +          i,j,k,nzap,from,to,ilo,jlo,klo,ihi,jhi,khi)
            end do
          end do
        end do
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine xxre4zap (i,j,k,nzap,from,to,ilo,jlo,klo,ihi,jhi,khi)
c
c ... recursive zapping routine (NEEDED FOR RE3ZAP !!!)
c
      include 'cavity.incl'
c
      integer nzap,ilo,jlo,klo,ihi,jhi,khi
c
      character*1 from,to
c
code ...
c
      call xxre3zap (i,j,k,nzap,from,to,ilo,jlo,klo,ihi,jhi,khi)
c
      return
      end
c
c ===========================================================================
c
      subroutine xxre5zap (istart,jstart,kstart,nzap,from,to,
     +                   mgrid1,mgrid2,mgrid3,mxy)
c
c ... recursive zapping routine (NEEDS RE6ZAP !!!)
c
      include 'cavity.incl'
c
      integer istart,jstart,kstart,inx,nzap,mgrid1,mgrid2,mgrid3,mxy
c
      character*1 from,to
c
code ...
c
      if (istart.le.1 .or. istart.ge.mgrid1) return
      if (jstart.le.1 .or. jstart.ge.mgrid2) return
      if (kstart.le.1 .or. kstart.ge.mgrid3) return
c
      inx = istart+(jstart-1)*mgrid1+(kstart-1)*mxy
      if (ibuff2 (inx) .eq. from) then
        ibuff2 (inx) = to
        nzap = nzap + 1
c
        do i=istart-1,istart+1
          do j=jstart-1,jstart+1
            do k=kstart-1,kstart+1
              call xxre6zap (
     +          i,j,k,nzap,from,to,mgrid1,mgrid2,mgrid3,mxy)
            end do
          end do
        end do
      end if
c
      return
      end
c
c ===========================================================================
c
      subroutine xxre6zap (i,j,k,nzap,from,to,
     +                       mgrid1,mgrid2,mgrid3,mxy)
c
c ... recursive zapping routine (NEEDED FOR RE5ZAP !!!)
c
      include 'cavity.incl'
c
      integer nzap,mgrid1,mgrid2,mgrid3,mxy
c
      character*1 from,to
c
code ...
c
      call xxre5zap (i,j,k,nzap,from,to,mgrid1,mgrid2,mgrid3,mxy)
c
      return
      end
c
c ===========================================================================
c
      subroutine set_border (mode,from,to,ibuff,ngrid)
c
c ... set border voxels
c     if mode = 0, then set all of them to value 'to'
c     if mode = 1, then only change those that are 'from' to 'to'
c
      implicit none
c
      integer mode,ngrid(3),i,j,k,j0,k0,nxy
c
      character*1 from,to,ibuff(*)
c
code ...
c
      nxy = ngrid(1)*ngrid(2)
c
      if (mode .eq. 1) goto 100
c
      i=1
      do j=1,ngrid(2)
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          ibuff (k0) = to
        end do
      end do
c
      i=ngrid(1)
      do j=1,ngrid(2)
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          ibuff (k0) = to
        end do
      end do
c
      do i=1,ngrid(1)
        j=1
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          ibuff (k0) = to
        end do
c
        j=ngrid(2)
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          ibuff (k0) = to
        end do
      end do
c
      do i=1,ngrid(1)
        do j=1,ngrid(2)
          j0 = (j-1)*ngrid(1) + i
          k=1
          k0 = (k-1)*nxy + j0
          ibuff (k0) = to
c
          k=ngrid(3)
          k0 = (k-1)*nxy + j0
          ibuff (k0) = to
        end do
      end do
c
      return
c
  100 continue
c
      i=1
      do j=1,ngrid(2)
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          if (ibuff(k0).eq.from) ibuff (k0) = to
        end do
      end do
c
      i=ngrid(1)
      do j=1,ngrid(2)
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          if (ibuff(k0).eq.from) ibuff (k0) = to
        end do
      end do
c
      do i=1,ngrid(1)
        j=1
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          if (ibuff(k0).eq.from) ibuff (k0) = to
        end do
c
        j=ngrid(2)
        j0 = (j-1)*ngrid(1) + i
        do k=1,ngrid(3)
          k0 = (k-1)*nxy + j0
          if (ibuff(k0).eq.from) ibuff (k0) = to
        end do
      end do
c
      do i=1,ngrid(1)
        do j=1,ngrid(2)
          j0 = (j-1)*ngrid(1) + i
          k=1
          k0 = (k-1)*nxy + j0
          if (ibuff(k0).eq.from) ibuff (k0) = to
c
          k=ngrid(3)
          k0 = (k-1)*nxy + j0
          if (ibuff(k0).eq.from) ibuff (k0) = to
        end do
      end do
c
      return
      end
c
c ===========================================================================
c
      logical function conver (oldval,newval,vtol,ptol,vdif,pdif)
c
c ... determine if convergence has occurred
c
      implicit none
c
      real oldval,newval,vtol,ptol,vdif,pdif
c
code ...
c
      vdif = abs (oldval - newval)
      pdif = 99999.999
      if (oldval .ne. 0.0) pdif = abs (100.0 * vdif / oldval)
c
      conver = (vdif .le. vtol .or. pdif .le. ptol)
c
      return
      end
c
c ===========================================================================
c
      subroutine surf_layer (buff,nx,ny,cavi,surf,prot,notp,np)
c
      implicit none
c
      integer nx,ny,i,j,k,ii,jj,kk,np
c
      character*1 cavi,surf,prot,notp,buff(*)
c
code ...
c
      do i=2,nx-1
        do j=2,ny-1
          k = (j-1)*nx + i
          if (buff(k) .eq. cavi) then
            do ii=i-1,i+1
              do jj=j-1,j+1
                kk = (jj-1)*nx + ii
                if (buff(kk) .eq. prot) then
                  buff (k) = surf
                  goto 738
                else if (buff(kk) .eq. notp) then
                  buff (k) = surf
                  goto 738
                end if
              end do
            end do
          end if
  738     continue
        end do
      end do
c
      np = 0
      do i=1,nx*ny
        if (buff(i) .eq. surf) np = np + 1
      end do
c
      return
      end
c
c ===========================================================================
c
      subroutine cont_layer (buff,nx,ny,mode,zoff,xco,yco,unit,
     +                       surf,done,temp)
c
c ... mode = 1 => Y/Z data
c     mode = 2 => X/Z data
c     mode = 3 => X/Y data
c
      implicit none
c
      real xco(*),yco(*),zoff
c
      integer nx,ny,i,j,k,mode,unit,length,if,jf,ig,jg,leng1
c
      character*1 surf,done,temp,buff(*)
      character line*80
c
code ...
c
c ... find a starting point labelled 'surf'
c
      do i=2,nx-1
        do j=2,ny-1
          k = (j-1)*nx + i
          if (buff(k) .eq. surf) then
c
c ... found the start of a new contour !
c
            if (mode .eq. 1) then
              write (line,'(a,3f10.3)') 'm ',zoff,xco(i),yco(j)
            else if (mode .eq. 2) then
              write (line,'(a,3f10.3)') 'm ',xco(i),zoff,yco(j)
            else if (mode .eq. 3) then
              write (line,'(a,3f10.3)') 'm ',xco(i),yco(j),zoff
            end if
            call pretty (line)
            write (unit,'(a)') line(1:leng1(line))
c
c ... now label starting point and trace contour in a
c     counter-clockwise fashion
c
            buff (k) = temp
            if = i
            jf = j
c
   10       continue
c
c ... (1) six o'clock
c
            ig = if
            jg = jf - 1
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (2) half past four
c
            ig = if + 1
            jg = jf - 1
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (3) three o'clock
c
            ig = if + 1
            jg = jf
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (4) half past one
c
            ig = if + 1
            jg = jf + 1
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (5) twelve o'clock
c
            ig = if
            jg = jf + 1
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (6) half past ten
c
            ig = if - 1
            jg = jf + 1
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (7) nine o'clock
c
            ig = if - 1
            jg = jf
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (8) half past seven
c
            ig = if - 1
            jg = jf - 1
            k = (jg-1)*nx + ig
            if (buff(k) .eq. temp) goto 20
            if (buff(k) .eq. surf) goto 30
c
c ... (9) no connections, just an isolated point
c
            ig = if
            jg = jf
            k = (jg-1)*nx + ig
            goto 20
c
c ... arrived at next point of contour
c
   30       continue
c
            buff (k) = done
c
            if (mode .eq. 1) then
              write (line,'(a,3f10.3)') 'l ',zoff,xco(ig),yco(jg)
            else if (mode .eq. 2) then
              write (line,'(a,3f10.3)') 'l ',xco(ig),zoff,yco(jg)
            else if (mode .eq. 3) then
              write (line,'(a,3f10.3)') 'l ',xco(ig),yco(jg),zoff
            end if
            call pretty (line)
            write (unit,'(a)') line(1:leng1(line))
c
            if = ig
            jf = jg
c
            goto 10
c
c ... arrived at start point of contour
c
   20       continue
c
            buff (k) = done
c
            if (mode .eq. 1) then
              write (line,'(a,3f10.3)') 'l ',zoff,xco(ig),yco(jg)
            else if (mode .eq. 2) then
              write (line,'(a,3f10.3)') 'l ',xco(ig),zoff,yco(jg)
            else if (mode .eq. 3) then
              write (line,'(a,3f10.3)') 'l ',xco(ig),yco(jg),zoff
            end if
            call pretty (line)
            write (unit,'(a)') line(1:leng1(line))
c
          end if
        end do
      end do
c
      return
      end
c
c ===========================================================================
c
      subroutine odl_con_surf (mgrid,xh,yh,zh)
c
      include 'cavity.incl'
c
      real xh(*),yh(*),zh(*)
c
      integer mgrid(3),mxy,j0,k0,ii,jj,kk,jj0,kk0
      integer mxyz,nhole,nhol2,length,leng1
c
      character flagit*1
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Connecting all surface points',
     +  '-----------------------------',' '
c
      mxy  = mgrid(1) * mgrid(2)
      mxyz = mxy      * mgrid(3)
c
c ... zap the inside (i.e., all points that are surrounded completely
c     by other cavity points); in this implementation, set a grid point
c     from 'cavi' to 'surf' if it has at least one protein or
c     different-cavity neighbour
c
      do i=2,mgrid(1)-1
        do j=2,mgrid(2)-1
          j0 = (j-1)*mgrid(1) + i
          do k=2,mgrid(3)-1
            k0 = (k-1)*mxy + j0
c
            if (ibuff2(k0) .eq. cavi) then
              do ii=i-1,i+1
                do jj=j-1,j+1
                  jj0 = (jj-1)*mgrid(1) + ii
                  do kk=k-1,k+1
                    kk0 = (kk-1)*mxy + jj0
c
                    if (ibuff2(kk0) .eq. prot) then
                      ibuff2 (k0) = surf
                      goto 738
                    end if
c
                    if (ibuff2(kk0) .eq. notp) then
                      ibuff2 (k0) = surf
                      goto 738
                    end if
c
                  end do
                end do
              end do
            end if
  738       continue
          end do
        end do
      end do
c
cc      write (line,'(a,i,a,i,a,i,a)')
cc     +  '(',mgrid(3),'(/1x,',mgrid(2),'(/',mgrid(1),'a1)))'
cc      call remspa (line)
cc      call textut (' format :',line)
cc      write (tty,fmt=line) (ibuff2(i),i=1,mxyz)
c
      nhole = 0
      nhol2 = 0
      do i=1,mxyz
        if (ibuff2(i) .eq. surf) nhole = nhole + 1
        if (ibuff2(i) .eq. cavi) nhol2 = nhol2 + 1
      end do
      call jvalut (' Nr of points on "cavity surface" :',1,nhole)
      call jvalut (' Nr of points inside cavity       :',1,nhol2)
      flagit = surf
      if (nhole .eq. 0) flagit = cavi
c
cc      call gkdcpu (total,user,sys)
cc      write (tty,'(1x,a,3f10.1)') '18 CPU total/user/sys :',
cc     +  total,user,sys
c
c ... now draw lines from all points marked 'flagit' to their neighbours
c     which are also marked 'flagit', then mark the point '4' to avoid
c     duplication of line elements
c
      do i=2,mgrid(1)-1
        do j=2,mgrid(2)-1
          j0 = (j-1)*mgrid(1) + i
          do k=2,mgrid(3)-1
            k0 = (k-1)*mxy + j0
c
            if (ibuff2(k0) .eq. flagit) then
              do ii=i-1,i+1
                do jj=j-1,j+1
                  jj0 = (jj-1)*mgrid(1) + ii
                  do kk=k-1,k+1
                    kk0 = (kk-1)*mxy + jj0
c
                    if (ii.eq.i) then
                      if (jj.eq.j) then
                        if (kk.eq.k) goto 736
                      end if
                    end if
c
                    if (ibuff2(kk0) .eq. flagit) then
                      write (line,'(a,3f10.3)')
     +  'm ',xh(i),yh(j),zh(k)
                      call pretty (line)
                      write (junit,'(a)') line(1:leng1(line))
                      write (line,'(a,3f10.3)')
     +  'l ',xh(ii),yh(jj),zh(kk)
                      call pretty (line)
                      write (junit,'(a)') line(1:leng1(line))
                    end if
c
  736               continue
                  end do
                end do
              end do
              ibuff2 (k0) = done
            end if
c
          end do
        end do
      end do
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '19 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Connecting all surface points done',
     +  '----------------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine odl_1sw_cont (mgrid,xh,yh,zh)
c
      include 'cavity.incl'
c
      real xh(*),yh(*),zh(*)
c
      integer mgrid(3),mxy,j0,k0,kk0,np
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Creating 1-sweep contour',
     +  '------------------------',' '
c
      mxy  = mgrid(1) * mgrid(2)
c
      do i=2,mgrid(1)-1
c
c ... loop over (j,k) layers; first extract 2D layer
c
        do j=1,mgrid(2)
          j0 = (j-1)*mgrid(1) + i
          do k=1,mgrid(3)
            k0  = (k-1)*mxy + j0
            kk0 = (k-1)*mgrid(2) + j
            ibuff3 (kk0) = ibuff2 (k0)
          end do
        end do
c
c ... mark surface points 'surf'
c
        call surf_layer (
     +    ibuff3,mgrid(2),mgrid(3),cavi,surf,prot,notp,np)
c
c ... trace the contour(s)
c
        if (np .gt. 0)
     +    call cont_layer (
     +      ibuff3,mgrid(2),mgrid(3),1,xh(i),yh,zh,junit,
     +      surf,done,temp)
c
      end do
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '24 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Creating 1-sweep contour done',
     +  '-----------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine odl_3sw_cont (mgrid,xh,yh,zh)
c
      include 'cavity.incl'
c
      real xh(*),yh(*),zh(*)
c
      integer mgrid(3),mxy,j0,k0,kk0,np
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Creating 3-sweep contour',
     +  '------------------------',' '
c
      mxy  = mgrid(1) * mgrid(2)
c
c ... (1) Y/Z layers
c
      do i=2,mgrid(1)-1
c
c ... loop over (j,k) layers; first extract 2D layer
c
        do j=1,mgrid(2)
          j0 = (j-1)*mgrid(1) + i
          do k=1,mgrid(3)
            k0  = (k-1)*mxy + j0
            kk0 = (k-1)*mgrid(2) + j
            ibuff3 (kk0) = ibuff2 (k0)
          end do
        end do
c
c ... mark surface points 'surf'
c
        call surf_layer (
     +    ibuff3,mgrid(2),mgrid(3),cavi,surf,prot,notp,np)
c
c ... trace the contour(s)
c
        if (np .gt. 0)
     +    call cont_layer (
     +      ibuff3,mgrid(2),mgrid(3),1,xh(i),yh,zh,junit,
     +      surf,done,temp)
c
      end do
c
c ... (2) X/Z layers
c
      do j=2,mgrid(2)-1
c
c ... loop over (i,k) layers; first extract 2D layer
c
        do i=1,mgrid(1)
          j0 = (j-1)*mgrid(1) + i
          do k=1,mgrid(3)
            k0  = (k-1)*mxy + j0
            kk0 = (k-1)*mgrid(1) + i
            ibuff3 (kk0) = ibuff2 (k0)
          end do
        end do
c
c ... mark surface points 'surf'
c
        call surf_layer (
     +    ibuff3,mgrid(1),mgrid(3),cavi,surf,prot,notp,np)
c
c ... trace the contour(s)
c
        if (np .gt. 0)
     +    call cont_layer (
     +      ibuff3,mgrid(1),mgrid(3),2,yh(j),xh,zh,junit,
     +      surf,done,temp)
c
      end do
c
c ... (3) X/Y layers
c
      do k=2,mgrid(3)-1
c
c ... loop over (i,j) layers; first extract 2D layer
c
        do i=1,mgrid(1)
          j0 = (k-1)*mxy +i 
          do j=1,mgrid(2)
            k0  = (j-1)*mgrid(1) + j0
            kk0 = (j-1)*mgrid(1) + i
            ibuff3 (kk0) = ibuff2 (k0)
          end do
        end do
c
c ... mark surface points 'surf'
c
        call surf_layer (
     +    ibuff3,mgrid(1),mgrid(2),cavi,surf,prot,notp,np)
c
c ... trace the contour(s)
c
        if (np .gt. 0)
     +    call cont_layer (
     +      ibuff3,mgrid(1),mgrid(2),3,zh(k),xh,yh,junit,
     +      surf,done,temp)
c
      end do
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '25 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Creating 3-sweep contour done',
     +  '-----------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine odl_dots (mgrid,xh,yh,zh)
c
      include 'cavity.incl'
c
      real xh(*),yh(*),zh(*)
c
      integer mgrid(3),mxy,j0,k0,mxyz,length,leng1
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Dotting all cavity points',
     +  '-------------------------',' '
c
      mxy  = mgrid(1) * mgrid(2)
      mxyz = mxy      * mgrid(3)
c
c ... now draw dots for all cavity points
c
      do i=2,mgrid(1)-1
        do j=2,mgrid(2)-1
          j0 = (j-1)*mgrid(1) + i
          do k=2,mgrid(3)-1
            k0 = (k-1)*mxy + j0
c
            if (ibuff2(k0) .eq. cavi) then
              write (line,'(a,3f10.3)') 'dot ',xh(i),yh(j),zh(k)
              call pretty (line)
              write (junit,'(a)') line(1:leng1(line))
            end if
c
          end do
        end do
      end do
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '26 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Dotting all cavity points done',
     +  '------------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine prot_surf
c
      include 'cavity.incl'
c
      real xh(maxgrd),yh(maxgrd),zh(maxgrd),xlo,xhi,ylo,yhi,zlo,zhi
c
      integer j0,k0,length,mgrid(3),mxy,mxyz,ii,jj,kk,jj0,kk0
      integer nsurf,nprot,nscale,leng1
c
      character myline*120
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Plotting all protein surface points',
     +  '-----------------------------------',' '
c
c      xlo = xmin - probe
c      xhi = xmax + probe
c      ylo = ymin - probe
c      yhi = ymax + probe
c      zlo = zmin - probe
c      zhi = zmax + probe
c
      xlo = xmin
      xhi = xmax
      ylo = ymin
      yhi = ymax
      zlo = zmin
      zhi = zmax
c
      call def_grd (dgrid,xlo,xhi,ylo,yhi,zlo,zhi,
     +              mgrid(1),mgrid(2),mgrid(3))
c
      call jvalut (' Protein plot grid :',3,mgrid)
c
      mxy  = mgrid(1)*mgrid(2)
      mxyz = mgrid(1)*mgrid(2)*mgrid(3)
c
      if (mxyz     .gt. maxbuf .or.
     +    mgrid(1) .gt. maxgrd .or.
     +    mgrid(2) .gt. maxgrd .or.
     +    mgrid(3) .gt. maxgrd ) then
c
        call errcon ('Grid too big for buffer')
        ierror = -1
        return
      end if
c
      do i=1,mgrid(1)
        xh (i) = xlo + float(i-1)*dgrid
      end do
      do i=1,mgrid(2)
        yh (i) = ylo + float(i-1)*dgrid
      end do
      do i=1,mgrid(3)
        zh (i) = zlo + float(i-1)*dgrid
      end do
c
      call set_up_grid (0,tty,1.0,probe,natoms,mgrid,dgrid,
     +  xlo,xhi,ylo,yhi,zlo,zhi,xh,yh,zh,mxyz,mxy,
     +  xd,yd,zd,dvdw,ibuff2,ierror,prot,notp,cavi,ltrace)
c
      nprot = 0
      do i=1,mxyz
        if (ibuff2(i) .eq. prot) nprot = nprot + 1
      end do
c
      call jvalut (' Nr of protein points :',1,nprot)
c
c ... DOT ODL FILE
c
      if (apdot .eq. 'D') then
c
        do i=2,mgrid(1)-1
          do j=2,mgrid(2)-1
            j0 = (j-1)*mgrid(1) + i
            do k=2,mgrid(3)-1
              k0 = (k-1)*mxy + j0
c
              if (ibuff2(k0) .eq. prot) then
                do ii=i-1,i+1
                  do jj=j-1,j+1
                    jj0 = (jj-1)*mgrid(1) + ii
                    do kk=k-1,k+1
                      kk0 = (kk-1)*mxy + jj0
c
                      if (ibuff2(kk0) .eq. notp) then
                        ibuff2 (k0) = surf
                        goto 738
                      end if
c
                    end do
                  end do
                end do
              end if
  738         continue
            end do
          end do
        end do
c
        nsurf = 0
        do i=1,mxyz
          if (ibuff2(i) .eq. surf) nsurf = nsurf + 1
        end do
c
        call jvalut (' Nr of protein-surface points :',1,nsurf)
c
        if (nsurf .le. 0) goto 999
c
        close (junit)
        call xopxua (junit,dfile,xinter(),ierror)
        if (ierror .ne. 0) goto 999
c
        write (junit,'(a)') 'begin dot_prot'
        write (junit,'(a)') 'colour 16799999'
c
        do i=2,mgrid(1)-1
          do j=2,mgrid(2)-1
            j0 = (j-1)*mgrid(1) + i
            do k=2,mgrid(3)-1
              k0 = (k-1)*mxy + j0
c
              if (ibuff2(k0) .eq. surf) then
                write (line,'(a,3f10.3)') 'dot ',xh(i),yh(j),zh(k)
                call pretty (line)
                write (junit,'(a)') line(1:leng1(line))
              end if
c
            end do
          end do
        end do
c
        write (junit,'(a)') 'end'
        close (junit)
c
c ... (NEW-) EZD FILE
c
      else
c
        do i=1,mxyz
          if (ibuff2(i) .eq. prot) ibuff2(i) = temp
        end do
c
        do i=1,mxyz
          if (ibuff2(i) .ne. temp) ibuff2(i) = '0'
        end do
c
        do i=1,mxyz
          if (ibuff2(i) .eq. temp) ibuff2(i) = '1'
        end do
c
        nsurf = 0
        do i=1,mxyz
          if (ibuff2(i) .eq. '1') nsurf = nsurf + 1
        end do
c
        call jvalut (' Nr of protein points :',1,nsurf)
c
        close (junit)
        call xopxua (junit,dfile,xinter(),ierror)
        if (ierror .ne. 0) goto 999
c
        nscale = 100
 6511   continue
        if (max(mgrid(1),mgrid(2),mgrid(3)) .gt. (nscale-20)) then
          nscale = nscale * 2
          goto 6511
        end if
c
        if (apdot .eq. 'O') then
c
          write (junit,7000) nint(xh(1)/dgrid),nint(yh(1)/dgrid),
     +      nint(zh(1)/dgrid)
          write (junit,7000) mgrid(1),mgrid(2),mgrid(3)
          write (junit,7000) nscale,nscale,nscale
          write (junit,7010) float(nscale)*dgrid,float(nscale)*dgrid,
     +      float(nscale)*dgrid,90.0,90.0,90.0
          write (junit,7020) 'MAP','(1X,75F1.0)',1.0
          write (junit,7030) (ibuff2(i),i=1,mxyz)
c
        else
c
          write (junit,'(a)') 'EZD_MAP'
c
          call stamp (myline)
          myline = '! '//myline
          write (junit,'(a)') myline(1:leng1(myline))
c
          write (myline,8010) 'CELL',float(nscale)*dgrid,
     +      float(nscale)*dgrid,float(nscale)*dgrid,90.0,90.0,90.0
          call pretty (myline)
          write (junit,'(a)') myline(1:leng1(myline))
c
          write (myline,8000) 'ORIGIN',nint(xh(1)/dgrid),
     +      nint(yh(1)/dgrid),nint(zh(1)/dgrid)
          call pretty (myline)
          write (junit,'(a)') myline(1:leng1(myline))
c
          write (myline,8000) 'EXTENT',mgrid(1),mgrid(2),mgrid(3)
          call pretty (myline)
          write (junit,'(a)') myline(1:leng1(myline))
c
          write (myline,8000) 'GRID',nscale,nscale,nscale
          call pretty (myline)
          write (junit,'(a)') myline(1:leng1(myline))
c
          write (myline,8010) 'SCALE',1.0
          call pretty (myline)
          write (junit,'(a)') myline(1:leng1(myline))
c
          write (junit,'(a)') 'MAP'
          write (junit,8030) (ibuff2(i),i=1,mxyz)
c
c ... 940715 - last card must be END
c
          write (junit,'(a)') 'END'
c
        end if
c
        close (junit)
c
 7000 format (8i5)
 7010 format (8f10.3)
 7020 format (a3,17x,a20,10x,1pe20.6)
 7030 format (1x,75a1)
c
 8000 format (a,1x,8i5)
 8010 format (a,1x,8f10.3)
 8030 format (39(1x,a1))
c
      end if
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '32 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Plotting all protein surface points done',
     +  '----------------------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine ezd_write (mgrid,xh,yh,zh)
c
      include 'cavity.incl'
c
      real xh(*),yh(*),zh(*)
c
      integer mgrid(3),mxyz,nscale
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Writing Old-EZD file for cavity',
     +  '-------------------------------',' '
c
      mxyz = mgrid(1) * mgrid(2) * mgrid(3)
c
      do i=1,mxyz
        if (ibuff2(i) .ne. cavi) ibuff2 (i) = '0'
      end do
c
      do i=1,mxyz
        if (ibuff2(i) .eq. cavi) ibuff2 (i) = '1'
      end do
c
      nscale = 100
 6511 continue
      if (max(mgrid(1),mgrid(2),mgrid(3)) .gt. (nscale-20)) then
        nscale = nscale * 2
        goto 6511
      end if
c
      write (junit,7000) nint(xh(1)/pgrid),nint(yh(1)/pgrid),
     +  nint(zh(1)/pgrid)
      write (junit,7000) mgrid(1),mgrid(2),mgrid(3)
      write (junit,7000) nscale,nscale,nscale
      write (junit,7010) float(nscale)*pgrid,float(nscale)*pgrid,
     +  float(nscale)*pgrid,90.0,90.0,90.0
      write (junit,7020) 'MAP','(1X,75F1.0)',1.0
      write (junit,7030) (ibuff2(i),i=1,mxyz)
c
 7000 format (8i5)
 7010 format (8f10.3)
 7020 format (a3,17x,a20,10x,1pe20.6)
 7030 format (1x,75a1)
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '33 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Writing Old-EZD file for cavity done',
     +  '------------------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c ===========================================================================
c
      subroutine sphere_volume (vol)
c
      implicit none
c
      real twopi,factor,third
c
      parameter (twopi  = 6.2831853071796)
      parameter (factor = 0.75/(0.50*twopi))
      parameter (third  = 1.0/3.0)
c
      real vol,rad
c
code ...
c
      rad = (factor * vol) ** third
      call rvalut (
     +  ' Volume corresponds to a sphere of radius (A) :',1,rad)
c
      return
      end
c
c ===========================================================================
c
      subroutine new_ezd_write (mgrid,xh,yh,zh)
c
      include 'cavity.incl'
c
      real xh(*),yh(*),zh(*)
c
      integer mgrid(3),mxyz,nscale,length,leng1
c
      character myline*120
c
code ...
c
      if (ltrace) write (tty,1200) ' ',
     +  'Writing New-EZD file for cavity',
     +  '-------------------------------',' '
c
      mxyz = mgrid(1) * mgrid(2) * mgrid(3)
c
      do i=1,mxyz
        if (ibuff2(i) .ne. cavi) ibuff2 (i) = '0'
      end do
c
      do i=1,mxyz
        if (ibuff2(i) .eq. cavi) ibuff2 (i) = '1'
      end do
c
      nscale = 100
 6511 continue
      if (max(mgrid(1),mgrid(2),mgrid(3)) .gt. (nscale-20)) then
        nscale = nscale * 2
        goto 6511
      end if
c
      write (junit,'(a)') 'EZD_MAP'
c
      call stamp (myline)
      myline = '! '//myline
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,8010) 'CELL',float(nscale)*pgrid,
     +  float(nscale)*pgrid,float(nscale)*pgrid,90.0,90.0,90.0
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,8000) 'ORIGIN',nint(xh(1)/pgrid),
     +  nint(yh(1)/pgrid),nint(zh(1)/pgrid)
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,8000) 'EXTENT',mgrid(1),mgrid(2),mgrid(3)
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,8000) 'GRID',nscale,nscale,nscale
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (myline,8010) 'SCALE',1.0
      call pretty (myline)
      write (junit,'(a)') myline(1:leng1(myline))
c
      write (junit,'(a)') 'MAP'
      write (junit,8030) (ibuff2(i),i=1,mxyz)
c
c ... 940715 - last card must be END
c
      write (junit,'(a)') 'END'
c
 8000 format (a,1x,8i5)
 8010 format (a,1x,8f10.3)
cxyz
c ... O 5.9.1 -> fix at SEVEN items per line
c
c 8030 format (39(1x,a1))
 8030 format (7(1x,a1))
cxyz
c
c ... end
c
  999 continue
      call gkdcpu (total,user,sys)
      write (tty,'(1x,a,3f10.1)') '34 CPU total/user/sys :',
     +  total,user,sys
c
      if (ltrace) write (tty,1200) ' ',
     +  'Writing New-EZD file for cavity done',
     +  '------------------------------------',' '
c
      ierror = 0
      return
c
 1200 format (99(1x,a:/))
c
      end
c
c
c
      subroutine add_probe (buff,mgrid,probe,hsize,cavi,temp,nex)
c
c ... add all points within the probe radius from the cavity
c     surface to the cavity itself
c
      implicit none
c
      integer maxdif
      parameter (maxdif=1000)
c
      real hsize,probe,pr2,dx,dy,dz,sqrs(-maxdif:maxdif)
c
      integer mgrid(3),dijk,i,j,k,ii,jj,kk,j0,k0,jj0,kk0
      integer i1,i2,j1,j2,k1,k2,nex,mxy
c
      character buff(*)*1,cavi*1,temp*1
c
code ...
c
c ... calculate range of points around each cavity point
c
      dijk = 1 + int ( probe / hsize )
c
c ... calculate square of probe radius in "grid units"
c
      pr2 = (probe**2) / (hsize**2)
c
      mxy = mgrid(1)*mgrid(2)
c
c ... pre-calculate I**2
c
      do i=0,maxdif
        sqrs (i) = i**2
        sqrs (-i) = sqrs (i)
      end do
c
      do i=1,mgrid(1)
        i1 = max (1, i-dijk)
        i2 = min (mgrid(1), i+dijk)
        do j=1,mgrid(2)
          j0 = (j-1)*mgrid(1) + i
          j1 = max (1, j-dijk)
          j2 = min (mgrid(2), j+dijk)
          do k=1,mgrid(3)
            k0 = (k-1)*mxy + j0
            if (buff(k0) .eq. cavi) then
              k1 = max (1, k-dijk)
              k2 = min (mgrid(3), k+dijk)
cc              print *,'I -',i,i1,i2
cc              print *,'J -',j,j1,j2
cc              print *,'K -',k,k1,k2
              do ii=i1,i2
                dx = sqrs (ii-i)
                if (dx .gt. pr2) goto 110
                do jj=j1,j2
                  dy = dx + sqrs (jj-j)
                  if (dy .gt. pr2) goto 120
                  jj0 = (jj-1)*mgrid(1) + ii
                  do kk=k1,k2
                    dz = dy + sqrs (kk-k)
                    if (dz .gt. pr2) goto 130
                    kk0 = (kk-1)*mxy + jj0
                    if (buff(kk0) .eq. cavi) goto 130
                    if (buff(kk0) .eq. temp) goto 130
                    buff(kk0) = temp
cc                    print *,'+ ',ii,jj,kk
  130               continue
                  end do
  120             continue
                end do
  110           continue
              end do
            end if
          end do
        end do
      end do
c
      nex = 0
      do i=1,mgrid(1)
        do j=1,mgrid(2)
          j0 = (j-1)*mgrid(1) + i
          do k=1,mgrid(3)
            k0 = (k-1)*mxy + j0
            if (buff(k0) .eq. temp) then
              nex = nex + 1
              buff(k0) = cavi
            end if
          end do
        end do
      end do
c
      return
      end
