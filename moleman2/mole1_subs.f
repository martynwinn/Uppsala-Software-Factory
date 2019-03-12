c
c ... mole1_subs.f - B-FACTOR/OCCUPANCY/CHAIN/DISTANCE/VRML
c                    subroutines for MOLEMAN2
c
c ... i.e., subroutines that *DO* include 'moleman2.incl'
c
c
      subroutine bq_stats (how,which)
c
      include 'moleman2.incl'
c
      real ave,sdv,xmin,xmax,xtot,xrms,xhav
c
      integer i,j,nat,nch,mode
c
      character*(*) how,which
      character now*1,seg*4
c
code ...
c
      call upcase (how)
      call remspa (how)
      if (how(1:1) .ne. 'T') how = 'C'
c
      mode = 1
      if (which(1:1) .eq. 'Q') mode = 2
c
 6000 format (8(1x,a10))
 6010 format (1x,a10,1x,i10,6(3x,f8.2))
c
c ... chain-based B-factor/Occupancy statistics
c
      if (how(1:1) .eq. 'C') then
c
        now = char(1+ichar(achain(1)))
        seg = '????'
        nch = 0
        nat = 0
        if (mode .eq. 1) then
          write (*,6000) 'Chain name','Atoms',
     +      'Bave','Bsdv','Bmin','Bmax','Brms','Bharm.av'
        else
          write (*,6000) 'Chain name','Atoms',
     +      'Qave','Qsdv','Qmin','Qmax','Qrms','Qharm.av'
        end if
        achain (natoms+1) = char(1+ichar(achain(natoms)))
c
        do i=1,natoms+1
          if (achain(i) .ne. now .or.
     +        inote(i)(7:10) .ne. seg) then
            if (nch.gt.0 .and. nat.gt.0) then
              call xstats (rbuf,nat,ave,sdv,xmin,xmax,xtot)
              call xstat2 (rbuf,nat,xrms,xhav)
              write (*,6010) (now//'<->'//seg),nat,ave,sdv,
     +          xmin,xmax,xrms,xhav
            end if
            if (i.gt.natoms) goto 1000
            nch = nch + 1
            nat = 0
            now = achain(i)
            seg = inote(i)(7:10)
          end if
          if (select(i)) then
            nat = nat + 1
            if (mode .eq. 1) then
              rbuf (nat) = batom(i)
            else
              rbuf (nat) = qatom(i)
            end if
          end if
        end do
 1000   continue
        call jvalut (' Nr of chains encountered  :',1,nch)
c
      else if (how(1:1) .eq. 'T') then
c
        if (mode .eq. 1) then
          write (*,6000) 'Type','Atoms',
     +      'Bave','Bsdv','Bmin','Bmax','Brms','Bharm.av'
        else
          write (*,6000) 'Type','Atoms',
     +      'Qave','Qsdv','Qmin','Qmax','Qrms','Qharm.av'
        end if
c
        do i=1,nrlcat+1
          if (icnt(i) .gt. 0) then
            nat = 0
            do j=1,natoms
              if (select(j)) then
                if (restyp(resptr(j)) .eq. i) then
                  nat = nat + 1
                  if (mode .eq. 1) then
                    rbuf (nat) = batom(j)
                  else
                    rbuf (nat) = qatom(j)
                  end if
                end if
              end if
            end do
            if (nat .gt. 0) then
              call xstats (rbuf,nat,ave,sdv,xmin,xmax,xtot)
              call xstat2 (rbuf,nat,xrms,xhav)
              write (*,6010) libcat(i),nat,ave,sdv,
     +          xmin,xmax,xrms,xhav
            end if
          end if
        end do
c
      end if
c
      return
      end
c
c
c
      subroutine bq_limit (which)
c
      include 'moleman2.incl'
c
      integer i,mode
c
      character*(*) which
c
code ...
c
      mode = 1
      if (which(1:1) .eq. 'Q') mode = 2
c
 6000 format (' Reset ',a,' to lie in range ',f8.2,' to ',f8.2)
c
      if (mode .eq. 1) then
        do i=1,natoms
          if (select(i)) then
            batom(i) = max (blimlo, min (blimhi, batom(i)))
          end if
        end do
        write (*,6000) 'B-factors',blimlo,blimhi
      else
        do i=1,natoms
          if (select(i)) then
            qatom(i) = max (qlimlo, min (qlimhi, qatom(i)))
          end if
        end do
        write (*,6000) 'occupancies',qlimlo,qlimhi
      end if
c
      return
      end
c
c
c
      subroutine bq_prod_plus (which,apr,apl)
c
      include 'moleman2.incl'
c
      real prod,plus
c
      integer i,mode,nok,ierr
c
      character*(*) which,apr,apl
c
code ...
c
      mode = 1
      if (which(1:1) .eq. 'Q') mode = 2
c
      call str2r (apr,prod,ierr)
      if (ierr .ne. 0) return
c
      call str2r (apl,plus,ierr)
      if (ierr .ne. 0) return
c
 6000 format (' New ',a,' = ',f10.4,' * Old ',a,' + ',f10.4)
c
      nok = 0
      if (mode .eq. 1) then
        write (*,6000) 'B',prod,'B',plus
        do i=1,natoms
          if (select(i)) then
            nok = nok + 1
            batom(i) = prod * batom(i) + plus
          end if
        end do
      else
        write (*,6000) 'Q',prod,'Q',plus
        do i=1,natoms
          if (select(i)) then
            nok = nok + 1
            qatom(i) = prod * qatom(i) + plus
          end if
        end do
      end if
c
      call jvalut (' Nr of atoms updated :',1,nok)
c
      return
      end
c
c
c
      subroutine bscale (blo,bhi)
c
      include 'moleman2.incl'
c
      real blo,bhi,xmin,xmax
c
      integer i,nok
c
code ...
c
      write (*,6000) blo,bhi
c
 6000 format (' Scale B-factors of selected atoms'/
     +        ' Requested range (A2) :',2f8.3)
 6010 format (' Current range (A2)   : ',2f8.3)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          if (nok .eq. 1) then
            xmin = batom (i)
            xmax = batom (i)
          else
            if (batom(i) .lt. xmin) xmin = batom (i)
            if (batom(i) .gt. xmax) xmax = batom (i)
          end if
        end if
      end do
c
      call jvalut (' Nr of selected atoms :',1,nok)
      if (nok .lt. 1) then
        call errcon ('At least two selected atoms needed for scaling')
        return
      end if
      write (*,6010) xmin,xmax
      if (xmin .ge. xmax) then
        call errcon ('Range is zero - cannot scale')
        return
      end if
c
      do i=1,natoms
        if (select(i)) then
          batom (i) = blo + (batom(i)-xmin)*(bhi-blo)/(xmax-xmin)
        end if
      end do
c
      call prompt (' B-factors scaled')
c
      return
      end
c
c
c
      subroutine bpseudo (how,par1,par2,par3,ierr)
c
      include 'moleman2.incl'
c
      real cog(3),rad,vol,xmin,xmax,xave
c
      integer ierr,i,nok,now,inow,j,nbr
c
      character*(*) how,par1,par2,par3
c
code ...
c
      ierr = -1
c
      if (how(1:1) .eq. 'X') then
        call prompt (' Pseudo-B = X = X-coordinate')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            batom (i) = xyz (1,i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'Y') then
        call prompt (' Pseudo-B = Y = Y-coordinate')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            batom (i) = xyz (2,i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'Z') then
        call prompt (' Pseudo-B = Z = Z-coordinate')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            batom (i) = xyz (3,i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'Q') then
        call prompt (' Pseudo-B = Q = Occupancy')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            batom (i) = qatom (i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'A') then
        call prompt (' Pseudo-B = A = Atom index (sequential)')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            batom (i) = float (i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'I') then
        call prompt (' Pseudo-B = I = Residue index (sequential)')
        nok = 0
        now = resptr(1) - 999
        inow = 0
        do i=1,natoms
          if (select(i)) then
            if (resptr(i) .ne. now) then
              inow = inow + 1
              now = resptr (i)
            end if
            batom (i) = float (inow)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'R') then
        call prompt (' Pseudo-B = R = Residue number (in sequence)')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            batom (i) = float (iresid (i))
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'G') then
        call prompt (' Pseudo-B = G = Distance from Centre-of-Gravity')
        nok = 0
        cog (1) = 0.0
        cog (2) = 0.0
        cog (3) = 0.0
        do i=1,natoms
          if (select(i)) then
            cog (1) = cog (1) + xyz (1,i)
            cog (2) = cog (2) + xyz (2,i)
            cog (3) = cog (3) + xyz (3,i)
            nok = nok + 1
          end if
        end do
        if (nok .lt. 1) then
          call errcon ('No atoms selected')
          return
        end if
        cog (1) = cog (1) / float (nok)
        cog (2) = cog (2) / float (nok)
        cog (3) = cog (3) / float (nok)
        call fvalut (' Centre-of-gravity :',3,cog)
        do i=1,natoms
          if (select(i)) then
            batom (i) = sqrt (
     +                    (xyz(1,i)-cog(1))*(xyz(1,i)-cog(1)) +
     +                    (xyz(2,i)-cog(2))*(xyz(2,i)-cog(2)) +
     +                    (xyz(3,i)-cog(3))*(xyz(3,i)-cog(3)) )
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'D') then
        call prompt (' Pseudo-B = D = Distance from a point')
        nok = 0
        call str2r (par1,cog(1),ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading X-coordinate of point')
          return
        end if
        call str2r (par2,cog(2),ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading Y-coordinate of point')
          return
        end if
        call str2r (par3,cog(3),ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading Z-coordinate of point')
          return
        end if
        call fvalut (' Point :',3,cog)
        do i=1,natoms
          if (select(i)) then
            batom (i) = sqrt (
     +                    (xyz(1,i)-cog(1))*(xyz(1,i)-cog(1)) +
     +                    (xyz(2,i)-cog(2))*(xyz(2,i)-cog(2)) +
     +                    (xyz(3,i)-cog(3))*(xyz(3,i)-cog(3)) )
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'H') then
        call prompt (' Pseudo-B = H = Halle B (LDM-prediction)')
        nok = 0
        call str2r (par1,rad,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading sphere radius')
          return
        end if
        call fvalut (' Radius (A) :',1,rad)
        if (rad .lt. 0.1) then
          call errcon ('Radius too small')
          return
        end if
c
        xave = 0.0
        do i=1,natoms
          if (select(i)) then
            nok = nok + 1
            xave = xave + batom (i)
          end if
        end do
        call jvalut (' Nr of atoms selected :',1,nok)
        if (nok .lt. 1) then
          call errcon ('No atoms selected')
          return
        end if
        xave = xave / float (nok)
        call rvalut (' Average experimental B (A2) :',1,xave)
c
        rad = rad*rad
        xmin = 0.0
        do i=1,natoms
          if (select(i)) then
            nbr = 0
            do j=1,natoms
              if (rad .ge.
     +             ( (xyz(1,i)-xyz(1,j))*(xyz(1,i)-xyz(1,j)) +
     +               (xyz(2,i)-xyz(2,j))*(xyz(2,i)-xyz(2,j)) +
     +               (xyz(3,i)-xyz(3,j))*(xyz(3,i)-xyz(3,j)) ) ) then
                nbr = nbr + 1
              end if
            end do
            batom (i) = twopi * twopi / float (nbr)
            xmin = xmin + batom (i)
          end if
        end do
c
        xmin = xmin / float (nok)
        call rvalut (' Average predicted B (A2) :',1,xmin)
        xmax = xave / xmin
        call rvalut (' ==> Scale factor :',1,xmax)
        do i=1,natoms
          if (select(i)) then
            batom (i) = xmax * batom (i)
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'N') then
        call prompt (' Pseudo-B = N = Nr of neighbour atoms in sphere')
        nok = 0
        call str2r (par1,rad,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading sphere radius')
          return
        end if
        call fvalut (' Radius (A) :',1,rad)
        if (rad .lt. 0.1) then
          call errcon ('Radius too small')
          return
        end if
        rad = rad*rad
        do i=1,natoms
          if (select(i)) then
            nbr = 0
            do j=1,natoms
              if (rad .ge.
     +             ( (xyz(1,i)-xyz(1,j))*(xyz(1,i)-xyz(1,j)) +
     +               (xyz(2,i)-xyz(2,j))*(xyz(2,i)-xyz(2,j)) +
     +               (xyz(3,i)-xyz(3,j))*(xyz(3,i)-xyz(3,j)) ) ) then
                nbr = nbr + 1
              end if
            end do
            batom (i) = float (nbr)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else if (how(1:1) .eq. 'C') then
        call prompt (' Pseudo-B = C = CX values (atomic convexity)')
        nok = 0
        call str2r (par1,rad,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading sphere radius')
          return
        end if
        call fvalut (' Radius (A) :',1,rad)
        if (rad .lt. 0.1) then
          call errcon ('Radius too small')
          return
        end if
        call str2r (par2,vol,ierr)
        if (ierr .ne. 0) then
          call errcon ('While reading atomic volume')
          return
        end if
        call fvalut (' Atomic volume (A3) :',1,vol)
        if (vol .lt. 0.1) then
          call errcon ('Atomic volume too small')
          return
        end if
        vol = 2.0 * twopi * (rad**3) / (3.0 * vol)
        rad = rad*rad
        do i=1,natoms
          if (select(i)) then
            nbr = 0
            do j=1,natoms
              if (rad .ge.
     +             ( (xyz(1,i)-xyz(1,j))*(xyz(1,i)-xyz(1,j)) +
     +               (xyz(2,i)-xyz(2,j))*(xyz(2,i)-xyz(2,j)) +
     +               (xyz(3,i)-xyz(3,j))*(xyz(3,i)-xyz(3,j)) ) ) then
                nbr = nbr + 1
              end if
            end do
            batom (i) = (vol / float (nbr)) - 1.0
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms processed :',1,nok)
c
      else
        call errcon ('Invalid property for Pseudo-B')
        return
      end if
c
c ... print statistics
c
      if (nok .le. 0) then
        call errcon ('No atoms selected !')
        return
      end if
c
      j = 0
      do i=1,natoms
        if (select(i)) then
          j = j + 1
          if (j .eq. 1) then
            xmin = batom (i)
            xmax = batom (i)
            xave = batom (i)
          else
            if (batom(i) .lt. xmin) xmin = batom (i)
            if (batom(i) .gt. xmax) xmax = batom (i)
            xave = xave + batom (i)
          end if
        end if
      end do
c
      write (*,*)
      xave = xave / float (j)
      call fvalut (' Minimum Pseudo-B :',1,xmin)
      if (xmin .le. 0.0) call prompt (
     +  ' Warning - non-positive values for Pseudo-B !')
      call fvalut (' Maximum Pseudo-B :',1,xmax)
      if (xmax .ge. 1000.0) call prompt (
     +  ' Warning - Pseudo-B values >= 1000 do not fit PDB format !')
      call fvalut (' Dynamic range    :',1,(xmax-xmin))
      call fvalut (' Average Pseudo-B :',1,xave)
      call prompt (' You may use BF LImit, BF PRod_plus and BF SCale')
      call prompt (' to modify the range and values of the Pseudo-Bs')
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine bgroup (how)
c
      include 'moleman2.incl'
c
      real sum,su2
c
      integer i,j,mode,nok,no2,nat,na2
c
      character*(*) how
c
code ...
c
      call upcase (how)
      call remspa (how)
      if (how(1:1) .eq. 'O') then
        mode = 1
        call prompt (' Group Bs for all selected atoms')
      else if (how(1:1) .eq. 'R') then
        mode = 2
        call prompt (' Group Bs per residue')
      else
        mode = 3
        call prompt (' Group Bs for main/side chain atoms')
      end if
c
      if (mode .eq. 1) then
        nok = 0
        sum = 0.0
        do i=1,natoms
          if (select(i)) then
            nok = nok + 1
            sum = sum + batom(i)
          end if
        end do
        call jvalut (' Nr of selected atoms :',1,nok)
        if (nok .lt. 1) return
        sum = sum / float(nok)
        call fvalut (' Average B-factor :',1,sum)
        do i=1,natoms
          if (select(i)) batom(i) = sum
        end do
c
      else if (mode .eq. 2) then
c
        nok = 0
        do i=1,nres
          nat = 0
          sum = 0.0
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              nat = nat + 1
              sum = sum + batom(j)
            end if
          end do
          if (nat .gt. 0) then
            sum = sum / float(nat)
            do j=atmptr(1,i),atmptr(2,i)
              if (select(j)) batom(j) = sum
            end do
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of residues grouped :',1,nok)
c
      else if (mode .eq. 3) then
c
        nok = 0
        no2 = 0
        do i=1,nres
          nat = 0
          sum = 0.0
          na2 = 0
          su2 = 0.0
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              if (lmain(j)) then
                nat = nat + 1
                sum = sum + batom(j)
              else
                na2 = na2 + 1
                su2 = su2 + batom(j)
              end if
            end if
          end do
          if (nat .gt. 0) then
            sum = sum / float(nat)
            do j=atmptr(1,i),atmptr(2,i)
              if (select(j) .and. lmain(j)) batom(j) = sum
            end do
            nok = nok + 1
          end if
          if (na2 .gt. 0) then
            su2 = su2 / float(na2)
            do j=atmptr(1,i),atmptr(2,i)
              if (select(j) .and. (.not.lmain(j)) ) batom(j) = su2
            end do
            no2 = no2 + 1
          end if
        end do
        call jvalut (' Nr of main-chain groups :',1,nok)
        call jvalut (' Nr of side-chain groups :',1,no2)
c
      end if
c
      return
      end
c
c
c
      subroutine bsmooth ()
c
      include 'moleman2.incl'
c
      real dist
c
      integer i,nok,j,k
c
code ...
c
      call prompt (' Smoothing Bs for bonded atoms')
      call fvalut (' Max bonded distance :',1,mxbond)
c
      nok = 0
      do i=1,natoms
        lbuf(i) = (select(i))
        if (lbuf(i)) nok = nok + 1
      end do
c
      call jvalut (' Nr of selected atoms :',1,nok)
      if (nok .lt. 2) then
        call errcon ('Not enough atoms selected')
        return
      end if
c
      do i=1,natoms
        if (lbuf(i)) then
          rbuf(i) = batom(i)
          k = 1
          do j=1,natoms
            if (i.eq.j) goto 100
            if (.not. lbuf(j)) goto 100
            if (abs(xyz(1,i)-xyz(1,j)) .gt. mxbond) goto 100
            if (abs(xyz(2,i)-xyz(2,j)) .gt. mxbond) goto 100
            if (abs(xyz(3,i)-xyz(3,j)) .gt. mxbond) goto 100
            if (dist(i,j,xyz) .gt. mxbond) goto 100
            k = k + 1
            rbuf(i) = rbuf(i) + batom(j)
  100       continue
          end do
          rbuf (i) = rbuf (i) / float(k)
        end if
      end do
c
      do i=1,natoms
        if (lbuf(i)) batom(i) = rbuf(i)
      end do
c
      return
      end
c
c
c
      subroutine bbond ()
c
      include 'moleman2.incl'
c
      real dist,sbo,sno,qq,db,smc,ssc,sms
c
      integer i,nok,j,nbo,non,nmc,nsc,nms
c
code ...
c
      call prompt (' Calculating RMS delta-B for (non-)bonded atoms')
      call fvalut (' Max bonded distance     :',1,mxbond)
      call fvalut (' Max non-bonded distance :',1,mxnonb)
c
      nok = 0
      do i=1,natoms
        lbuf(i) = (select(i))
        if (lbuf(i)) nok = nok + 1
      end do
c
      call jvalut (' Nr of selected atoms :',1,nok)
      if (nok .lt. 2) then
        call errcon ('Not enough atoms selected')
        return
      end if
c
      nbo = 0
      non = 0
      nmc = 0
      nsc = 0
      nms = 0
      sbo = 0.0
      sno = 0.0
      smc = 0.0
      ssc = 0.0
      sms = 0.0
c
      do i=1,natoms-1
        if (lbuf(i)) then
          do j=i+1,natoms
            if (i.eq.j) goto 100
            if (.not. lbuf(j)) goto 100
            if (abs(xyz(1,i)-xyz(1,j)) .gt. mxnonb) goto 100
            if (abs(xyz(2,i)-xyz(2,j)) .gt. mxnonb) goto 100
            if (abs(xyz(3,i)-xyz(3,j)) .gt. mxnonb) goto 100
            qq = dist(i,j,xyz)
            if (qq .gt. mxnonb) goto 100
            db = (batom(i)-batom(j))**2
            if (qq .le. mxbond) then
              nbo = nbo + 1
              sbo = sbo + db
              if (lmain(i) .and. lmain(j)) then
                nmc = nmc + 1
                smc = smc + db
              else if ((.not.lmain(i)).and.(.not. lmain(j))) then
                nsc = nsc + 1
                ssc = ssc + db
              else
                nms = nms + 1
                sms = sms + db
              end if
              goto 100
            end if
            non = non + 1
            sno = sno + db
  100       continue
          end do
        end if
      end do
c
      call ivalut (  ' Nr of bonds (all atoms)   :',1,nbo)
      if (nbo .gt. 0) then
        call fvalut ('        RMS delta-B (A**2) :',1,
     +    sqrt(sbo/float(nbo)))
      end if
c
      call ivalut (  ' Nr of bonds (main-main)   :',1,nmc)
      if (nmc .gt. 0) then
        call fvalut ('        RMS delta-B (A**2) :',1,
     +    sqrt(smc/float(nmc)))
      end if
c
      call ivalut (  ' Nr of bonds (side-side)   :',1,nsc)
      if (nsc .gt. 0) then
        call fvalut ('        RMS delta-B (A**2) :',1,
     +    sqrt(ssc/float(nsc)))
      end if
c
      call ivalut (  ' Nr of bonds (main-side)   :',1,nms)
      if (nms .gt. 0) then
        call fvalut ('        RMS delta-B (A**2) :',1,
     +    sqrt(sms/float(nms)))
      end if
c
      call ivalut (  ' Nr of non-bonded contacts :',1,non)
      if (non .gt. 0) then
        call fvalut ('        RMS delta-B (A**2) :',1,
     +    sqrt(sno/float(non)))
      end if
c
      return
      end
c
c
c
      subroutine bodb (iunit,file,molnam)
c
      include 'moleman2.incl'
c
      integer i,iunit,ierr,length
c
      character*(*) file,molnam
      character line*256,myfmt*20
c
code ...
c
      call upcase (molnam)
c
      call textut (' Save B-values in ODB file :',file)
      call textut (' Molecule name in O :',molnam)
c
      call xopxua (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      write (iunit,'(a1,1x,a)') '!',line(1:length(line))
c
      myfmt = '(10F8.3)'
      write (iunit,'(2a,i10,2a)',err=999) molnam(1:length(molnam)),
     +  '_ATOM_B     R    ',natoms,' ',myfmt(1:length(myfmt))
c
      write (iunit,myfmt,err=999) (batom(i),i=1,natoms)
c
      call prompt (' File written')
c
      goto 1000
c
  999 continue
      call errcon ('While writing ODB file')
c
 1000 continue
      close (iunit)
c
      return
      end
c
c
c
      subroutine bplot (iunit,file,how)
c
      include 'moleman2.incl'
c
c ... size of distance bins for radial B-factor plots
c
      real bin
      parameter (bin = 2.0)
c
      real cog(3)
      real sum,ave,sdv,xmin,xmax,xtot,d,bmax
c
      integer i,iunit,mode,j,nok,ndo,ierr,maxbin,k,kmax
c
      character*(*) file,how
      character line*256
c
code ...
c
c ... max nr of distance bins
c
      maxbin = 100
c
      mode = 1
      call upcase (how)
      call remspa (how)
      if (how(1:1) .eq. 'R') mode = 2
c
      call xopxua (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      if (mode .eq. 1) then
c
c ... average B per residue of selected atoms
c
        ndo = 0
        do i=1,nres
          nok = 0
          sum = 0
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              nok = nok + 1
              sum = sum + batom(j)
            end if
          end do
          if (nok .gt. 0) then
            ndo = ndo + 1
            rbuf (ndo) = sum / float(nok)
          end if
        end do
        call jvalut (' Residues with selected atoms :',1,ndo)
        if (ndo .lt. 3) then
          call errcon ( 'Fewer than three residues')
          goto 1000
        end if
        call xstats (rbuf,ndo,ave,sdv,xmin,xmax,xtot)
c
        call stamp (line)
        line = 'REMARK '//line
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Average temperature-factor per residue'//
     +    ' for selected atoms'
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Generated from PDB file '//pdbfil
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Selected atoms: '//selstr
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'REMARK Min ',xmin,' Max ',xmax
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'REMARK Ave ',ave,' Sdv ',sdv
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'NPOINT',ndo
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'COLOUR 4',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,'XLABEL Residue',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,'YLABEL Average B-factor',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'XYVIEW',0,ndo+1,0,xmax+0.01*(xmax-xmin)
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'XLIMIT 1.0 1.0',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,'YVALUE *',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        write (iunit,'(10f8.2)',err=1000) (rbuf(i),i=1,ndo)
c
        call putlin (iunit,'END',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        close (iunit)
        return
c
c ... radial B-factor plot
c
      else
c
        do i=1,maxbin
          ibuf (i) = 0
          rbuf (i) = 0.0
        end do
c
        do i=1,3
          cog(i) = 0.0
        end do
c
        nok = 0
        do i=1,natoms
          if (select(i)) then
            nok = nok + 1
            do j=1,3
              cog(j) = cog(j) + xyz(j,i)
            end do
          end if
        end do
c
        call jvalut (' Selected atoms :',1,nok)
        if (nok .lt. 10) then
          call errcon ( 'Fewer than ten atoms')
          goto 1000
        end if
c
        do i=1,3
          cog(i) = cog(i) / float (nok)
        end do
        call fvalut (' Centre-of-gravity :',3,cog)
c
        kmax = 0
        do i=1,natoms
          if (select(i)) then
            d = 0.0
            do j=1,3
              d = d + (xyz(j,i)-cog(j))**2
            end do
            d = sqrt(d)
            k = max (1, min (maxbin,int(d/bin)+1))
            if (k.gt.kmax) kmax = k
            ibuf (k) = ibuf (k) + 1
            rbuf (k) = rbuf (k) + batom(i)
ccc            if (k.eq.1) call print_atom(i)
          end if
        end do
c
 6000 format (' Shell ',f6.1,' - ',f6.1,' A - ',i6,' atoms; <B> =',
     +  f6.2,' A**2')
c
        bmax = 0.0
        do i=1,maxbin
          if (ibuf(i) .gt. 0) then
            rbuf(i) = rbuf(i) / float(ibuf(i))
            write (*,6000) bin*(i-1),bin*i,ibuf(i),rbuf(i)
            if (rbuf(i) .gt. bmax) bmax = rbuf(i)
          end if
        end do
c
        call stamp (line)
        line = 'REMARK '//line
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Radial temperature-factor plot'//
     +    ' for selected atoms'
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Generated from PDB file '//pdbfil
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Selected atoms: '//selstr
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'NPOINT',kmax
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'REMARK Bin size (A) : ',bin
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'COLOUR 4',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,
     +    'XLABEL Distance from Centre-of-Gravity (A)',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,
     +    'YLABEL Average B in bin (A**2)',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'XYVIEW',-1,bin*float(kmax+1),0,1.01*bmax
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'XLIMIT 0',bin
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'YVALUE *',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        write (iunit,'(10f8.2)',err=1000) (rbuf(i),i=1,kmax)
c
        call putlin (iunit,'END',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        close (iunit)
        return
c
      end if
c
 1000 continue
      call errcon ('While writing plot file')
      close (iunit)
      return
c
      end
c
c
c
      subroutine qplot (iunit,file,how)
c
      include 'moleman2.incl'
c
c ... size of distance bins for radial occupancy plots
c
      real bin
      parameter (bin = 2.0)
c
      real cog(3)
      real sum,ave,sdv,xmin,xmax,xtot,d,bmax
c
      integer i,iunit,mode,j,nok,ndo,ierr,maxbin,k,kmax
c
      character*(*) file,how
      character line*256
c
code ...
c
c ... max nr of distance bins
c
      maxbin = 100
c
      mode = 1
      call upcase (how)
      call remspa (how)
      if (how(1:1) .eq. 'R') mode = 2
c
      call xopxua (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      if (mode .eq. 1) then
c
c ... average Q per residue of selected atoms
c
        ndo = 0
        do i=1,nres
          nok = 0
          sum = 0
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              nok = nok + 1
              sum = sum + qatom(j)
            end if
          end do
          if (nok .gt. 0) then
            ndo = ndo + 1
            rbuf (ndo) = sum / float(nok)
          end if
        end do
        call jvalut (' Residues with selected atoms :',1,ndo)
        if (ndo .lt. 3) then
          call errcon ( 'Fewer than three residues')
          goto 1000
        end if
        call xstats (rbuf,ndo,ave,sdv,xmin,xmax,xtot)
c
        call stamp (line)
        line = 'REMARK '//line
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Average occupancy per residue'//
     +    ' for selected atoms'
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Generated from PDB file '//pdbfil
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Selected atoms: '//selstr
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'REMARK Min ',xmin,' Max ',xmax
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'REMARK Ave ',ave,' Sdv ',sdv
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'NPOINT',ndo
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'COLOUR 4',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,'XLABEL Residue',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,'YLABEL Average B-factor',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'XYVIEW',0,ndo+1,0,xmax+0.01*(xmax-xmin)
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'XLIMIT 1.0 1.0',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,'YVALUE *',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        write (iunit,'(10f8.2)',err=1000) (rbuf(i),i=1,ndo)
c
        call putlin (iunit,'END',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        close (iunit)
        return
c
c ... radial occupancy plot
c
      else
c
        do i=1,maxbin
          ibuf (i) = 0
          rbuf (i) = 0.0
        end do
c
        do i=1,3
          cog(i) = 0.0
        end do
c
        nok = 0
        do i=1,natoms
          if (select(i)) then
            nok = nok + 1
            do j=1,3
              cog(j) = cog(j) + xyz(j,i)
            end do
          end if
        end do
c
        call jvalut (' Selected atoms :',1,nok)
        if (nok .lt. 10) then
          call errcon ( 'Fewer than ten atoms')
          goto 1000
        end if
c
        do i=1,3
          cog(i) = cog(i) / float (nok)
        end do
        call fvalut (' Centre-of-gravity :',3,cog)
c
        kmax = 0
        do i=1,natoms
          if (select(i)) then
            d = 0.0
            do j=1,3
              d = d + (xyz(j,i)-cog(j))**2
            end do
            d = sqrt(d)
            k = max (1, min (maxbin,int(d/bin)+1))
            if (k.gt.kmax) kmax = k
            ibuf (k) = ibuf (k) + 1
            rbuf (k) = rbuf (k) + qatom(i)
ccc            if (k.eq.1) call print_atom(i)
          end if
        end do
c
 6000 format (' Shell ',f6.1,' - ',f6.1,' A - ',i6,
     +        ' atoms; <OCC> =',f6.2)
c
        bmax = 0.0
        do i=1,maxbin
          if (ibuf(i) .gt. 0) then
            rbuf(i) = rbuf(i) / float(ibuf(i))
            write (*,6000) bin*(i-1),bin*i,ibuf(i),rbuf(i)
            if (rbuf(i) .gt. bmax) bmax = rbuf(i)
          end if
        end do
c
        call stamp (line)
        line = 'REMARK '//line
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Radial occupancy plot'//
     +    ' for selected atoms'
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Generated from PDB file '//pdbfil
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        line = 'REMARK Selected atoms: '//selstr
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'NPOINT',kmax
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'REMARK Bin size (A) : ',bin
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'COLOUR 4',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,
     +    'XLABEL Distance from Centre-of-Gravity (A)',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        call putlin (iunit,
     +    'YLABEL Average B in bin (A**2)',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'XYVIEW',-1,bin*float(kmax+1),0,1.01*bmax
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        write (line,*) 'XLIMIT 0',bin
        call pretty (line)
        call putlin (iunit,line,.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        call putlin (iunit,'YVALUE *',.false.,ierr)
        if (ierr .ne. 0) goto 1000
        write (iunit,'(10f8.2)',err=1000) (rbuf(i),i=1,kmax)
c
        call putlin (iunit,'END',.false.,ierr)
        if (ierr .ne. 0) goto 1000
c
        close (iunit)
        return
c
      end if
c
 1000 continue
      call errcon ('While writing plot file')
      close (iunit)
      return
c
      end
c
c
c
      subroutine chain_auto (imode)
c
c ... mode = 0 -> automatic
c            1 -> ask user
c
      include 'moleman2.incl'
c
      real dist
c
      integer i,j,imode,mode,iz,ia,nseg,i1,i2,i3,i4
c
      logical lnews
c
      character nsegid*4,asksid*4,askchn*1,reason*80
c
code ...
c
      mode = imode
      if (mode .ne. 1) mode = 0
c
      ia = ichar('A')
      iz = ichar('Z')
c
      i1 = ia
      i2 = ia
      i3 = ia
      i4 = ia - 1
c
      nseg = 0
c
      do i=1,nres
        lnews = .false.
c
c ... first residue -> new segment
c
        if (i .eq. 1) then
          reason = 'First residue'
          goto 1000
        end if
c
c ... different type (e.g., PROTein -> WATEr) -> new segment
c
        if (restyp(i) .ne. restyp(i-1)) then
          reason = 'Different type of residue: '//
     +      libcat(restyp(i))
          goto 1000
        end if
c
c ... if both PROTein, check CA-CA distance
c
        if (restyp(i) .eq. iprot) then
          if (captr(i) .le. 0) then
            call errcon ('Missing CA atom !!!')
            reason = 'Missing CA atom = ERROR !!!'
            goto 1000
          end if
          if (dist(captr(i),captr(i-1),xyz) .gt. mxcaca) then
            reason = 'CA-CA distance to previous residue too long'
            goto 1000
          end if
        end if
c
c ... if both WATEr -> never a new segment
c
        if (restyp(i) .eq. iwate) then
          goto 2000
        end if
c
c ... if other, check if sequential residue numbers
c
        if ( iresid(atmptr(1,i)) .ne.
     +       (1+iresid(atmptr(1,i-1))) ) then
          reason = 'Break in sequential residue numbering'
          goto 1000
        end if
c
c ... else, no new segment
c
        goto 2000
c
c ... new segment if here
c
 1000   continue
        nseg = nseg + 1
        call print_res(i,1)
        call jvalut (' New segment :',1,nseg)
        call textut (' Because :',reason)
c
        i4 = i4 + 1
        if (i4 .gt. iz) then
          i4 = ia
          i3 = i3 + 1
          if (i3 .gt. iz) then
            i3 = ia
            i2 = i2 + 1
            if (i2 .gt. iz) then
              i2 = ia
              i1 = i1 + 1
              if (i1 .gt. iz) then
                call errcon ('Too many segments (bloody hell !)')
                i1 = ia
              end if
            end if
          end if
        end if
c
        nsegid = char(i1)//char(i2)//char(i3)//char(i4)
        asksid = nsegid
        askchn = nsegid(4:4)
        if (mode .eq. 1) then
          call textin (' Segment ID (4 char) ?',asksid)
          call upcase (asksid)
          askchn = asksid(4:4)
          call textin (' Chain name (1 char) ?',askchn)
          call upcase (askchn)
          i1 = max (ia, min (iz, ichar(asksid(1:1)) ))
          i2 = max (ia, min (iz, ichar(asksid(2:2)) ))
          i3 = max (ia, min (iz, ichar(asksid(3:3)) ))
          i4 = max (ia, min (iz, ichar(asksid(4:4)) ))
        else
          call textut (' Segment ID (4 char) :',asksid)
          call textut (' Chain name (1 char) :',askchn)
        end if
c
c ... set segment and chain
c
 2000   continue
        do j=atmptr(1,i),atmptr(2,i)
          inote(j)(7:10) = asksid
          achain(j) = askchn
        end do
c
      end do
c
      write (*,*)
      call jvalut (' Total number of segment/chains found :',1,nseg)
c
      return
      end
c
c
c
      subroutine chain_rename (old,new)
c
      include 'moleman2.incl'
c
      integer i,nok
c
      character old*1,new*1
c
code ...
c
      call upcase (old)
      call upcase (new)
      write (*,6000) old,new
 6000 format (' Replace chain name |',a1,'| by |',a1,'|')
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          if (achain(i) .eq. old .or. old .eq. '*') then
            achain (i) = new
            nok = nok + 1
          end if
        end if
      end do
c
      call jvalut (' Nr of chain names replaced :',1,nok)
c
      return
      end
c
c
c
      subroutine chain_reseg (old,new)
c
      include 'moleman2.incl'
c
      integer i,nok
c
      character old*4,new*4
c
code ...
c
      call upcase (old)
      call upcase (new)
      write (*,6000) old,new
 6000 format (' Replace segment id |',a4,'| by |',a4,'|')
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          if (inote(i)(7:10) .eq. old .or. old .eq. '*   ') then
            inote(i)(7:10) = new
            nok = nok + 1
          end if
        end if
      end do
c
      call jvalut (' Nr of segment ids replaced :',1,nok)
c
      return
      end
c
c
c
      subroutine chain_select (chn,seg)
c
      include 'moleman2.incl'
c
      integer i,no1,no2
c
      character chn*1,seg*4
c
code ...
c
      call upcase (chn)
      call upcase (seg)
      write (*,6000) chn,seg
 6000 format (' Give selected atoms chain name |',a1,
     +  '| and segment id |',a4,'|')
c
      no1 = 0
      no2 = 0
      do i=1,natoms
        if (select(i)) then
          if (chn .ne. '=') then
            achain (i) = chn
            no1 = no1 + 1
          end if
          if (seg .ne. '=') then
            inote(i)(7:10) = seg
            no2 = no2 + 1
          end if
        end if
      end do
c
      call jvalut (' Nr of chain names replaced :',1,no1)
      call jvalut (' Nr of segment ids replaced :',1,no2)
c
      return
      end
c
c
c
      subroutine chain_from_seg (how)
c
      include 'moleman2.incl'
c
      integer i
c
      logical lask
c
      character how*(*),chn*1,seg*4
c
code ...
c
      call remspa (how)
      call upcase (how)
      if (how(1:2) .ne. 'AS') how = 'AUto'
      lask = (how(1:2) .eq. 'AS')
c
      chn = '?'
      seg = '????'
c
      do i=1,natoms
        if (inote(i)(7:10) .ne. seg) then
          seg = inote(i)(7:10)
          chn = seg(4:4)
          call print_res (resptr(i),1)
          if (lask) call textin (' Chain name ?',chn)
          call upcase (chn)
          call textut (' New chain name :',chn)
        end if
        achain (i) = chn
      end do
c
      return
      end
c
c
c
      subroutine chain_to_seg (how)
c
      include 'moleman2.incl'
c
      integer i
c
      logical lask
c
      character how*(*),chn*1,seg*4
c
code ...
c
      call remspa (how)
      call upcase (how)
      if (how(1:2) .ne. 'AS') how = 'AUto'
      lask = (how(1:2) .eq. 'AS')
c
      chn = '?'
      seg = '????'
c
      do i=1,natoms
        if (achain(i) .ne. chn) then
          chn = achain(i)
          seg = chn//chn//chn//chn
          call print_res (resptr(i),1)
          if (lask) call textin (' Segment id ?',seg)
          call upcase (seg)
          call textut (' New segment id :',seg)
        end if
        inote(i)(7:10) = seg
      end do
c
      return
      end
c
c
c
      subroutine chain_ot2 ()
c
      include 'moleman2.incl'
c
      real dist,q(3)
c
      integer i,j,ii,nat,nok
c
code ...
c
      call prompt (' Looking for protein chain ends ...')
c
      nok = 0
      do i=1,nres
        if (restyp(i) .eq. iprot) then
          if (i .eq. nres) goto 1000
          if (restyp(i+1) .ne. iprot) goto 1000
          if (captr(i) .le. 0 .or. captr(i+1) .le. 0) goto 2000
          if (dist(captr(i),captr(i+1),xyz) .gt. mxcaca) goto 1000
        end if
        goto 2000
c
c ... suggest them
c
 1000   continue
        call print_res (i,1)
        nat = atmptr(2,i) - atmptr(1,i) + 1
        ii = atmptr(1,i)
        call suggot (nat,atmnam(ii),xyz(1,ii),j,q,.true.)
        nok = nok + 1
c
 2000   continue
c
      end do
c
      if (nok .gt. 0) write (*,*)
      call jvalut (' Nr of chain ends found :',1,nok)
c
      return
      end
c
c
c
      subroutine dist_plot (file,f4)
c
      include 'moleman2.incl'
c
      real dist,ave,sdv,xmin,xmax,xtot
c
      integer i,j,f4,ierr,np,nd,leng1
c
      character file*(*),line*256
c
code ...
c
      call textut (' Distance plot file   :',file)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      if (nselect .lt. 3) then
        call errcon ('Not enough atoms')
        return
      end if
c
      np = nselect*nselect
      if (np .gt. maxbuf) then
        call errcon ('Too many atoms selected for 2D contour plot')
        call jvalut (' Requested points :',1,np)
        call jvalut (' Max buffer size  :',1,maxbuf)
        return
      end if
c
 5000 format (a6,1x,a)
 5010 format (a6,1x,12i6)
 5020 format (a6,1x,6f12.4)
c
      call xopxua (f4,file,linter,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      write (f4,5000,err=1000) 'REMARK',
     +    'Distance plot'
      write (f4,5000,err=1000) 'REMARK',line(1:leng1(line))
c
      line = 'REMARK Generated from PDB file '//pdbfil
      call putlin (f4,line,.false.,ierr)
      if (ierr .ne. 0) goto 1000
c
      line = 'REMARK Selected atoms: '//selstr
      call putlin (f4,line,.false.,ierr)
      if (ierr .ne. 0) goto 1000
c
      write (f4,5000,err=1000) 'REMARK'
      write (f4,5000,err=1000) 'XLABEL','Selected atom'
      write (f4,5000,err=1000) 'YLABEL','Selected atom'
      write (f4,5010,err=1000) 'NLEVEL',7
      write (f4,5000,err=1000) 'LEVELS'
      write (f4,5000,err=1000) '6 7 8 ','9 10 11 12'
      write (f4,5000,err=1000) 'COLOUR'
      write (f4,5000,err=1000) '1 1 5 ','5 2 6 4'
      write (f4,5010,err=1000) 'XPOINT',nselect
      write (f4,5010,err=1000) 'YPOINT',nselect
      write (f4,5010,err=1000) 'XLIMIT',1,nselect
      write (f4,5010,err=1000) 'YLIMIT',1,nselect
      write (f4,5000,err=1000) 'ZVALUE','*'
c
      nd = 0
      do i=1,natoms
        if (select(i)) then
          do j=1,natoms
            if (select(j)) then
              nd = nd + 1
              rbuf (nd) = dist (i,j,xyz)
            end if
          end do
        end if
      end do
c
      call jvalut (' Nr of distances :',1,nd)
      call xstats (rbuf(1),nd,ave,sdv,xmin,xmax,xtot)
      call fvalut (' Average distance :',1,ave)
      call fvalut (' St. deviation    :',1,sdv)
      call fvalut (' Maximum distance :',1,xmax)
c
      write (f4,'(7f10.1)',err=1000) (rbuf(i),i=1,nd)
      write (f4,5000,err=1000) 'END   '
c
      close (f4)
      call prompt (' Plot file written')
      return
c
c ... error
c
 1000 continue
      call errcon ('While writing file')
      close (f4)
c
      return
      end
c
c
c
      subroutine conv_plot (file,ncolen,acotyp,f4)
c
      include 'moleman2.incl'
c
      integer maxcon
      parameter (maxcon=10000)
c
      real myxyz(3,maxcon),rt(12)
      real ave,sdv,xmin,xmax,xtot,rmsd
c
      integer i,j,f4,ierr,np,nd,leng1,ncolen,na
c
      character file*(*),line*256,acotyp*4
c
code ...
c
      call textut (' Convolution plot file :',file)
      call jvalut (' Fragment length :',1,ncolen)
      if (ncolen .lt. 3) then
        call errcon ('Fragment length too small')
        return
      end if
      call textut (' Atom type :',acotyp)
c
      na = 0
      do i=1,natoms
        if (select(i)) then
          if (atmnam(i) .eq. acotyp) then
            na = na + 1
            if (na .le. maxcon) then
              myxyz(1,na) = xyz(1,i)
              myxyz(2,na) = xyz(2,i)
              myxyz(3,na) = xyz(3,i)
            else
              call errcon ('Too many atoms; rest skipped')
              na = maxcon
              goto 1234
            end if
          end if
        end if
      end do
 1234 continue
      call jvalut (' Nr of selected atoms :',1,na)
c
      if (na .lt. (ncolen+1)) then
        call errcon ('Not enough atoms')
        return
      end if
c
      np = na*na
      if (np .gt. maxbuf) then
        call errcon ('Too many atoms selected for 2D contour plot')
        call jvalut (' Requested points :',1,np)
        call jvalut (' Max buffer size  :',1,maxbuf)
        return
      end if
c
 5000 format (a6,1x,a)
 5010 format (a6,1x,12i6)
 5015 format (a6,1x,a,1x,i8)
 5020 format (a6,1x,6f12.4)
 5025 format (a6,1x,a,1x,f8.2)
c
      call xopxua (f4,file,linter,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      write (f4,5000,err=1000) 'REMARK',
     +    'Convolution plot'
      write (f4,5015,err=1000) 'REMARK',
     +    'Fragment length',ncolen
      write (f4,5000,err=1000) 'REMARK',line(1:leng1(line))
c
      line = 'REMARK Generated from PDB file '//pdbfil
      call putlin (f4,line,.false.,ierr)
      if (ierr .ne. 0) goto 1000
c
      write (f4,5000,err=1000) 'REMARK'
      write (f4,5000,err=1000) 'XLABEL','Selected atom'
      write (f4,5000,err=1000) 'YLABEL','Selected atom'
      write (f4,5010,err=1000) 'NLEVEL',5
      write (f4,5000,err=1000) 'LEVELS'
c
      if (ncolen .le. 7) then
        write (f4,5000,err=1000) ' ','0.2 0.4 0.6 0.8 1.0'
      else if (ncolen .le. 12) then
        write (f4,5000,err=1000) ' ','0.3 0.6 0.9 1.2 1.5'
      else
        write (f4,5000,err=1000) ' ','0.4 0.8 1.2 1.6 2.0'
      end if
c
      write (f4,5000,err=1000) 'COLOUR'
      write (f4,5000,err=1000) ' ','1 5 2 6 4'
      write (f4,5010,err=1000) 'XPOINT',na-ncolen+1
      write (f4,5010,err=1000) 'YPOINT',na-ncolen+1
      write (f4,5010,err=1000) 'XLIMIT',1,na-ncolen+1
      write (f4,5010,err=1000) 'YLIMIT',1,na-ncolen+1
c
      nd = 0
      do i=1,na-ncolen+1
        do j=1,na-ncolen+1
          nd = nd + 1
          call lsqgjk (myxyz(1,i),myxyz(1,j),ncolen,rmsd,rt,ierr)
          if (ierr .ne. 0) rmsd = 99.99
          rbuf (nd) = rmsd
        end do
      end do
c
      call jvalut (' Nr of RMSDs   :',1,nd)
      call xstats (rbuf(1),nd,ave,sdv,xmin,xmax,xtot)
      call fvalut (' Average RMSD  :',1,ave)
      call fvalut (' St. deviation :',1,sdv)
      call fvalut (' Maximum RMSD  :',1,xmax)
c
      write (f4,5025,err=1000) 'REMARK',
     +    'Average RMSD (A)',ave
      write (f4,5025,err=1000) 'REMARK',
     +    'St. dev. (A)',sdv
c
      write (f4,5000,err=1000) 'ZVALUE','*'
      write (f4,'(10f7.2)',err=1000) (rbuf(i),i=1,nd)
      write (f4,5000,err=1000) 'END   '
c
      close (f4)
      call prompt (' Plot file written')
      return
c
c ... error
c
 1000 continue
      call errcon ('While writing file')
      close (f4)
c
      return
      end
c
c
c
      subroutine dist_dist (par)
c
      include 'moleman2.incl'
c
      real dist,qq,bin
c
      integer i,j,ierr,isum,nd,nb,ii
c
      character par*(*)
c
code ...
c
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      if (nselect .lt. 3) then
        call errcon ('Not enough atoms')
        return
      end if
c
      call str2r (par,bin,ierr)
      if (ierr .ne. 0) return
      call fvalut (' Bin size :',1,bin)
      if (bin .lt. 0.1) then
        call errcon ('Invalid bin size')
        return
      end if
c
      nb = min ( maxbuf, nint(1000.0/bin) )
      do i=1,nb
        ibuf (i) = 0
      end do
c
      nd = 0
      do i=1,natoms-1
        if (select(i)) then
          do j=i+1,natoms
            if (select(j)) then
              nd = nd + 1
              qq = dist (i,j,xyz)
              ii = min (nb, 1 + int (qq/bin))
              ibuf (ii) = ibuf (ii) + 1
            end if
          end do
        end if
      end do
c
      call jvalut (' Nr of distances :',1,nd)
      isum = 0
      write (*,6000) 'lower','upper','# dist','# cumul',
     +  '% dist','% cumul'
      do i=1,nb
        if (ibuf(i) .gt. 0) then
          isum = isum + ibuf(i)
          write (*,6010) (i-1)*bin,i*bin,ibuf(i),isum,
     +      100.0*float(ibuf(i))/float(nd),
     +      100.0*float(isum)/float(nd)
        end if
      end do
c
 6000 format (2(1x,a8),2(1x,a10),2(1x,a8))
 6010 format (2(1x,f8.1),2(1x,i10),2(1x,f8.3))
c
      return
      end
c
c
c
      subroutine dist_short (par)
c
      include 'moleman2.incl'
c
      real dist,qq,cut,q1,q2
c
      integer i,j,ierr,nd,nb,k
c
      character par*(*)
c
code ...
c
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      if (nselect .lt. 3) then
        call errcon ('Not enough atoms')
        return
      end if
c
      call str2r (par,cut,ierr)
      if (ierr .ne. 0) return
      call fvalut (' Cut-off (A) :',1,cut)
      call fvalut (' Ditto, bonded atoms (A) :',1,mxbond)
      if (cut .le. mxbond) then
        call errcon ('Invalid cut-off')
        return
      end if
c
      nd = 0
      nb = 0
      do i=1,natoms-1
        if (select(i)) then
          do j=i+1,natoms
            if (select(j)) then
c              if (iresid(i) .eq. iresid(j)) goto 1000
c              if (iresid(i) .eq. (iresid(j)+1) ) goto 1000
c              if (iresid(i) .eq. (iresid(j)-1) ) goto 1000
              nd = nd + 1
              if (abs(xyz(1,i)-xyz(1,j)) .gt. cut) goto 1000
              if (abs(xyz(2,i)-xyz(2,j)) .gt. cut) goto 1000
              if (abs(xyz(3,i)-xyz(3,j)) .gt. cut) goto 1000
              qq = dist (i,j,xyz)
              if (qq .gt. mxbond .and. qq .le. cut) then
c
c ... check if 1-3 interaction
c
                do k=atmptr(1,resptr(i)),atmptr(2,resptr(i))
                  if (k .ne. i) then
                    q1 = dist (j,k,xyz)
                    if (q1 .le. mxbond) then
                      q2 = dist (i,k,xyz)
                      if (q2 .le. mxbond) goto 1000
                    end if
                  end if
                end do
c
                if (iresid(i) .ne. iresid(j)) then
                  do k=atmptr(1,resptr(j)),atmptr(2,resptr(j))
                    if (k .ne. j) then
                      q1 = dist (i,k,xyz)
                      if (q1 .le. mxbond) then
                        q2 = dist (j,k,xyz)
                        if (q2 .le. mxbond) goto 1000
                      end if
                    end if
                  end do
                end if
c
c ... looks like a short contact (or an S=S link ;-)
c
                nb = nb + 1
                write (*,6000) qq
                call print_atom (i)
                call print_atom (j)
              end if
 1000         continue
            end if
          end do
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of distances      :',1,nd)
      call jvalut (' Nr of short contacts :',1,nb)
c
 6000 format (/' Short contact: ',f8.2,' A')
c
      return
      end
c
c
c
      subroutine dist_chain (chna,chnb,par)
c
      include 'moleman2.incl'
c
      real dist,qq,cut
c
      integer i,j,ierr,nb
c
      character chna*(*),chnb*(*),par*(*)
c
code ...
c
      call upcase (chna)
      call upcase (chnb)
      call textut (' Chain 1 :',chna(1:1))
      call textut (' Chain 2 :',chnb(1:1))
c
      call str2r (par,cut,ierr)
      if (ierr .ne. 0) return
      call fvalut (' Cut-off (A) :',1,cut)
c
      nb = 0
      do i=1,natoms
        if (achain(i) .eq. chna(1:1)) nb = nb + 1
      end do
      call jvalut (' Nr of atoms in chain 1 :',1,nb)
      if (nb .lt. 1) return
c
      nb = 0
      do i=1,natoms
        if (achain(i) .eq. chnb(1:1)) nb = nb + 1
      end do
      call jvalut (' Nr of atoms in chain 2 :',1,nb)
      if (nb .lt. 1) return
c
      nb = 0
      do i=1,natoms
        if (achain(i) .eq. chna(1:1)) then
          do j=1,natoms
            if (achain(j) .eq. chnb(1:1)) then
              if (abs(xyz(1,i)-xyz(1,j)) .gt. cut) goto 1000
              if (abs(xyz(2,i)-xyz(2,j)) .gt. cut) goto 1000
              if (abs(xyz(3,i)-xyz(3,j)) .gt. cut) goto 1000
              qq = dist (i,j,xyz)
              if (qq .le. cut) then
                nb = nb + 1
                write (*,6000) qq
                call print_atom (i)
                call print_atom (j)
              end if
 1000         continue
            end if
          end do
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of contacts :',1,nb)
c
 6000 format (/' Contact: ',f8.2,' A')
c
      return
      end
c
c
c
      subroutine dist_list (par1,par2)
c
      include 'moleman2.incl'
c
      real dist,qq,cutlo,cuthi,sum,xmin,xmax
c
      integer i,j,ierr,nd,nb
c
      character par1*(*),par2*(*)
c
code ...
c
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      if (nselect .lt. 3) then
        call errcon ('Not enough atoms')
        return
      end if
c
      call str2r (par1,cutlo,ierr)
      if (ierr .ne. 0) return
      call str2r (par2,cuthi,ierr)
      if (ierr .ne. 0) return
      call rlohi (cutlo,cuthi)
      call fvalut (' Lower cut-off :',1,cutlo)
      call fvalut (' Upper cut-off :',1,cuthi)
c
      nd = 0
      nb = 0
      sum = 0.0
c
      do i=1,natoms-1
        if (select(i)) then
          do j=i+1,natoms
            if (select(j)) then
              nd = nd + 1
              if (abs(xyz(1,i)-xyz(1,j)) .gt. cuthi) goto 1000
              if (abs(xyz(2,i)-xyz(2,j)) .gt. cuthi) goto 1000
              if (abs(xyz(3,i)-xyz(3,j)) .gt. cuthi) goto 1000
              qq = dist (i,j,xyz)
              if (qq .ge. cutlo .and. qq .le. cuthi) then
                nb = nb + 1
                write (*,6000) qq
                call print_atom (i)
                call print_atom (j)
                sum = sum + qq
                if (nb .eq. 1) then
                  xmin = qq
                  xmax = qq
                else
                  if (qq .lt. xmin) xmin = qq
                  if (qq .gt. xmax) xmax = qq
                end if
              end if
 1000         continue
            end if
          end do
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of distances checked :',1,nd)
      call jvalut (' Nr of distances listed  :',1,nb)
c
      if (nb .gt. 0) then
        sum = sum / float(nb)
        call fvalut (' Minimum distance :',1,xmin)
        call fvalut (' Maximum distance :',1,xmax)
        call fvalut (' Average distance :',1,sum)
      end if
c
 6000 format (/' Distance: ',f8.2,' A between:')
c
      return
      end
c
c
c
      subroutine dist_select (par,towhom)
c
      include 'moleman2.incl'
c
      real dist,qq,cut
c
      integer i,j,ierr,nd,nb,mode
c
      character par*(*),towhom*(*)
c
code ...
c
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      if (nselect .lt. 1) then
        call errcon ('No atoms selected')
        return
      end if
c
      call str2r (par,cut,ierr)
      if (ierr .ne. 0) return
      call fvalut (' Cut-off (A) :',1,cut)
      if (cut .le. 0.1) then
        call errcon ('Invalid cut-off')
        return
      end if
c
      call upcase (towhom)
      if (towhom(1:2) .eq. 'AL') then
        mode = 1
        call prompt (' Distances to ALl atoms')
      else if (towhom(1:2) .eq. 'SE') then
        mode = 2
        call prompt (' Distances to SElected atoms')
      else if (towhom(1:2) .eq. 'UN') then
        mode = 3
        call prompt (' Distances to UNselected atoms')
      else
        mode = 1
        call errcon ('Invalid mode; assuming AL')
        call prompt (' Distances to ALl atoms')
      end if
c
      nd = 0
      nb = 0
      do i=1,natoms
        if (select(i)) then
          do j=1,natoms
            if (i .eq. j) goto 1000
            if (mode .eq. 2 .and. (.not. select(j))) goto 1000
            if (mode .eq. 3 .and. select(j)) goto 1000
            nd = nd + 1
            if (abs(xyz(1,i)-xyz(1,j)) .gt. cut) goto 1000
            if (abs(xyz(2,i)-xyz(2,j)) .gt. cut) goto 1000
            if (abs(xyz(3,i)-xyz(3,j)) .gt. cut) goto 1000
            qq = dist (i,j,xyz)
            if (qq .le. cut) then
              nb = nb + 1
              write (*,6000) qq
              call print_atom (i)
              call print_atom (j)
            end if
 1000       continue
          end do
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of distances :',1,nd)
      call jvalut (' Nr listed       :',1,nb)
c
 6000 format (/' Distance: ',f8.2,' A')
c
      return
      end
c
c
c
      subroutine vtrace (vrdist,cavrml)
c
      include 'moleman2.incl'
c
      real vrdist,dist
c
      integer nca,ntot,k,i
c
      character cavrml*4
c
code ...
c
      nca = 0
      ntot = 0
      k = -1
c
      do i=1,natoms
        if (select(i)) then
          if (atmnam(i) .eq. cavrml) then
c
            if (nca .gt. 0) then
              if (dist(i,k,xyz) .gt. vrdist) then
                call xvrml_polyline (nca,rbuf)
                nca = 0
                k = -1
              end if
            end if
c
            nca = nca + 1
            ntot = ntot + 1
            rbuf(3*(nca-1)+1) = xyz(1,i)
            rbuf(3*(nca-1)+2) = xyz(2,i)
            rbuf(3*(nca-1)+3) = xyz(3,i)
            k = i
          end if
        end if
      end do
c
      if (nca .gt. 0) then
        call xvrml_polyline (nca,rbuf)
      end if
c
      call jvalut (' Nr of central atoms written :',1,ntot)
c
      return
      end
c
c
c
      subroutine vcpk (lmono,xrad)
c
      include 'moleman2.incl'
c
      real xrad
c
      integer ntot,i
c
      logical lmono
c
code ...
c
      ntot = 0
c
      do i=1,natoms
        if (select(i)) then
          ntot = ntot + 1
          rbuf(3*(ntot-1)+1) = xyz(1,i)
          rbuf(3*(ntot-1)+2) = xyz(2,i)
          rbuf(3*(ntot-1)+3) = xyz(3,i)
          rbuf(3*maxatm+ntot) = xrad
ccc          if (cvbrad(i).gt.0.0) rbuf(3*maxatm+ntot) = cvbrad(i)
          if (.not. lmono) then
            ibuf (ntot) = colour (i)
          end if
        end if
      end do
c
      if (ntot .gt. 0) then
        if (lmono) then
          call xvrml_cpk (ntot,rbuf,rbuf(3*maxatm+1))
        else
          call xvrml_cpk_col (ntot,rbuf,rbuf(3*maxatm+1),ibuf)
        end if
      end if
c
      call jvalut (' Nr of atoms written :',1,ntot)
c
      return
      end
c
c
c
      subroutine vstick (lball,lcyl,lwire,lmono,strad)
c
      include 'moleman2.incl'
c
      real cut,dist,strad
c
      integer ntot,i,j,k,ip,jp,maxbnd,nbo
c
      logical lmono,lball,lcyl,lwire
c
code ...
c
      if (.not. (lball.or.lcyl.or.lwire)) return
c
      ntot = 0
      nbo = 0
      maxbnd = 4
c
c ... save current radii; get real covalent bond radii
c
      do i=1,natoms
        rbuf (4*maxatm+i) = cvbrad (i)
      end do
      call pdb_chem_charge (.false.,.false.)
c
c ... which atoms to draw ?
c
      do i=1,natoms
        if (select(i)) then
          ntot = ntot + 1
          rbuf(3*(ntot-1)+1) = xyz(1,i)
          rbuf(3*(ntot-1)+2) = xyz(2,i)
          rbuf(3*(ntot-1)+3) = xyz(3,i)
          rbuf(3*maxatm+ntot) = 1.0
          if (rbuf (4*maxatm+i).gt.0.0)
     +      rbuf(3*maxatm+ntot) = rbuf (4*maxatm+i)
          jbuf (ntot) = i
          if (.not. lmono) then
            jbuf (maxatm+ntot) = colour (i)
          end if
        end if
      end do
c
      if (ntot .gt. 0) then
c
        do i=1,ntot
          ibuf (i) = 0
          do j=1,maxbnd
            ibuf (maxatm+(i-1)*maxbnd+j) = 0
          end do
        end do
c
c ... generate neighbour lists
c
        do i=1,ntot-1
          ip = jbuf(i)
          do j=i+1,ntot
            jp = jbuf(j)
            cut = cvbrad(ip) + cvbrad(jp) + 0.45
            if (dist(ip,jp,xyz) .le. cut) then
              nbo = nbo + 1
              if (ibuf(i) .lt. maxbnd) then
                ibuf (i) = ibuf (i) + 1
                ibuf (maxatm+(i-1)*maxbnd+ibuf(i)) = j
              else if (ibuf(j) .lt. maxbnd) then
                ibuf (j) = ibuf (j) + 1
                ibuf (maxatm+(j-1)*maxbnd+ibuf(j)) = i
              else
                call errcon ('Too many bonds !!!')
              end if
            end if
          end do
        end do
c
c ... flag atoms without any bonds at all (draw as 3D crosses)
c
        do i=1,ntot
          if (ibuf(i) .eq. 0) then
            do j=1,ntot
              if (i.ne.j .and. ibuf(j) .gt. 0) then
                do k=1,ibuf(j)
                  if (i.eq.ibuf(maxatm+(j-1)*maxbnd+k)) then
                    ibuf(i) = -1
                    goto 123
                  end if
                end do
              end if
            end do
  123       continue
          end if
        end do
c
        ip = 0
        jp = 0
        do i=1,ntot
          if (ibuf(i) .eq. 0) ip = ip + 1
          if (ibuf(i) .lt. 0) jp = jp + 1
        end do
        call jvalut (' Atoms without bonds :',1,ip)
        call jvalut (' Terminal atoms      :',1,jp)
c
        call jvalut (' Nr of atoms :',1,ntot)
        call jvalut (' Nr of bonds :',1,nbo)
c
        if (lball) then
          if (lmono) then
            call xvrml_cpk (ntot,rbuf,rbuf(3*maxatm+1))
          else
            call xvrml_cpk_col (ntot,rbuf,rbuf(3*maxatm+1),
     +             jbuf(maxatm+1))
          end if
        end if
c
        if (lcyl) then
          if (lmono) then
            call xvrml_stick (ntot,rbuf,maxbnd,ibuf(1),
     +             ibuf(maxatm+1),strad)
          else
            call xvrml_stick_col (ntot,rbuf,maxbnd,ibuf(1),
     +             ibuf(maxatm+1),strad,jbuf(maxatm+1))
          end if
        end if
c
        if (lwire) then
          if (lmono) then
            call xvrml_wire (ntot,rbuf,maxbnd,ibuf(1),
     +             ibuf(maxatm+1))
          else
            call xvrml_wire_col (ntot,rbuf,maxbnd,ibuf(1),
     +             ibuf(maxatm+1),jbuf(maxatm+1))
          end if
        end if
c
      else
        call errcon ('No atoms selected !')
      end if
c
c ... restore user's current radii
c
      do i=1,natoms
        cvbrad (i) = rbuf (4*maxatm+i)
      end do
c
      return
      end
c
c
c
      subroutine vturd (lmono,vrdist,cavrml,fatness)
c
      include 'moleman2.incl'
c
      real vrdist,dist,fatness
c
      integer nca,ntot,k,i
c
      logical lmono
c
      character cavrml*4
c
code ...
c
      nca = 0
      ntot = 0
      k = -1
c
      do i=1,natoms
        if (select(i)) then
          if (atmnam(i) .eq. cavrml) then
c
            if (nca .gt. 0) then
              if (dist(i,k,xyz) .gt. vrdist) then
                ibuf (nca) = 0
                if (lmono) then
                  call xvrml_cpk (nca,rbuf,rbuf(3*maxatm+1))
                  call xvrml_stick (nca,rbuf,1,ibuf(1),
     +               ibuf(maxatm+1),fatness)
                else
                  call xvrml_cpk_col (nca,rbuf,rbuf(3*maxatm+1),jbuf)
                  call xvrml_stick_col (nca,rbuf,1,ibuf(1),
     +               ibuf(maxatm+1),fatness,jbuf)
                end if
                nca = 0
                k = -1
              end if
            end if
c
            nca = nca + 1
            ntot = ntot + 1
            rbuf (3*(nca-1)+1) = xyz(1,i)
            rbuf (3*(nca-1)+2) = xyz(2,i)
            rbuf (3*(nca-1)+3) = xyz(3,i)
            rbuf (3*maxatm+nca) = fatness
            ibuf (nca) = 1
            ibuf (maxatm+nca) = nca + 1
            jbuf (nca) = colour (i)
            k = i
          end if
        end if
      end do
c
      if (nca .gt. 0) then
        ibuf (nca) = 0
        if (lmono) then
          call xvrml_cpk (nca,rbuf,rbuf(3*maxatm+1))
          call xvrml_stick (nca,rbuf,1,ibuf(1),
     +      ibuf(maxatm+1),fatness)
        else
          call xvrml_cpk_col (nca,rbuf,rbuf(3*maxatm+1),jbuf)
          call xvrml_stick_col (nca,rbuf,1,ibuf(1),
     +      ibuf(maxatm+1),fatness,jbuf)
        end if
      end if
c
      call jvalut (' Nr of central atoms written :',1,ntot)
c
      return
      end
c
c
c
      subroutine vcramp (how)
c
      include 'moleman2.incl'
c
      real r,g,b,h1,h2,h3,v1,v2,s1,s2,delta,xmin,xmax,h
c
      integer ntot,i,k,idum
c
      character how*2
c
code ...
c
      ntot = 0
c
      if (how .eq. 'BF') then
        call prompt (' Colour-ramping by B-factor ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = batom (i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'OC') then
        call prompt (' Colour-ramping by occupancy ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = qatom (i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'MA') then
        call prompt (' Colour-ramping by atomic mass ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = atmass (i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'RA') then
        call prompt (' Colour-ramping by atomic radius ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = cvbrad (i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'EL') then
        call prompt (' Colour-ramping by element number ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = element (i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'AT') then
        call prompt (' Colour-ramping by atom number ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = atomnr (i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'X ') then
        call prompt (' Colour-ramping by X-coordinate ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = xyz (1,i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'Y ') then
        call prompt (' Colour-ramping by Y-coordinate ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = xyz (2,i)
            ibuf (ntot) = i
          end if
        end do
      else if (how .eq. 'Z ') then
        call prompt (' Colour-ramping by Z-coordinate ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = xyz (3,i)
            ibuf (ntot) = i
          end if
        end do
      else
        call prompt (' Colour-ramping by residue number ...')
        do i=1,natoms
          if (select(i)) then
            ntot = ntot + 1
            rbuf (ntot) = iresid (i)
            ibuf (ntot) = i
          end if
        end do
      end if
c
      call jvalut (' Nr of atoms selected :',1,ntot)
      if (ntot .le. 1) then
        call errcon ('Not enough atoms for ramping')
        return
      end if
      xmin = rbuf (1)
      xmax = rbuf (1)
      do i=1,ntot
        if (rbuf(i) .lt. xmin) xmin = rbuf (i)
        if (rbuf(i) .gt. xmax) xmax = rbuf (i)
      end do
      call rvalut (' Lowest  value found :',1,xmin)
      call rvalut (' Highest value found :',1,xmax)
      if (xmin .ge. xmax) then
        call errcon ('Invalid range of values for ramping')
        return
      end if
c
      call rgbhsv (0.0,0.0,1.0,h1,s1,v1)
      call rgbhsv (1.0,0.0,0.0,h2,s2,v2)
      if (h1 .eq. h2) h2 = h2 + 180.0
      h3 = (h2 - h1) / 100.0
      delta = (h2-h1) / (xmax-xmin)
c
      do i=1,ntot
        k = int(float(int((rbuf(i)-xmin)*delta/h3))*h3 + h1) 
        h = float(k)
        if (h2 .gt. h1) then
          h = max (h1, min (h2, h))
        else
          h = max (h2, min (h1, h))
        end if
        call fix360 (h)
        call hsvrgb (h,1.0,1.0,r,g,b)
        call xvrml_encode_rgb (r,g,b,idum)
        colour (ibuf(i)) = idum
c
ccc        write (*,'(i6,1pe12.4,i15,i8)') i,rbuf(i),idum,ibuf(i)
c
      end do
c
      call prompt (' Ramping done')
c
      return
      end
