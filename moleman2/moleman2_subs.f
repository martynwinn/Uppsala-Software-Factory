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
c
c ... mole2_gen.f - GENERAL subroutines for MOLEMAN2
c
c ... i.e., those that do not include 'moleman2.incl'
c
      subroutine detaja (aresid,achain,iresid,insert)
c
c ... deduce chain id and residue nr and insert from an Alwyn-type
c     character residue id [i.e., go from (A6) to (A1,I4,A1)]
c
c     ARESID = 'AA132'  -> ACHAIN = 'A' IRESID =  132 INSERT = 'A'
c              'W43'    ->          'W'            43          ' '
c              '191'    ->          ' '           191          ' '
c              '12A'    ->          'A'            12          ' '
c              '34AA'   ->          'A'            34          'A'
c              'X1292A' ->          'X'          1292          'A'
c
c ... Gerard Kleywegt @ 951108
c
      implicit none
c
      integer iresid,ll,length
c
      character aresid*6,achain*1,insert*1,esline*128
c
code ...
c
c ... remove spaces and left-justify text
c
      call remspa (aresid)
      ll = length (aresid)
c
c ... first guess: no chain id/no insert (legal PDB)
c
      achain = ' '
      insert = ' '
      read (aresid(1:6),*,err=10) iresid
      if (iresid .gt. 9999) then
        write (aresid,'(i6)') iresid
        achain = aresid(1:2)
        read (aresid(3:6),*) iresid
      end if
      return
c
c ... second guess: single-character chain id (legal PDB)
c
   10 continue
      achain = aresid(1:1)
      insert = ' '
      read (aresid(2:6),*,err=20) iresid
      return
c
c ... third guess: two-character chain id
c     (this is an extension to "pure PDB")
c
   20 continue
      achain = aresid(1:1)
      insert = aresid(2:2)
      read (aresid(3:6),*,err=30) iresid
      return
c
c ... fourth guess: one-character ID APPENDED to residue nr
c
   30 continue
      achain = aresid (ll:ll)
      insert = ' '
      read (aresid(1:ll-1),*,err=40) iresid
      return
c
c ... fifth guess: two-character ID APPENDED to residue nr
c     (this is an extension to "pure PDB")
c
   40 continue
      achain = aresid (ll-1:ll-1)
      insert = aresid (ll:ll)
      read (aresid(1:ll-2),*,err=50) iresid
      return
c
c ... final guess: "A123B"-type
c     (this is an extension to "pure PDB")
c
   50 continue
      achain = aresid (1:1)
      insert = aresid(ll:ll)
      read (aresid(2:ll-1),*,err=60) iresid
      return
c
c ... if here, then PDB convention grossly violated
c
   60 continue
      esline = 'While converting Alwyn-type residue id ('
     +         // aresid // ')'
      call errcon (esline)
      achain = '?'
      insert = ' '
      iresid = 9999
c
      return
      end
c
c
c
      subroutine suggot (nat,nam,xyz,ip,newxyz,lprint)
c
      implicit NONE
c
      real codist,ocoang
      parameter (codist=1.23, ocoang=122.5)
c
      integer nat
c
      real xyz(3,nat),otxyz(3,5),otint(3,5),newxyz(3)
      real dist,angle,tangle
c
      integer otref(3,5),iot(5)
      integer i,j,ip
c
      logical lprint
c
      character nam(nat)*4
c
code ...
c
      do i=1,5
        iot(i) = 0
        do j=1,3
          otxyz (j,i) = 0.0
          otint (j,i) = 0.0
          otref (j,i) = 0
        end do
      end do
c
      do i=1,nat
        if (nam(i) .eq. ' N  ') then
          iot (1) = i
        else if (nam(i) .eq. ' CA ') then
          iot (2) = i
        else if (nam(i) .eq. ' C  ') then
          iot (3) = i
        else if (nam(i) .eq. ' O  ') then
          iot (4) = i
        else if (nam(i) .eq. ' OT1') then
          iot (4) = i
        else if (nam(i) .eq. ' OT2') then
          iot (5) = i
        else if (nam(i) .eq. ' OTX') then
          iot (5) = i
        else if (nam(i) .eq. ' OXT') then
          iot (5) = i
        end if
      end do
c
      if (iot(1) .le. 0) then
        if (lprint) call errcon (' Atom " N  " not found')
        return
      end if
c
      if (iot(2) .le. 0) then
        if (lprint) call errcon (' Atom " CA " not found')
        return
      end if
c
      if (iot(3) .le. 0) then
        if (lprint) call errcon (' Atom " C  " not found')
        return
      end if
c
      j = 3
      if (iot(4) .gt. 0) j = 4
      if (j .eq. 4 .and. iot(5) .gt. 0) j = 5
      do i=1,j
        otxyz (1,i) = xyz(1,iot(i))
        otxyz (2,i) = xyz(2,iot(i))
        otxyz (3,i) = xyz(3,iot(i))
      end do
c
      otint (1,2) = dist(1,2,otxyz)
      otref (1,2) = 1
cc        call fvalut (' Distance N-CA          =',1,otint(1,2))
c
      otint (1,3) = dist (2,3,otxyz)
      otref (1,3) = 2
cc        call fvalut (' Distance CA-C          =',1,otint(1,3))
      otint (2,3) = angle (1,2,3,otxyz)
      otref (2,3) = 1
cc        call fvalut (' Angle N-CA-C           =',1,otint(2,3))
c
      otint (1,4) = codist
      otref (1,4) = 3
cc        call fvalut (' Distance C-OT1         =',1,
cc     +    dist(3,4,otxyz))
cc        call fvalut (' Reset to               =',1,otint(1,4))
c ... keep CA-COO flat
      otint (2,4) = 0.5 * (360.0 - ocoang)
      otref (2,4) = 2
cc        call fvalut (' Angle CA-C-OT1         =',1,
cc     +    angle(2,3,4,otxyz))
cc        call fvalut (' Reset to               =',1,otint(2,4))
      otint (3,4) = tangle (1,2,3,4,otxyz)
      otref (3,4) = 1
cc        call fvalut (' Dihedral N-CA-C-OT1    =',1,otint(3,4))
c
      otint (1,5) = codist
      otref (1,5) = 3
cc        call fvalut (' Distance C-OT2         =',1,otint(1,5))
      otint (2,5) = ocoang
      otref (2,5) = 4
cc        call fvalut (' Angle OT1-C-OT2        =',1,otint(2,5))
      otint (3,5) = 180.0
      otref (3,5) = 2
cc        call fvalut (' Dihedral CA-OT1-C-OT2  =',1,otint(3,5))
c
cc        call fvalut (' Int coords :',15,otint)
cc        call ivalut (' Ref coords :',15,otref)
cc        call fvalut (' XYZ coords :',15,otxyz)
c
      call c2car4 (5,otxyz,otint,otref)
c
cc        write (*,*)
cc        write (*,6374) ' N  ',xf(iot(1)),yf(iot(1)),zf(iot(1)),
cc     +    (otxyz(i,1),i=1,3)
cc        write (*,6374) ' CA ',xf(iot(2)),yf(iot(2)),zf(iot(2)),
cc     +    (otxyz(i,2),i=1,3)
cc        write (*,6374) ' C  ',xf(iot(3)),yf(iot(3)),zf(iot(3)),
cc     +    (otxyz(i,3),i=1,3)
c
      if (lprint) then
        if (iot(4) .gt. 0) then
          write (*,6374) ' OT1',(xyz(i,iot(4)),i=1,3),
     +      (otxyz(i,4),i=1,3)
        else
          write (*,6374) ' OT1',0.0,0.0,0.0,(otxyz(i,4),i=1,3)
        end if
c
        if (iot(5) .gt. 0) then
          write (*,6374) ' OT2',(xyz(i,iot(5)),i=1,3),
     +      (otxyz(i,5),i=1,3)
        else
          write (*,6374) ' OT2',0.0,0.0,0.0,(otxyz(i,5),i=1,3)
        end if
c
cc      write (*,*)
cc      write (*,*) 'Check geometry of carboxylate group :'
cc      call fvalut (' Dist       C-OT1 =',1,dist(3,4,otxyz))
cc      call fvalut (' Dist       C-OT2 =',1,dist(3,5,otxyz))
cc      call fvalut (' Dist     OT1-OT2 =',1,dist(4,5,otxyz))
cc      call fvalut (' Angle  OT1-C-OT2 =',1,angle(4,3,5,otxyz))
cc      call fvalut (' Angle   OT1-C-CA =',1,angle(4,3,2,otxyz))
cc      call fvalut (' Angle   OT2-C-CA =',1,angle(5,3,2,otxyz))
cc      call fvalut (' Dih   OT1-C-CA-N =',1,tangle(4,3,2,1,otxyz))
cc      call fvalut (' Dih   OT2-C-CA-N =',1,tangle(5,3,2,1,otxyz))
cc      call fvalut (' Dih CA-OT1-C-OT2 =',1,tangle(2,4,3,5,otxyz))
c
        call prompt (' You must add/edit OT1/OT2 yourself !')
      end if
c
      ip = iot(5)
      newxyz(1) = otxyz(1,5)
      newxyz(2) = otxyz(2,5)
      newxyz(3) = otxyz(3,5)
c
cc        write (*,*)
c
 6374 format (' ==> ',a4,' NOW : ',3f8.3,' SUGGESTED : ',3f8.3)
c
      return
      end
c
c
c
      subroutine find_type (nat,resnam,atmnam,myres,myatm,
     +  nok,iptr,lprint)
c
      implicit none
c
      integer nat,nok,i
      integer iptr(nat)
c
      logical lprint
c
      character resnam(nat)*3,atmnam(nat)*4,myres*3,myatm*4
c
code ...
c
      if (lprint) write (*,6000) myres,myatm
 6000 format (' Looking for ',a3,'-',a4,' atoms ...')
c
      nok = 0
c
      do i=1,nat
        if (resnam(i) .eq. myres) then
          if (atmnam(i) .eq. myatm) then
            nok = nok + 1
            iptr (nok) = i
            if (lprint) call print_atom (i)
          end if
        end if
      end do
c
      return
      end
c
c
c
      subroutine odicts (myres,iopt,f1,napr,otxyz,line,bndcut,
     +  maxtor,dismat,valtor,deftor,isbond,okbond,afftor,tornam,
     +  atmnam,dummy)
c
      integer f1,napr,iopt,maxtor
c
      real otxyz(3,napr),dismat(napr,napr),valtor(maxtor)
      real dist,tangle,dummy,xx,bndcut,xdmin
c
      integer deftor(4,maxtor)
      integer i,j,k,l,nl,i1,j1,j2,ncnt,idmin
      integer length,leng1
c
      logical isbond(napr,napr),okbond(napr,napr)
      logical afftor(napr,maxtor)
      logical lhydro
c
      character tornam(maxtor)*6,atmnam(napr)*4
      character line*(*)
      character line2*256,junky*256,myres*3
c
code ...
c
      call fvalut (' Cut-off distance for bonded atoms :',
     +      1,bndcut)
      call fvalut (' Special torsion-angle tolerance   :',
     +      1,dummy)
c
cccc      call textut (' Residue type :',myres)
cccc      call textut (' Atom types   :',line)
c
      if (iopt .eq. 3 .or. iopt .eq. 4) goto 2313
c
c ... rs-fit or rsr datablock file
c
      call stamp (line2)
      write (f1,'(a1,1x,a)') '!',line2(1:leng1(line2))
      write (f1,'(a1)') '!'
      if (iopt .eq.1) then        
        write (f1,'(a,1x,a)') '! RS-fit datablock for',myres
      else if (iopt .eq. 2) then
        write (f1,'(a,1x,a)') '! RSR datablock for',myres
      end if
      write (f1,'(a1)') '!'
c
      if (iopt .eq.1) then        
        junky = 'rsfit_'//myres
      else if (iopt .eq. 2) then
        junky = 'RSR_dict_'//myres
      end if
      call remspa (junky)
      call textut (' Datablock name :',junky)
      nl = 1 + ((napr*5)/71)
      write (line2,'(a,1x,a,1x,i3,1x,i3)')
     +  junky(1:leng1(junky)),'T',nl,70
      call pretty (line2)
      write (f1,'(a)') line2(1:leng1(line2))
c
      i1 = 0
      do i=1,napr,14
        j = min (napr,i+13)
        write (f1,'(a)') line(1+5*(i-1):4+5*(j-1))
        i1 = i1 + 1
      end do
c
      if (i1 .lt. nl) then
        do i=i1+1,nl
          write (f1,*)
        end do
      end if
c
      close (f1)
      call prompt (' Datablock written')
      return
c
c ... connect file or torsion datablock
c
 2313 continue
c
      nl = 0
      do i=1,napr
        do j=i+1,napr
          isbond (i,j) = .false.
          isbond (j,i) = .false.
          okbond (i,j) = .false.
          okbond (j,i) = .false.
c
          dismat (i,j) = dist (i,j,otxyz)
          dismat (j,i) = dismat (i,j)
c
          if (dismat(i,j) .le. bndcut) then
            if (.not. (lhydro(atmnam(i)) .and.
     +                 lhydro(atmnam(j)) )) then
              isbond (i,j) = .true.
              isbond (j,i) = .true.
              nl = nl + 1
            end if
          end if
        end do
        isbond (i,i) = .false.
        okbond (i,i) = .false.
        dismat (i,i) = 0.0
      end do
      call jvalut (' Nr of bonds :',1,nl)
c
c ... 981014 - hydrogens only bonded to nearest neighbour !
c
      do i=1,napr
        if (lhydro(atmnam(i))) then
          xdmin = 999.9
          idmin = -1
          do j=1,napr
            if (isbond(i,j)) then
              if (dismat(i,j) .lt. xdmin) then
                xdmin = dismat(i,j)
                idmin = j
              end if
            end if
          end do
          if (idmin .gt. 0) then
            do j=1,napr
              if (j.ne.idmin .and. isbond(i,j)) then
                isbond (i,j) = .false.
                isbond (j,i) = .false.
                nl = nl - 1
              end if
            end do
          end if
        end if
      end do
c
      call jvalut (' Nr of bonds :',1,nl)
c
      if (iopt .eq. 4) then
        if (nl .le. 0) then
          call errcon (' No bonds means no torsions')
          return
        end if
        goto 2452
      end if
c
c ... connect file
c
      j2 = 0
      write (f1,'(a3)') myres
      do i=1,napr,10
        j = min (napr,i+9)
        write (f1,'(a4,1x,a)') 'ATOM',line(1+5*(i-1):4+5*(j-1))
        j2 = j2 + 1
      end do
c
      if (nl .gt. 0) then
c
        j1 = 0
c
 2336   continue
c
        i1 = 0
c
        do i=1,napr
          do j=1,napr
            if ( isbond(i,j) .and. (.not. okbond(i,j)) ) then
              j1 = j1 + 1
              if (j1 .eq. 1) then
                line2 = 'CONNECT - '//line(1+5*(i-1):4+5*(i-1))//
     +                  ' '//line(1+5*(j-1):4+5*(j-1))
              else
                line2 = 'CONNECT '//line(1+5*(i-1):4+5*(i-1))//
     +                  ' '//line(1+5*(j-1):4+5*(j-1))
              end if
              i1 = 2
              k = j
              okbond (i,j) = .true.
              okbond (j,i) = .true.
              goto 2339
            end if
          end do
        end do
        goto 2341
c
 2339   continue
        if (i1 .ge. 12) goto 2343
        do i=1,napr
          if (isbond(k,i) .and. (.not. okbond(k,i)) ) then
            line2 = line2(1:leng1(line2))//' '//
     +              line(1+5*(i-1):4+5*(i-1))
            i1 = i1 + 1
            okbond (i,k) = .true.
            okbond (k,i) = .true.
            k = i
            goto 2339
          end if
        end do
c
 2343   continue
        if (j1 .eq. 1) then
          line2 = line2(1:leng1(line2))//' +'
        end if
        write (f1,'(a)') line2(1:leng1(line2))
        goto 2336
c
 2341   continue
c
c ... check orphan atoms without any bonds
c
        do i=1,napr
          do j=1,napr
            if (isbond(i,j)) goto 2424
          end do
          write (f1,'(a7,1x,a,1x,a)') 'CONNECT',
     +      line(1+5*(i-1):4+5*(i-1)),line(1+5*(i-1):4+5*(i-1))
          j1 = j1 + 1
 2424     continue
        end do
c
      else
c
c ... no bonds at all; just "connect" individual atoms with themselves
c
        do i=1,napr
          write (f1,'(a7,1x,a,1x,a)') 'CONNECT',
     +      line(1+5*(i-1):4+5*(i-1)),line(1+5*(i-1):4+5*(i-1))
        end do
        j1 = napr
      end if
c
      close (f1)
      call jvalut (' Nr of ATOM records :',1,j2)
      call jvalut (' Nr of CONNECT ,,   :',1,j1)
      call prompt (' Connect file written (append to all.dat)')
      return
c
c ... torsion datablock
c
 2452 continue
c
      ncnt = 0
      do i=1,napr
        do j=1,napr
          if (i.eq.j) goto 2465
          if (dismat(i,j) .gt. bndcut) goto 2465
          do k=1,napr
            if (k.eq.j .or. k.eq.i) goto 2466
            if (dismat(i,k) .gt. bndcut) goto 2466
            do l=1,napr
              if (l.eq.j .or. l.eq.i .or. l.eq.k) goto 2467
              if (dismat(l,j) .gt. bndcut) goto 2467
c
c ... we have found a dihedral
c
              xx = tangle(k,i,j,l,otxyz)
              call fixang (xx)
              write (*,'(/1x,a8,4(1x,a4),1x,f8.2)')
     +          'DIHEDRAL',atmnam(k),atmnam(i),
     +          atmnam(j),atmnam(l),xx
c
c ... skip if torsion near 0 or +/- 180 degrees
c
              if (abs (xx) .le. dummy) then
                call prompt (' Skip -> torsion ~ 0')
                goto 2467
              else if (abs(xx-180) .le. dummy .or.
     +                 abs(xx+180) .le. dummy) then
                call prompt (' Skip -> torsion ~ 180')
                goto 2467
              end if
c
c ... get all atoms affected by the torsion
c
              do i1=1,napr
                do j1=1,napr
                  okbond (i1,j1) = .false.
                end do
              end do
c
c ... flag all atoms connected to atom "J", except "I", as affected
c
              do i1=1,napr
                if (i1.ne.i) then
                  if (isbond(j,i1)) okbond (i1,i1) = .true.
                end if
              end do
c
 2514         continue
              j2 = 0
              do i1=1,napr
                if (okbond(i1,i1)) then
                  do j1=1,napr
                    if (j1.ne.j) then
                      if (isbond(i1,j1) .and.
     +                    (.not. okbond(j1,j1))) then
                         okbond (j1,j1) = .true.
                         j2 = j2 + 1
                      end if
                    end if
                  end do
                end if
              end do
c
c ... if I or K are flagged, we probably have a ring torsion
c
              if (okbond(i,i) .or. okbond(k,k)) then
                call prompt (' Skip -> ring torsion')
                goto 2467
              end if
c
              if (j2 .gt. 0) goto 2514
c
              j2 = 1
              line2 = atmnam (l)
              do i1=1,napr
                if (i1 .ne. l .and. okbond(i1,i1)) then
                  line2 = line2(1:leng1(line2))//' '//
     +                    atmnam(i1)
                  j2 = j2 + 1
                end if
              end do
c
c ... if almost all atoms affected, better to use it the other
c     way around
c
              call textut (' Affected atoms :',line2)
              if (j2 .ge. (napr-4) ) then
                call prompt (' Skip -> too many affected atoms')
                goto 2467
              end if
c
c ... check if it is a permutation of a previous torsion
c
              if (ncnt .ge. 1) then
                do i1=1,ncnt
                  if (
     +         (i.eq.deftor(2,i1) .and. j.eq.deftor(3,i1))) then
cccc     +           (j.eq.deftor(2,i1) .and. i.eq.deftor(3,i1))) then
                    do j1=1,napr
                      if (okbond(j1,j1) .neqv. afftor(j1,i1))
     +                  goto 2532
                    end do
                    call textut (' Permutation of :',tornam(i1))
                    call prompt (' Skip -> permutation')
                    goto 2467
 2532               continue
                  end if
                end do
              end if
c
              if (ncnt .ge. maxtor) then
                call jvalut (' Max nr of torsions :',1,maxtor)
                call errcon (' Too many torsion angles !!!')
                goto 2467
              end if
c
c ... we have a new torsion !!!
c
              ncnt = ncnt + 1
              write (tornam(ncnt),'(a3,i3)') 'TOR',ncnt
              call remspa (tornam(ncnt))
              call textut (' OKAY ==> new torsion :',tornam(ncnt))
c
              deftor (1,ncnt) = k
              deftor (2,ncnt) = i
              deftor (3,ncnt) = j
              deftor (4,ncnt) = l
c
              valtor (ncnt) = float(nint(xx))
c
              do i1=1,napr
                afftor(i1,ncnt) = okbond(i1,i1)
              end do
c
 2467         continue
            end do
 2466       continue
          end do
 2465     continue
        end do
      end do
c
      write (*,*)
      call ivalut (' Nr of unique rotatable torsions :',1,ncnt)
c
c ... now write the file
c
      write (f1,'(a,1x,a)') 'RESIDUE',myres
      nl = 1
c
      do i=1,ncnt
        write (line2,'(a,1x,a,1x,f8.0,4(1x,a4))') 'TORSION',
     +    tornam(i),valtor(i),atmnam(deftor(1,i)),
     +    atmnam(deftor(2,i)),atmnam(deftor(3,i)),
     +    atmnam(deftor(4,i))
        call pretty (line2)
c
        do j=1,napr
          if (afftor(j,i)) then
            if (length(line2) .ge. 65) then
              write (f1,'(a,1x,a1)') line2(1:leng1(line2)),'\\'
              nl = nl + 1
              line2 = '   '//atmnam(j)
            else
              line2 = line2(1:leng1(line2))//' '//atmnam(j)
            end if
            call pretty (line2)
          end if
        end do
        write (f1,'(a)') line2(1:leng1(line2))
        nl = nl + 1
      end do
c
      close (f1)
      call jvalut (' Nr of lines written :',1,nl)
      call prompt (' Torsion file written (append to torsion.o)')
      return
c
      end
c
c
c
      subroutine salpha (x)
c
c ... SALPHA - set up penta-residue ALPHA helix CA coordinates
c     (from O's alpha template)
c
      real x(3,5)
c
code ...
c
      x(1,1) = 3.633
      x(2,1) = 15.082
      x(3,1) = 31.410
c
      x(1,2) = 3.847
      x(2,2) = 16.332
      x(3,2) = 27.802
c
      x(1,3) = 4.017
      x(2,3) = 12.736
      x(3,3) = 26.460
c
      x(1,4) = 0.906
      x(2,4) = 11.724
      x(3,4) = 28.417
c
      x(1,5) = -1.006
      x(2,5) = 14.798
      x(3,5) = 27.217
c
      return
      end
c
c
c
      subroutine sbeta (x)
c
c ... SBETA - set up penta-residue BETA strand CA coordinates
c     (from O's beta template)
c
      real x(3,5)
c
code ...
c
      x(1,1) = 43.840
      x(2,1) = 28.467
      x(3,1) = -6.991
c
      x(1,2) = 45.504
      x(2,2) = 26.414
      x(3,2) = -4.284
c
      x(1,3) = 43.734
      x(2,3) = 23.423
      x(3,3) = -2.720
c
      x(1,4) = 44.607
      x(2,4) = 21.760
      x(3,4) = 0.604
c
      x(1,5) = 43.505
      x(2,5) = 18.374
      x(3,5) = 1.905
c
      return
      end
c
c
c
      subroutine slefth (x)
c
c ... SLEFTH - set up tetra-residue LEFT-Handed helix CA coordinates
c     (from GJK's lefth template)
c
      real x(3,4)
c
code ...
c
      x(1,1) = 0.000
      x(2,1) = 0.000
      x(3,1) = 0.000
c
      x(1,2) = -0.972
      x(2,2) = 0.114
      x(3,2) = -3.703
c
      x(1,3) = -4.115
      x(2,3) = 2.184
      x(3,3) = -2.993
c
      x(1,4) = -5.086
      x(2,4) = 0.000
      x(3,4) = 0.000
c
      return
      end
c
c
c
      subroutine psnini (iunit,psfile,prognm)
c
c ... initialise PostScript Duarte-Pyle plot
c
      implicit none
c
      real phi,psi
c
      integer iunit,i
c
      character psfile*(*),prognm*(*)
      character labx*80,laby*80
c
code ...
c
      call xps_init ()
      call xps_open (iunit,psfile,prognm)
      call xps_scale (0.0,360.0,0.0,360.0)
      call xps_stroke ()
c
      call xps_ps_comment ('Boxes for core regions')
      call xps_light_box (150.0,190.0,0.0,360.0)
      call xps_light_box (0.0,360.0,190.0,260.0)
c
      call xps_move (0.,0.)
      call xps_draw (360.,0.)
      call xps_draw (360.,360.)
      call xps_draw (0.,360.)
      call xps_draw (0.,0.)
c
c --- ETA axis ticks
c
      call xps_ps_comment ('ETA axis ticks')
      do 300 i=1,11
        phi = i*30.
        call xps_move (phi,0.)
        call xps_draw (phi,5.)
        call xps_move (phi,360.)
        call xps_draw (phi,355.)
300   continue
c
c --- THETA axis ticks
c
      call xps_ps_comment ('THETA axis ticks')
      do 310 i=1,11
        psi = i*30.
        call xps_move (0.,psi)
        call xps_draw (5.,psi)
        call xps_move (360.,psi)
        call xps_draw (355.,psi)
310   continue
c
      labx = 'ETA mapped to [0,360>'
      laby = 'THETA mapped to [0,360>'
      call xps_label (labx,laby)
      call xps_ps_comment ('Duarte-Pyle plot initialisation done')
c
      return
      end
c
c ... mole2_subs.f - SPECIFIC subroutines for MOLEMAN2
c
c ... i.e., subroutines that *DO* include 'moleman2.incl'
c
c
      subroutine readdb (iunit,ierr)
c
      include 'moleman2.incl'
c
      integer maxopt
      parameter (maxopt=50)
c
      integer nline,iunit,ierr,ii,jj,i,j,k,nopt,kk,length
c
      logical lmore,lmch
c
      character line*256,optpar(maxopt)*80,quote*1
c
      parameter (quote='''')
c
code ...
c
      nline  = 0
      nalias = 0
      nmrtyp = 0
      nmatyp = 0
      nlring = 0
      ierr = 0
c
      libcat(1)='PROT'
      libcat(2)='NUCL'
      libcat(3)='WATE'
      libcat(4)='META'
      libcat(5)='INOR'
      libcat(6)='CARB'
      libcat(7)='ORGA'
      libcat(8)='HETE'
c
 6000 format (A,A,A,A,A)
c
   10 continue
      read (iunit,6000,err=8000,end=9000) line
      nline = nline + 1
      if (line(1:1) .eq. '!') goto 10
      call upcase (line)
      if (line(1:3) .eq. 'REM') then
        call textut (' >',line(4:))
        goto 10
      end if
      if (line(1:3) .eq. 'RES') goto 20
      if (line(1:3) .eq. 'RNG') goto 120
      call errcon ('Invalid keyword (expected !, REM, RES, RNG)')
      call jvalut (' Line number    :',1,nline)
      call textut (' Offending line :',line)
      ierr = -1
      return
c
c ... new residue type
c
   20 continue
      if (nmrtyp .ge. mxrtyp) then
        call errcon ('Too many residue definitions; skipping rest')
        call jvalut (' Max allowed    :',1,mxrtyp)
        call jvalut (' Line number    :',1,nline)
        call textut (' Offending line :',line)
        goto 9000
      end if
      nmrtyp = nmrtyp + 1
      ii = nmrtyp
      nmrptr (1,ii) = -1
      nmrptr (2,ii) = -1
      nmaptr (1,ii) = -1
      nmaptr (2,ii) = -1
      lrname (ii) = line(5:7)
      lrdesc (ii) = line(8:)
      call pretty (lrdesc(ii))
      lrtype (ii) = ihete
      libolc (ii) = '?'
      jj = nmatyp
      kk = nalias
c
   30 continue
      read (iunit,6000,err=8000,end=9000) line
      nline = nline + 1
      if (line(1:1) .eq. '!') goto 30
      call upcase (line)
      if (line(1:3) .eq. 'REM') then
        call textut (' >',line(4:))
        goto 30
      end if
c
      if (line(1:3) .eq. 'TYP') then
        do i=1,nrlcat
          if (line(5:8) .eq. libcat(i)) lrtype (ii) = i
        end do
c
      else if (line(1:3) .eq. 'OLC') then
        libolc (ii) = line(5:5)
c
      else if (line(1:3) .eq. 'AKA') then
        call extrop (line(5:),nopt,maxopt,optpar,ierr)
        if (ierr .ne. 0) then
          call errcon ('While parsing AKA line')
          call jvalut (' Line number    :',1,nline)
          call textut (' Offending line :',line)
          ierr = -1
          return
        end if
        if ( (nalias+nopt) .gt. mxrtal) then
          call errcon ('Too many AKA aliases; rest skipped')
          call jvalut (' Max allowed    :',1,mxrtal)
          call jvalut (' Line number    :',1,nline)
          call textut (' Offending line :',line)
        else
          do i=1,nopt
            nalias = nalias + 1
            lalias (nalias) = optpar(i)(1:3)
          end do
        end if
c
      else if (line(1:3) .eq. 'MCH' .or.
     +         line(1:3) .eq. 'SCH') then
c
        lmore = .false.
        lmch  = (line(1:3) .eq. 'MCH')
        line  = line (4:)
   40   continue
        j = length (line)
        if (line(j:j) .eq. '-') then
          lmore = .true.
          line(j:) = ' '
        end if
   50   continue
        j = index (line,quote)
        if (j .le. 0) goto 60
        k = j + index (line(j+1:),quote)
        if (k .le. 0) k = j + 5
        if (nmatyp .ge. mxrtat) then
          call errcon ('Too many atom types')
          call jvalut (' Max allowed    :',1,mxrtat)
          call jvalut (' Line number    :',1,nline)
          call textut (' Offending line :',line)
          ierr = -1
          return
        end if
        nmatyp = nmatyp + 1
        lratom (nmatyp) = line(j+1:j+4)
        ismain (nmatyp) = lmch
ccc          call textut (' >',line)
ccc          print *,j,k,' |',lratom (nmatyp),'| '
        line = line(k+1:)
        goto 50
c
   60   continue
c
        if (lmore) then
          lmore = .false.
          read (iunit,6000,err=8000,end=9000) line
          nline = nline + 1
          goto 40
        end if
c
      else if (line(1:3) .eq. 'END') then
c
        if (nmatyp .eq. jj) then
          call errcon ('Residue type has no atoms')
          call textut (' Type :',lrname(ii))
          call textut (' Name :',lrdesc(ii))
          nmrtyp = nmrtyp - 1
          goto 10
        end if
c
        nmrptr (1,ii) = jj + 1
        nmrptr (2,ii) = nmatyp
c
        if (nalias .ne. kk) then
          nmaptr (1,ii) = kk + 1
          nmaptr (2,ii) = nalias
        end if
c
        goto 10
c
      else
c
        call errcon ('Invalid line in residue definition')
        call jvalut (' Line number    :',1,nline)
        call textut (' Offending line :',line)
        ierr = -1
        return
c
      end if
c
      goto 30
c
c ... RING DEFINITION
c
  120 continue
      if (nlring .ge. mxring) then
        call errcon ('Too many ring definitions')
        call jvalut (' Max allowed    :',1,mxring)
        call jvalut (' Line number    :',1,nline)
        call textut (' Offending line :',line)
        ierr = -1
        return
      end if
      nlring = nlring + 1
      rngres (nlring) = line(5:7)
      read (line(9:10),'(i2)',err=8000,end=8000) rngnat (nlring)
      if (rngnat(nlring) .gt. 10) then
        call errcon ('Max 10 ring atoms; rest ignored !')
        rngnat(nlring) = 10
      end if
      read (line(12:),'(10(a4,1x))')
     +  (rngatm(j,nlring),j=1,rngnat(nlring))
c
      goto 10
c
 8000 continue
      call errcon ('While reading file')
      ierr = -1
      return
c
 9000 continue
      write (*,*)
      call jvalut (' Lines read       :',1,nline)
      call jvalut (' Residue types    :',1,nmrtyp)
      call jvalut (' Atom types       :',1,nmatyp)
      call jvalut (' Aliases          :',1,nalias)
      call jvalut (' Ring definitions :',1,nlring)
c
      if (nmrtyp .lt. 10) then
        call errcon ('Fewer than 10 defined residue types')
        ierr = -1
        return
      end if
c
      call prompt ('0First and last residue types:')
      call tellib (lrname(1),i,.true.)
      call tellib (lrname(nmrtyp),i,.true.)
c
c ... check integrity
c
      call prompt ('0Check integrity:')
      lmore = .false.
      do i=1,nmrtyp
        do j=1,nmrtyp
          if (i.ne.j .and. lrname(i) .eq. lrname(j)) then
            lmore = .true.
            write (*,6200) i,lrname(i),j,lrname(j)
          end if
          if (nmaptr(1,j) .gt. 0 .and. i.ne.j) then
            do k=nmaptr(1,j),nmaptr(2,j)
              if (lalias(k) .eq. lrname(i)) then
                lmore = .true.
                write (*,6200) i,lrname(i),j,lalias(k)
              end if
            end do
          end if
        end do
      end do
      if (lmore) then
        call errcon ('Non-unique residue names/aliases')
      else
        call prompt ('There are no name or alias conflicts')
      end if
c
 6200 format (' WARNING - name or alias conflict: ',i5,' = ',
     +  a3,' and ',i5,' = ',a3)
c
      call prompt ('0Count types:')
      do i=1,nrlcat+1
        icnt(i) = 0
      end do
      do i=1,nmrtyp
        icnt(lrtype(i)) = icnt(lrtype(i)) + 1
      end do
      call jvalut (' Nr of amino acid residue types :',
     +  1,icnt(iprot))
      call jvalut (' Nr of nucleic acid types       :',
     +  1,icnt(inucl))
      call jvalut (' Nr of water types              :',
     +  1,icnt(iwate))
      call jvalut (' Nr of metal types              :',
     +  1,icnt(imeta))
      call jvalut (' Nr of inorganic types          :',
     +  1,icnt(iinor))
      call jvalut (' Nr of carbohydrate types       :',
     +  1,icnt(icarb))
      call jvalut (' Nr of organic compound types   :',
     +  1,icnt(iorga))
      call jvalut (' Nr of other compound types     :',
     +  1,icnt(ihete))
c
      if (icnt(iprot) .lt. 20) call errstp (
     +  'READDB - Not enough amino acid types')
      if (icnt(iwate) .lt. 1) call errstp (
     +  'READDB - No water type defined')
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine tellib (myname,ires,lprint)
c
      include 'moleman2.incl'
c
      integer i,j,k,ires,leng1
c
      logical lprint
c
      character myname*(*),name*3
c
code ...
c
      name = myname
      call upcase (name)
c
c ... first: check normal names
c
      do i=1,nmrtyp
        if (name .eq. lrname(i)) then
          ires = i
          if (lprint) then
            write (*,6100) i,lrname(i),libolc(i),
     +       libcat(lrtype(i)),
     +       lrdesc(i)(1:leng1(lrdesc(i)))
            if (nmaptr(1,i) .gt. 0) then
              write (*,6110) (lalias(j),j=nmaptr(1,i),nmaptr(2,i))
            end if
            write (*,6120)
     +        (lratom(j),ismain(j),j=nmrptr(1,i),nmrptr(2,i))
          end if
          return
        end if
      end do
c
c ... second: check aliases
c
      do i=1,nmrtyp
        if (nmaptr(1,i) .gt. 0) then
          do k=nmaptr(1,i),nmaptr(2,i)
            if (name .eq. lalias(k)) then
              ires = i
              if (lprint) then
                write (*,6100) i,lrname(i),libolc(i),
     +           libcat(lrtype(i)),
     +           lrdesc(i)(1:leng1(lrdesc(i)))
                if (nmaptr(1,i) .gt. 0) then
                  write (*,6110) (lalias(j),j=nmaptr(1,i),nmaptr(2,i))
                end if
                write (*,6120)
     +            (lratom(j),ismain(j),j=nmrptr(1,i),nmrptr(2,i))
              end if
              return
            end if
          end do
        end if
      end do
c
      ires = -1
      if (lprint) call textut (' Residue type not found :',name)
c
 6100 format (/' Residue # ',i5,' = ',a3,1x,a1,' (',a4,') = ',a)
 6110 format (10(' A.k.a. ',10(1x,a3),:,/))
 6120 format (10(' Atoms  ',6('  |',a4,'| (',l1,')',:),:,/))
c
      return
      end
c
c
c
      subroutine print_atom (kk)
c
      include 'moleman2.incl'
c
      integer kk,j,leng1
c
code ...
c
 6000 format (' ATOM ',i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a)
      write (*,6000) atomnr(kk),atmnam(kk),altloc(kk),
     +    resnam(kk),achain(kk),iresid(kk),insert(kk),
     +    (xyz(j,kk),j=1,3),qatom(kk),batom(kk),
     +    inote(kk)(1:leng1(inote(kk)))
c
      return
      end
c
c
c
      subroutine print_res (i,n)
c
c ... print N blank lines and then the name of the residue with pointer I
c
      include 'moleman2.incl'
c
      integer i,kk,n,j,length
c
      character mystr*80,ustr*(*)
c
code ...
c
      mystr = 'RESIDUE'
      goto 100
c
      entry string_res (i,n,ustr)
c
      mystr = ustr
c
  100 continue
c
      kk = atmptr(1,i)
      if (n .gt. 0) then
        do j=1,min(100,n)
          write (*,*)
        end do
      end if
c
 6000 format (1x,a,1x,a1,a3,1x,a1,i4,a1,1x,a4)
      write (*,6000) mystr(1:length(mystr)),altloc(kk),resnam(kk),
     +  achain(kk),iresid(kk),insert(kk),inote(kk)(7:10)
c
      return
      end
c
c
c
      subroutine do_select (how,what,which)
c
      include 'moleman2.incl'
c
      integer i,j
c
      logical lhydro
c
      character*(*) how,what,which
      character line*132
c
code ...
c
ccc      call textut (' How   :',how)
ccc      call textut (' What  :',what)
ccc      call textut (' Which :',which)
c
      if (how(1:2) .eq. 'AL') then
        call prompt (' Select ALL atoms')
        do i=1,natoms
          select(i) = .true.
        end do
        nselect = natoms
        selstr = 'ALL |'
c
      else if (how(1:1) .eq. '?') then
        goto 1000
c
      else if (how(1:2) .eq. 'NO') then
        call prompt (' Select NO atoms')
        do i=1,natoms
          select(i) = .false.
        end do
        nselect = 0
        selstr = 'NONE |'
c
      else if (how(1:2) .eq. 'HY') then
        call prompt (' Select HYDROGEN atoms')
        nselect = 0
        do i=1,natoms
          select(i) = lhydro(atmnam(i))
          if (select(i)) nselect = nselect + 1
        end do
        selstr = 'HYDROGEN |'
c
      else if (how(1:2) .eq. 'EX') then
        call prompt (' Select NON-HYDROGEN atoms')
        nselect = 0
        do i=1,natoms
          select(i) = (.not. lhydro(atmnam(i)))
          if (select(i)) nselect = nselect + 1
        end do
        selstr = 'NON-HYDROGEN |'
c
      else if (how(1:2) .eq. 'NE') then
        call prompt (' NEGATE atom selection')
        nselect = 0
        do i=1,natoms
          select(i) = (.not. select(i))
          if (select(i)) nselect = nselect + 1
        end do
        call appstr (selstr,' | NEGATE')
c
      else if (how(1:2) .eq. 'AN' .or.
     +         how(1:2) .eq. 'OR' .or.
     +         how(1:2) .eq. 'BU') then
c
        if (how(1:2) .eq. 'AN') then
          call prompt (' AND atom selection')
        else if (how(1:2) .eq. 'BU') then
          call prompt (' BUtnot atom selection')
        else if (how(1:2) .eq. 'OR') then
          call prompt (' OR atom selection')
        end if
        call textut (' With atoms for which :',what)
        call textut (' Equals :',which)
c
        if (what(1:2) .eq. 'CH') then
          what = 'CHain'
          do i=1,natoms
            lbuf (i) = (achain(i) .eq. which(1:1))
          end do
        else if (what(1:2) .eq. 'SE') then
          what = 'SEgid'
          do i=1,natoms
            lbuf (i) = (inote(i)(7:10) .eq. which(1:4))
          end do
        else if (what(1:2) .eq. 'RE') then
          what = 'REsidue'
          do i=1,natoms
            lbuf (i) = (resnam(i) .eq. which(1:3))
          end do
        else if (what(1:2) .eq. 'AT') then
          what = 'ATom'
          do i=1,natoms
            lbuf (i) = (atmnam(i) .eq. which(1:4))
          end do
        else if (what(1:2) .eq. 'AL') then
          what = 'ALtloc'
          do i=1,natoms
            lbuf (i) = (altloc(i) .eq. which(1:1))
          end do
        else if (what(1:2) .eq. 'AN') then
          what = 'ANisou'
          if (which(1:1) .eq. 'F') then
            do i=1,natoms
              lbuf (i) = (.not. laniso(i))
            end do
          else if (which(1:1) .eq. 'T') then
            do i=1,natoms
              lbuf (i) = (laniso(i))
            end do
          else
            call errcon ('ANiso must be True or False')
          end if
        else if (what(1:2) .eq. 'TY') then
          what = 'TYpe'
          j = -1
          do i=1,nrlcat+1
            if (libcat(i) .eq. which(1:4)) j=i
          end do
          if (j .le. 0) then
            call errcon ('Invalid type')
            call asciut (' Must be one of :',nrlcat+1,libcat)
            goto 1000
          end if
          which = libcat (j)
          do i=1,natoms
            lbuf (i) = (restyp(resptr(i)) .eq. j)
          end do
        else if (what(1:2) .eq. 'ST') then
          what = 'STruc_sec'
          j = -1
          if (which(1:2) .eq. 'LO') j=0
          if (which(1:2) .eq. 'TU') j=0
          if (which(1:2) .eq. 'AL') j=1
          if (which(1:2) .eq. 'BE') j=2
          if (which(1:2) .eq. 'LE') j=3
          if (which(1:2) .eq. 'NO') j=-1
          if (j .le. 0) then
            call errcon ('Invalid secondary structure type')
            call prompt (' Must be one of : LOop, TUrn, ALpha,')
            call prompt ('     BEta, LEft-handed, NOn-protein')
            goto 1000
          end if
          which = ssenam(j)
          do i=1,natoms
            lbuf (i) = (sstype(resptr(i)) .eq. j)
          end do
        else if (what(1:2) .eq. 'CL') then
          what = 'CLass'
          if (which(1:1) .eq. 'M') then
            which = 'Main'
            do i=1,natoms
              lbuf (i) = lmain (i)
            end do
          else if (which(1:1) .eq. 'S') then
            which = 'Side'
            do i=1,natoms
              lbuf (i) = (.not. lmain (i))
            end do
          else
            call errcon ('Invalid class')
            call prompt (' Must be Main or Side chain')
            goto 1000
          end if
        else
          call errcon ('Invalid AND/OR/BUTNOT option')
          goto 1000
        end if
c
        if (how(1:2) .eq. 'AN') then
          nselect = 0
          do i=1,natoms
            select(i) = (select(i) .and. lbuf(i))
            if (select(i)) nselect = nselect + 1
          end do
          line = ' AND '//what(1:6)//' = '//which(1:6)//' | '
          call pretty (line(2:))
          call appstr (selstr,line)
        else if (how(1:2) .eq. 'OR') then
          nselect = 0
          do i=1,natoms
            select(i) = (select(i) .or. lbuf(i))
            if (select(i)) nselect = nselect + 1
          end do
          line = ' OR '//what(1:6)//' = '//which(1:6)//' | '
          call pretty (line(2:))
          call appstr (selstr,line)
        else if (how(1:2) .eq. 'BU') then
          nselect = 0
          do i=1,natoms
            select(i) = (select(i) .and. (.not.lbuf(i)))
            if (select(i)) nselect = nselect + 1
          end do
          line = ' BUTNOT '//what(1:6)//' = '//which(1:6)//' | '
          call pretty (line(2:))
          call appstr (selstr,line)
        end if
c
      else
        call errcon ('Invalid SElection command')
      end if
c
 1000 continue
ccc      print *,'okay'
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine do_numeric_select (andor,what,lo,hi)
c
      include 'moleman2.incl'
c
      real rlo,rhi
c
      integer i,ilo,ihi,ierr,length,leng1,j
c
      character*(*) andor,what,lo,hi
      character line*132
c
code ...
c
      call remspa (andor)
      call remspa (what)
      call remspa (lo)
      call remspa (hi)
c
      if (andor(1:1) .eq. 'A') then
        line = ' AND '
      else if (andor(1:1) .eq. 'B') then
        line = ' BUTNOT '
      else
        line = ' OR '
      end if
c
      do i=1,natoms
        lbuf (i) = .false.
      end do
c
      if (what(1:1) .eq. 'R') then
        call str2i (lo,ilo,ierr)
        if (ierr .ne. 0) return
        call str2i (hi,ihi,ierr)
        if (ierr .ne. 0) return
        call ilohi (ilo,ihi)
        call appstr (line,' Residue_nr ')
        write (line(length(line)+1:),*) ilo,ihi
        do i=1,natoms
          lbuf (i) = (iresid(i) .ge. ilo .and.
     +                iresid(i) .le. ihi)
        end do
c
      else if (what(1:1) .eq. 'S') then
        call str2i (lo,ilo,ierr)
        if (ierr .ne. 0) return
        call str2i (hi,ihi,ierr)
        if (ierr .ne. 0) return
        call ilohi (ilo,ihi)
        call appstr (line,' Sec_struc ')
        write (line(length(line)+1:),*) ilo,ihi
        do i=1,natoms
          j = resptr (i)
          lbuf (i) = (sstype(j) .ge. ilo .and.
     +                sstype(j) .le. ihi)
        end do
c
      else if (what(1:1) .eq. 'E') then
        call str2i (lo,ilo,ierr)
        if (ierr .ne. 0) return
        call str2i (hi,ihi,ierr)
        if (ierr .ne. 0) return
        call ilohi (ilo,ihi)
        call appstr (line,' Element ')
        write (line(length(line)+1:),*) ilo,ihi
        do i=1,natoms
          lbuf (i) = (element(i) .ge. ilo .and.
     +                element(i) .le. ihi)
        end do
c
      else if (what(1:1) .eq. 'B') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' B-factor ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (batom(i) .ge. rlo .and.
     +                batom(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'O') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Occupancy ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (qatom(i) .ge. rlo .and.
     +                qatom(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'M') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Mass ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (atmass(i) .ge. rlo .and.
     +                atmass(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'C') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Cov_radius ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (cvbrad(i) .ge. rlo .and.
     +                cvbrad(i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'X') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' X-coordinate ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (xyz(1,i) .ge. rlo .and.
     +                xyz(1,i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'Y') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Y-coordinate ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (xyz(2,i) .ge. rlo .and.
     +                xyz(2,i) .le. rhi)
        end do
c
      else if (what(1:1) .eq. 'Z') then
        call str2r (lo,rlo,ierr)
        if (ierr .ne. 0) return
        call str2r (hi,rhi,ierr)
        if (ierr .ne. 0) return
        call rlohi (rlo,rhi)
        call appstr (line,' Z-coordinate ')
        write (line(length(line)+1:),*) rlo,rhi
        do i=1,natoms
          lbuf (i) = (xyz(3,i) .ge. rlo .and.
     +                xyz(3,i) .le. rhi)
        end do
c
      else
        call errcon ('Invalid numeric property selected')
        return
      end if
c
      call pretty (line(2:))
      call textut (' Select Numeric :',line)
      line = line(1:leng1(line))//' |'
      call appstr (selstr,line)
c
      if (andor(1:1) .eq. 'A') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      else if (andor(1:1) .eq. 'B') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. (.not. lbuf(i)) )
          if (select(i)) nselect = nselect + 1
        end do
      else
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .or. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      end if
c
 1000 continue
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine do_point_select (andor,mx,my,mz,lo,hi)
c
      include 'moleman2.incl'
c
      real rlo,rhi,xxx,rxyz(3)
      real distce
c
      integer i,ierr,leng1
c
      character*(*) andor,mx,my,mz,lo,hi
      character line*132
c
code ...
c
      call remspa (andor)
      call remspa (mx)
      call remspa (my)
      call remspa (mz)
      call remspa (lo)
      call remspa (hi)
c
      if (andor(1:1) .eq. 'A') then
        line = ' AND '
      else if (andor(1:1) .eq. 'B') then
        line = ' BUTNOT '
      else
        line = ' OR '
      end if
c
      call appstr (line,' POint_distance ')
      call pretty (line(2:))
      line = line(1:leng1(line))//' |'
      call appstr (selstr,line)
c
      do i=1,natoms
        lbuf (i) = .false.
      end do
c
      call str2r (mx,rxyz(1),ierr)
      if (ierr .ne. 0) return
      call str2r (my,rxyz(2),ierr)
      if (ierr .ne. 0) return
      call str2r (mz,rxyz(3),ierr)
      if (ierr .ne. 0) return
      call str2r (lo,rlo,ierr)
      if (ierr .ne. 0) return
      call str2r (hi,rhi,ierr)
      if (ierr .ne. 0) return
      call rlohi (rlo,rhi)
c
      call fvalut (' Point (A) :',3,rxyz)
      call fvalut (' Minimum distance (A) :',1,rlo)
      call fvalut (' Maximum distance (A) :',1,rhi)
c
      do i=1,natoms
        xxx = distce (xyz(1,i),rxyz)
        lbuf (i) = (xxx .ge. rlo .and. xxx .le. rhi)
      end do
c
      if (andor(1:1) .eq. 'A') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      else if (andor(1:1) .eq. 'B') then
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .and. (.not. lbuf(i)) )
          if (select(i)) nselect = nselect + 1
        end do
      else
        nselect = 0
        do i=1,natoms
          select(i) = (select(i) .or. lbuf(i))
          if (select(i)) nselect = nselect + 1
        end do
      end if
c
 1000 continue
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine select_near (lo,hi)
c
      include 'moleman2.incl'
c
      real dist,qq,rlo,rhi
c
      integer i,j,ierr,leng1
c
      character*(*) lo,hi
      character line*132
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
      call remspa (lo)
      call remspa (hi)
c
      do i=1,natoms
        lbuf (i) = .false.
      end do
c
      call str2r (lo,rlo,ierr)
      if (ierr .ne. 0) return
      call str2r (hi,rhi,ierr)
      if (ierr .ne. 0) return
      call rlohi (rlo,rhi)
c
      line = ' DIstance '
      write (line(leng1(line)+1:),'(2f6.2)') rlo,rhi
      line = line(1:leng1(line))//' |'
      call appstr (selstr,line)
c
      do i=1,natoms
        if (select(i)) then
          do j=1,natoms
            if (i .eq. j) goto 1000
            if (select(j)) goto 1000
            if (abs(xyz(1,i)-xyz(1,j)) .gt. rhi) goto 1000
            if (abs(xyz(2,i)-xyz(2,j)) .gt. rhi) goto 1000
            if (abs(xyz(3,i)-xyz(3,j)) .gt. rhi) goto 1000
            qq = dist (i,j,xyz)
            if (qq .ge. rlo .and. qq .le. rhi) then
              lbuf (j) = .true.
            end if
 1000       continue
          end do
        end if
      end do
c
      nselect = 0
      do i=1,natoms
        select(i) = (select(i) .or. lbuf(i))
        if (select(i)) nselect = nselect + 1
      end do
c
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine select_by_residue ()
c
      include 'moleman2.incl'
c
      integer i,j,k,nok
c
code ...
c
      nok = 0
      call prompt (' Select by residue')
      do i=1,nres
        do j=atmptr(1,i),atmptr(2,i)
          if (select(j)) then
            do k=atmptr(1,i),atmptr(2,i)
              select (k) = .true.
              nok = nok + 1
            end do
            goto 10
          end if
        end do
   10   continue
      end do
c
      call appstr (selstr,' BY_residue | ')
c
      nselect = nok
      call textut (' Selection history :',selstr)
      call jvalut (' Nr of selected atoms :',1,nselect)
c
      return
      end
c
c
c
      subroutine list_select (which)
c
      include 'moleman2.incl'
c
      integer i,j,nok
c
      character which*(*)
c
code ...
c
      call upcase (which)
      call remspa (which)
c
      if (which(1:1) .ne. 'A') then
        call prompt (' List first selected atom of every residue')
        nok = 0
        do i=1,nres
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              call print_atom (j)
              nok = nok + 1
              goto 10
            end if
          end do
   10     continue
        end do
        call jvalut (' Nr of residues listed :',1,nok)
c
      else
        call prompt (' List all selected atoms')
        nok = 0
        do i=1,natoms
          if (select(i)) then
            call print_atom (i)
            nok = nok + 1
          end if
        end do
        call jvalut (' Nr of atoms listed :',1,nok)
      end if
c
      return
      end
c
c
c
      subroutine pdb_chem_charge (lprint,lcol)
c
      include 'moleman2.incl'
c
      real mass,radius
c
      integer i,j,ne,nz,np,nn
c
      logical lhydro,lprint,lcol
c
      character try*2,chem*2,charge*2,line*80
c
code ...
c
      if (lprint) call prompt (
     +  ' Deriving chemical name and charge ...')
c
      nz = 0
      nn = 0
      np = 0
      ne = 0
c
      do i=1,natoms
c
c ... is it hydrogen ?
c
        if (lhydro(atmnam(i))) then
          chem = ' H'
          call elinfo (chem,line,j,mass,radius,.false.)
          goto 50
        end if
c
c ... try first two characters of atom name
c
        try = atmnam(i)(1:2)
        call elinfo (try,line,j,mass,radius,.false.)
        if (j .ge. 1) then
          chem = try
          goto 50
        end if
c
c ... try space + second char of atom name
c
        try = ' '//atmnam(i)(2:2)
        call elinfo (try,line,j,mass,radius,.false.)
        if (j .ge. 1) then
          chem = try
          goto 50
        end if
c
c ... try space + first char of atom name
c
        try = ' '//atmnam(i)(1:1)
        call elinfo (try,line,j,mass,radius,.false.)
        if (j .ge. 1) then
          chem = try
          goto 50
        end if
c
c ... unknown element
c
        ne = ne + 1
        call errcon ('Unknown chemical element; assume Carbon')
        call print_atom (i)
        chem = ' C'
        call elinfo (chem,line,j,mass,radius,.false.)
c
c ... try to read charge from atom name
c
   50   continue
c
        element (i) = j
        atmass  (i) = mass
        cvbrad  (i) = 1.0
        if (radius .gt. 0.0) cvbrad (i) = radius
c
        if (lcol) then
          if (element(i) .eq. 1) then
            call xvrml_encode_rgb (1.0,0.0,1.0,colour(i))
          else if (element(i) .eq. 6) then
            call xvrml_encode_rgb (1.0,1.0,0.0,colour(i))
          else if (element(i) .eq. 7) then
            call xvrml_encode_rgb (0.0,0.0,1.0,colour(i))
          else if (element(i) .eq. 8) then
            call xvrml_encode_rgb (1.0,0.0,0.0,colour(i))
          else if (element(i) .eq. 16) then
            call xvrml_encode_rgb (0.0,1.0,0.0,colour(i))
          else if (element(i) .eq. 15) then
            call xvrml_encode_rgb (1.0,0.5,0.0,colour(i))
          else
            call xvrml_col_index (20+3*element(i),colour(i))
          end if
        end if
c
        if (atmnam(i)(3:4) .eq. '++') then
          charge = '+2'
        else if (atmnam(i)(3:4) .eq. '--') then
          charge = '-2'
        else if (atmnam(i)(3:3) .eq. '+') then
          charge = atmnam(i)(3:4)
        else if (atmnam(i)(3:3) .eq. '-') then
          charge = atmnam(i)(3:4)
        else if (atmnam(i)(4:4) .eq. '+') then
          charge = '+'//atmnam(i)(3:3)
        else if (atmnam(i)(4:4) .eq. '-') then
          charge = '-'//atmnam(i)(3:3)
        else
          charge = ' 0'
        end if
c
        if (charge(1:1) .eq. '+' .or. charge(1:1) .eq. '-') then
          if (charge(2:2) .eq. ' ') then
            charge (2:2) = '1'
          else if (charge(2:2) .lt. '1' .or.
     +             charge(2:2) .gt. '9') then
            charge = ' 0'
          end if
        end if
c
        if (lprint) then
          inote (i)(11:12) = chem
          inote (i)(13:14) = charge
        end if
c
        if (charge .eq. ' 0') then
          nz = nz + 1
        else if (charge(1:1) .eq. '+') then
          np = np + 1
        else if (charge(1:1) .eq. '-') then
          nn = nn + 1
        end if
c
      end do
c
      if (lprint) then
        call jvalut (' Nr of atoms processed    :',1,natoms)
        call jvalut (' Unknown chemical element :',1,ne)
        call jvalut (' Nr of positive atoms     :',1,np)
        call jvalut (' Nr of negative atoms     :',1,nn)
        call jvalut (' Nr uncharged or unknown  :',1,nz)
      end if
c
      return
      end
c
c
c
      subroutine pdb_number (first)
c
      include 'moleman2.incl'
c
      integer i,j,nok,ierr,i1,inow
c
      logical select_res
c
      character*(*) first
c
code ...
c
      call str2i (first,i1,ierr)
      if (ierr .ne. 0) return
c
      call jvalut (' Renumber selected residues starting at :',1,i1)
      nok = 0
      inow = i1
      do i=1,nres
        if (select_res(i)) then
          do j=atmptr(1,i),atmptr(2,i)
            iresid(j) = inow
          end do
          nok = nok + 1
          inow = inow + 1
        end if
      end do
      call jvalut (' Nr of last changed residue :',1,(inow-1))
      call jvalut (' Nr of residues changed :',1,nok)
c
      return
      end
c
c
c
      subroutine pdb_name (which,old,new)
c
      include 'moleman2.incl'
c
      integer i,j,nok,length
c
      logical select_res
c
      character*(*) which,old,new
c
code ...
c
      call upcase (which)
      call upcase (old)
      call upcase (new)
      call remspa (which)
c
      if (length(which) .lt. 0) then
        call errcon ('Invalid selection')
        return
      else if (length(old) .lt. 0) then
        call errcon ('Invalid value for old name')
        return
      else if (length(new) .lt. 0) then
        call errcon ('Invalid value for new name')
        return
      end if
c
      if (which(1:1) .ne. 'R') which(1:1) = 'A'
c
 6000 format (' Replace atom name |',a4,'| by |',a4,'|')
 6010 format (' Replace residue name |',a3,'| by |',a3,'|')
c
      if (which(1:1) .eq. 'A') then
        write (*,6000) old,new
        nok = 0
        do i=1,natoms
          if (select(i)) then
            if (atmnam(i) .eq. old(1:4)) then
              nok = nok + 1
              atmnam (i) = new(1:4)
            end if
          end if
        end do
        call jvalut (' Nr of atom names changed :',1,nok)
      else
        write (*,6010) old,new
        nok = 0
        do i=1,nres
          if (resnam(atmptr(1,i)) .eq. old(1:3)) then
            if (select_res(i)) then
              do j=atmptr(1,i),atmptr(2,i)
                resnam(j) = new(1:3)
              end do
              nok = nok + 1
            end if
          end if
        end do
        call jvalut (' Nr of residues changed :',1,nok)
      end if
c
      return
      end
c
c
c
      subroutine pdb_remark (what,which)
c
      include 'moleman2.incl'
c
      integer i,ii,ierr,leng1
c
      character what*(*),which*(*)
c
code ...
c
      call upcase (what)
      call remspa (what)
      if (what(1:1) .ne. 'D' .and. what(1:1) .ne. 'A') what = 'L'
c
 6000 format (1x,i5,': ',a)
c
      if (what(1:1) .eq. 'L') then
        call prompt (' List REMARK records')
        call jvalut (' Nr of REMARK records :',1,nrem)
        if (nrem .gt. 0) then
          do i=1,nrem
            write (*,6000) i,remark(i)(1:leng1(remark(i)))
          end do
        end if
c
      else if (what(1:1) .eq. 'D') then
        if (which(1:1) .eq. '*') then
          call prompt (' Delete all REMARK records')
          nrem = 0
        else
          call str2i (which,ii,ierr)
          if (ierr .ne. 0) return
          if (ii .le. 0 .or. ii .gt. nrem) then
            call errcon ('Invalid REMARK number')
            return
          end if
          call jvalut (' Delete REMARK record nr:',1,ii)
          if (ii .eq. nrem) then
            nrem = nrem - 1
          else
            do i=ii,nrem-1
              remark(i) = remark(i+1)
            end do
            nrem = nrem - 1
          end if
        end if
c
      else if (what(1:1) .eq. 'A') then
        call textut (' Add REMARK record :',which)
        if (nrem .lt. maxcom) then
          nrem = nrem + 1
          remark (nrem) = 'REMARK '//which(1:leng1(which))
          write (*,6000) nrem,remark(nrem)(1:leng1(remark(nrem)))
        else
          call errcon ('Too many REMARK records')
          call jvalut (' Maximum :',1,maxcom)
        end if
      end if
c
      return
      end
c
c
c
      subroutine pdb_farout ()
c
      include 'moleman2.incl'
c
      integer i,i1,i2,nah,nbs,ityp,ifirst,nn
c
      logical busy
c
      character line*128
c
code ...
c
      call prompt (
     +  ' Generating quick-n-dirty HELIX and SHEET records ...')
c
      call pdb_remark ('A',
     +  'YASSPA quick-n-dirty HELIX and SHEET records')
c
      nah = 0
      nbs = 0
c
      busy = .false.
      ityp = 0
c
      do i=1,nres
c
c ... start of a new SSE ?
c
        if (.not. busy .and. sstype(i) .gt. 0) then
          ityp = sstype (i)
          busy = .true.
          ifirst = i
          goto 50
        end if
c
c ... end of an SSE ?
c
        if (busy .and.
     +      (sstype(i+1) .ne. ityp .or. i .eq. nres)) then
          nn = i - ifirst + 1
          i1 = atmptr (1,ifirst)
          i2 = atmptr (1,i)
          if (nn .ge. 3) then
            if (ityp .eq. 1) then
              nah = nah + 1
              write (line,6000) 'HELIX ',nah,nah,
     +          resnam(i1),achain(i1),iresid(i1),insert(i1),
     +          resnam(i2),achain(i2),iresid(i2),insert(i2),
     +          1
            else
              nbs = nbs + 1
              write (line,6100) 'SHEET ',nbs,nbs,1,
     +          resnam(i1),achain(i1),iresid(i1),insert(i1),
     +          resnam(i2),achain(i2),iresid(i2),insert(i2),
     +          0
            end if
c
            call textut (' Record :',line)
            if (nother .lt. maxcom) then
              nother = nother + 1
              other (nother) = line
            else
              call errcon ('Too many records')
              call jvalut (' Maximum :',1,maxcom)
            end if
c
          end if
          busy = .false.
          ityp = 0
          goto 50
        end if
c
   50   continue
      end do
c
 6000 format (a6,1x,i3,1x,i3,2(1x,a3,1x,a1,1x,i4,a1),i2)
 6100 format (a6,1x,i3,1x,i3,i2,2(1x,a3,1x,a1,i4,a1),i2)
c
cATOM    230  CA  ARG    29      23.498  24.524  11.798  1.00 18.43      1CBS 452
cATOM     67  CA  MET A   9      -2.009  19.060  -6.222  1.00 33.52      1CBR 370
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
cHELIX    1   1 PHE     15  LEU     22  1
cHELIX    1   1 GLY A   41  MET A   43  5                                1FSS 170
cSHEET    6   6 1 GLU    70  GLU    72  0) 
cSHEET    1   A 3 LEU A   7  THR A  10  0                                1FSS 194
cSHEET    2   A 3 GLY A  13  MET A  16 -1  N  VAL A  15   O  VAL A   8   1FSS 195
cSHEET    3   A 3 VAL A  57  ASN A  59  1  N  TRP A  58   O  LYS A  14   1FSS 196
c
      call jvalut (' Nr of HELIX records :',1,nah)
      call jvalut (' Nr of SHEET records :',1,nbs)
c
      return
      end
c
c
c
      subroutine pdb_seqres ()
c
      include 'moleman2.incl'
c
      integer maxchn
      parameter (maxchn = 100)
c
      integer i,i1,i2,i3,i4,j,k,ndone
c
      logical rsdone(maxres)
c
      character seqres(maxres)*3
      character line*128,chdone*(maxchn)
c
code ...
c
      call prompt (
     +  ' Generating quick-n-dirty SEQRES records ...')
c
      call pdb_remark ('A',
     +  'Quick-n-dirty SEQRES records')
c
      if (nres .lt. 3) then
        call errcon ('Fewer than 3 residues')
        return
      end if
c
      chdone = ' '
      ndone = 0
      do i=1,nres
        rsdone (i) = .true.
        if (restyp (i) .eq. iprot .or.
     +      restyp (i) .eq. inucl) then
          rsdone (i) = .false.
          i1 = atmptr (1,i)
          if (ndone .lt. 1) then
            ndone = 1
            chdone = achain (i1)
          else
            if (index(chdone,achain(i1)) .le. 0) then
              if (ndone .lt. maxchn) then
                ndone = ndone + 1
                chdone(ndone:ndone) = achain(i1)
              else
                call errcon ('Too many chain IDs - rest skipped')
              end if
            end if
          end if
        end if
      end do
c
      if (ndone .lt. 1) then
        call errcon ('No chain IDs found ???')
        return
      end if
c
      call textut (' Chain IDs found :',chdone)
c
      do j=1,ndone
c
        i2 = 0
        do i=1,nres
          if (.not. rsdone(i)) then
            i1 = atmptr (1,i)
            if (achain(i1) .eq. chdone(j:j)) then
              i2 = i2 + 1
              seqres (i2) = resres (i)
            end if
          end if
        end do
c
        write (*,*)
        call textut (' Chain ID :',chdone(j:j))
        call jvalut (' Nr of residues :',1,i2)
c
        if (i2 .gt. 0) then
c
          i4 = 0
          do i=1,i2,13
            i3 = min(i2,i+12)
            i4 = i4 + 1
            write (line,6000) i4,chdone(j:j),i2,(seqres(k),k=i,i3)
c
ccc            write (*,'(1x,a)') line(1:length(line))
c
            if (nother .lt. maxcom) then
              nother = nother + 1
              other (nother) = line
            else
              call errcon ('Too many records')
              call jvalut (' Maximum :',1,maxcom)
              return
            end if
          end do
          call jvalut (' Nr of SEQRES records :',1,i4)
c
        end if
c
      end do
c
 6000 format ('SEQRES',i4,1x,a1,i5,1x,13(1x,a3))
c
      return
      end
c
c
c
      subroutine pdb_ssbond (mywhat)
c
      include 'moleman2.incl'
c
      real dist,dd
c
      integer i,nok,j,ii,jj,ncys,leng1
c
      character mywhat*(*),what*10
c
code ...
c
      what = mywhat
      call upcase (what)
      call remspa (what)
      if (what(1:1) .ne. 'D' .and. what(1:1) .ne. 'G') what = 'L'
c
      if (what(1:1) .eq. 'L') then
        call prompt (' List SSBOND records')
        nok = 0
        if (nother .gt. 0) then
          do i=1,nother
            if (other(i)(1:6) .eq. 'SSBOND') then
              write (*,'(1x,a)') other(i)(1:leng1(other(i)))
              nok = nok + 1
            end if
          end do
        end if
        call jvalut (' Nr of SSBOND records listed :',1,nok)
c
      else if (what(1:1) .eq. 'D') then
        call prompt (' Delete SSBOND records')
        nok = 0
        if (nother .gt. 0) then
          j = 0
          do i=1,nother
            if (other(i)(1:6) .eq. 'SSBOND') then
              nok = nok + 1
            else
              j = j + 1
              if (i .ne. j) other(j) = other(i)
            end if
          end do
          nother = j
        end if
        call jvalut (' Nr of SSBOND records deleted :',1,nok)
c
      else if (what(1:1) .eq. 'G') then
        call prompt (' Generate SSBOND records')
        nok = 0
        call find_type (natoms,resnam,atmnam,'CYS',' SG ',
     +    ncys,ibuf,.true.)
        if (ncys .lt. 2) then
          call prompt (' No disulfide links')
          return
        end if
c
        call fvalut (' Max SG-SG distance for link :',1,mxcyss)
        nok = 0
        do i=1,ncys-1
          ii = ibuf(i)
          do j=i+1,ncys
            jj = ibuf(j)
            dd = dist (ii,jj,xyz)
            if (dd .le. mxcyss) then
              if (nother .eq. maxcom) then
                call errcon ('No room for more records')
                call jvalut (' Maximum :',1,maxcom)
                return
              end if
              nok = nok + 1
              nother = nother + 1
              write (other(nother),
     + '(a6,1x,i3,1x,a3,a2,1x,i4,a1,3x,a3,a2,1x,i4,a1,4x,a,f6.2,a)')
     +            'SSBOND',nok,'CYS',achain(ii),iresid(ii),
     +            insert(ii),'CYS',achain(jj),iresid(jj),
     +            insert(jj),'S-S = ',dd,' A'
              write (*,'(1x,a)') other(nother)(1:leng1(other(nother)))
            end if
          end do
        end do
c
        call jvalut (' Nr of SSBOND records generated :',1,nok)
c
      end if
c
      return
      end
c
c
c
      subroutine pdb_hetero (what)
c
      include 'moleman2.incl'
c
      integer i,j,n1,n2
c
      character what*(*)
c
code ...
c
      call upcase (what)
      call remspa (what)
      if (what(1:1) .ne. 'A') what = 'D'
c
      if (what(1:1) .eq. 'A') then
        do i=1,natoms
          lhet (i) = .false.
        end do
        call prompt (' All atoms set to type ATOM')
      else if (what(1:1) .eq. 'D') then
        n1 = 0
        n2 = 0
        call prompt (' Deducing ATOM/HETATM types ...')
        do i=1,natoms
          j = restyp(resptr(i))
          if (j .eq. iprot .or. j .eq. inucl) then
            lhet (i) = .false.
            n1 = n1 + 1
          else
            lhet (i) = .true.
            n2 = n2 + 1
          end if
        end do
        call jvalut (' Nr set to ATOM   :',1,n1)
        call jvalut (' Nr set to HETATM :',1,n2)
      end if
c
      return
      end
c
c
c
      subroutine mc_check (iunit,file,what,which)
c
      include 'moleman2.incl'
c
      real dist,tangle,omega,plan,phi,psi,dx,dy,ptsize,chir
      real lefthh(4)
c
      integer coregn(37,37)
      integer i,nok,j,iphi,ipsi,nend,nerr,ngly,nout,nyes,ndaa
      integer ncis,nraar,nbent,iunit,mode,length,ierr,kk,leng1
      integer nrleft,i1left,inleft
c
      logical lstart,lend,lomega,lcis,lraar,lplan,lbent,lpsi,lphi
      logical lpos,lrama,lp,label,ldaa,lchir
c
      character*(*) file,what,which
      character line*256,mychn*1
c
      data lefthh /30.0,130.0,-50.0,100.0/
c
code ...
c
      call defcor (coregn)
c
      mode = 0
      if (length(file) .gt. 0) then
        call upcase (what)
        call remspa (what)
        if (what(1:1) .eq. 'T') then
          mode = 1
        else if (what(1:1) .eq. 'L') then
          mode = 2
          label = .true.
        else
          mode = 2
          label = .false.
        end if
      end if
c
      mychn = which(1:1)
      call upcase (mychn)
c
      if (mode .eq. 1) then
c
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
c
      else if (mode .eq. 2) then
c
        dx = 2.0
        dy = 2.0
        ptsize = 10.0
c
        call psrini (iunit,file,prognm)
c
      end if
c
      do i=1,5*nres
        ibuf (i) = -1
      end do
c
      do i=1,2*nres
        rbuf (i) = -999.99
      end do
c
      write (*,*)
      call textut (' Chain ID to check (* = all) :',mychn)
c
      nok = 0
      do i=1,nres
        lbuf (i) = (restyp(i).eq.iprot)
        if (mychn .ne. '*') then
          lbuf(i) = ( lbuf(i) .and.
     +                (achain(atmptr(1,i)) .eq. mychn) )
        end if
        if (lbuf(i)) then
          nok = nok + 1
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' N  ') then
              ibuf (i) = j
            else if (atmnam(j) .eq. ' CA ') then
              ibuf (nres+i) = j
            else if (atmnam(j) .eq. ' C  ') then
              ibuf (2*nres+i) = j
            else if (atmnam(j) .eq. ' O  ') then
              ibuf (3*nres+i) = j
            else if (atmnam(j) .eq. ' CB ') then
              ibuf (4*nres+i) = j
            end if
          end do
        end if
      end do
c
      call jvalut (' Nr of residues to check :',1,nok)
      if (nok .lt. 3) then
        call errcon ('I don''t call this a protein ...')
        call xps_delete ()
        call prompt (' PostScript file empty and deleted')
        return
      end if
c
      nend = 0
      nerr = 0
      ngly = 0
      ndaa = 0
      nout = 0
      nyes = 0
      ncis = 0
      nraar = 0
      nbent = 0
c
      nrleft = 0
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,4(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,4(1x,f10.2))
 6010 format (/
     +  ' Listing of Phi, Psi, Omega and Planarity (pseudo) torsion'/
     +  ' angles for PDB file ',a,' chain ',a1/
     + /' PHI = C(i-1) - N(i) - CA(i) - C(i) - usually negative'/
     +  ' PSI = N(i) - CA(i) - C(i) - N(i+1)'/
     +  ' OMEGA = CA(i) - C(i) - N(i+1) - CA (i+1) - cis=0, trans=180'/
     +  ' PLANARITY = C(i) - CA(i) - N(i+1) - O(i) - planar < 5'/
     + /' An entry of "-999.9" means that the torsion angle could'/
     +  ' not be calculated (for terminal residues and residues with'/
     +  ' missing atoms).'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil)),mychn
        write (iunit,6100,err=9000)
     +    'Residue','Phi','Psi','Omega','Planarity'
        write (iunit,6100,err=9000)
     +    '-------','---','---','-----','---------'
      end if
c
      do i=1,nres
c
        if (.not. lbuf(i)) then
          nrleft = 0
          goto 1000
        end if
c
        if (captr(i) .le. 0) then
          call print_res (i,1)
          call prompt (' ERROR - missing CA atom')
          nerr = nerr + 1
          nrleft = 0
          goto 1000
        end if
c
        lstart = .false.
        lend = .false.
        lomega = .false.
        lcis = .false.
        lraar = .false.
        lchir = .false.
        ldaa = .false.
        lbent = .false.
        lphi = .false.
        lpsi = .false.
        lpos = .false.
        lrama = .false.
c
        omega = -999.9
        plan = -999.9
        phi = -999.9
        psi = -999.9
c
        if (i .eq. 1) then
          lstart = .true.
        else if (captr(i-1) .gt. 0) then
          lstart =  (dist(captr(i-1),captr(i),xyz) .gt. mxcaca)
        else
          lstart = .true.
        end if
c
        if (i .eq. nres) then
          lend = .true.
        else  if (captr(i+1) .gt. 0) then
          lend =  (dist(captr(i+1),captr(i),xyz) .gt. mxcaca)
        else
          lend = .true.
        end if
c
        if (.not. lend) then
c
c ... OMEGA = CA(i) - C(i) - N(i+1) - CA (i+1)
c
          lomega = (ibuf(nres+i) .gt. 0 .and.
     +              ibuf(2*nres+i) .gt. 0 .and.
     +              ibuf(i+1) .gt. 0 .and.
     +              ibuf(nres+i+1) .gt. 0)
          if (lomega) then
            omega = tangle (ibuf(nres+i),ibuf(2*nres+i),
     +                ibuf(i+1),ibuf(nres+i+1),xyz)
            lcis = (abs(omega) .le. 30)
            lraar = (abs(omega) .gt. 30 .and.
     +               abs(omega) .lt. 150)
          end if
c
c ... PLANARITY = C(i) - CA(i) - N(i+1) - O(i)
c
          lplan = (ibuf(2*nres+i) .gt. 0 .and.
     +             ibuf(nres+i) .gt. 0 .and.
     +             ibuf(i+1) .gt. 0 .and.
     +             ibuf(3*nres+i) .gt. 0)
          if (lplan) then
            plan = tangle (ibuf(2*nres+i),ibuf(nres+i),
     +               ibuf(i+1),ibuf(3*nres+i),xyz)
            lbent = (abs(plan) .gt. 5.0)
          end if
c
c ... PSI = N(i) - CA(i) - C(i) - N(i+1)
c
          lpsi = (ibuf(i) .gt. 0 .and.
     +            ibuf(nres+i) .gt. 0 .and.
     +            ibuf(2*nres+i) .gt. 0 .and.
     +            ibuf(i+1) .gt. 0)
          if (lpsi) then
            psi = tangle (ibuf(i),ibuf(nres+i),
     +              ibuf(2*nres+i),ibuf(i+1),xyz)
          end if
        end if
c
c ... PHI = C(i-1) - N(i) - CA(i) - C(i)
c
        if (.not. lstart) then
          lphi = (ibuf(2*nres+i-1) .gt. 0 .and.
     +            ibuf(i) .gt. 0 .and.
     +            ibuf(nres+i) .gt. 0 .and.
     +            ibuf(2*nres+i) .gt. 0)
          if (lphi) then
            phi = tangle (ibuf(2*nres+i-1),ibuf(i),
     +              ibuf(nres+i),ibuf(2*nres+i),xyz)
            if (resnam(atmptr(1,i)) .ne. 'GLY') then
              lpos = (phi .ge. 0.0)
            end if
          end if
        end if
c
c ... CHIR = CA - N - C - CB
c
        lchir = (ibuf(nres+i) .gt. 0 .and.
     +           ibuf(i) .gt. 0 .and.
     +           ibuf(2*nres+i) .gt. 0 .and.
     +           ibuf(4*nres+i) .gt. 0)
        if (lchir) then
          chir = tangle (ibuf(nres+i),ibuf(i),
     +             ibuf(2*nres+i),ibuf(4*nres+i),xyz)
          ldaa = (chir .lt. 0.0)
        end if
c
        if (lphi .and. lpsi) then
          rbuf (i) = phi
          rbuf (nres+i) = psi
          if (resnam(atmptr(1,i)) .ne. 'GLY') then
            if (.not. ldaa) then
              iphi = int ( (180.0 + phi) / 10.0 ) + 1
              ipsi = int ( (180.0 + psi) / 10.0 ) + 1
              lrama = (coregn(iphi,ipsi) .eq. 1)
            else
              iphi = int ( (180.0 - phi) / 10.0 ) + 1
              ipsi = int ( (180.0 - psi) / 10.0 ) + 1
              lrama = (coregn(iphi,ipsi) .eq. 1)
              iphi = int ( (180.0 + phi) / 10.0 ) + 1
              ipsi = int ( (180.0 + psi) / 10.0 ) + 1
            end if
          end if
        end if
c
        lp = .false.
c
        if (lstart) then
          call print_res(i,2)
          lp = .true.
          call prompt (' <<<<< Start of new chain')
          if (mode .eq. 1) write (iunit,*,err=9000)
        end if
c
        if (lend) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' >>>>> End of chain')
        end if
c
        if (.not. lend) then
          if (.not. lomega) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate OMEGA')
          else
            if (lcis) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,f8.1,a,a3)')
     +          ' Cis-peptide; omega = ',omega,'; next residue is ',
     +          resnam(atmptr(1,i+1))
              call prompt (line)
            end if
            if (lraar) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - strange omega = ',omega
              call prompt (line)
            end if
          end if
          if (.not. lplan) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate planarity')
          else
            if (lbent) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - non-planar peptide; planarity = ',plan
              call prompt (line)
            end if
          end if
          if (.not. lpsi) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate PSI')
          end if
        end if
c
        if (.not. lstart) then
          if (.not. lphi) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate PHI')
          else
            if (resnam(atmptr(1,i)) .ne. 'GLY') then
              if (lpos) then
                if (.not. lp) call print_res(i,1)
                lp = .true.
                write (line,'(a,f8.1)')
     +            ' Warning - positive PHI = ',phi
                call prompt (line)
              end if
            end if
          end if
        end if
c
        if (lchir) then
          if (ldaa) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            write (line,'(a,f8.1)')
     +        ' Warning - D-amino acid; improper = ',chir
            call prompt (line)
          end if
        end if
c
        if (lphi .and. lpsi) then
          if (resnam(atmptr(1,i)) .ne. 'GLY') then
c
c ... plus for all non-glycines
c
            if (mode .eq. 2) then
              if (.not. ldaa) then
                call xps_colour (0)
                call xps_symbol (1,phi-dx,phi+dx,psi-dy,psi+dy)
              else
                call xps_colour (1)
                call xps_symbol (4,phi-dx,phi+dx,psi-dy,psi+dy)
              end if
            end if
            if (.not. lrama) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              write (line,'(a,2f8.1)')
     +          ' Warning - unusual PHI,PSI combination = ',phi,psi
              call prompt (line)
c
c ... add cross for outliers and optional label with residue name etc.
c
              if (mode .eq. 2) then
                call xps_symbol (2,phi-dx,phi+dx,psi-dy,psi+dy)
                if (label) then
                  kk = atmptr(1,i)
                  write (line,6000) resnam(kk),achain(kk),iresid(kk),
     +              insert(kk),inote(kk)(7:10)
                  call remspa(line)
                  call xps_colour (4)
                  call xps_text (phi+dx+0.5,psi-dy-0.5,ptsize,line)
                end if
              end if
            end if
          else
c
c ... square for glycines
c
            if (mode .eq. 2) then
              call xps_colour (0)
              call xps_symbol (0,phi-dx,phi+dx,psi-dy,psi+dy)
            end if
          end if
        end if
c
        if (lstart .or. lend) then
          nend = nend + 1
        else if (.not. (lphi .and. lpsi)) then
          nerr = nerr + 1
        else if (resnam(atmptr(1,i)) .eq. 'GLY') then
          ngly = ngly + 1
        else if (.not. lrama) then
          nout = nout + 1
        else
          nyes = nyes + 1
        end if
        if (lcis) ncis = ncis + 1
        if (lraar) nraar = nraar + 1
        if (lbent) nbent = nbent + 1
        if (ldaa) ndaa = ndaa + 1
c
        if (mode .eq. 1) then
          kk = atmptr(1,i)
          write (iunit,6200,err=9000)
     +      resnam(kk),achain(kk),iresid(kk),
     +      insert(kk),inote(kk)(7:10),phi,psi,omega,plan
        end if
c
c ... check for left-handed helices
c
        if (lphi .and. lpsi) then
          if (phi .ge. lefthh(1) .and. phi .le. lefthh(2) .and.
     +        psi .ge. lefthh(3) .and. psi .le. lefthh(4)) then
            nrleft = nrleft + 1
            if (nrleft .eq. 1) i1left = i
            inleft = i
            goto 1000
          end if
        end if
c
        if (nrleft .ge. 4) then
          write (*,*)
          call jvalut (
     +      ' Left-handed helix - nr of residues :',1,nrleft)
          call string_res (i1left,0,'Left-handed helix first :')
          call string_res (inleft,0,'Left-handed helix last  :')
        end if
        nrleft = 0
c
 1000   continue
c
      end do
c
      write (*,*)
      call textut (' Chain ID checked (* = all)   :',mychn)
      call jvalut (' Total nr of residues checked :',1,nok)
      call jvalut (' D-amino acids                :',1,ndaa)
      call jvalut (' Cis-peptide bonds            :',1,ncis)
      call jvalut (' Unusual OMEGA values         :',1,nraar)
      call jvalut (' Non-planar peptide bonds     :',1,nbent)
      call jvalut (' Start/end residues           :',1,nend)
      call jvalut (' Problems (missing atoms ?)   :',1,nerr)
      call jvalut (' Glycine residues             :',1,ngly)
      call jvalut (' Remaining residues in Ramach :',1,(nout+nyes))
      call jvalut ('   In core regions            :',1,nyes)
      call jvalut ('   Outliers                   :',1,nout)
      if ( (nout+nyes) .gt. 0) then
        call fvalut ('   Outlier percentage         :',1,
     +    100.0*float(nout)/float(nout+nyes))
      end if
      call prompt (' An average <= 2.0 A model has ~0-5% outliers')
      call prompt (' See: Kleywegt, G.J. and Jones, T.A. (1996).')
      call prompt ('      Structure 4, 1395-1400.')
      if (ndaa .gt. 0) then
        call prompt (' For D-amino acids, -Phi and -Psi were used')
      end if
c
      if (mode .eq. 1) then
c
        write (iunit,*)
        close (iunit)
c
      else if (mode .eq. 2) then
c
        line = ' PDB file : '//pdbfil
        call xps_legend (line)
c
        write (line,*)
     +    ' Glycines (open squares):',ngly,
     +    ' ; Start/end residues :',nend
        call xps_legend (line)
c
        write (line,*)
     +    ' D-amino acids :',ndaa,
     +    ' ; Residues with missing atoms :',nerr
        call xps_legend (line)
c
        write (line,*) ' Residues in Ramachandran plot checked :',
     +    (nout+nyes),' out of ',nok
        call xps_legend (line)
c
        write (line,*) ' In core regions (plus signs):',nyes,
     +    ' ; Outliers (asterisks):',nout
        call xps_legend (line)
c
        if ( (nout+nyes) .gt. 0) then
          write (line,'(a,1x,f8.1)')
     +      ' Percentage outliers:',
     +      100.0*float(nout)/float(nout+nyes)
          call xps_legend (line)
        end if
c
        line = ' An average <= 2.0 A model has ~0-5% outliers'
        call xps_legend (line)
c
        line = ' See: Kleywegt, G.J. and Jones, T.A. (1996). '//
     +         'Structure 4, 1395-1400.'
        call xps_legend (line)
c
        if (ndaa .gt. 0) then
          line = ' For D-amino acids, -Phi and -Psi were used !'
          call xps_legend (line)
        end if
c
        if ((nout+nyes) .gt. 0) then
          call xps_close ()
          call prompt (' PostScript file created')
        else
          call xps_delete ()
          call prompt (' PostScript file empty and deleted')
        end if
c
      end if
c
      return
c
 9000 continue
      call errcon ('While writing text file')
c
      return
      end
c
c
c
      subroutine sc_check (iunit,file)
c
      include 'moleman2.incl'
c
      real chi(5)
      real tangle,chir,cideal,ctoler,xbad
c
      integer iptr(9)
      integer i,nok,j,iunit,mode,length,ierr,k,nchi,kk
      integer ndaa,nbad,nswap,leng1
c
      logical llaa,lbadc,lp,lswap
c
      character*(*) file
      character line*256
c
code ...
c
      cideal = 34.0
      ctoler = tortol
c
      mode = 0
      if (length(file) .gt. 0) then
        mode = 1
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
      end if
c
      call fvalut (' Ideal chirality improper :',1,cideal)
      call fvalut (' Ideal flatness improper  :',1,0.0)
      call fvalut (' Tolerance                :',1,ctoler)
      write (*,*)
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,6(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,6(1x,f10.2))
 6010 format (/
     +  ' Listing of Chi1-5 and Chirality (pseudo) torsion'/
     +  ' angles for PDB file ',a/
     + /' CHIRAL = CA    - N     - C     - CB (+34=L-aa; -34=D-aa)'/
     +  ' CHI-1  = N     - CA    - CB    - ?G(1)'/
     +  ' CHI-2  = CA    - CB    - ?G(1) - ?D(1)'/
     +  ' CHI-3  = CB    - ?G(1) - ?D(1) - ?E(1)'/
     +  ' CHI-4  = ?G(1) - ?D(1) - ?E(1) - ?Z(1)'/
     +  ' CHI-5  = ?D(1) - ?E(1) - ?Z(1) - ?H(1)'/
     + /' An entry of "-999.9" means that the torsion angle could'/
     +  ' not be calculated (e.g., for residues with missing atoms).'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil))
        write (iunit,6100,err=9000)
     +    'Residue','CA-chir','Chi-1','Chi-2','Chi-3','Chi-4','Chi-5'
        write (iunit,6100,err=9000)
     +    '-------','-------','-----','-----','-----','-----','-----'
      end if
c
      nok = 0
      ndaa = 0
      nbad = 0
      nswap = 0
c
      do i=1,nres
        if (restyp(i) .ne. iprot) goto 100
        if (resnam(atmptr(1,i)) .eq. 'GLY') goto 100
c
        do k=1,9
          iptr(k) = -1
        end do
c
        do k=1,5
          chi(k) = -999.9
        end do
c
        llaa = .true.
        lbadc = .false.
        lswap = .false.
c
        do j=atmptr(1,i),atmptr(2,i)
          if (atmnam(j) .eq. ' C  ') then
            iptr(1) = j
          else if (atmnam(j) .eq. ' N  ') then
            iptr(2) = j
          else if (atmnam(j) .eq. ' CA ') then
            iptr(3) = j
          else if (atmnam(j) .eq. ' CB ') then
            iptr(4) = j
          else if (atmnam(j)(3:4) .eq. 'G ' .or.
     +             atmnam(j)(3:4) .eq. 'G1') then
            iptr(5) = j
          else if (atmnam(j)(3:4) .eq. 'D ' .or.
     +             atmnam(j)(3:4) .eq. 'D1') then
            iptr(6) = j
          else if (atmnam(j)(3:4) .eq. 'E ' .or.
     +             atmnam(j)(3:4) .eq. 'E1') then
            iptr(7) = j
          else if (atmnam(j)(3:4) .eq. 'Z ' .or.
     +             atmnam(j)(3:4) .eq. 'Z1') then
            iptr(8) = j
          else if (atmnam(j)(3:4) .eq. 'H ' .or.
     +             atmnam(j)(3:4) .eq. 'H1') then
            iptr(9) = j
          end if
        end do
c
c ... CA chirality
c
        if (iptr(1).le.0 .or. iptr(2).le.0 .or.
     +      iptr(3).le.0 .or. iptr(4).le.0) then
          call print_res (i,1)
          call prompt (' ERROR - missing atom(s)')
          goto 100
        end if
c
        chir = tangle (iptr(3),iptr(2),iptr(1),iptr(4),xyz)
        if (chir .lt. 0) llaa = .false.
        if (llaa) then
          lbadc = ( abs(chir - cideal) .gt. ctoler)
        else
          lbadc = ( abs(-chir - cideal) .gt. ctoler)
        end if
c
        nchi = 0
        if (iptr(5) .le. 0) goto 90
        chi(1) = tangle (iptr(2),iptr(3),iptr(4),iptr(5),xyz)
        nchi = nchi + 1
c
        if (iptr(6) .le. 0) goto 90
        chi(2) = tangle (iptr(3),iptr(4),iptr(5),iptr(6),xyz)
        nchi = nchi + 1
c
        if (resnam(atmptr(1,i)) .eq. 'ASP' .or.
     +      resnam(atmptr(1,i)) .eq. 'PHE' .or.
     +      resnam(atmptr(1,i)) .eq. 'TYR') then
          lswap = (abs(chi(2)) .gt. 90.0)
          xbad = chi(2)
        end if
c
        if (iptr(7) .le. 0) goto 90
        chi(3) = tangle (iptr(4),iptr(5),iptr(6),iptr(7),xyz)
        nchi = nchi + 1
c
        if (resnam(atmptr(1,i)) .eq. 'GLU') then
          lswap = (abs(chi(3)) .gt. 90.0)
          xbad = chi(3)
        end if
c
        if (iptr(8) .le. 0) goto 90
        chi(4) = tangle (iptr(5),iptr(6),iptr(7),iptr(8),xyz)
        nchi = nchi + 1
c
        if (iptr(9) .le. 0) goto 90
        chi(5) = tangle (iptr(6),iptr(7),iptr(8),iptr(9),xyz)
        nchi = nchi + 1
c
        if (resnam(atmptr(1,i)) .eq. 'ARG') then
          lswap = (abs(chi(5)) .gt. 90.0)
          xbad = chi(5)
        end if
c
   90   continue
        nok = nok + 1
        lp = .false.
c
        if (.not. llaa) then
          if (.not. lp) call print_res (i,1)
          lp = .true.
          call prompt (' Warning - D-amino acid !')
          ndaa = ndaa + 1
        end if
c
        if (lbadc) then
          if (.not. lp) call print_res (i,1)
          lp = .true.
          write (line,'(a,f8.1)')
     +      ' Warning - poor CA chirality improper = ',chir
          call prompt (line)
          nbad = nbad + 1
        end if
c
        if (lswap) then
          if (.not. lp) call print_res (i,1)
          lp = .true.
          write (line,'(a,f8.1)')
     +      ' Warning - sidechain 1/2 atom names incorrect; Chi = ',
     +      xbad
          call prompt (line)
          nswap = nswap + 1
        end if
c
c ... residue-specific checks can go here
c
c ... Ile - CB chirality
c
        if (resnam(atmptr(1,i)) .eq. 'ILE') then
          call get_geom (i,' CB ',' CG1',' CG2',' CA ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' ERROR - wrong CB chirality; improper = ',chir
              call prompt (line)
            end if
          end if
c
c ... Leu, Val, Thr - CG / CB tetrahedral
c
        else if (resnam(atmptr(1,i)) .eq. 'LEU') then
          call get_geom (i,' CG ',' CD2',' CD1',' CB ',chir,ierr)
ccc      print *,'LEU ',chir,ierr
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG tetrahedral; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'VAL') then
          call get_geom (i,' CB ',' CG2',' CG1',' CA ',chir,ierr)
ccc      print *,'VAL ',chir,ierr
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CB tetrahedral; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'THR') then
          call get_geom (i,' CB ',' OG1',' CG2',' CA ',chir,ierr)
ccc      print *,'THR ',chir,ierr
          if (ierr .eq. 0) then
            if (abs(chir-cideal) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CB tetrahedral; improper = ',chir
              call prompt (line)
            end if
          end if
        end if
c
c ... sp2-carbon flatness
c
        if (resnam(atmptr(1,i)) .eq. 'ARG') then
          call get_geom (i,' CZ ',' NH1',' NH2',' NE ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CZ flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'ASN') then
          call get_geom (i,' CG ',' OD1',' ND2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'ASP') then
          call get_geom (i,' CG ',' OD1',' OD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'GLN') then
          call get_geom (i,' CD ',' OE1',' NE2',' CG ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CD flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'GLU') then
          call get_geom (i,' CD ',' OE1',' OE2',' CG ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CD flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'HIS') then
          call get_geom (i,' CG ',' ND1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'PHE') then
          call get_geom (i,' CG ',' CD1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'TYR') then
          call get_geom (i,' CG ',' CD1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        else if (resnam(atmptr(1,i)) .eq. 'TRP') then
          call get_geom (i,' CG ',' CD1',' CD2',' CB ',chir,ierr)
          if (ierr .eq. 0) then
            if (abs(chir) .gt. ctoler) then
              if (.not. lp) call print_res (i,1)
              lp = .true.
              write (line,'(a,f8.1)')
     +          ' Warning - poor CG flatness; improper = ',chir
              call prompt (line)
            end if
          end if
        end if
c
c ... ring flatness here ?
c
c
c ... favourable chi1-4 torsion angles here ?
c
c
c ... favourable chi1-2 rotamers here ?
c
c
c ... add to text file
c
        if (mode .eq. 1) then
          kk = atmptr(1,i)
          if (nchi .le. 0) then
            write (iunit,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),chir
          else
            if (resnam(atmptr(1,i)) .eq. 'TYR' .or.
     +          resnam(atmptr(1,i)) .eq. 'PHE' .or.
     +          resnam(atmptr(1,i)) .eq. 'HIS' .or.
     +          resnam(atmptr(1,i)) .eq. 'TRP') then
              nchi = min (2, nchi)
            end if
            write (iunit,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),chir,(chi(j),j=1,nchi)
          end if
        end if
c
  100   continue
      end do
c
 8000 continue
      write (*,*)
      call jvalut (' Nr of residues checked    :',1,nok)
      call jvalut (' Nr of D-amino acids       :',1,ndaa)
      call jvalut (' Nr with poor CA chirality :',1,nbad)
      call jvalut (' Nr wrong 1/2 atom names   :',1,nswap)
c
      if (mode .eq. 1) close (iunit)
c
      return
c
 9000 continue
      call errcon ('While writing file')
      return
c
      end
c
c
c
      subroutine get_geom (ires,a1,a2,a3,a4,val,ierr)
c
      include 'moleman2.incl'
c
      real dist,angle,tangle,val
c
      integer ip(4),ires,ierr,ndo,nok,i,j
c
      character*4 a1,a2,a3,a4,ats(4)
c
code ...
c
      ierr = -1
      if (a1 .eq. ' ' .or. a2 .eq. ' ') return
c
      ats(1) = a1
      ats(2) = a2
      ats(3) = a3
      ats(4) = a4
c
      ndo = 4
      do i=1,4
        if (ats(i) .eq. ' ') then
          ndo = i-1
          goto 10
        end if
      end do
c
   10 continue
      if (ndo .lt. 2) return
c
      nok = 0
      do i=1,ndo
        ip(i) = -1
        do j=atmptr(1,ires),atmptr(2,ires)
          if (atmnam(j) .eq. ats(i)) then
            ip(i) = j
            nok = nok + 1
            goto 20
          end if
        end do
   20   continue
      end do
c
      if (nok .ne. ndo) return
c
      if (ndo .eq. 2) then
        val = dist (ip(1),ip(2),xyz)
      else if (ndo .eq. 3) then
        val = angle (ip(1),ip(2),ip(3),xyz)
      else if (ndo .eq. 4) then
        val = tangle (ip(1),ip(2),ip(3),ip(4),xyz)
      end if
c
      ierr = 0
c
      return
      end
c
c
c
      subroutine ca_check (iunit,file,what,which)
c
      include 'moleman2.incl'
c
      real dist,angle,tangle,dx,dy,ptsize,phi,psi,cadist,cators
      real cangle,ave,sdv,xmin,xmax,xtot,x,y
c
      integer caqual(0:60,-60:60)
      integer lstat(maxres),dcnts(5),atcnts(0:3)
      integer i,nok,j,iunit,mode,length,ierr,kk,i1,natok,i2,ii
      integer ncheck,niso,leng1
c
      logical lp,label,lend,lgly,lbadd,lbadat,lstart
c
      character*(*) file,what,which
      character line*256,labx*40,laby*40,mychn*1
c
code ...
c
      call initca (caqual(0,-60))
c
c ... bug fix
c
      do i=0,60
        caqual (i,-60) = max (caqual (i,-59),caqual (i,-60),
     +                        caqual (i,60))
      end do
c
      do i=1,5
        dcnts(i) = 0
      end do
c
      do i=0,3
        atcnts(i) = 0
      end do
c
      mode = 0
      if (length(file) .gt. 0) then
        call upcase (what)
        call remspa (what)
        if (what(1:1) .eq. 'T') then
          mode = 1
        else if (what(1:1) .eq. 'L') then
          mode = 2
          label = .true.
        else
          mode = 2
          label = .false.
        end if
      end if
c
      mychn = which(1:1)
      call upcase (mychn)
c
      if (mode .eq. 1) then
c
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
c
      else if (mode .eq. 2) then
c
        dx = 1.0
        dy = 2.0
        ptsize = 10.0
c
        call xps_init ()
        call xps_open (iunit,file,'MOLEMAN2')
        call xps_scale (0.0,180.0,-180.0,180.0)
        call xps_stroke ()
        call xps_ps_comment ('Set up CA-Ramachandran plot')
c
        do i=0,60
          do j=-60,60
            phi = float(i)*3.0
            psi = float(j)*3.0
            if (caqual(i,j) .eq. 1) then
              call xps_dark_box (phi,phi+3.0,psi,psi+3.0)
            else if (caqual(i,j) .eq. 2) then
              call xps_grey_box (phi,phi+3.0,psi,psi+3.0)
            else if (caqual(i,j) .eq. 3) then
              call xps_light_box (phi,phi+3.0,psi,psi+3.0)
            end if
          end do
        end do
c
        call xps_move (  0.,-180.)
        call xps_draw (180.,-180.)
        call xps_draw (180., 180.)
        call xps_draw (  0., 180.)
        call xps_draw (  0.,-180.)
c
c ---	phi axis ticks
c
        do 300 i=1,5
          phi = i*30.
          call xps_move (phi,-180.)
          call xps_draw (phi,-175.)
          call xps_move (phi, 180.)
          call xps_draw (phi, 175.)
300     continue
c
c ---	psi axis ticks
c
        do 310 i=1,11
          psi = -180+i*30.
          call xps_move (  0.,psi)
          call xps_draw (  3.,psi)
          call xps_move (180.,psi)
          call xps_draw (177.,psi)
310     continue
c
        labx = 'CA-CA-CA angle mapped to [0,180>'
        laby = 'CA-CA-CA-CA dihedral mapped to [-180,180>'
        call xps_label (labx,laby)
        call xps_ps_comment ('CA-Ramachandran initialisation done')
c
      end if
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,4(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,4(1x,f10.2))
 6010 format (/
     +  ' Listing for PDB file ',a,' chain ',a1/
     + /' CA(i)   - CA(i+1) distance (2.9=cis, 3.8=trans)'/
     +  ' CA(i-1) - CA(i) - CA(i+1) angle'/
     +  ' CA(i-1) - CA(i) - CA(i+1) - CA(i+2) pseudo-torsion angle'/
     + /' An entry of "-999.9" means that the value could not be'/
     +  ' calculated (for near-terminal residues)'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil)),mychn
        write (iunit,6100,err=9000)
     +    'Residue','Distance','Angle','Torsion'
        write (iunit,6100,err=9000)
     +    '-------','--------','-----','-------'
      end if
c
      nok = 0
      natok = 0
      ncheck = 0
      niso = 0
c
      do i=1,nres
c
        lstat (i) = 0
c
        if (mychn .ne. '*') then
          if (achain(atmptr(1,i)) .ne. mychn) then
            lstat (i) = -10
            goto 1000
          end if
        end if
c
        if (captr(i) .le. 0) then
          lstat (i) = -10
          goto 1000
        end if
c
        lstart = .false.
        lend = .false.
        lgly = (resnam(captr(i)) .eq. 'GLY')
        lbadd = .false.
        lbadat = .false.
c
        cadist = -999.9
        cangle = -999.9
        cators = -999.9
c
        if (i .eq. 1) then
          lstart = .true.
        else if (captr(i-1) .gt. 0) then
          lstart = (dist(captr(i-1),captr(i),xyz) .gt. mxcaca)
        else
          lstart = .true.
        end if
c
        if (i .eq. nres) then
          lend = .true.
        else if (captr(i+1) .gt. 0) then
          lend = (dist(captr(i+1),captr(i),xyz) .gt. mxcaca)
        else
          lend = .true.
        end if
c
        if (lend) then
c
          lstat (i) = -20
c
        else
c
          nok = nok + 1
          cadist = dist(captr(i+1),captr(i),xyz)
          rbuf (nok) = cadist
          if (cadist .le. 2.80) then
            dcnts(1)=dcnts(1)+1
            lbadd = .true.
          else if (cadist .le. 3.00) then
            dcnts(2)=dcnts(2)+1
          else if (cadist .le. 3.70) then
            dcnts(3)=dcnts(3)+1
            lbadd = .true.
          else if (cadist .le. 3.90) then
            dcnts(4)=dcnts(4)+1
          else
            dcnts(5)=dcnts(5)+1
            lbadd = .true.
          end if
c
          if (lstart) goto 100
c
          if (captr(i+2) .le. 0) goto 100
          if (dist(captr(i+1),captr(i+2),xyz) .gt. mxcaca) goto 100
c
          cangle = angle (captr(i-1),captr(i),captr(i+1),xyz)
c
          if (captr(i+3) .le. 0) goto 100
          if (dist(captr(i+2),captr(i+3),xyz) .gt. mxcaca) goto 100
c
          natok = natok + 1
          cators = tangle (captr(i-1),captr(i),captr(i+1),
     +                     captr(i+2),xyz)
c
          if (.not. lgly) then
            ncheck = ncheck + 1
            i1 = max (0, min (60, nint (cangle/3.0) ) )
            i2 = max (-60, min (60, nint (cators/3.0) ) )
            ii = caqual (i1,i2)
            atcnts(ii) = atcnts(ii) + 1
            lbadat = (ii .eq. 0 .and. (.not. lgly))
          end if
c
          if (mode .eq. 2) then
            if (lgly) then
              call xps_colour (0)
              call xps_symbol (0,cangle-dx,cangle+dx,
     +          cators-dy,cators+dy)
            else
              if (lbadat) then
                call xps_colour (4)
                call xps_symbol (1,cangle-dx,cangle+dx,
     +            cators-dy,cators+dy)
                call xps_symbol (2,cangle-dx,cangle+dx,
     +            cators-dy,cators+dy)
                if (label) then
                  kk = atmptr(1,i)
                  write (line,6000) resnam(kk),achain(kk),iresid(kk),
     +              insert(kk),inote(kk)(7:10)
                  call remspa(line)
                  call xps_text (cangle+dx+0.5,cators-dy-0.5,
     +                           ptsize,line)
                end if
              else
                call xps_colour (0)
                call xps_symbol (1,cangle-dx,cangle+dx,
     +            cators-dy,cators+dy)
              end if
            end if
          end if
c
        end if
c
  100   continue
c
        lp = .false.
c
        if (lstart) then
          call print_res(i,2)
          lp = .true.
          call prompt (' <<<<< Start of new chain')
          if (mode .eq. 1) write (iunit,*,err=9000)
        end if
c
        if (lstart .and. lend) then
          call prompt (' Warning - isolated CA atom !!!')
          niso = niso + 1
        end if
c
        if (lend) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' >>>>> End of chain')
        end if
c
        if (lbadd) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          write (line,'(a,f8.2)') ' Warning - unusual CA-CA distance ',
     +      cadist
          call prompt (line)
        end if
c
        if (lbadat) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' Warning - unusual CA geometry !')
          lstat (i) = 1
        end if
c
        if (mode .eq. 1) then
          kk = atmptr(1,i)
          write (iunit,6200,err=9000)
     +      resnam(kk),achain(kk),iresid(kk),
     +      insert(kk),inote(kk)(7:10),cadist,cangle,cators
        end if
c
 1000   continue
      end do
c
      write (*,*)
      call textut (' Chain ID checked (* = all) :',mychn)
      call jvalut (' Nr of residues checked     :',1,nres)
      call jvalut (' Nr of isolated CAs         :',1,niso)
      call jvalut (' Nr with CA-CA distance     :',1,nok)
      call jvalut (' Nr with CA ang/torsion     :',1,natok)
      call jvalut ('   Ditto, non-Gly           :',1,ncheck)
c
      if (nok .gt. 2) then
        call prompt ('0CA-CA distances:')
        call xstats (rbuf,nok,ave,sdv,xmin,xmax,xtot)
        write (*,7000) nok,ave,sdv,xmin,xmax
c
        x = 100.0*float(dcnts(1))/float(nok)
        y = abs(x-0.008)/0.079
        write (*,7010) 'Short (<= 2.8 A) ',
     +      dcnts(1),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(2))/float(nok)
        y = abs (x-0.240)/0.458
        write (*,7010) 'CIS   (<= 3.0 A) ',
     +      dcnts(2),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(3))/float(nok)
        y = abs (x-1.517)/3.492
        write (*,7010) 'Poor  (<= 3.7 A) ',
     +      dcnts(3),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(4))/float(nok)
        y = abs (x-96.818)/7.044
        write (*,7010) 'TRANS (<= 3.9 A) ',
     +      dcnts(4),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
        x = 100.0*float(dcnts(5))/float(nok)
        y = abs (x-1.416)/4.004
        write (*,7010) 'Long  (>  3.9 A) ',
     +      dcnts(5),x,y
        if (y .gt. 3.0) call prompt (
     +      ' >>> WARNING - 3 SIGMA deviant !')
c
      end if
c
      if (ncheck .gt. 2) then
        call prompt ('0CA geometry:')
c
        x=100.0*float(atcnts(3))/float(ncheck)
        y=abs(x-72.8)/8.9
        write (line,5010) 'CORE',atcnts(3),x,y
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        if (y .gt. 3.0) then
          line = ' >>> WARNING - 3 SIGMA deviant !'
          call prompt (line)
        end if
c
        write (line,5000) 'Additional',atcnts(2),
     +    100.0*float(atcnts(2))/float(ncheck)
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
c
        write (line,5000) 'Generous',atcnts(1),
     +    100.0*float(atcnts(1))/float(ncheck)
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
c
        x=100.0*float(atcnts(0))/float(ncheck)
        y=abs(x-3.1)/2.2
        write (line,5010) 'DISALLOWED',atcnts(0),x,y
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        if (y .gt. 3.0) then
          line = ' >>> WARNING - 3 SIGMA deviant !'
          call prompt (line)
        end if
c
        line = ' For <= 2.0 A structures :'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Core       -  7.1 % area - 72.8 % (8.9) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Additional -  5.2 % area - 12.8 % (4.0) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Generous   - 15.0 % area - 11.3 % (4.6) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
        line = ' Disallowed - 72.6 % area -  3.1 % (2.2) residues'
        if (mode .eq. 2) call xps_legend (line)
        call prompt (line)
c
      end if
c
c ... check for stretches of poor residues
c
      if (atcnts(0) .ge. 2) then
c
        i = 0
 3241   continue
        i = i + 1
        if (i .ge. nres) goto 3242
c
c ... 5 out of 8 sequential ?
c
          if (i .le. (nres-7)) then
            j = lstat(i) + lstat(i+1) + lstat(i+2) +
     +          lstat(i+3) + lstat(i+4) + lstat(i+5) +
     +          lstat(i+6) + lstat(i+7)
            if (j .ge. 5) then
              write (*,*)
              call prompt (
     +  ' WARNING !!!  At least five of eight sequential residues')
              call prompt (
     +  '              have poor CA geometry, starting at:')
              call print_res (i,0)
              i = i + 7
              goto 3241
            end if
          end if
c
c ... 3 out of 5 sequential ?
c
          if (i .le. (nres-4)) then
            j = lstat(i) + lstat(i+1) + lstat(i+2) +
     +          lstat(i+3) + lstat(i+4)
            if (j .ge. 3) then
              write (*,*)
              call prompt (
     +  ' WARNING !!!  At least three of five sequential residues')
              call prompt (
     +  '              have poor CA geometry, starting at:')
              call print_res (i,0)
              i = i + 4
              goto 3241
            end if
          end if
c
c ... 2 out of 3 sequential ?
c
          if (i .le. (nres-2)) then
            j = lstat(i) + lstat(i+1) + lstat(i+2)
            if (j .ge. 2) then
              write (*,*)
              call prompt (
     +  ' WARNING !!!  At least two of three sequential residues')
              call prompt (
     +  '              have poor CA geometry, starting at:')
              call print_res (i,0)
              i = i + 2
              goto 3241
            end if
          end if
c
        goto 3241
c
 3242   continue
c
      end if
c
 5000 format (1x,a10,' res : ',i6,' = ',f8.2,' %')
 5010 format (1x,a10,' res : ',i6,' = ',f8.2,' % = ',f4.1,
     +  ' SIGMA from mean')
c
 7000 format (1x,i7,' CA-CA distances'/
     +  ' Average CA-CA distance = ',f8.3,' Sigma = ',f8.3/
     +  ' Minimum CA-CA distance = ',f8.3,' Maxim = ',f8.3)
 7010 format (1x,a,' CA-CA dists : ',i6,' res = ',f8.2,' % = ',
     +  f4.1,' SIGMA from mean')
c
 8000 continue
c
      if (mode .eq. 1) then
        close (iunit)
      else if (mode .eq. 2) then
        if (natok .gt. 0) then
          call xps_close ()
          call prompt (' PostScript file created')
        else
          call xps_delete ()
          call prompt (' PostScript file empty and deleted')
        end if
      end if
c
      return
c
 9000 continue
      call errcon ('While writing file')
      return
c
      end
c
c
c
      subroutine mole_stats ()
c
      include 'moleman2.incl'
c
      real xn,yn,zn,ave,sdv,xmin,xmax,rog,xtot,xrms,xhave
      real xdim,ydim,zdim
c
      integer i,nh,ns,na
c
      logical lhydro
c
code ...
c
      call jvalut (' Nr of atoms    :',1,natoms)
      call jvalut (' Nr of residues :',1,nres)
      if (natoms .le. 0) return
c
      write (*,*)
      call jvalut (' Nr of amino acid residues :',
     +  1,icnt(iprot))
      call jvalut (' Nr of nucleic acids       :',
     +  1,icnt(inucl))
      call jvalut (' Nr of waters              :',
     +  1,icnt(iwate))
      call jvalut (' Nr of metals              :',
     +  1,icnt(imeta))
      call jvalut (' Nr of inorganics          :',
     +  1,icnt(iinor))
      call jvalut (' Nr of carbohydrates       :',
     +  1,icnt(icarb))
      call jvalut (' Nr of organic compounds   :',
     +  1,icnt(iorga))
      call jvalut (' Nr of other compounds     :',
     +  1,icnt(ihete))
c
 8702 format (' ',a10,6a11)
 8704 format (' ',a10,6f11.3)
c
      nh = 0
      ns = 0
      na = 0
      do i=1,natoms
        if (select(i)) then
          ns = ns + 1
          rbuf (ns) = xyz(1,i)
          rbuf (natoms+ns) = xyz(2,i)
          rbuf (2*natoms+ns) = xyz(3,i)
          rbuf (3*natoms+ns) = batom(i)
          rbuf (4*natoms+ns) = qatom(i)
          if (lhydro(atmnam(i))) nh = nh + 1
          if (laniso(i)) na = na + 1
        end if
      end do
c
      write (*,*)
      call jvalut (' Nr of selected atoms :',1,ns)
      call jvalut ('      Ditto, hydrogen :',1,nh)
      call jvalut ('      Ditto, ANISOU   :',1,na)
      if (ns .le. 0) return
c
      write (*,*)
      write (*,8702) 'Item','Average','St.Dev','Min','Max',
     +  'RMS','Harm.ave.'
      write (*,8702) '----','-------','------','---','---',
     +  '---','---------'
c
      call xstats (rbuf(1),ns,ave,sdv,xmin,xmax,xtot)
      write (*,8704) 'X-coord',ave,sdv,xmin,xmax
      xn = ave
      xdim = xmax - xmin
      call xstats (rbuf(natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      write (*,8704) 'Y-coord',ave,sdv,xmin,xmax
      yn = ave
      ydim = xmax - xmin
      call xstats (rbuf(2*natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      write (*,8704) 'Z-coord',ave,sdv,xmin,xmax
      zn = ave
      zdim = xmax - xmin
c
      call xstats (rbuf(3*natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      call xstat2 (rbuf(3*natoms+1),ns,xrms,xhave)
      write (*,8704) 'B-factor',ave,sdv,xmin,xmax,xrms,xhave
      if (xmin .lt. 2.0) then
        call prompt (' Warning - there are B-factors < 2.0 A**2 !')
      end if
      if (xmax .gt. 100.0) then
        call prompt (' Warning - there are B-factors > 100 A**2 !')
      end if
c
      call xstats (rbuf(4*natoms+1),ns,ave,sdv,xmin,xmax,xtot)
      call xstat2 (rbuf(4*natoms+1),ns,xrms,xhave)
      write (*,8704) 'Occpncy',ave,sdv,xmin,xmax,xrms,xhave
      if (xmin .lt. 0.0) then
        call prompt (' Warning - there are occupancies < 0 !')
      end if
      if (xmax .gt. 1.0) then
        call prompt (' Warning - there are occupancies > 1 !')
      end if
c
      rog = 0.0
      do i=1,ns
        rog = rog + (rbuf(i)-xn)**2 +
     +              (rbuf(natoms+i)-yn)**2 +
     +              (rbuf(2*natoms+i)-zn)**2
      end do
      rog = max (0.0, rog / float(ns))
      rog = sqrt (rog)
      write (*,6010) rog
 6010 format (/' The radius of gyration is ',f6.1,' A')
c
      write (*,6000) xdim,ydim,zdim
 6000 format (/' Range of X, Y, and Z coordinates: ',
     +  f6.1,' A * ',f6.1,' A * ',f6.1,' A'/
     +  ' If you have used XYz ALign_inertia_axes, these numbers'/
     +  ' give you an indication of the dimensions of the selected'/
     +  ' molecule (or set of atoms).')
c
c ... more (centre-of-mass; sum of masses ?)
c
      return
      end
c
c
c
      subroutine geom_multi (mytyp)
c
      include 'moleman2.incl'
c
      real dist,angle,tangle,xave,xsdv,xmin,xmax,xtot,cosave,sinave
c
      integer iptr(maxcpy,maxapr)
      integer i,j,nok,ne,k,l,m,ilib,nlib,nn
c
      logical lbond(maxapr,maxapr)
      logical ldum,select_res
c
      character mytyp*3,text*11
c
code ...
c
      call textut (' Multiple copy geometry for :',mytyp)
      ilib = 0
      do i=1,nmrtyp
        if (lrname(i) .eq. mytyp) then
          ilib = i
          goto 10
        end if
      end do
   10 continue
      if (ilib .le. 0) then
        call errcon ('Residue type not in library')
        return
      end if
      nlib = nmrptr(2,ilib) - nmrptr(1,ilib) + 1
      call jvalut (' Nr of atoms :',1,nlib)
      call asciut (' Atoms :',nlib,lratom(nmrptr(1,ilib)))
      if (nlib .gt. maxapr) then
        call errcon ('Residue type has too many atoms')
        call jvalut (' Maximum :',1,maxapr)
        return
      end if
      if (nlib .le. 1) then
        call errcon ('Residue type has fewer than two atoms')
        return
      end if
c
      call prompt (' Looking for selected residues ...')
      nok = 0
      do i=1,nres
        if (resnam(atmptr(1,i)) .eq. mytyp) then
          if (select_res(i)) then
            nok = nok + 1
            if (nok .gt. maxcpy) then
              call errcon ('Too many copies')
              call jvalut (' Maximum :',1,maxcpy)
              nok = maxcpy
              goto 40
            end if
            goto 20
          end if
          goto 30
c
   20     continue
          call print_res (i,0)
          do j=1,nlib
            iptr (nok,j) = 0
ccc         print *,' ... ',lratom(nmrptr(1,ilib)+j-1)
            do k=atmptr(1,i),atmptr(2,i)
              if (lratom(nmrptr(1,ilib)+j-1) .eq.
     +            atmnam(k)) then
                iptr (nok,j) = k
ccc             print *,'     ',atmnam(k)
              end if
            end do
          end do
c
   30     continue
        end if
      end do
   40 continue
c
      call jvalut (' Nr of copies found :',1,nok)
      if (nok .lt. 2) then
        call errcon ('Fewer than two copies found')
        return
      end if
c
c ... two atoms are bonded if they are in bonding distance
c     in at least one of the copies
c
      ne = 0
      do i=1,nlib-1
        do j=i+1,nlib
          ldum = .false.
          do k=1,nok
            if (iptr(k,i) .gt. 0 .and.
     +          iptr(k,j) .gt. 0) then
              if (dist(iptr(k,i),iptr(k,j),xyz) .le. mxbond) then
                ldum = .true.
                ne = ne + 1
                goto 50
              end if
            end if
          end do
   50     continue
          lbond (i,j) = ldum
          lbond (j,i) = ldum
        end do
      end do
      call jvalut (' Nr of bonds :',1,ne)
      if (ne .lt. 1) return
c
      call fvalut (' Bond distance range large if >',1,largeb)
      call fvalut (' Bond angle    range large if >',1,largea)
c
c ... do all bonds
c
      write (*,6000) mxbond
      do i=1,nlib-1
        do j=i+1,nlib
          if (.not. lbond(i,j)) goto 60
          nn = 0
          do k=1,nok
            if (iptr(k,i) .gt. 0 .and.
     +          iptr(k,j) .gt. 0) then
              nn = nn + 1
              rbuf (nn) = dist(iptr(k,i),iptr(k,j),xyz)
            end if
          end do
          if (nn .le. 0) then
            xave = 0.0
            xsdv = 0.0
            xmin = 0.0
            xmax = 0.0
          else if (nn .eq. 1) then
            xave = rbuf (1)
            xsdv = 0.0
            xmin = rbuf (1)
            xmax = rbuf (1)
          else
            call xstats (rbuf(1),nn,xave,xsdv,xmin,xmax,xtot)
          end if
          text = ' '
          if ( (xmax-xmin) .gt. largeb) text='Large range'
          write (*,6010) lratom(nmrptr(1,ilib)+i-1),
     +      lratom(nmrptr(1,ilib)+j-1),nn,xave,xsdv,xmin,xmax,text
   60     continue
        end do
      end do
c
 6000 format (/' Bonded distances with cut-off : ',f8.3,' A'/
     +         ' ==========================================')
 6010 format (1x,a4,' - ',a4,' # ',i3,' Ave, Sdv, Min, Max ',
     +  4f8.3,1x,a)
c
      if (ne .lt. 2) return
c
c ... do all angles
c
      write (*,7000)
      do i=1,nlib
        do j=1,nlib-1
          if (.not. lbond(i,j)) goto 70
          do k=j+1,nlib
            if (.not. lbond(i,k)) goto 80
            nn = 0
            do l=1,nok
              if (iptr(l,i) .gt. 0 .and.
     +            iptr(l,j) .gt. 0 .and.
     +            iptr(l,k) .gt. 0) then
                nn = nn + 1
                rbuf (nn) = angle (iptr(l,j),iptr(l,i),iptr(l,k),xyz)
                rbuf (maxatm+nn) = dist (iptr(l,j),iptr(l,k),xyz)
              end if
            end do
c
            if (nn .le. 0) then
              xave = 0.0
              xsdv = 0.0
              xmin = 0.0
              xmax = 0.0
            else if (nn .eq. 1) then
              xave = rbuf (1)
              xsdv = 0.0
              xmin = rbuf (1)
              xmax = rbuf (1)
            else
              call xstats (rbuf(1),nn,xave,xsdv,xmin,xmax,xtot)
            end if
            text = ' '
            if ( (xmax-xmin) .gt. largea) text='Large range'
            write (*,7010) lratom(nmrptr(1,ilib)+j-1),
     +        lratom(nmrptr(1,ilib)+i-1),lratom(nmrptr(1,ilib)+k-1),
     +        nn,xave,xsdv,xmin,xmax,text
c
            if (nn .le. 0) then
              xave = 0.0
              xsdv = 0.0
              xmin = 0.0
              xmax = 0.0
            else if (nn .eq. 1) then
              xave = rbuf (maxatm+1)
              xsdv = 0.0
              xmin = rbuf (maxatm+1)
              xmax = rbuf (maxatm+1)
            else
              call xstats (rbuf(maxatm+1),nn,xave,xsdv,xmin,xmax,xtot)
            end if
            write (*,7020) nn,xave,xsdv,xmin,xmax
c
   80       continue
          end do
   70     continue
        end do
      end do
c
 7000 format (/' Angles and 1-3 angle distances'/
     +         ' ==============================')
 7010 format (1x,a4,' - ',a4,' - ',a4,' Angle     : ',i3,4f8.2,1x,a)
 7020 format (19x,' 1-3 Dist  : ',i3,4f8.3,1x,a)
c
      if (ne .lt. 3) return
c
c ... do all dihedrals
c
      write (*,8000)
      do i=1,nlib-1
        do j=i+1,nlib
          if (.not. lbond(i,j)) goto 100
          do k=1,nlib
            if (k.eq.i .or. k.eq.j) goto 110
            if (.not. lbond(k,i)) goto 110
            do l=1,nlib
              if (l.eq.i .or. l.eq.j .or. l.eq.k) goto 120
              if (.not. lbond(l,j)) goto 120
              nn = 0
              do m=1,nok
                if (iptr(m,i) .gt. 0 .and.
     +              iptr(m,j) .gt. 0 .and.
     +              iptr(m,k) .gt. 0 .and.
     +              iptr(m,l) .gt. 0) then
                  nn = nn + 1
                  rbuf (nn) = tangle (iptr(m,k),iptr(m,i),
     +                        iptr(m,j),iptr(m,l),xyz)
c
ccc                  if (rbuf(nn) .le. -140.) rbuf(nn)=rbuf(nn)+360.
ccc                  if (rbuf(nn) .gt.  240.) rbuf(nn)=rbuf(nn)-360.
c
c  (The proper way to average angles is to calculate ATAN2 ( <SIN>, <COS> )
c
                  rbuf (maxatm+nn) = dist (iptr(m,k),iptr(m,l),xyz)
                end if
              end do
c
              if (nn .le. 0) then
                xave = 0.0
                xsdv = 0.0
                xmin = 0.0
                xmax = 0.0
              else if (nn .eq. 1) then
                xave = rbuf (1)
                xsdv = 0.0
                xmin = rbuf (1)
                xmax = rbuf (1)
              else
c
                call xstats (rbuf(1),nn,xave,xsdv,xmin,xmax,xtot)
c
                cosave = 0.0
                sinave = 0.0
                do m=1,nn
                  cosave = cosave + cos (degtor*rbuf(m))
                  sinave = sinave + sin (degtor*rbuf(m))
                end do
                cosave = cosave / float (nn)
                sinave = sinave / float (nn)
                xave = rtodeg * atan2 (sinave,cosave)
              end if
c
ccc              text = ' '
ccc              if ( (xmax-xmin) .gt. larged) text='Large range'
c
              write (*,8010) lratom(nmrptr(1,ilib)+k-1),
     +          lratom(nmrptr(1,ilib)+i-1),lratom(nmrptr(1,ilib)+j-1),
     +          lratom(nmrptr(1,ilib)+l-1),nn,xave,xsdv,xmin,xmax
c
              if (nn .le. 0) then
                xave = 0.0
                xsdv = 0.0
                xmin = 0.0
                xmax = 0.0
              else if (nn .eq. 1) then
                xave = rbuf (maxatm+1)
                xsdv = 0.0
                xmin = rbuf (maxatm+1)
                xmax = rbuf (maxatm+1)
              else
                call xstats (rbuf(maxatm+1),nn,xave,xsdv,xmin,xmax,xtot)
              end if
              write (*,8020) nn,xave,xsdv,xmin,xmax
c
  120         continue
            end do
  110       continue
          end do
  100     continue
        end do
      end do
c
 8000 format (/' Dihedrals and 1-4 torsion distances'/
     +         ' ===================================')
 8010 format (1x,a4,' - ',a4,' - ',a4,' - ',a4,' Dihedral : ',
     +  i3,4f8.2,1x,a)
 8020 format (26x,' 1-4 Dist : ',i3,4f8.3,1x,a)
c
      return
      end
c
c
c
      subroutine geom_select ()
c
      include 'moleman2.incl'
c
      real dist,angle,tangle
c
      integer i,j,nok,maxp,np,nd,ne,ii,jj,k,kk,l,ll
c
code ...
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          ibuf (nok) = i
        end if
      end do
      call jvalut (' Nr of selected atoms :',1,nok)
c
      maxp = maxbuf
      np = int (sqrt(float(maxp)))
      if (nok .gt. np) then
        call errcon ('Too many atoms selected')
        call jvalut (' Maximum allowed :',1,np)
        return
      end if
c
      nd = 0
      do i=1,nok-1
        ii = ibuf(i)
        do j=i+1,nok
          jj = ibuf(j)
          nd = (i-1)*nok + j
          rbuf (nd) = dist (ii,jj,xyz)
          ne = (j-1)*nok + i
          rbuf (ne) = rbuf (nd)
        end do
      end do
c
c ... bonded distances
c
      write (*,6000) mxbond
      ne = 0
      do i=1,nok-1
        ii = ibuf(i)
        do j=i+1,nok
          jj = ibuf(j)
          nd = (i-1)*nok + j
          if (rbuf(nd) .le. mxbond) then
            ne = ne + 1
            write (*,6010) atmnam(ii),achain(ii),iresid(ii),insert(ii),
     +        atmnam(jj),achain(jj),iresid(jj),insert(jj),rbuf(nd)
          end if
        end do
      end do
      call jvalut (' Nr of bonded distances :',1,ne)
 6000 format (/' Bonded distances with cut-off : ',f8.3,' A'/
     +         ' ==========================================')
 6010 format (1x,(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' = ',f8.3,' A')
c
      if (ne .lt. 2) return
c
c ... angles
c
      write (*,7000)
      ne = 0
      do i=1,nok
        ii = ibuf (i)
        do j=1,nok-1
          jj = ibuf (j)
          nd = (i-1)*nok + j
          if (rbuf(nd) .le. mxbond .and. (ii.ne.jj)) then
            do k=j+1,nok
              kk = ibuf (k)
              if (kk .ne. ii .and. kk .ne. jj) then
                nd = (i-1)*nok + k
                if (rbuf(nd) .le. mxbond) then
                  ne = ne + 1
                  write (*,7010)
     +        atmnam(jj),achain(jj),iresid(jj),insert(jj),
     +        atmnam(ii),achain(ii),iresid(ii),insert(ii),
     +        atmnam(kk),achain(kk),iresid(kk),insert(kk),
     +        angle(jj,ii,kk,xyz),dist(jj,kk,xyz)
                end if
              end if
            end do
          end if
        end do
      end do
      call jvalut (' Nr of angles :',1,ne)
 7000 format (/' Angles and 1-3 angle distances'/
     +         ' ==============================')
 7010 format (1x,(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' - ',(a4,' [',a1,i4,a1,']'),' = ',f8.3,' deg = ',f8.3,' A')
c
c ... dihedrals
c
      write (*,8000)
      ne = 0
      do i=1,nok-1
        ii = ibuf (i)
        do j=i+1,nok
          jj = ibuf (j)
          nd = (i-1)*nok + j
          if (rbuf(nd) .le. mxbond .and. (ii.ne.jj)) then
            do k=1,nok
              kk = ibuf (k)
              if (kk .ne. ii .and. kk .ne. jj) then
                nd = (i-1)*nok + k
                if (rbuf(nd) .le. mxbond) then
                  do l=1,nok
                    ll = ibuf (l)
                    if (ll .ne. ii .and. ll .ne. jj .and.
     +                  ll .ne. kk) then
                      nd = (j-1)*nok + l
                      if (rbuf(nd) .le. mxbond) then
                        ne = ne + 1
                        write (*,8010)
     +        atmnam(kk),achain(kk),iresid(kk),insert(kk),
     +        atmnam(ii),achain(ii),iresid(ii),insert(ii),
     +        atmnam(jj),achain(jj),iresid(jj),insert(jj),
     +        atmnam(ll),achain(ll),iresid(ll),insert(ll),
     +        tangle(kk,ii,jj,ll,xyz),dist(kk,ll,xyz)
                      end if
                    end if
                  end do
                end if
              end if
            end do
          end if
        end do
      end do
      call jvalut (' Nr of dihedrals :',1,ne)
 8000 format (/' Dihedrals and 1-4 torsion distances'/
     +         ' ===================================')
 8010 format (1x,(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' - ',(a4,' [',a1,i4,a1,']'),' - ',(a4,' [',a1,i4,a1,']'),
     +  ' = ',f8.3,' deg = ',f8.3,' A')
c
      return
      end
c
c
c
      logical function select_res (ires)
c
      include 'moleman2.incl'
c
      integer i,ires
c
code ...
c
      select_res = .false.
c
      if (ires .gt. 0 .and. ires .le. nres) then
        do i=atmptr(1,ires),atmptr(2,ires)
          if (select(i)) then
            select_res = .true.
            return
          end if
        end do
      end if
c
      return
      end
c
c
c
      subroutine na_duarte (iunit,file,what,which)
c
      include 'moleman2.incl'
c
      real dist,tangle,dx,dy,ptsize,eta,theta
c
      integer i,nok,j,nend,nerr,nyes,isym,icol
      integer iunit,mode,length,ierr,kk,leng1
      integer na,nc,ng,nt,nu,nx,inu
c
      logical lstart,lend,leta,ltheta,lp,label
c
      character*(*) file,what,which
      character typ1lc(6)*1
      character line*256,mychn*1
c
      data typ1lc /'A','G','C','T','U','X'/
c
code ...
c
      mode = 0
      if (length(file) .gt. 0) then
        call upcase (what)
        call remspa (what)
        if (what(1:1) .eq. 'T') then
          mode = 1
        else if (what(1:1) .eq. 'A') then
          mode = 2
          label = .true.
        else
          mode = 2
          label = .false.
        end if
      end if
c
      mychn = which(1:1)
      call upcase (mychn)
c
      if (mode .eq. 1) then
c
        call xopxua (iunit,file,linter,ierr)
        if (ierr .ne. 0) return
c
      else if (mode .eq. 2) then
c
        dx = 2.0
        dy = 2.0
        ptsize = 10.0
c
        call psnini (iunit,file,prognm)
c
      end if
c
      do i=1,2*nres
        ibuf (i) = -1
      end do
c
      do i=1,2*nres
        rbuf (i) = -999.99
      end do
c
      write (*,*)
      call textut (' Chain ID to check (* = all) :',mychn)
c
      nok = 0
      do i=1,nres
        lbuf (i) = (restyp(i).eq.inucl)
        if (mychn .ne. '*') then
          lbuf(i) = ( lbuf(i) .and.
     +                (achain(atmptr(1,i)) .eq. mychn) )
        end if
        if (lbuf(i)) then
          nok = nok + 1
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' P  ') then
              ibuf (i) = j
            else if (atmnam(j) .eq. ' C4*') then
              ibuf (nres+i) = j
            else if (atmnam(j) .eq. ' C4''') then
              ibuf (nres+i) = j
            end if
          end do
        end if
      end do
c
      call jvalut (' Nr of residues to check :',1,nok)
      if (nok .lt. 3) then
        call errcon ('Fewer than 3 nucleic acid residues ...')
        call xps_delete ()
        call prompt (' PostScript file empty and deleted')
        return
      end if
c
      nend = 0
      nerr = 0
      nyes = 0
c
      na = 0
      nc = 0
      ng = 0
      nt = 0
      nu = 0
      nx = 0
c
 6000 format (a3,'-',a1,i4,a1,'-',a4)
 6100 format (1x,a15,4(1x,a10))
 6200 format (1x,a3,'-',a1,i4,a1,'-',a4,2(1x,f10.2),' [',a1,']')
 6010 format (/
     +  ' Listing of Eta and Theta pseudo-torsion'/
     +  ' angles for PDB file ',a,' chain ',a1/
     + /' ETA   = C4*(i-1) - P(i) - C4*(i) - P(i+1)'/
     +  ' THETA = P(i) - C4*(i) - P(i+1) - C4*(i+1)'/
     + /' An entry of "-999.9" means that the pseudo-torsion could'/
     +  ' not be calculated (for terminal residues and residues with'/
     +  ' missing atoms).'/)
c
      if (mode .eq. 1) then
        call stamp (line)
        write (iunit,'(1x,a)') line(1:leng1(line))
        write (iunit,6010,err=9000) pdbfil(1:leng1(pdbfil)),mychn
        write (iunit,6100,err=9000)
     +    'Residue','Eta','Theta'
        write (iunit,6100,err=9000)
     +    '-------','---','-----'
      end if
c
      do i=1,nres
c
        if (.not. lbuf(i)) goto 1000
c
        if (ibuf(i) .le. 0) goto 1000
c
        lstart = .false.
        lend = .false.
        leta = .false.
        ltheta = .false.
c
        eta = -999.9
        theta = -999.9
c
        if (i .eq. 1) then
          lstart = .true.
        else if (ibuf(i-1) .gt. 0) then
          lstart =  (dist(ibuf(i-1),ibuf(i),xyz) .gt. mxpp)
        else
          lstart = .true.
        end if
c
        if (i .eq. nres) then
          lend = .true.
        else  if (ibuf(i+1) .gt. 0) then
          lend =  (dist(ibuf(i+1),ibuf(i),xyz) .gt. mxpp)
        else
          lend = .true.
        end if
c
        if (.not. lend) then
c
c ... THETA = P(i) - C4*(i) - P(i+1) - C4*(i+1)
c
          ltheta = (ibuf(i) .gt. 0 .and.
     +              ibuf(nres+i) .gt. 0 .and.
     +              ibuf(i+1) .gt. 0 .and.
     +              ibuf(nres+i+1) .gt. 0)
          if (ltheta) then
            theta = tangle (ibuf(i),ibuf(nres+i),
     +              ibuf(i+1),ibuf(nres+i+1),xyz)
            call fix360 (theta)
          end if
c
c ... ETA   = C4*(i-1) - P(i) - C4*(i) - P(i+1)
c
          if (.not. lstart) then
            leta = (ibuf(nres+i-1) .gt. 0 .and.
     +              ibuf(i) .gt. 0 .and.
     +              ibuf(nres+i) .gt. 0 .and.
     +              ibuf(i+1) .gt. 0)
            if (leta) then
              eta = tangle (ibuf(nres+i-1),ibuf(i),
     +              ibuf(nres+i),ibuf(i+1),xyz)
              call fix360 (eta)
            end if
          end if
c
        end if
c
        if (ltheta .and. leta) then
          rbuf (i) = eta
          rbuf (nres+i) = theta
        end if
c
        lp = .false.
c
        if (lstart) then
          call print_res(i,2)
          lp = .true.
          call prompt (' <<<<< Start of new chain')
          if (mode .eq. 1) write (iunit,*,err=9000)
        end if
c
        if (lend) then
          if (.not. lp) call print_res(i,1)
          lp = .true.
          call prompt (' >>>>> End of chain')
        end if
c
        if (.not. lend) then
          if (.not. ltheta) then
            if (.not. lp) call print_res(i,1)
            lp = .true.
            call prompt (' ERROR - could not calculate THETA')
          end if
          if (.not. lstart) then
            if (.not. leta) then
              if (.not. lp) call print_res(i,1)
              lp = .true.
              call prompt (' ERROR - could not calculate ETA')
            end if
          end if
        end if
c
        if (lstart .or. lend) then
          nend = nend + 1
        else if (.not. (leta .and. ltheta)) then
          nerr = nerr + 1
        else
          nyes = nyes + 1
        end if
c
        if (ltheta .and. leta) then
c
          call nuctyp (resnam(atmptr(1,i)),inu)
c
          if (inu .eq. 1) then
            isym = 0
            icol = 0
            na = na + 1
          else if (inu .eq. 3) then
            isym = 1
            icol = 1
            nc = nc + 1
          else if (inu .eq. 2) then
            isym = 2
            icol = 6
            ng = ng + 1
          else if (inu .eq. 4) then
            isym = 4
            icol = 4
            nt = nt + 1
          else if (inu .eq. 5) then
            isym = 5
            icol = 5
            nu = nu + 1
          else
            isym = 3
            icol = 2
            nx = nx + 1
            inu = 6
          end if
c
          kk = atmptr(1,i)
          write (*,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),eta,theta,typ1lc(inu)
          if (mode .eq. 0) goto 1000
c
          if (mode .eq. 1) then
            write (iunit,6200,err=9000)
     +        resnam(kk),achain(kk),iresid(kk),
     +        insert(kk),inote(kk)(7:10),eta,theta,typ1lc(inu)
            goto 1000
          end if
c
          if (.not. label) then
            call xps_colour (0)
            call xps_symbol (1,eta-dx,eta+dx,theta-dy,theta+dy)
          else
            call xps_colour (icol)
            call xps_symbol (isym,eta-dx,eta+dx,theta-dy,theta+dy)
          end if
        end if
c
 1000   continue
      end do
c
      write (*,*)
      call textut (' Chain ID checked (* = all)   :',mychn)
      call jvalut (' Total nr of residues checked :',1,nok)
      call jvalut (' Start/end residues           :',1,nend)
      call jvalut (' Problems (missing atoms ?)   :',1,nerr)
      call jvalut (' Remaining residues in plot   :',1,nyes)
      call jvalut (' ... Number of A :',1,na)
      call jvalut (' ... Number of C :',1,nc)
      call jvalut (' ... Number of G :',1,ng)
      call jvalut (' ... Number of T :',1,nt)
      call jvalut (' ... Number of U :',1,nu)
      call jvalut (' ... Others      :',1,nx)
c
      if (mode .eq. 2) then
        line = ' PDB file : '//pdbfil
        call xps_legend (line)
c
        write (line,*) ' Residues in plot :',nyes
        call xps_legend (line)
c
        if (label) then
          line = ' A=black box, C=red plus, G=cyan cross,'//
     +      ' T=blue diamond, U=magenta triangle, other=green Z'
          call xps_legend (line)
        end if
c
        line = ' Shaded area has ETA in [150,190]'//
     +    ' and/or THETA in [190,260]'
        call xps_legend (line)
c
        line = ' See: CM Duarte & AM Pyle (1998). '//
     +         'J. Mol. Biol. 284, 1465-1478.'
        call xps_legend (line)
c
        if (mode .eq. 1) then
          close (iunit)
        else if (mode .eq. 2) then
          if (nyes .gt. 0) then
            call xps_close ()
            call prompt (' PostScript file created')
          else
            call xps_delete ()
            call prompt (' PostScript file empty and deleted')
          end if
        end if
      end if
c
      return
c
 9000 continue
      call errcon ('While writing text file')
c
      return
      end
c
c
c
      subroutine pdb_sanity_check ()
c
      include 'moleman2.incl'
c
c ... max nr of allowed alternative conformations for any atom
c
      integer maxcon
      parameter (maxcon=20)
c
      real sumocc(maxapr)
      real distce,dd
c
      integer basptr(maxapr),numcon(maxapr),conptr(maxcon,maxapr)
      integer i,j,k,m,n,nat,nuniq
c
code ...
c
      do i=1,nres
        nat = atmptr(2,i) - atmptr(1,i) + 1
        if (nat .gt. maxapr) then
          write (*,*)
          call errcon ('Too many atoms in residue')
          call jvalut (' Maximum allowed :',1,maxapr)
          call print_res (i,1)
          goto 900
        end if
        do j=1,nat
          basptr (j) = 0
          numcon (j) = 0
          sumocc (j) = 0.0
        end do
        nuniq = 1
        basptr (1) = atmptr(1,i)
        numcon (1) = 1
        conptr (1,1) = atmptr(1,i)
        sumocc (1) = qatom (atmptr(1,i))
        do j=atmptr(1,i)+1,atmptr(2,i)
          do k=1,nuniq
            if (atmnam(j) .eq. atmnam(basptr(k))) then
              numcon (k) = numcon (k) + 1
              conptr (numcon(k),k) = j
              sumocc (k) = sumocc (k) + qatom (j)
              goto 800
            end if
          end do
          nuniq = nuniq + 1
          basptr (nuniq) = j
          numcon (nuniq) = 1
          conptr (1,nuniq) = j
          sumocc (nuniq) = qatom (j)
  800     continue
        end do
c
        do k=1,nuniq
c
c ... check if occupancies sum to 1.0 (not neccessary for atoms on symmetry
c     axes etc.)
c
          j = nint (100*sumocc(k))
          if (j .ne. 100) then
            write (*,6000) numcon(k),sumocc(k)
            do m=1,numcon(k)
              call print_atom (conptr(m,k))
            end do
          end if
c
c ... if only one alt. loc., then normally the flag should be blank
c
          if (numcon(k) .eq. 1) then
            if (altloc(basptr(k)) .ne. ' ') then
              write (*,6010) altloc(basptr(k))
              call print_atom (basptr(k))
            end if
          end if
c
c ... if more than 3 alt. loc., issue a warning
c
          if (numcon(k) .gt. 3) then
            write (*,6015) numcon(k)
            do m=1,numcon(k)
              call print_atom (conptr(m,k))
            end do
          end if
c
c ... if more than one alt. loc., none of the flags should be blank
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)
              if (altloc(conptr(m,k)) .eq. ' ') then
                write (*,6020) numcon(k)
                call print_atom (conptr(m,k))
              end if
            end do
          end if
c
c ... if more than one alt. loc., none of the flags should be duplicated
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)-1
              do n=m+1,numcon(k)
                if (altloc(conptr(m,k)) .eq. altloc(conptr(n,k))) then
                  write (*,6030) altloc(conptr(m,k))
                  call print_atom (conptr(m,k))
                  call print_atom (conptr(n,k))
                end if
              end do
            end do
          end if
c
c ... if more than one alt. loc., the flags ought to be alphabetic
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)-1
              do n=m+1,numcon(k)
                if (ichar(altloc(conptr(n,k))) .lt.
     +              ichar(altloc(conptr(n,k)))) then
                  write (*,6040)
                  call print_atom (conptr(m,k))
                  call print_atom (conptr(n,k))
                end if
              end do
            end do
          end if
c
c ... if more than one alt. loc., they should not be too close in space
c
          if (numcon(k) .gt. 1) then
            do m=1,numcon(k)-1
              do n=m+1,numcon(k)
                dd = distce (xyz(1,conptr(m,k)),xyz(1,conptr(n,k)))
                if (dd .le. 0.05) then
                  write (*,6050) dd
                else if (dd .le. 0.2) then
                  write (*,6052) dd
                else if (dd .le. 1.0) then
                  write (*,6054) dd
                else
                  goto 700
                end if
                call print_atom (conptr(m,k))
                call print_atom (conptr(n,k))
  700           continue
              end do
            end do
          end if
c
        end do
c
  900   continue
      end do
c
 6000 format (/
     +  ' WARNING - OCCUPANCIES DO NOT SUM TO 1.00'/
     +  '           for the following atom, the occupancies of the ',i3/
     +  '           alternate locations add up to ',f5.2,' instead of'/
     +  '           1.00 (this can be okay if you are sure that the'/
     +  '           atom has partial occupancy, or if it lies in a'/
     +  '           special position, such as on a twofold axis):')
c
 6010 format (/
     +  ' WARNING - SINGLE LOCATION WITH NON-BLANK LABEL'/
     +  '           the following atom has only one location, but its'/
     +  '           alternate location label is "',a1,'" instead of'/
     +  '           blank (this can be okay in some cases, e.g. for'/
     +  '           waters that only interact with one of a set of'/
     +  '           alternative conformations of an amino acid'/
     +  '           residue):')
c
 6015 format (/
     +  ' WARNING - MORE THAN 3 ALTERNATE LOCATIONS'/
     +  '           the following atom has ',i3,' alternate locations'/
     +  '           (you may want to verify that these are supported'/
     +  '           by the electron density):')
c
 6020 format (/
     +  ' ERROR   - ALTERNATE LOCATION WITH BLANK LABEL'/
     +  '           the following atom occupies one of ',i3/
     +  '           alternate locations and should therefore have a'/
     +  '           non-blank label):')
c
 6030 format (/
     +  ' ERROR   - ALTERNATE LOCATIONS WITH IDENTICAL LABELS'/
     +  '           the following atoms have identical alternate'/
     +  '           location labels ("',a1,'"):')
c
 6040 format (/
     +  ' NOTE    - LABELS NOT IN ALPHABETICAL ORDER'/
     +  '           the following two atoms have alternate location'/
     +  '           labels that are not in alphabetical order (this'/
     +  '           is not strictly necessary either, but may help'/
     +  '           avoid problems with programs that assume that'/
     +  '           labels will be called A, B, C, etc.):')
c
 6050 format (/
     +  ' ERROR   - ALTERNATE LOCATIONS IN IDENTICAL POSITIONS'/
     +  '           the following two atoms have alternate locations'/
     +  '           but their distance is only ',f5.2,' A, suggesting'/
     +  '           that they are in identical positions in space and'/
     +  '           should be merged into a single location:')
c
 6052 format (/
     +  ' WARNING - ALTERNATE LOCATIONS IN ALMOST IDENTICAL POSITIONS'/
     +  '           the following two atoms have alternate locations'/
     +  '           but their distance is only ',f5.2,' A, suggesting'/
     +  '           that they are in almost identical positions in'/
     +  '           space and could be merged into a single location:')
c
 6054 format (/
     +  ' NOTE    - ALTERNATE LOCATIONS IN SIMILAR POSITIONS'/
     +  '           the following two atoms have alternate locations'/
     +  '           and their distance is ',f5.2,' A, suggesting that'/
     +  '           they are in similar positions in space and could'/
     +  '           perhaps be merged into a single location:')
c
c1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c         1         2         3         4         5         6         7         8         9        10
c
      return
      end
c
c
c
      subroutine rel_cont_order (rcocut)
c
      include 'moleman2.incl'
c
      real rcocut,rcosqr,dd,rco
c
      integer i,j,k,l,nrok,ncon,sumsep
c
      logical select_res
c
code ...
c
      call prompt (
     +  ' Calculate relative contact order of selected atoms')
      call fvalut (' Cut-off contact distance (A) :',1,rcocut)
c
      rcosqr = rcocut * rcocut
      nrok = 0
      ncon = 0
      sumsep = 0
c
      do i=1,nres-1
        if (select_res(i)) then
          nrok = nrok + 1
          do j=i+1,nres
            if (select_res(j)) then
              do k=atmptr(1,i),atmptr(2,i)
                if (.not. select(k)) goto 300
                do l=atmptr(1,j),atmptr(2,j)
                  if (.not. select(l)) goto 100
                  dd = (xyz(1,k)-xyz(1,l))**2
                  if (dd .gt. rcosqr) goto 100
                  dd = dd + (xyz(2,k)-xyz(2,l))**2
                  if (dd .gt. rcosqr) goto 100
                  dd = dd + (xyz(3,k)-xyz(3,l))**2
                  if (dd .gt. rcosqr) goto 100
c
                  ncon = ncon + 1
c
c ... assume sequence separation is |res_nr(atom_l)-res_nr(atom_k)|
c
                  sumsep = sumsep + abs(iresid(l) - iresid(k))
c
                  goto 200
c
  100             continue
                end do
  300           continue
              end do
  200         continue
            end if
          end do
        end if
      end do
c
      if (select_res(nres)) nrok = nrok + 1
      call jvalut (' Nr of selected residues :',1,nrok)
      call jvalut (' Nr of contacting pairs  :',1,ncon)
      call jvalut (' Sum of separations      :',1,sumsep)
c
      rco = float(sumsep)/(float(nrok)*float(ncon))
      call fvalut (' Relative contact order  :',1,rco)
c
      return
      end
c
c ... mole3_subs.f - XYZ & ONO & SQUENCE subroutines for MOLEMAN2
c
c ... i.e., subroutines that *DO* include 'moleman2.incl'
c
c
      subroutine frac_cart (mymode)
c
      include 'moleman2.incl'
c
      real mat(9),dummy(3)
c
      integer mode,mymode,i,nok
c
code ...
c
      mode = mymode
      if (mode .ne. 1) mode = 0
c
      if (mode .eq. 0) then
        call prompt (' Fractionalise coordinates')
        call orthog (cell,mat,1)
      else
        call prompt (' Orthogonalise coordinates')
        call orthog (cell,mat,0)
      end if
c
      call fvalut (' Cell axes   :',3,cell)
      call fvalut (' Cell angles :',3,cell(4))
      call print_rot (mat)
c
ccc      call matejd (mat)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          call mulmtx (mat,xyz(1,i),dummy,3,3,1)
          if (nok .eq. 1) then
            call print_atom (i)
            call fvalut (' Before :',3,xyz(1,i))
            call fvalut (' After  :',3,dummy)
          end if
          xyz(1,i) = dummy(1)
          xyz(2,i) = dummy(2)
          xyz(3,i) = dummy(3)
        end if
      end do
c
      call jvalut (' Nr of transformed atoms :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_ranrot ()
c
      include 'moleman2.incl'
c
      real mat(9),dummy(3)
c
      integer i,nok
c
code ...
c
      call ranrot (mat)
      call prompt ('0Random rotation matrix:')
      call print_rot (mat)
      call matejd (mat)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          call mulmtx (mat,xyz(1,i),dummy,3,3,1)
          if (nok .eq. 1) then
            call print_atom (i)
            call fvalut (' Before :',3,xyz(1,i))
            call fvalut (' After  :',3,dummy)
          end if
          xyz(1,i) = dummy(1)
          xyz(2,i) = dummy(2)
          xyz(3,i) = dummy(3)
        end if
      end do
c
      call jvalut (' Nr of transformed atoms :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_axis_rotate (xhow,xangle)
c
      include 'moleman2.incl'
c
      real rotmat(3,3),dumxyz(3),dx,det3,detx,angle
c
      integer i,j,nok,ierr
c
      character how*1,xhow*(*),xangle*(*)
c
code ...
c
      how = xhow
      call upcase (how)
c
      call str2r (xangle,angle,ierr)
      if (ierr .ne. 0) then
        call errcon ('Invalid angle')
        return
      end if
c
      call prompt ('0Rotate around X, Y or Z axis')
      call textut (' Rotation axis :',how)
      call fvalut (' Angle (deg)   :',1,angle)
c
      do i=1,3
        do j=1,3
          rotmat (i,j) = 0.0
        end do
        rotmat (i,i) = 1.0
      end do
c
      dx = degtor * angle
c
      if (how .eq. 'X') then
        rotmat (2,2) = cos(dx)
        rotmat (3,3) = rotmat (2,2)
        rotmat (3,2) = sin(dx)
        rotmat (2,3) = -rotmat (3,2)
ccc        write (*,'(a12,3(3f13.7,:,/,12x))')
ccc     +     ' X Matrix : ',((rotmat(i,j),j=1,3),i=1,3)
c
      else if (how .eq. 'Y') then
        rotmat (1,1) = cos(dx)
        rotmat (3,3) = rotmat (1,1)
        rotmat (3,1) = sin(dx)
        rotmat (1,3) = -rotmat (3,1)
ccc        write (*,'(a12,3(3f13.7,:,/,12x))')
ccc     +     ' Y Matrix : ',((rotmat(i,j),j=1,3),i=1,3)
c
      else if (how .eq. 'Z') then
        rotmat (2,2) = cos(dx)
        rotmat (1,1) = rotmat (2,2)
        rotmat (1,2) = sin(dx)
        rotmat (2,1) = -rotmat (1,2)
ccc        write (*,'(a12,3(3f13.7,:,/,12x))')
ccc     +     ' Z Matrix : ',((rotmat(i,j),j=1,3),i=1,3)
c
      else
        call errcon ('Unrecognised axis (not one of X, Y, Z) !')
        return
      end if
c
      detx = det3 (rotmat)
      call rvalut (' Determinant of rotation matrix :',1,detx)
c
      call print_rot (rotmat)
      call matejd (rotmat)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          call mulmtx (rotmat,xyz(1,i),dumxyz,3,3,1)
          if (nok .eq. 1) then
            call print_atom (i)
            call fvalut (' Before :',3,xyz(1,i))
            call fvalut (' After  :',3,dumxyz)
          end if
          xyz(1,i) = dumxyz(1)
          xyz(2,i) = dumxyz(2)
          xyz(3,i) = dumxyz(3)
        end if
      end do
c
      call jvalut (' Nr of transformed atoms :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_rotate (how,a1,a2,a3)
c
      include 'moleman2.incl'
c
      real mat(9),dummy(3),det3,det
c
      integer i,nok,ierr
c
      character*(*) how,a1,a2,a3
c
code ...
c
      call str2r (a1,dummy(1),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a2,dummy(2),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a3,dummy(3),ierr)
      if (ierr .ne. 0) return
c
      if (how(1:2) .eq. 'CP') then
        call fvalut (' CCP4 Polar angles :',3,dummy)
        call ccppol (dummy,mat)
      else if (how(1:2) .eq. 'ME') then
        call fvalut (' Merlot Euler angles :',3,dummy)
        call mereul (dummy,mat)
      else if (how(1:2) .eq. 'MP') then
        call fvalut (' Merlot Polar angles :',3,dummy)
        call merpol (dummy,mat)
      else if (how(1:2) .eq. 'XP') then
        call fvalut (' X-PLOR Polar angles :',3,dummy)
        call xplpol (dummy,mat)
      else if (how(1:2) .eq. 'XL') then
        call fvalut (' X-PLOR Lattman Euler angles :',3,dummy)
        call eulatt (dummy,mat)
      else
        call fvalut (' CCP4 Euler angles :',3,dummy)
        call ccpeul (dummy,mat)
      end if
c
      call print_rot (mat)
      call matejd (mat)
      det = det3 (mat)
      write (*,'(1x,a,f15.6)') 'Determinant :',det
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          call mulmtx (mat,xyz(1,i),dummy,3,3,1)
          if (nok .eq. 1) then
            call print_atom (i)
            call fvalut (' Before :',3,xyz(1,i))
            call fvalut (' After  :',3,dummy)
          end if
          xyz(1,i) = dummy(1)
          xyz(2,i) = dummy(2)
          xyz(3,i) = dummy(3)
        end if
      end do
c
      call jvalut (' Nr of transformed atoms :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_trans (a1,a2,a3)
c
      include 'moleman2.incl'
c
      real dummy(3)
      real dd
c
      integer i,nok,ierr
c
      character*(*) a1,a2,a3
c
code ...
c
      call str2r (a1,dummy(1),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a2,dummy(2),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a3,dummy(3),ierr)
      if (ierr .ne. 0) return
c
      dd = sqrt (dummy(1)**2 + dummy(2)**2 + dummy(3)**2)
      call fvalut (' Translation :',3,dummy)
      call fvalut (' Length of translation vector :',1,dd)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          xyz(1,i) = xyz(1,i) + dummy(1)
          xyz(2,i) = xyz(2,i) + dummy(2)
          xyz(3,i) = xyz(3,i) + dummy(3)
        end if
      end do
c
      call jvalut (' Nr of transformed atoms :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_perturb (a1,a2,a3,a4,a5)
c
      include 'moleman2.incl'
c
      real dummy(5),rmses(6),avers(5),shift(5)
      real dp
c
      integer i,nok,ierr,j
c
      character*(*) a1,a2,a3,a4,a5
c
code ...
c
      call str2r (a1,dummy(1),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a2,dummy(2),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a3,dummy(3),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a4,dummy(4),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a5,dummy(5),ierr)
      if (ierr .ne. 0) return
c
      do i=1,5
        dummy(i) = max(0.0,abs(dummy(i)))
        rmses(i) = 0.0
        avers(i) = 0.0
      end do
      rmses (6) = 0.0
      call rvalut (' Max allowed shifts X/Y/Z/B/Q :',5,dummy)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          do j=1,5
            shift (j) = 0.0
            if (dummy(j) .gt. 0.0) then
              call gkrand (shift(j),-dummy(j),dummy(j),0)
            end if
            avers (j) = avers (j) + shift (j)
            rmses (j) = rmses (j) + (shift(j)*shift(j))
          end do
          dp = (shift(1)**2) + (shift(2)**2) + (shift(3)**2)
          rmses (6) = rmses (6) + sqrt(dp)
          do j=1,3
            xyz(j,i) = xyz(j,i) + shift(j)
          end do
          batom (i) = batom (i) + shift(4)
          qatom (i) = qatom (i) + shift(5)
        end if
      end do
c
      call jvalut (' Nr of perturbed atoms :',1,nok)
      if (nok .lt. 1) return
c
      do i=1,5
        avers(i) = avers(i) / float(nok)
        rmses(i) = sqrt (rmses(i) / float(nok))
      end do
      rmses (6) = rmses(6)/float(nok)
c
      call rvalut (' Average     shifts X/Y/Z/B/Q :',5,avers)
      call rvalut (' RMS         shifts X/Y/Z/B/Q :',5,rmses)
      call rvalut (' AVERAGE POSITIONAL SHIFT (A) :',1,rmses(6))
c
      return
      end
c
c
c
      subroutine xyz_mirror (plane,which)
c
      include 'moleman2.incl'
c
      real dummy
c
      integer i,nok,ierr,j
c
      character*(*) plane,which
      character line*80
c
code ...
c
      call upcase (plane)
      call remspa (plane)
      if (plane(1:1) .ne. 'Y' .and. plane(1:1) .ne. 'Z') then
        plane = 'X'
      end if
c
      call str2r (which,dummy,ierr)
      if (ierr .ne. 0) return
c
      write (line,'(3a,f10.3)') ' Mirror selected atoms in plane: ',
     +  plane(1:1),' = ',dummy
      call pretty (line(2:))
      call prompt (line)
c
      j = 1
      if (plane(1:1) .eq. 'Y') then
        j = 2
      else if (plane(1:1) .eq. 'Z') then
        j = 3
      end if
c
c ... mirror: say in X=3 plane, then 4 -> 2
c     -> equation: Xnew = Xmir - (Xold - Xmir) = 2 * Xmir - Xold
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          xyz(j,i) = 2.0*dummy - xyz(j,i)
        end if
      end do
c
      call jvalut (' Nr of atoms mirrored :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_invert (a1,a2,a3)
c
      include 'moleman2.incl'
c
      real dummy(3)
c
      integer i,nok,ierr
c
      character*(*) a1,a2,a3
c
code ...
c
      call str2r (a1,dummy(1),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a2,dummy(2),ierr)
      if (ierr .ne. 0) return
c
      call str2r (a3,dummy(3),ierr)
      if (ierr .ne. 0) return
c
      call fvalut (' Invert through the point :',3,dummy)
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          xyz(1,i) = 2.0*dummy(1) - xyz(1,i)
          xyz(2,i) = 2.0*dummy(2) - xyz(2,i)
          xyz(3,i) = 2.0*dummy(3) - xyz(3,i)
        end if
      end do
c
      call jvalut (' Nr of atoms inverted :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_user_rt (elms,ldist)
c
      include 'moleman2.incl'
c
      real rtatob(12),dum(3)
c
      integer i,nok,ierr
c
      logical ldist
c
      character*(*) elms(12)
c
code ...
c
      do i=1,12
        call str2r (elms(i),rtatob(i),ierr)
        if (ierr .ne. 0) return
      end do
c
      call anancs (1,rtatob,.true.,ierr)
      if (ierr .ne. 0) then
        if (ldist) then
          call prompt (' Operator applied at your own risk !')
        else
          call errcon ('Operation aborted')
          return
        end if
      end if
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          call vecrtv (xyz(1,i),dum,1,rtatob(1),rtatob(10))
          xyz (1,i) = dum (1)
          xyz (2,i) = dum (2)
          xyz (3,i) = dum (3)
        end if
      end do
c
      call jvalut (' Nr of atoms transformed :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_o_rt (iunit,file)
c
      include 'moleman2.incl'
c
      real rtatob(12),dum(3)
c
      integer i,nok,ierr,iunit
c
      character*(*) file
c
code ...
c
      i = 0
      call rdoncs (iunit,file,i,1,rtatob,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading LSQ operator file')
        return
      end if
      if (i .ne. 1) then
        call errcon ('No LSQ operator read')
        return
      end if
c
      call anancs (1,rtatob,.false.,ierr)
      if (ierr .ne. 0) then
        call errcon ('Operation aborted')
        return
      end if
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          call vecrtv (xyz(1,i),dum,1,rtatob(1),rtatob(10))
          xyz (1,i) = dum (1)
          xyz (2,i) = dum (2)
          xyz (3,i) = dum (3)
        end if
      end do
c
      call jvalut (' Nr of atoms transformed :',1,nok)
c
      return
      end
c
c
c
      subroutine xyz_origin ()
c
      include 'moleman2.incl'
c
      real dum(3)
c
      integer i,j,nok
c
code ...
c
      call prompt (' Moving CofG of selected atoms to (0,0,0)')
      nok = 0
      do i=1,3
        dum(i) = 0.0
      end do
      do j=1,natoms
        if (select(j)) then
          nok = nok + 1
          dum(1)=dum(1) + xyz(1,j)
          dum(2)=dum(2) + xyz(2,j)
          dum(3)=dum(3) + xyz(3,j)
        end if
      end do
c
      call jvalut (' Nr of selected atoms :',1,nok)
      if (nok .lt. 1) return
      do i=1,3
        dum(i) = dum(i) / float(nok)
      end do
      call fvalut (' Centre-of-Gravity :',3,dum)
c
      do j=1,natoms
        if (select(j)) then
          xyz(1,j)=xyz(1,j) - dum(1)
          xyz(2,j)=xyz(2,j) - dum(2)
          xyz(3,j)=xyz(3,j) - dum(3)
        end if
      end do
      call prompt (' CofG now at (0,0,0)')
c
      return
      end
c
c
c
      subroutine xyz_mass_origin ()
c
      include 'moleman2.incl'
c
      real dum(3),sumass,mass
c
      integer i,j,nok,nerr
c
code ...
c
      call prompt (' Moving CofG of selected atoms to (0,0,0)')
      call prompt (' Using atom masses to calculate CofG')
      nok = 0
      do i=1,3
        dum(i) = 0.0
      end do
      sumass = 0.0
      nerr = 0
      do j=1,natoms
        if (select(j)) then
          nok = nok + 1
          call getelm (atmnam(j),inote(j)(11:12),i,mass)
          if (i .le. 0) then
            nerr = nerr + 1
            call errcon ('Unknown chemical element')
            call print_atom (j)
          else
            dum(1)=dum(1) + mass*xyz(1,j)
            dum(2)=dum(2) + mass*xyz(2,j)
            dum(3)=dum(3) + mass*xyz(3,j)
            sumass = sumass + mass
          end if
        end if
      end do
c
      call jvalut (' Nr of selected atoms   :',1,nok)
      call jvalut (' Nr of unknown elements :',1,nerr)
      call rvalut (' Sum of masses          :',1,sumass)
      if (nok .lt. 1 .or. nerr .gt. 0) return
c
      do i=1,3
        dum(i) = dum(i) / sumass
      end do
      call fvalut (' Centre-of-Mass :',3,dum)
c
      do j=1,natoms
        if (select(j)) then
          xyz(1,j)=xyz(1,j) - dum(1)
          xyz(2,j)=xyz(2,j) - dum(2)
          xyz(3,j)=xyz(3,j) - dum(3)
        end if
      end do
      call prompt (' CofG now at (0,0,0)')
c
      return
      end
c
c
c
      subroutine xyz_inertia ()
c
      include 'moleman2.incl'
c
      real a(3,3),d(3),v(3,3),rtx(12)
      real q11,q22,q33,q12,q13,q23,rms
      real det3
c
      integer i,nok,j,k
c
code ...
c
      call xyz_origin ()
c
      q11=0.0
      q22=0.0
      q33=0.0
      q12=0.0
      q13=0.0
      q23=0.0
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          q11=q11+xyz(1,i)*xyz(1,i)
          q22=q22+xyz(2,i)*xyz(2,i)
          q33=q33+xyz(3,i)*xyz(3,i)
          q12=q12+xyz(1,i)*xyz(2,i)
          q13=q13+xyz(1,i)*xyz(3,i)
          q23=q23+xyz(2,i)*xyz(3,i)
        end if
      end do
      if (nok .lt. 3) then
        call errcon ('Fewer than 3 atoms')
        return
      end if
      a(1,1)=q11
      a(1,2)=q12
      a(2,1)=q12
      a(1,3)=q13
      a(3,1)=q13
      a(2,2)=q22
      a(2,3)=q23
      a(3,2)=q23
      a(3,3)=q33
c
      call jacobi (a,3,3,d,v,i,j)
      if (j .ne. 0) then
        call errcon ('Fatal error in eigenvalue determination')
        return
      end if
c
      call eigsrt (d,v,3,3)
c
      k = 0
 1234 continue
      do i=1,3
        write (*,6000) i,d(i),(v(j,i),j=1,3)
      end do
 6000 format (' Eigen value ',i1,' = ',f12.1,' Vector : ',3f10.6)
c
      rms = det3(v)
      call fvalut (' Determinant :',1,rms)
c
      if (rms .lt. 0.0) then
        if (k.eq.0) then
          call errcon (
     +      'Negative determinant; change hand of inertia axes')
          do i=1,3
            do j=1,3
              v(i,j) = -v(i,j)
            end do
          end do
          k = 1
          goto 1234
        else
          call errcon ('Cannot fix negative determinant !!!???')
          return
        end if
      end if
c
c ... get the operator to align the inertia axes with the x, y and z axes
c     (this is the inverse (=transpose) of the matrix of eigenvectors V)
c
      do i=1,3
        do j=1,3
          a(i,j) = 0.0
        end do
        a (i,i) = 1.0
      end do
c
      call lsqgjk (a,v,3,rms,rtx,j)
      if (j .ne. 0) then
        call errcon ('Fatal error in LSQ superpositioning')
        return
      end if
c
      call anancs (1,rtx,.true.,j)
ccc      call fvalut (' RMSD :',1,rms)
c
      do i=1,natoms
        if (select(i)) then
          call vecrtv (xyz(1,i),xyz(1,i),1,rtx(1),rtx(10))
        end if
      end do
c
      return
      end
c
c
c
      subroutine ono_dicts (what,which,file,number,iunit)
c
      include 'moleman2.incl'
c
      integer length,ires,ierr,iopt,i,j,napr,iunit,leng1
c
      character myfile*256
      character line*(6*maxapr),tornam(maxtor)*6
      character*(*) what,which,file,number
c
code ...
c
      call upcase (which)
c
ccc      call remspa (which)
c
      if (what(1:2) .eq. 'RS') then
        what = 'RSR_dict'
        iopt = 2
      else if (what(1:2) .eq. 'FI') then
        what = 'RS_fit'
        iopt = 1
      else if (what(1:2) .eq. 'TO') then
        what = 'Torsion'
        iopt = 4
      else if (what(1:2) .eq. 'CO') then
        what = 'Connect'
        iopt = 3
      end if
c
      if (length(number) .gt. 0) then
        call str2i (number,ires,ierr)
        if (ierr .ne. 0) return
        do i=1,nres
          j = atmptr(1,i)
          if (iresid(j) .eq. ires) then
            if (resnam(j) .eq. which(1:3)) then
              ires = i
              goto 10
            else
              call errcon ('Selected residue of wrong type')
              call textut (' Type requested :',which)
              call prompt (' Found :')
              call print_res (i,0)
              return
            end if
          end if
        end do
        call errcon ('Selected residue number not found')
        return
      else
        do i=1,nres
          j = atmptr(1,i)
          if (resnam(j) .eq. which(1:3)) then
            ires = i
            goto 10
          end if
        end do
        call errcon ('Residue type not found')
        call textut (' Type requested :',which)
        return
      end if
c
   10 continue
      call textut (' ONO file type :',what)
      call textut (' Residue type  :',which)
      call prompt (' Residue used  :')
      call print_res (ires,0)
c
      napr = atmptr(2,ires) - atmptr(1,ires) + 1
      call jvalut (' Atoms in residue :',1,napr)
      if (napr .gt. maxapr) then
        call errcon ('Too many atoms in this residue')
        call jvalut (' Maximum  :',1,maxapr)
        return
      end if
c
      j = 0
      line = ' '
      do i=atmptr(1,ires),atmptr(2,ires)
        call print_atom (i)
        j = j + 1
        ibuf (j) = j
        rbuf ( (j-1)*3 + 1 ) = xyz (1,i)
        rbuf ( (j-1)*3 + 2 ) = xyz (2,i)
        rbuf ( (j-1)*3 + 3 ) = xyz (3,i)
        line (1+5*(j-1):4+5*(j-1)) = atmnam(i)
      end do
c
      if (length(file) .le. 0) then
        if (iopt .eq. 1 .or. iopt .eq. 2) then
          myfile = which(1:3)//'_'//
     +             what(1:leng1(what))//'.odb'
        else
          myfile = which(1:3)//'_'//
     +             what(1:leng1(what))//'.dat'
        end if
        call locase (myfile)
        call remspa (myfile)
      else
        myfile = file
        call remspa (myfile)
      end if
      call textut (' Output file   :',myfile)
c
      call xopxua (iunit,myfile,linter,ierr)
      if (ierr .ne. 0) return
c
c ... use our buffer arrays:
c
c (1) RBUF - 1 -> 3*NAPR = OTXYZ
c            3*NAPR + 1 -> 3*NAPR + NAPR*NAPR = DISMAT
c            3*NAPR + NAPR*NAPR + 1 -> .. + MAXTOR = VALTOR
c
c (2) IBUF - 1 -> 4*MAXTOR = DEFTOR
c
c (3) LBUF - 1 -> NAPR*NAPR = ISBOND
c            NAPR*NAPR + 1 -> 2*NAPR*NAPR = OKBOND
c            2*NAPR*NAPR + 1 -> 2*NAPR*NAPR + MAXTOR = AFFTOR

      call odicts (which(1:3),iopt,iunit,napr,rbuf(1),line,mxbond,
     +  maxtor,rbuf(3*napr+1),rbuf(3*napr+(napr*napr)+1),ibuf(1),
     +  lbuf(1),lbuf(napr*napr+1),lbuf(2*napr*napr+1),tornam,
     +  atmnam(atmptr(1,ires)),tortol)
c
      return
      end
c
c
c
      subroutine ono_disulfide (omol,file,myhow,objnam,iunit)
c
      include 'moleman2.incl'
c
      real strad,dd,dist
c
      integer length,ierr,i,j,iunit,ncys,nss,ii,jj,leng1
c
      character*(*) file,omol,myhow,objnam
      character s1name*6,s2name*6,line*256,how*1
c
code ...
c
      strad = 0.2
c
      call upcase (omol)
      call remspa (omol)
      if (length(omol) .lt. 1) then
        call errcon ('No O molecule name supplied')
        return
      end if
      call textut (' O molecule name :',omol)
c
      if (length(file) .lt. 1) then
        call errcon ('No ODL file name supplied')
        return
      end if
      call textut (' ODL file name :',file)
c
      how = myhow
      call upcase (how)
      if (how .ne. 'L') then
        how = 'S'
        call prompt (' Draw as Sticks')
      else
        call prompt (' Draw as Lines')
      end if
c
      if (length(objnam) .lt. 1) then
        call errcon ('No ODL object name supplied')
        return
      end if
      call textut (' ODL object name :',objnam)
c
c ... find disulfides (if any)
c
      call prompt ('0Looking for disulfides ...')
      call find_type (natoms,resnam,atmnam,'CYS',' SG ',
     +  ncys,ibuf,.true.)
c
 6040 format (' Disulfide # ',i4,1x,i5,1x,a4,' <-> ',i5,1x,a4,
     +  ' @ ',f6.2,' A')
 6600 format (20(a,1x))
 6610 format (7(a8,1x),f8.2)
c
      call jvalut (' Nr of CYS SG atoms :',1,ncys)
      if (ncys .lt. 2) return
c
      if (ncys .gt. 1) then
        call fvalut (' Max SG-SG distance for link :',1,mxcyss)
        nss = 0
        do i=1,ncys-1
          ii = ibuf(i)
          write (s1name,'(a1,i4,a1)') achain(ii),iresid(ii),
     +      insert(ii)
          call remspa (s1name)
          do j=i+1,ncys
            jj = ibuf(j)
            dd = dist (ii,jj,xyz)
            if (dd .le. mxcyss) then
              nss = nss + 1
              write (s2name,'(a1,i4,a1)') achain(jj),iresid(jj),
     +          insert(jj)
              call remspa (s2name)
              write (*,6040) nss,iresid(ii),inote(ii)(7:10),
     +          iresid(jj),inote(jj)(7:10),dd
              if (nss .eq. 1) then
                call xopxua (iunit,file,linter,ierr)
                if (ierr .ne. 0) return
                write (iunit,6600) 'begin',objnam
                write (iunit,6600) '  colour green'
              end if
c
              if (how .eq. 'S') then
                write (line,6610) 'stick',omol,s1name,atmnam(ii),
     +            omol,s2name,atmnam(jj),strad
                call pretty (line)
                call locase (line)
                write (iunit,6600) ('  '//line(1:leng1(line)))
              else
                write (line,6600) 'move_atom',omol,s1name,atmnam(ii)
                call pretty (line)
                call locase (line)
                write (iunit,6600) ('  '//line(1:leng1(line)))
                write (line,6600) 'line_atom',omol,s2name,atmnam(jj)
                call pretty (line)
                call locase (line)
                write (iunit,6600) ('  '//line(1:leng1(line)))
              end if
            end if
          end do
        end do
        call jvalut (' Nr of disulfides :',1,nss)
      end if
c
      if (nss .gt. 0) then
        write (iunit,6600) 'end_object'
        close (iunit)
        call prompt (' ODL file written')
      end if
c
      return
      end
c
c
c
      subroutine ono_water_fit (file,f9)
c
      include 'moleman2.incl'
c
      integer ierr,i,j,f9,ns,leng1
c
      character*(*) file
      character line*256,s1name*6,s2name*6
c
code ...
c
      call textut (' File name :',file)
c
      call xopxua (f9,file,linter,ierr)
      if (ierr .ne. 0) return
c
 5110 format (20(a,1x))
c
      call stamp (line)
      write (f9,5110) '!',line(1:leng1(line))
c
      write (f9,5110) '!'
      write (f9,5110) 'bell Message Set up for water fitting ...'
      write (f9,5110) '!'
      write (f9,5110) 'symbol mymol # Molecule name ? #'
      write (f9,5110) 'mol $mymol'
      write (f9,5110) '!'
      write (f9,5110) 'symbol mymap # Map file name ? #'
      write (f9,5110) 'map_file $mymap'
      write (f9,5110) 'rsr_map $mymap'
      write (f9,5110) '!'
      write (f9,5110) 'rsr_setup'
      write (f9,5110) 'yes'
      write (f9,5110) 'no'
      write (f9,5110) 'conv'
      write (f9,5110) 'RFAC'
      write (f9,5110) 'yes'
      write (f9,5110) ';'
      write (f9,5110) ';'
      write (f9,5110) '20.0'
      write (f9,5110) '3.5'
      write (f9,5110) ';'
      write (f9,5110) '3'
      write (f9,5110) '10.0'
      write (f9,5110) '!'
      write (f9,5110)
     +    'bell Message Fitting waters of mol $mymol',
     +    'in map $mymap ...'
c
      ns = 0
      do i=1,nres
        if (restyp(i) .eq. iwate) then
          j = atmptr(1,i)
          ns = ns + 1
          write (s1name,'(a1,i4,a1)') achain(j),iresid(j),
     +      insert(j)
          call remspa (s1name)
          if (ns .eq. 1) s2name = s1name
          write (f9,5110) '!'
          write (f9,5110) 'rs_fit',s1name,';'
          write (f9,5110) 'rsr_rigid',s1name,'; yes'
        end if
      end do
c
      call jvalut (' Nr of waters to be fitted :',1,ns)
c
      if (ns .gt. 0) then
        write (f9,5110) '!'
        write (f9,5110) 'Message Calculating new RS-fit values'
        write (f9,5110) '! copy_db ${mymol}//_residue_rsprefit',
     +    '${mymol}//_residue_rsfit'
        write (f9,5110) 'rs_fit',s2name,s1name
      end if
      write (f9,5110) '!'
      write (f9,5110) 'bell Message Done'
c
      close (f9)
      call prompt (' O macro written')
c
      return
      end
c
c
c
      subroutine ono_fix_hydro ()
c
      include 'moleman2.incl'
c
      integer i,j
c
      logical hd21,hd22,hh11,hh12,hh21,hh22
c
code ...
c
      call prompt (' Fixing hydrogen names for Asn and Arg')
c
      do i=1,nres
        if (resnam(atmptr(1,i)) .eq. 'ASN') then
          hd21 = .false.
          hd22 = .false.
          call print_res (i,1)
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' HD2') then
              if (.not. hd21) then
                atmnam (j) = 'HD21'
                hd21 = .true.
              else if (.not. hd22) then
                atmnam (j) = 'HD22'
                hd22 = .true.
              else
                call errcon ('Too many HD2* atoms ...')
                atmnam (j) = 'HD2?'
              end if
              call print_atom (j)
            end if
          end do
c
        else if (resnam(atmptr(1,i)) .eq. 'ARG') then
          hh11 = .false.
          hh12 = .false.
          hh21 = .false.
          hh22 = .false.
          call print_res (i,1)
          do j=atmptr(1,i),atmptr(2,i)
            if (atmnam(j) .eq. ' HH1') then
              if (.not. hh11) then
                atmnam (j) = 'HH11'
                hh11 = .true.
              else if (.not. hh12) then
                atmnam (j) = 'HH12'
                hh12 = .true.
              else
                call errcon ('Too many HH1* atoms ...')
                atmnam (j) = 'HH1?'
              end if
              call print_atom (j)
            else if (atmnam(j) .eq. ' HH2') then
              if (.not. hh21) then
                atmnam (j) = 'HH21'
                hh21 = .true.
              else if (.not. hh22) then
                atmnam (j) = 'HH22'
                hh22 = .true.
              else
                call errcon ('Too many HH2* atoms ...')
                atmnam (j) = 'HH2?'
              end if
              call print_atom (j)
            end if
          end do
        end if
      end do
c
      return
      end
c
c
c
      subroutine ono_oops (omol,file,iunit)
c
      include 'moleman2.incl'
c
      integer length,ierr,iunit,lo,ii,leng1
c
      character*(*) file,omol
      character s1name*6,s2name*6,line*256,ocom*80
c
code ...
c
      call upcase (omol)
      call remspa (omol)
      if (length(omol) .lt. 1) then
        call errcon ('No O molecule name supplied')
        return
      end if
c
      call textut (' File name :',file)
c
      call xopxua (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
 6000 format (20(a,1x))
c
      call stamp (line)
      write (iunit,6000) '!',line(1:leng1(line))
      write (iunit,6000) '! File:',file(1:leng1(file))
      write (iunit,6000) '!'
c
      lo = length(omol)
      write (iunit,6000) 'mol',omol(1:lo)
      write (iunit,6000) '! read stereo_chem.odb'
      write (iunit,6000) '! refi_init',omol(1:lo)
      write (iunit,6000) '! refi_gen',omol(1:lo),';'
c
      write (iunit,6000) '!'
      write (iunit,6000) 'read omac/message.odb'
c
      write (iunit,6000) '!'
      write (iunit,6000) 'read omac/lefth.odb'
      write (iunit,6000) 'yasspa',omol(1:lo),'lefth 0.3'
      write (iunit,6000) 'yasspa',omol(1:lo),'alpha 0.5'
      write (iunit,6000) 'yasspa',omol(1:lo),'beta 0.8'
c
      ocom = omol(1:lo)//'_residue_name'
      write (iunit,6000) '! write',ocom(1:leng1(ocom)),
     +  'resnam.o ;'
c
      ocom = omol(1:lo)//'_residue_type'
      write (iunit,6000) '! write',ocom(1:leng1(ocom)),
     +  'restyp.o ;'
c
      ii = atmptr(1,1)
      write (s1name,'(a1,i4,a1)') achain(ii),iresid(ii),
     +      insert(ii)
      call remspa (s1name)
c
      ii = atmptr(1,nres)
      write (s2name,'(a1,i4,a1)') achain(ii),iresid(ii),
     +      insert(ii)
      call remspa (s2name)
c
      write (iunit,6000) '!'
      write (iunit,6000) 'pep_flip',omol(1:lo),s1name,s2name
      ocom = omol(1:lo)//'_residue_pepflip'
      write (iunit,6000) 'write',ocom(1:leng1(ocom)),
     +  'pepflip.o ;'
c
      write (iunit,6000) '!'
      write (iunit,6000) 'rsc_fit',omol(1:lo),s1name,s2name
      ocom = omol(1:lo)//'_residue_rsc'
      write (iunit,6000) 'write',ocom(1:leng1(ocom)),
     +  'rsc.o ;'
c
      ocom = omol(1:lo)//'_residue_rsfit'
c
      write (iunit,6000) '!'
      write (iunit,6000) 'symbol mymap # Map file name ? #'
      write (iunit,6000) 'map_file $mymap'
      write (iunit,6000) 'rsr_map $mymap ;'
c
      write (iunit,6000) '!'
      write (iunit,6000) '! calculate RSCC values'
      write (iunit,6000) '!'
      write (iunit,6000) 'rsr_setup'
      write (iunit,6000) 'yes'
      write (iunit,6000) 'no'
      write (iunit,6000) 'conv'
      write (iunit,6000) 'RSCC'
      write (iunit,6000) 'yes'
      write (iunit,6000) ';'
      write (iunit,6000) ';'
      write (iunit,6000) '20.0'
      write (iunit,6000) '3.5'
      write (iunit,6000) '1.04 0.9'
      write (iunit,6000) '3'
      write (iunit,6000) '10.0'
      write (iunit,6000) '!'
      write (iunit,6000) 'db_kill',ocom(1:leng1(ocom))
      write (iunit,6000) 'read odat/rsfit_all.o'
      write (iunit,6000) 'rs_fit',omol(1:lo),s1name,s2name
      write (iunit,6000) 'write',ocom(1:leng1(ocom)),
     +  'rsfit_all.o ;'
      write (iunit,6000) '!'
      write (iunit,6000) '! db_kill',ocom(1:leng1(ocom))
      write (iunit,6000) '! read odat/rsfit_mc.o'
      write (iunit,6000) '! rs_fit',omol(1:lo),s1name,s2name
      write (iunit,6000) '! write',ocom(1:leng1(ocom)),
     +  'rsfit_mc.o ;'
      write (iunit,6000) '!'
      write (iunit,6000) '! db_kill',ocom(1:leng1(ocom))
      write (iunit,6000) '! read odat/rsfit_sc.o'
      write (iunit,6000) '! rs_fit',omol(1:lo),s1name,s2name
      write (iunit,6000) '! write',ocom(1:leng1(ocom)),
     +  'rsfit_sc.o ;'
c
      write (iunit,6000) '!'
      write (iunit,6000) '! calculate RSR values'
      write (iunit,6000) '!'
      write (iunit,6000) 'rsr_setup'
      write (iunit,6000) 'yes'
      write (iunit,6000) 'no'
      write (iunit,6000) 'conv'
      write (iunit,6000) 'RFAC'
      write (iunit,6000) 'yes'
      write (iunit,6000) ';'
      write (iunit,6000) ';'
      write (iunit,6000) '20.0'
      write (iunit,6000) '3.5'
      write (iunit,6000) '1.04 0.9'
      write (iunit,6000) '3'
      write (iunit,6000) '10.0'
      write (iunit,6000) '!'
      write (iunit,6000) 'db_kill',ocom(1:leng1(ocom))
      write (iunit,6000) 'read odat/rsfit_all.o'
      write (iunit,6000) 'rs_fit',omol(1:lo),s1name,s2name
      write (iunit,6000) 'write',ocom(1:leng1(ocom)),
     +  'rsrfac_all.o ;'
      write (iunit,6000) '!'
      write (iunit,6000) 'bell message Done'
      write (iunit,6000) 'print Done'
      write (iunit,6000) '!'
c
      call prompt (' O macro written')
      close (iunit)
c
      return
      end
c
c
c
      subroutine sequence_list (how)
c
      include 'moleman2.incl'
c
      integer i
c
      character how*(*)
c
code ...
c
      call upcase (how)
      call remspa (how)
c
      if (how(1:1) .ne. '1' .and. how(1:1) .ne. 'F') how(1:1) = '3'
c
 6000 format (6(1x,10a1))
 6100 format (15(1x,a3))
c
      if (how(1:1) .eq. '1') then
        call prompt ('0Sequence in 1-letter code:')
        write (*,6000) (onelc(i),i=1,nres)
      else if (how(1:1) .eq. '3') then
        call prompt ('0Sequence in 3-letter code:')
        write (*,6100) (resnam(atmptr(1,i)),i=1,nres)
      else
        do i=1,nres
          call print_res (i,0)
        end do
      end if
c
      return
      end
c
c
c
      subroutine sequence_pir (f8,file,name,title)
c
      include 'moleman2.incl'
c
      integer f8,i,ierr,leng1,now
c
      logical select_res
c
      character file*(*),name*(*),title*(*),line*128
      character selpir(maxres)*1
c
code ...
c
      call xopxua (f8,file,linter,ierr)
      if (ierr .ne. 0) return
c
      call stamp (line)
      write (f8,'(a1,1x,a)') '!',line(1:leng1(line))
      write (f8,*)
      write (f8,'(a4,a)') '>P1;',name(1:leng1(name))
      write (f8,'(a)') title(1:leng1(title))
c
      now = 0
      do i=1,nres
        if (select_res(i)) then
          now = now + 1
          selpir (now) = onelc(i)
        end if
      end do
c
      call jvalut (' Nr of selected residues  :',1,now)
      call jvalut (' Total number of residues :',1,nres)
c
      if (now .gt. 0) then
c        write (f8,'(6(10a1,1x))') (onelc(i),i=1,nres)
        write (f8,'(6(10a1,1x))') (selpir(i),i=1,now)
        write (*,'(6(10a1,1x))') (selpir(i),i=1,now)
      else
        call errcon ('No residues selected')
      end if
c
      write (f8,'(a1)') '*'
c
      close (f8)
      call prompt (' PIR file written ...')
c
      return
      end
c
c
c
      subroutine sequence_glyco ()
c
      include 'moleman2.incl'
c
      integer i,nok
c
code ...
c
      call prompt (' Looking for potential N-glycosylation sites')
      call prompt (' Consensus: Asn - not Pro - Ser/Thr')
c
      nok = 0
      do i=1,nres-2
        if (resnam(atmptr(1,i)) .eq. 'ASN') then
          if (iresid(atmptr(1,i+2)) .eq. (2+iresid(atmptr(1,i)))) then
            if (resnam(atmptr(1,i+1)) .ne. 'PRO') then
              if (resnam(atmptr(1,i+2)) .eq. 'SER' .or.
     +            resnam(atmptr(1,i+2)) .eq. 'THR') then
                call print_res (i,1)
                call print_res (i+1,0)
                call print_res (i+2,0)
                nok = nok + 1
              end if
            end if
          end if
        end if
      end do
c
      call jvalut (' Potential sites found :',1,nok)
c
      return
      end
c
c
c
      subroutine sequence_motif (motif)
c
      include 'moleman2.incl'
c
      integer i,nok,length,lm,j,k
c
      character motif*(*)
c
code ...
c
      call upcase (motif)
      call remspa (motif)
      call textut (' Look for sequence motif :',motif)
      lm = length(motif)
      if (lm .lt. 1) return
c
      nok = 0
      do i=1,nres-lm+1
        do j=1,lm
          k = i + j - 1
          if (.not. ( motif(j:j) .eq. '?' .or.
     +                motif(j:j) .eq. onelc(k) ) ) goto 1000
cccc          if (j.eq.1) write (*,*) (onelc(k),k=i,i+lm-1)
        end do
        if ( (iresid(atmptr(1,i))+lm-1) .ne.
     +       iresid(atmptr(1,i+lm-1)) ) goto 1000
c
        call print_res (i,1)
        do j=i+1,i+lm-1
          call print_res (j,0)
        end do
        nok = nok + 1
 1000   continue
      end do
      call jvalut (' Nr of hits :',1,nok)
c
      return
      end
c
c
c
      subroutine sequence_count ()
c
      include 'moleman2.incl'
c
      integer i,nok,leng1
c
      logical select_res
c
code ...
c
      nok = 0
      do i=1,nres
        lbuf (i) = select_res (i)
        if (lbuf(i)) nok = nok + 1
      end do
c
      call jvalut (' Nr of selected residues :',1,nok)
      if (nok .lt. 1) return
c
      do i=1,nmrtyp + 1
        ibuf (i) = 0
      end do
c
      do i=1,nres
        if (lbuf(i)) then
ccc          call print_res (i,0)
ccc          print *,' >>> ',typptr(i)
          if (typptr(i) .gt. 0) then
            ibuf (typptr(i)) = ibuf (typptr(i)) + 1
          else
            call errcon ('Residue not in library :')
            call print_res (i,0)
            ibuf (nmrtyp + 1) = ibuf (nmrtyp + 1) + 1
          end if
        end if
      end do
c
      if (ibuf(nmrtyp+1) .gt .0) write (*,*)
      do i=1,nmrtyp
        if (ibuf(i) .gt. 0) then
          write (*,6000) ibuf(i),libolc(i),lrname(i),
     +      lrdesc(i)(1:leng1(lrdesc(i)))
        end if
      end do
      call jvalut (' Residues not in library :',1,ibuf(nmrtyp+1))
c
 6000 format (1x,i6,' |',a1,'| ',a3,' = ',a)
c
      return
      end
c
c
c
      subroutine sequence_ex_280 ()
c
      include 'moleman2.incl'
c
      real e
c
      integer i,j,nok,ny,nw,nc
c
      logical select_res
c
code ...
c
      call prompt (' Reference: Gill & Von Hippel,')
      call prompt (' Anal. Biochem. 182: 319-326 (1989).')
      call prompt (' e280 = (1280*TYR + 5690*TRP + 120*CYS)')
      call prompt (' accuracy roughly 5%; units M-1 cm-1')
c
      nok = 0
      ny = 0
      nw = 0
      nc = 0
      do i=1,nres
        lbuf (i) = select_res (i)
        if (lbuf(i)) then
          nok = nok + 1
          j=atmptr(1,i)
          if (resnam(j) .eq. 'TYR') then
            ny = ny + 1
          else if (resnam(j) .eq. 'TRP') then
            nw = nw + 1
          else if (resnam(j) .eq. 'CYS') then
            nc = nc + 1
          end if
        end if
      end do
c
      call jvalut (' Nr of selected residues :',1,nok)
      if (nok .lt. 1) return
      call ivalut (' Nr of TYR      :',1,ny)
      call ivalut (' Nr of TRP      :',1,nw)
      call ivalut (' Nr of CYS      :',1,nc)
c
      e = 1280.0*ny + 5690.0*nw + 120.0*nc
      call rvalut (' e280 (M-1 cm-1) ~',1,e)
c
      return
      end
c
c
c
      subroutine ono_ls_angles (i1,j1)
c
      include 'moleman2.incl'
c
      real xlen(maxpln)
      real cosa,ang,x
c
      integer i1,j1,ilo,jlo,ihi,jhi,i,j
c
code ...
c
      if (i1 .le. 0 .or. i1 .gt. maxpln) then
        ilo = 1
        ihi = maxpln
      else
        ilo = i1
        ihi = i1
      end if
c
      if (j1 .le. 0 .or. j1 .gt. maxpln) then
        jlo = 1
        jhi = maxpln
      else
        jlo = j1
        jhi = j1
      end if
c
      do i=1,maxpln
        xlen (i) = vplane(1,i)*vplane(1,i) +
     +             vplane(2,i)*vplane(2,i) +
     +             vplane(3,i)*vplane(3,i)
        xlen (i) = sqrt(xlen(i))
      end do
c
      do i=ilo,ihi
        if (.not. lplane(i)) goto 200
        do j=jlo,jhi
          if (.not. lplane(j)) goto 100
          if (i .eq. j) goto 100
c
          x = vplane(1,i)*vplane(1,j) +
     +        vplane(2,i)*vplane(2,j) +
     +        vplane(3,i)*vplane(3,j)
          cosa = x / sqrt(xlen(i)*xlen(j))
          ang = rtodeg * acos (cosa)
          write (*,6000) i,nplane(i),j,nplane(j),cosa,ang,180.0-ang
c
  100     continue
        end do
  200   continue
      end do
c
 6000 format (/' Plane # ',i3,' [',a,'] and plane # ',i3,' [',a,']'/
     +         ' COS (angle) = ',f8.5,' -> Angle = ',f8.2,' = ',f8.2)
c
      return
      end
c
c
c
      subroutine ono_cell (iunit,odlfil,mode,col)
c
      include 'moleman2.incl'
c
      real a(3,3),x(3),flo(3),fhi(3)
c
      integer iunit,ierr,leng1,i
c
      logical lsolid
c
      character odlfil*(*),mode*(*),col*(*)
c
code ...
c
c ... generate ODL file to draw the cell
c
      call xopxua (iunit,odlfil,linter,ierr)
      if (ierr .ne. 0) return
c
      lsolid = (mode(1:1) .ne. 'L')
c
c --- A = frac to ang, B=Ang to frac
c
      call orthog (cell, a, 0)
c
      do i=1,3
        flo(i) = 0.0
        fhi (i) = 1.0
      end do
c
      call fvalut (' Origin (fract.) :',3,flo)
      call fvalut (' Top    (fract.) :',3,fhi)
c
      call mulmtx (a, flo, x, 3, 3, 1)
      call fvalut (' Origin (cart.)  :',3,x)
      call mulmtx (a, fhi, x, 3, 3, 1)
      call fvalut (' Top    (cart.)  :',3,x)
c
c ... now write the ODL object
c
      write (iunit,6000) 'begin celly'
      if (lsolid) write (iunit,6000) 'mode solid'
      write (iunit,6010) 'colour ',col(1:leng1(col))
c
      if (lsolid) then
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'poly 4'
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) ' ',(x(i),i=1,3)
c
        write (iunit,6000) 'mode line'
c
      else
c
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
c
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
c
        call convec (flo(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (flo(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),flo(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (fhi(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),flo(3),x,a)
        write (iunit,6000) 'm ',(x(i),i=1,3)
        call convec (flo(1),fhi(2),fhi(3),x,a)
        write (iunit,6000) 'l ',(x(i),i=1,3)
c
      end if
c
      write (iunit,6000) 'end_object'
c
 6000 format (a,3f10.3)
 6010 format (9a)
c
      close (iunit)
      write (*,*) 'ODL file written'
c
      return
      end
c
c
c
      subroutine ono_rings (iunit,odlfil,obj,col)
c
      include 'moleman2.incl'
c
      real rngxyz(3,10)
c
      integer i,nok,j,k,iunit,ierr,length,leng1,nhits
      integer kres,k1,k2,l,m
c
      character odlfil*(*),obj*(*),col*(*),line*128
c
code ...
c
      if (nlring .lt. 1) then
        call errcon ('Library file contained no RNG definitions !')
        return
      end if
c
c ... generate ODL file to draw the rings
c
      call xopxua (iunit,odlfil,linter,ierr)
      if (ierr .ne. 0) return
c
 7000 format (20(a,1x))
 7010 format (2x,3(f10.2,1x))
 7020 format (a,i10)
c
      call stamp (line)
      write (iunit,7000,err=8000) '!',line(1:leng1(line))
      write (iunit,7000,err=8000) 'begin',obj(1:length(obj))
      write (iunit,7000,err=8000) '  mode solid'
      write (iunit,7000,err=8000) '  colour',col(1:length(col))
c
      nok = 0
c
      do i=1,nlring
        call find_type (natoms,resnam,atmnam,rngres(i),
     +                  rngatm(1,i),nhits,ibuf,.false.)
        if (nhits .gt. 0) then
          do j=1,nhits
            k = ibuf(j)
            kres = resptr(k)
            k1 = atmptr(1,kres)
            k2 = atmptr(2,kres)
            do l=k1,k2
              if (select(l)) goto 100
            end do
            goto 200
c
c ... at least one atom in residue selected
c
  100       continue
            rngxyz (1,1) = xyz (1,k)
            rngxyz (2,1) = xyz (2,k)
            rngxyz (3,1) = xyz (3,k)
            do m=2,rngnat(i)
              do l=k1,k2
                if (rngatm(m,i) .eq. atmnam(l)) then
                  rngxyz (1,m) = xyz (1,l)
                  rngxyz (2,m) = xyz (2,l)
                  rngxyz (3,m) = xyz (3,l)
                  goto 110
                end if
              end do
              goto 200
  110         continue
            end do
c
            write (iunit,7020,err=8000) '  poly ',rngnat(i)
            do m=1,rngnat(i)
              write (iunit,7010,err=8000) (rngxyz(l,m),l=1,3)
            end do
            nok = nok + 1
c
  200       continue
          end do
        end if
      end do
c
      write (iunit,7000,err=8000) '  mode line'
      write (iunit,7000,err=8000) 'end_object'
c
      call jvalut (' Nr of rings found :',1,nok)
      close (iunit)
      return
c
 8000 continue
      call errcon ('While writing ODL file')
      close (iunit)
c
      return
      end
c
c
c
      subroutine xyz_ls_plane (iunit,odlfil,obj,col,border,
     +                         pnum,pnam,mymode)
c
      include 'moleman2.incl'
c
      real a(3,3),c(3),d(3),e(3),v(3,3),h(3,3),rtx(12),circum(3)
      real drawem(5),q(3),inout(3),axvec(3)
      real q11,q22,q33,q12,q13,q23,rms,x1,x2,x3
      real xmin,xmax,ymin,ymax,zmin,zmax,border,xdum
      real det3,distce
c
      integer i,nok,j,k,iunit,ierr,length,leng1,pnum
c
      logical lodl,liax
c
      character odlfil*(*),obj*(*),col*(*),line*128,pnam*(*)
      character mode*1,mymode*(*)
c
code ...
c
      lodl = (iunit .gt. 0 .and. length(odlfil).gt.0)
      mode = mymode
      call upcase (mode)
      liax = (mode .eq. 'I')
c
      nok = 0
      do i=1,natoms
        if (select(i)) then
          nok = nok + 1
          ibuf (nok) = i
          j = 3*(nok-1)
          rbuf (j+1) = xyz(1,i)
          rbuf (j+2) = xyz(2,i)
          rbuf (j+3) = xyz(3,i)
        end if
      end do
      call jvalut (' Nr of selected atoms :',1,nok)
      if (nok .lt. 3) then
        call errcon ('Fewer than 3 atoms')
        return
      end if
c
      do i=1,3
        c(i) = 0.0
      end do
      xmin = rbuf(1)
      xmax = rbuf(1)
      ymin = rbuf(2)
      ymax = rbuf(2)
      zmin = rbuf(3)
      zmax = rbuf(3)
      do i=1,nok
        j = 3*(i-1)
        c(1) = c(1) + rbuf(j+1)
        c(2) = c(2) + rbuf(j+2)
        c(3) = c(3) + rbuf(j+3)
        xmin = min (xmin, rbuf(j+1))
        xmax = max (xmax, rbuf(j+1))
        ymin = min (ymin, rbuf(j+2))
        ymax = max (ymax, rbuf(j+2))
        zmin = min (zmin, rbuf(j+3))
        zmax = max (zmax, rbuf(j+3))
      end do
c
      do i=1,3
        c(i) = c(i)/float(nok)
      end do
      call fvalut (' Centre of Gravity :',3,c)
      do i=1,nok
        j = 3*(i-1)
        rbuf (j+1) = rbuf(j+1) - c(1)
        rbuf (j+2) = rbuf(j+2) - c(2)
        rbuf (j+3) = rbuf(j+3) - c(3)
      end do
c
      q11=0.0
      q22=0.0
      q33=0.0
      q12=0.0
      q13=0.0
      q23=0.0
      do i=1,nok
        j = 3*(i-1)
        x1 = rbuf(j+1)
        x2 = rbuf(j+2)
        x3 = rbuf(j+3)
        q11=q11+x1*x1
        q22=q22+x2*x2
        q33=q33+x3*x3
        q12=q12+x1*x2
        q13=q13+x1*x3
        q23=q23+x2*x3
      end do
c
      a(1,1)=q11
      a(1,2)=q12
      a(2,1)=q12
      a(1,3)=q13
      a(3,1)=q13
      a(2,2)=q22
      a(2,3)=q23
      a(3,2)=q23
      a(3,3)=q33
c
      call jacobi (a,3,3,d,v,i,j)
      if (j .ne. 0) then
        call errcon ('Fatal error in eigenvalue determination')
        return
      end if
c
      call eigsrt (d,v,3,3)
c
      k = 0
 1234 continue
      do i=1,3
        write (*,6000) i,d(i),(v(j,i),j=1,3)
      end do
 6000 format (' Eigen value ',i1,' = ',f12.1,' Vector : ',3f10.6)
c
      rms = det3(v)
      call fvalut (' Determinant :',1,rms)
c
      if (rms .lt. 0.0) then
        if (k.eq.0) then
          call errcon (
     +      'Negative determinant; change hand of inertia axes')
          do i=1,3
            do j=1,3
              v(i,j) = -v(i,j)
            end do
          end do
          k = 1
          goto 1234
        else
          call errcon ('Cannot fix negative determinant !!!???')
          return
        end if
      end if
c
c ... print plane
c
      call prompt (' Eigenvector #3 defines the least-squares plane')
      x1 = v(1,3)*c(1) + v(2,3)*c(2) + v(3,3)*c(3)
      write (*,6100) v(1,3),v(2,3),v(3,3),x1
 6100 format (' Equation: ',f10.6,' X + ',f10.6,' Y + ',f10.6,
     +  ' Z = ',f10.6)
c
c ... darw inertia axes
c
      if (liax) then
        call xopxua (iunit,odlfil,linter,ierr)
        if (ierr .ne. 0) return
c
        if (length(obj) .lt. 1) obj = '_iaxes'
        if (length(col) .lt. 1) col = 'cyan'
c
        call textut (' Creating ODL file   :',odlfil)
        call textut (' Object name         :',obj)
c
        if (border .gt. 0.001) then
ccc          call fvalut (' Axis length (A)     :',1,border)
          do i=1,3
            d(i) = border
          end do
        else if (border .lt. -0.001) then
          border = -1.0 * border
ccc          call fvalut (' Longest axis (A)    :',1,border)
          xdum = 1.0 / d(1)
          do i=1,3
            d(i) = d(i) * border * xdum
          end do
        else
          d(1) = 3.0
          d(2) = 2.0
          d(3) = 1.0
        end if
        call fvalut (' Axis lengths (A)    :',3,d)
c
        call stamp (line)
        write (iunit,7000,err=8000) '!',line(1:leng1(line))
        write (iunit,7000,err=8000) 'begin',obj(1:length(obj))
        write (iunit,7000,err=8000) '  colour',col(1:length(col))
c
        do i=1,3
          axvec (1) = c(1) + d(i) * v(1,i)
          axvec (2) = c(2) + d(i) * v(2,i)
          axvec (3) = c(3) + d(i) * v(3,i)
          write (iunit,7020,err=8000) '  m ',(c(j),j=1,3)
          write (iunit,7020,err=8000) '  l ',(axvec(j),j=1,3)
        end do
c
        write (iunit,7000,err=8000) 'end_object'
c
        call prompt (' ODL file written')
        close (iunit)
        return
c
      end if
c
      if (lodl .and. pnum .gt. 0) then
        pnum = min (pnum,maxpln)
        lplane (pnum) = .true.
        if (pnam .ne. ' ') nplane (pnum) = pnam
        vplane (1,pnum) = v(1,3)
        vplane (2,pnum) = v(2,3)
        vplane (3,pnum) = v(3,3)
        line = ' Stored plane "'//nplane(pnum)//'" as number :'
        call ivalut (line,1,pnum)
      end if
c
c ... get the operator to align the inertia axes with the x, y and z axes
c     (this is the inverse (=transpose) of the matrix of eigenvectors V)
c
      do i=1,3
        do j=1,3
          a(i,j) = 0.0
        end do
        a (i,i) = 1.0
      end do
c
      call lsqgjk (a,v,3,rms,rtx,j)
      if (j .ne. 0) then
        call errcon ('Fatal error in LSQ superpositioning')
        return
      end if
c
ccc      call anancs (1,rtx,.true.,j)
ccc      call fvalut (' RMSD :',1,rms)
c
      do i=1,nok
        j = 3*(i-1)
        call vecrtv (rbuf(j+1),rbuf(j+1),1,rtx(1),rtx(10))
      end do
c
      if (lodl) write (*,*)
      x3 = 0.0
      do i=1,nok
        j = 3*(i-1)
        k = ibuf(i)
        if (lodl) then
          write (*,6200) atomnr(k),atmnam(k),altloc(k),
     +      resnam(k),achain(k),iresid(k),insert(k),rbuf(j+3)
        end if
        x3 = x3 + rbuf(j+3)**2
      end do
 6200 format (' ATOM ',i5,1x,a4,a1,a3,1x,a1,i4,a1,' Distance ',f8.3)
      x3 = sqrt ( x3 / float(nok) )
      if (lodl) write (*,*)
      call fvalut (' RMSD to plane :',1,x3)
c
      if (.not. lodl) return
c
ccc      if (v(3,3) .eq. 0.0) then
ccc        call errcon ('Z coeff. of plane is zero; sorry !')
ccc      end if
c
c ... generate ODL file to draw the plane
c
      call xopxua (iunit,odlfil,linter,ierr)
      if (ierr .ne. 0) return
c
      if (length(obj) .lt. 1) obj = '_plane'
      if (length(col) .lt. 1) col = 'cyan'
c
 7000 format (20(a,1x))
 7010 format (2x,3(f10.2,1x))
 7020 format (a,1x,3(f10.2,1x))
c
ccc      border = 1.0
c
      call textut (' Creating ODL file   :',odlfil)
      call textut (' Object name         :',obj)
      call textut (' Object colour       :',col)
      call fvalut (' Border around atoms :',1,border)
c
      xmin = xmin - border
      xmax = xmax + border
      ymin = ymin - border
      ymax = ymax + border
      zmin = zmin - border
      zmax = zmax + border
c
      x1 = v(1,3)*c(1) + v(2,3)*c(2) + v(3,3)*c(3)
c
      call stamp (line)
      write (iunit,7000,err=8000) '!',line(1:leng1(line))
      write (iunit,7000,err=8000) 'begin',obj(1:length(obj))
      write (iunit,7000,err=8000) '  mode solid'
      write (iunit,7000,err=8000) '  colour',col(1:length(col))
      write (iunit,7000,err=8000) '  poly 4'
c
c ... find out which way gives the least elongated plane
c
      if (v(3,3) .eq. 0) then
        circum (1) = 1.0E30
      else
        c(1) = xmin
        c(2) = ymin
        c(3) = (x1 - v(1,3)*xmin - v(2,3)*ymin) / v(3,3)
        d(1) = xmin
        d(2) = ymax
        d(3) = (x1 - v(1,3)*xmin - v(2,3)*ymax) / v(3,3)
        e(1) = xmax
        e(2) = ymax
        e(3) = (x1 - v(1,3)*xmax - v(2,3)*ymax) / v(3,3)
        circum (1) = 2.0 * (distce(c,d) + distce(d,e))
        call gkjifz (c,d,e,v,h(1,1),inout(1))
      end if
c
      if (v(2,3) .eq. 0) then
        circum (2) = 1.0E30
      else
        c(1) = xmin
        c(3) = zmin
        c(2) = (x1 - v(1,3)*xmin - v(3,3)*zmin) / v(2,3)
        d(1) = xmin
        d(3) = zmax
        d(2) = (x1 - v(1,3)*xmin - v(3,3)*zmax) / v(2,3)
        e(1) = xmax
        e(3) = zmax
        e(2) = (x1 - v(1,3)*xmax - v(3,3)*zmax) / v(2,3)
        circum (2) = 2.0 * (distce(c,d) + distce(d,e))
        call gkjifz (c,d,e,v,h(1,2),inout(2))
      end if
c
      if (v(1,3) .eq. 0) then
        circum (3) = 1.0E30
      else
        c(2) = ymin
        c(3) = zmin
        c(1) = (x1 - v(2,3)*ymin - v(3,3)*zmin) / v(1,3)
        d(2) = ymin
        d(3) = zmax
        d(1) = (x1 - v(2,3)*ymin - v(3,3)*zmax) / v(1,3)
        e(2) = ymax
        e(3) = zmax
        e(1) = (x1 - v(2,3)*ymax - v(3,3)*zmax) / v(1,3)
        circum (3) = 2.0 * (distce(c,d) + distce(d,e))
        call gkjifz (c,d,e,v,h(1,3),inout(3))
      end if
c
      q(1) = -0.01 * v(1,3)
      q(2) = -0.01 * v(2,3)
      q(3) = -0.01 * v(3,3)
      q11 = sqrt (q(1)*q(1) + q(2)*q(2) + q(3)*q(3))
      q(1) = 0.5 * q(1)/q11
      q(2) = 0.5 * q(2)/q11
      q(3) = 0.5 * q(3)/q11
c
      if ( circum(1) .le. circum(2) .and.
     +     circum(1) .le. circum(3) ) then
c
ccc        print *,' 11111111111111111111111111111111111 ',inout(1)
        do i=1,3
          q(i) = -h(i,1)
        end do
        q11 = sqrt (q(1)*q(1) + q(2)*q(2) + q(3)*q(3))
        q(1) = 0.5 * q(1)/q11
        q(2) = 0.5 * q(2)/q11
        q(3) = 0.5 * q(3)/q11
c
        x2 = (x1 - v(1,3)*xmin - v(2,3)*ymin) / v(3,3)
        drawem (1) = x2
        write (iunit,7010,err=8000) xmin,ymin,drawem(1)
        x2 = (x1 - v(1,3)*xmin - v(2,3)*ymax) / v(3,3)
        drawem (2) = x2
        write (iunit,7010,err=8000) xmin,ymax,drawem(2)
        x2 = (x1 - v(1,3)*xmax - v(2,3)*ymax) / v(3,3)
        drawem (3) = x2
        write (iunit,7010,err=8000) xmax,ymax,drawem(3)
        x2 = (x1 - v(1,3)*xmax - v(2,3)*ymin) / v(3,3)
        drawem (4) = x2
        write (iunit,7010,err=8000) xmax,ymin,drawem(4)
ccc        x2 = (x1 - v(1,3)*xmin - v(2,3)*ymin) / v(3,3)
ccc        drawem (5) = x2
ccc        write (iunit,7010,err=8000) xmin,ymin,drawem(5)
c
c        write (iunit,7000,err=8000) '  colour',col(1:length(col))
c        write (iunit,7000,err=8000) '  poly 5'
c        write (iunit,7010,err=8000)
c     +    xmin+q(1),ymin+q(2),drawem(5)+q(3)
c        write (iunit,7010,err=8000)
c     +    xmax+q(1),ymin+q(2),drawem(4)+q(3)
c        write (iunit,7010,err=8000)
c     +    xmax+q(1),ymax+q(2),drawem(3)+q(3)
c        write (iunit,7010,err=8000)
c     +    xmin+q(1),ymax+q(2),drawem(2)+q(3)
c        write (iunit,7010,err=8000)
c     +    xmin+q(1),ymin+q(2),drawem(1)+q(3)
c
      else if ( circum(2) .le. circum(1) .and.
     +          circum(2) .le. circum(3) ) then
c
ccc        print *,' 22222222222222222222222222222222222 ',inout(2)
        do i=1,3
          q(i) = -h(i,2)
        end do
        q11 = sqrt (q(1)*q(1) + q(2)*q(2) + q(3)*q(3))
        q(1) = 0.5 * q(1)/q11
        q(2) = 0.5 * q(2)/q11
        q(3) = 0.5 * q(3)/q11
c
        x2 = (x1 - v(1,3)*xmin - v(3,3)*zmin) / v(2,3)
        drawem (1) = x2
        write (iunit,7010,err=8000) xmin,drawem(1),zmin
        x2 = (x1 - v(1,3)*xmin - v(3,3)*zmax) / v(2,3)
        drawem (2) = x2
        write (iunit,7010,err=8000) xmin,drawem(2),zmax
        x2 = (x1 - v(1,3)*xmax - v(3,3)*zmax) / v(2,3)
        drawem (3) = x2
        write (iunit,7010,err=8000) xmax,drawem(3),zmax
        x2 = (x1 - v(1,3)*xmax - v(3,3)*zmin) / v(2,3)
        drawem (4) = x2
        write (iunit,7010,err=8000) xmax,drawem(4),zmin
ccc        x2 = (x1 - v(1,3)*xmin - v(3,3)*zmin) / v(2,3)
ccc        drawem (5) = x2
ccc        write (iunit,7010,err=8000) xmin,drawem(5),zmin
c
c        write (iunit,7000,err=8000) '  colour',col(1:length(col))
c        write (iunit,7000,err=8000) '  poly 5'
c        write (iunit,7010,err=8000)
c     +    xmin+q(1),drawem(5)+q(2),zmin+q(3)
c        write (iunit,7010,err=8000)
c     +    xmax+q(1),drawem(4)+q(2),zmin+q(3)
c        write (iunit,7010,err=8000)
c     +    xmax+q(1),drawem(3)+q(2),zmax+q(3)
c        write (iunit,7010,err=8000)
c     +    xmin+q(1),drawem(2)+q(2),zmax+q(3)
c        write (iunit,7010,err=8000)
c     +    xmin+q(1),drawem(1)+q(2),zmin+q(3)
c
      else
c
ccc        print *,' 3333333333333333333333333333333333333 ',inout(3)
        do i=1,3
          q(i) = -h(i,2)
        end do
        q11 = sqrt (q(1)*q(1) + q(2)*q(2) + q(3)*q(3))
        q(1) = 0.5 * q(1)/q11
        q(2) = 0.5 * q(2)/q11
        q(3) = 0.5 * q(3)/q11
c
        x2 = (x1 - v(2,3)*ymin - v(3,3)*zmin) / v(1,3)
        drawem (1) = x2
        write (iunit,7010,err=8000) drawem(1),ymin,zmin
        x2 = (x1 - v(2,3)*ymin - v(3,3)*zmax) / v(1,3)
        drawem (2) = x2
        write (iunit,7010,err=8000) drawem(2),ymin,zmax
        x2 = (x1 - v(2,3)*ymax - v(3,3)*zmax) / v(1,3)
        drawem (3) = x2
        write (iunit,7010,err=8000) drawem(3),ymax,zmax
        x2 = (x1 - v(2,3)*ymax - v(3,3)*zmin) / v(1,3)
        drawem (4) = x2
        write (iunit,7010,err=8000) drawem(4),ymax,zmin
ccc        x2 = (x1 - v(2,3)*ymin - v(3,3)*zmin) / v(1,3)
ccc        drawem (5) = x2
ccc        write (iunit,7010,err=8000) drawem(5),ymin,zmin
c
c        write (iunit,7000,err=8000) '  colour',col(1:length(col))
c        write (iunit,7000,err=8000) '  poly 5'
c        write (iunit,7010,err=8000)
c     +    drawem(5)+q(1),ymin+q(2),zmin+q(3)
c        write (iunit,7010,err=8000)
c     +    drawem(4)+q(1),ymax+q(2),zmin+q(3)
c        write (iunit,7010,err=8000)
c     +    drawem(3)+q(1),ymax+q(2),zmax+q(3)
c        write (iunit,7010,err=8000)
c     +    drawem(2)+q(1),ymin+q(2),zmax+q(3)
c        write (iunit,7010,err=8000)
c     +    drawem(1)+q(1),ymin+q(2),zmin+q(3)
c
      end if
c
      write (iunit,7000,err=8000) '  mode line'
      write (iunit,7000,err=8000) 'end_object'
c
      call prompt (' ODL file written')
      close (iunit)
      return
c
 8000 continue
      call errcon ('While writing ODL file')
      close (iunit)
c
      return
      end
c
c
c
      subroutine gkjifz (c,d,e,v,h,x)
c
      implicit none
c
      real c(3),d(3),e(3),v(3),f(3),g(3),h(3),x,y,z
c
      integer i
c
code ...
c
      do i=1,3
        f(i) = d(i) - c(i)
        g(i) = e(i) - d(i)
      end do
      h(1) = f(2)*g(3) - f(3)*g(2)
      h(2) = f(1)*g(3) - f(3)*g(1)
      h(3) = f(1)*g(2) - f(2)*g(1)
      x = 0.0
      y = 0.0
      z = 0.0
      do i=1,3
        x = x + v(i)*h(i)
        y = y + v(i)*v(i)
        z = z + h(i)*h(i)
      end do
      x = x / sqrt(y*z)
c
      return
      end
c
c
c
      subroutine auto_bones (type,coords)
c
      include 'moleman2.incl'
c
      real start(3),end(3)
c
      integer ierr,num
c
      character type*(*),coords(6)*(*)
c
code ...
c
      if (type(1:1) .ne. 'A' .and. type(1:1) .ne. 'B') then
        call errcon ('Invalid type')
        return
      end if
c
      call str2r (coords(1),start(1),ierr)
      if (ierr .ne. 0) return
      call str2r (coords(2),start(2),ierr)
      if (ierr .ne. 0) return
      call str2r (coords(3),start(3),ierr)
      if (ierr .ne. 0) return
      call str2r (coords(4),end(1),ierr)
      if (ierr .ne. 0) return
      call str2r (coords(5),end(2),ierr)
      if (ierr .ne. 0) return
      call str2r (coords(6),end(3),ierr)
      if (ierr .ne. 0) return
c
      num = -1
c
      call do_spink (start,end,type,num)
c
      return
      end
c
c
c
      subroutine auto_spink (type,cnres)
c
      include 'moleman2.incl'
c
      real start(3),end(3)
c
      integer i,ierr,num
c
      character type*(*),cnres*(*)
c
code ...
c
      if (type(1:1) .ne. 'A' .and. type(1:1) .ne. 'B') then
        call errcon ('Invalid type')
        return
      end if
c
      do i=1,3
        start(i) = 0.0
        end (i) = 0.0
      end do
c
      call str2i (cnres,num,ierr)
      if (ierr .ne. 0) return
c
      call reslen (type(1:1),num,end(1))
c
ccc      if (type(1:1) .eq. 'A') then
ccc        end (1) = 1.41 * float(num-1)
ccc      else
ccc        end (1) = 3.30 * float(num-1)
ccc      end if
c
      call do_spink (start,end,type,num)
c
      return
      end
c
c
c
      subroutine do_spink (start,end,type,num)
c
      include 'moleman2.incl'
c
      real start(3),end(3)
c
      integer num,iat,nold,i
c
      character type*1,ich*1
c
code ...
c
      nold = natoms
c
c ... skip one in residue numbering
c
      if (nold .gt. 0) then
        iat = iresid(natoms) + 1
      else
        iat = 0
      end if
c
      call auto_alpha_beta (type,start,end,num,
     +  natoms,maxatm,atmnam,iresid,iat,xyz,rbuf,ibuf)
c
      if (natoms .le. nold) return
c
      if (nold .gt. 0) then
        iat = atomnr(nold)
        ich = 'Z'
        if (ichar(achain(nold)) .ge. ichar('A') .and.
     +      ichar(achain(nold)) .le. ichar('Y'))
     +    ich = char(ichar(achain(nold))+1)
      else
        iat = 0
        ich = 'Z'
      end if
c
      do i=nold+1,natoms
        iat = iat + 1
        atomnr (i) = iat
        altloc (i) = ' '
        resnam (i) = 'ALA'
        achain (i) = ich
        insert (i) = ' '
        qatom  (i) = 1.0
        batom  (i) = 20.0
        inote  (i) = ' '
        resptr (i) = -1
        lmain  (i) = .false.
        lhet   (i) = .false.
        weight (i) = 1.0
        select (i) = .true.
        laniso (i) = .false.
      end do
c
      call book_keep (.false.)
c
      return
      end
c
c
c
      subroutine auto_sse (iunit,file)
c
      include 'moleman2.incl'
c
      real start(3),end(3)
c
      integer num,iat,nold,i,ierr,length,iunit,ires,nstart
c
      character type*1,ich*1,line*256,file*(*)
      character name*40,res1*6,res2*6,ins*1
c
code ...
c
      if (length(file) .lt. 1) then
        call errcon ('No filename provided')
        return
      end if
c
      call xopxua (iunit,file,linter,ierr)
      if (ierr .ne. 0) return
c
      nstart = natoms
c
   10 continue
      read (iunit,'(a)',end=1000,err=9000) line
      call upcase (line(1:5))
      if (line(1:5) .eq. 'ALPHA') goto 100
      if (line(1:4) .eq. 'BETA') goto 110
      if (line(1:1) .ne. '!') call textut (' >',line)
      goto 10
c
  100 continue
      type = 'A'
      goto 200
c
  110 continue
      type = 'B'
c
  200 continue
      write (*,*)
      call upcase (line)
      call textut (' >',line)
      read (line(6:),*,end=9000,err=9000) name,res1,res2,num,
     +  (start(i),i=1,3),(end(i),i=1,3)
      call detaja (res2,ich,i,ins)
      call detaja (res1,ich,ires,ins)
      if (ires .gt. i) then
        call prompt ('Swapping start and end')
        call rswap (start(1),end(1))
        call rswap (start(2),end(2))
        call rswap (start(3),end(3))
      end if
      nold = natoms
c
      call auto_alpha_beta (type,start,end,num,
     +  natoms,maxatm,atmnam,iresid,ires-1,xyz,rbuf,ibuf)
c
      if (natoms .le. nold) goto 10
c
      do i=nold+1,natoms
        iat = iat + 1
        atomnr (i) = iat
        altloc (i) = ' '
        resnam (i) = 'ALA'
        achain (i) = ich
        insert (i) = ins
        qatom  (i) = 1.0
        batom  (i) = 20.0
        inote  (i) = ' '
        resptr (i) = -1
        lmain  (i) = .false.
        lhet   (i) = .false.
        weight (i) = 1.0
        select (i) = .true.
        laniso (i) = .false.
      end do
c
      goto 10
c
c ... all done
c
 1000 continue
      close (iunit)
      if (natoms .gt. nstart) then
        call jvalut (' New atoms :',1,(natoms-nstart+1))
        call book_keep (.false.)
      end if
c
      return
c
 9000 continue
      call errcon ('While reading SSE file')
      goto 1000
c
      end
c
c
c
      subroutine auto_alpha_beta (type,start,end,num,
     +  natoms,maxatm,atmnam,iresid,iat,xyz,xxyyzz,ibuff)
c
      implicit none
c
      integer num,natoms,maxatm,iat
      integer ibuff(3,*),iresid(*)
c
      real start(3),end(3),x1(3),x2(3)
      real try(3,4),act(3,4),rt(12),rms
      real xyz(3,*),xxyyzz(3,*)
      real d,d1,d2,q1,q2,qq
      real dist
ccc      real det3
c
      integer i,j,ires,nn,ierr
c
      character type*1
      character atmnam(*)*4
c
code ...
c
c ... get desired length of SSE
c
      d = sqrt ( (start(1)-end(1))**2 + (start(2)-end(2))**2 + 
     +           (start(3)-end(3))**2 )
      call fvalut (' Length (A) :',1,d)
c
      if (num .eq. 1) then
        call errcon ('Fewer than two residues')
        return
      end if
c
c ... derive nr of residues
c
      if (num .lt. 1) then
        call reslen (type,num,d)
      end if
c
c ... helix rises 1.41 A per residue
c
ccc      if (num .lt. 1 .and. type .eq. 'A') then
ccc        num = nint(d/1.41)
ccc      else if (num .lt. 1) then
c
c ... strand rises 3.30 A per residue
c
ccc        num = nint(d/3.30)
ccc      end if
c
      call fvalut (' Start CA coordinates :',3,start)
      call fvalut (' End   CA coordinates :',3,end)
      call jvalut (' Nr of residues  :',1,num)
      if (num .lt. 2) then
        call errcon ('Fewer than 2 residues to generate')
        return
      end if
      call jvalut (' Nr of new atoms :',1,5*num)
c
c ... can't get too many atoms
c
      if ( (natoms+5*num) .gt. maxatm) then
        call jvalut (' Old nr of atoms :',1,natoms)
        call jvalut (' New nr of atoms :',1,natoms+5*num)
        call jvalut (' Max nr of atoms :',1,maxatm)
        call errcon ('Too many atoms - sorry !')
        return
      end if
c
      nn = 5*num
      ires = iat
c
      if (type .eq. 'A') then
c
c ... ideal helix geometry (internal coordinates)
c
cBAD      11  N   ALA     3       1.317 117.032 -50.524  1.00 20.00     8     7     6
cBAD      12  CA  ALA     3       1.483 121.842 177.987  1.00 20.00    11     8     7
cBAD      13  C   ALA     3       1.534 111.564 -53.333  1.00 20.00    12    11     8
cBAD      14  O   ALA     3       1.270 119.208 133.221  1.00 20.00    13    12    11
cBAD      15  CB  ALA     3       1.549 109.843 178.252  1.00 20.00    12    11     8
c
        do i=1,nn-1,5
c
          ires = ires + 1
c
c ... N
c
          j = i
          xxyyzz(1,j) = 1.333
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 117.1
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = -50.0
          ibuff (3,j) = max (0,j-5)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' N  '
c
c ... CA
c
          j = i+1
          xxyyzz(1,j) = 1.480
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 121.5
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 180.0
          ibuff (3,j) = max (0,j-5)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' CA '
c
c ... C
c
          j = i+2
          xxyyzz(1,j) = 1.520
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 110.0
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = -60.0
          ibuff (3,j) = max (0,j-5)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' C  '
c
c ... O
c
          j = i+3
          xxyyzz(1,j) = 1.245
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 120.0
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = 133.0
          ibuff (3,j) = max (0,j-3)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' O  '
c
c ... CB
c
          j = i+4
          xxyyzz(1,j) = 1.548
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 110.8
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 180.0
          ibuff (3,j) = max (0,j-7)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' CB '
c
c ... CB first residue different
c
          if (j .eq. 5) then
            xxyyzz(1,j) = 2.568
            ibuff (1,j) = 3
            xxyyzz(2,j) = 32.68
            ibuff (2,j) = 2
            xxyyzz(3,j) = -123.0
            ibuff (3,j) = 1
          end if
c
        end do
c
      else
c
c
c ... ideal strand geometry (internal coordinates)
c
cBAD      11  N   ALA     3       1.305 115.737 147.052  1.00 20.00     8     7     6
cBAD      12  CA  ALA     3       1.466 123.431-179.357  1.00 20.00    11     8     7
cBAD      13  C   ALA     3       1.532 108.628-131.194  1.00 20.00    12    11     8
cBAD      14  O   ALA     3       1.282 115.511 -52.183  1.00 20.00    13    12    11
cBAD      15  CB  ALA     3       1.551 108.713 105.733  1.00 20.00    12    11     8
c
        do i=1,nn-1,5
c
          ires = ires + 1
c
c ... N
c
          j = i
          xxyyzz(1,j) = 1.32
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 115.7
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 120.0
          ibuff (3,j) = max (0,j-5)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' N  '
c
c ... CA
c
          j = i+1
          xxyyzz(1,j) = 1.466
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 121.9
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 180.0
          ibuff (3,j) = max (0,j-5)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' CA '
c
c ... C
c
          j = i+2
          xxyyzz(1,j) = 1.530
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 108.6
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = -125.0
          ibuff (3,j) = max (0,j-5)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' C  '
c
c ... O
c
          j = i+3
          xxyyzz(1,j) = 1.255
          ibuff (1,j) = max (0,j-1)
          xxyyzz(2,j) = 120.0
          ibuff (2,j) = max (0,j-2)
          xxyyzz(3,j) = -40.0
          ibuff (3,j) = max (0,j-3)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' O  '
c
c ... CB
c
          j = i+4
          xxyyzz(1,j) = 1.548
          ibuff (1,j) = max (0,j-3)
          xxyyzz(2,j) = 111.3
          ibuff (2,j) = max (0,j-4)
          xxyyzz(3,j) = 108.9
          ibuff (3,j) = max (0,j-7)
          iresid (natoms+j) = ires
          atmnam (natoms+j) = ' CB '
c
c ... CB first residue different
c
          if (j .eq. 5) then
            xxyyzz(1,j) = 2.487
            ibuff (1,j) = 3
            xxyyzz(2,j) = 35.99
            ibuff (2,j) = 2
            xxyyzz(3,j) = -121.3
            ibuff (3,j) = 1
          end if
c
        end do
c
      end if
c
c ... convert to Cartesian coordinates
c
      call c2cart (nn,xyz(1,natoms+1),xxyyzz,ibuff)
c
c ... get actual distance between first and last CA
c
      d = dist (natoms+2,natoms+(num-1)*5+2,xyz)
      call fvalut (' Actual distance terminal CAs (A) :',1,d)
c
c ... replace "END" by a point on the line START -> END which has
c     distance D to START and is closest to the old END
c
c ... line START -> END has equation:
c
c         ( Xstart )       (Xend - Xstart)
c ... x = ( Ystart ) + L * (Yend - Ystart)
c         ( Zstart )       (Zend - Zstart)
c
c ... there are two points on this line which have distance D to
c     the point N; these occur for values of L = +/- SQRT (P),
c     with P = D**2 / { (Xend - Xstart)**2 + ... }
c
c ... try both roots and keep the one which is closest to the old
c     point END
c
      qq = d*d / ( (start(1)-end(1))**2 + (start(2)-end(2))**2 +
     +             (start(3)-end(3))**2 )
c
c ... + root
c
      q1 = sqrt(qq)
      x1(1) = start(1) + q1*(end(1)-start(1))
      x1(2) = start(2) + q1*(end(2)-start(2))
      x1(3) = start(3) + q1*(end(3)-start(3))
      d1 = sqrt ( (x1(1)-end(1))**2 + (x1(2)-end(2))**2 +
     +            (x1(3)-end(3))**2 )
c
c ... - root
c
      q2 = -q1
      x2(1) = start(1) + q2*(end(1)-start(1))
      x2(2) = start(2) + q2*(end(2)-start(2))
      x2(3) = start(3) + q2*(end(3)-start(3))
      d2 = sqrt ( (x2(1)-end(1))**2 + (x2(2)-end(2))**2 +
     +            (x2(3)-end(3))**2 )
c
ccc      print *,' Q1, D1 = ',q1,d1
ccc      print *,' Q2, D2 = ',q2,d2
c
      if (d1 .lt. d2) then
        try(1,3) = x1(1)
        try(2,3) = x1(2)
        try(3,3) = x1(3)
      else
        try(1,3) = x2(1)
        try(2,3) = x2(2)
        try(3,3) = x2(3)
      end if
c
ccc      call fvalut (' Actual end point (CA) :',3,try(1,3))
c
      try(1,1) = start(1)
      try(2,1) = start(2)
      try(3,1) = start(3)
c
c      try(1,4) = 0.5*(try(1,1)+try(1,3))
c      try(2,4) = 0.5*(try(2,1)+try(2,3))
c      try(3,4) = 0.5*(try(3,1)+try(3,3))
c
      try(1,4) = (try(1,1)+2.0*try(1,3))/3.0
      try(2,4) = (try(2,1)+2.0*try(2,3))/3.0
      try(3,4) = (try(3,1)+2.0*try(3,3))/3.0
c
      act(1,1) = xyz(1,natoms+2)
      act(2,1) = xyz(2,natoms+2)
      act(3,1) = xyz(3,natoms+2)
c
      act(1,3) = xyz(1,natoms+(num-1)*5+2)
      act(2,3) = xyz(2,natoms+(num-1)*5+2)
      act(3,3) = xyz(3,natoms+(num-1)*5+2)
c
c      act(1,4) = 0.5*(act(1,1)+act(1,3))
c      act(2,4) = 0.5*(act(2,1)+act(2,3))
c      act(3,4) = 0.5*(act(3,1)+act(3,3))
c
      act(1,4) = (act(1,1)+2.0*act(1,3))/3.0
      act(2,4) = (act(2,1)+2.0*act(2,3))/3.0
      act(3,4) = (act(3,1)+2.0*act(3,3))/3.0
c
c      q1 = dist(1,3,act)
c      q2 = dist(1,3,try)
c      print *,' ... ACT ',q1
c      print *,' ... TRY ',q2
c
c ... add some "random noise" to the user's desired positions
c     to make sure that the LSQ routine works (otherwise we
c     have three points exactly co-linear -> LSQ fails)
c
c ... 960414 - actually, this gives problems if there are few
c     residues (instabilities: determinant of operator <> 1.0)
c     the following should work provided neither of the two
c     axes has direction vector (1 1 1) ...
c
c      do i=1,3
c        try (i,2) = try(i,2) + 0.1
c        act (i,2) = act(i,2) + 0.1
c      end do
c
c     It is even better to find a point in the plane that bisects
c     the start->end line at a given distance from that line (for
c     both TRY and ACT):
c
c     D = mid-point start->end = (end+start)/2  = (a b c)
c     direction vector start->end = End - Start = (p q r)
c     the midway perpendicular plane is: px + qy + rz = (ap+bq+cr)
c
c     one line in the plane through D perpendicular to the
c     line start->end is: x = lambda.(k l m) + (a b c)
c     select the direction such that k=-q, l=-p, m=0
c     the point on this line which lies g A from point D
c     follows from: dist(x-d) = g
c     --> lambda = +/- g / SQRT (1 + 1 + m**2)
c     take positive root (arbitrary)
c     ==> gives coordinates for our third point in the least-squares
c     superpositioning
c     ==> point has: x = a + lambda*(-q)
c                    y = b + lambda*(p)
c                    z = c
c
c ... in fact the problem may have been that I took the mid-point
c     between start and end which may yield problems; now, instead,
c     we use the point on the start->end line that lies 2/3 along
c     the way from start to end;
c     the distance g from this point is taken equal to the distance
c     between start and end to get three well-separated points for
c     the leasts-squares
c     indeed, now the RMSD for the three points after LSQ is of
c     the order of 10**-6, and the determinant is always +1.0
c     (the distance from this point to start is g.sqrt(13/9) and
c     to end it is g.sqrt(10/9); so they're all of the same order
c     of magnitude)
c
      q1 = d / sqrt ( (try(1,3)-try(1,1))**2 +
     +                   (try(2,3)-try(2,1))**2)
c      print *,' try lambda = ',q1
      try (1,2) = try (1,4) + (try(2,1)-try(2,3))*q1
      try (2,2) = try (2,4) + (try(1,3)-try(1,1))*q1
      try (3,2) = try (3,4)
c
      q1 = d / sqrt ( (act(1,3)-act(1,1))**2 +
     +                   (act(2,3)-act(2,1))**2)
c      print *,' act lambda = ',q1
      act (1,2) = act (1,4) + (act(2,1)-act(2,3))*q1
      act (2,2) = act (2,4) + (act(1,3)-act(1,1))*q1
      act (3,2) = act (3,4)
c
c      q1 = 1.0 / sqrt (
c     +  (xyz(1,natoms+(num-1)*5+2)-xyz(1,natoms+2))**2 +
c     +  (xyz(2,natoms+(num-1)*5+2)-xyz(2,natoms+2))**2 )
c      act (1,2) = act (1,4) +
c     +  (xyz(2,natoms+2)-xyz(2,natoms+(num-1)*5+2))*q1
c      act (2,2) = act (2,4) +
c     +  (xyz(1,natoms+(num-1)*5+2)-xyz(1,natoms+2))*q1
c      act (3,2) = act (3,4)
c
c ... check that both distances are correct:
c
c      print *,' Dact = ',dist (2,4,act),dist (1,2,act),
c     +  dist (3,2,act)
c      print *,' Dtry = ',dist (2,4,try),dist (1,2,try),
c     +  dist (3,2,try)
c
c      do i=1,3       
c        call fvalut (' Before :',3,try(1,i))
c        do j=1,3
c          print *,' I,J = ',i,j
ccc          act(j,i)=act(j,i)-act(j,1)+try(j,1)
c          call gkrand (qq,-0.05,0.05,0)
c          try(j,i)=try(j,i)+qq
c          try(j,i)=try(j,i) - x1(j)
ccc          try(j,i)=10.0*try(j,i)
ccc          act(j,i)=10.0*act(j,i)
c        end do
c        call fvalut (' After  :',3,try(1,i))
c      end do
c      call fvalut (' TRY 1 :',3,try(1,1))
c      call fvalut (' ACT 1 :',3,act(1,1))
c      call fvalut (' TRY 2 :',3,try(1,2))
c      call fvalut (' ACT 2 :',3,act(1,2))
c      call fvalut (' TRY 3 :',3,try(1,3))
c      call fvalut (' ACT 3 :',3,act(1,3))
c
c ... get operator to position SSE according to the user's
c     target coordinates
c
      call lsqgjk (try,act,3,rms,rt,ierr)
c
c      print *,'RMS = ',rms,' - IERR = ',ierr
c
c      call anancs (1,rt,.true.,i)
c
c      q1 = det3 (rt)
c      if (abs(q1-1.0) .gt. 0.00001) then
c        call errcon ('Non-unity determinant')
ccc        goto 6321
c      end if
c
ccc      rt(10) = 0.0
ccc      rt(11) = 0.0
ccc      rt(12) = 0.0
c
      call vecrtv (xyz(1,natoms+1),xyz(1,natoms+1),
     +  nn,rt(1),rt(10))
c
ccc      print *,'FIRST = ',xyz(1,natoms+1),xyz(2,natoms+1),
ccc     +  xyz(3,natoms+1)
c      do i=1,nn
c        do j=1,3
c          xyz(j,natoms+i) = xyz(j,natoms+i) + try(j,2)
c        end do
c      end do
c
      call fvalut (' Actual N-term CA :',3,xyz(1,natoms+2))
      call fvalut (' Actual C-term CA :',3,xyz(1,natoms+(num-1)*5+2))
c
      q1 = 0.0
      do i=1,3
        q1 = q1 + (xyz(i,natoms+2)-start(i))**2 +
     +            (xyz(i,natoms+(num-1)*5+2)-end(i))**2
      end do
      q1 = sqrt (0.5 * q1)
      call fvalut (' RMSD (A) :',1,q1)
c
      natoms = natoms + nn
c
      return
      end
c
c
c
      subroutine reslen (type,num,leng)
c
      implicit none
c
c ... distance covered per residue for very very long SSEs
c
      real arise,brise
      parameter (arise=1.41, brise=3.30)
c
      real alen(10),blen(10)
      real leng,dist,q
c
      integer num,i,inow
c
      character type*1
c
      data alen /-1.0,  3.83,  5.54,  5.09,  5.75,
     +            8.20, 9.63, 10.09, 11.48, 13.52/
c
      data blen /-1.0,   3.81,  6.61, 10.09, 13.22,
     +           16.63, 19.83, 23.21, 26.44, 29.80/
c
code ...
c
      if (type .eq. 'A') then
        if (num .ge. 2 .and. num .le. 10) then
          leng = alen(num)
          return
        else if (num .gt. 10) then
          leng = arise * float(num)
          return
        else if (num .le. 1) then
          if (leng .gt. alen(9)) then
            num = nint (leng/arise)
            return
          else
            inow = 2
            dist = abs(leng-alen(2))
            do i=3,10
              q = abs(leng-alen(i))
              if (q .le. dist) then
                inow = i
                dist = q
              end if
            end do
            num = inow
            return
          end if
        end if
      else if (type .eq. 'B') then
        if (num .ge. 2 .and. num .le. 10) then
          leng = blen(num)
          return
        else if (num .gt. 10) then
          leng = brise * float(num)
          return
        else if (num .le. 1) then
          if (leng .gt. blen(9)) then
            num = nint (leng/brise)
            return
          else
            inow = 2
            dist = abs(leng-blen(2))
            do i=3,10
              q = abs(leng-blen(i))
              if (q .le. dist) then
                inow = i
                dist = q
              end if
            end do
            num = inow
            return
          end if
        end if
      end if
c
c ... if here, error
c
      return
      end
c
c
c
      subroutine ono_molray_trace (iunit,trafil,mode,xmult,neafar,
     +                             forrev,trares,tratom,bufxyz)
c
      include 'moleman2.incl'
c
      real bufxyz(3,*)
      real cog(3),vec(3),x1(3),x2(3),xprev(3),rtemp(3,4)
      real distce,xmult,dx1,dx2
c
      integer i,nok,j,k,iunit,ierr,length,ndone,iuse
c
      logical xinter,lnear,lfor
c
      character trafil*(*),mode*(*),neafar*(*),forrev*(*)
      character trares*3,tratom*4
c
code ...
c
      write (*,*)
      if (length(trafil) .lt. 1) then
        call errcon ('No output file name provided')
        return
      end if
      call textut (' Opening output file :',trafil)
      call xopxua (iunit,trafil,xinter(),ierr)
      if (ierr .ne. 0) return
c
      call upcase (mode)
      call textut (' Mode :',mode)
      if (mode(1:4) .ne. 'PEPT') then
        call errcon ('Unknown mode')
        return
      end if
c
      call fvalut (' Multiplier for normal vector :',1,xmult)
c
      call upcase (neafar)
      call textut (' Near/Far :',neafar)
      if (neafar(1:3) .ne. 'NEA' .and.
     +    neafar(1:3) .ne. 'FAR') then
        call errcon ('Unknown nearest/farthest mode')
        return
      end if
      lnear = (neafar(1:3) .eq. 'NEA')
c
      call upcase (forrev)
      call textut (' Forward/reverse :',forrev)
      if (forrev(1:3) .ne. 'FOR' .and.
     +    forrev(1:3) .ne. 'REV') then
        call errcon ('Unknown forward/reverse mode')
        return
      end if
      lfor = (forrev(1:3) .eq. 'FOR')
c
      call upcase (trares)
      call textut (' Trace residue type :',trares)
      call upcase (tratom)
      call textut (' Trace atom type :',tratom)
c
      nok = 0
      ndone = 0
c
c ... mode PEPT (only defined for proteins)
c
      if (mode(1:4) .ne. 'PEPT') goto 2000
c
      do i=1,3*nres
        ibuf (i) = -1
      end do
c
      nok = 0
      do i=1,nres
        lbuf (i) = (restyp(i).eq.iprot)
        if (lbuf(i)) then
          k = 0
          do j=atmptr(1,i),atmptr(2,i)
            if (select(j)) then
              if (atmnam(j) .eq. ' N  ') then
                ibuf (i) = j
                k = k + 1
              else if (atmnam(j) .eq. ' CA ') then
                ibuf (nres+i) = j
                k = k + 1
              else if (atmnam(j) .eq. ' C  ') then
                ibuf (2*nres+i) = j
                k = k + 1
              end if
            end if
          end do
          if (k .gt. 0) nok = nok + 1
        end if
      end do
c
      call jvalut (' Nr of residues to check :',1,nok)
      if (nok .lt. 3) then
        call errcon ('I don''t call this a protein ...')
        return
      end if
c
      ndone = 0
c
      do i=1,nres-1
c
        if (ibuf(nres+i) .gt. 0 .and.
     +      ibuf(2*nres+i) .gt. 0 .and.
     +      ibuf(i+1) .gt. 0 .and.
     +      ibuf(nres+i+1) .gt. 0) then
c
          ndone = ndone + 1
c
          do j=1,3
            rtemp (j,1) = xyz(j,ibuf(nres+i))
            rtemp (j,2) = xyz(j,ibuf(2*nres+i))
            rtemp (j,3) = xyz(j,ibuf(1+i))
            rtemp (j,4) = xyz(j,ibuf(nres+1+i))
          end do
c
          call lsplan (4,rtemp,cog,vec)
c
          do j=1,3
            x1(j) = cog(j) + xmult*vec(j)
            x2(j) = cog(j) - xmult*vec(j)
          end do
c
          iuse = 1
          if (ndone .gt. 1) then
            dx1 = distce (xprev,x1)
            dx2 = distce (xprev,x2)
            if (lnear) then
              if (dx2 .lt. dx1) iuse = 2
            else
              if (dx2 .gt. dx1) iuse = 2
            end if
ccc            write (*,'(2x,i8,2f8.2,i3)') ndone,dx1,dx2,iuse
          end if
c
          if (iuse .eq. 1) then
ccc            write (iunit,6000) ndone,tratom,trares,ndone,x1
            do j=1,3
              bufxyz (j,ndone) = x1(j)
              xprev (j) = x1 (j)
            end do
          else
ccc            write (iunit,6000) ndone,tratom,trares,ndone,x2
            do j=1,3
              bufxyz (j,ndone) = x2(j)
              xprev (j) = x2 (j)
            end do
          end if
c
 6000 format ('ATOM  ',i5,1x,a4,1x,a3,1x,i5,4x,3f8.3,
     +  '  1.00 20.00   8')
c
cATOM    163  CB  LYS    20      62.798   1.269   0.832  1.00 20.00   6
c1234567890123456789012345678901234567890123456789012345678901234567890
c
        end if
      end do
c
      goto 9000
c
c ... next mode (if required in the future)
c
 2000 continue
      return
c
c ... write PDB file
c
 9000 continue
c
      call jvalut (' Nr of trace atoms found :',1,ndone)
      if (ndone .lt. 1) return
c
      if (lfor) then
        do i=1,ndone
          write (iunit,6000) i,tratom,trares,i,
     +      bufxyz(1,i),bufxyz(2,i),bufxyz(3,i)
        end do
      else
        i = 0
        do j=ndone,1,-1
          i = i + 1
          write (iunit,6000) i,tratom,trares,i,
     +      bufxyz(1,j),bufxyz(2,j),bufxyz(3,j)
        end do
      end if
c
      write (iunit,'(a6)') 'END   '
c
      call prompt (' PDB file written')
c
      return
      end
c
c
c
      subroutine lsplan (nok,rbuf,cog,vec)
c
      implicit none
c
      integer nok
c
      real rbuf(3*nok),cog(3),vec(3)
      real a(3,3),c(3),d(3),v(3,3)
      real q11,q22,q33,q12,q13,q23,rms,x1,x2,x3
      real xmin,xmax,ymin,ymax,zmin,zmax,xmag
      real det3
c
      integer i,j,k
c
code ...
c
      do i=1,3
        c(i) = 0.0
      end do
      xmin = rbuf(1)
      xmax = rbuf(1)
      ymin = rbuf(2)
      ymax = rbuf(2)
      zmin = rbuf(3)
      zmax = rbuf(3)
      do i=1,nok
        j = 3*(i-1)
        c(1) = c(1) + rbuf(j+1)
        c(2) = c(2) + rbuf(j+2)
        c(3) = c(3) + rbuf(j+3)
        xmin = min (xmin, rbuf(j+1))
        xmax = max (xmax, rbuf(j+1))
        ymin = min (ymin, rbuf(j+2))
        ymax = max (ymax, rbuf(j+2))
        zmin = min (zmin, rbuf(j+3))
        zmax = max (zmax, rbuf(j+3))
      end do
c
      do i=1,3
        c(i) = c(i)/float(nok)
        cog (i) = c(i)
      end do
ccc      call fvalut (' Centre of Gravity :',3,c)
      do i=1,nok
        j = 3*(i-1)
        rbuf (j+1) = rbuf(j+1) - c(1)
        rbuf (j+2) = rbuf(j+2) - c(2)
        rbuf (j+3) = rbuf(j+3) - c(3)
      end do
c
      q11=0.0
      q22=0.0
      q33=0.0
      q12=0.0
      q13=0.0
      q23=0.0
      do i=1,nok
        j = 3*(i-1)
        x1 = rbuf(j+1)
        x2 = rbuf(j+2)
        x3 = rbuf(j+3)
        q11=q11+x1*x1
        q22=q22+x2*x2
        q33=q33+x3*x3
        q12=q12+x1*x2
        q13=q13+x1*x3
        q23=q23+x2*x3
      end do
c
      a(1,1)=q11
      a(1,2)=q12
      a(2,1)=q12
      a(1,3)=q13
      a(3,1)=q13
      a(2,2)=q22
      a(2,3)=q23
      a(3,2)=q23
      a(3,3)=q33
c
      call jacobi (a,3,3,d,v,i,j)
      if (j .ne. 0) then
        call errcon ('Fatal error in eigenvalue determination')
        return
      end if
c
      call eigsrt (d,v,3,3)
c
      k = 0
 1234 continue
ccc      do i=1,3
ccc        write (*,6000) i,d(i),(v(j,i),j=1,3)
ccc      end do
 6000 format (' Eigen value ',i1,' = ',f12.1,' Vector : ',3f10.6)
c
      rms = det3(v)
ccc      call fvalut (' Determinant :',1,rms)
c
      if (rms .lt. 0.0) then
        if (k.eq.0) then
ccc          call errcon (
ccc     +      'Negative determinant; change hand of inertia axes')
          do i=1,3
            do j=1,3
              v(i,j) = -v(i,j)
            end do
          end do
          k = 1
          goto 1234
        else
          call errcon ('Cannot fix negative determinant !!!???')
          return
        end if
      end if
c
c ... return unit normal vector
c
      xmag = 0.0
      do i=1,3
        vec (i) = v(i,3)
        xmag = xmag + (vec(i)*vec(i))
      end do
      xmag = sqrt(xmag)
      do i=1,3
        vec (i) = vec (i) / xmag
      end do
c
ccc      call fvalut (' Unit normal vector :',3,vec)
c
      return
      end
