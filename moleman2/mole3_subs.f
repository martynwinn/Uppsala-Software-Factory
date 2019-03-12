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
