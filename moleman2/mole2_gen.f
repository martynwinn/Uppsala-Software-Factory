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
