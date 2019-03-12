      program seaman
c
c ... SEAMAN - SEArch-model MANipulator for Molecular Replacement
c
c ... Gerard Kleywegt @ 950622
c
c ... 0.1 @ 950622 - first version
c
      include 'seaman.incl'
c
      integer iunit
      parameter (iunit=11)
c
      real blim(2),qlim(2)
      real total,user,sys,xave,xsdv,xmin,xmax,xtot,factor
c
      integer zone(2),selzon(2)
      integer i,j,k,ierr,nn,seresi,length
c
      logical xinter,linter
c
      character line*256,comand*4,file*128,saved*128,selchn*1
      character change1*3,change2*3,seltyp*3
c
code ...
c
      call gkinit (prognm,vers)
c
      call jvalut (' Max nr of residue types :',1,maxtyp)
      call jvalut (' Max nr of type atoms    :',1,maxatp)
      call jvalut (' Max nr of residues      :',1,maxres)
      call jvalut (' Max nr of atoms         :',1,maxatm)
      call jvalut (' Max nr of chains/models :',1,maxchn)
      write (*,*)
c
      linter = xinter()
      file = 'm1.pdb'
      saved = 'not_saved_yet'
      natoms = 0
      comand = 'READ'
c
      blim (1) = 1.0
      blim (2) = 60.0
      qlim (1) = 0.99
      qlim (2) = 1.01
      zone (1) = 9998
      zone (2) = 9999
      selzon (1) = 0
      selzon (2) = 0
      seresi = 1
      selchn = '*'
      change1= 'VAL'
      change2= 'ALA'
      seltyp = 'LYS'
      defq = 1.0
      defb = 20.0
      factor = 1.0
c
c ... read library
c
      line = 'seaman.lib'
      call gklibf (line)
c
      write (*,*)
      call textin (' Name of library file ?',line)
      call xopxoa (iunit,line,linter,ierr)
      if (ierr .ne. 0) then
        call errcon ('While opening SEAMAN library file')
        goto 9000
      end if
c
      call prompt ('0Reading library ...')
      call readdb (iunit,ierr)
      if (ierr .ne. 0) then
        call errcon ('While reading SEAMAN library file')
        goto 9000
      end if
c
c ... formats
c
 6000 format (/
     +  ' SEAMAN options :'//
     +  ' ? (list options)                    ',
     +              ' ! (comment)'/
     +  ' QUIT program                        ',
     +              ' $ issue shell command'/
     + /' READ molecule from PDB file         ',
     +              ' WRITe molecule to PDB file'/
     +  ' LIST molecule statistics            ',
     +              ' SEQUence of molecule'/
     + /' BFACtor delete                      ',
     +              ' OCCUpancy delete'/
     +  ' ISOLated atoms delete               ',
     +              ' ZONE delete'/
     +  ' LOOPs delete                        ',
     +              ' TURNs delete'/
     + /' GLYCine zone                        ',
     +              ' ALANine zone'/
     +  ' SERIne zone                         ',
     +              ' RESIdue type to Gly/Ala/Ser'/
     +  ' MINImalist substitutions            ',
     +              ' MUTAte replace rotamer'/
     +  ' ROTAmer residue                     ',
     +              ' '/
     + )
c
c --- LIST OPTIONS
c
  100 continue
      write (*,6000)
c
c --- MAIN COMMAND LOOP
c
   10 continue
      call gkdcpu (total,user,sys)
      if (total .ge. 1.0)
     +  write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      call textin (' Command ?',comand)
      if (length(comand) .lt. 1) goto 10
c
      call upcase (comand)
c
      if (comand(1:1) .eq. '!') goto  10
c
c ... system command
c
      if (comand(1:1) .eq. '$') then
        line = ' '
        call textin (' System command ?',line)
        call gksys (line)
        goto 10
      end if
c
      if (comand(1:1) .eq. '?') goto 100
c
c ... QUIT
c
      if (comand .eq. 'QUIT') then
c
        goto 9000
c
c ... WRITE
c
      else if (comand .eq. 'WRIT') then
c
        call textin (' PDB file name ?',saved)
c
        call xopxua (iunit,saved,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening PDB file')
          goto 10
        end if
c
        call putpdb (iunit,ierr)
c
        goto 10
c
c ... READ
c
      else if (comand .eq. 'READ') then
c
        call textin (' PDB file name ?',file)
c
        call xopxoa (iunit,file,linter,ierr)
        if (ierr .ne. 0) then
          call errcon ('While opening PDB file')
          goto 10
        end if
c
ccc        call fvalin (' Sensitivity for P-SEA ?',1,factor)
c
        call getpdb (iunit,ierr)
        saved = 'not_saved_yet'
        if (ierr .ne. 0) then
          call errcon ('While reading PDB file')
          natoms = 0
          goto 10
        end if
c
        call xstats (batom,natoms,defb,xsdv,xmin,xmax,xtot)
c
        call bookkp (ierr)
        if (ierr .ne. 0) then
          call errcon ('Sorry')
          natoms = 0
          goto 10
        end if
c
        call yasspa ()
ccc        call psea (factor)
        call surfer ()
c
        goto 10
c
      end if
c
      if (natoms .lt. 1) then
        call errcon ('No molecule read yet')
        goto 10
      end if
c
 6010 format (1x,a12,' : Ave, Sdv : ',2f8.2/
     +           13x,'   Min, Max : ',2f8.2)
c
 6020 format (/' Type      : ',10(1x,a3,3x))
 6021 format ( ' Name      : ',10(1x,a6))
 6022 format ( ' Structure : ',10(1x,a6))
 6023 format ( ' Surface ? : ',10(1x,l1,5x))
 6024 format ( ' <Bresi>   : ',10(1x,f5.1,1x))
 6025 format ( ' # atoms   : ',10(1x,i2,4x))
c
c ... LIST
c
      if (comand .eq. 'LIST') then
c
        write (*,*)
        call textut (' Model read from file :',file)
        call textut (' Last saved to file   :',saved)
        call jvalut (' Nr of atoms in molecule     :',1,natoms)
        call jvalut (' Nr of chains or models      :',1,nchain)
        call jvalut (' Nr of residues              :',1,nres)
c
        call xstats (batom,natoms,xave,xsdv,xmin,xmax,xtot)
        write (*,6010) 'B-factors',xave,xsdv,xmin,xmax
c
        call xstats (qatom,natoms,xave,xsdv,xmin,xmax,xtot)
        write (*,6010) 'Occupancies',xave,xsdv,xmin,xmax
c
        write (*,*)
c
        goto 10
c
c ... SEQUence
c
      else if (comand .eq. 'SEQU') then
c
        do i=1,nres
          j=resptr(2,i)-resptr(1,i)+1
          natres(i)=j
          if (j .gt. 0) then
            call xstats (batom(resptr(1,i)),j,
     +                   xave,xsdv,xmin,xmax,xtot)
            avbres(i)=xave
          else
            avbres(i)=-999.999
          end if
        end do
c
        do i=1,nres,10
          j=min(i+9,nres)
          write (*,6020) (resnam(resptr(1,k)),k=i,j)
          write (*,6021) (namres(k),k=i,j)
          write (*,6022) (struct(k),k=i,j)
          write (*,6024) (avbres(k),k=i,j)
          write (*,6025) (natres(k),k=i,j)
ccc          write (*,6023) (surfac(k),k=i,j)
        end do
c
        write (*,*)
c
        goto 10
c
c ... BFACtor delete
c
      else if (comand .eq. 'BFAC') then
c
        call prompt (' Provide B-factor limits for atoms to KEEP')
        call fvalin (' Limits ?',2,blim)
        call bfadel (blim(1),blim(2))
c
        goto 10
c
c ... OCCUpancy delete
c
      else if (comand .eq. 'OCCU') then
c
        call prompt (' Provide occupancy limits for atoms to KEEP')
        call fvalin (' Limits ?',2,qlim)
        call occdel (qlim(1),qlim(2))
c
        goto 10
c
c ... ZONE delete
c
      else if (comand .eq. 'ZONE') then
c
        call selecm (selzon,selchn,nn)
        if (nn .gt. 0) call zondel ()
c
        goto 10
c
c ... ISOLated atoms delete
c
      else if (comand .eq. 'ISOL') then
c
        call isodel ()
c
        goto 10
c
c ... LOOPs delete
c
      else if (comand .eq. 'LOOP') then
c
        call loodel ()
c
        goto 10
c
c ... TURNs delete
c
      else if (comand .eq. 'TURN') then
c
        call turdel ()
c
        goto 10
c
c ... GLYCine zone
c
      else if (comand .eq. 'GLYC') then
c
        call selecm (selzon,selchn,nn)
        if (nn .gt. 0) call polgly ()
c
        goto 10
c
c ... ALANine zone
c
      else if (comand .eq. 'ALAN') then
c
        call prompt (' Glycines left intact')
        call selecm (selzon,selchn,nn)
        if (nn .gt. 0) call polala ()
c
        goto 10
c
c ... SERIne zone
c
      else if (comand .eq. 'SERI') then
c
        call prompt (' Glycines and alanines left intact')
        call selecm (selzon,selchn,nn)
        if (nn .gt. 0) call polser ()
c
        goto 10
c
c ... RESIdue type to Gly/Ala/Ser
c
      else if (comand .eq. 'RESI') then
c
        call prompt (' Glycines and alanines left intact')
        call selecr (seltyp,nn)
        if (nn .le. 0) goto 10
c
        call textin (' Change to Gly/Ala/Ser ?',change2)
        call upcase (change2)
c
        if (change2(1:1) .eq. 'G') then
          call polgly ()
        else if (change2(1:1) .eq. 'A') then
          call polala ()
        else if (change2(1:1) .eq. 'S') then
          call polser ()
        else
          call errcon ('Invalid residue type')
        end if
c
        goto 10
c
c ... MINImalist substitutions
c
      else if (comand .eq. 'MINI') then
c
        call selec1 (seresi,selchn,nn)
        if (nn .le. 0) goto 10
c
        call textin (' Change to type ?',change1)
        call upcase (change1)
        if (change1.eq.'GLY' .or. change1.eq.'ALA' .or.
     +      change1.eq.'SER') then
          call errcon (' Cannot change to Gly/Ala/Ser')
          call prompt (' Use option GLYC/ALAN/SERI to do that')
          goto 10
        end if
c
        do i=1,nres
          j=resptr(1,i)
          if (selres(i)) then
            call textut (' MINI :',namres(i))
            call minima (i,change1)
          end if
        end do
c
        goto 10
c
c ... MUTAte replace
c
      else if (comand .eq. 'MUTA') then
c
        call selec1 (seresi,selchn,nn)
        if (nn .le. 0) goto 10
c
        call textin (' Change to type ?',change1)
        call upcase (change1)
        if (change1.eq.'GLY') then
          call errcon (' Cannot change to Gly')
          call prompt (' Use option GLYC to do that')
          goto 10
        end if
c
        do i=1,nres
          j=resptr(1,i)
          if (selres(i)) then
            call textut (' MUTA :',namres(i))
            call mutate (i,change1)
          end if
        end do
c
        goto 10
c
c ... ROTAmer residue
c
      else if (comand .eq. 'ROTA') then
c
        call selec1 (seresi,selchn,nn)
        if (nn .le. 0) goto 10
c
        do i=1,nres
          j=resptr(1,i)
          if (selres(i)) then
            call textut (' ROTA :',namres(i))
            call mutate (i,resnam(j))
          end if
        end do
c
        goto 10
c
      end if
c
c ... invalid command
c
      call errcon ('Command not recognised')
      call textut (' Command :',comand)
c
      goto 10
c
c ... all done
c
 9000 continue
      call gkquit ()
c
      end
