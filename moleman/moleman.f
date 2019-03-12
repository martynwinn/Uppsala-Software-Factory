      program moleman
c
c ... manipulate PDB files - the most dreadful piece of Fortran
c     code I've ever written
c
c ... Gerard Kleywegt @ 920330
c
c     Department of Molecular Biology
c     BMC, University of Uppsala
c     Uppsala, Sweden
c
c --- 1992 ---
c ... Modified 920424,27
c ... Modified 920504,05
c ... Modified 920630
c ... Modified 920702,07,10,15,24,29
c ... Modified 920810,21,24,25,27
c ... Modified 920901,04,07,11
c ... Modified 921022
c ... Modified 921102
c --- 1993 ---
c ... Modified 930118
c ... Modified 930202,05,06,11,12,14,22
c ... Modified 930305
c ... Modified 930415
c ... Modified 930519
c ... Modified 930607,11,16,18,23,29
c ... Modified 930805,23
c ... Modified 930929,30
c ... Modified 931006,28
c ... Modified 931109
c --- 1994 ---
c ... Modified 940214
c ... Modified 940307,16,23
c ... Modified 940415
c ... Modified 940525
c ... Modified 940803,12,14,15,16,17,29
c ... Modified 940906
c ... Modified 941219,23,30
c --- 1995 ---
c ... Modified 950101,02,24
c ... Modified 950201,07,15,23,27
c ... Modified 950311,30,31
c ... Modified 950412,13,17,21
c ... Modified 950505,06,30
c ... Modified 950609,14,21,30
c ... Modified 950705,15,23
c ... Modified 950829
c ... Modified 950922,23,28,29
c ... Modified 951030
c ... Modified 951102,10
c --- 1996 ---
c ... Modified 960208,09
c ... Modified 960405
c ... Modified 960712
c ... Modified 961121
c ... Modified 961218
c --- 1999 ---
c ... Modified 991230
c
c REMARK   2 RESOLUTION. 2.6  ANGSTROMS.                                  1GUH  17
c REMARK   2 RESOLUTION. NOT APPLICABLE.  SEE REMARK 4.                   1TRX  32
c
      implicit none
c
      character prognm*12,version*12
      parameter (prognm = 'MOLEMAN', version = '041001/7.4')
c
      integer f1,f2,f3,f4,f5,f7,f8,f9
      parameter (f1=11,f2=12,f3=13,f4=14,f5=15,f7=17,f8=18,f9=19)
c
      include 'mole_dim.incl'
c
      real xf(maxatm),yf(maxatm),zf(maxatm),batom(maxatm),qatom(maxatm)
      real bfacts(2),qfacts(2),xsym(4),ysym(4),zsym(4),xn,yn,zn
      real bsum,bave,bsid,basc,btotal,shifts(5),rmses(6),avers(5)
      real xshf(3),dx,dp,bold,ave,sdv,xmin,xmax,xtot,xbuff(maxbuf)
      real blim(2),qlim(2),ymin,zmin,dmin,dij,rotmat(3,4),trans(3)
      real otxyz(3,maxapr),otint(3,5),dist,angle,codist,ocoang,tangle
      real dismat(maxapr,maxapr),bndcut,allxyz(3,maxatm),dummy
      real cell(6),otof(3,3),ftoo(3,3),ftrans(3),xfra(3),xcar(3)
      real rotang(3),phipsi(2,maxres),ca2fa(12),bstats(7),qstats(7)
      real cog(3),b1(maxatm),b2(maxatm),xxyyzz(3,maxatm),nobcut
      real ymax,xf1,xf2,yf1,yf2,strad,ssdist,xx,valtor(maxtor)
      real volume,badcut,start(3),end(3),total,user,sys,q,rog,som
      real bmax(7),radius
c
      integer natoms,atomnr(maxatm),iresid(maxatm),ierr,oldres,ires
      integer insert(maxatm),i,nlines,kk,length,mave,inum,nresat,nres
      integer iold,id,nats,j,nsum,nsid,nmulti,nalone,ncopies,inew
      integer icopies(maxatm),idaddy,iat,nline,lenbas,nwrit,iseed
      integer ibuff(maxbuf),resptr(maxres),nwater,lowest,jold,jnew
      integer icom,iz1,iz2,iz3,ndum,izone(2),ldum,iot(maxapr)
      integer otref(3,maxapr),igeo,napr,k,l,jres,ncnt,sorted(maxatm)
      integer iref,jref,itra,irot,bqcnts(7),nowres,i1,j1,j2,leng1
      integer n1,ca1,cb1,sg1,n2,ca2,cb2,sg2,inb,iflag,ns,nss
      integer ntyp,typcnt(maxtyp),deftor(4,maxtor),nhskip,zmol
      integer caqual(0:60,-60:60)
c
      logical mainch,baveok(maxatm),lkeeph,lalwyn,lforce,ldummy
      logical llabel,lnorem,lhydro,xinter,lrema,lwater,lhread,ldump
      logical hd21,hd22,hh11,hh12,hh21,hh22,lstrip,lswapm,ldorem
      logical lecho
      logical isbond(maxapr,maxapr),okbond(maxapr,maxapr)
      logical afftor(maxapr,maxtor)
c
      character file1*80,file2*80,file3*80,base*40,myres*3
      character line*1024,key*6,option*40,nowchn*2,which*1
      character atmnam(maxatm)*4,altloc(maxatm)*1,keephs*1
      character resnam(maxatm)*3,achain(maxatm)*2,newchn*2
      character inote(maxatm)*40,aresid(maxatm)*6,aforce*1
      character borq*1,file4*80,file5*80,file8*80,newchx*4
      character amain*1,resold*3,resnew*3,file6*80,alabel*1
      character awater*3,commnt*72,anorem*1
      character whichn*2,file7*80,names(0:maxapr)*4,xnote*40
      character amacid*200,onelet*200,stars*80,pirnam*20
      character pirtit*80,pircod(maxres)*1,arema*1,astrip*1
      character answer*1,mplane*1,ligand*200,aswapm*1,junk80*80
      character omol*6,flagfn*80,s1name*7,s2name*7,pres*3
      character typnam(maxtyp)*3,tornam(maxtor)*6
      character spgrp*11
c
      data bfacts/1.0,0.0/, blim /2.0, 99.9/, trans /3*0.0/
      data qfacts/1.0,0.0/, qlim /0.0, 1.0/, ftrans /3*0.0/
      data xsym /1.0,0.0,0.0,0.0/, ysym /0.0,1.0,0.0,0.0/
      data zsym /0.0,0.0,1.0,0.0/, cell /3*1.0,3*90.0/
      data shifts /0.1,0.1,0.1,0.0,0.0/, rotang /3*0.0/
      data ca2fa /12*0.0/
      data newchn /'MC'/, keephs /'N'/, aforce /'Y'/
      data mave /1/, iseed/0/, borq /'B'/, ires,jres /1,2/
      data resold /'SOL'/, resnew /'HOH'/, alabel /'N'/
      data awater /'HOH'/, iz1,iz2,iz3 /1,1,1/, igeo/1/
      data izone /0,0/, commnt /'MoleMan PDB file'/
      data anorem /'Y'/, codist /1.23/, ocoang /122.5/
      data bndcut /2.0/, whichn/'**'/, itra,irot /1,1/
      data pirnam /'SEQUENCE'/, pirtit /'Any title'/
      data arema /'Y'/, astrip /'N'/, mplane /'X'/, inb /1/
      data ligand /'???'/, aswapm /'N'/, nobcut /3.6/
      data iflag /1/, omol/'M1'/, strad /0.2/, ssdist /2.5/
      data nowchn /' A'/, badcut /2.4/, start /3*0.0/
      data end /3*10.0/, line /' '/
c
      equivalence (xbuff(1),ibuff(1))
c
code ...
c
      call gkinit (prognm,version)
c
c ... flag CAQUAL as uninitialised
c
      caqual (0,0) = -1
c
      call jvalut (' Maximum number of atoms     :',1,maxatm)
      call jvalut (' Maximum number of residues  :',1,maxres)
      call jvalut (' Max nr of atoms per residue :',1,maxapr)
      call jvalut (' Maximum buffer (words)      :',1,maxbuf)
      call jvalut (' Max nr of non-ATOM records  :',1,maxcom)
c
      amacid = 'ALA ARG ASN ASP CYS GLN GLU GLY HIS '//
     +         'ILE LEU LYS MET PHE PRO SER THR TRP '//
     +         'TYR VAL CPR ASX GLX UNK CYH CSS PCA'
      onelet = 'A   R   N   D   C   Q   E   G   H   '//
     +         'I   L   K   M   F   P   S   T   W   '//
     +         'Y   V   J   B   Z   X   C   C   Q   '
c
      call pretty (amacid)
      call pretty (onelet)
c
      natoms = 0
      nrem = 0
      option = 'READ_pdb_file'
      file1 = 'in.pdb'
      file2 = 'out.pdb'
      base  = 'out'
      flagfn = 'atom_flag.odb'
c
      which = 'A'
      amain = 'A'
      newchx = 'AAAA'
c
      spgrp = 'P 1'
      zmol = 1
c
      file4 = 'atom_b.plt'
      file5 = 'resi_b.plt'
      file6 = 'ca_dist.pl2'
      file7 = 'export.bad'
      file8 = 'sequence.pir'
c
      write (*,*)
      call gkrand (dx,0.0,0.0,-1)
c
 4711 format (a6,a)
 4713 format (i5,1x,a4,a1,a3,a2,i4,a1,3x,3f8.3,2f6.2,a40)
 4715 format (i5,1x,a4,a1,a3,   a6,a1,3x,3f8.3,2f6.2,a40)
 6666 format (/' *** ERROR *** ',a/)
c
   10 continue
      write (*,'(99(a/:))') ' ',
     +  ' ===> INPUT/OUTPUT/CONTROL:',
     +  ' ?                         !',
     +  ' READ_pdb_file             NO_H_read',
     +  ' ALWYn_format_read         APPEnd_pdb_file',
     +  ' WRITe_pdb_file            SPLIt_pdb_file',            
     +  ' EXPOrt_bad_file           IMPOrt_bad_file',
     +  ' SAME_export               QUIT',
     +  ' HELIx_generate            STRAnd_generate',
     +  ' DUMP_pdb_file             ECHO_pdb_file',
     +  ' ',
     +  ' ===> MISCELLANEOUS:',
     +  ' REMArk_etc_cards          CRYStal_PDB_card',
     +  ' TALLy_residues            COUNt_elements',
     +  ' GLYCo_sites               EXTInction_280',
     +  ' WATEr_sort                PIR_sequence_file',
     +  ' SSBOnd_records',
     +  ' ',
     +  ' ===> O FILES:',
     +  ' RSR_datablock             CONNect_file',
     +  ' TORSion_datablock         RSFIt_datablock',
     +  ' DISUlfide_ODL_file        FIT_water_macro',
     +  ' FLAG_colours',
     +  ' ',
     +  ' ===> ANALYSIS:',
     +  ' STATistics                CA_Distance_plot',
     +  ' PLOT_Bs_or_Qs             LIST_residue',
     +  ' GEOMetry_list             SEQUence_list',
     +  ' RAMAchandran_plot         BALAsubramanian_plot',
     +  ' CA_Ramachandran_plot      CACA_distances',
     +  ' PLANar_peptides           BURIed_charges',
     +  ' DISTance_distribution     SHORt_contacts',
     +  ' RADIal_B_plot             CHI_list',
     +  ' ',
     +  ' ===> TEMPERATURE FACTORS & OCCUPANCIES:',
     +  ' LIMIt_B_and_Q             AVERage_temp_factors',
     +  ' TEMP_factors_set          OCCUpancies_set',
     +  ' SMOOth_Bs                 B_Q_statistics',
     +  ' BONDed_Bs                 NONBonded_Bs',
     +  ' ',
     +  ' ===> ATOMS, RESIDUES, CHAINS & SEGMENTS:',
     +  ' O2XHydrogens              SUGGest_OT2',
     +  ' CHECk_nomenclature        CORRect_nomenclature',
     +  ' CHAIn_name                XPLOr_ids',
     +  ' FROM_chain_to_XID         XID_to_chain',
     +  ' AUTO_chain_segid          ASK_auto_chain_segid',
     +  ' RENUmber_atoms            ALTEr_residue_name',
     +  ' RESIdu_renumber           ZONE_renumber',
     +  ' ATRE_number',
     +  ' ',
     +  ' ===> MANIPULATION OPTIONS:',
     +  ' FRACtional_to_cartesian   CARTesian_to_fractional',
     +  ' ROTAte_molecule           TRANslate_molecule',
     +  ' APPLy_random_rotation     RANDom_shifts',
     +  ' ORIGin_move               MIRRor_zone',
     +  ' INVErt_zone',
     +  ' '
c
   11 continue
c
      call gkdcpu (total,user,sys)
      if (total .ge. 1.0)
     +  write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      call textin (' Option ?',option)
      if (length(option) .lt. 1) goto 11
      if (option(1:1) .eq. '!') goto 11
      if (option(1:1) .eq. '?') goto 10
      call upcase (option(1:4))
c
      if (option(1:4) .eq. 'READ' .or. 
     +    option(1:4) .eq. 'NO_H' .or. 
     +    option(1:4) .eq. 'APPE' .or. 
     +    option(1:4) .eq. 'ALWY' .or. 
     +    option(1:4) .eq. 'ECHO') then
c
        lhread = (option(1:4) .ne. 'NO_H')
c
        lecho = (option(1:4) .eq. 'ECHO')
c
        write (*,*)
        call textin (' Input PDB file ?',file1)
        close (f1)
        call xopxoa (f1,file1,xinter(),ierr)
        if (ierr .ne. 0) goto 999
c
c ... find out if it's Alwyn's idea of a PDB file (i.e.,
c     chain id as a prefix to the residue nr) ...
c
        lalwyn = .false.
c
        if (option(1:4) .eq. 'ALWY') goto 310
c
  300   continue
ccc        line = ' '
        read (f1,4711,end=390,err=997) key,line
c
ccc        call upcase (key)
ccc        call upcase (line)
ccc        print *,'|',key,'|',line(1:length(line)),'|'
c
        if (key .ne. 'ATOM  ' .and.
     +      key .ne. 'HETATM') goto 300
c
        read (line,4713,err=310,end=310) atomnr(maxatm),
     +    atmnam(maxatm),altloc(maxatm),resnam(maxatm),
     +    achain(maxatm),iresid(maxatm),insert(maxatm),
     +    xf(maxatm),yf(maxatm),zf(maxatm),qatom(maxatm),
     +    batom(maxatm),inote(maxatm)
c
        goto 390
c
  310   continue
        lalwyn = .true.
        write (*,*)
     +    'OnO !!!  It is an "Alwyn-Jones-PDB-file" ...'
c
  390   continue
        rewind (f1)
        if (option(1:4) .eq. 'READ' .or.
     +      option(1:4) .eq. 'ALWY' .or.
     +      option(1:4) .eq. 'NO_H' .or.
     +      option(1:4) .eq. 'ECHO') then
          natoms = 0
          nrem = 0
          do i=1,3
            cell (i) = 1.0
            cell (i+3) = 90.0
          end do
        end if
        nlines = 0
        nhskip = 0
        spgrp = 'P 1'
        zmol = 1
c
   20   continue
          read (f1,4711,end=21,err=997) key,line
          call upcase (key)
          nlines = nlines + 1
c
          if (lecho) call textut (' >',line)
c
          if (key .eq. 'END   ') goto 20
c
          if (key(1:5) .eq. 'ORIGX') goto 20
          if (key(1:5) .eq. 'SCALE') goto 20
c
          if (key(1:5) .eq. 'CRYST') then
            read (line,*) (cell(i),i=1,6)
            spgrp = line(50:60)
            zmol = 1
            if (length(line) .ge. 64) then
              read (line(61:64),*,err=20) zmol
            end if
            goto 20
          end if
c
          if (key .eq. 'REMARK' .or.
     +        key .eq. 'HEADER' .or.
     +        key .eq. 'AUTHOR' .or.
     +        key .eq. 'COMPND' .or.
     +        key .eq. 'JRNL  ' .or.
     +        key .eq. 'SOURCE' .or.
     +        key .eq. 'REVDAT' .or.
     +        key .eq. 'SEQRES' .or.
     +        key .eq. 'HELIX ' .or.
     +        key .eq. 'SHEET ' .or.
     +        key .eq. 'SSBOND' .or.
     +        key .eq. 'FTNOTE' .or.
     +        key .eq. 'MTRIX1' .or.
     +        key .eq. 'MTRIX2' .or.
     +        key .eq. 'MTRIX3' .or.
     +        key .eq. 'HET   ' .or.
     +        key .eq. 'FORMUL' ) then
            if (nrem .lt. maxcom) then
              nrem = nrem + 1
              remark (nrem) = key//line
              if (nrem .eq. maxcom) call errcon (
     +          'WARNING - max nr of non-ATOM records read')
            end if
            goto 20
          end if
c
          if (key .ne. 'ATOM  ' .and.
     +        key .ne. 'HETATM') goto 20
c
          natoms = natoms + 1
          call upcase (line)
c
          if (natoms .gt. maxatm) then
            write (*,6666) 'Too many atoms -- rest skipped'
            call jvalut (' Max nr of atoms =',1,maxatm)
            natoms = natoms - 1
            goto 21
          end if
c
          kk = natoms
          if (lalwyn) then
            read (line,4715,err=5372) atomnr(kk),atmnam(kk),altloc(kk),
     +        resnam(kk),aresid(kk),insert(kk),xf(kk),
     +        yf(kk),zf(kk),qatom(kk),batom(kk),inote(kk)
            call detaj (aresid(kk),achain(kk),iresid(kk))
          else
            read (line,4713,err=5372) atomnr(kk),atmnam(kk),altloc(kk),
     +        resnam(kk),achain(kk),iresid(kk),insert(kk),xf(kk),
     +        yf(kk),zf(kk),qatom(kk),batom(kk),inote(kk)
          end if
c
          if (.not. lhread) then
            if (lhydro(atmnam(kk))) then
              natoms = natoms - 1
              nhskip = nhskip + 1
            end if
          end if
c
        goto 20
c
c ... read from string error
c
 5372   continue
        call errcon ('While reading PDB file')
        call jvalut (' Line nr :',1,nlines)
        call jvalut (' Atom nr :',1,natoms)
        call textut (' Line :',line)
        goto 11
c
   21   continue
        call jvalut (' Number of lines read :',1,nlines)
        if (.not. lhread) then
          call jvalut (' Hydrogens skipped    :',1,nhskip)
        end if
        call jvalut (' Number of atoms now  :',1,natoms)
        newchn = achain(1)
        goto 11
c
      else if (option(1:4) .eq. 'QUIT') then
c
        goto 9999
c
      else if (option(1:4) .eq. 'HELI' .or.
     +         option(1:4) .eq. 'STRA') then
c
        if (natoms .gt. 0) then
          call prompt (
     +      ' WARNING - This will ERASE your current molecule !')
          answer = 'N'
          call textin (' Are you sure ?',answer)
          call upcase (answer)
          if (answer .ne. 'Y') goto 11
        end if
c
        call prompt (' Provide the coordinates of the N- and the C-')
        call prompt (' terminal CA atoms (approximately).')
        call prompt (' If you only want to generate X residues,')
        call prompt (' enter 0,0,0 for the N-terminus and:')
        call prompt (' HELIX : 1.46*X,0,0 for the C-terminus,')
        call prompt (' STRAND: 3.32*X,0,0 for the C-terminus')
        call fvalin (' Start coordinates (N-term) ?',3,start)
        call fvalin (' End   coordinates (C-term) ?',3,end)
c
        if (option(1:4) .eq. 'HELI') then
          call helstr ('A',start,end,natoms,atomnr,atmnam,altloc,
     +      resnam,achain,iresid,insert,xf,yf,zf,allxyz,qatom,
     +      batom,inote,xxyyzz,ibuff,maxatm)
        else
          call helstr ('B',start,end,natoms,atomnr,atmnam,altloc,
     +      resnam,achain,iresid,insert,xf,yf,zf,allxyz,qatom,
     +      batom,inote,xxyyzz,ibuff,maxatm)
        end if
c
        call prompt (' Use Move_zone and RSR_rigid in O to')
        call prompt (' position the structure element into the')
        call prompt (' density; then use Mutate_* and Lego_side_ch')
        call prompt (' for the correct residue types and rotamers.')
c
        goto 11
c
      endif
c
      if (natoms .le. 0 .and. option(1:4) .ne. 'IMPO') then
        write (*,6666) 'No molecule read yet'
        goto 11
      end if
c
      if (option(1:4) .eq. 'TEMP') then
c
        call fvalin (' TFnew = A * TFold + B.  Enter A, B :',2,bfacts)
        call ivalin (
     +    ' Residue range to apply (0 0 = all molecule) ?',2,izone)
c
        ldum = 0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 2011
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 2011
          batom (i) = bfacts(1)*batom(i) + bfacts(2)
          ldum = ldum + 1
 2011     continue
        end do
        call jvalut (' Nr of temperature factors updated :',1,ldum)
c
      else if (option(1:4) .eq. 'REMA') then
c
        if (nrem .gt. 0) then
          do i=1,nrem
            write (*,'(a,i3,a,a)') ' #',i,' > ',
     +        (remark(i)(1:leng1(remark(i))))
          end do
        else
          call prompt (' No REMARK etc. cards in memory')
        end if
c
c ... Gerard's TEST option (least-cluttered view)
c
      else if (option(1:4) .eq. 'TEST') then
c
        call minclu (natoms,xf,yf,zf,allxyz)
c
      else if (option(1:4) .eq. 'CRYS') then
c
        call fvalin (' Unit-cell constants ?',6,cell)
        call celvol (cell,volume)
        call rvalut (' Unit-cell volume (A3) :',1,volume)
        if (length(spgrp) .eq. 0) spgrp = 'P 1'
        call upcase (spgrp)
        call textin (' Spacegroup ?',spgrp)
        call upcase (spgrp)
        zmol = max (1,zmol)
        call ivalin (' Value of Z ?',1,zmol)
        zmol = max (1,zmol)
c
      else if (option(1:4) .eq. 'PIR_') then
c
        call textin (' Amino acid residue names ?',amacid)
        call pretty (amacid)
        call upcase (amacid)
        call textin (' One letter codes         ?',onelet)
        call pretty (onelet)
        call upcase (onelet)
c
        call textin (' PIR output file name     ?',file8)
        close (f8)
        call xopxua (f8,file8,xinter(),ierr)
        if (ierr .ne. 0) goto 11
c
        call textin (' PIR sequence name        ?',pirnam)
        call textin (' PIR sequence title       ?',pirtit)
c
        call stamp (line)
        write (f8,'(a1,1x,a)') '!',line(1:leng1(line))
        write (f8,*)
        write (f8,'(a4,a)') '>P1;',pirnam(1:leng1(pirnam))
        write (f8,'(a)') pirtit(1:leng1(pirtit))
c
        nres = 0
        nowres = -1
        ndum = 0
c
        do i=1,natoms
c
c ... a new residue ?
c
          if (iresid(i) .ne. nowres) then
            nres = nres + 1
            j = index (amacid,resnam(i))
            if (j .le. 0) then
              pircod (nres) = ' '
            else
              j = 1 + 2*(j/4)
              pircod (nres) = onelet(j:j)
              ndum = ndum + 1
            end if
            nowres = iresid (i)
          end if
c
        end do
c
        write (f8,'(6(10a1,1x))') (pircod(i),i=1,nres)
        write (f8,'(a1)') '*'
c
        close (f8)
        call jvalut (' Number of residues encountered :',1,nres)
        call jvalut (' Recognised one-letter codes    :',1,ndum)
c
      else if (option(1:4) .eq. 'TALL') then
c
        nres = 0
        nowres = -1
c
        do i=1,maxtyp
          typcnt (i) = 0
          typnam (i) = '   '
        end do
        ntyp = 1
        typnam (1) = 'H2O'
c
        do i=1,natoms
          if (iresid(i) .ne. nowres) then
            nres = nres + 1
            nowres = iresid (i)
            if (lwater(resnam(i))) then
              typcnt (1) = typcnt (1) + 1
              goto 2052
            end if
            do j=1,ntyp
              if (resnam(i) .eq. typnam(j)) then
                typcnt (j) = typcnt (j) + 1
                goto 2052
              end if
            end do
            ntyp = ntyp + 1
            typnam (ntyp) = resnam(i)
            typcnt (ntyp) = 1
 2052       continue
          end if
        end do
c
        if (ntyp .gt. 1) then
          do i=2,ntyp
            dummy = 100.0 * float(typcnt(i)) / float(nres)
            write (*,2357) typnam(i),typcnt(i),dummy
          end do
        end if
c
 2357 format (' Type ',a3,' Number = ',i6,' ( ',f6.2,' %)')
c
        if (typcnt(1) .gt. 0) then
          call jvalut (' Nr of waters found   :',1,typcnt(1))
          call prompt (' (Using the most usual names for waters)')
        else
          call prompt (' There are *NO* waters')
        end if
c
        call jvalut (' Total nr of residues :',1,nres)
c
      else if (option(1:4) .eq. 'GLYC') then
c
        nres = 0
        nowres = -1
c
        do i=1,natoms
          if (iresid(i) .ne. nowres) then
            nres = nres + 1
            resptr (nres) = i
            nowres = iresid (i)
          end if
        end do
c
        call jvalut (' Nr of residues found :',1,nres)
        if (nres .lt. 3) goto 11
c
        ndum = 0
        do i=2,nres-2
          j = resptr (i)
          j2 = resptr (i+2)
          j1 = resptr (i+1)
          if (iresid(j2) .eq. (iresid(j)+2) .and.
     +        resnam(j) .eq. 'ASN') then
            if (resnam(j2) .eq. 'SER' .or.
     +          resnam(j2) .eq. 'THR') then
              if (resnam(j1) .eq. 'PRO') then
                write (*,6514) resnam(j),iresid(j),
     +            resnam(j2),iresid(j2)
              else
                ndum = ndum + 1
                write (*,6515) ndum,resnam(j),iresid(j),
     +            resnam(j2),iresid(j2)
              end if
            end if
          end if
        end do
c
 6514 format (' *NOT*  an N-glycosylation site       : ',
     +        a3,i5,' -    PRO    - ',a3,i5)
 6515 format (' Potential N-glycosylation site # ',i3,' : ',
     +        a3,i5,' - (not PRO) - ',a3,i5)
c
        call jvalut (' Potential glycosylation sites :',1,ndum)
c
      else if (option(1:4) .eq. 'EXTI') then
c
        call extinc (natoms,atmnam,resnam)
c
      else if (option(1:4) .eq. 'COUN') then
c
        call counte (natoms,atmnam)
c
      else if (option(1:4) .eq. 'SMOO') then
c
        call fvalin (' Cut-off distance for bonded atoms ?',1,bndcut)
        call cpcoor (natoms,xf,yf,zf,xxyyzz)
        call smoobs (natoms,bndcut,batom,b1,b2,xxyyzz,atmnam)
c
      else if (option(1:4) .eq. 'BOND') then
c
        call fvalin (' Cut-off distance for bonded atoms ?',1,bndcut)
        call cpcoor (natoms,xf,yf,zf,xxyyzz)
        call bondbs (natoms,bndcut,batom,b1,b2,xxyyzz,atmnam)
c
      else if (option(1:4) .eq. 'SHOR') then
c
        call fvalin (' Cut-off distance for SHORT contacts ?',1,badcut)
        call cpcoor (natoms,xf,yf,zf,xxyyzz)
        call shorty (natoms,badcut,xxyyzz,
     +    atmnam,resnam,iresid,achain)
c
      else if (option(1:4) .eq. 'NONB') then
c
        write (*,'(99(1x,a:/))')
     +    'RMS delta-B of non-bonded atoms',
     +    'Select one of the following options :',
     +    ' 1 = all atoms',
     +    ' 2 = N and O only (H-bonds & salt links)',
     +    ' 3 = all atoms, except N and O'
        call ivalin (' Option (1-3) ?',1,inb)
        inb = max (1,min(3,inb))
c
        call fvalin (' Cut-off distance for bonded atoms ?',1,bndcut)
        call fvalin (' Cut-off for NON-bonded atoms      ?',1,nobcut)
        call cpcoor (natoms,xf,yf,zf,xxyyzz)
        call nonbbs (natoms,bndcut,nobcut,batom,b1,b2,xxyyzz,
     +    atmnam,inb,baveok)
c
      else if (option(1:4) .eq. 'B_Q_') then
c
        call textin (' Amino acid residue names     ?',amacid)
        call pretty (amacid)
        call upcase (amacid)
c
c        call textin (' Residue name for waters      ?',awater)
c        call remspa (awater)
c        call upcase (awater)
c
        call textin (' Names of ligands/substrates  ?',ligand)
        call pretty (ligand)
        call upcase (ligand)
c
        call textin (' Which chain (** = all)       ?',whichn)
        call upcase (whichn)
c
        call textin (
     +    ' Include HYDROGEN atoms (Y/N) ?',keephs)
        call upcase (keephs)
        lkeeph = (keephs .eq. 'Y')
c
        do i=1,7
          bstats (i) =  0.0
          bmax   (i) = -1.0
          qstats (i) =  0.0
          bqcnts (i) =  0
        end do
c
        do i=1,natoms
c
          if (.not. lkeeph) then
            if (lhydro(atmnam(i))) goto 6592
          end if
c
          if (whichn .ne. '**') then
            if (achain (i) .ne. whichn) goto 6592
          end if
c
          bqcnts (5) = bqcnts (5) + 1
          bstats (5) = bstats (5) + batom(i)
          qstats (5) = qstats (5) + qatom(i)
          bmax   (5) = max (bmax(5), batom(i))
c
          if (lwater(resnam(i))) then
            bqcnts (4) = bqcnts (4) + 1
            bstats (4) = bstats (4) + batom(i)
            qstats (4) = qstats (4) + qatom(i)
            bmax   (4) = max (bmax(4), batom(i))
          else if (index(ligand,resnam(i)) .gt. 0) then
            bqcnts (7) = bqcnts (7) + 1
            bstats (7) = bstats (7) + batom(i)
            qstats (7) = qstats (7) + qatom(i)
            bmax   (7) = max (bmax(7), batom(i))
          else if (index(amacid,resnam(i)) .gt. 0) then
            bqcnts (6) = bqcnts (6) + 1
            bstats (6) = bstats (6) + batom(i)
            qstats (6) = qstats (6) + qatom(i)
            bmax   (6) = max (bmax(6), batom(i))
            if (mainch(atmnam(i))) then
              bqcnts (1) = bqcnts (1) + 1
              bstats (1) = bstats (1) + batom(i)
              qstats (1) = qstats (1) + qatom(i)
              bmax   (1) = max (bmax(1), batom(i))
            else
              bqcnts (2) = bqcnts (2) + 1
              bstats (2) = bstats (2) + batom(i)
              qstats (2) = qstats (2) + qatom(i)
              bmax   (2) = max (bmax(2), batom(i))
            end if
          else
            bqcnts (3) = bqcnts (3) + 1
            bstats (3) = bstats (3) + batom(i)
            qstats (3) = qstats (3) + qatom(i)
            bmax   (3) = max (bmax(3), batom(i))
          end if
 6592     continue
        end do
c
        do i=1,7
          if (bqcnts(i) .gt. 0) then
            bstats(i)=bstats(i)/float(bqcnts(i))
            qstats(i)=qstats(i)/float(bqcnts(i))
          else
            bstats (i) = 0.0
            qstats (i) = 0.0
            bqcnts (i) = 0
            bmax   (i) = 0.0
          end if
        end do
c
        write (*,*)
        call textut (' B & Q statistics for chain :',whichn)
        write (*,6752) 'Atom type         ',
     +    'Number','Average B','Maximum B','Average Q'
        write (*,6753) 'Protein main chain',
     +    bqcnts(1),bstats(1),bmax(1),qstats(1)
        write (*,6753) 'Protein side chain',
     +    bqcnts(2),bstats(2),bmax(2),qstats(2)
        write (*,6753) 'Protein all atoms ',
     +    bqcnts(6),bstats(6),bmax(6),qstats(6)
        write (*,6753) 'Ligand/substrate  ',
     +    bqcnts(7),bstats(7),bmax(7),qstats(7)
        write (*,6753) 'Water molecules   ',
     +    bqcnts(4),bstats(4),bmax(4),qstats(4)
        write (*,6753) 'Other entities    ',
     +    bqcnts(3),bstats(3),bmax(3),qstats(3)
        write (*,6753) 'All atoms         ',
     +    bqcnts(5),bstats(5),bmax(5),qstats(5)
c
        if (bqcnts(5) .gt. 0) then
          call r5bq (whichn,lkeeph,bqcnts,bstats,bmax)
        end if
c
        if (lkeeph) call prompt (' INCLUDES HYDROGEN ATOMS !!!')
        write (*,*)
c
 6752 format (/1x,a,1x,a6,3(1x,a10))
 6753 format (1x,a,1x,i6,3(1x,f10.3))
c
      else if (option(1:4) .eq. 'SEQU') then
c
        jres = ires
        ires = 0
        nwrit = 0
        do i=1,natoms
          if (iresid(i) .ne. ires) then
            ires = iresid(i)
            write (line,4713) i,atmnam(i),altloc(i),
     +        resnam(i),achain(i),iresid(i),insert(i),xf(i),
     +        yf(i),zf(i),qatom(i),batom(i),inote(i)
            write (*,'(a)') line(1:leng1(line))
            nwrit = nwrit + 1
          end if
        end do
c
        ires = jres
        call jvalut (' Nr of residues found :',1,nwrit)
c
      else if (option(1:4) .eq. 'LIST') then
c
        call ivalin (' Which residue number ?',1,ires)
        if (ires .le. 0) then
          call errcon ('Residue number must be positive !')
          goto 11
        end if
c
        nwrit = 0
        do i=1,natoms
          if (iresid(i) .eq. ires) then
            write (line,4713) i,atmnam(i),altloc(i),
     +        resnam(i),achain(i),iresid(i),insert(i),xf(i),
     +        yf(i),zf(i),qatom(i),batom(i),inote(i)
            write (*,'(a)') line(1:leng1(line))
            nwrit = nwrit + 1
          end if
        end do
c
        if (nwrit .gt. 0) then
          call ivalut (' Nr of atoms listed :',1,nwrit)
        else
          call errcon (' Residue does not exist')
        end if
        write (*,*)
c
      else if (option(1:4) .eq. 'RADI') then
c
        line = 'b_radial.plt'
        call textin (' Plot file ?',line)
        if (length(line).le.0) goto 11
c
        newchn = ' A'
        call textin (' Which chain (2 characters !) ?',newchn)
c
        nwrit = 0
        do i=1,natoms
          if (.not. lhydro(atmnam(i))) then
            if (achain(i) .eq. newchn) then
              nwrit = nwrit + 1
              allxyz(1,nwrit)=xf(i)
              allxyz(2,nwrit)=yf(i)
              allxyz(3,nwrit)=zf(i)
              xbuff(nwrit)=batom(i)
            end if
          end if
        end do
        call jvalut (' Nr of atoms selected (no Hs) :',1,nwrit)
        if (nwrit .le. 3) goto 11
c
        call bradpl (f1,line,nwrit,allxyz,xbuff)
c
      else if (option(1:4) .eq. 'EXPO') then
c
        call badman (
     +    option,f7,file7,natoms,maxapr,maxbuf,iresid,atmnam,
     +    iot,names,sorted,allxyz,xf,yf,zf,altloc,resnam,achain,
     +    insert,qatom,batom,ibuff)
c
      else if (option(1:4) .eq. 'SAME') then
c
        call badman (
     +    option,f7,file7,natoms,maxapr,maxbuf,iresid,atmnam,
     +    iot,names,sorted,allxyz,xf,yf,zf,altloc,resnam,achain,
     +    insert,qatom,batom,ibuff)
c
      else if (option(1:4) .eq. 'IMPO') then
c
        write (*,*)
        call textin (
     +    ' Input Bonds/Angles/Dihedrals file ?',file7)
        close (f7)
        call xopxoa (f7,file7,xinter(),ierr)
        if (ierr .ne. 0) goto 899
c
        natoms = 0
        nlines = 0
        nrem   = 0
c
  920   continue
          read (f7,4711,end=921,err=897) key,line
          call upcase (key)
          nlines = nlines + 1
c
          if (key .ne. 'BAD   ') goto 920
c
          natoms = natoms + 1
          call upcase (line)
c
          if (natoms .gt. maxatm) then
            write (*,6666) 'Too many atoms -- rest skipped'
            call jvalut (' Max nr of atoms =',1,maxatm)
            natoms = natoms - 1
            goto 921
          end if
c
          kk = natoms
          iref = 3*(kk-1)
          jref = 3*(maxatm+kk-1)
          read (line,4713) atomnr(kk),atmnam(kk),altloc(kk),
     +      resnam(kk),achain(kk),iresid(kk),insert(kk),
     +      xbuff(iref+1),xbuff(iref+2),xbuff(iref+3),
     +      qatom(kk),batom(kk),xnote
          read (xnote,*) ibuff(jref+1),ibuff(jref+2),ibuff(jref+3)
          inote (kk) = ' '
c
        goto 920
c
  921   continue
c
        call jvalut (' Nr of lines read :',1,nlines)
        call jvalut (' Nr of atoms read :',1,natoms)
        write (*,*) '... Converting to Cartesian coordinates ...'
c
        call c2cart (natoms,allxyz,xbuff(1),ibuff(3*maxatm+1))
        do i=1,natoms
          xf (i) = allxyz (1,i)
          yf (i) = allxyz (2,i)
          zf (i) = allxyz (3,i)
        end do
c
      else if (option(1:4) .eq. 'PLAN') then
c
        write (*,*)
        write (*,*)
     +    'Checking peptide planarity ...'
        write (*,*)
     +    'Improper Ci-CAi-Ni+1-Oi should be ZERO degrees'
        write (*,*)
     +    'The nr of * indicates deviation from planarity'
c
        write (*,*)
        do i=1,8
          iot (i) = 0
        end do
        napr = 0
        iold = 0
c
        do k=1,natoms
          if (atmnam(k) .eq. ' N  ') then
            napr = napr + 1
            iot (5) = k
            otxyz (1,5) = xf (k)
            otxyz (2,5) = yf (k)
            otxyz (3,5) = zf (k)
          else if (atmnam(k) .eq. ' CA ') then
            napr = napr + 1
            iot (6) = k
            otxyz (1,6) = xf (k)
            otxyz (2,6) = yf (k)
            otxyz (3,6) = zf (k)
          else if (atmnam(k) .eq. ' C  ') then
            napr = napr + 1
            iot (7) = k
            otxyz (1,7) = xf (k)
            otxyz (2,7) = yf (k)
            otxyz (3,7) = zf (k)
          else if (atmnam(k) .eq. ' O  ') then
            napr = napr + 1
            iot (8) = k
            otxyz (1,8) = xf (k)
            otxyz (2,8) = yf (k)
            otxyz (3,8) = zf (k)
          else
            goto 8270
          end if
c
          if (iot(5) .gt. 0 .and. iot(6) .gt. 0 .and.
     +        iot(7) .gt. 0 .and. iot(8) .gt. 0) then
c
            if (napr .gt. 4) then
              if (dist(2,6,otxyz) .ge. 6.0) goto 8253
              iold = iold + 1
              i = iot (4)
              dummy = tangle(3,2,5,4,otxyz)
              write (line,4713) i,atmnam(i),altloc(i),
     +          resnam(i),achain(i),iresid(i),insert(i),
     +          dummy
c
              if (abs(dummy) .gt. 0.001) then
                i = min (25, int(abs(dummy/0.5)))
                if (i.gt.0) then
                  write (stars,'(1x,79a1)') ('*',j=1,i)
                  call appstr (line,stars)
                end if
              end if
c
              dummy = abs (tangle (2,3,5,6,otxyz))
              if (dummy. le. 30.0) then
                write (*,'(a,f8.1)')
     +   ' *** Next residue has a CIS-peptide !  OMEGA = ',dummy
              else if (dummy .lt. 150.0) then
                write (*,'(a,f8.1)')
     +   ' *** Next residue is distorted      !  OMEGA = ',dummy
              end if
c
              write (*,'(a)') line(1:leng1(line))
 8253         continue
            end if
c
            do i=1,4
              do j=1,3
                otxyz (j,i) = otxyz (j,i+4)
              end do
              iot (i)   = iot (i+4)
              iot (i+4) = 0
            end do
c
          end if
c
 8270     continue
        end do
c
        write (*,*)
        call jvalut (' Nr of peptide planes :',1,iold)
        write (*,*)
c
      else if (option(1:4) .eq. 'FIT_') then
c
        write (*,*)
        line = 'fit_waters.omac'
        call textin (' Name of water-fitting O macro ?',line)
        close (f9)
        call xopxua (f9,line,xinter(),ierr)
        if (ierr .ne. 0) goto 11
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
        write (f9,5110) '# RS-fit RFAC or RSCC ? #'
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
        do i=1,natoms
          if (.not. lhydro(atmnam(i))) then
            if (lwater(resnam(i))) then
              ns = ns + 1
              write (s1name,'(a2,i4,a1)') achain(i),iresid(i),
     +          insert(i)
              call remspa (s1name)
              if (ns .eq. 1) s2name = s1name
              j = length(s1name)
              write (f9,5110) '!'
              write (f9,5110) 'rs_fit',s1name(1:j),';'
              write (f9,5110) 'rsr_rigid',s1name(1:j),'; yes'
            end if
          end if
        end do
c
        if (ns .gt. 0) then
          write (f9,5110) '!'
          write (f9,5110) 'Message Calculating new RS-fit values'
          write (f9,5110) '! copy_db ${mymol}//_residue_rsprefit',
     +      '${mymol}//_residue_rsfit'
          write (f9,5110) 'rs_fit',s2name,s1name
        end if
        write (f9,5110) '!'
        write (f9,5110) 'bell Message Done'
c
        close (f9)
        call prompt (' O macro written')
c
        call jvalut (' Nr of waters to be fitted :',1,ns)
c
 5110   format (20(a,1x))
c
      else if (option(1:4) .eq. 'SSBO') then
c
        call chkrem (ldorem,'SSBOND ')
        if (ldorem) then
          call errcon ('SSBOND records already present')
          goto 11
        end if
c
        call fvalin (' Max SS bond length (A) ?',1,ssdist)
        write (*,*)
c
        ns = 0
        do i=1,natoms
          if (resnam(i) .eq. 'CYS') then
            if (atmnam(i) .eq. ' SG ') then
              write (line,4713) i,atmnam(i),altloc(i),
     +          resnam(i),achain(i),iresid(i),insert(i),
     +          xf(i),yf(i),zf(i),qatom(i),batom(i),inote(i)
              write (*,'(a)') line(1:leng1(line))
              ns = ns + 1
              ibuff (ns) = i
              xxyyzz (1,ns) = xf (i)
              xxyyzz (2,ns) = yf (i)
              xxyyzz (3,ns) = zf (i)
            end if
          end if
        end do
c
        write (*,*)
        call jvalut (' Nr of CYS SG atoms found :',1,ns)
        if (ns .le. 1) goto 11
c
        write (*,*)
c
        nss = 0
        do i=1,ns-1
          i1 = ibuff (i)
          do j=i+1,ns
            j1 = ibuff (j)
            xx = dist (i,j,xxyyzz)
            if (xx .le. ssdist) then
              nss = nss + 1
              if (nrem .lt. maxcom) then
                nrem = nrem + 1
                write (remark(nrem),
     + '(a6,1x,i3,1x,a3,a2,1x,i4,a1,3x,a3,a2,1x,i4,a1,4x,a,f6.2,a)')
     +            'SSBOND',nss,'CYS',achain(i1),iresid(i1),
     +            insert(i1),'CYS',achain(j1),iresid(j1),
     +            insert(j1),'S-S = ',xx,' A'
              else
                call errcon ('Too many remark etc. records')
              end if
              write (*,'(1x,a)')
     +          remark(nrem)(1:leng1(remark(nrem)))
            end if
          end do
        end do
c
        write (*,*)
        call jvalut (' Nr of disulfide bridges :',1,nss)
c
      else if (option(1:4) .eq. 'DISU') then
c
        write (*,*)
        line = 'ss.odl'
        call textin (' Name of disulfide ODL file ?',line)
        call xopxua (f9,line,xinter(),ierr)
        if (ierr .ne. 0) goto 11
c
        call textin (' Name of the molecule in O ?',omol)
        call fvalin (' Stick radius (A) ?',1,strad)
        call fvalin (' Max SS bond length (A) ?',1,ssdist)
        write (*,*)
c
        ns = 0
        do i=1,natoms
          if (resnam(i) .eq. 'CYS') then
            if (atmnam(i) .eq. ' SG ') then
              write (line,4713) i,atmnam(i),altloc(i),
     +          resnam(i),achain(i),iresid(i),insert(i),
     +          xf(i),yf(i),zf(i),qatom(i),batom(i),inote(i)
              write (*,'(a)') line(1:leng1(line))
              ns = ns + 1
              ibuff (ns) = i
              xxyyzz (1,ns) = xf (i)
              xxyyzz (2,ns) = yf (i)
              xxyyzz (3,ns) = zf (i)
            end if
          end if
        end do
c
        write (*,*)
        call jvalut (' Nr of CYS SG atoms found :',1,ns)
        if (ns .le. 1) goto 11
c
        write (*,*)
        write (f9,6600) 'begin ssbond'
        write (f9,6600) '  colour green'
c
        nss = 0
        do i=1,ns-1
          i1 = ibuff (i)
          write (s1name,'(a2,i4,a1)') achain(i1),iresid(i1),
     +      insert(i1)
          call remspa (s1name)
          do j=i+1,ns
            j1 = ibuff (j)
            xx = dist (i,j,xxyyzz)
            if (xx .le. ssdist) then
              write (s2name,'(a2,i4,a1)') achain(j1),iresid(j1),
     +          insert(j1)
              call remspa (s2name)
              write (line,6610) 'stick',omol,s1name,atmnam(i1),
     +          omol,s2name,atmnam(j1),strad
              call pretty (line)
              call locase (line)
              write (f9,6600) ('  '//line(1:leng1(line)))
              nss = nss + 1
              call textut ('S=S bond :',line)
              call fvalut (' SG-SG distance (A) :',1,xx)
            end if
          end do
        end do
c
        write (*,*)
        call jvalut (' Nr of disulfide bridges :',1,nss)
        write (f9,6600) 'end_object'
        close (f9)
c
 6600 format (a)
 6610 format (7(a8,1x),f8.2)
c
      else if (option(1:4) .eq. 'RSFI' .or.
     +         option(1:4) .eq. 'RSR_' .or.
     +         option(1:4) .eq. 'CONN' .or.
     +         option(1:4) .eq. 'TORS') then
c
        call odicts (
     +    option,natoms,ires,bndcut,iot,iresid,atmnam,xf,yf,zf,
     +    altloc,resnam,achain,insert,qatom,batom,inote,otxyz,
     +    f1,isbond,okbond,dismat,deftor,valtor,afftor,tornam)
c
      else if (option(1:4) .eq. 'GEOM') then
c
        write (*,*)
        write (*,*) 'Select an option:'
        write (*,*) '1 - analyse one residue'
        write (*,*) '2 - analyse backbone (N-CA-C)'
        write (*,*) '3 - analyse CAs'
        write (*,*) '4 - analyse two residues (eg, disulfide)'
        write (*,*) '5 - analyse zone of residues'
        write (*,*) '6 - analyse disulfide'
        call ivalin (' Option (1-6; 0 to abort) ?',1,igeo)
        if (igeo .lt. 1 .or. igeo .gt. 6) goto 11
c
        if (igeo.eq.1 .or. igeo.eq.4 .or.
     +      igeo.eq.5 .or. igeo.eq.6) then
c
c ... analyse one residue
c
          if (igeo.eq.1) then
c
            call ivalin (' Which residue number ?',1,ires)
            if (ires .le. 0) then
              call errcon ('Residue number must be positive !')
              goto 11
            end if
c
c ... analyse two residues
c
          else if (igeo.eq.4 .or. igeo .eq.6) then
c
            if (igeo .eq. 6) then
              ires = 0
              nwrit = 0
              do i=1,natoms
                if (iresid(i) .ne. ires) then
                  if (resnam(i) .eq. 'CYS' .or.
     +                resnam(i) .eq. 'CYH' ) then
                    jres = ires
                    ires = iresid(i)
                    write (line,4713) i,atmnam(i),altloc(i),
     +                resnam(i),achain(i),iresid(i),insert(i),
     +                xf(i),yf(i),zf(i),qatom(i),batom(i),inote(i)
                    write (*,'(a)') line(1:leng1(line))
                    nwrit = nwrit + 1
                  end if
                end if
              end do
              call ivalut (' Nr of CYS/CYH residues :',1,nwrit)
              if (nwrit .lt. 2) then
                call errcon ('No disulphides')
                goto 11
              end if
              nwrit = ires
              ires = jres
              jres =nwrit
            end if
c
            call ivalin (' First  residue number ?',1,ires)
            if (ires .le. 0) then
              call errcon ('Residue number must be positive !')
              goto 11
            end if
c
            call ivalin (' Second residue number ?',1,jres)
            if (jres .le. 0) then
              call errcon ('Residue number must be positive !')
              goto 11
            end if
c
c ... analyse zone of residues
c
          else if (igeo.eq.5) then
c
            call ivalin (' First residue number ?',1,ires)
            if (ires .le. 0) then
              call errcon ('Residue number must be positive !')
              goto 11
            end if
c
            call ivalin (' Last  residue number ?',1,jres)
            if (jres .le. 0) then
              call errcon ('Residue number must be positive !')
              goto 11
            end if
c
            call ilohi (ires,jres)
c
          end if
c
          if (igeo .eq. 6) then
             n1 = 0
             n2 = 0
             ca1 = 0
             ca2 = 0
             cb1 = 0
             cb2 = 0
             sg1 = 0
             sg2 = 0
          else
            call fvalin (' Cut-off distance for bonded atoms ?',
     +        1,bndcut)
            bndcut = max (0.1, bndcut)
          end if
c
          write (*,*)
          napr = 0
          do i=1,natoms
            if ( (igeo.eq.1 .and. iresid(i) .eq. ires) .or.
     +           (igeo.eq.4 .and. (iresid(i) .eq. ires .or.
     +                             iresid(i) .eq. jres)) .or.
     +           (igeo.eq.6 .and. (iresid(i) .eq. ires .or.
     +                             iresid(i) .eq. jres)) .or.
     +           (igeo.eq.5 .and. (iresid(i) .ge. ires .and.
     +                             iresid(i) .le. jres)) ) then
              if (napr .ge. maxapr) then
                call errcon ('Too many atoms')
                call ivalut (' Maximum :',1,maxapr)
                goto 8362
              end if
c
              if (igeo .eq. 6) then
                if (atmnam(i) .eq. ' N  ') then
                  if (iresid(i) .eq. ires) then
                    n1 = napr + 1
                  else
                    n2 = napr + 1
                  end if
                else if (atmnam(i) .eq. ' CA ') then
                  if (iresid(i) .eq. ires) then
                    ca1 = napr + 1
                  else
                    ca2 = napr + 1
                  end if
                else if (atmnam(i) .eq. ' CB ') then
                  if (iresid(i) .eq. ires) then
                    cb1 = napr + 1
                  else
                    cb2 = napr + 1
                  end if
                else if (atmnam(i) .eq. ' SG ') then
                  if (iresid(i) .eq. ires) then
                    sg1 = napr + 1
                  else
                    sg2 = napr + 1
                  end if
                else
                  goto 6501
                end if
              end if
c
              napr = napr + 1
              iot (napr) = i
              otxyz (1,napr) = xf (i)
              otxyz (2,napr) = yf (i)
              otxyz (3,napr) = zf (i)
              write (line,4713) i,atmnam(i),altloc(i),
     +          resnam(i),achain(i),iresid(i),insert(i),xf(i),
     +          yf(i),zf(i),qatom(i),batom(i),inote(i)
              write (*,'(a)') line(1:leng1(line))
c
 6501         continue
            end if
          end do              
c
 8362     continue
          write (*,*)
          call ivalut (' Nr of atoms found :',1,napr)
          if (napr .le. 0) then
            call errcon ('No atoms means nothing to do')
            goto 11
          end if
c
          if (igeo .eq. 6) then
c
            if (n1 .lt. 1 .or. n2 .lt. 1 .or.
     +          ca1 .lt. 1 .or. ca2 .lt. 1 .or.
     +          cb1 .lt. 1 .or. cb2 .lt. 1 .or.
     +          sg1 .lt. 1 .or. sg2 .lt. 1) then
              call errcon ('Missing atoms')
              goto 11
            end if
c
            if (dist(sg1,sg2,otxyz) .gt. 3.0) call prompt (
     +        ' WARNING >>> S-S distance > 3.0 A !')
c
 6512 format (5(1x,a5))
 6513 format (2i6,3f6.2)
            write (*,*)
            write (*,6512) '  1  ','  2  ','SG-SG',
     +        'CB-CB','CA-CA'
            write (*,6512) ('=====',i=1,5)
            write (*,6513) ires,jres,dist(sg1,sg2,otxyz),
     +        dist(cb1,cb2,otxyz),dist(ca1,ca2,otxyz)
c
 6522 format (2(1x,a5),6(1x,a11))
 6523 format (2i6,6(3x,f6.1,3x))
            write (*,*)
            write (*,6522) '  1  ','  2  ','SG1-SG2-CB2',
     +        'SG2-SG1-CB1','SG1-CB1-CA1','SG2-CB2-CA2',
     +        'CB1-CA1-N1','CB2-CA2-N2'
            write (*,6522) '=====','=====',
     +        ('===========',i=1,6)
            write (*,6523) ires,jres,angle(sg1,sg2,cb2,otxyz),
     +        angle(sg2,sg1,cb1,otxyz),angle(sg1,cb1,ca1,otxyz),
     +        angle(sg2,cb2,ca2,otxyz),angle(cb1,ca1,n1,otxyz),
     +        angle(cb2,ca2,n2,otxyz)
c
 6532 format (2(1x,a5),5(1x,a15))
 6533 format (2i6,5(5x,f6.1,5x))
            write (*,*)
            write (*,6532) '  1  ','  2  ','N1-CA1-CB1-SG1',
     +        'CA1-CB1-SG1-SG2','CB1-SG1-SG2-CB2',
     +        'N2-CA2-CB2-SG2','CA2-CB2-SG2-SG1'
            write (*,6532) '=====','=====',
     +        ('===============',i=1,5)
            write (*,6533) ires,jres,
     +        tangle(n1 ,ca1,cb1,sg1,otxyz),
     +        tangle(ca1,cb1,sg1,sg2,otxyz),
     +        tangle(cb1,sg1,sg2,cb2,otxyz),
     +        tangle(n2 ,ca2,cb2,sg2,otxyz),
     +        tangle(ca2,cb2,sg2,sg1,otxyz)
c
            write (*,*)
            goto 11
          end if
c
          write (*,'(/a)') ' *** Distances ***'
          do i=1,napr
            do j=i+1,napr
              dismat (i,j) = dist (i,j,otxyz)
              dismat (j,i) = dismat (i,j)
            end do
            dismat (i,i) = 0.0
          end do
c
          answer = 'N'
          call textin (
     +      ' Do you want to see the complete distance matrix ?',
     +      answer)
          call upcase (answer)
          if (answer .eq. 'Y') then
            write (*,*) 'Complete distance matrix'
            write (*,'(10x,20(1x,a4))') (atmnam(iot(i)),i=1,napr)
            do i=1,napr
              write (*,'(1x,a4,1x,a4,20(f5.2))') 'DIST',
     +          atmnam(iot(i)),(dismat(i,j),j=1,napr)
            end do
            write (*,*)
          end if
c
          answer = 'N'
          call textin (
     +      ' Do you want to see the list of bonded distances ?',
     +      answer)
          call upcase (answer)
          if (answer .eq. 'Y') then
            write (*,*) 'List of bonded distances'
            ncnt = 0
            do i=1,napr-1
              do j=i+1,napr
                if (dismat(i,j) .le. bndcut) then
                   write (*,'(1x,a4,2(1x,a4,i4),1x,f8.2)')
     +               'BOND',atmnam(iot(i)),iresid(iot(i)),
     +               atmnam(iot(j)),iresid(iot(j)),
     +               dismat(i,j)
                   ncnt = ncnt + 1
                 end if
               end do
             end do
             call ivalut (' Nr of bonded distances :',1,ncnt)
          end if
c
          write (*,'(/a)') ' *** Bond angles ***'
          ncnt = 0
          do i=1,napr
            do j=1,napr-1
              if (i.eq.j) goto 8363
              if (dismat(i,j) .gt. bndcut) goto 8363
              do k=j+1,napr
                if (k.eq.j .or. k.eq.i) goto 8364
                if (dismat(i,k) .gt. bndcut) goto 8364
                write (*,'(1x,a5,3(1x,a4,i4),1x,f8.2)')
     +            'ANGLE',atmnam(iot(j)),iresid(iot(j)),
     +            atmnam(iot(i)),iresid(iot(i)),
     +            atmnam(iot(k)),iresid(iot(k)),
     +            angle(j,i,k,otxyz)
                  ncnt = ncnt + 1
 8364           continue
              end do
 8363         continue
            end do
          end do
          call ivalut (' Nr of bond angles :',1,ncnt)
c
          write (*,'(/a)') ' *** Dihedral angles ***'
          ncnt = 0
          do i=1,napr
            do j=1,napr
              if (i.eq.j) goto 8365
              if (dismat(i,j) .gt. bndcut) goto 8365
              do k=1,napr-1
                if (k.eq.j .or. k.eq.i) goto 8366
                if (dismat(i,k) .gt. bndcut) goto 8366
                do l=k+1,napr
                  if (l.eq.j .or. l.eq.i .or. l.eq.k) goto 8367
                  if (dismat(l,j) .gt. bndcut) goto 8367
                  write (*,'(1x,a8,4(1x,a4,i4),1x,f8.2)')
     +              'DIHEDRAL',atmnam(iot(k)),iresid(iot(k)),
     +              atmnam(iot(i)),iresid(iot(i)),
     +              atmnam(iot(j)),iresid(iot(j)),
     +              atmnam(iot(l)),iresid(iot(l)),
     +              tangle(k,i,j,l,otxyz)
                  ncnt = ncnt + 1
 8367             continue
                end do
 8366           continue
              end do
 8365         continue
            end do
          end do
          call ivalut (' Nr of 1-4 dihedrals :',1,ncnt)
c
        else if (igeo .eq. 2) then
c
c ... analyse N-CA-C backbone
c
          write (*,*)
          do i=1,6
            iot (i) = 0
          end do
          napr = 0
c
          do k=1,natoms
            if (atmnam(k) .eq. ' N  ') then
              napr = napr + 1
              iot (4) = k
              otxyz (1,4) = xf (k)
              otxyz (2,4) = yf (k)
              otxyz (3,4) = zf (k)
            else if (atmnam(k) .eq. ' CA ') then
              napr = napr + 1
              iot (5) = k
              otxyz (1,5) = xf (k)
              otxyz (2,5) = yf (k)
              otxyz (3,5) = zf (k)
            else if (atmnam(k) .eq. ' C  ') then
              napr = napr + 1
              iot (6) = k
              otxyz (1,6) = xf (k)
              otxyz (2,6) = yf (k)
              otxyz (3,6) = zf (k)
            else
              goto 8370
            end if
c
            if (iot(4) .gt. 0 .and. iot(5) .gt. 0 .and.
     +          iot(6) .gt. 0) then
              if (napr .gt. 3) then
                do j=4,6
                  i = iot (j)
                  write (line,4713) i,atmnam(i),altloc(i),
     +              resnam(i),achain(i),iresid(i),insert(i),
     +              dist(j-1,j,otxyz),angle(j-2,j-1,j,otxyz),
     +              tangle(j-3,j-2,j-1,j,otxyz)
                  write (*,'(a)') line(1:leng1(line))
                end do
              else
                j = 4
                i = iot (j)
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            0.0,0.0,0.0
                write (*,'(a)') line(1:leng1(line))
c
                j = 5
                i = iot (j)
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            dist(j-1,j,otxyz),0.0,0.0
                write (*,'(a)') line(1:leng1(line))
c
                j = 6
                i = iot (j)
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            dist(j-1,j,otxyz),angle(j-2,j-1,j,otxyz),0.0
                write (*,'(a)') line(1:leng1(line))
              end if
c
              do i=1,3
                do j=1,3
                  otxyz (j,i) = otxyz (j,i+3)
                end do
                iot (i)   = iot (i+3)
                iot (i+3) = 0
              end do
c
            end if
c
 8370       continue
          end do
c
          write (*,*)
          call ivalut (' Nr of N-CA-C atoms found :',1,napr)
c
        else if (igeo .eq. 3) then
c
c ... analyse CA backbone
c
          write (*,*)
          napr = 0
c
          do i=1,natoms
            if (atmnam(i) .eq. ' CA ') then
              napr = napr + 1
              j = 3*napr - 2
              xbuff (j)   = xf(i)
              xbuff (j+1) = yf(i)
              xbuff (j+2) = zf(i)
              j = napr
              if (napr .gt. 3) then
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            dist(j-1,j,xbuff),angle(j-2,j-1,j,xbuff),
     +            tangle(j-3,j-2,j-1,j,xbuff)
                write (*,'(a)') line(1:leng1(line))
              else if (napr .eq. 1) then
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            0.0,0.0,0.0
                write (*,'(a)') line(1:leng1(line))
              else if (napr .eq. 2) then
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            dist(j-1,j,xbuff),0.0,0.0
                write (*,'(a)') line(1:leng1(line))
              else if (napr .eq. 3) then
                write (line,4713) i,atmnam(i),altloc(i),
     +            resnam(i),achain(i),iresid(i),insert(i),
     +            dist(j-1,j,xbuff),angle(j-2,j-1,j,xbuff),0.0
                write (*,'(a)') line(1:leng1(line))
              end if
            end if
          end do
c
          write (*,*)
          call ivalut (' Nr of CA atoms found :',1,napr)
c
        end if
c
        write (*,*)
c
      else if (option(1:4) .eq. 'DIST') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        call disdis (natoms,allxyz,atmnam)
c
      else if (option(1:4) .eq. 'CHEC') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        ldummy = .false.
        call chknom (natoms,allxyz,atmnam,resnam,iresid,ldummy)
c
      else if (option(1:4) .eq. 'CORR') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        ldummy = .true.
        call chknom (natoms,allxyz,atmnam,resnam,iresid,ldummy)
c
      else if (option(1:4) .eq. 'CHI_') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        call calchi (natoms,allxyz,atmnam,resnam,iresid)
c
      else if (option(1:4) .eq. 'BURI') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        call buried (natoms,allxyz,atmnam,resnam,iresid)
c
      else if (option(1:4) .eq. 'RAMA') then
c
        do i=1,natoms
          allxyz (1,10+i) = xf (i)
          allxyz (2,10+i) = yf (i)
          allxyz (3,10+i) = zf (i)
        end do
c
        call calcpp (natoms,allxyz,atmnam,iresid,
     +    nres,resptr,phipsi)
c
        call dorama (natoms,iresid,nres,resptr,phipsi,resnam)
c
      else if (option(1:4) .eq. 'BALA') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        call calcpp (natoms,allxyz,atmnam,iresid,
     +    nres,resptr,phipsi)
c
        call dobala (natoms,iresid,nres,resptr,phipsi,resnam)
c
      else if (option(1:4) .eq. 'CACA') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        call calcca (natoms,allxyz,atmnam,iresid,
     +    nres,resptr,phipsi,1)
c
      else if (option(1:4) .eq. 'CA_R') then
c
        do i=1,natoms
          allxyz (1,i) = xf (i)
          allxyz (2,i) = yf (i)
          allxyz (3,i) = zf (i)
        end do
c
        call calcca (natoms,allxyz,atmnam,iresid,
     +    nres,resptr,phipsi,0)
c
        if (caqual(0,0) .lt. 0) call initca (caqual(0,-60))
c
        call docara (natoms,iresid,nres,resptr,phipsi,resnam,
     +               caqual(0,-60))
c
      else if (option(1:4) .eq. 'SUGG') then
c
        call ivalin (' Which residue number ?',1,ires)
        if (ires .le. 0) then
          call errcon ('Residue number must be positive !')
          goto 11
        end if
c
        do i=1,4
          iot (i) = 0
        end do
c
c ... find CA, C and O
c
        do i=1,natoms
          if (iresid(i) .eq. ires) then
            if (atmnam(i) .eq. ' N  ') then
              iot (1) = i
              write (*,*) '... found N'
            else if (atmnam(i) .eq. ' CA ') then
              iot (2) = i
              write (*,*) '... found CA'
            else if (atmnam(i) .eq. ' C  ') then
              iot (3) = i
              write (*,*) '... found C'
            else if (atmnam(i) .eq. ' O  ') then
              iot (4) = i
              write (*,*) '... found O'
            else if (atmnam(i) .eq. ' OT1') then
              iot (4) = i
              write (*,*) '... found OT1'
            end if
          end if
          if (iot(1) .gt. 0 .and. iot(2) .gt. 0 .and.
     +        iot(3) .gt. 0 .and. iot(4) .gt. 0) goto 9201
        end do
c
        if (iot(1) .le. 0) call errcon (' N not found')
        if (iot(2) .le. 0) call errcon (' CA not found')
        if (iot(3) .le. 0) call errcon (' C not found')
        if (iot(4) .le. 0) call errcon (' O(T1) not found')
        goto 11
c
c ... construct arrays
c
 9201   continue
        do i=1,4
          otxyz (1,i) = xf(iot(i))
          otxyz (2,i) = yf(iot(i))
          otxyz (3,i) = zf(iot(i))
        end do
c
        do i=1,5
          do j=1,3
            otint (j,i) = 0.0
            otref (j,i) = 0
          end do
        end do
c
        call fvalut (' Using C-OT distance (A)   :',1,codist)
        call fvalut (' Using OT-C-OT angle (deg) :',1,ocoang)
c
        otint (1,2) = dist (1,2,otxyz)
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
        call fvalut (' Distance C-OT1         =',1,
     +    dist(3,4,otxyz))
        call fvalut (' Reset to               =',1,otint(1,4))
c ... keep CA-COO flat
        otint (2,4) = 0.5 * (360.0 - ocoang)
        otref (2,4) = 2
        call fvalut (' Angle CA-C-OT1         =',1,
     +    angle(2,3,4,otxyz))
        call fvalut (' Reset to               =',1,otint(2,4))
        otint (3,4) = tangle (1,2,3,4,otxyz)
        otref (3,4) = 1
        call fvalut (' Dihedral N-CA-C-OT1    =',1,otint(3,4))
c
        otint (1,5) = codist
        otref (1,5) = 3
        call fvalut (' Distance C-OT2         =',1,otint(1,5))
        otint (2,5) = ocoang
        otref (2,5) = 4
        call fvalut (' Angle OT1-C-OT2        =',1,otint(2,5))
        otint (3,5) = 180.0
        otref (3,5) = 2
        call fvalut (' Dihedral CA-OT1-C-OT2  =',1,otint(3,5))
c
cc        call fvalut (' Int coords :',15,otint)
cc        call ivalut (' Ref coords :',15,otref)
cc        call fvalut (' XYZ coords :',15,otxyz)
c
        call c2car4 (5,otxyz,otint,otref)
c
        write (*,*)
cc        write (*,6374) ' N  ',xf(iot(1)),yf(iot(1)),zf(iot(1)),
cc     +    (otxyz(i,1),i=1,3)
cc        write (*,6374) ' CA ',xf(iot(2)),yf(iot(2)),zf(iot(2)),
cc     +    (otxyz(i,2),i=1,3)
cc        write (*,6374) ' C  ',xf(iot(3)),yf(iot(3)),zf(iot(3)),
cc     +    (otxyz(i,3),i=1,3)
        write (*,6374) ' OT1',xf(iot(4)),yf(iot(4)),zf(iot(4)),
     +    (otxyz(i,4),i=1,3)
        write (*,6374) ' OT2',0.0,0.0,0.0,(otxyz(i,5),i=1,3)
        write (*,*)
c
        write (*,*) 'Check geometry of carboxylate group :'
        call fvalut (' Dist       C-OT1 =',1,dist(3,4,otxyz))
        call fvalut (' Dist       C-OT2 =',1,dist(3,5,otxyz))
        call fvalut (' Dist     OT1-OT2 =',1,dist(4,5,otxyz))
        call fvalut (' Angle  OT1-C-OT2 =',1,angle(4,3,5,otxyz))
        call fvalut (' Angle   OT1-C-CA =',1,angle(4,3,2,otxyz))
        call fvalut (' Angle   OT2-C-CA =',1,angle(5,3,2,otxyz))
        call fvalut (' Dih   OT1-C-CA-N =',1,tangle(4,3,2,1,otxyz))
        call fvalut (' Dih   OT2-C-CA-N =',1,tangle(5,3,2,1,otxyz))
        call fvalut (' Dih CA-OT1-C-OT2 =',1,tangle(2,4,3,5,otxyz))
c
        write (*,*) ' ==> YOU MUST ADD/EDIT OT1/OT2 YOURSELF !!!'
        write (*,*)
c
 6374 format (' ==> ',a4,' NOW : ',3f8.3,' SUGGESTED : ',3f8.3)
c
      else if (option(1:4) .eq. 'WATE') then
c
        write (*,*) 'Waters can be called H2O, HOH, WAT, SOL, etc.'
        call textin (' What residue TYPE are your waters ?',awater)
        call upcase (awater)
c
        write (*,*) 'You may either swap XYZ,B,Q only, OR'
        write (*,*) 'also swap residue numbers etc.; in the'
        write (*,*) 'latter case, your waters may *NOT* be'
        write (*,*) 'sorted by residue ID afterwards !!!'
        call textin (' Swap residue numbers etc. ?',aswapm)
        call upcase (aswapm)
        lswapm = (aswapm .eq. 'Y')
c
        nwater = 0
        do i=1,natoms
          if (resnam(i) .eq. awater) then
            nwater = nwater + 1
            resptr (nwater) = i
            if (nwater .eq. 1) then
              xmin = xf (i)
              ymin = yf (i)
              zmin = zf (i)
              lowest = 1
            else if (xf(i) .lt. xmin .and.
     +               yf(i) .lt. ymin .and.
     +               zf(i) .lt. zmin) then
              lowest = nwater
              xmin = xf (i)
              ymin = yf (i)
              zmin = zf (i)
            end if
          end if
        end do
        call jvalut (' Nr of water molecules found :',1,nwater)
        if (nwater .le. 2) goto 11
c
        call jvalut (
     +    ' Water molecule with minimum coordinates :',1,lowest)
        iold = 1
        inew = lowest
c
 7383   continue
        if (iold .ne. inew) then
          jold = resptr (iold)
          jnew = resptr (inew)
          print *,' ... swap waters ',iresid(jold),' and ',iresid(jnew)
c
c ... swap XYZBQ
c
          call rswap (xf(jold),xf(jnew))
          call rswap (yf(jold),yf(jnew))
          call rswap (zf(jold),zf(jnew))
          call rswap (batom(jold),batom(jnew))
          call rswap (qatom(jold),qatom(jnew))
c
c ... swap other stuff as well ?
c
          if (lswapm) then
            call iswap (iresid(jold),iresid(jnew))
            call iswap (insert(jold),insert(jnew))
            write (line,'(a1,a2,a40)') altloc(jold),
     +        achain(jold),inote(jold)
            write (junk80,'(a1,a2,a40)') altloc(jnew),
     +        achain(jnew),inote(jnew)
            read (junk80,'(a1,a2,a40)') altloc(jold),
     +        achain(jold),inote(jold)
            read (line,'(a1,a2,a40)') altloc(jnew),
     +        achain(jnew),inote(jnew)
          end if
c
        end if
c
        icom = iold
        iold = iold + 1
        if (iold .ge. nwater) goto 7484
        inew = iold
        jold = resptr (icom)
        jnew = resptr (inew)
        dmin = (xf(jold)-xf(jnew))**2 + (yf(jold)-yf(jnew))**2 + 
     +         (zf(jold)-zf(jnew))**2
        if (iold+1 .le. nwater) then
          do i=iold+1,nwater
            jnew = resptr(i)
            dij  = (xf(jold)-xf(jnew))**2 + (yf(jold)-yf(jnew))**2 + 
     +             (zf(jold)-zf(jnew))**2
            if (dij .lt. dmin) then
              dmin = dij
              inew = i
            end if
          end do
        end if
c
        jnew = resptr (inew)
        dmin = sqrt(dmin)
        call jvalut (' Closest water is :',1,iresid(jnew))
        call fvalut (' Distance (A)     :',1,dmin)
        if (dmin .lt. 2.5) write (*,*) ' ==> CLOSER THAN 2.5 A !!!'
c
        goto 7383
c
 7484   continue
        call jvalut (' Nr of waters sorted :',1,nwater)
c
      else if (option(1:4) .eq. 'INVE') then
c
        call ivalin (
     +    ' Residue range to invert (0 0 = all molecule) ?',2,izone)
        ldum = 0
        cog (1) = 0.0
        cog (2) = 0.0
        cog (3) = 0.0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 6672
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 6672
          cog(1) = cog(1) + xf(i)
          xf (i) = -xf(i)
          cog(2) = cog(2) + yf(i)
          yf (i) = -yf(i)
          cog(3) = cog(3) + zf(i)
          zf (i) = -zf(i)
          ldum = ldum + 1
 6672     continue
        end do
        call jvalut (' Nr of atoms affected :',1,ldum)
        if (ldum .gt. 0) then
          cog (1) = 2.0 * cog(1) / float(ldum)
          cog (2) = 2.0 * cog(2) / float(ldum)
          cog (3) = 2.0 * cog(3) / float(ldum)
          call fvalut (' Translation to apply :',3,cog)
          do i=1,natoms
            if (izone(1) .ne. 0 .and.
     +        iresid(i) .lt. izone(1)) goto 6674
            if (izone(2) .ne. 0 .and.
     +        iresid(i) .gt. izone(2)) goto 6674
            xf (i) = xf(i) + cog(1)
            yf (i) = yf(i) + cog(2)
            zf (i) = zf(i) + cog(3)
 6674       continue
          end do
          call jvalut (' Nr of atoms translated :',1,ldum)
        end if
c
      else if (option(1:4) .eq. 'MIRR') then
c
        call textin (
     +    ' Mirror X, Y, or Z (X/Y/Z) ?',mplane)
        call upcase (mplane)
        if (mplane .ne. 'Y' .and. mplane .ne. 'Z') mplane = 'X'
c
        call ivalin (
     +    ' Residue range to mirror (0 0 = all molecule) ?',2,izone)
        ldum = 0
        cog (1) = 0.0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 6472
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 6472
          if (mplane .eq. 'X') then
            cog(1) = cog(1) + xf(i)
            xf (i) = -xf(i)
          else if (mplane .eq. 'Y') then
            cog(1) = cog(1) + yf(i)
            yf (i) = -yf(i)
          else if (mplane .eq. 'Z') then
            cog(1) = cog(1) + zf(i)
            zf (i) = -zf(i)
          end if
          ldum = ldum + 1
 6472     continue
        end do
        call jvalut (' Nr of atoms affected :',1,ldum)
        if (ldum .gt. 0) then
          cog (1) = 2.0 * cog(1) / float(ldum)
          call fvalut (' Translation to apply :',1,cog)
          do i=1,natoms
            if (izone(1) .ne. 0 .and.
     +        iresid(i) .lt. izone(1)) goto 6474
            if (izone(2) .ne. 0 .and.
     +        iresid(i) .gt. izone(2)) goto 6474
            if (mplane .eq. 'X') then
              xf (i) = xf(i) + cog(1)
            else if (mplane .eq. 'Y') then
              yf (i) = yf(i) + cog(1)
            else if (mplane .eq. 'Z') then
              zf (i) = zf(i) + cog(1)
            end if
 6474       continue
          end do
          call jvalut (' Nr of atoms translated :',1,ldum)
        end if
c
      else if (option(1:4) .eq. 'ORIG') then
c
        call ivalin (
     +    ' Residue range to move (0 0 = all molecule) ?',2,izone)
        ldum = 0
        cog (1) = 0.0
        cog (2) = 0.0
        cog (3) = 0.0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 6572
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 6572
          cog(1) = cog(1) + xf(i)
          cog(2) = cog(2) + yf(i)
          cog(3) = cog(3) + zf(i)
          ldum = ldum + 1
 6572     continue
        end do
        call jvalut (' Nr of atoms affected :',1,ldum)
        if (ldum .gt. 0) then
          cog (1) = -cog(1)/float(ldum)
          cog (2) = -cog(2)/float(ldum)
          cog (3) = -cog(3)/float(ldum)
          call fvalut (' Translation to apply :',3,cog)
          do i=1,natoms
            if (izone(1) .ne. 0 .and.
     +        iresid(i) .lt. izone(1)) goto 6574
            if (izone(2) .ne. 0 .and.
     +        iresid(i) .gt. izone(2)) goto 6574
            xf(i) = xf(i) + cog(1)
            yf(i) = yf(i) + cog(2)
            zf(i) = zf(i) + cog(3)
 6574       continue
          end do
          call jvalut (' Nr of atoms translated :',1,ldum)
        end if
c
      else if (option(1:4) .eq. 'LIMI') then
c
        call fvalin (
     +    ' Enter MIN and MAX temperature factor :',2,blim)
        call fvalin (
     +    ' Enter MIN and MAX occupancy          :',2,qlim)
        call ivalin (
     +    ' Residue range to apply (0 0 = all molecule) ?',2,izone)
        ldum = 0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 2021
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 2021
          batom (i) = max (blim(1), min (batom(i), blim(2)))
          qatom (i) = max (qlim(1), min (qatom(i), qlim(2)))
          ldum = ldum + 1
 2021     continue
        end do
        call jvalut (' Nr of atoms updated :',1,ldum)
c
      else if (option(1:4) .eq. 'ALTE') then
c
        call textin (' OLD residue name ?',resold)
        call textin (' NEW residue name ?',resnew)
c
        iold = 0
        do i=1,natoms
          if (resnam(i) .eq. resold) then
            iold = iold + 1
            resnam (i) = resnew
          end if
        end do
        call jvalut (' Nr of atoms affected :',1,iold)
c
      else if (option(1:4) .eq. 'XPLO') then
c
        newchx = ' '
        call textin (' X-PLOR chain label (4 characters) ?',newchx)
        call ivalin (
     +    ' Residue range to apply (0 0 = all molecule) ?',2,izone)
        ldum = 0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 2031
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 2031
          inote (i) (7:10) = newchx
          ldum = ldum + 1
 2031     continue
        end do
        call jvalut (' Nr of X-PLOR chain labels updated :',1,ldum)
c
      else if (option(1:4) .eq. 'AUTO') then
c
        call prompt (' Generating chain and segids ...')
c
        j = 0
        iold = -999
        pres = '???'
c
        do i=1,natoms
          if (iresid (i) .ne. iold) then
            if (iresid(i) .eq. (iold+1)) then
              iold = iresid (i)
              pres = resnam (i)
            else if (lwater(pres) .and.
     +               lwater(resnam(i))) then
              iold = iresid (i)
              pres = resnam (i)
            else
              j = j + 1
              if (j .eq. 27) then
                call errcon ('More than 26 chains; rest = "Z" !')
              end if
              nowchn = ' '//char (ichar('A') - 1 + min(j,26))
              newchx (1:1) = nowchn(2:2)
              newchx (2:2) = nowchn(2:2)
              newchx (3:3) = nowchn(2:2)
              newchx (4:4) = nowchn(2:2)
              write (*,'(1x,a,a,a,a,a,i4)') 'New chain ',
     +          nowchn,', segid ',newchx,' @ residue ',
     +          iresid(i)
              iold = iresid (i)
              pres = resnam (i)
            end if
          end if
          achain (i) = nowchn
          inote (i) (7:10) = newchx
        end do
        call jvalut (' Nr of segments found :',1,j)
c
c
      else if (option(1:4) .eq. 'ASK_') then
c
        call prompt (' Generating chain and segids ...')
c
        j = 0
        iold = -999
        pres = '???'
c
        do i=1,natoms
          if (iresid (i) .ne. iold) then
            if (iresid(i) .eq. (iold+1)) then
              iold = iresid (i)
              pres = resnam (i)
            else if (lwater(pres) .and.
     +               lwater(resnam(i))) then
              iold = iresid (i)
              pres = resnam (i)
            else
              j = j + 1
              j = max (1, min (26, j))
              nowchn = ' '//char (ichar('A') - 1 + min(j,26))
              newchx (1:1) = nowchn(2:2)
              newchx (2:2) = nowchn(2:2)
              newchx (3:3) = nowchn(2:2)
              newchx (4:4) = nowchn(2:2)
c
              write (line,4713) i,atmnam(i),altloc(i),
     +          resnam(i),achain(i),iresid(i),insert(i),xf(i),
     +          yf(i),zf(i),qatom(i),batom(i),inote(i)
              write (*,'(/a)') line(1:leng1(line))
c
              call textin (' New chain name ?',nowchn)
              call upcase (nowchn)
              newchx (1:1) = nowchn(2:2)
              newchx (2:2) = nowchn(2:2)
              newchx (3:3) = nowchn(2:2)
              newchx (4:4) = nowchn(2:2)
              call textin (' New SEGId name ?',newchx)
              call upcase (newchx)
c
              write (*,'(1x,a,a,a,a,a,i4)') 'New chain ',
     +          nowchn,', segid ',newchx,' @ residue ',
     +          iresid(i)
              iold = iresid (i)
              pres = resnam (i)
            end if
          end if
          achain (i) = nowchn
          inote (i) (7:10) = newchx
        end do
        write (*,*)
        call jvalut (' Nr of segments found :',1,j)
c
      else if (option(1:4) .eq. 'FROM') then
c
        nowchn = achain (1)
        call textin (' Chain label*2 to translate ?',nowchn)
        newchx = nowchn//nowchn
        call remspa (newchx)
        if (length(newchx).eq.2) newchx=newchx(1:2)//newchx(1:2)
        call upcase (newchx)
        call textin (' Corresponding X-PLOR chain label*4 ?',newchx)
c
        iold = 0
        do i=1,natoms
          if (achain(i) .eq. nowchn) then
            inote (i) (7:10) = newchx
            iold = iold + 1
          end if
        end do
        call jvalut (' Nr of chain labels translated :',1,iold)
c
      else if (option(1:4) .eq. 'XID_') then
c
        call textin (' X-PLOR chain label*4 to translate ?',newchx)
        nowchn = ' '//newchx(1:1)
        call upcase (nowchn)
        call textin (' Corresponding chain label*2 ?',nowchn)
c
        iold = 0
        do i=1,natoms
          if (inote (i) (7:10) .eq. newchx) then
            achain (i) = nowchn
            iold = iold + 1
          end if
        end do
        call jvalut (' Nr of X-PLOR labels translated :',1,iold)
c
      else if (option(1:4) .eq. 'AVER') then
c
        write (*,'(99(a/:))')
     +    ' Valid options are:',
     +    ' 1. Average over all atoms (i.e., compute Boverall)',
     +    ' 2. Average per residue over all atoms',
     +    ' 3. Average per residue, separately for main and side-chain',
     +    ' 4. Average corresponding atoms in different chains'
c
        call ivalin (' Option ?',1,mave)
c
        if (mave. lt. 1 .or. mave .gt. 4) then
          call ivalut (' ERROR - Invalid option :',1,mave)
          goto 11
        end if
c
        if (mave .eq. 1) then
c
          bsum = 0.0
          do i=1,natoms
            bsum = bsum + batom (i)
          end do
c
          bave = bsum / float(natoms)
          do i=1,natoms
            batom (i) = bave
          end do
c
          call jvalut (' Nr of temperature factors updated :',1,natoms)
          call fvalut (' Average temperature factor :',1,bave)
c
        else if (mave .eq. 2) then
c
          iold = 1
          bsum = batom (1)
          id   = iresid (1)
c
          do i=2,natoms
            if (iresid(i) .eq. id) then
              bsum = bsum + batom (i)
            else
              nats = i - iold
              bave = bsum / float (nats)
              do j = iold, i-1
                batom (j) = bave
              end do
              write (*,6837) id,nats,bave
              iold = i
              bsum = batom (i)
              id   = iresid (i)
            end if
          end do
c
          nats = natoms - iold
          bave = bsum / float (nats)
          do j = iold, natoms
            batom (j) = bave
          end do
          write (*,6837) id,nats,bave
c
          call jvalut (' Nr of temperature factors updated :',1,natoms)
c
        else if (mave .eq. 3) then
c
          iold = 1
          bsum = 0.0
          bsid = 0.0
          nsum = 0
          nsid = 0
          id   = iresid (1)
          if (mainch(atmnam(1))) then
            bsum = bsum + batom (1)
            nsum = nsum + 1
          else
            bsid = bsid + batom (1)
            nsid = nsid + 1
          end if
c
          do i=2,natoms
            if (iresid(i) .eq. id) then
              if (mainch(atmnam(i))) then
                bsum = bsum + batom (i)
                nsum = nsum + 1
              else
                bsid = bsid + batom (i)
                nsid = nsid + 1
              end if
            else
              nats = i - iold
              bave = 0.0
              basc = 0.0
              if (nsum .gt. 0) bave = bsum / float (nsum)
              if (nsid .gt. 0) basc = bsid / float (nsid)
c
              do j = iold, i-1
                if (mainch(atmnam(j))) then
                  batom (j) = bave
                else
                  batom (j) = basc
                end if
              end do
              write (*,6839) id,nats,bave,nsum,basc,nsid
c
              iold = i
              bsum = 0.0
              bsid = 0.0
              nsum = 0
              nsid = 0
              id   = iresid (i)
              if (mainch(atmnam(i))) then
                bsum = bsum + batom (i)
                nsum = nsum + 1
              else
                bsid = bsid + batom (i)
                nsid = nsid + 1
              end if
            end if
          end do
c
          nats = nsum + nsid
          bave = 0.0
          basc = 0.0
          if (nsum .gt. 0) bave = bsum / float (nsum)
          if (nsid .gt. 0) basc = bsid / float (nsid)
c
          do j = iold, natoms
            if (mainch(atmnam(j))) then
              batom (j) = bave
            else
              batom (j) = basc
            end if
          end do
          write (*,6839) id,nats,bave,nsum,basc,nsid
c
          call jvalut (' Nr of temperature factors updated :',1,natoms)
c
        else if (mave .eq. 4) then
c
          do i=1,natoms
            baveok (i) = .false.
          end do
c
          nalone = 0
          nmulti = 0
c
          do i=1,natoms
            if (.not. baveok(i)) then
c
              ncopies = 1
              icopies (1) = i
              idaddy = i
              btotal = batom (i)
c
              if (i .lt. natoms) then
                do j=i+1,natoms
                  if (.not. baveok(j)) then
                    if (iresid(j) .eq. iresid(idaddy)) then
                      if (resnam(j) .eq. resnam (idaddy) .and.
     +                    atmnam(j) .eq. atmnam (idaddy)) then
                        ncopies = ncopies + 1
                        icopies (ncopies) = j
                        btotal = btotal + batom (j)
                      end if
                    end if
                  end if
                end do
              end if
c
              bave = btotal / float (ncopies)
c
              do j=1,ncopies
                iat = icopies (j)
                batom (iat)  = bave
                baveok (iat) = .true.
              end do
c
              if (ncopies .eq. 1) then
                nalone = nalone + 1
              else
                nmulti = nmulti + 1
              end if
c
            end if
          end do
c
          call jvalut (' Nr of atoms averaged       :',1,nmulti)
          call jvalut (' Nr of atoms occurring once :',1,nalone)
c
        end if
c
 6837 format ('  Res ',i5,5x,'Nr_atoms ',i2,5x,'Bave ',f10.2)
 6839 format ('  Res ',i5,5x,'Nr_atoms ',i2,5x,'Bave-MC ',f10.2,' (',i2,
     +  ')',5x,'Bave-SC ',f10.2,' (',i2,')')
c
      else if (option(1:4) .eq. 'OCCU') then
c
        call fvalin (' OFnew = A * OFold + B.  Enter A, B :',2,qfacts)
        call ivalin (
     +    ' Residue range to apply (0 0 = all molecule) ?',2,izone)
        ldum = 0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 2041
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 2041
          qatom (i) = qfacts(1)*qatom(i) + qfacts(2)
          ldum = ldum + 1
 2041     continue
        end do
        call jvalut (' Nr of occupancy factors updated :',1,ldum)
c
      else if (option(1:4) .eq. 'FLAG') then
c
        call prompt (' Select one of the following:')
        call prompt (' 1) Dutch/Italian etc. flag colours')
        call prompt (' 2) Swedish/Danish etc. flag colours')
        call ivalin (' Which type of flag ?',1,iflag)
        iflag = max(1,min(2,iflag))
c
        if (iflag .eq. 1) then
          call xstats (yf,natoms,ave,sdv,xmin,xmax,xtot)
          xf1 = xmin + (xmax-xmin)/3.0
          xf2 = xmax - (xmax-xmin)/3.0
          do i=1,natoms
            if (yf(i) .ge. xf2) then
              ibuff(i)=1
            else if (yf(i) .le. xf1) then
              ibuff(i)=3
            else
              ibuff(i)=2
            end if
          end do
        else if (iflag .eq. 2) then
          call xstats (xf,natoms,ave,sdv,xmin,xmax,xtot)
          call xstats (yf,natoms,ave,sdv,ymin,ymax,xtot)
          xf1 = (xmin+xmax)/3.0 - (xmax-xmin)/8.0
          xf2 = (xmin+xmax)/3.0 + (xmax-xmin)/8.0
          yf1 = (ymin+ymax)/2.0 - (ymax-ymin)/8.0
          yf2 = (ymin+ymax)/2.0 + (ymax-ymin)/8.0
          do i=1,natoms
            ibuff (i) = 1
            if (xf(i).ge.xf1 .and. xf(i).le.xf2) ibuff(i)=2
            if (yf(i).ge.yf1 .and. yf(i).le.yf2) ibuff(i)=2
          end do
        end if
c
        call textin (' Molecule name in O ?',omol)
        flagfn = omol//'_atom_flag.odb'
        call remspa (flagfn)
        call locase (flagfn)
        call textin (' O datablock file ?',flagfn)
        call xopxua (f4,flagfn,xinter(),ierr)
        if (ierr .ne. 0) goto 11
        junk80 = omol//'_atom_flag'
        call remspa (junk80)
        call locase (junk80)
        write (line,'(a,a,i10,a)') junk80(1:leng1(junk80)),
     +    ' I ',natoms,' (80i1)'
        call pretty (line)
        write (f4,'(a)') line(1:leng1(line))
        write (f4,'(80i1)') (ibuff(i),i=1,natoms)
        close (f4)
        call prompt (' Datablock written; use Paint_case in O')
c
      else if (option(1:4) .eq. 'RESI') then
c
        iold = 0
        inew = 1
        call ivalin (' First NEW residue number ?',1,inew)
        inew = inew - 1
c
        do i=1,natoms
          if (iresid (i) .ne. iold) then
            iold = iresid (i)
            inew = inew + 1
            iresid (i) = inew
c           print *,' changed ',iold,' to ',inew,' atom ',i
          else
            iresid (i) = inew
          end if
        end do
c
        call jvalut (' Last residue number :',1,inew)
c
      else if (option(1:4) .eq. 'O2XH') then
c
        call prompt (' Fixing hydrogen names for Asn & Arg')
c
        iold = -1
        inew = 0
c
        do i=1,natoms
          if (iresid (i) .ne. iold) then
            iold = iresid (i)
            hd21 = .false.
            hd22 = .false.
            hh11 = .false.
            hh12 = .false.
            hh21 = .false.
            hh22 = .false.
          end if
c
          if (resnam(i) .eq. 'ASN') then
            if (atmnam(i) .eq. ' HD2') then
              if (.not. hd21) then
                atmnam (i) = 'HD21'
                hd21 = .true.
              else if (.not. hd22) then
                atmnam (i) = 'HD22'
                hd22 = .true.
              else
                call errcon ('Too many HD2 atoms ...')
                atmnam (i) = 'HD2?'
              end if
              write (*,6503) resnam(i),iresid(i),
     +            ' HD2',atmnam(i)
              inew = inew + 1
            end if
          else if (resnam(i) .eq. 'ARG') then
            if (atmnam(i) .eq. ' HH1') then
              if (.not. hh11) then
                atmnam (i) = 'HH11'
                hh11 = .true.
              else if (.not. hh12) then
                atmnam (i) = 'HH12'
                hh12 = .true.
              else
                call errcon ('Too many HH1 atoms ...')
                atmnam (i) = 'HH1?'
              end if
              write (*,6503) resnam(i),iresid(i),
     +            ' HH1',atmnam(i)
              inew = inew + 1
            else if (atmnam(i) .eq. ' HH2') then
              if (.not. hh21) then
                atmnam (i) = 'HH21'
                hh21 = .true.
              else if (.not. hh22) then
                atmnam (i) = 'HH22'
                hh22 = .true.
              else
                call errcon ('Too many HH2 atoms ...')
                atmnam (i) = 'HH2?'
              end if
              write (*,6503) resnam(i),iresid(i),
     +            ' HH2',atmnam(i)
              inew = inew + 1
            end if
          end if
        end do
c
        call jvalut (' Nr of hydrogens renamed :',1,inew)
c
 6503 format (' Renamed ',a3,i4,1x,a4,' to ',a4)
c
      else if (option(1:4) .eq. 'ZONE') then
c
        call ivalin (' First residue to renumber    ?',1,iz1)
        call ivalin (' Last residue to renumber     ?',1,iz2)
        call ivalin (' NEW number for first residue ?',1,iz3)
c
        ndum = 0
        iold = 0
        inew = iz3 - 1
c
        do i=1,natoms
          if (iresid (i) .ge. iz1 .and.
     +        iresid (i) .le. iz2) then
            ndum = ndum + 1
            if (iresid (i) .ne. iold) then
              iold = iresid (i)
              inew = inew + 1
              iresid (i) = inew
c              print *,' changed ',iold,' to ',inew,' atom ',i
            else
              iresid (i) = inew
            end if
          end if
        end do
c
        call jvalut (' Last residue number :',1,inew)
        call jvalut (' Nr of atoms changed :',1,ndum)
c
      else if (option(1:4) .eq. 'RENU') then
c
        inew = 1
        call ivalin (' First NEW atom number ?',1,inew)
        inew = inew - 1
c
        do i=1,natoms
          atomnr (i) = inew + i
        end do
c
        call jvalut (' Nr of atom numbers updated :',1,natoms)
        call jvalut (' Nr of last atom            :',1,
     +    atomnr(natoms))
c
      else if (option(1:4) .eq. 'ATRE') then
c
        inew = 1
        call ivalin (' First NEW atom/residue number ?',1,inew)
        inew = inew - 1
c
        do i=1,natoms
          atomnr (i) = inew + i
          iresid (i) = inew + i
        end do
c
        call jvalut (' Nr of atoms/residues updated :',1,natoms)
        call jvalut (' Nr of last atom/residue      :',1,
     +    atomnr(natoms))
c
      else if (option(1:4) .eq. 'STAT') then
c
 8702 format (' ',a10,4a15)
 8704 format (' ',a10,4f15.3)
c
        if (natoms .gt. 0) then
          call jvalut (' Nr of atom numbers in memory :',1,natoms)
          write (*,*)
          write (*,8702) 'Item','Average','St.Dev','Min','Max'
          write (*,8702) '----','-------','------','---','---'
c
          call xstats (xf,natoms,ave,sdv,xmin,xmax,xtot)
          write (*,8704) 'X-coord',ave,sdv,xmin,xmax
          xn = ave
          call xstats (yf,natoms,ave,sdv,xmin,xmax,xtot)
          write (*,8704) 'Y-coord',ave,sdv,xmin,xmax
          yn = ave
          call xstats (zf,natoms,ave,sdv,xmin,xmax,xtot)
          write (*,8704) 'Z-coord',ave,sdv,xmin,xmax
          zn = ave
          call xstats (batom,natoms,ave,sdv,xmin,xmax,xtot)
          write (*,8704) 'B-factor',ave,sdv,xmin,xmax
          call xstats (qatom,natoms,ave,sdv,xmin,xmax,xtot)
          write (*,8704) 'Occpncy',ave,sdv,xmin,xmax
          write (*,*)
c
c ... radius of gyration
c
          rog = 0.0
          do i=1,natoms
            q = (xf(i)-xn)**2 + (yf(i)-yn)**2 + (zf(i)-zn)**2
            rog = rog + q
          end do
          rog = rog / float(natoms)
          rog = sqrt (rog)
          write (*,'(1x,a,f12.2)') 'Radius of gyration (A) : ',rog
c
c ... centre-of-mass
c
          xn = 0.0
          yn = 0.0
          zn = 0.0
          som = 0.0
          do i=1,natoms
            if (lhydro(atmnam(i))) then
              junk80 = ' H'
            else
              junk80 = atmnam(i)(1:2)
            end if
            call elinfo (junk80(1:2),junk80(3:),j,q,radius,.false.)
            if (j .le. 0) q = 12.0
            xn = xn + xf(i)*q
            yn = yn + yf(i)*q
            zn = zn + zf(i)*q
            som = som + q
          end do
          xn = xn/som
          yn = yn/som
          zn = zn/som
          write (*,'(1x,a,f21.3)')  'Sum of masses  : ',som
          write (*,'(1x,a,3f10.2)') 'Centre-of-mass : ',xn,yn,zn
        else
          call errcon ('No atoms read yet !')
        end if
c
      else if (option(1:4) .eq. 'ROTA') then
c
        write (*,*) 'Options:'
        write (*,*) '1 - enter your own rotation matrix'
        write (*,*) '2 - CCP4 (ALMN) Eulerian angles'
        write (*,*) '3 - CCP4 (ALMN) Omega/Phi/Kappa'
        write (*,*) '4 - MERLOT Eulerian angles'
        write (*,*) '5 - MERLOT Polar angles'
        write (*,*) 'Options 2 - 5 programmed by Rams'
        call ivalin (' Option ?',1,irot)
c
        if (irot .eq. 1) then
          call fvalin (
     +      ' Xnew = A*Xold + B*Yold + C*Zold.  Enter A-C :',3,xsym)
          call fvalin (
     +      ' Ynew = D*Xold + E*Yold + F*Zold.  Enter D-F :',3,ysym)
          call fvalin (
     +      ' Znew = G*Xold + H*Yold + I*Zold.  Enter G-I :',3,zsym)
c
          do i=1,natoms
            xn = xsym(1)*xf(i) + xsym(2)*yf(i) + xsym(3)*zf(i)
            yn = ysym(1)*xf(i) + ysym(2)*yf(i) + ysym(3)*zf(i)
            zn = zsym(1)*xf(i) + zsym(2)*yf(i) + zsym(3)*zf(i)
            xf (i) = xn
            yf (i) = yn
            zf (i) = zn
          end do
          call jvalut (' Nr of atoms rotated :',1,natoms)
c
        else if (irot .ge. 2 .and. irot .le. 5) then
c
          do i=1,4
            do j=1,3
              rotmat (j,i) = 0.0
            end do
          end do
c
          if (irot .eq. 2) then
            call fvalin (' Euler angles ?',3,rotang)
            call ccpeul (rotang,rotmat)
          else if (irot .eq. 3) then
            call fvalin (' Omega, Phi, Kappa ?',3,rotang)
            call ccppol (rotang,rotmat)
          else if (irot .eq. 4) then
            call fvalin (' Alpha, Beta, Gamma ?',3,rotang)
            call mereul (rotang,rotmat)
          else if (irot .eq. 5) then
            call fvalin (' Phi, Psi, Kappa ?',3,rotang)
            call merpol (rotang,rotmat)
          endif
c
ccc          call matana (rotmat)
c
          call anancs (1,rotmat,.true.,ierr)
c
          do i=1,natoms
            xcar (1) = xf(i)
            xcar (2) = yf(i)
            xcar (3) = zf(i)
            do j=1,3
              xfra (j) = 0.0
              do k=1,3
                xfra (j) = xfra(j) + xcar(k)*rotmat(j,k)
              end do
            end do
            xf (i) = xfra (1)
            yf (i) = xfra (2)
            zf (i) = xfra (3)
          end do
          call jvalut (' Nr of atoms rotated :',1,natoms)
c
        else
          call errcon ('Invalid option')
        end if
c
      else if (option(1:4) .eq. 'CART') then
c
        call fvalin (' Unit-cell constants ?',6,cell)
        call orthog (cell,otof,1)
        call orthog (cell,ftoo,-1) 
        call rvalut (' Matrix :',3,otof(1,1))
        call rvalut (' Matrix :',3,otof(1,2))
        call rvalut (' Matrix :',3,otof(1,3))
        do i=1,natoms
          xfra (1) = xf(i)
          xfra (2) = yf(i)
          xfra (3) = zf(i)
          do j=1,3
            xcar (j)= 0
            do k=1,3
              xcar (j) = xcar (j) + xfra(k)*otof(j,k)
            end do
          end do
          xf (i) = xcar (1)
          yf (i) = xcar (2)
          zf (i) = xcar (3)
        end do
        call jvalut (' Nr of atoms converted :',1,natoms)
c
      else if (option(1:4) .eq. 'FRAC') then
c
        call fvalin (' Unit-cell constants ?',6,cell)
        call orthog (cell,otof,1)
        call orthog (cell,ftoo,-1)
        call rvalut (' Matrix :',3,ftoo(1,1))
        call rvalut (' Matrix :',3,ftoo(1,2))
        call rvalut (' Matrix :',3,ftoo(1,3))
        do i=1,natoms
          xcar (1) = xf(i)
          xcar (2) = yf(i)
          xcar (3) = zf(i)
          do j=1,3
            xfra (j)= 0
            do k=1,3
              xfra (j) = xfra (j) + xcar(k)*ftoo(j,k)
            end do
          end do
          xf (i) = xfra (1)
          yf (i) = xfra (2)
          zf (i) = xfra (3)
        end do
        call jvalut (' Nr of atoms converted :',1,natoms)
c
      else if (option(1:4) .eq. 'TRAN') then
c
        call ivalin (' 1 = Cartesian, 2 = Fractional. Option ?',1,itra)
c
        if (itra .eq. 1) then
          call fvalin (' Translation vector ?',3,trans)
          do i=1,natoms
            xf(i) = xf(i) + trans(1)
            yf(i) = yf(i) + trans(2)
            zf(i) = zf(i) + trans(3)
          end do
          call jvalut (' Nr of atoms translated :',1,natoms)
c
        else if (itra .eq. 2) then
          call fvalin (' Unit cell constants    ?',6,cell)
          call fvalin (' Fractional translation ?',3,ftrans)
          call orthog (cell,otof,1)
          call orthog (cell,ftoo,-1) 
          do i=1,natoms
            xcar (1) = xf(i)
            xcar (2) = yf(i)
            xcar (3) = zf(i)
            do j=1,3
              xfra (j) = ftrans(j)
              do k=1,3
                xfra (j) = xfra(j) + xcar(k)*otof(j,k)
              end do
            end do
            do j=1,3
              xcar (j)= 0
              do k=1,3
                xcar (j) = xcar (j) + xfra(k)*ftoo(j,k)
              end do
            end do
            xf (i) = xcar (1)
            yf (i) = xcar (2)
            zf (i) = xcar (3)
          end do
          call jvalut (' Nr of atoms translated :',1,natoms)
c
        else
          call errcon ('Invalid option')
        end if
c
      else if (option(1:4) .eq. 'APPL') then
c
        call jvalin (' Seed for random number generator     ?',1,iseed)
        call gkrand (dx,0.0,0.0,iseed)
        iseed = 0
c
        call ranrot (rotmat)
c
        do i=1,natoms
          xn = rotmat(1,1)*xf(i) + rotmat(1,2)*yf(i) + rotmat(1,3)*zf(i)
          yn = rotmat(2,1)*xf(i) + rotmat(2,2)*yf(i) + rotmat(2,3)*zf(i)
          zn = rotmat(3,1)*xf(i) + rotmat(3,2)*yf(i) + rotmat(3,3)*zf(i)
          xf (i) = xn
          yf (i) = yn
          zf (i) = zn
        end do
        call jvalut (' Nr of atoms rotated :',1,natoms)
c
      else if (option(1:4) .eq. 'CHAI') then
c
        newchn = ' '
        call textin (' Chain label (2 characters) ?',newchn)
        call ivalin (
     +    ' Residue range to apply (0 0 = all molecule) ?',2,izone)
        ldum = 0
        do i=1,natoms
          if (izone(1) .ne. 0 .and.
     +      iresid(i) .lt. izone(1)) goto 2051
          if (izone(2) .ne. 0 .and.
     +      iresid(i) .gt. izone(2)) goto 2051
          achain (i) = newchn
          ldum = ldum + 1
 2051     continue
        end do
        call jvalut (' Nr of chain labels updated :',1,natoms)
c
      else if (option(1:4) .eq. 'WRIT' .or.
     +         option(1:4) .eq. 'DUMP') then
c
        ldump = (option(1:4) .eq. 'DUMP')
c
        write (*,*)
        file2 = 'out.pdb'
        call textin (
     +    ' Output PDB file                               ?',file2)
        close (f2)
        call xopxna (f2,file2,xinter(),ierr)
        if (ierr .ne. 0) goto 998
c
        if (.not. ldump) call textin (
     +    ' REMARK at start of file ?',commnt)
c
        if (.not. ldump) then
          call textin (
     +      ' Copy all REMARK, HEADER etc. cards from input ?',arema)
          call upcase (arema)
        end if
        lrema = (arema .eq. 'Y')
c
        if (.not. ldump) then
          call textin (
     +      ' Which chain to write (** = any and all)       ?',whichn)
        end if
c
        if (.not. ldump) call ivalin (
     +    ' Residue range to write (0 0 = all molecule)   ?',2,izone)
c
        if (.not. ldump) then
          write (*,*) 'You may output All atoms, only Main-chain atoms,'
          write (*,*) 'a Poly-alanine (Gly intact), a poly-Serine,'
          write (*,*) '(Gly and Ala intact) or a poly-Glycine'
          call textin (
     +      ' Which option do you want (All/M/P/S/G)        ?',which)
          call upcase (which)
          if (which .ne. 'M' .and. which .ne. 'P' .and.
     +        which .ne. 'S' .and. which .ne. 'G') which = 'A'
          if (which .eq. 'P') write (*,*)
     +      'NOTE: GLYcines are left intact !!!'
          if (which .eq. 'S') then
            write (*,*)
     +       'NOTE: GLYcines and ALAnines are left intact !!!'
            write (*,*)
     +        'NOTE: atoms CG and CG1 are renamed to OG !!!'
          end if
        end if
c
        if (.not. ldump) call textin (
     +    ' Write HYDROGEN atoms (Y/N)                    ?',keephs)
        call upcase (keephs)
        lkeeph = (keephs .eq. 'Y')
c
        if (.not. ldump) call textin (
     +    ' Force consecutive atom numbering (Y/N)        ?',aforce)
        call upcase (aforce)
        lforce = (aforce .eq. 'Y')
c
        if (.not. ldump) then
          write (*,*) 'X-PLOR needs OT1 and OT2, but O hates them'
          write (*,*) 'If your file contains OT1/2 you may either'
          write (*,*) 'keep them, or replace them by O/OXT'
          call textin (
     +      ' Write X-PLOR OT1/2 ? (Y/N)                     ?',astrip)
          call upcase (astrip)
        end if
        lstrip = (astrip .eq. 'Y')
c
        if (.not. ldump) then
          call fvalut (' Cell :',6,cell)
          write (*,*) 'CCP4 requires CRYST, SCALE and ORIGX cards'
          write (*,*) 'X-PLOR does not like them at all'
          write (*,*) 'Therefore: reply Y for CCP4 and N for X-PLOR :'
          call textin (
     +      ' Write CRYST, SCALE, ORIGX cards (Y/N)         ?',anorem)
          call upcase (anorem)
        end if
        lnorem = (anorem .eq. 'Y')
c
        key = 'ATOM  '
        nwrit = 0
c
        call stamp (line)
        line = ' '//line
        write (f2,4711,err=996) 'REMARK',line(1:leng1(line))
c
        if (length(commnt) .gt. 0) then
          write (f2,4711,err=996) 'REMARK',
     +      (' '//commnt(1:leng1(commnt)))
        end if
c
        if (lrema .and. nrem .gt. 0) then
          do kk=1,nrem
            ldum = length(remark(kk))
            if (ldum .gt. 0) then
              write (f2,'(a)',err=996) remark(kk)(1:ldum)
            end if
 5326       continue
          end do
        end if
c
c ... write CRYST1 etc. cards
c
        if (lnorem) then
c
          write (f2,'(a6,3f9.3,3f7.2,1x,a11,i4)',err=996) 
     +      'CRYST1',(cell(i),i=1,6),spgrp,zmol
c
          write (f2,'(a)',err=996)
     +      'ORIGX1      1.000000  0.000000  0.000000        0.00000'
          write (f2,'(a)',err=996)
     +      'ORIGX2      0.000000  1.000000  0.000000        0.00000'
          write (f2,'(a)',err=996)
     +      'ORIGX3      0.000000  0.000000  1.000000        0.00000'
c
          call orthog (cell,ca2fa,1)
          write (f2,'(a6,4x,3f10.6,f15.5)',err=996)
     +      'SCALE1',(ca2fa(i),i=1,10,3)
          write (f2,'(a6,4x,3f10.6,f15.5)',err=996)
     +      'SCALE2',(ca2fa(i),i=2,11,3)
          write (f2,'(a6,4x,3f10.6,f15.5)',err=996)
     +      'SCALE3',(ca2fa(i),i=3,12,3)
c
        end if
c
        do kk=1,natoms
c
c ...screen chain ID
c
          if (whichn .eq. '**') goto 7352
          if (achain(kk) .ne. whichn) goto 6511
c
c ... screen residue numbers
c
 7352     continue
          if (izone(1) .eq. 0 .and. izone(2) .eq. 0) then
            goto 6262
          else if (izone(1) .eq. 0) then
            if (iresid(kk) .le. izone(2)) goto 6262
            goto 6511
          else if (izone(2) .eq. 0) then
            if (iresid(kk) .ge. izone(1)) goto 6262
            goto 6511
          else
            if (iresid(kk) .ge. izone(1) .and.
     +          iresid(kk) .le. izone(2)) goto 6262
            goto 6511
          end if
c
 6262     continue
c
c ... screen atom types
c
          myres = resnam(kk)
c *** ALL ATOMS
          if (which .eq. 'A') goto 6265
c *** MAIN-CHAIN ATOMS ONLY
          if (which .eq. 'M') then
            if (mainch(atmnam(kk))) goto 6265
          end if
c *** POLY-GLYCINE
          if (which .eq. 'G') then
            myres = 'GLY'
            if (mainch(atmnam(kk))) goto 6265
          end if
c *** POLY-ALANINE
          if (which .eq. 'P') then
            if (myres .eq. 'GLY') then
              if (mainch(atmnam(kk))) goto 6265
            else
              myres = 'ALA'
              if (mainch(atmnam(kk)) .or.
     +            atmnam(kk) .eq. ' CB ') goto 6265
            end if
          end if
c *** POLY-SERINE
          if (which .eq. 'S') then
            if (myres .eq. 'GLY') then
              if (mainch(atmnam(kk))) goto 6265
            else if (myres .eq. 'ALA') then
              if (mainch(atmnam(kk)) .or.
     +            atmnam(kk) .eq. ' CB ') goto 6265
            else
              myres = 'SER'
              if (mainch(atmnam(kk)) .or.
     +            atmnam(kk) .eq. ' CB ') goto 6265
              if (atmnam(kk) .eq. ' CG ' .or.
     +            atmnam(kk) .eq. ' OG ' .or.
     +            atmnam(kk) .eq. ' SG ' .or.
     +            atmnam(kk) .eq. ' CG1' .or.
     +            atmnam(kk) .eq. ' OG1') then
                atmnam (kk) = ' OG '
                goto 6265
              end if
            end if
          end if
          goto 6511
c
 6265     continue
c
c ... screen OT1/OT2
c
          if (.not. lstrip) then
            if (atmnam(kk) .eq. ' OT1') atmnam (kk) = ' O  '
            if (atmnam(kk) .eq. ' OT2') atmnam (kk) = ' OXT'
          end if
c
c ... screen hydrogens
c
          if (.not. lkeeph) then
            if (lhydro(atmnam(kk))) goto 6511
          end if
c
          inum = atomnr(kk)
          if (lforce) inum = nwrit + 1
c
          write (line,4713) inum,atmnam(kk),altloc(kk),
     +      myres,achain(kk),iresid(kk),insert(kk),xf(kk),
     +      yf(kk),zf(kk),qatom(kk),batom(kk),inote(kk)
          write (f2,4711,err=996) key(1:6),(line(1:leng1(line)))
          nwrit = nwrit + 1
c
 6511     continue
        end do
c
        write (f2,4711,err=996) 'END   '
c
        call jvalut (' Nr of atoms written :',1,nwrit)
c
c ... CLOSE THE FILE HERE to make sure it is written out
c     completely !!!
c
        close (f2)
c
      else if (option(1:4) .eq. 'SPLI') then
c
        write (*,*)
        call textin (' Basename of PDB files ?',base)
        nowchn = '$#'
        key = 'ATOM  '
        nline = 0
        lenbas = length (base)
c
        do kk=1,natoms
c
          if (kk .eq. 1 .or. achain(kk) .ne. nowchn) then
c
            if (kk.gt.1) write (f2,4711,err=996) 'END   '
            close (f2)
            file3 = base(1:lenbas)//achain(kk)//'.pdb'
            call remspa (file3)
            call locase (file3)
            if (kk.gt.1) call jvalut (
     +        ' Nr of atoms written to it :',1,nline)
            call textut (' New chain id :',achain(kk))
            call textut (' New pdb file :',file3)
            call xopxua (f2,file3,xinter(),ierr)
            if (ierr .ne. 0) goto 998
            nowchn = achain (kk)
            nline = 0
c
          end if
c
          write (line,4713) atomnr(kk),atmnam(kk),altloc(kk),
     +      resnam(kk),achain(kk),iresid(kk),insert(kk),xf(kk),
     +      yf(kk),zf(kk),qatom(kk),batom(kk),inote(kk)
          write (f2,4711,err=996) key(1:6),(line(1:leng1(line)))
c
          nline = nline + 1
c
        end do
c
        write (f2,4711,err=996) 'END   '
        close (f2)
c
        call jvalut (' Nr of atoms written to it :',1,nline)
        call jvalut (' Nr of atoms written in core :',1,natoms)
c
      else if (option(1:4) .eq. 'RAND') then
c
        write (*,*)
        call rvalin (' Max magnitude for shifts (X/Y/Z/B/Q) ?',
     +    5,shifts)
c
        do i=1,5
          shifts (i) = max (0.0, shifts(i))
          rmses (i) = 0.0
          avers (i) = 0.0
        end do
        rmses (6) = 0.0
        call rvalut (' Maximum shifts :',5,shifts)
c
        call jvalin (' Seed for random number generator     ?',1,iseed)
        call gkrand (dx,0.0,0.0,iseed)
        iseed = 0
c
        do i=1,natoms
c
          dp = 0.0
          do j=1,3
            call gkrand (dx,-shifts(j),shifts(j),0)
            avers (j) = avers (j) + dx
            rmses (j) = rmses (j) + dx**2
            dp = dp + dx**2
            xshf (j) = dx
          end do
          xf (i) = xf (i) + xshf (1)
          yf (i) = yf (i) + xshf (2)
          zf (i) = zf (i) + xshf (3)
          rmses (6) = rmses (6) + sqrt(dp)
c
          if (shifts(4) .gt. 0.0) then
            call gkrand (dx,-shifts(4),shifts(4),0)
            bold = batom (i)
            batom (i) = batom (i) + dx
            batom (i) = max (2.0, batom(i))
            dx = batom(i) - bold
            avers (4) = avers (4) + dx
            rmses (4) = rmses (4) + dx**2
          end if
c
          if (shifts(5) .gt. 0.0) then
            call gkrand (dx,-shifts(5),shifts(5),0)
            bold = qatom (i)
            qatom (i) = qatom (i) + dx
            qatom (i) = min (1.0, max (0.0, qatom(i)))
            dx = qatom(i) - bold
            avers (5) = avers (5) + dx
            rmses (5) = rmses (5) + dx**2
          end if
        end do
c
        do i=1,5
          avers(i) = avers(i)/float(natoms)
          rmses(i) = sqrt (rmses(i)/float(natoms))
        end do
        rmses (6) = rmses(6)/float(natoms)
c
        call rvalut (' Max allowed shifts X/Y/Z/B/Q :',5,shifts)
        call rvalut (' Average     shifts X/Y/Z/B/Q :',5,avers)
        call rvalut (' RMS         shifts X/Y/Z/B/Q :',5,rmses)
        call rvalut (' AVERAGE POSITIONAL SHIFT     :',1,rmses(6))
c
      else if (option(1:4) .eq. 'CA_D') then
c
        write (*,*)
c
        call textin (
     +    ' Name for C-alpha distance plot file ?',file6)
        close (f4)
        call xopxua (f4,file6,xinter(),ierr)
        if (ierr .ne. 0) goto 11
c
        iold = 0
        do i=1,natoms
          if (atmnam(i) .eq. ' CA ') iold = iold + 1
        end do
c
        call jvalut (' Nr of C-alpha atoms  :',1,iold)
        if (iold .lt. 3) then
          call errcon ('Not enough C-alpha atoms')
          close (f4)
          goto 11
        end if
c
        if (iold*iold .gt. maxbuf) then
          call errcon ('Plot data does not fit in buffer')
          close (f4)
          goto 11
        end if
c
        call stamp (line)
c
        write (f4,5000) 'REMARK',
     +    ' C-alpha distance plot'
        write (f4,5000) 'REMARK',(' '//line(1:leng1(line)))
        write (f4,5000) 'REMARK',
     +    (' From PDB file '//file1(1:leng1(file1)))
        write (f4,5000) 'REMARK'
        write (f4,5000) 'XLABEL','Residue'
        write (f4,5000) 'YLABEL','Residue'
        write (f4,5010) 'NLEVEL',7
        write (f4,5000) 'LEVELS'
        write (f4,5000) '6 7 8 ','9 10 11 12'
        write (f4,5000) 'COLOUR'
        write (f4,5000) '1 1 5 ','5 2 6 4'
        write (f4,5010) 'XPOINT',iold
        write (f4,5010) 'YPOINT',iold
        write (f4,5010) 'XLIMIT',1,iold
        write (f4,5010) 'YLIMIT',1,iold
        write (f4,5000) 'ZVALUE','*'
c
        iold = 0
        do i=1,natoms
          if (atmnam(i) .eq. ' CA ') then
            do j=1,natoms
              if (atmnam(j) .eq. ' CA ') then
                iold = iold + 1
                xbuff (iold) = sqrt ( (xf(i)-xf(j))**2 +
     +            (yf(i)-yf(j))**2 + (zf(i)-zf(j))**2 )
              end if
            end do
          end if
        end do
c
        write (f4,'(7f10.1)') (xbuff(i),i=1,iold)
        write (f4,5000) 'END   '
c
        close (f4)
        call jvalut (' Nr of points written :',1,iold)
c
      else if (option(1:4) .eq. 'PLOT') then
c
        write (*,*)
        call textin (' Make plot file for Bs or Qs ?',borq)
        call upcase (borq)
        if (borq .ne. 'Q') borq = 'B'
c
        close (f4)
        close (f5)
c
        call textin (' Filename for per_atom plot    ?',file4)
        call xopxua (f4,file4,xinter(),ierr)
        if (ierr .ne. 0) goto 11
c
        call textin (' Filename for per_residue plot ?',file5)
        call xopxua (f5,file5,xinter(),ierr)
        if (ierr .ne. 0) then
          close (f4)
          goto 11
        end if
c
        write (*,'(99(1x,a,:,/))') ' ',
     +    'You may plot the following for each residue:',
     +    'R = RMS B/Q over all atoms / average over molecule',
     +    'A = average B/Q for all atoms',
     +    'M = average B/Q for main-chain atoms',
     +    'S = average B/Q for side-chain atoms'
        call textin (' Option (R/A/M/S) ?',amain)
        call upcase (amain)
        if (index ('MSAR',amain) .le. 0) amain = 'A'
c
        call textin (
     +    ' Write atom/residue labels to file (Y/N)    ?',alabel)
        call upcase (alabel)
        llabel = (alabel .eq. 'Y')
c
        call prompt (' WARNING - if there are hydrogen atoms they')
        call prompt ('           will be included !')
c
        if (borq .eq. 'B') then
          call xstats (batom,natoms,ave,sdv,xmin,xmax,xtot)
        else
          call xstats (qatom,natoms,ave,sdv,xmin,xmax,xtot)
        end if
c
        call stamp (line)
c
 5000 format (a6,1x,a)
 5010 format (a6,1x,12i6)
 5020 format (a6,1x,6f12.4)
c
c ... do per_atom plot first
c
        if (borq .eq. 'B') then
          write (f4,5000) 'REMARK',
     +      ' Per-atom temperature factor plot'
        else
          write (f4,5000) 'REMARK',
     +      ' Per-atom occupancy plot'
        end if
c
        write (f4,5000) 'REMARK',(' '//line(1:leng1(line)))
        write (f4,5000) 'REMARK',
     +    (' From PDB file '//file1(1:leng1(file1)))
        write (f4,5000) 'REMARK'
        write (f4,5010) 'NPOINT',natoms
        write (f4,5010) 'COLOUR',4
        write (f4,5000) 'XLABEL','Atom'
c
        if (borq .eq. 'B') then
          write (f4,5000) 'YLABEL','Temperature factor'
c
          call xstats (batom,natoms,ave,sdv,xmin,xmax,xtot)
          write (f4,5020) 'XYVIEW',0.0,float(natoms+1),
     +      0.0,(1.05*xmax)
        else
          write (f4,5000) 'YLABEL','Occupancy'
          write (f4,5020) 'XYVIEW',0.0,float(natoms+1),
     +      0.0,1.0
        end if
c
        write (f4,5020) 'XLIMIT',1.0,1.0
        write (f4,5000) 'YVALUE','(8f10.4)'
c
        if (borq .eq. 'B') then
          write (f4,'(8f10.4)') (batom(i),i=1,natoms)
        else
          write (f4,'(8f10.4)') (qatom(i),i=1,natoms)
        end if
c
        if (llabel) then
          write (f4,5000) 'LABELS','(1x,a)'
          do i=1,natoms
            write (line,'(4a,i6,2a)') achain(i),'_',resnam(i),
     +        '_',iresid(i),'_',atmnam(i)
            call remspa (line)
            write (f4,'(1x,a)') line(1:leng1(line))
          end do
        end if
c
        write (f4,5000) 'END   '
        close (f4)
c
c ... now per_residue plot
c
        oldres = -1
        xmax = -999.999
c
        do i=1,natoms
c
cc          if (amain .eq. 'M') then
cc            if (.not. mainch(atmnam(i))) goto 7373
cc          else if (amain .eq. 'S') then
cc            if (mainch(atmnam(i))) goto 7373
cc          end if
c
          ires = iresid (i)
c
          if (oldres .eq. -1) then
c
c ... First residue
c
            nres = 1
            oldres = ires
            nresat = 0
            resptr (nres) = i
            xbuff (1) = 0.0
c
            if (amain .eq. 'M') then
              if (.not. mainch(atmnam(i))) goto 7373
            else if (amain .eq. 'S') then
              if (mainch(atmnam(i))) goto 7373
            end if
c
            nresat = nresat + 1
            if (borq .eq. 'B') then
              if (amain .ne. 'R') then
                xbuff (1) = batom(i)
              else
                xbuff (1) = batom(i)**2
              end if
            else
              if (amain .ne. 'R') then
                xbuff (1) = qatom(i)
              else
                xbuff (1) = qatom(i)**2
              end if
            end if
c
          else if (ires .eq. oldres) then
c
c ... Same residue as before
c
            if (amain .eq. 'M') then
              if (.not. mainch(atmnam(i))) goto 7373
            else if (amain .eq. 'S') then
              if (mainch(atmnam(i))) goto 7373
            end if
c
            nresat = nresat + 1
            if (borq .eq. 'B') then
              if (amain .ne. 'R') then
                xbuff (nres) = xbuff (nres) + batom(i)
              else
                xbuff (nres) = xbuff (nres) + batom(i)**2
              end if
            else
              if (amain .ne. 'R') then
                xbuff (nres) = xbuff (nres) + qatom(i)
              else
                xbuff (nres) = xbuff (nres) + qatom(i)**2
              end if
            end if
c
          else
c
c ... New residue
c
            if (nresat .gt. 0) then
              if (amain .ne. 'R') then
                xbuff (nres) = xbuff (nres) / float (nresat)
              else
                xbuff (nres) = sqrt (xbuff (nres) / float (nresat))
                xbuff (nres) = xbuff (nres) / ave
              end if
              xmax = max (xmax,xbuff (nres))
            else
              xbuff (nres) = 0.0
            end if
c
            nres = nres + 1
            oldres = ires
            nresat = 0
            resptr (nres) = i
            xbuff (nres) = 0.0
c
            if (amain .eq. 'M') then
              if (.not. mainch(atmnam(i))) goto 7373
            else if (amain .eq. 'S') then
              if (mainch(atmnam(i))) goto 7373
            end if
c
            nresat = nresat + 1
            if (borq .eq. 'B') then
              if (amain .ne. 'R') then
                xbuff (nres) = batom(i)
              else
                xbuff (nres) = batom(i)**2
              end if
            else
              if (amain .ne. 'R') then
                xbuff (nres) = qatom(i)
              else
                xbuff (nres) = qatom(i)**2
              end if
            end if
c
          end if
c
 7373     continue
c
        end do
c
        if (nresat .gt. 0) then
          if (amain .ne. 'R') then
            xbuff (nres) = xbuff (nres) / float (nresat)
          else
            xbuff (nres) = sqrt (xbuff (nres) / float (nresat))
            xbuff (nres) = xbuff (nres) / ave
          end if
          xmax = max (xmax,xbuff (nres))
        else
          xbuff (nres) = 0.0
        end if
c
        if (borq .eq. 'B') then
          write (f5,5000) 'REMARK',
     +      ' Per-residue temperature factor plot'
        else
          write (f5,5000) 'REMARK',
     +      ' Per-residue occupancy plot'
        end if
c
        call stamp (line)
        write (f5,5000) 'REMARK',(' '//line(1:leng1(line)))
        write (f5,5000) 'REMARK',
     +    (' From PDB file '//file1(1:leng1(file1)))
c
        if (amain .eq. 'M') then
          write (f5,5000) 'REMARK',
     +      ' Averaged over main-chain atoms'
        else if (amain .eq. 'S') then
          write (f5,5000) 'REMARK',
     +      ' Averaged over side-chain atoms'
        else if (amain .eq. 'R') then
          write (f5,5000) 'REMARK',
     +      ' RMS value over all atoms / average over molecule'
        else if (amain .eq. 'A') then
          write (f5,5000) 'REMARK',
     +      ' Averaged over all atoms'
        end if
c
        write (f5,5000) 'REMARK'
        write (f5,5010) 'NPOINT',nres
        write (f5,5010) 'COLOUR',4
        write (f5,5000) 'XLABEL','Residue'
c
        if (borq .eq. 'B') then
          write (f5,5000) 'YLABEL','Temperature factor'
          write (f5,5020) 'XYVIEW',0.0,float(nres+1),0.0,(1.05*xmax)
        else
          write (f5,5000) 'YLABEL','Occupancy'
          write (f5,5020) 'XYVIEW',0.0,float(nres+1),0.0,1.0
        end if
c
        write (f5,5020) 'XLIMIT',1.0,1.0
        write (f5,5000) 'YVALUE','(8f10.4)'
        write (f5,'(8f10.4)') (xbuff(i),i=1,nres)
c
        if (llabel) then
          write (f5,5000) 'LABELS','(1x,a)'
          do i=1,nres
            j = resptr (i)
            write (line,'(4a,i6)') achain(j),'_',resnam(j),
     +        '_',iresid(j)
            call remspa (line)
            write (f5,'(1x,a)') line(1:leng1(line))
          end do
        end if
c
        write (f5,5000) 'END   '
        close (f5)
c
c ... INVALID OPTION IF HERE
c
      else
c
        call textut (' ERROR - Invalid option :',option)
c
      end if
c
      goto 11
c
c ... all done
c
 9999 continue
      call gkquit
c
      stop
c
  999 write (*,*) ' *** ERROR - while opening input PDB file'
      goto 11
c
  998 write (*,*) ' *** ERROR - while opening output PDB file'
      goto 11
c
  997 write (*,*) ' *** ERROR - while reading input PDB file'
      goto 11
c
  996 write (*,*) ' *** ERROR - while writing output PDB file'
      goto 11
c
  899 write (*,*) ' *** ERROR - while opening input BAD file'
      goto 11
c
  897 write (*,*) ' *** ERROR - while reading input BAD file'
      goto 11
c
ccc  896 write (*,*) ' *** ERROR - while writing output BAD file'
ccc      goto 11
c
      end
