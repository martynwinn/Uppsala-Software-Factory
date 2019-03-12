c
c --------------------------------------------------------------------------
c
      subroutine detaj (aresid,achain,iresid)
c
c ... deduce chain id and residue nr from an Alwyn-type of
c     character residue id [i.e., go from (A6) to (A2,I4)]
c
c     example: ARESID = 'AA132'  -> ACHAIN = 'AA' IRESID =  132
c                       'W43'    ->          ' W'            43
c                       '191'    ->          '  '           191
c                       '12A'    ->          ' A'            12
c                       '1292X2' ->          'X2'          1292
c
c ... Gerard Kleywegt
c
      implicit none
c
      integer iresid,ll,length
c
      character aresid*6,achain*2,esline*128
c
code ...
c
c ... remove spaces and left-justify text
c
      call remspa (aresid)
      ll = length (aresid)
c
c ... first guess: no chain id (legal PDB)
c
      achain = '  '
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
      achain = ' '//aresid(1:1)
      read (aresid(2:6),*,err=20) iresid
      goto 200
c
c ... third guess: two-character chain id
c     (this is an extension to "pure PDB")
c
   20 continue
      achain = aresid(1:2)
      read (aresid(3:6),*,err=30) iresid
      return
c
c ... fourth guess: one-character ID APPENDED to residue nr
c
   30 continue
      achain = aresid (ll:ll)
      read (aresid(1:ll-1),*,err=40) iresid
      goto 200
c
c ... fifth guess: two-character ID APPENDED to residue nr
c     (this is an extension to "pure PDB")
c
   40 continue
      achain = aresid (ll-1:ll)
      read (aresid(1:ll-2),*,err=50) iresid
      return
c
c ... final guess: "A123B"-type
c     (this is an extension to "pure PDB")
c
   50 continue
      achain = aresid (1:1)//aresid(ll:ll)
      read (aresid(2:ll-1),*,err=60) iresid
      return
c
c ... if here, then PDB convention grossly violated
c
   60 continue
      esline = 'While converting Alwyn-type residue id: '
     +         //aresid
      call errcon (esline)
      achain = '??'
      iresid = 9999
      return
c
c ... case of mixed alpha-numeric chain ID
c
  200 continue
      if (iresid .gt. 9999) then
        write (aresid,'(a1,i5)') achain(2:2),iresid
        achain = aresid(1:2)
        read (aresid(3:6),*) iresid
      end if
c
      return
      end
