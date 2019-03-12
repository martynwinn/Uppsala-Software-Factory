c
c ========================================================================
c
      logical function lhydro (atmnam)
c
      implicit none
c
      integer i1,i2,ia,iz
c
      character atmnam*(*), mynam*4
c
code ...
c
      mynam = atmnam
      call upcase (mynam)
c
c ... default PDB names: ' H' in columns 1 and 2
c
      if (mynam (1:2) .eq. ' H') goto 6511
c
      if (mynam (1:2) .eq. 'H ') goto 6511
c
      i1 = ichar (mynam(1:1))
      i2 = ichar (mynam(2:2))
      ia = ichar ('A')
      iz = ichar ('Z')
c
c ... allow names like 1H, 'H, H5, etc.
c
      if (mynam (1:1) .eq. 'H' .and.
     +    (i2 .lt. ia .or. i2 .gt. iz)) goto 6511
c
      if (mynam (2:2) .eq. 'H' .and.
     +    (i1 .lt. ia .or. i1 .gt. iz)) goto 6511
c
c ... allow for X-PLOR idiosyncrasies (Arg-HH, Asn-HD, Gln-HE)
c
      if (mynam (1:2) .eq. 'HH') goto 6511
      if (mynam (1:2) .eq. 'HD') goto 6511
      if (mynam (1:2) .eq. 'HE') goto 6511
c
      if (mynam (1:2) .eq. 'HM') goto 6511
c
c ... allow things like " 1H " etc.
c
      if (mynam(1:1) .eq. ' ' .and. mynam(3:3) .eq. 'H' .and.
     +    (i2 .lt. ia .or. i2 .gt. iz)) goto 6511
c
c ... other NMR inventions (' Qxx' pseudo-atoms)
c
      if (mynam (1:2) .eq. ' Q') goto 6511
c
c ... now check for Deuterium ...
c
c ... default PDB names: ' D' in columns 1 and 2
c
      if (mynam (1:2) .eq. ' D') goto 6511
c
      if (mynam (1:2) .eq. 'D ') goto 6511
c
      lhydro = .false.
      return
c
 6511 continue
      lhydro = .true.
c
      return
      end
