c
c
c
      logical function lsatom (atom1,atom2)
c
c ... LSATOM - returns .TRUE. if same atom name (modulo * or ')
c              (used for nucleic acids)
c
c ... Gerard J Kleywegt @ 001212
c
      implicit none
c
      character*1 quote,star
      parameter (quote = '''', star='*')
c
      character atom1*4,atom2*4
c
code ...
c
      if (atom1 .eq. atom2) goto 1000
c
      if (atom1(1:3) .eq. atom2(1:3)) then
        if (atom1(4:4) .eq. quote .and.
     +      atom2(4:4) .eq. star) goto 1000
        if (atom2(4:4) .eq. quote .and.
     +      atom1(4:4) .eq. star) goto 1000
      end if
c
      lsatom = .false.
      return
c
 1000 continue
      lsatom = .true.
c
      return
      end
