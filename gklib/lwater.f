c
c ==========================================================================
c
      logical function lwater (restyp)
c
c ... LWATER - returns .true. if residue type indicates WATER
c
      implicit none
c
      character rnam*3,restyp*(*)
c
code ...
c
      rnam = restyp
      call upcase (rnam)
c
      lwater = (rnam.eq.'WAT' .or. rnam.eq.'HOH' .or.
     +          rnam.eq.'SOL' .or. rnam.eq.'OHH' .or.
     +          rnam.eq.'HHO' .or. rnam.eq.'H2O' .or.
     +          rnam.eq.'OH2' .or. rnam.eq.'H3O' .or.
     +          rnam.eq.'OH3' .or. rnam.eq.'EAU' .or.
     +          rnam.eq.'DOD' .or. rnam.eq.'TIP' .or.
     +          rnam.eq.'D2O' .or. rnam.eq.'ODD')
c
      return
      end
