c
c
c
      logical function mainuc (atname)
c
c ... MAINUC (atom_name) - returns .true. if atom_name is a
c                          backbone atom of a DNA or RNA molecule
c                          (i.e., sugar or phosphate)
c
c ... Gerard Kleywegt @ 20001208
c
      implicit none
c
      character*1 quote,star
      parameter (quote = '''', star='*')
c
      character atname*(*),myat*4
c
code ...
c
      myat = atname
c
c ... allow for old and new names for the phosphate oxygens
c
      if (myat .eq. ' P  ' .or.
     +    myat .eq. ' O1P' .or.
     +    myat .eq. ' O2P' .or.
     +    myat .eq. ' OP1' .or.
     +    myat .eq. ' OP2') goto 9000
c
c ... allow for both old names with * and new ones with '
c
      if (myat(4:4) .eq. quote .or.
     +    myat(4:4) .eq. star) then
        if (myat(1:3) .eq. ' O5' .or.
     +      myat(1:3) .eq. ' C5' .or.
     +      myat(1:3) .eq. ' C4' .or.
     +      myat(1:3) .eq. ' O4' .or.
     +      myat(1:3) .eq. ' C3' .or.
     +      myat(1:3) .eq. ' O3' .or.
     +      myat(1:3) .eq. ' C2' .or.
     +      myat(1:3) .eq. ' O2' .or.
     +      myat(1:3) .eq. ' C1') goto 9000
      end if
c
      mainuc = .false.
      return
c
 9000 continue
      mainuc = .true.
      return
      end
