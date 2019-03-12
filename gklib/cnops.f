c
c=====================================================
c
      logical function cnops (name)
c
      implicit none
c
      character name*4,two*2
c
code ...
c
      two = name(1:2)
      call upcase (two)
c
      cnops = (two.eq.' C' .or. two.eq.' N' .or.
     +         two.eq.' O' .or. two.eq.' P' .or.
     +         two.eq.' S')
c
      return
      end
