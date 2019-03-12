c
c ===========================================================================
c
      real function det3 (rm)
c
c ... calculate determinant of matrix RM(3,3)
c
      implicit none
c
      real rm(3,3)
c
c
code ...
c
      det3 =   rm(1,1) * (rm(2,2)*rm(3,3) - rm(3,2)*rm(2,3))
     +       - rm(2,1) * (rm(1,2)*rm(3,3) - rm(3,2)*rm(1,3))
     +       + rm(3,1) * (rm(1,2)*rm(2,3) - rm(2,2)*rm(1,3))
c
      return
      end
