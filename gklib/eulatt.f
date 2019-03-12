c
c ===========================================================================
c
      subroutine eulatt (ang,mat)
c
c ... convert Lattman angles THETA+, THETA2, THETA-
c     into a rotation matrix
c
      double precision twopi,degtor
      parameter (twopi=6.2831853071796)
      parameter (degtor=twopi/360.0)
c
      real ang(3),mat(3,3)
c
      double precision tp,t2,tm,cp,sp,c2,s2,cm,sm,cc1,cc2,ttp,ttm,dhalf
c
code ...
c
      dhalf = 0.5D0
c
      tp = ang(1)*degtor
      t2 = ang(2)*degtor
      tm = ang(3)*degtor
c
      cp = cos(tp)
      sp = sin(tp)
      c2 = cos(t2)
      s2 = sin(t2)
      cm = cos(tm)
      sm = sin(tm)
c
      cc1 = 1.0D0 + c2
      cc2 = 1.0D0 - c2
c
      ttp = dhalf*(tp+tm)
      ttm = dhalf*(tp-tm)
c
      mat (1,1) = dhalf * ( cc2*cm + cc1*cp )
      mat (2,1) = dhalf * ( cc2*sm - cc1*sp )
      mat (3,1) = sin(ttp)*s2
c
      mat (1,2) = dhalf * ( cc2*sm + cc1*sp )
      mat (2,2) = dhalf * ( cc1*cp - cc2*cm )
      mat (3,2) = -cos(ttp)*s2
c
      mat (1,3) = sin(ttm)*s2
      mat (2,3) = cos(ttm)*s2
      mat (3,3) = c2
c
      return
      end
