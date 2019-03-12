c
c
c
      subroutine celvol (cell,volume)
c
c ... CELVOL - calculate unit-cell volume
c
      implicit NONE
c
      real twopi,rtodeg,degtor
      parameter (twopi  = 6.2831853071796)
      parameter (rtodeg = 360.0 / twopi)
      parameter (degtor = twopi / 360.0)
c
      real cell(6),volume,ca,cb,cg
c
code ...
c
      ca = cos(cell(4)*degtor)
      cb = cos(cell(5)*degtor)
      cg = cos(cell(6)*degtor)
c
      volume = cell(1)*cell(2)*cell(3) *
     +  sqrt (1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg)
c
      return
      end
