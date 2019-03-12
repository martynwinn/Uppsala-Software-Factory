c
c --------------------------------------------------------------------------
c
      logical function mainch (atname)
c
c ... return .TRUE. if atom name (char. 1-4) are recognised
c     as a main-chain atom's name (incl XPLOR's OT1/OT2)
c     950330 - include OTX and OXT and PROLSQ OT
c
c ... Gerard Kleywegt
c
      implicit none
c
      character atname*(*)
c
code ...
c
      mainch = (atname (1:4) .eq. ' C  ' .or.
     +          atname (1:4) .eq. ' O  ' .or.
     +          atname (1:4) .eq. ' OT1' .or.
     +          atname (1:4) .eq. ' OT2' .or.
     +          atname (1:4) .eq. ' OTX' .or.
     +          atname (1:4) .eq. ' OXT' .or.
     +          atname (1:4) .eq. ' OT ' .or.
     +          atname (1:4) .eq. ' CA ' .or.
     +          atname (1:4) .eq. ' N  ')
c
      return
      end
