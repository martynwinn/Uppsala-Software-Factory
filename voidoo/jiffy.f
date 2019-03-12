      program jiffy
c
      real twopi
      parameter (twopi=6.2831853071796)
c
code ...
c
      call gkinit ('JIFFY','961213')
c
      vol = 30.0
      third = 1.0 / 3.0
c
   10 continue
      print *
      call fvalin (' Volume ?',1,vol)
      if (vol .le. 0.0) goto 20
      rad = (3.0 * vol) / (2.0 * twopi)
      rad = rad ** third
      call fvalut (' Radius :',1,rad)
      goto 10
c
   20 continue
      end
