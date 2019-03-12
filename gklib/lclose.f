c
c=========================================================================
c
      integer function lclose (xlook4,rarray,ndata,tolmax,differ)
c
      implicit none
c
c     INTEGER FUNCTION LCLOSE (XLOOK4,RARRAY,NDATA,TOLMAX,DIFFER)
c
c ... return index of the element of the array RARRAY (1..NDATA)
c     that is closest to the value of XLOOK4, but not further than
c     TOLMAX; return the index as LCLOSE and the absolute difference
c     in DIFFER
c
c ... ERRORS : LCLOSE = -1 : no element within tolerance
c                       -2 : NDATA < 1
c                       -3 : TOLMAX < 0
c
c ... Gerard J Kleywegt @ 900515
c
      integer nill
      parameter (nill=0)
c
      real    rarray(*),xlook4,differ,tolmax,xdif
c
      integer ndata,itemp,ind
c
code ...
c
      if (ndata.le.nill) then
        lclose = -2
        return
      endif
c
      if (tolmax.lt.0.0) then
        lclose = -3
        return
      endif
c
      itemp  = nill
      differ = abs (rarray(1)-xlook4) + 1.0
c
      do 10 ind = 1, ndata
        xdif = abs (rarray(ind)-xlook4)
        if (xdif.le.tolmax.and.xdif.lt.differ) then
          itemp = ind
          differ = xdif
        endif
   10 continue
c
      if (itemp.le.nill) then
        lclose = -1
      else
        lclose = itemp
      endif
c
      return
      end
