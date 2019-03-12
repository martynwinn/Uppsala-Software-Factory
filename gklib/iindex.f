c
c=========================================================================
c
      integer function iindex (ix, itoler, ndata, idata)
c
      implicit none
c
c=========================================================================
c
c FUNCTION: IINDEX
c ----------------
c
c VERSION : 881031/890126
c PURPOSE : check whether an integer number IX is equal to one of the NDATA
c           elements of the integer array IDATA, given a tolerance ITOLER
c           if there exists an IDATA (M) such that |IX - IDATA(M)| <= ITOLER
c           then IINDEX is set to M, otherwise IINDEX = -1
c           if NDATA < 1 then IINDEX = -2
c           if ITOLER < 0 then IINDEX = -3
c
c CALL    : IINDEX(IX,ITOLER,NDATA,IDATA)
c
c AUTHOR  : G J Kleywegt
c           Laboratory of Organic Chemistry
c           Padualaan 8 , 3584 CH Utrecht, The Netherlands
c
c=========================================================================
c
      integer nill
      parameter (nill=0)
c
      integer idata(*), i, itoler, ix, ndata
c
code ...
c
      if (ndata.lt.1) then
        iindex = -2
        return
      endif
c
      if (itoler.lt.nill) then
        iindex = -3
        return
      endif
c
c ... first try a direct hit
c
      do 10 i=1,ndata
        if (idata(i).eq.ix) goto 20
  10  continue
c
c ... if this fails, invoke the tolerance (if any)
c
      if (itoler.gt.nill) then
        do 15 i=1,ndata
          if (iabs(idata(i)-ix).le.itoler) goto 20
  15    continue
      end if
c
c ... ix does not occur in idata
c
      iindex = -1
      return
c
c ... ix is the i-th element of idata
c
  20  iindex = i
      return
      end
