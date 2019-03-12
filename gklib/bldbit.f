c
c
c
      subroutine bldbit (map, e1, e2, e3, orgn, grid, cell,
     $                   min, max, rt, ct, bit, x, errcod)
c
c ---	If we were unable to place point in asymm. unit, get the
c	grid points around it within 4x4x4
c ---	this is needed for grid point between edge of extent
c	and next integral point.
c
      implicit none
c
      integer e1, e2, e3, orgn(3), ct, errcod
      real map(e1, e2, e3), grid(3), cell(6), min(3), max(3)
      real rt(12,1), bit(4,4,4), x(3)
c
      integer i, ig(3), ix(3), j, k, l
      real x1(3)
c
code ...
c
      do 110 l=1,3
110     ig(l) = (x(l)*cell(l)/grid(l)-1.)
c
      do 100 i=1,4
      do 100 j=1,4
      do 100 k=1,4
        ix(1) = ig(1)+ i-1
        ix(2) = ig(2)+ j-1
        ix(3) = ig(3)+ k-1
c
c ---	ix is the grid value of the bitarray (i,j,k)
c ---	convert to fractionals again and off we go
c
        do 120 l=1,3
120       x1(l) = float(ix(l))*grid(l)/cell(l)
c
        call frcsym (x1, min, max, rt, ct, errcod)
        if (errcod .ne. 0) return
c
c ---	this should now be in the volume and on a grid point
c
          do 130 l=1,3
130         ix(l) = nint(x1(l)*cell(l)/grid(l))- orgn(l)+ 1
c
          bit(i,j,k) = map(ix(1), ix(2), ix(3))
c
100   continue
      errcod = 0
c
      return
      end
