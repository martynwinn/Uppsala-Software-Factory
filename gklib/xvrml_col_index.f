c
c
c
      subroutine xvrml_col_index (i,j)
c
      include 'xvrml.incl'
c
      integer i,j,k
c
code ...
c
      k = max (1, min (numcol, i))
      j = colnum (k)
ccc      print *,i,j,k,colnam(k)
c
      return
      end
