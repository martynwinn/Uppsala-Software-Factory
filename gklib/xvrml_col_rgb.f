c
c
c
      subroutine xvrml_col_rgb (i,r,g,b)
c
c ... used by PostScript routines
c
      include 'xvrml.incl'
c
      real r,g,b
c
      integer i,j,k
c
code ...
c
      if (i .gt. numcol) then
        call xvrml_decode_rgb (r,g,b,i)
      else
        k = max (1, min (numcol, i))
        call xvrml_decode_rgb (r,g,b,colnum(k))
      end if
c
      return
      end
