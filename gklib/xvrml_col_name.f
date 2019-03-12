c
c
c
      subroutine xvrml_col_name (col,r,g,b,icol)
c
      include 'xvrml.incl'
c
      real r,g,b
c
      integer i,icol,ll,length
c
      character col*(*),mycol*21
c
code ...
c
      mycol = col
      call locase (mycol)
      call subchr (mycol,' ','_',i)
      ll = length (mycol)
      if (ll .le. 1) then
        icol = colnum (1)
        call xvrml_decode_rgb (r,g,b,icol)
        return
      end if
c
      do i=1,numcol
        if (mycol .eq. colnam(i)) then
          icol = colnum (i)
          call xvrml_decode_rgb (r,g,b,icol)
          return
        end if
      end do
c
      do i=1,numcol
        if (index (colnam(i),mycol(1:ll)) .gt. 0) then
          icol = colnum (i)
          call xvrml_decode_rgb (r,g,b,icol)
          return
        end if
      end do
c
      icol = colnum (1)
      call xvrml_decode_rgb (r,g,b,icol)
c
      return
      end
