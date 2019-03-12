c
c
c
      subroutine xvrml_rgb_name (colnam,r,g,b)
c
      implicit none
c
      real r,g,b
c
      integer j
c
      character colnam*(*),mycol*80
c
code ...
c
      mycol = colnam
c
      read (mycol,*,err=10) r,g,b
      call xvrml_encode_rgb (r,g,b,j)
      call xvrml_decode_rgb (r,g,b,j)
      return
c
   10 continue
      call xvrml_col_name (mycol,r,g,b,j)
      return
c
      end
