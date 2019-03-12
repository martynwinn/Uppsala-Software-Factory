c
c
c
      subroutine xvrml_col_list ()
c
      include 'xvrml.incl'
c
      real r,g,b
c
      integer i
c
code ...
c
      call jvalut (' Nr of colours :',1,numcol)
      do i=1,numcol
        call xvrml_decode_rgb (r,g,b,colnum(i))
        write (*,6000) i,colnam(i),colnum(i),r,g,b
      end do
c
 6000 format (' # ',i3,' (',a,') = ',i12,' RGB ',3f8.3)
c
      return
      end
