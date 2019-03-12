c
c
c
      subroutine xvrml_encode_rgb (r,g,b,irgb)
c
c ... encode RGB value triplet into one integer
c
      implicit none
c
      real r,g,b
c
      integer irgb,i,j,k
c
code ...
c
      i = max(-500.0, min (500.0,1000.0*(r-0.5)))
      j = max(-500.0, min (500.0,1000.0*(g-0.5)))
      k = max(-500.0, min (500.0,1000.0*(b-0.5)))
c
      call packin (i,j,k,0,irgb)
c
      return
      end
