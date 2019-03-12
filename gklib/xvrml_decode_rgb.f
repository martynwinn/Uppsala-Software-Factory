c
c
c
      subroutine xvrml_decode_rgb (r,g,b,irgb)
c
c ... decode RGB value triplet from integer
c
      implicit none
c
      real r,g,b
c
      integer irgb,i,j,k,ipart
c
code ...
c
      call packut (i,j,k,ipart,irgb)
c
      r = max (0.0, min (1.0, float(i+500)/1000.0))
      g = max (0.0, min (1.0, float(j+500)/1000.0))
      b = max (0.0, min (1.0, float(k+500)/1000.0))
c
      return
      end
