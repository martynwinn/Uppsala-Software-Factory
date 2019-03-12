c
c ... xalwyn.f - subroutines taken from alwyn
c
c ===========================================================================
c
      subroutine packin (H,K,L,PART,REFL)
c
      implicit none
c
      INTEGER H,K,L,PART,REFL,I1,I2,I3,I4
c
      parameter (i1 = 512)
      parameter (i2 = 2*i1)
      parameter (i3 = i2*i2)
      parameter (i4 = i3*i2)
c
ccc      DATA I1,I2,I3,I4/256,512,262144,134217728/
c
      REFL = (H+I1) + (K+I1)*I2 + (L+I1)*I3 + PART*I4
c
      RETURN
      END
