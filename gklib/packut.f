c
c ===========================================================================
c
      subroutine packut (H,K,L,PART,REFL)
c
      implicit none
c
      INTEGER H,K,L,PART,REFL,I1,I2,I3,I4,I5,I6,I7
c
      parameter (i1 = 512)
      parameter (i2 = 2*i1 - 1)
      parameter (i3 = 2*i1)
      parameter (i4 = i2*i3)
      parameter (i5 = i3*i3)
      parameter (i6 = i4*i3)
      parameter (i7 = i5*i3)
c
cc      DATA I1,I2,I3,I4,I5,I6,I7
cc     $	/256,511,512,261632,262144,133955584,134217728/
c
      H = IAND(REFL,I2)    - I1
      K = IAND(REFL,I4)/I3 - I1
      L = IAND(REFL,I6)/I5 - I1
      PART = REFL/I7
c
      RETURN
      END
