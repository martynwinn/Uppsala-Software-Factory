c
c ... xhp_graf.f - HPGL graphics routines
c
c ... 930305 - Gerard Kleywegt - converted Alwyn's routines
c
      SUBROUTINE xhp_init (LUNHP)
c
C ---	INITIATE PLOTTER
C
      include 'xhp_graf.incl'
c
      integer lunhp
c
      BYTE L(3),I(3),a(6),b(7),c(10)
      character*1 L_c(3),I_c(3),a_c(6), b_c(7), c_c(10)
      equivalence(a,a_c), (b,b_c), (c,c_c),(L,L_c) ! mrh
c
!mrh
      DATA L_c/'', '.', '('/
      DATA I_c /'I','N',';'/
c
! mrh
      data a_c/'','.','@',';','2',':'/
      data b_c/'','.','N',';','1','9',':'/
      data c_c/'','.','I','8','1',';',';','1','7',':'/

      L(1) = 27
      a(1) = 27
      b(1) = 27
      c(1) = 27

C
C ---	LOAD IN DEVICE,LUN AND USER->PLOTTER SPACE DEFAULTS
c
      LUN=LUNHP
c
      XSCALE(1)=1.
      YSCALE(1)=1.
      XSCALE(2)=0.
      YSCALE(2)=0.
c
      WRITE(LUN,10)L
      write(lun,10)a
      write(lun,10)b
      write(lun,10)c
      WRITE(LUN,10)I
c
      RETURN
c
  10  FORMAT(1X,10A1)
c
      END
