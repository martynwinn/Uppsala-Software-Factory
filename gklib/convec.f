c
c
c
      subroutine convec (f1,f2,f3,x,a)
c
c ... convert cartesian to fractionals
c
      implicit NONE
c
      real f1,f2,f3,x(3),a(3,3),t(3)
c
code ...
c
      t(1) = f1
      t(2) = f2
      t(3) = f3
      call mulmtx (a,t,x,3,3,1)
c
      return
      end
