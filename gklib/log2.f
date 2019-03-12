c
c
c
      real function log2 (a)
c
c ... This function computes the base 2 logarithm of A
c     http://www.cs.utah.edu/~hamlet/release/lessons/fortran08/examples/user2.f#user2f
c
      implicit none
c
      real a
c
code ...
c
      log2 = log10(a)/log10(2.0)
c
      return
      end
