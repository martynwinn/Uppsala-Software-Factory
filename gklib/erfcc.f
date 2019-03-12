c
c
c
      real function erfcc (x)
c
c ... from "Numerical Recipes", p. 164 (chapter 6.2)
c
      real x,z,t
c
code ...
c
      z=abs(x)
      t=1.0/(1.0+0.5*z)
      erfcc = t*exp(-z*z - 1.26551223 + t*(1.00002368 +
     +  t*(0.37409196 + t*(0.09678418 + t*(-0.18628806 +
     +  t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
     +  t*(-0.82215223 + t*0.17087277)))))))))
      if (x.lt.0.0) erfcc = 2.0 - erfcc
c
      return
      end
