c
c
c
      real function betai (a,b,x)
c
c ... BETAI - Numerical Recipes, p. 167
c
      implicit none
c
      real a,b,x
      real bt
      real gammln,betacf
c
code ...
c
      if (x .lt. 0.0 .or. x .gt. 1.0) then
        betai = 0.0
        call errcon ('Invalid argument in BETAI !')
        return
      end if
c
      if (x .eq. 0.0 .or. x .eq. 1.0) then
        bt = 0.0
      else
        bt = exp(gammln(a+b)-gammln(a)-gammln(b) +
     +           a*alog(x)+b*alog(1.0-x))
      end if
c
      if (x .lt. (a+1.0)/(a+b+2.0)) then
        betai = bt * betacf(a,b,x)/a
        return
      else
        betai = 1.0 - bt*betacf(b,a,1.0-x)/b
        return
      end if
c
      return
      end
