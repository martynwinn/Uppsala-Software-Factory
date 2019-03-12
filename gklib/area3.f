c
c
c
      real function area3 (x1,x2,x3)
c
c ... AREA3 = calculate area of triangle using Hero's formula
c     calculate the three sides a, b, c
c     let s = (a+b+c) / 2
c     then area = sqrt (s*(s-a)*(s-b)*(s-c))
c
c ... Gerard Kleywegt @ 2001-07-27
c
      implicit none
c
      real x1(3),x2(3),x3(3)
c
      real a,b,c,s,distce
c
code ...
c
      a = distce(x1,x2)
      b = distce(x1,x3)
      c = distce(x2,x3)
      s = 0.5 * (a+b+c)
      area3 = sqrt (s*(s-a)*(s-b)*(s-c))
c
      return
      end
