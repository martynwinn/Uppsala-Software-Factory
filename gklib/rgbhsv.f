c
c
c
      subroutine rgbhsv (r, g, b, h, s, v)
c
c  Translate a colour given in the (R,G,B) triplet into hue, saturation,
c  and value (intensity) as required by the PS300.
c  See Foley & Van Dam p. 615.
c
      implicit none
c
      real r, g, b, h, s, v
c
c When ------- Who ---------------- What -------------------------------
c 10-May-1990  Morten Kjeldgaard    Written, in Dallas.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      real rgbmax, rgbmin, q, rc, gc, bc
c
code ...
c
      rgbmax = max (r,g,b)
      rgbmin = min (r,g,b)
      q = rgbmax - rgbmin
      v = rgbmax
      if (rgbmax .ne. 0.0) then
        s = q/rgbmax
      else
        s = 0.0
      end if
      if (s .eq. 0.0) then
        h = 0.0
      else
c ---    rc measures the distance of color from red
        rc = (rgbmax - r)/q
        gc = (rgbmax - g)/q
        bc = (rgbmax - b)/q
        if (r .eq. rgbmax) then
c ---       resulting color between yellow and magenta
          h = bc - gc
        else if (g .eq. rgbmax) then
c ---       resulting color between cyan and yellow
          h = 2.0 + rc - bc
        else if (b .eq. rgbmax) then
c ---       resulting color between magenta and cyan
          h = 4.0 + gc - rc
        end if
c ---      convert to degrees
        h = h * 60.0 + 120.0
        if (h .lt. 0.0) h = h + 360.0
        if (h .gt. 360.0) h = h - 360.0
      end if
c
      return
      end
