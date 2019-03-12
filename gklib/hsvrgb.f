c
c
c
      subroutine hsvrgb (qh, s, v, r, g, b)
c
c  Translate a colour given in hue, saturation and value (intensity) as
c  given by the PS300 into an (R,G,B) triplet.
c  This routine is not used in O, it is just here because I wrote it 
c  for the conversion of the program to rgb colours.
c  See Foley & Van Dam p. 615.
c
      implicit none
c
      real qh, s, v, r, g, b
c
c When ------- Who ---------------- What -------------------------------
c 11-May-1990  Morten Kjeldgaard    Written, in Dallas.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      real f, p, q, t, h
      integer i
c
code ...
c
      h = qh
      if (s .eq. 0.0) then
        r = v
        g = v
        b = v
        return
      end if
      h = h - 120.0
      if (h .lt. 0.0) h = h + 360.0
      if (h .ge. 360.0) h = h - 360.0
      h = h/60.
      i = int(h)
      f = h - float(i)
      p = v * (1.0 - s)
      q = v * (1.0 - (s*f))
      t = v * (1.0 - (s*(1.0 - f)))
      if (i .eq. 0) then
        r = v
        g = t
        b = p
      else if (i .eq. 1) then
        r = q
        g = v
        b = p
      else if (i .eq. 2) then
        r = p
        g = v
        b = t
      else if (i .eq. 3) then
        r = p
        g = q
        b = v
      else if (i .eq. 4) then
        r = t
        g = p
        b = v
      else if (i .eq. 5) then
        r = v
        g = p
        b = q
      end if
c
      return
      end
