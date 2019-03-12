c
c ================================================================
c
      subroutine xps_mapol (rx,ry,zx,zy)
c
c ... convert polar to Cartesian coordinates if necessary
c
      include 'xps.incl'
c
      real rx,ry,zx,zy,qx,qy
c
code ...
c
      if (pspola) then
        qx = -poxmin + max (poxmin, min (rx, poxmax))
        qy = ry
        call fix360 (qy)
        zx = qx * cos (qy*degtor)
        zy = qx * sin (qy*degtor)
      else
        zx = rx
        zy = ry
      end if
c
      return
      end
