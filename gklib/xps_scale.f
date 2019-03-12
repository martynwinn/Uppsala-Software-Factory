c
c ================================================================
c
      subroutine xps_scale (rxmin,rxmax,rymin,rymax)
c
c ... define scale
c
      include 'xps.incl'
c
      real rxmin,rxmax,rymin,rymax
c
code ...
c
      if (.not. psopen) then
        write (*,*) 'ERROR - No PostScript file is open'
        return
      end if
c
      psxsca = (psxmax-psxmin)/(rxmax-rxmin)
      psysca = (psymax-psymin)/(rymax-rymin)
      psrxin = rxmin
      psryin = rymin
      psrxax = rxmax
      psryax = rymax
c
      return
      end
